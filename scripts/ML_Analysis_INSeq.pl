#!/usr/bin/perl -w
#--------------------------------------------------------------------------#
#  Author: Meng Wu, the Gordon Lab, Washington University in St. Louis     #
#                                                                          #
#  File: ML_Analysis_INSeq.pl						   #
#  Descriptions: Fitness index calculation based on EM		           #
#									   #
#  Usage:								   #
# 	see Usage in the code	                                           #
#					                                   #
#  Contact: mengwu@wustl.edu						   #
#--------------------------------------------------------------------------#

use strict;
use warnings;
use Getopt::Long;
#use Statistics::Distributions;

if (!@ARGV){
    print "Usage:$0 <logratio file> <counts file>\n";
    exit;
}


my $qlogratio;
my $file = shift@ARGV;
my $countfile=shift @ARGV;
open(IN, $file) || die "cannot open the logratio file: $file\n";
my $total;
my $header = <IN>;
chomp $header;
my @fileNames = split/\s+/, $header;
shift @fileNames;

while (<IN>){
   chomp;
   my @tmp=split /\s+/;
   my $geneID = shift@tmp;
   for(my $n=0; $n<= $#tmp; $n++)
   {
	$qlogratio -> {$fileNames[$n]} -> {$geneID} = $tmp[$n];  
   }
   $total ++;
}
close IN;


SAMPLE: foreach my $sample (sort keys %{$qlogratio}){
my $uA;
my $uB;
my $u_null;
my $stdev_null;
my $q = 0.5;
my $stdev;
my $iteration = 1000;
my $i;
my $diff;
my @populationtypes = qw(aa bb);
my %output;
my $solve;
my %likelihood_aa;
my %likelihood_bb;
my %likelihood;
 
my $out=$sample.".geno";
open(OUT, ">$out") || die "cannot open $out\n";


my @data=();

foreach my $geneKey(sort keys %{$qlogratio -> {$sample}})
{
	push @data, $qlogratio->{$sample}->{$geneKey};
}

my @meansd = &quantiles_std(@data);
	$uA = $meansd[0];
	$uB = $meansd[2];
 	$stdev = $meansd[3];
        ($u_null,$stdev_null)=&stddev(@data);


for ($i=0;$i<=$iteration;$i++)
{
	my $groups;
	my $population;
	my @diffs;
	my $Nq;
        my $Nstdev;
        my $mean;
	my %num;

	my @parameters = ($uA,$uB,$q, $stdev);
	my @newparam;

	for(my $j=0; $j<= $#populationtypes; $j++)
	{
		$mean -> {$populationtypes[$j]} = $parameters[$j];
		$num{$populationtypes[$j]} = 0;
	}

	foreach my $k (sort keys%{$qlogratio -> {$sample}})
	{
		($population,$likelihood_aa{$k},$likelihood_bb{$k},$likelihood{$k})=&E($qlogratio->{$sample}->{$k}, @parameters);  		
		$output{$k} = $population;
                next SAMPLE if ($likelihood_aa{$k} eq "NA");
		push @{$groups->{$population}},$qlogratio->{$sample} ->{$k};
	}


	foreach my $populationtype (@populationtypes)
	{
		if(defined(@{$groups->{$populationtype}}))
		{
			my ($g_ave, $g_stdev, $g_num) = &ave_var(@{$groups->{$populationtype}});
			$mean -> {$populationtype} = $g_ave;		
			$Nstdev += $g_stdev;
			$num{$populationtype} = $g_num;
		}
		push @newparam, $mean -> {$populationtype};
	}


	$Nstdev /= $total;
	$Nstdev = sqrt($Nstdev);				
	
	my $Nbb = $num{$populationtypes[1]};

	$Nq=$Nbb/$total;			

	push @newparam, $Nq;
	push @newparam, $Nstdev;

	for(my $index=0; $index <= $#parameters; $index++)
	{
		push @diffs, abs($parameters[$index] - $newparam[$index]);
	}

	$diff=&max(@diffs);

#	print SUM "$i\t$uA\t$uB\t$q\t$stdev\t$diff\n";
	  
	if ($diff<10e-10)
	{
		$solve++;
		if($solve > 2)
		{
			last;
		}	       					
	}
	else
	{			 				
		$uA=$newparam[0];	     			 
	      	$uB=$newparam[1];
	      	$q=$newparam[2];
	      	$stdev=$newparam[3];
	}
}

$stdev_null=$stdev;
my $likelihood_null_all=0;
my $likelihood_alter_all=0;
my $likelihood_null;
foreach my $k (sort keys%{$qlogratio -> {$sample}})
        {
            $likelihood_null->{$k}=&likelihood($qlogratio->{$sample}->{$k},$u_null,$stdev_null);
            $likelihood_null_all+=$likelihood_null->{$k};
            $likelihood_alter_all+=$likelihood{$k};
        }

my $LR;  
$LR=-2*$likelihood_null_all+2*$likelihood_alter_all;


print OUT "ID\tlogratio\tGeno\tLikelihood_aa\tLikelihood_bb\tLikelihood\tLikelihood_null\n";
foreach my $key (sort {$a cmp $b} keys %output)
{
	my $idkey = $key;
	my $tmpqlogratio = $qlogratio -> {$sample} -> {$idkey};
	print OUT "$idkey\t$tmpqlogratio\t$output{$key}\t$likelihood_aa{$key}\t$likelihood_bb{$key}\t$likelihood{$key}\t$likelihood_null->{$key}\n";
}

close OUT;


foreach my $populationtype (@populationtypes)
{
	my $outgenofile = $sample."-".$populationtype.".out";
	open(GEN, ">$outgenofile") || die "cannot open $outgenofile\n";

	open(IN2, $out) || die "cannot open $out\n";
        <IN2>;
	while(<IN2>)
	{
		chomp;
		if($_ =~ /$populationtype/)
		{
			print GEN "$_\n";
		}
	}
	close IN2;
	close GEN;
}


my $Rcode = "R.script";

open(R, ">$Rcode") || die "cannot open $Rcode\n";

print R "
 counts_all<-read.table(\"$countfile\",sep=\"\\t\",header=T,row.names=1)
 refcount<-counts_all[,1]
 samplecount<-counts_all[[\"$sample\"]]
 data_all<-read.table(\"$sample.geno\",sep=\"\\t\",header=T)
 all <- data_all\$logratio
 pall <-  hist(all,breaks=60)
 data <- read.table(\"$sample-aa.out\", sep=\"\\t\")
 aa<-data\$V2
 c1 <- ceiling((max(aa) - min(aa))/diff(pall\$mids[1:2]))
 p1 <- hist(aa,breaks=c1)
 data <- read.table(\"$sample-bb.out\", sep=\"\\t\")
 bb<-data\$V2
 c3 <- ceiling((max(bb) - min(bb))/diff(pall\$mids[1:2]))
 p3 <- hist(bb,breaks=c3)
 lb <- min(min(aa), min(bb)) - 0.5
 ub <- max(max(aa), max(bb)) + 0.5
 pdf(\"$sample-MLAnalysisPlot.pdf\")
 hist(all,breaks=60, main=\"$sample\", xlab=\"\")
 xfit<-seq(lb,ub,length=400)
 yfit<-dnorm(xfit,mean=mean(aa),sd=sd(aa))
 yfit <- yfit*diff(p1\$mids[1:2])*length(aa)
 lines(xfit, yfit, col=\"red\", lwd=2)
 xfit<-seq(lb,ub,length=400)
 yfit<-dnorm(xfit,mean=mean(bb),sd=sd(bb))
 yfit <- yfit*diff(p3\$mids[1:2])*length(bb)
 lines(xfit, yfit, col=\"green\", lwd=2)
 zscore<-(all-mean(bb))/sd(bb)
 pvalue<-1-pnorm(abs(zscore))
 padj<-p.adjust(pvalue,method=\"fdr\")
 total<-cbind(data_all,zscore,pvalue,padj,refcount,samplecount)
 write.table(total,\"$sample-pvalue.txt\",sep=\"\\t\")
 dev.off()
 q()
";

system ("R < $Rcode --vanilla");

my $tmpfile = "Rplots.pdf";
unlink $tmpfile;
#unlink $Rcode;

}

sub E
{
        use constant PI    => 4 * atan2(1, 1);
	my ($E_pheno, $E_uAA, $E_uBB, $E_q, $E_stdev) = @_;
	my $egeno;
        my $lnaa;
        my $lnbb;
        my $ln;

        if (($E_q!=1)&&($E_q!=0)) {
	  $lnaa = -log($E_stdev) - 0.5*log(2*PI) - 0.5*(($E_pheno-$E_uAA)**2)/($E_stdev**2) + log(1-$E_q);
	  $lnbb = -log($E_stdev) - 0.5*log(2*PI) - 0.5*(($E_pheno-$E_uBB)**2)/($E_stdev**2) + log($E_q);
       
	  if ($lnaa>$lnbb)
	  {
		$egeno="aa";
	  } 
	  else 
	  {
      		$egeno="bb";
   	  }
          if ((exp($lnaa)+exp($lnbb))!=0){
            $ln=log(exp($lnaa)+exp($lnbb));
          }else{
            $ln="AA";
          }
   	  return $egeno,$lnaa,$lnbb,$ln;
          }else {
              print "Only one distribution can be estimated\n";
              $egeno="bb";
              $lnaa="NA";
              $lnbb="NA";
              $ln="AA"
          }
        
}


sub quantiles_std
{
        my @vector = @_;
        my %p;
        my %pkey;
        my $index = 1;
	my $totalnum = $#vector + 1;

        foreach my $key (sort{$a <=> $b} @vector)
        {
                my $quantile = ($index-1)/($totalnum-1);	
                $p{$key} = $quantile;
                $pkey{$quantile} = $key;
                $index ++;
        }

        my @quans = (0.05, 0.5, 0.95);			
        my @lower;
        my @upper;
        my @quant;

        foreach my $q (@quans)
        {
                foreach my $order(values %p)
                {
                        if($order <= $q)
                        {
                                push @lower, $order;
                        }
                	if($order > $q)
                        {
                                push @upper, $order;
                        }
                }

                my $qlb = &max(@lower);				
                my $qub = &min(@upper);
                my $f1;

                if($qub eq $qlb)
                {
                        $f1 = $qub;
                }
                else
                {
                        $f1 = ($q - $qlb)/($qub - $qlb);
                }

                my $vqt = (1 - $f1)*$pkey{$qlb} + $f1*$pkey{$qub};	
                push @quant, $vqt;
        }
	my @tmpstd = &ave_var(@vector);				
	my $t_std = sqrt($tmpstd[1]/$tmpstd[2])/3;
	push @quant, $t_std;
	
        return @quant;
}


sub max 
{
  	my @array=@_;
  	my $max;
  	$max=$array[0];

  	foreach my $item (@array)
	{
     		if($max < $item)
		{
	  		$max=$item;
     		}
  	} 
  	return $max;
}	

sub min
{
        my @array = @_;
        my $min = $array[0];
        foreach my $item (@array)
        {
                if($min >= $item)
                {
                        $min = $item;
                }
        }
        return $min;
}


sub ave_var
{
        my @array = @_;
        my $sum=0;
        my $mean=0;
        for (my $m=0; $m <= $#array; $m++)
        {
                $sum += $array[$m];
        }
        my $number = $#array+1;
        $mean=$sum/$number;

        my $var = 0;
        for (my $j=0;$j<=$#array;$j++)
        {
                $var=$var+($array[$j]-$mean)**2;
        }

        return $mean,$var, $number;
}

sub likelihood {
   use constant PI    => 4 * atan2(1, 1);
   my $sub_ind =shift @_;
   my $sub_u =shift @_;
   my $sub_stdev =shift @_;
   my $ln = -log($sub_stdev) - 0.5*log(2*PI) - 0.5*(($sub_ind-$sub_u)**2)/($sub_stdev**2); 
   return $ln;
}

sub stddev {
   my @data=@_;
   my $number;
   my $average;
   my $vector_length=scalar(@data);
   my $sum=0;
   for (my $sub_i=0;$sub_i<$vector_length;$sub_i++){
      $sum += $data[$sub_i];
   }
   my $mean=$sum/($vector_length);
   my $variancesum=0;
   for (my $sub_i=0;$sub_i<$vector_length;$sub_i++){
      $variancesum=$variancesum +($data[$sub_i]-$mean)**2;
   }
   my $stdev = ($variancesum/($vector_length))**.5;
   return ($mean,$stdev);
}

