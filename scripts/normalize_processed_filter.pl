#!usr/bin/perl -w -s

#11/2/09
#per-sample count per million normalization on processed mapping output data
#06/12/2013
#
#input file format:
#SampleID \t chromosome \t Coordinate \t L_read \t R_read \t Total_reads \n
#04/09/2015
#Adjusted the cutoff to 3


if (!@ARGV){
  print "\nUsage: perl $0 <input.txt>\n";
  exit;
}


open IN, "$ARGV[0]";
open OUT, ">$ARGV[0]_filter_cpm.txt";

my $logfile;
if ($ARGV[0]=~/(\w+\.scarf)/){
   $logfile=$1.".log";
}

open LOG,">>$logfile";
print LOG "SampleID\tRawReadsAfterFiltered\tCoverage\tScale factor\n";
#print OUT "SampleID\tLocation\tCPM-normalized count\n";

#read file into 2D hash
my $data; # $data->{sampleID}->{location} = count
my $sample;

while (my $line = <IN>){
  chomp $line;
  my @temp = split /\s+/, $line;
  $data->{$temp[0]}->{$temp[1]} = $temp[4];
}


  my @sites;
  my @counts;
  my @norm_counts;
  foreach my $chromosome (keys %$data){ 
    foreach my $location (keys %{$data->{$chromosome}}){  #for each location within that sampleID
      if ($data->{$chromosome}->{$location}>=3){
        push (@sites, $location);
        push (@counts, $data->{$chromosome}->{$location});
      }
    }
  }
  my $sum=0;
  my $coverage=scalar @sites;
  for (my $i=0; $i<scalar @sites; $i++){
    $sum = $sum+$counts[$i];
  }
  my $scale_factor = 1000000/$sum;
  print LOG "$ARGV[0]\t$sum\t$coverage\t$scale_factor\n";
  foreach my $chromosome (keys %$data){ 
    foreach my $location (keys %{$data->{$chromosome}}){  #for each location within that sampleID
              if ($data->{$chromosome}->{$location}>5){
                  print OUT"$chromosome\t$location\t".$data->{$chromosome}->{$location}*$scale_factor."\n";
              }
    }
  }

close LOG;
