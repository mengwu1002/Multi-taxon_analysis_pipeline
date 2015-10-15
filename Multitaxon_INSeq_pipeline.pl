#!/usr/bin/perl -w

#--------------------------------------------------------------------------#
#  Author: Meng Wu, the Gordon Lab, Washington University in St. Louis     #                                                                   
#                                                                          #
#  File:   Multitaxon_INSeq_pipeline.pl                                    #
#  Description: This is the analysis pipeline script 
#  Date:   2014-09-05                                                      #                  
#  Version: 1.30                                                           #
#                                                                          # 
#  Usage:                                                                  #
#       see Usage in the code                                              #
#                                                                          #
#  Contact: mengwu@wustl.edu                                               #
#--------------------------------------------------------------------------#



use strict;
use warnings;
use Getopt::Long;
use IO::File;
use FindBin qw($Bin);
use IO::Uncompress::AnyUncompress;  # can read RFC 1950, RFC 1951, gzip, zip, bzip2, lzop, lzf

if (!@ARGV){
   print "Usage: perl $0 -i <the raw reads file>  -m <Barcodes mapping file>  -s <indexed genome name> [-index <index seq reads file>  -d <length_disrupt_percent (max=1)> [-operon -c <operon_probability_cutoff (max=1)>] [-arrayed] [\n";
   print "Required argument: ";
   print "-i gives the input raw reads file, -m gives the mapping file with the barcode sequence and name for each sample in the format as <barcode>\t<Sample name>, -s gives the name of the indexed genome that the reads should map to";
   print "Optional argument: \n";
   print "-d gives the region of the gene in which insertions are expected to disrupt gene function. The default is 1, which means when insertion falls anywhere in the gene (100%), the gene function will be disrupted. Setting the -d argument to 0.9, for example, would exclude insertions in the distal 10% of the gene when calculating the total number of reads/insertions for that gene. -operon (no argument) specifies that putative downstream (polar) insertions should be calculated based on a user-provided operon probability file. -c is the cutoff for operon probability, default is 1, which means only when the probability of two genes being in an operon is equal to 1 (100%), the polar effect will be considered for the downstream gene. Setting the -c argument to 0.8, for example, will calculate a polar effect for the downstream gene if the probability of the genes being in an operon is at least 0.8 (80%). -arrayed (no argument) is the option for the arrayed library. Refer to README for more information\n";
   exit;
}

my $datasource;   #### The raw reads file which will be analyzed
my $mismatch1=0;  #### Mismatches allowed in finding the transposon sequence
my $mismatch2=0;  #### Mismatches allowed when mapping the sequences to reference using bowtie
my $mapping;      #### Mapping file name
my $poolname="INSEQ_experiment";    #### The prefix of all the analyzed files
my $index;                          #### The name of the index fold in the indexes
my $operon=0;                       #### The option if analyze the sequence using operon information
my $cutoff=100;                     #### The cutoff for operon probability, default is 100
my $disruption=100;                 #### The disruptable percentage of genes. Only when transposon was inserted into this proximal region of gene, the gene is considered interrupted since the insertion in the distal region of a gene may not affect the function of gene at all. The default is 100.
my $arrayed=0;     #### The option if the data is from an arrayed library
my $path=$Bin;     #### Get the path for this analysis package
my $indexfile;

GetOptions (
  "index=s"=>\$indexfile,
  "i|input=s"=>\$datasource,
  "s|species=s"=>\$index,
  "m|mappingfile=s"=> \$mapping,
  "operon"=>\$operon,
   "d|disruption=s"=>\$disruption,
   "c|cutoff=i"=>\$cutoff,
   "arrayed"=>\$arrayed,
   "mismatch1=i"=>\$mismatch1,
   "mismatch2=i"=>\$mismatch2,
);

print "test $mismatch1\n";
my $bowtie_path="";

open BT,"$path/config.txt" || die "Error, please check configuration file as config.txt in the package";
while (<BT>){
    chomp;
    if ($_=~m/\=\"(.*)\"/){
         $bowtie_path=$1;
         if ($1 eq ""){
            die "please specify where the bowtie directory is\n";
         }
    }
}


my $outputdir=`pwd`;
chomp $outputdir;
#print "$outputdir \n";

my $input=$poolname.".scarf";

`rm $input`;

`ln -s $datasource $input`;

mkdir "bcsortedseqs";     
mkdir "results";

# Step 1
# Use the mapping file to assign each read to a barcode and store the output file as inputfile_assigned.txt", and store the statistics in the log file

my $barcode_assigned=$input."_assigned.txt";
my $logfile=$input.".log";

open CODES,$mapping || die "Error: can't open the mapping file, check the README for more details about mapping file\n";
my @codes;
my %codes_hash; #$codes_hash{$code} = sampleID
my %codes_number; #$codes_number{$code} = sample number (within a given sequencing lane)
my $sample_count=0;
my %has_barcode; #$has_barcode{$code} = number of reads with barcode $code 
my $total=0;
my $total_mapped=0;
my $code_length;
my $unbarcoded="INSEQ_without_BC.txt";
my $libraries;

while (my $line=<CODES>){
  chomp $line;
  my @codes_array = split (/\s+/, $line);
  $code_length=length ($codes_array[0]);   
  $codes_hash{$codes_array[0]}=$codes_array[1];
  $codes_number{$codes_array[0]}=$sample_count;
  $sample_count++;
  @{$libraries->{$codes_array[1]}}=split (/\,/,$codes_array[2]);
}
close CODES;

my $nonbc_count=0;
open NONBC,">$unbarcoded";

my $fh_input=new IO::Uncompress::AnyUncompress $input or die "can't open $input: $!\n";
my $fh_index=new IO::Uncompress::AnyUncompress $indexfile or die "can't open $indexfile: $!\n";

my $seq_length;

open OUT,">$barcode_assigned";
while (my $lineF=<$fh_input>) {
  chomp $lineF;
  $total++;
  my $seq1=<$fh_input>;
  chomp $seq1;
  my $qualN=<$fh_input>;
  my $qualS=<$fh_input>;
  <$fh_index>;
  my $indexseq=<$fh_index>;
  chomp $indexseq;
  <$fh_index>;
  <$fh_index>; 
  my $seq_barcode = substr ($seq1, 0, $code_length,"");
  if (exists $codes_hash{$seq_barcode}){ #if first bases of read are one of the barcodes
    $total_mapped++;
    unless (exists $has_barcode{$seq_barcode}) {
      $has_barcode{$seq_barcode}=1;
    } else {
      $has_barcode{$seq_barcode}++;
    }
    $seq_length=length ($seq1);
    my $seq = $seq1.$indexseq;
    print OUT ">".$codes_hash{$seq_barcode}.":".$seq_barcode."\n".$seq."\n";
  }else{
    $nonbc_count++;
    my $header="seq_".$nonbc_count;
    print NONBC ">$header\n$seq1\n";
  }
}
close OUT;

my $index_position=$seq_length+1;
my $filehandlehash;
foreach my $a (keys %codes_hash) {
        if ($has_barcode{$a}){
           foreach my $library (@{$libraries->{$codes_hash{$a}}}){
              my $fh = IO::File->new(">bcsortedseqs\/$input\_$codes_hash{$a}\_$library\.fas");
              $filehandlehash->{$codes_hash{$a}}->{$library} = $fh;
           }
        }
}

open LOG,">$logfile";
print LOG "Number of total reads in $input is $total\n";
print LOG "Total mapped $total_mapped\t".100*$total_mapped/$total." %\n";
print LOG "Barcode\tSample\tReads\tPercent\n";
foreach my $key (sort {$codes_number{$a} <=> $codes_number{$b}} keys %codes_number){
  if ($has_barcode{$key}){
    my $pct = 100*$has_barcode{$key}/$total;
    print LOG "$key\t$codes_hash{$key}\t$has_barcode{$key}\t$pct\n";
  }
}

#Step 2: Trimmed reads to remove transposon, append 16bp reads with a 5'N 
my $species;
my $position;

$species->{"ATCG"}="Bt7330";
$species->{"TCGA"}="Bt7330";
$species->{"ACGT"}="Bvulgatus";
$species->{"TACG"}="Bvulgatus";
$species->{"CGAT"}="Bovatus";
$species->{"GATC"}="Bovatus";
$species->{"GTAC"}="Buniformis";
$species->{"CGTA"}="Buniformis";
$species->{"TCAG"}="BWH2";
$species->{"ACTG"}="BWH2";
$species->{"TCGT"}="BtVPI";
$position->{"ATCG"}="Forward";
$position->{"TCGA"}="Reverse";
$position->{"ACGT"}="Forward";
$position->{"TACG"}="Reverse";
$position->{"CGAT"}="Forward";
$position->{"GATC"}="Reverse";
$position->{"GTAC"}="Forward";
$position->{"CGTA"}="Reverse";
$position->{"TCAG"}="Forward";
$position->{"ACTG"}="Reverse";
$position->{"TCGT"}="Both";

open IN,$barcode_assigned;

my @tn_array=split (//,"ACAGGTTG");
$total=0;
my $count=0;
my $exact_match_count=0;
my $mismatch1_count=0;
my $mismatch2_count=0;
my $header; 
my %tn_match; #$tn_match{sampleID} = number of reads that have TN sequence at allowable #mismatches
my @header; 
my $bc;
my $number;
my $length;
my $NON_TN="INSEQ_nonTN";
my $count_nonTN=0;
my $innerbarcode;
open NONTN,">$NON_TN";
my $bag;
my $inner_count;

while (my $line =<IN>){ #go through each line
  chomp $line;
  if ($line =~m/^>/){ #if a header row
    $header = $line;
    @header = split (/:/,$line);
    $bc=$header[0];
    $bc=~s/>//;
    unless (exists $tn_match{$header[0]}){
      $tn_match{$header[0]}=0;
    }
  } else { #if a sequence row
    my $seq = $line;
    my $pos1_match=0;
    my $pos2_match=0;
    #does read contain transposon at bp 20 or 21?
    my $tn_pos1=substr($seq, 16, 8);
    my $tn_pos2=substr($seq, 17, 8);
    if ($tn_pos1 eq 'ACAGGTTG') {
      $pos1_match=1;
      $exact_match_count++;
      $count++;
    } else {
      if ($tn_pos2 eq 'ACAGGTTG') {
        $pos2_match=1;
        $exact_match_count++;
        $count++;
      }
    } 
    if ($mismatch1 == 1) { #if 1bp tn mismatches allowed
      my $pos1_score=0;
      my $pos2_score=0;
      my @pos1_array = split (//,$tn_pos1);
      my @pos2_array = split (//,$tn_pos2);
      for (my $n=0; $n<8;$n++){
        if ($pos1_array[$n] eq $tn_array[$n]) {
          $pos1_score++;
        }
        if ($pos2_array[$n] eq $tn_array[$n]) {
          $pos2_score++;
        }
      }
      if ($pos1_score == 7) {
        $pos1_match=1;
        $mismatch1_count++;
        $count++;
      } 
      if ($pos2_score == 7) {
        $pos2_match=1;
        $mismatch1_count++;
        $count++;
      }
    }
    if ($mismatch1 == 2) { #if 2bp tn mismatches allowed
      my $pos1_score=0;
      my $pos2_score=0;
      my @pos1_array = split (//,$tn_pos1);
      my @pos2_array = split (//,$tn_pos2);
      for (my $n=0; $n<8;$n++){
        if ($pos1_array[$n] eq $tn_array[$n]) {
          $pos1_score++;
        }
        if ($pos2_array[$n] eq $tn_array[$n]) {
          $pos2_score++;
        }
      }
      if ($pos1_score == 7) {
        $pos1_match=1;
        $mismatch1_count++;
        $count++;
      } 
      if ($pos2_score == 7) {
        $pos2_match=1;
        $mismatch1_count++;
        $count++;
      }
      if ($pos1_score == 6) {
        $pos1_match=1;
        $mismatch2_count++;
        $count++;
      } 
      if ($pos2_score == 6) {
        $pos2_match=1;
        $mismatch2_count++;
        $count++;
      }
    }
    
     #if match found
    if ($pos1_match==1){
      #trim transposon sequence  
      my $trimmed_seq = substr ($seq, 0, 16);
      $tn_match{$header[0]}++;
      $innerbarcode=substr($seq,$index_position,4);
    #  print "$bc\t$innerbarcode\t$species->{$innerbarcode}\t$position->{$innerbarcode}\n";
      if ($species->{$innerbarcode}){
        push @{$bag->{$bc}->{$species->{$innerbarcode}}},$trimmed_seq;
        $inner_count->{$bc}->{$species->{$innerbarcode}}->{$position->{$innerbarcode}}++;
      }
    }
    elsif ($pos2_match==1){
      #trim transposon sequence
      my $trimmed_seq = substr ($seq, 0, 17);
      $tn_match{$header[0]}++;
      $innerbarcode=substr($seq,$index_position,4);
      if ($species->{$innerbarcode}){
         push @{$bag->{$bc}->{$species->{$innerbarcode}}},$trimmed_seq;
         $inner_count->{$bc}->{$species->{$innerbarcode}}->{$position->{$innerbarcode}}++;
      }
    }else{
      $count_nonTN++;
      $header=$header."_".$count_nonTN;
      print NONTN "$header\n$seq\n";
    }
  }
}
close OUT;

foreach my $sample (sort keys %$bag){   
       foreach my $library (@{$libraries->{$sample}}){
           my $fh=$filehandlehash->{$sample}->{$library};
             foreach my $read (@{$bag->{$sample}->{$library}}){
             $header=$sample.":".$library;
             print $fh ">$header\n";
             print $fh "$read\n";
           }
      }
}

my $non_TN_percent=sprintf ("%.1f",$count_nonTN/$total_mapped*100);
print LOG "In the trimming process, $count_nonTN ($non_TN_percent%) are not trimmed, Here is the statistics for each sample\n";
print LOG "Sample\tTrimmed\tPercentage\n";
foreach my $key (sort {$codes_number{$a} <=> $codes_number{$b}} keys %codes_number){
  if ($has_barcode{$key}){
    my $trimkey=">".$codes_hash{$key};
    my $percentage=sprintf("%.1f",$tn_match{$trimkey}/$has_barcode{$key}*100);
    print LOG "$codes_hash{$key}\t$has_barcode{$key}\t$tn_match{$trimkey}\t$percentage\n";
  }
}

foreach my $sample (sort keys %{$bag}){
        foreach my $code (sort keys %$species){
                  if ($inner_count->{$sample}->{$species->{$code}}->{$position->{$code}}){
                  print LOG "$sample\t$species->{$code}\t$position->{$code}\t$inner_count->{$sample}->{$species->{$code}}->{$position->{$code}}\n";
                  }
        }
}

close LOG;

my $clean_up="clean_up.sh";
open CL,">$clean_up";
print CL "rm $barcode_assigned\n",
         "rm mappingjobs_*\n",
         "rm -r bcsortedseqs\n";

#Create job files for each sample 
chdir "bcsortedseqs";
my @sorted_seq=`ls *.fas`;
chdir $outputdir;

foreach my $sorted_seq (@sorted_seq) {
          chomp $sorted_seq;
          if ($sorted_seq=~m/$input\_(\w+)\_(\w+)\.fas/){
             my $sample=$1;
             my $index=$2;
             open (TASKSFORBC, ">$outputdir\/mappingjobs\_$sorted_seq\.job") || die "Error: Can't create $outputdir\/mappingjobs\_$sorted_seq\.job\n\n";
             my @ptt=`ls $path/indexes/$index/*.ptt`;
             print TASKSFORBC  "$bowtie_path/bowtie -m 1 --best --strata -a --fullref -n $mismatch2 -l 17  $path\/indexes\/$index\/$index -f bcsortedseqs\/$sorted_seq results\/$sorted_seq\.bowtiemap\n";
             print TASKSFORBC  "perl $path/scripts/process_bowtie_output.pl results\/$sorted_seq\.bowtiemap \n";
             foreach my $ptt (@ptt){
               chomp $ptt;
               if ($ptt =~m/$path\/indexes\/$index\/(.*)\.ptt/){
               print TASKSFORBC  "perl $path/scripts/normalize_processed_filter.pl results\/$sorted_seq\.bowtiemap_processed.txt_$1\n";
               if (!$arrayed){
                  if ($operon==0){
                      print TASKSFORBC  "perl $path/scripts/map_genes.pl $ptt results\/$sorted_seq\.bowtiemap_processed.txt_$1_filter_cpm.txt $disruption\n";
                  }elsif ($operon==1){
                      print TASKSFORBC "perl $path/scripts/map_genes_operon.pl $ptt $path/indexes\/$index\/$1.operons  $disruption $cutoff results\/$sorted_seq\.bowtiemap_processed.txt_$1_filter_cpm.txt \n";
                  }
               }else{
               }  
              }
         }
            close TASKSFORBC;
       }
         system("qsub -l h_vmem=4G mappingjobs\_$sorted_seq\.job");
}





