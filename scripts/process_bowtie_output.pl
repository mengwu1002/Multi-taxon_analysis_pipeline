#!/usr/bin/perl
use warnings;
use strict;

#--------------------------------------------------------------------------#
#  Author: Meng Wu, the Gordon Lab, Washington University in St. Louis     #
#                                                                          #
#  File:   process_bowtie_output.pl                                        #
#  Date:   2014-09-05                                                      #
#  Version: 1.30                                                           #
#                                                                          #
#  Usage:                                                                  #
#       see Usage in the code                                              #
#                                                                          #
#  Contact: mengwu@wustl.edu                                               #
#--------------------------------------------------------------------------#


#Take the bowtie output with the parameter as "-m 1 -l 17 -n 1 -a --best --strata"

if (!@ARGV) {
  print "\nUsage: perl $0 <bowtie_output_file> \n";
  exit;
}

open IN, "$ARGV[0]";



my $total_count=0;
my $count; 

my $insertion_F;
my $insertion_R;
my %sites;     #$sites{sampleID} = number of insertion locations for that sample
my $chromosome;

while (my $line=<IN>) {
  if ($line !~/^#/){
    my @data = split (/\t/, $line);
    $total_count++;
    $chromosome=$data[2];
    my $length=length $data[4];
    my $site;
    $count->{$chromosome}->{$data[2]}++;
    if ($data[1] eq "+"){ #set insertion location of "T" in "TA" top-strand dinucleotide
	$site = $data[3]+$length-1;
      #add read to count of forward reads for that site in that sample
	$insertion_F->{$chromosome}->{$site}++;
    }
    if ($data[1] eq "-"){
	$site=$data[3]+1;
        $insertion_R->{$chromosome}->{$site}++;
   }
  }
}

my $coverage;
my $sites;
my $sample;



foreach my $key (keys %$insertion_F){ #for each sample in %insertion_F
  if ($key=~m/(.*)/){
   open OUT, ">$ARGV[0]_processed.txt_".$1;
   foreach my $location (keys %{$insertion_F->{$key}}) {
    push @{$sites{$key}},$location;
    print OUT ">$key\t$location\t$insertion_F->{$key}->{$location}\t";
    my $sum = $insertion_F->{$key}->{$location};
    if (exists ($insertion_R->{$key}->{$location})){
      print OUT "$insertion_R->{$key}->{$location}\t";
      $sum = $sum + $insertion_R->{$key}->{$location};
      delete $insertion_R->{$key}->{$location};
    } else {
      print OUT "0\t";
    }
    print OUT "$sum\n";  
   }
  }
} 
foreach my $key (keys %$insertion_R){
 if ($key=~m/(.*)/){
   open OUT, ">>$ARGV[0]_processed.txt_".$1;
   foreach my $location (keys %{$insertion_R->{$key}}) {
   push @{$sites{$key}},$location;
   print OUT ">$key\t$location\t0\t$insertion_R->{$key}->{$location}\t$insertion_R->{$key}->{$location}\n";
   }
 }
}


