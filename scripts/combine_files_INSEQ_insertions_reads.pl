#!/usr/bin/perl -w
#--------------------------------------------------------------------------#
#  Author: Meng Wu, the Gordon Lab, Washington University in St. Louis     #                                                                   
#                                                                          # 
#  Usage:                                                                  #
#       see Usage in the code                                              #
#                                                                          #
#  Contact: mengwu@wustl.edu                                               #
#--------------------------------------------------------------------------#

if (!@ARGV){
   print "Usage: perl $0 <the output files from the analysis pipeline, the mapped insertions binned to each gene>\n";
   exit;
}

my @experiment=@ARGV;
my $sample;
my @sample;
my $ratio;
my $insertion;

foreach my $file(@experiment){
        chomp $file;
        print "Get the files\n";
        if ($file=~m/(.*)\_filter\_cpm\.txt/){
          $sample=$1;
          push @sample,$sample;
          open FILE,$file;
          <FILE>;
          while (<FILE>){
             chomp;
             @temp=split /\t/,$_;
             $ratio->{$temp[0]}->{$sample}=$temp[2];
             $insertion->{$temp[0]}->{$sample}=$temp[1];
         }
    }
}

my $out1="summary_reads.txt";
my $out2="summary_insertions.txt";
open OUT1,">$out1";
open OUT2,">$out2";
sort @sample;
my @header;
push @header,"Gene";
foreach my $sample (@sample){
  #  if (defined $candidate->{$sample}){
           push @header,$sample;
  #  }
}
my $line=join "\t",@header;
print OUT1 "$line\n";
print OUT2 "$line\n";
shift @header;
foreach my $gene(sort keys %$ratio){
    my @out1;
    my @out2;
    push @out1,$gene;
    push @out2,$gene;
    foreach my $sample (@sample){
          push @out1,$ratio->{$gene}->{$sample};
          push @out2,$insertion->{$gene}->{$sample};
    }
    my $line1=join "\t",@out1;
    my $line2=join "\t",@out2;
    print OUT1 "$line1\n";
    print OUT2 "$line2\n";
}

close OUT1;
close OUT2;

