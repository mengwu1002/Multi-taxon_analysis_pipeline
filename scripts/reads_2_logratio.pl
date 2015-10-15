#!/usr/bin/perl

#--------------------------------------------------------------------------#
#  Author: Meng Wu, the Gordon Lab, Washington University in St. Louis     #                                                                   
#                                                                          #
#  Usage:                                                                  #
#       see Usage in the code                                              #
#                                                                          #
#  Contact: mengwu@wustl.edu                                               #
#--------------------------------------------------------------------------#

if (!@ARGV){
   print "Usage: perl $0 <Input reads> <Summary reads> <Meta file>\n";
   exit;
}

use FindBin qw($Bin);
my $path=$Bin;

my $ref=shift @ARGV;
my $rawreads=shift @ARGV;
my $metafile=shift @ARGV;

my $Rcode="R.script";
open (R,">$Rcode")|| die "cannot open $Rcode\n";

print R"
   source(\"$Bin/Reads_logtransform_correlation_pseudo.R\")
   Log_corr(\"$ref\",\"$rawreads\",\"$metafile\")
   q()
";

system ("R <$Rcode --vanilla");
