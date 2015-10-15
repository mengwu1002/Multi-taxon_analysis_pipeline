#!usr/bin/perl -w -s


#--------------------------------------------------------------------------#
#  Author: Meng Wu, the Gordon Lab, Washington University in St. Louis     #
#                                                                          #
#  File:   map_genes.pl                                                    #
#  Date:   2014-09-05                                                      #
#  Version: 1.30                                                           #
#                                                                          #
#  Usage:                                                                  #
#       see Usage in the code                                              #
#                                                                          #
#  Contact: mengwu@wustl.edu                                               #
#--------------------------------------------------------------------------#

use POSIX qw(log10);

if (@ARGV !=5){
  print "\nUsage: perl map_genes.pl <ptt> <operons_file> <length_disrupt_percent (max=1)> <operon_probability_cutoff (max=1)> <infile.txt>\n";
  exit;
}

open PTT, $ARGV[0];
open OPS, $ARGV[1];

my $percent = $ARGV[2]; #percent of gene considered "hit"
my $operon_percent=$ARGV[3]; #operon probability score cutoff

##########
# Step 1: process_ptt
##########
# this module delimits genes based on ptt file and % of gene considered "hit"
# input is (.ptt file from NCBI, gene_disrupt_percent)
# output is a series of arrays:
#print "processing ptt file\n";
my @left_end;
my @right_end;
my @gene_number;
my @annotation;
my @gene_length;
my @strand;
my @left_end_percent;
my @right_end_percent;
my @gene_length_percent;
my @unique_left;
my @unique_right;
my @shared_left;
my @shared_right;
my @unique_length;
my @shared_length;

my $distal=1-$percent;
#skip past header info in .ptt file
my $ptt_line=<PTT>;
$ptt_line=<PTT>;
$ptt_line=<PTT>;
#go through the ptt file and read info into arrays
my $line_number=0;
while ($ptt_line=<PTT>){  
  chomp $ptt_line;
  my @temp_array = split (/\t/, $ptt_line);  #split ptt line by tabs
  my $gene = $temp_array[5];
#  $gene =~ s/_//g;
  $gene_number[$line_number]=$gene; 
  $strand[$line_number]=$temp_array[1];
  $annotation[$line_number]=$temp_array[8];
  my @coordinates_array = split (/\../, $temp_array[0]);
  $left_end[$line_number]=$coordinates_array[0];
  $right_end[$line_number]=$coordinates_array[1];
  $gene_length[$line_number]=abs($right_end[$line_number]-$left_end[$line_number]);
  $line_number++;
}
my $number_of_genes=scalar @gene_number;
#trim genes to exclude distal x%
for (my $i=0;$i<$number_of_genes; $i++) {
  $gene_length_percent[$i]=int (($gene_length[$i]*$percent)+.5);
  if ($strand[$i] eq "+") { #if on positive strand
    $left_end_percent[$i]=$left_end[$i];
    $right_end_percent[$i]=$left_end[$i]+$gene_length_percent[$i]; 
  }
  if ($strand[$i] eq "-") {
    $right_end_percent[$i]=$right_end[$i];
    $left_end_percent[$i]=$right_end[$i]-$gene_length_percent[$i];
  }
}
# delimit unique and shared regions of first X% of genes
for (my $j=0;$j<$number_of_genes;$j++) {#go through all genes
  my $set_start=0;
  my $set_end=0;
  if ($j==0){ # if first gene, set no overlap from last gene
    $unique_left[$j]=$left_end_percent[$j];
    $set_start=1;
  }
  if ($j==($number_of_genes-1)){ #if last gene
    $unique_right[$j]=$right_end_percent[$j]; # set no overlap from first gene
    $set_end=1;
    $shared_left[$j]="none";
    $shared_right[$j]="none";
    $shared_length[$j]="none";
  }
  #assign unique left end of gene
  if ($set_start==0){ # if not first or last gene
    if ($left_end_percent[$j]>$right_end_percent[$j-1]) { #if left_edge is after previous gene right_edge
      $unique_left[$j]=$left_end_percent[$j];      #unique left_edge is left_end_percent
    }
    if ($left_end_percent[$j]<=$right_end_percent[$j-1]) { #if left_edge is before previous gene right_edge
      $unique_left[$j]=$right_end_percent[$j-1]+1;  #unique left_edge is (previous gene right_end_percent)+1 
    }
  }
  #assign unique right end of gene
  if ($set_end==0){
    if ($right_end_percent[$j]<$left_end_percent[$j+1]) { #if right_edge is before next gene left_edge
      $unique_right[$j]=$right_end_percent[$j]; #unique right_edge is right_end_percent
      $shared_left[$j]="none";   #no shared region for gene j
      $shared_right[$j]="none";
      $shared_length[$j]="none";
    }
    if ($right_end_percent[$j]>=$left_end_percent[$j+1]) { #if right_edge is after next gene left_edge
      $unique_right[$j]=$left_end_percent[$j+1]-1; #unique right_edge is (next gene left_end_percent)-1
      $shared_left[$j]=$left_end_percent[$j+1]; 
      $shared_right[$j]=$right_end_percent[$j];
      $shared_length[$j]=$shared_right[$j]-$shared_left[$j];
    }
  } 
}

##########
# Step 2: establish operon boundaries
##########
#print "establishing operon boundaries\n";
#OPS file should be tab-delim two-column
#column 1 : gene number
#column 2 : probability (max 1) that this gene is in an operon w/ next gene
#example: Gene1    .5
#50% probability that Gene1 is in an operon w/ Gene2

my @operon_prob;
my @operon_name;
my @operon_left;
my @operon_right;
my @operon_strand;
my @operon_left_gene;
my @operon_right_gene;
my @operon_left_gene_number;
my @operon_right_gene_number;
my $name;
my $in_operon=0;
my $count=0;
while (my $line = <OPS>){
  chomp $line;
  my @temp_array = split (/\t/, $line);
  if ($temp_array[0] ne $gene_number[$count]){
    print "ptt and operon files not lining up: $gene_number[$count]\t$temp_array[0]\n";
    exit;
  }
  $operon_prob[$count]=$temp_array[1];
  $count++;
} 

for (my $i=0; $i<$number_of_genes-1; $i++){ #go through each gene
  my $assigned=0;
  if (($in_operon==0) && ($operon_prob[$i]>=$operon_percent)) { #if starting an operon
    $in_operon=1;
    push (@operon_left, $left_end[$i]);
    push (@operon_left_gene_number, $i);
    push (@operon_strand, $strand[$i]);
    $name = $gene_number[$i];
    $assigned=1;
  }
  if (($assigned==0) && ($in_operon==1) && ($operon_prob[$i]<$operon_percent)){ #if ending an operon
    $in_operon=0;
    push (@operon_right, $right_end[$i]);
    push (@operon_right_gene_number, $i);
    $name = $name."-".$gene_number[$i];
    push (@operon_name, $name);
    $assigned=1;
  }
  if (($assigned==0) && ($in_operon==0) && ($operon_prob[$i]<$operon_percent)){ #if one-gene operon
    push (@operon_left, $left_end[$i]);
    push (@operon_right, $right_end[$i]);
    push (@operon_left_gene_number, $i);
    push (@operon_right_gene_number, $i);
    push (@operon_name, $gene_number[$i]."-".$gene_number[$i]);
    push (@operon_strand, $strand[$i]);
    $assigned=1;
  }
}
my $number_of_operons=scalar @operon_name;

########
# Step 3: Map insertion counts to genes and operons
########
my %sampleIDs; #$sampleIDs{sampleID}=1 if present in file
my %genes_hit_direct; #$genes_hit_direct{coordinate} = array of genes hit directly by coordinate
my %genes_hit_polar; #$genes_hit_polar{coordinate} = array of genes in downstream operon from coordinate
my %gene_hitcount_direct; #sample-specific value for how many reads mapping directly to that gene
my %gene_hitcount_polar; #sample-specific value for how many reads mapping directly or upstream of that gene
my %gene_unique_sites; #sample-specific value for number of unique insertion locations in that gene
my %seen_coordinate; #$seen_coordinate{coordinate} = 0 if new, 1 if previously mapped
my %hitcount;

open IN, "$ARGV[4]";

while (my $line=<IN>) {
  chomp $line;
  if ($line =~m/^>/){
    my @temp = split (/\t/, $line);
    my $sampleID = $temp[0];
    my $coordinate = $temp[1]; 
    my $norm_count = $temp[2];
    #print "sample ID $sampleID coordinate $coordinate count $norm_count\n";
    #if sample is new, add it to list of samples
    unless (exists $sampleIDs{$sampleID}){
      $sampleIDs{$sampleID}=1;
    }
    #if coordinate is new, map it to genes and operons
    unless (exists ($seen_coordinate{$coordinate})){
      $seen_coordinate{$coordinate}=1;
      ########
      # identify genes hit directly by coordinate
      ########
      my $found=0;
      my $found_operon=0;
      for (my $check=0;$check<$number_of_genes; $check++){ #go through each pair in the left-right arrays
	if (($left_end_percent[$check]<=$coordinate) && ($coordinate<=$right_end_percent[$check])){ #if hit is in proximal portion of gene
	  $found=1;
	  $found_operon=1;
	  push (@{$genes_hit_direct{$coordinate}}, $gene_number[$check]);  #allows for single coordinate to hit multiple genes
	  #print "$coordinate mapped to gene $gene_number[$check]\n";
	  push (@{$genes_hit_polar{$coordinate}}, $gene_number[$check]);
	} #end gene-hit section 
      } #end go through genes section 
      if ($found==0){
	@{$genes_hit_direct{$coordinate}}[0]="none";
      } #end not-found in gene section
      ##########
      # identify genes hit by polar effect
      ##########
      for (my $check_operon=0;$check_operon<$number_of_operons;$check_operon++) { #go through each operon
	if (($operon_left[$check_operon]<=$coordinate) && ($coordinate <=$operon_right[$check_operon])) {
	  if ($operon_strand[$check_operon] eq "+") {
	    for (my $b=$operon_left_gene_number[$check_operon]; $b<=$operon_right_gene_number[$check_operon]; $b++){ #go through each gene in the operon
	    if ($left_end[$b]>$coordinate){ #if gene starts downstream of coordinate
	      #add this gene to the list of genes hit polarly by this coordinate
	      $found_operon=1;
	      push (@{$genes_hit_polar{$coordinate}}, $gene_number[$b]);
	      #print "$coordinate mapped (polar) to gene $gene_number[$b]\n";
	    } #end gene-hit-polar section
	  } #end go through genes in operon section 
	  } # end positive-strand operon section
	  if ($operon_strand[$check_operon] eq "-"){
	    for (my $j=$operon_right_gene_number[$check_operon]; $j>=$operon_left_gene_number[$check_operon];$j--){
	      if ($right_end[$j]<$coordinate){ #if gene starts downstream of coordinate
		$found_operon=1;
		push (@{$genes_hit_polar{$coordinate}}, $gene_number[$j]);
		#print "$coordinate mapped (polar) to gene $gene_number[$j]\n";
	      } #end gene-hit-polar section
	    } #end go through genes in operon section
	  } # end minus-strand operon section
	} # end coordinate in operon section
      } # end of operon section
      if ($found_operon==0){
	@{$genes_hit_polar{$coordinate}}[0]="none";
      } # end not-found-operon section
    } #end "new coordinate" loop 
   

    ########
    # assign input-specific data to coordinate 
    ########
    $hitcount{$coordinate}->{$sampleID}=$norm_count;
    unless ($genes_hit_direct{$coordinate}[0] eq "none") { 
      for (my $k=0; $k< scalar @{$genes_hit_direct{$coordinate}}; $k++) { #go through each gene in the array referenced by $genes_hit_direct{$coordinate} and add hitcount for this coordinate/input
	$gene_hitcount_direct{$genes_hit_direct{$coordinate}[$k]}->{$sampleID}=0 if !defined $gene_hitcount_direct{$genes_hit_direct{$coordinate}[$k]}->{$sampleID};
        $gene_hitcount_polar{$genes_hit_polar{$coordinate}[$k]}->{$sampleID}=0 if !defined $gene_hitcount_polar{$genes_hit_polar{$coordinate}[$k]}->{$sampleID};
	$gene_unique_sites{$genes_hit_direct{$coordinate}[$k]}->{$sampleID}=0 if !defined $gene_unique_sites{$genes_hit_direct{$coordinate}[$k]}->{$sampleID};
	$gene_hitcount_direct{$genes_hit_direct{$coordinate}[$k]}->{$sampleID}=$gene_hitcount_direct{$genes_hit_direct{$coordinate}[$k]}->{$sampleID}+$hitcount{$coordinate}->{$sampleID};
	$gene_hitcount_polar{$genes_hit_polar{$coordinate}[$k]}->{$sampleID}=$gene_hitcount_polar{$genes_hit_polar{$coordinate}[$k]}->{$sampleID}+$hitcount{$coordinate}->{$sampleID};
	$gene_unique_sites{$genes_hit_direct{$coordinate}[$k]}->{$sampleID}++;
      } #end add direct hits
    } #end coordinate hits direct genes
    unless ($genes_hit_polar{$coordinate}[0] eq "none") {
      for (my $l=0; $l<scalar @{$genes_hit_polar{$coordinate}}; $l++){ #go through each gene in the array referenced by $genes_hit_polar{$coordinate} and add hitcount for this coordinate/input
	$gene_hitcount_polar{$genes_hit_polar{$coordinate}[$l]}->{$sampleID}=0 if !defined $gene_hitcount_polar{$genes_hit_polar{$coordinate}[$l]}->{$sampleID};
	$gene_hitcount_polar{$genes_hit_polar{$coordinate}[$l]}->{$sampleID}=$gene_hitcount_polar{$genes_hit_polar{$coordinate}[$l]}->{$sampleID}+$hitcount{$coordinate}->{$sampleID};
      } #end add polar hits
    } #end coordinate hits polar genes 
  } #end "if line starts with ">" loop
} #end of "go through each line of input" loop

########
# Produce output files
########
my @columns;
foreach my $key (sort (keys %sampleIDs)){
  push (@columns, $key);
}

my $output1=$ARGV[4]."_genes_direct.txt";
my $output2=$ARGV[4]."_genes_polar.txt";
my $output3=$ARGV[4]."_genes_sites.txt";
my $output4=$ARGV[4]."_insertions.txt";


open OUT, ">$output1";
open OUT2, ">$output2";
open OUT3, ">$output3";
open OUT4, ">$output4";

# write output files
# print header info
print OUT "GeneID\t";
print OUT2 "GeneID\t";
print OUT3 "GeneID\t";
print OUT4 "Coordinate\tGenes\t";

for (my $i=0; $i<scalar @columns; $i++){
  print OUT "$columns[$i]\t";
  print OUT2 "$columns[$i]\t";
  print OUT3 "$columns[$i]\t";
  print OUT4 "$columns[$i]\t";
}
print OUT "\n";
print OUT2 "\n";
print OUT3 "\n";
print OUT4 "\n";


for (my $j=0; $j<$number_of_genes; $j++){ #go through each gene
  print OUT "$gene_number[$j]\t";
  print OUT2 "$gene_number[$j]\t";
  print OUT3 "$gene_number[$j]\t";
  for (my $k=0; $k<scalar @columns; $k++){
    $gene_hitcount_direct{$gene_number[$j]}->{$columns[$k]}=0 if !defined $gene_hitcount_direct{$gene_number[$j]}->{$columns[$k]};
    $gene_hitcount_polar{$gene_number[$j]}->{$columns[$k]}=0 if !defined $gene_hitcount_polar{$gene_number[$j]}->{$columns[$k]};
    $gene_unique_sites{$gene_number[$j]}->{$columns[$k]}=0 if !defined $gene_unique_sites{$gene_number[$j]}->{$columns[$k]};
    print OUT "$gene_hitcount_direct{$gene_number[$j]}->{$columns[$k]}\t";
    print OUT2 "$gene_hitcount_polar{$gene_number[$j]}->{$columns[$k]}\t";
    print OUT3 "$gene_unique_sites{$gene_number[$j]}->{$columns[$k]}\t";
  }
  print OUT "$annotation[$j]\n";
  print OUT2 "$annotation[$j]\n";
  print OUT3 "$annotation[$j]\n";
}


#print out coordinate data to insertions.txt file
foreach my $key (keys %hitcount){ #go through each key (coordinate) in %hitcount
  print OUT4 "$key\t";
  for (my $k=0; $k< scalar @{$genes_hit_direct{$key}}; $k++) { #go through each gene in the array referenced by $genes_hit_direct{$key}
    print OUT4 "$genes_hit_direct{$key}[$k]";
    unless ($k==((scalar @{$genes_hit_direct{$key}})-1)){
      print OUT4 ":";
    }
  }
  print OUT4 "\t";
  for (my $j=0; $j<scalar @columns; $j++){
    $hitcount{$key}->{$columns[$j]} = 0 if !defined $hitcount{$key}->{$columns[$j]};
    print OUT4 "$hitcount{$key}->{$columns[$j]}\t";
  }
  print OUT4 "\n";
}


