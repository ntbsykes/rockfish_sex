#!/bin/perl
use warnings;
use strict;

#Open position file & use faidx
open(POSITIONS,"meta/outgroup_amh_pos.tsv");
while(<POSITIONS>){
  chomp;
  my ($species,$chr,$start,$end,$direction) = split(/\s/);
  if ($direction eq "f") {
    open(SAMTOOLS,"bin/samtools-1.13/samtools faidx /home/ntbsykes/projects/def-gowens/gowens/sebastes/aleutianus_ragtag_genomes/$species/*/ragtag_output/ragtag.scaffolds.fasta $chr:$start-$end |");
  }else{
    open(SAMTOOLS,"bin/samtools-1.13/samtools faidx -i /home/ntbsykes/projects/def-gowens/gowens/sebastes/aleutianus_ragtag_genomes/$species/*/ragtag_output/ragtag.scaffolds.fasta $chr:$start-$end |");
  }
  while(my $line = <SAMTOOLS>){
    if ($line =~ m/^>/){
      print ">$species.$chr\n";
    } else {
      print $line;
    }
  }
  close(SAMTOOLS);
}
close(POSITIONS);
