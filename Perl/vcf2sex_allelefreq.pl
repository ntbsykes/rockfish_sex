#!/bin/perl

#This tests if there is a difference in allele frequency between sexes.
use strict;
use warnings;

my $pop_file = $ARGV[0];
my $vcf = $ARGV[1];

my %pop;
open POP, $pop_file;
while(<POP>){
  chomp;
  my @a = split(/\t/,$_);
  $pop{$a[0]} = $a[1];
}
close POP;

my %sex;
print "Contig\tPosition\tderived_m\tderived_f\tgenotyped_m\tgenotyped_f\tchi_square";

open VCF, $vcf;
while(<VCF>){
  chomp;
  my %geno;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^##/){next;}
  if ($_ =~ m/^#/){
    foreach my $i (9..$#a){
      $sex{$i} = $pop{$a[$i]};
    }
  }else{
    my $Contig = $a[0];
    my $Position = $a[1];
    my @alt = split(/,/,$a[4]);
    if ($alt[1]){next;} #Skip multi-allelic sites
    foreach my $i (9..$#a){
      if (($a[$i] eq '.') or ($a[$i] eq './.')){
      }else{
        my @infos = split(/:/,$a[$i]);
        my @genotypes = split(/\//,$infos[0]);
        foreach my $j (0..1){
          $geno{$sex{$i}}{$genotypes[$j]}++;
        }
      }
    }
    #Now calculate the two proportion test.
    my @sexes = ("M","F");
    foreach my $sex (@sexes){
      foreach my $j (0..1){
        unless($geno{$sex}{$j}){
          $geno{$sex}{$j} = 0;
        }
      }
    }
    my $total_m = ($geno{"M"}{1} + $geno{"M"}{0});
    my $total_f = ($geno{"F"}{1} + $geno{"F"}{0});
    if ($total_f eq 0){
      next;
    }
    if ($total_m eq 0){
      next;
    }
    my $derived_m = $geno{"M"}{1};
    my $derived_f = $geno{"F"}{1};
    my $n = $total_m + $total_f;

    my $chi_bot = (($geno{"M"}{0} + $geno{"F"}{0})*($derived_m + $derived_f)*($geno{"M"}{0} + $derived_m)*($geno{"F"}{0} + $derived_f));
    my $chi_top = ($n*(abs($geno{"M"}{0}*$derived_f - $geno{"F"}{0}*$derived_m) - $n/2)**2);
    my $chi_square;
    if ($chi_bot > 0){
      $chi_square = $chi_top / $chi_bot;
    }else{
      $chi_square = "NA";
    }
    print "\n$Contig\t$Position\t$derived_m\t$derived_f\t$total_m\t$total_f\t$chi_square";

  }
}
