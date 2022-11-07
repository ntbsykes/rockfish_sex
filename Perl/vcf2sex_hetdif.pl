#!/bin/perl

#This tests if there is a difference in heterozygosity between sexes.
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
print "Contig\tPosition\thet_m\thet_f\tgenotyped_m\tgenotyped_f\tchi_square";

open VCF, $vcf;
while(<VCF>){
  chomp;
  my %het;
  my %present;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^##/){next;}
  if ($_ =~ m/^#/){
    foreach my $i (9..$#a){
      $sex{$i} = $pop{$a[$i]};
    }
  }else{
    my $Contig = $a[0];
    my $Position = $a[1];
    foreach my $i (9..$#a){
      if (($a[$i] eq '.') or ($a[$i] eq './.')){
      }else{
        $present{$sex{$i}}++;
        my @infos = split(/:/,$a[$i]);
        my @genotypes = split(/\//,$infos[0]);
        if ($genotypes[0] ne $genotypes[1]){
          $het{$sex{$i}}++;
        }
      }
    }
    #Now calculate the chi-square test.
    unless($present{"M"}){
      $present{"M"} = 0;
      next;
    }
    unless($present{"F"}){
      $present{"F"} = 0;
      next;
    }
    unless($het{"M"}){
      $het{"M"} = 0;
    }
    unless($het{"F"}){
      $het{"F"} = 0;
    }
    my $hom_m = $present{"M"} - $het{"M"};
    my $hom_f = $present{"F"} - $het{"F"};
    my $het_m = $het{"M"};
    my $het_f = $het{"F"};
    my $present_m = $present{"M"};
    my $present_f = $present{"F"};
    my $n = $present{"M"} + $present{"F"};
    my $chi_bot = (($het{"M"} + $het{"F"})*($hom_m + $hom_f)*($het{"M"} + $hom_m)*($het{"F"} + $hom_f));
    my $chi_top = ($n*(abs($het{"M"}*$hom_f - $het{"F"}*$hom_m) - $n/2)**2);
    my $chi_square;
    if ($chi_bot > 0){
      $chi_square = $chi_top / $chi_bot;
    }else{
      $chi_square = "NA";
    }
    print "\n$Contig\t$Position\t$het_m\t$het_f\t$present_m\t$present_f\t$chi_square";
  }
}
