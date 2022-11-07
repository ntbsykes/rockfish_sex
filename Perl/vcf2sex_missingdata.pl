#!/bin/perl

#This tests if there is a difference in missing data between sexes.
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
print "Contig\tPosition\tpresent_m\tpresent_f\tbias\tchi_square";

open VCF, $vcf;
while(<VCF>){
  chomp;
  my %total;
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
        $total{$sex{$i}}++;
      }else{
        $total{$sex{$i}}++;
        $present{$sex{$i}}++;
      }
    }
    #Now calculate the chi-square test.
    unless($present{"M"}){
      $present{"M"} = 0;
    }
    unless($present{"F"}){
      $present{"F"} = 0;
    }
    unless($total{"M"}){
      $total{"M"} = 0;
    }
    unless($total{"F"}){
      $total{"F"} = 0;
    }
    my $absent_m = $total{"M"} - $present{"M"};
    my $absent_f = $total{"F"} - $present{"F"};
    my $p_m = $present{"M"}/$total{"M"};
    my $p_f = $present{"F"}/$total{"F"};
    my $bias = $p_m - $p_f;
    my $n = $total{"M"} + $total{"F"};
    my $chi_bot = (($present{"M"} + $present{"F"})*($absent_m + $absent_f)*($present{"M"} + $absent_m)*($present{"F"} + $absent_f));
    my $chi_top = ($n*(abs($present{"M"}*$absent_f - $present{"F"}*$absent_m) - $n/2)**2);
    my $chi_square;
    if ($chi_bot > 0){
      $chi_square = $chi_top / $chi_bot;
    }else{
      $chi_square = "NA";
    }
    print "\n$Contig\t$Position\t$p_m\t$p_f\t$bias\t$chi_square";
  }
}
