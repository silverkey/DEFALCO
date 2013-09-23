#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my @file = glob('*.bam');
my $href = {};
my $tab = {};

mkdir('BAM');

open(TAB,">targets.txt");
open(SCRIPT,">SCRIPT");

foreach my $file(@file) {
  print "Found BAM file $file\n";
  my @section = split("\_",$file);
  print "Relevant name section: $section[1]\n\n";
  push(@{$href->{$section[1]}},$file);
  my $sample = "$section[1]";
  $sample =~ s/\-\d//;
  $tab->{"$section[1]\t$sample\n"} ++;
}

print Dumper $href;

foreach my $section(keys %$href) {
  my $command = "samtools merge BAM\/$section\.bam";
  foreach my $file(@{$href->{$section}}) {
    $command .= " $file";
  }
  print SCRIPT $command."\n\n";
  print $command."\n\n";
  system($command);
  print "\n\nDONE $section\n\n";
}

print "\n\nDONE ALL!!!\n\n";

foreach my $row(keys %$tab) {
  print TAB $row;
}

print  "THE END!";
