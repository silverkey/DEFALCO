#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my @file = glob('*.bam');
my $href = {};
my $tab = {};

mkdir('ANALYSIS');
system("cp index_BAM_in_folder.pl ANALYSIS/");

open(TAB,">ANALYSIS/targets.txt");
open(TAB2,">ANALYSIS/mapping.txt");

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
  my $c = 1;
  foreach my $file(@{$href->{$section}}) {
    print "Copy $file to ANALYSIS/$section\_$c\.bam\n";
    system("cp $file ANALYSIS/$section\_$c.bam");
    print TAB2 "$file\t$section\t$section\_$c\.bam\n";
    $c++;
  }
}

print "\n\nIndexing BAM files\n\n";

chdir("ANALYSIS");

system("perl index_BAM_in_folder.pl > index.bam.folder.log");

foreach my $row(keys %$tab) {
  print TAB $row;
}

print  "THE END!";
