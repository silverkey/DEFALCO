#!/usr/bin/perl
use strict;
use warnings;

my @bam = glob('*.bam');

foreach my $bam(@bam) {
  my $command = "samtools index $bam";
  print "$command\n";
  system($command);
}

