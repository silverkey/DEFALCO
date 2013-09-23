#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

#--------------------------------------------------------------------------
# THIS SCRIPT DO THE FOLLOWING:
# 1) ENTER INTO THE FOLDER WITH THE BAM TO ANALYZE
# 2) CHANGE THE HEADER OF THE BAM TO MAKE THEM COMPATIBLE WITH ENSEMBL GTF
# 3) ENTER INTO THE FOLDER WITH THE NEW BAM
# 4) CONVERT EACH BAM INTO SAM 
# 5) RUN HTSEQ WITH THE GIVEN PARAMETER AGAINST THE GIVEN GFF
#--------------------------------------------------------------------------

#------------------------------------------
# PARAMETERS ASSOCIATE TO FILES AND FOLDERS
#------------------------------------------
my $bamdir = '.'; # the folder containing the bam
my $gff = '/home/remo/ANALYSIS/sandro_defalco/BAM/BAM/Homo_sapiens.GRCh37.72.gtf';
my $dir = 'BAM_NEW_HEADERS';
my $headerdir = 'HEADERS';
#------------------------------------------

#------------------------------------------
# PARAMETERS ASSOCIATED TO HTSEQ ANALYSES
#------------------------------------------
my $mode = 'union';
my $stranded = 'no';
my $type = 'exon';
my $idattr = 'gene_id';
#------------------------------------------

#------------------------------------------
# LET'S THE ANALYSES BEGIN...
#------------------------------------------
chdir($bamdir);
mkdir($dir) or die "\nCannot create directory $dir!!!\n\n";
mkdir($headerdir) or die "\nCannot create directory $headerdir!!!\n\n";

my @file = glob('*.bam');

# THE FOLLOWING SECTION IS TO CHANGE THE HEADER TO MAKE IT COMPATIBLE
# THE ENSEMBL CHROMOSOME NAMING THAT IS ALSO SHORTER .... 
foreach my $file(@file) {
  my $header = "HEADERS_$file";
  my $newheader = "NEW_HEADERS_$file";
  my $gethead = "samtools view -H $file > $header";
  exec_command($gethead);
  # CHANGE THE HEADERS
  change_header($header,$newheader);
  my $rehead = "samtools reheader $newheader $file > $dir\/$file";
  exec_command($rehead);
  my $cleanup = "mv HEADERS_* HEADERS/; mv NEW_HEADERS_* HEADERS/";
  exec_command($cleanup);
}

# NOW WE RUN HTSEQ ON EACH OF THE BAM FILES IN PARALLEL
chdir($dir) or die "\nCannot change directory $dir!!!\n\n";
my $pm = new Parallel::ForkManager(20);
$pm->run_on_finish(
  sub {
    my($pid,$exit_code) = @_;
    print "** Just got out of the pool ".
          "with PID $pid and exit code: $exit_code\n";
  }
);
foreach my $file(@file) {
  my $pid = $pm->start and next;
  my $sam = make_sam_name_from_bam($file);
  my $bam2sam = "samtools view $file > $sam";
  exec_command($bam2sam);
  my $htseq = "htseq-count --mode=$mode --stranded=$stranded --type=$type --idattr=$idattr $sam $gff > COUNTS_$file 2>ERR_$file";
  exec_command($htseq);
  $pm->finish;
}

# EXEC SYSTEM CALL AND SAFELY CHECK FOR ERRORS
sub exec_command {
  my $command = shift;
  print "\nLAUNCHING SYSTEM CALL:\n\t$command\n";
  system($command);
  die "ERROR using command:\n\t$command\:\n\t$!" unless $? == 0;
  print "DONE!\n";
}

# HERE WE MAKE THE CHANGES IN THE CHROMOSOME NAMES WORKING ON THE
# HEADER OF THE BAM FILES
sub change_header {
  my $file = shift;
  my $newfile = shift;
  open(IN,$file);
  open(OUT,">$newfile");
  while(my $row = <IN>) {
    if($row =~ /^\@SQ/) {
      # TAKE OUT "chr"
      $row =~ s/SN:chr/SN:/;
      # TAKE OUT ".fa"
      $row =~ s/(SN:\S+)\.fa/$1/;
      # CHANGE "M" in "MT"
      $row =~ s/SN:M/SN:MT/ unless $row =~ /SN:MT/;
    }
    print OUT $row;
  }
}

# CHANGE THE EXTENSION OF THE BAM BEING SURE IT IS THERE
# OTHERWISE WE SIMPLY ADD .SAM EXTENSION
sub make_sam_name_from_bam {
  my $bam = shift;
  $bam =~ s/\.bam$//;
  $bam .= ".sam";
  return $bam;
}

__END__

EXAMPLE OUT BAM HEADERS AS OUTPUT FROM samtools view -H....

@PG	ID:illumina_export2sam.pl	VN:2.3.1	CL:/g/solexa/bin/illumina/CASAVA_v1.8.2/bin/illumina_export2sam.pl --read1=lane213s003331_export.txt
@SQ	SN:chr10.fa	LN:135534747
@SQ	SN:chr11.fa	LN:135006516
@SQ	SN:chr12.fa	LN:133851895
@SQ	SN:chr13.fa	LN:115169878
@SQ	SN:chr14.fa	LN:107349540
@SQ	SN:chr15.fa	LN:102531392
@SQ	SN:chr16.fa	LN:90354753
@SQ	SN:chr17.fa	LN:81195210
@SQ	SN:chr18.fa	LN:78077248
@SQ	SN:chr19.fa	LN:59128983
@SQ	SN:chr1.fa	LN:249250621
@SQ	SN:chr20.fa	LN:63025520
@SQ	SN:chr21.fa	LN:48129895
@SQ	SN:chr22.fa	LN:51304566
@SQ	SN:chr2.fa	LN:243199373
@SQ	SN:chr3.fa	LN:198022430
@SQ	SN:chr4.fa	LN:191154276
@SQ	SN:chr5.fa	LN:180915260
@SQ	SN:chr6.fa	LN:171115067
@SQ	SN:chr7.fa	LN:159138663
@SQ	SN:chr8.fa	LN:146364022
@SQ	SN:chr9.fa	LN:141213431
@SQ	SN:chrM.fa	LN:16571
@SQ	SN:chrX.fa	LN:155270560
@SQ	SN:chrY.fa	LN:59373566

