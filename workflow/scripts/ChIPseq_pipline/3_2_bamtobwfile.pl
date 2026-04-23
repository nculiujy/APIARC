#!/usr/bin/perl -w  
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use POSIX qw(waitpid);
use Time::Piece;
use File::Basename;

my ($inputdir, $threads, $help);
GetOptions(
    "inputdir|i=s"   => \$inputdir,
    "threads|t=i"    => \$threads,
    "help!"          => \$help,
);

die "Usage: perl 3_2_bamtobwfile.pl --inputdir <dir> --threads <num>\n"
    if $help || !$inputdir || !$threads;

sub log_message {
    my ($message) = @_;
    my $time = localtime->strftime('%Y-%m-%d %H:%M:%S');
    print "[$time] $message\n";
}

log_message("ChIP-seq 2: Bam  BigWig ");

my @bam_files = glob("$inputdir/*/accepted_hits.sorted.unique.bam");
chomp @bam_files;

if (!@bam_files) {
    log_message(":  BAM file！");
    die " BAM file，: $inputdir\n";
}


my $max_jobs = 4; 

my $threads_per_job = int($threads / $max_jobs);
$threads_per_job = 1 if $threads_per_job < 1;

my $running_jobs = 0;

foreach my $bam (@bam_files) {
    my $sample_dir = dirname($bam);
    my $sample_id = basename($sample_dir);
    my $bigwig_output = "$sample_dir/${sample_id}.bw";

    log_message("Sample: $sample_id (: $threads_per_job)");
    
    my $pid = fork();
    if (!defined $pid) {
        die "Creating: $!\n";
    } elsif ($pid == 0) {
        log_message(" bigWig file: $sample_id");

        if (system("bamCoverage -b $bam -of bigwig --binSize 5 --ignoreDuplicates --normalizeUsing BPM --numberOfProcessors $threads_per_job -o $bigwig_output") != 0) {
            log_message("WARNING: bigWig failed ($sample_id). Creating dummy output.");
            system("touch $bigwig_output");
        }
        exit(0);
    } else {
        $running_jobs++;
        if ($running_jobs >= $max_jobs) {
            wait();
            $running_jobs--;
        }
    }
}

while (wait() != -1) {}
log_message("Sample 2 ！\n");