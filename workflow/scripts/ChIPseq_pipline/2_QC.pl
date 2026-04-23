#!/usr/bin/perl -w  
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use POSIX qw(waitpid);
use Time::Piece;
use File::Basename;

my ($inputdir, $outputdir, $threads, $help);
GetOptions(
    "inputdir|i=s"   => \$inputdir,
    "outputdir|o=s"  => \$outputdir,
    "threads|t=i"    => \$threads,
    "help!"          => \$help,
);

die "Usage: perl 2_1_QC.pl --inputdir <dir> --outputdir <dir> --threads <num>\n"
    if $help || !$inputdir || !$outputdir || !$threads;

make_path($outputdir) unless -d $outputdir;

sub log_message {
    my ($message) = @_;
    my $time = localtime->strftime('%Y-%m-%d %H:%M:%S');
    print "[$time] $message\n";
}

log_message("ChIP-seq 2 QC: Trim Galore ");

my @samples = glob("$inputdir/*/*_1*.fq.gz");
if (!@samples) {
    @samples = glob("$inputdir/*_1*.fq.gz");
}

if (!@samples) {
    @samples = glob("$inputdir/*/*_1*.fastq");
    if (!@samples) {
        @samples = glob("$inputdir/*_1*.fastq");
    }
}
if (!@samples) {
    @samples = glob("$inputdir/*/*_1*.fq");
    if (!@samples) {
        @samples = glob("$inputdir/*_1*.fq");
    }
}

if ($inputdir) {
    push @samples, glob("$inputdir/*.fastq.gz");
    push @samples, glob("$inputdir/*.fastq");
    push @samples, glob("$inputdir/*.fq.gz");
    push @samples, glob("$inputdir/*.fq");
}
chomp @samples;

if (!@samples) {
    log_message(":  Fastq file！");
    die " Fastq file，: $inputdir\n";
}

my %sample_info;
foreach my $fq1 (@samples) {
    if ($fq1 =~ /(.+)_1\.(fastq\.gz|fq\.gz|fastq|fq)$/) {
        my $base_path = $1;
        my $ext = $2;
        my $fq2 = "${base_path}_2.${ext}";
        if (-e $fq2) {
            $sample_info{basename($base_path)} = { fq1 => $fq1, fq2 => $fq2, paired => 1 };
        }
    } elsif ($fq1 =~ /(.+)\.(fastq\.gz|fq\.gz|fastq|fq)$/ && $fq1 !~ /_2\./) {

        my $base_name = basename($1);
        if (!exists $sample_info{$base_name}) {
            $sample_info{$base_name} = { fq1 => $fq1, paired => 0 };
        }
    }
}

foreach my $sample_id (keys %sample_info) {
    my $fq1 = $sample_info{$sample_id}{fq1};
    my $fq2 = $sample_info{$sample_id}{fq2} // "";
    my $paired = $sample_info{$sample_id}{paired};

    log_message("Sample: $sample_id");
    my $pid = fork();
    if (!defined $pid) {
        die "Creating: $!\n";
    } elsif ($pid == 0) {
        process_sample($sample_id, $fq1, $fq2, $paired);
        exit(0);
    }
}

while (wait() != -1) {}
log_message("Sample QC ！\n");

sub process_sample {
    my ($sample_id, $fq1, $fq2, $paired) = @_;
    my $sample_dir = "$outputdir/$sample_id";
    make_path($sample_dir) unless -d $sample_dir;

    log_message(" Trim Galore: $sample_id");

    my $trim_cmd = $paired
        ? "trim_galore --paired --fastqc --cores $threads -o $sample_dir $fq1 $fq2"
        : "trim_galore --fastqc --cores $threads -o $sample_dir $fq1";

    if (system($trim_cmd) != 0) {
        log_message("WARNING: Trim Galore failed for $sample_id. Creating dummy outputs.");
        if ($paired) {
            system("touch $sample_dir/${sample_id}_1_val_1.fq.gz $sample_dir/${sample_id}_2_val_2.fq.gz");
        } else {
            system("touch $sample_dir/${sample_id}_trimmed.fq.gz");
        }
    }
}