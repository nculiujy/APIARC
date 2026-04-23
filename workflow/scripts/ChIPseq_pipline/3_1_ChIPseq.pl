#!/usr/bin/perl -w  
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use POSIX qw(waitpid);
use Time::Piece;
use File::Basename;

my ($inputdir, $outputdir, $indexdir, $picarddir, $threads, $help, $single_end);
GetOptions(
    "inputdir|i=s"   => \$inputdir,
    "outputdir|o=s"  => \$outputdir,
    "indexdir|x=s"   => \$indexdir,
    "picarddir|p=s"  => \$picarddir,
    "threads|t=i"    => \$threads,
    "single_end|s!"  => \$single_end,
    "help!"          => \$help,
);

die "Usage: perl 3_1_ChIPseq.pl --inputdir <dir> --outputdir <dir> --indexdir <dir> --picarddir <dir> --threads <num> [--single_end]\n"
    if $help || !$inputdir || !$outputdir || !$indexdir || !$picarddir || !$threads;

make_path($outputdir) unless -d $outputdir;

sub log_message {
    my ($message) = @_;
    my $time = localtime->strftime('%Y-%m-%d %H:%M:%S');
    print "[$time] $message\n";
}

log_message("ChIP-seq 1:  ");


my @samples = glob("$inputdir/*/*_1*.fq.gz");
if (!@samples) {
    @samples = glob("$inputdir/*_1*.fq.gz");
}

if (!@samples) {
    @samples = glob("$inputdir/*/*_1*.fq");
    if (!@samples) {
        @samples = glob("$inputdir/*_1*.fq");
    }
}
if (!@samples) {
    @samples = glob("$inputdir/*/*_1*.fastq");
    if (!@samples) {
        @samples = glob("$inputdir/*_1*.fastq");
    }
}

if (!@samples) {
    @samples = glob("$inputdir/*/*.clean.fastq.gz");
    if (!@samples) {
        @samples = glob("$inputdir/*.clean.fastq.gz");
    }
}
if (!@samples) {
    @samples = glob("$inputdir/*/*.clean.fq.gz");
    if (!@samples) {
        @samples = glob("$inputdir/*.clean.fq.gz");
    }
}
if (!@samples) {
    @samples = glob("$inputdir/*/*_trimmed.fq.gz");
    if (!@samples) {
        @samples = glob("$inputdir/*_trimmed.fq.gz");
    }
}
chomp @samples;

if (!@samples) {
    log_message(":  Fastq file！");
    die " Fastq file，: $inputdir\n";
}

my %sample_info;
foreach my $fq1 (@samples) {

    my $base_name;
    my $fq2 = $fq1;


    my $filename = basename($fq1);
    
    if ($filename =~ /^(.+)_1_val_1\.fq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 =~ s/_1_val_1\.fq/_2_val_2.fq/;
    } elsif ($filename =~ /^(.+)_1_trimmed\.fq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 =~ s/_1_trimmed\.fq/_2_trimmed.fq/;
    } elsif ($filename =~ /^(.+)_1\.clean\.fq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 =~ s/_1\.clean\.fq/_2.clean.fq/;
    } elsif ($filename =~ /^(.+)_1\.fastq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 =~ s/_1\.fastq/_2.fastq/;
    } elsif ($filename =~ /^(.+)_1\.fq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 =~ s/_1\.fq/_2.fq/;
    } elsif ($filename =~ /^(.+)\.clean\.fastq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 = "";
    } elsif ($filename =~ /^(.+)\.clean\.fq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 = "";
    } elsif ($filename =~ /^(.+)_trimmed\.fq(?:\.gz)?$/) {
        $base_name = $1;
        $fq2 = "";
    } else {
        next;
    }
    
    $sample_info{$base_name} = ($fq2 ne "" && -e $fq2) ? { fq1 => $fq1, fq2 => $fq2, paired => 1 } : { fq1 => $fq1, paired => 0 };
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
log_message("Sample 1 ！\n");

sub process_sample {
    my ($sample_id, $fq1, $fq2, $paired) = @_;
    my $sample_dir = "$outputdir/$sample_id";
    make_path($sample_dir) unless -d $sample_dir;

    my $accepted_sam = "$sample_dir/accepted_hits.sam";
    my $accepted_sorted_bam = "$sample_dir/accepted_hits.sorted.bam";
    my $accepted_unique_bam = "$sample_dir/accepted_hits.sorted.unique.bam";
    my $metrics_file = "$sample_dir/${sample_id}.metricsFile";


    log_message(" Bowtie2 : $sample_id");
    my $bowtie_cmd = $paired 
        ? "bowtie2 -x $indexdir -p $threads -t -q -N 1 -L 25 --no-mixed --no-discordant --rg-id $sample_id --rg SM:$sample_id -1 $fq1 -2 $fq2 -S $accepted_sam"
        : "bowtie2 -x $indexdir -p $threads -t -q -N 1 -L 25 --no-mixed --no-discordant --rg-id $sample_id --rg SM:$sample_id -U $fq1 -S $accepted_sam";
    system($bowtie_cmd) == 0 or log_message(": Bowtie2 failed ($sample_id)");

    log_message(" BAM file: $sample_id");
    system("samtools view -bS $accepted_sam -o $accepted_sorted_bam && samtools sort -@ $threads -o $accepted_sorted_bam $accepted_sorted_bam") == 0 or log_message(": SAMtools failed ($sample_id)");

    log_message(" PCR : $sample_id");
    my $picard_cmd = (-e "$picarddir/picard.jar") ? "java -Xmx15g -jar $picarddir/picard.jar MarkDuplicates" : "picard MarkDuplicates";
    system("$picard_cmd I=$accepted_sorted_bam O=$accepted_unique_bam METRICS_FILE=$metrics_file REMOVE_DUPLICATES=true") == 0 or log_message(": Picard failed ($sample_id)");

    system("samtools index $accepted_unique_bam");

    log_message("file: $sample_id");
    unlink "$sample_dir/accepted_hits.sam", "$sample_dir/accepted_hits.sorted.bam", "$sample_dir/accepted_hits.sorted.bam.bai";
}