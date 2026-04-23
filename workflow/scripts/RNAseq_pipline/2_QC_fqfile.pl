#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(mkpath);
use File::Basename;
use File::Spec;
use File::Find;


my ($outputdir, $inputdir, $threads, $max_parallel);

sub parse_command_line {
    GetOptions(
        "outputdir|o=s" => \$outputdir,
        "inputdir|i=s"  => \$inputdir,
        "threads=i"     => \$threads,
        "maxparallel=i" => \$max_parallel,
    ) or die "：$!";

    die "\n" unless defined $outputdir;
    die "\n" unless defined $inputdir;
}

sub set_default_values {
    $max_parallel  ||= 5;
    $threads       ||= 4;
}

parse_command_line();
set_default_values();

mkpath($outputdir) unless -d $outputdir;


my @samples;
find(sub {
    return unless -f $_;
    return unless /\.(fastq|fq|fastq\.gz|fq\.gz)$/;
    push @samples, $File::Find::name;
}, $inputdir);

if (scalar @samples == 0) {
    print "No fastq files found in: $inputdir. Skipping QC...\n";
    mkpath($outputdir) unless -d $outputdir;
    my $finished_flag = "$outputdir/QC_finished.txt";
    open(my $FH, '>', $finished_flag) or die "Creatingfile: $!";
    print $FH "All QC jobs finished (no files).\n";
    close($FH);
    exit 0;
}

my @active_processes;
my $qc_report_dir = "$outputdir/qc_reports";
mkpath($qc_report_dir) unless -d $qc_report_dir;

my %processed_samples;

foreach my $file (@samples) {
    my ($sample_id, $read_type);
    

    my $rel_path = File::Spec->abs2rel($file, $inputdir);
    my $rel_dir = dirname($rel_path);
    
    my $clean_rel_dir = $rel_dir;
    $clean_rel_dir =~ s{/rawdata/}{/}g;
    $clean_rel_dir =~ s{^rawdata/}{}g;
    $clean_rel_dir =~ s{/rawdata$}{}g;
    
    if ($file =~ /(.+)_1\.(fastq|fq|fastq\.gz|fq\.gz)$/) {
        $sample_id = basename($1);
        $read_type = 'paired';
    } elsif ($file =~ /(.+)_2\.(fastq|fq|fastq\.gz|fq\.gz)$/) {
        $sample_id = basename($1);
        $read_type = 'paired';
    } elsif ($file =~ /(.+)\.(fastq|fq|fastq\.gz|fq\.gz)$/) {
        $sample_id = basename($1);
        $read_type = 'single';
    } else {
        next;
    }
    

    my $unique_id = "$clean_rel_dir/$sample_id";
    next if $processed_samples{$unique_id};
    $processed_samples{$unique_id} = 1;

    print "Starting QC for sample: $unique_id ($read_type)\n";
    
    my $target_out_dir = "$outputdir/$clean_rel_dir";
    mkpath($target_out_dir) unless -d $target_out_dir;


    my $report_subdir = "$qc_report_dir/$clean_rel_dir";
    mkpath($report_subdir) unless -d $report_subdir;
    my $html_report = "$report_subdir/${sample_id}_fastp.html";
    my $json_report = "$report_subdir/${sample_id}_fastp.json";

    my $pid = fork();
    if (!defined $pid) {
        die "Creating: $!";
    } elsif ($pid == 0) {
        my $cmd;
        if ($read_type eq 'paired') {

            my $fq1 = (glob("$inputdir/$rel_dir/${sample_id}_1.fastq*") || glob("$inputdir/$rel_dir/${sample_id}_1.fq*"))[0];
            my $fq2 = (glob("$inputdir/$rel_dir/${sample_id}_2.fastq*") || glob("$inputdir/$rel_dir/${sample_id}_2.fq*"))[0];
            

            if (!$fq1 || !-e $fq1) { $fq1 = "$inputdir/$rel_dir/${sample_id}_1.fastq.gz" }
            if (!$fq1 || !-e $fq1) { $fq1 = "$inputdir/$rel_dir/${sample_id}_1.fastq" }
            if (!$fq2 || !-e $fq2) { $fq2 = "$inputdir/$rel_dir/${sample_id}_2.fastq.gz" }
            if (!$fq2 || !-e $fq2) { $fq2 = "$inputdir/$rel_dir/${sample_id}_2.fastq" }

            my $fq1_out = "$target_out_dir/${sample_id}_1.clean.fastq.gz";
            my $fq2_out = "$target_out_dir/${sample_id}_2.clean.fastq.gz";
            
            if (-e $fq1 && -e $fq2) {
                $cmd = "fastp -i $fq1 -I $fq2 -o $fq1_out -O $fq2_out -h $html_report -j $json_report --thread $threads";
                if (system($cmd) != 0) {
                    print "WARNING: fastp failed for $sample_id. Creating dummy outputs.\n";
                    system("touch $fq1_out $fq2_out");
                }
            } else {
                die "Sample $sample_id corresponding to fq file ($fq1  $fq2 does not exist)\n";
            }
        } else {

            my $fq = (glob("$inputdir/$rel_dir/${sample_id}.fastq*") || glob("$inputdir/$rel_dir/${sample_id}.fq*"))[0];
            
            if (!$fq || !-e $fq) { $fq = "$inputdir/$rel_dir/${sample_id}.fastq.gz" }
            if (!$fq || !-e $fq) { $fq = "$inputdir/$rel_dir/${sample_id}.fastq" }

            my $fq_out = "$target_out_dir/${sample_id}.clean.fastq.gz";
            
            if (-e $fq) {
                $cmd = "fastp -i $fq -o $fq_out -h $html_report -j $json_report --thread $threads";
                if (system($cmd) != 0) {
                    print "WARNING: fastp failed for $sample_id. Creating dummy outputs.\n";
                    system("touch $fq_out");
                }
            } else {
                die "Sample $sample_id corresponding to fq file ($fq does not exist)\n";
            }
        }
        
        exit(0);
    } else {
        push @active_processes, $pid;
    }

    while (@active_processes >= $max_parallel) {
        my $child = waitpid(-1, 0);
        @active_processes = grep { $_ != $child } @active_processes;
    }
}

while (wait() != -1) { }

my $finished_flag = "$outputdir/QC_finished.txt";
open(my $FH, '>', $finished_flag) or die "Creatingfile: $!";
print $FH "All QC jobs finished.\n";
close($FH);
