#!/usr/bin/perl
BEGIN { $ENV{PATH} = "/home/jyliu/miniconda3/envs/ChIPseq_Pipline/bin:$ENV{PATH}"; }
use Getopt::Long;
use File::Path qw(make_path);
use File::Basename;
use Time::Piece;


my ($metadata_file, $bamdir, $outdir, $picture_dir, $genome, $qval, $threads, $cutoff_type, $norm_method, $help);
GetOptions(
    'metadata=s'    => \$metadata_file,
    'bamdir=s'      => \$bamdir,
    'outdir=s'      => \$outdir,
    'picdir=s'      => \$picture_dir,
    'genome=s'      => \$genome,
    'qval=f'        => \$qval,
    'cutoff_type=s' => \$cutoff_type,
    'norm_method=s' => \$norm_method,
    'threads=i'     => \$threads,
    'help!'         => \$help,
) or die "。\n";


if ($help || !$metadata_file || !$bamdir || !$outdir || !$picture_dir || !$genome) {
    die "Usage: perl 4_1_peakcalling.pl --metadata <file> --bamdir <dir> --outdir <dir> --picdir <dir> --genome <str> --qval <float> --cutoff_type <qvalue|pvalue> --norm_method <BPM|TPM|none> --threads <num>\n";
}

$qval //= 0.05;
$cutoff_type //= "qvalue";
$norm_method //= "none";
$threads //= 8;

make_path($outdir) unless -d $outdir;
make_path($picture_dir) unless -d $picture_dir;

sub log_message {
    my ($message) = @_;
    my $time = localtime->strftime('%Y-%m-%d %H:%M:%S');
    print "[$time] $message\n";
}

log_message("ChIP-seq 4: Peak Calling ");
log_message(": cutoff_type=$cutoff_type, value=$qval, normalization=$norm_method");



open(my $IN, '<', $metadata_file) or die "file $metadata_file: $!";
my $header = <$IN>; # Skip header

my %group_samples;

while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;

    $line =~ s/^\s+|\s+$//g;
    $line =~ s/"//g;
    
    my ($ip_sample, $input_sample, $ip_name) = split(/,/, $line);
    next unless $ip_sample && $ip_name;


    my $group = $ip_name;
    $group =~ s/_rep\d+$//;

    push @{$group_samples{$group}{"IP_samples"}}, $ip_sample;
    push @{$group_samples{$group}{"Input_samples"}}, $input_sample;
}
close($IN);

my $processed_count = 0;

foreach my $group (keys %group_samples) {
    my $treat_srrs_ref = $group_samples{$group}{"IP_samples"};
    my $control_srrs_ref = $group_samples{$group}{"Input_samples"};
    

    my %seen_treat;
    my @treat_srrs = grep { !$seen_treat{$_}++ } @$treat_srrs_ref;
    
    my %seen_control;
    my @control_srrs = grep { !$seen_control{$_}++ } @$control_srrs_ref;
    

    my @treat_bams;
    foreach my $t_srr (@treat_srrs) {
        my $bam = "$bamdir/$t_srr/accepted_hits.sorted.unique.bam";
        push @treat_bams, $bam if -e $bam;
        log_message(":  BAM does not exist: $bam") unless -e $bam;
    }
    
    my @control_bams;
    foreach my $c_srr (@control_srrs) {
        my $bam = "$bamdir/$c_srr/accepted_hits.sorted.unique.bam";
        push @control_bams, $bam if -e $bam;
        log_message(":  BAM does not exist: $bam") unless -e $bam;
    }
    
    if (!@treat_bams || !@control_bams) {
        log_message(":  $group  BAM file，...");
        next;
    }

    $processed_count++;
    my $pair_name = "${group}_vs_Input";
    log_message(">>> : $pair_name | genome: $genome");

    my $sample_dir = "$outdir/$pair_name";
    make_path($sample_dir);

    #### ===[Step 1: MACS2 peak calling]===
    log_message(">>> [MACS2]  Peak Calling: $pair_name");
    
    my $treat_bams_str = join(" ", @treat_bams);
    my $control_bams_str = join(" ", @control_bams);
    

    my $cutoff_arg = ($cutoff_type eq "pvalue") ? "-p $qval" : "-q $qval";
    


    my $extra_args = "";
    if ($norm_method eq "BPM" || $norm_method eq "RPKM") {
        $extra_args .= " --SPMR";
    }
    
    my $macs2_genome = $genome;
    if (lc($genome) eq "tair") {
        $macs2_genome = "1.2e8";
    } elsif (lc($genome) eq "homo") {
        $macs2_genome = "hs";
    } elsif (lc($genome) eq "mm") {
        $macs2_genome = "mm";
    }

    my $macs2_cmd = "macs2 callpeak -t $treat_bams_str -c $control_bams_str -g $macs2_genome -n $pair_name --keep-dup all $cutoff_arg $extra_args --outdir $sample_dir";
    log_message(">>> [MACS2] Cmd: $macs2_cmd");
    if (system($macs2_cmd) != 0) {
        log_message("WARNING: MACS2 failed ($pair_name). Creating dummy outputs.");
        system("touch $sample_dir/${pair_name}_peaks.narrowPeak");
        system("touch $sample_dir/${pair_name}_peaks.xls");
        system("touch $sample_dir/${pair_name}_summits.bed");
    }


    my $narrow   = "$sample_dir/${pair_name}_peaks.narrowPeak";
    my $bed      = "$sample_dir/${pair_name}_peaks.bed";
    my $tab_bed  = "$sample_dir/${pair_name}_peaks_tab.bed";

    if (-e $narrow) {
        system("awk '{print \$1, \$2, \$3, \$4, \$5}' $narrow > $bed");
        system("sed 's/ \\+/\t/g' $bed > $tab_bed");
    } else {
        log_message(":  MACS2 file $narrow");
    }
}

if ($processed_count == 0) {
    log_message(": SampleCannot find BAM file。Peakcalling failed！");
    die " BAM file。\n";
}

log_message("Sample Peak Calling ！\n");