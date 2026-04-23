#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(mkpath);
use YAML 'LoadFile';

my ($motifdir, $outputdir, $tfmapfile, $motifmapfile, $geneoutput, $groupfile, $bwfile, $outpic, $help, $beddir, $groupname);

GetOptions(
    "motifdir=s"      => \$motifdir,
    "tfmapfile=s"     => \$tfmapfile,
    "motifmapfile=s"  => \$motifmapfile,
    "geneoutput=s"    => \$geneoutput,
    "outputdir=s"     => \$outputdir,
    "groupfile=s"     => \$groupfile,
    "bwfile=s"        => \$bwfile,
    "outpic=s"        => \$outpic,
    "beddir=s"        => \$beddir,
    "groupname=s"     => \$groupname,
    "help!"           => \$help,
);

exit 1 if $help or !($motifdir && $outputdir && $groupfile && $bwfile);

my %tf_files = (
    'homo' => "workflow/resources/Motif_tf_anno/Homo_sapiens_TF.txt",
    'mm' => "workflow/resources/Motif_tf_anno/Mus_musculus_TF.txt",
);

my $group_config = LoadFile($groupfile);
my $genome = $group_config->{species} || "mm";

if (!$tfmapfile) {
    exit 1 unless exists $tf_files{$genome};
    $tfmapfile = $tf_files{$genome};
}

my @types = ("up", "down");
$outpic //= $outputdir;
mkpath($outpic, 0755) unless -e $outpic;

foreach my $type (@types) {

    my $posfile_pattern = "$beddir/*_${type}_peaks.bed";
    my @posfiles = glob($posfile_pattern);
    next unless @posfiles;
    my $motif_result_dir = "${outputdir}/Motif_result_${type}";
    unless (-e $motif_result_dir) {
        mkpath($motif_result_dir, 0755) or die "$motif_result_dir\n";
    }
    foreach my $posfile (@posfiles) {
        my $prefix = "${groupname}_Enhancer_both_${type}_genes";
        my $fastafile = "$outputdir/${prefix}_target_seqs.fasta";
        extract_fasta_by_genome($groupfile, $posfile, $fastafile);
        motif_enrichment($fastafile, $motif_result_dir, $motifdir);
        my $networkfile = "$outputdir/${groupname}_TF_motif_network_${type}.txt";
        build_motif_network_with_significance($motif_result_dir, $networkfile);
        if ($tfmapfile && $motifmapfile) {
            my $geneout = $geneoutput // "$outputdir/${groupname}_TF_gene_mapping_${type}.txt";
            build_tf_gene_mapping($tfmapfile, $motifmapfile, $networkfile, $geneout);
        }
        process_bam_bw_files($groupfile, $bwfile, $posfile, $prefix, $outpic, $outputdir) if ($groupfile && $bwfile);
    }
}

sub extract_fasta_by_genome {
    my ($groupfile, $posfile, $fastafile) = @_;
    my $group_config = LoadFile($groupfile);
    my %genome_files = (
        'mm' => "workflow/resources/genome/GRCm38.p6.genome.fa.gz",
        'homo' => "workflow/resources/genome/GRCh38.p12.genome.fa.gz",
    );
    my $genome = $group_config->{species} || "mm";
    exit 1 unless exists $genome_files{$genome};
    my $genomefile = $genome_files{$genome};
    open(GENOME, "gzip -dc $genomefile|") or die "$genomefile: $!";
    my (%chromosome, $chromname);
    while (<GENOME>) {
        chomp;
        if (/^>(chr\S*)/) {
            $chromname = $1;
            $chromosome{$chromname} = '';
        } else {
            $chromosome{$chromname} .= $_;
        }
    }
    close GENOME;
    my %posHash;
    open(ERNA, "<$posfile") or die "$posfile: $!";
    while (<ERNA>) {
        chomp;
        my @fieldValues = split /\s+/,$_;
        $posHash{$fieldValues[3]} = [@fieldValues[0,1,2]];
    }
    close ERNA;
    open(FA, ">$fastafile") or die "$fastafile: $!";
    foreach my $eRNAid (keys %posHash) {
        my ($chr, $start, $end) = @{ $posHash{$eRNAid} };
        my $geneseq = substr($chromosome{$chr}, $start, $end - $start + 1);
        print FA ">$eRNAid\n";
        my $strnum = length($geneseq);
        for (my $i = 0; $i < $strnum; $i += 50) {
            print FA substr($geneseq, $i, 50), "\n";
        }
    }
    close FA;
}

sub process_bam_bw_files {
    my ($groupfile, $bwfile, $posfile, $prefix, $outpic, $outputdir) = @_;
    my $bw_config = LoadFile($bwfile);
    my @bw_files = ();
    foreach my $sample_name (keys %$bw_config) {
        if (exists $bw_config->{$sample_name}->{bw1}) {
            push @bw_files, $bw_config->{$sample_name}->{bw1};
        }
    }
    my $samples_label = join(" ", keys %$bw_config);
    exit 1 unless scalar(@bw_files) == scalar(split(" ", $samples_label));
    my @soft_colors = (
        "#66c2a5", "#8da0cb", "#fc8d62", "#e78ac3", "#ffd92f", "#a6d854", "#d9d9d9",
        "#b3b3b3", "#ccebc5", "#fdb462", "#bebada", "#fb8072", "#80b1d3", "#b3de69"
    );
    my @plot_colors = map { $soft_colors[$_ % @soft_colors] } (0..$#bw_files);
    my $color_str = join(" ", @plot_colors);
    my $matrix_file = "$outputdir/${prefix}_matrix.gz";
    my $signal_file = "$outputdir/${prefix}_signal.tab";
    my $plot_file = "$outpic/${prefix}_plotProfile.pdf";
    unless (-e $matrix_file) {
        my $computeMatrix_cmd = "computeMatrix reference-point -R $posfile -S @bw_files -b 1000 -a 1000 --binSize 10 -p 10 --samplesLabel $samples_label -out $matrix_file --outFileNameMatrix $signal_file";
        system($computeMatrix_cmd);
    }
    if (-e $matrix_file) {
        my $plot_cmd = "plotProfile -m $matrix_file -out $plot_file --plotType fill --perGroup --colors $color_str --plotTitle ${prefix}_from_peaks --samplesLabel $samples_label --plotHeight 6 --plotWidth 15";
        system($plot_cmd);
    }
}

sub motif_enrichment {
    my ($fastafile, $outputdir, $motifdir) = @_;
    my $motif_result_dir = "$outputdir/Motif_result";
    unless(-e $motif_result_dir){
        mkpath($motif_result_dir, 0755) or die "$motif_result_dir\n";
    }
    my @motiffiles = `find $motifdir -name "*meme"`;
    chomp @motiffiles;
    foreach my $motiffile (@motiffiles){
        $motiffile =~ /.*\/(.*)\.meme/;
        my $motifname = $1;
        my $outdir = "$motif_result_dir/$motifname";
        unless(-e $outdir){
            mkpath($outdir, 0755) or die "$outdir\n";
        }
        my $ame_cmd = "ame --control --shuffle-- --oc $outdir $fastafile $motiffile";
        my $out = system($ame_cmd);
    }
}

sub build_motif_network_with_significance {
    my ($motif_result_dir, $outfile) = @_;
    open(OUT, ">$outfile") or die "$outfile: $!\n";
    my @motifseqfiles = `find $motif_result_dir -type f -name "sequences.tsv"`;
    chomp @motifseqfiles;
    foreach my $motifseqfile (@motifseqfiles){
        (my $amefile = $motifseqfile) =~ s/sequences/ame/;
        next unless -e $amefile;
        open(my $AMETSV, "<", $amefile) or die "$amefile: $!\n";
        my ($stat, $adj_pvalue, $motif_Name) = ('No', '', '');
        while (<$AMETSV>) {
            chomp;
            next if /^rank/;
            my @fields = split /\t/;
            if(defined $fields[6] && $fields[6] < 0.05) {
                $stat = 'Yes';
            }
            $adj_pvalue = defined $fields[6] ? $fields[6] : '';
            $motif_Name  = defined $fields[3] ? $fields[3] : '';
            last;
        }
        close $AMETSV;
        open(my $MS, "<", $motifseqfile) or die "$motifseqfile: $!\n";
        while(<$MS>){
            chomp;
            my @fieldValues = split /\t/;
            if(defined $fieldValues[6] && $fieldValues[6] eq "tp" && $fieldValues[3] !~ /shuf/) {
                print OUT join("\t",
                    $fieldValues[1],   # eRNA_ID
                    $motif_Name,       # motif_name
                    $fieldValues[3],   # motif_name
                    $fieldValues[5],   # score or other field
                    $stat,             # Yes/No
                    $adj_pvalue
                ), "\n";
            }
        }
        close $MS;
    }
    close OUT;
}

sub build_tf_gene_mapping {
    my ($tfmapfile, $motifmapfile, $tfeRNAfile, $outputfile) = @_;
    my (%tfmapHash, %tfnameHash);
    open(my $MAP, "<", $tfmapfile) or die "$tfmapfile: $!\n";
    while (<$MAP>) {
        chomp;
        my @fieldValues = split /\t/;
        next unless @fieldValues >= 3;
        if ($fieldValues[0] eq 'Homo_sapiens' or $fieldValues[0] eq 'Mus_musculus') {
            $tfmapHash{uc($fieldValues[1])} = $fieldValues[2];
            $tfnameHash{$fieldValues[2]} = uc($fieldValues[1]);
        }
    }
    close $MAP;
    my %motifmapHash;
    open(my $MOTIF, "<", $motifmapfile) or die "$motifmapfile: $!\n";
    while (<$MOTIF>) {
        chomp;
        my @fieldValues = split /\t/;
        if ($fieldValues[1] =~ /\::/) {
            my @parts = split /::/, $fieldValues[1];
            foreach my $part (@parts) {
                $motifmapHash{uc($part)} = $fieldValues[0];
            }
        } else {
            $motifmapHash{uc($fieldValues[1])} = $fieldValues[0];
        }
    }
    close $MOTIF;
    my %tfmotifHash;
    foreach my $tfname (keys %tfmapHash) {
        if (exists $motifmapHash{$tfname}) {
            $tfmotifHash{$motifmapHash{$tfname}} = $tfmapHash{$tfname};
        }
    }
    my %genetfHash;
    open(my $TE, "<", $tfeRNAfile) or die "$tfeRNAfile: $!\n";
    while (<$TE>) {
        chomp;
        my @fieldValues = split /\t/;
        if (exists $tfmotifHash{$fieldValues[0]}) {
            $genetfHash{$fieldValues[1]}{ $tfmotifHash{$fieldValues[0]} } = 1;
        }
    }
    close $TE;
    open(my $OUT, ">", $outputfile) or die "$outputfile: $!\n";
    foreach my $gene (keys %genetfHash) {
        my $tfs = join(", ", sort keys %{ $genetfHash{$gene} });
        print $OUT "$gene\t$tfs\n";
    }
    close $OUT;
}