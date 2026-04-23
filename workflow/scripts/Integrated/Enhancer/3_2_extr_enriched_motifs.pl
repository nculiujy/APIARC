#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my ($tfmapfile, $motifdir, $outfile, $yamlfile, $help);

GetOptions(
    "tfmapfile=s" => \$tfmapfile,
    "motifdir=s"  => \$motifdir,
    "outfile=s"   => \$outfile,
    "yamlfile=s"  => \$yamlfile,
    "help!"       => \$help,
);

exit 1 if $help;
die "Missing motifdir or yamlfile" if !$motifdir or !$yamlfile;

$motifdir = abs_path($motifdir) or die "Cannot resolve abs_path for $motifdir: $!";
my $motifdir_name = basename($motifdir); 
my $parent_dir = dirname($motifdir);  

if (!$tfmapfile) {
    my $genome = "mm";
    open(my $IN, "<", $yamlfile) or die "Cannot open $yamlfile: $!";
    while(<$IN>) {
        if (/^species:\s*['"]?(\w+)['"]?/) {
            $genome = lc($1);
            last;
        }
    }
    close $IN;
    my %tf_files = (
        'homo' => "workflow/resources/Motif_tf_anno/Homo_sapiens_TF.txt",
        'mm' => "workflow/resources/Motif_tf_anno/Mus_musculus_TF.txt",
    );
    if (!exists $tf_files{$genome}) {
        die "Unknown genome species: $genome\n";
    }
    $tfmapfile = $tf_files{$genome};
}

if (!defined $outfile) {
    $outfile = "$parent_dir/${motifdir_name}_list.txt";
}

my %tfmapHash;
open(MAP, "<$tfmapfile") or die "$tfmapfile: $!\n";
while (<MAP>) {
    chomp;
    if (/^Mus_musculus/) {
        my @fieldValues = split /\t/;
        $tfmapHash{uc($fieldValues[1])} = $fieldValues[2];
    }
}
close MAP;

my @enfiles = `find $motifdir -name "ame.tsv"`;
open(OUT, ">$outfile") or die "$outfile: $!\n";
foreach my $enfile (@enfiles) {
    chomp $enfile;
    next if $enfile =~ /^\s*$/;
    open(EN, "<$enfile") or die "$enfile: $!\n";
    while (<EN>) {
        chomp;
        next if /^rank/;
        my @fieldValues = split /\t/;
        if (defined $fieldValues[6] && $fieldValues[6] < 1) {
            if (exists $tfmapHash{uc($fieldValues[3])}) {
                my $tfid   = $tfmapHash{uc($fieldValues[3])};
                my $tfname = uc($fieldValues[3]);
                print OUT "$tfid\t$tfname\t$fieldValues[6]\n";
            }
        }
    }
    close EN;
}
close OUT;
