#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($enhfile, $genefile, $outfile, $randbed, $help);
GetOptions(
    "enhfile=s"   => \$enhfile,
    "genefile=s"  => \$genefile,
    "outfile=s"   => \$outfile,
    "randbed=s"   => \$randbed,
    "help!"       => \$help,
);

exit(1) if $help or not ($enhfile and $genefile and $outfile and $randbed);

sub min {
    my $mn = shift;
    for (@_) { $mn = $_ if ($_ < $mn); }
    return $mn;
}

sub read_enh {
    my ($file) = @_;
    my %hash;
    open my $fh, '<', $file or die "Cannot open enhfile $file: $!\n";
    while (<$fh>) {
        chomp; next if /^\s*$/;
        my @f = split /\t/;
        $hash{$f[3]} = [$f[0], $f[1]+0, $f[2]+0];
    }
    close $fh;
    return %hash;
}

sub read_gene {
    my ($file) = @_;
    my %hash;
    open my $fh, '<', $file or die "Cannot open genefile $file: $!\n";
    while (<$fh>) {
        chomp; next if /^\s*$/;
        my @f = split /\t/;
        $hash{$f[0]}{$f[3]} = [$f[1]+0, $f[2]+0];
    }
    close $fh;
    return %hash;
}

my %enhHash = read_enh($enhfile);
my %geneposHash = read_gene($genefile);

open my $out_fh,  '>', $outfile or die "Cannot write outfile $outfile: $!\n";
open my $rand_fh, '>', $randbed  or die "Cannot write randbed $randbed: $!\n";

for my $enhid (keys %enhHash) {
    my ($enhchr, $enhstart, $enhend) = @{ $enhHash{$enhid} };
    next unless exists $geneposHash{$enhchr};
    for my $geneid (keys %{ $geneposHash{$enhchr} }) {
        my ($genestart, $geneend) = @{ $geneposHash{$enhchr}{$geneid} };
        my $dist = min(abs($genestart - $enhend), abs($geneend - $enhstart));
        if ($dist < 1_000_000) {
            print $out_fh "$enhchr\t$enhstart\t$enhend\t$enhid\t$geneid\t$dist\n";
            my $elength = $enhend - $enhstart;
            next if $elength <= 0;
            my $win_start = $genestart - 1_000_000;
            $win_start = 0 if $win_start < 0;
            my $win_end = $geneend + 1_000_000;
            my $max_start = $win_end - $elength;
            next if $max_start <= $win_start;
            my $rand_start = $win_start + int(rand($max_start - $win_start + 1));
            my $rand_end = $rand_start + $elength;
            print $rand_fh "$enhchr\t$rand_start\t$rand_end\tRAND_${enhid}_$geneid\n";
        }
    }
}
close $out_fh;
close $rand_fh;