#!/usr/bin/perl
use strict;
use warnings;

die "\nUsage:\n\tperl $0 <sample_info> <raw_dir>\n\n" unless(@ARGV == 2);

my ($list, $indir) = @ARGV;

open(IN, "<$list") || die "$!";
while(<IN>){
	chomp;
	next if(/^SampleID/);
	my @tmp = split /\t/, $_;
	my $id = $tmp[0];
	my $name = $tmp[1];
	unless(-e "$indir/neg/$name.raw"){
		`mv $indir/neg/$id.raw $indir/neg/$name.raw`;
		`mv $indir/pos/$id.raw $indir/pos/$name.raw`;	
	}
}
close(IN);
