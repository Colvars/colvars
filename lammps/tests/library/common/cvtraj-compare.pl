#!/usr/bin/perl -w
# Tool to validate and compare two LAMMPS data files
# with "inexact" floating point comparisons
# July 2013 by Axel Kohlmeyer <akohlmey@gmail.com>

use strict;
use warnings;

my $version = 'v0.1';

# delta for floating point comparisons.
my $small = 1.0e-7;

########################################################################
# floating point comparison
sub floatdiff {
  my ($n1,$n2,$rel) = @_;

  my $diff = abs($n1-$n2);
  my $avg = (abs($n1)+abs($n2))*0.5;
  return 0 if ($avg == 0.0);
  if ($rel) {
#    print "relative difference: ",$diff/$avg," vs. $small\n";
    return 0 if ($diff/$avg < $small);
  } else {
#    print "absolute difference: ",$diff," vs. $small\n";
    return 0 if ($diff < $small);
  }
  return 1;
}

########################################################################
# compare colvars trajectory file
sub compare_files {
    my ($fp1,$fp2) = @_;
    my (@l1,@l2,$i,$line,$num,$a,$b);
    $line = 0;

    while (! (eof $fp1 && eof $fp2)) {
        die "Truncated file." if (eof $fp1 || eof $fp2);

        ++$line;
        # the first line has headers with column assignment.
        $_ = <$fp1>;
        chomp;
        @l1 = split;

        $_ = <$fp2>;
        chomp;
        @l2 = split;

        # skip over comment lines
        next if (($l1[0] eq '#') || ($l2[0] eq '#'));

        # first column is step number
        die "Step number mismatch" if ($l1[0] != $l2[0]);

        # the rest is floating point. do inexact bsolute compare
        $num = $#l1;
        for ($i=1; $i <= $num; ++$i) {
            $a = $l1[$i];
            $b = $l2[$i];
            die "Data mismatch in line $line column $i: $a vs. $b"
                if floatdiff($a,$b,0);
        }
    }
}

########################################################################
# main program

my ($f1,$f2);

if ($#ARGV < 1) {
  die "usage $0 <file 1> <file 2>";
}

print "\nColvars trajectory file validation tool. $version\n\n";

# open files
open($f1, '<', $ARGV[0]) or die $!;
print "opened file 1: $ARGV[0]\n";
open($f2, '<', $ARGV[1]) or die $!;
print "opened file 2: $ARGV[1]\n";

compare_files($f1,$f2);

print "Files $ARGV[0] and $ARGV[1] match\n\n";

close $f1;
close $f2;

exit 0;

