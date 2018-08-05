#!/usr/bin/perl -w
use strict;
use FastaUtils qw(readFasta);

my $reffasta = shift;
my $outfile  = shift;
my $hmmout   = shift;

$reffasta =~ /(.*)\.\S+/;
my $refbase = $1;

my $seq = readFasta($reffasta);
print "Finished reading the genome file\n";

open (OUT, ">$outfile") or die $!;
my $counter = 0;
open (HMM, "<$hmmout") or die $!;
while (<HMM>){
  chomp $_;
  next if (substr($_,0,1) eq "#");
  my @a = split (/\s+/, $_);
  if ($a[5] - $a[4] >= 500){
    my $start = ($a[11] eq "+") ? $a[8] - 501 : $a[9] - 501;
    my $end   = ($a[11] eq "+") ? $a[9] + 500 : $a[8] + 500;
    $start = 0 if ($start < 0);
    $end   = $a[10] if ($end > $a[10]);
    print OUT ">$a[0]:".($start+1)."-$end\n";
    print OUT uc(substr($$seq{$a[0]}{'sequence'}, $start, $end - $start))."\n";
    $counter++;
  }
}
close HMM;
close OUT;

exit;
