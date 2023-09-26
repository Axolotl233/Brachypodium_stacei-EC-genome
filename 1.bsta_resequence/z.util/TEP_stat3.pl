#! perl

use warnings;
use strict;

open IN,'<',"Bsta.v2.liffoff.mcscanx.gff";
my %h;
while(<IN>){
    chomp;
    my @l = split/\t/;
    $l[2] = $l[2] - 2000;
    $l[3] = $l[3] + 2000;
    $h{$l[0]}{$l[1]} = [$l[2],$l[3]];
}
close IN;

open IN,'<',shift;
while(<IN>){
    chomp;
    my @l = split/\t/;
    my @gs;
    for my $g(sort {$a cmp $b} keys %{$h{$l[0]}}){
        my @t = @{$h{$l[0]}{$g}};
        if($t[1] > $l[1] && $l[2] > $t[0]){
            push @gs , $g;
        }
    }
    if (scalar @gs == 0){
        @gs = qw /NA/;
    }
    print "$_\t",(join ",",@gs),"\n";
}
