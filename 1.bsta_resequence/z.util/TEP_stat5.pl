#! perl

use warnings;
use strict;

my %r;
open IN,'<',"z.sample.list3";
while(<IN>){
    chomp;
    my @l = split/\t/;
    next if $l[1] eq "R";
    $r{$l[0]} = $l[1];
}
close IN;
my %r2;
my $c = 0;
my @d = "sample";
for my $f (@ARGV){
    open I,'<',"$f";
    while(<I>){
        chomp;
        my %t;
        my @l = split/\t/;
        my @n = split/,/,$l[6];
        $t{$_} = 1 for @n;
        $c += 1;
        push @d,$c;
        for my $s (keys %r){
            if(exists $t{$s}){
	push @{$r2{$s}} , 1;
            }else{
	push @{$r2{$s}} , 0;
            }
        }
    }
    close I;
}
push @d,"pop";
print join",",@d;
print "\n";
for my $k (sort {$a cmp $b} keys %r){
    my @t = @{$r2{$k}};
    print "$k,";
    print join",",@t;
    print ",$r{$k}\n";
}
