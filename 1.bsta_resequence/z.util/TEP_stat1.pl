#! perl

use warnings;
use strict;

my %h;
while(<>){
    chomp;
    my @l = split/\t/;
    my @t2 = split/,/,$l[-1];
    my $n = scalar @t2;
    $h{$n} += 1;
}
print "$_\t$h{$_}\n" for sort{$a <=>$b} keys %h;
