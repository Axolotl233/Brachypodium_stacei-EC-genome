#! perl

use warnings;
use strict;

while (<>){
    chomp;
    my @l = split/\t/;
    (my $n,my $m) = (split/,/,$l[5])[0,1];
    my $t = abs($n-$m);
    if($t > 0.4){
        print join"\t",@l[0,1,2,5,7];
        print "\n";
    }
}
