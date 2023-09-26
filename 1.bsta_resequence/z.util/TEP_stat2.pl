#! perl

my %h;
open IN,'<',shift;
readline IN;
my $c = 0;
while(<IN>){
    $c += 1;
    my @l = split/\t/;
    my @j = split/,/,$l[4];
    for my $b (@j){
        my @t = split/\//,$b;
        $h{$t[0]} += 1;
    }
}

for my $k(sort {$a cmp $b} keys %h){
    my $r = ($h{$k})/$c;
    $r = $r *100;
    #my $r = "a";
    print "$k\t$h{$k}\t$r\n";
}
