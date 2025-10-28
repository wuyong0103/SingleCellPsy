#!/usr/bin/perl
my %exp;
my $header = "CellType";
my $age;
my %ha;
open IN, $ARGV[0] or die;
while(<IN>){
    chomp;
    my @a = split(/\t/, $_);
    if($.==1){
        $age = $a[1];
    }
    else{
        $a[0] =~ /(.*)\-(.*)/;
        my $dis = $1;
        if(!exists $ha{$dis}){
            $header = $header."\t".$dis."_".$age;
            $ha{$dis} = 0;
        }
        my $cell = $2;
        $exp{$cell} = $exp{$cell}."\t".$a[1];
    }

}
close IN;

print "$header\n";
foreach $key (sort keys %exp){
    print "$key$exp{$key}\n";
}
