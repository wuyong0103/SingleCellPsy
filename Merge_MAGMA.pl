#!/usr/bin/perl

#read the cell type into a hash
my %id;
open IN, $ARGV[0] or die;
while(<IN>){
    chomp;
    $id{$_} = 0;
}
close IN;
print "CellType\t0";

open IN, $ARGV[1] or die;
while(<IN>){
    chomp;
    my $file = $_.".gsa.out";
    print "\t$_";
    open FILE, $file or die;
    my %temp;
    while(<FILE>){
        my @a = split(/\s+/, $_);
        $temp{$a[0]} = $a[6];
    }
    foreach $key (keys %id){
        if(exists $temp{$key}){
            $id{$key} = $id{$key}."\t".$temp{$key};
        }
        else{
            $id{$key} = $id{$key}."\tNA"
        }
    }
    close FILE;
}
close IN;
print "\n";

#output the result
foreach $key (sort keys %id){
    print "$key\t$id{$key}\n";
}
