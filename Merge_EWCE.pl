#!/usr/bin/perl
#read the cell type into a hash
my %idp;
my %idfc;
my %idsd;
my %idq;
open IN, $ARGV[0] or die;
while(<IN>){
    chomp;
    $idp{$_} = "";
    $idfc{$_} = "";
    $idsd{$_} = "";
    $idq{$_} = "";
}
close IN;

$po = $ARGV[2]."_EWCE_".$ARGV[3]."_p.txt";
$fco = $ARGV[2]."_EWCE_".$ARGV[3]."_fc.txt";
$sdo = $ARGV[2]."_EWCE_".$ARGV[3]."_sd.txt";
$qo = $ARGV[2]."_EWCE_".$ARGV[3]."_q.txt";

open PO, ">$po" or die;
open FCO, ">$fco" or die;
open SDO, ">$sdo" or die;
open QO, ">$qo" or die;
print PO "CellType";
print FCO "CellType";
print SDO "CellType";
print QO "CellType";

open IN, $ARGV[1] or die;
while(<IN>){
    chomp;
    m/result\/[A-Z]*?_[a-zA-Z]*?_(.*)_lvl[12]\.csv/;
    print PO "\t$1";
    print FCO "\t$1";
    print SDO "\t$1";
    print QO "\t$1";
    open FILE, $_ or die "can not open file:$_\n";
    my %p;
    my %fc;
    my %sd;
    my %q;
    while(<FILE>){
        my @a = split(/\s+/, $_);
        $p{$a[0]} = $a[1];
        $fc{$a[0]} = $a[2];
        $sd{$a[0]} = $a[3];
        $q{$a[0]} = $a[4];
    }
    foreach $key (keys %idp){
        if(exists $p{$key}){
            $idp{$key} = $idp{$key}."\t".$p{$key};
            $idfc{$key} = $idfc{$key}."\t".$fc{$key};
            $idsd{$key} = $idsd{$key}."\t".$sd{$key};
            $idq{$key} = $idq{$key}."\t".$q{$key};
        }
        else{
            $idp{$key} = $idp{$key}."\tNA";
            $idfc{$key} = $idfc{$key}."\tNA";
            $idsd{$key} = $idsd{$key}."\tNA";
            $idq{$key} = $idq{$key}."\tNA";
        }
    }
    close FILE;
}
close IN;
print PO "\n";
print FCO "\n";
print SDO "\n";
print QO "\n";

#output the result
foreach $key (sort keys %idp){
    print PO $ARGV[2]."-".$key.$idp{$key}."\n";
    print FCO $ARGV[2]."-".$key.$idfc{$key}."\n";
    print SDO $ARGV[2]."-".$key.$idsd{$key}."\n";
    print QO $ARGV[2]."-".$key.$idq{$key}."\n";
}
close PO;
close FCO;
close SDO;
close QO;

