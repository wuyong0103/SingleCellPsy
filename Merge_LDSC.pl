#!/usr/bin/perl

system("sed '1d' /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/Mean/LDSC/Bed/LDSC_cell_types_GWAS_pvalues.txt | cut -f 1 | sort -u >Disease.txt");
my %ct;
#/home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/subtype.txt
open IN, $ARGV[0] or die;
while(<IN>){
    chomp;
    $ct{$_}=0;
}
close IN;

my %dc;
open DS, "Disease.txt" or die;
while(<DS>){
    chomp;
    foreach $key (keys %ct){
       $dc_key = $_."-$key";
       $dc{$dc_key} = $dc_key;
    }
}
close DS;

my $header = "Dis-CT";
#/home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/filename.txt
open IN, $ARGV[1] or die;
while(<IN>){
    chomp;
    my %exp;
    $header = $header."\t$_";
    my $file = "/home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/".$_."/LDSC/Bed/LDSC_cell_types_GWAS_pvalues.txt";
    open LDSC, $file or die;
    while(<LDSC>){
        chomp;
        my @a = split(/\t/, $_);
        my $key = $a[0]."-".$a[1];
        if(exists $exp{$key}){
            $exp{$key} = $exp{$key}."\t".$a[5];
        }
        else{
            $exp{$key} = $a[5];
        }
    }
    close LDSC;
    foreach $key (sort keys %dc){
        if(exists $exp{$key}){
            $dc{$key} = $dc{$key}."\t".$exp{$key};
        }
        else{
            $dc{$key} = $dc{$key}."\tNA";
        }
    }
}
close IN;

print "$header\n";
foreach $key (sort keys %dc){
    print $dc{$key}."\n";
}
