#!/bin/bash
#Written by Yong Wu
#wuyong0103@126.com


#Create Seurat object for each study
#Run script in each directory
Rscript Bakken2021_Seurat.R
Rscript Batiuk2022_Seurat.R
Rscript Caglayan2023_Seurat.R
Rscript Clarence2025_Seurat.R
Rscript Emani2024_Seurat.R
Rscript Fan2020_Seurat.R
Rscript Huuki-Myers2024_Seurat.R
Rscript Jorstad2023_seurat.R
Rscript Maitra2023_Seurat.R
Rscript Mannens2024_Seurat.R
Rscript Morabito2021_Seurat.R
Rscript Pfisterer2020_Seurat.R
Rscript Pineda2024_Seurat.R
Rscript Schirmer2019_Seurat.R
Rscript Steyn2024_Seurat.R
Rscript Tran2021_Seurat.R
Rscript Velmeshev2023_Seurat.R
Rscript Wang2025_Seurat.R
Rscript Zhu2023_Seurat.R

#Integrage all dataset
Rscript 01_Integrate.R
Rscript 02_ExNeu.R
Rscript 03_InNeu.R

#Downsample
04_Downsample.R

#Create EWCE object
05_EwceCtd.R

cd /home/lilab/wuyong/project/scRNA/data
cat Bakken2021Nature/meta_Bakken2021.txt Batiuk2022SciAdv/meta_Batiuk2022.txt Caglayan2023Nature/meta_Caglayan2023.txt Cameron2023BioPsy/meta_Cameron2023.txt Clarence2025NatGenet/meta_Clarence2025.txt Emani2024Science/meta_Emani2024.txt Emani2024Science/meta_Ma2022.txt Fan2020SciAdv/meta_Fan2020.txt Gerstner2025SciAdv/meta_Gerstner2025.txt Huuki-Myers2024Science/meta_Huuki-Myers2024.txt Hwang2025Nature/meta_Hwang2025.txt Jorstad2023Science/meta_Jorstad2023.txt Maitra2023NatCommun/meta_Maitra2023.txt Mannens2024Nature/meta_Mannens2024.txt Morabito2021NatNeu/meta_Morabito2021.txt Pfisterer2020NC/meta_Pfisterer2020.txt Pineda2024Cell/meta_Pineda2024.txt Polioudakis2019Neuron/meta_Polioudakis2019.txt Schirmer2019Nature/meta_Schirmer2019.txt Smith2021PNAS/meta_Smith2021_8m.txt Smith2021PNAS/meta_Smith2021_pcw22.txt Steyn2024NatGenet/meta_Steyn2024.txt Tran2021Neuron/meta_Tran2021.txt Velmeshev2023Science/meta_Velmeshev2023.txt Wang2025Nature/meta_Wang2025.txt Zhu2023SciAdv/meta_Zhu2023.txt | awk 'NR==1 || $1!="barcode"' > ../integrate/meta_AllStudy.txt

cd /home/lilab/wuyong/project/scRNA/integrate
perl -ane 'if($.==1){chomp; print "$_\tstage\n";}else{chomp; if($F[3]=~/y$/){$F[3]=~s/y//; if($F[3]<=3){print "$_\tinfant\n"}elsif($F[3]<=6){print "$_\tkid\n"}elsif($F[3]<=12){print "$_\tchild\n"}elsif($F[3]<=18){print "$_\tjuvenile\n"}elsif($F[3]<=40){print "$_\tyouth\n"}elsif($F[3]<=60){print "$_\tmidlife\n"}else{print "$_\telder\n"}}elsif($F[3]=~/pcw$/){$F[3]=~s/pcw$//; if($F[3]<13){print "$_\tfirsttrim\n"}elsif($F[3]<28){print "$_\tsecondtrim\n"}else{print "$_\tthirdtrim\n"}}}' meta_AllStudy.txt >meta_AllStudy_stage.txt

cd /home/lilab/wuyong/project/scRNA/integrate
cut -f 7 Integrate-Meta.tsv | sed '1d' | sort -u | awk '$1!="Unknown"' | sed 's/-/./g' >type.txt
cut -f 8 Integrate-Meta.tsv | sed '1d' | sort -u | awk '$1!="Unknown"' | sed 's/-/./g' >subtype.txt


#=====================================================================================================
cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
mkdir Mean
mv Mean.csv Mean/Mean.csv
cd Mean
Rscript ../get_top10_gene.R Mean

cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
mkdir Mean_SubCellType
mv Mean_SubCellType.csv Mean_SubCellType/Mean_SubCellType.csv
cd Mean_SubCellType
Rscript ../get_top10_gene.R Mean_SubCellType

cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
for i in firsttrim secondtrim thirdtrim infant kid child juvenile youth midlife elder;
do
    mkdir  Mean_${i};
    mv stage_lineage_${i}_SumMean.csv Mean_${i}/Mean_${i}.csv;
    cd Mean_${i}/;
    Rscript ../get_top10_gene.R Mean_${i};
    cd ../;

    mkdir Mean_SubCellType_${i};
    mv stage_sub_lineage_${i}_SumMean.csv Mean_SubCellType_${i}/Mean_SubCellType_${i}.csv;
    cd Mean_SubCellType_${i}/;
    Rscript ../get_top10_gene.R Mean_SubCellType_${i};
    cd ../;
done

find ./Mean* -type d | sed 's/\.\///' | grep -v '/' >filename.txt


#=====================================================================================================
#let's do MAGMA
cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
mkdir MAGMA-result; cd MAGMA-result
for k in `cat ../filename.txt`
do
    for i in BD2025OConnel SCZ2022PGC3
    do
        magma --gene-results /home/lilab/wuyong/data/GWAS/${i}/${i}.annotated_35kbup_10_down.genes.raw --set-annot /home/lilab/wuyong/project/scRNA/integrate/downsample500/${k}/MAGMA/top10.txt --out ${i}-${k}
    done
done

cd /home/lilab/wuyong/project/scRNA/integrate/downsample500/MAGMA-result
for i in firsttrim secondtrim thirdtrim infant kid child juvenile youth midlife elder
do
    sed -i "s/^${i}\.//" *_${i}.gsa.out
done

#Manage the MAGMA results
cd /home/lilab/wuyong/project/scRNA/integrate/downsample500/MAGMA-result
for i in BD2025OConnel SCZ2022PGC3
do
    awk '$0~/SubCellType/{print "'${i}'""-"$0}' ../filename.txt | perl ../Merge_MAGMA.pl ../subtype.txt - | cut -f 1,3- >${i}_SubCellType_MAGMA.txt
    awk '$0!~/SubCellType/{print "'${i}'""-"$0}' ../filename.txt | perl ../Merge_MAGMA.pl ../type.txt - | cut -f 1,3- >${i}_CellType_MAGMA.txt
done

find *SubCellType_MAGMA.txt | perl -ae '/(.*)_SubCellType_MAGMA.txt/; my $dis=$1."-"; open DIS, $_ or die; while(<DIS>){chomp; if($.==1){s/\Q$dis\E//g; print $_."\n";}else{print $dis.$_."\n";}}close DIS;' | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > MAGMA_Subtype.txt
find *_CellType_MAGMA.txt | perl -ae '/(.*)_CellType_MAGMA.txt/; my $dis=$1."-"; open DIS, $_ or die; while(<DIS>){chomp; if($.==1){s/\Q$dis\E//g; print $_."\n";}else{print $dis.$_."\n";}}close DIS;' | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > MAGMA_Type.txt
cut -f 1,2 MAGMA_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Subtype.txt
cut -f 1,2 MAGMA_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Type.txt


#=====================================================================================================
#let's do LDSC
cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
for i in firsttrim secondtrim thirdtrim infant kid child juvenile youth midlife elder
do
    cd Mean_${i}/LDSC/Bed
    rename "s/^${i}\.//" ${i}.*
    cd ../../../Mean_SubCellType_${i}/LDSC/Bed
    rename "s/^${i}\.//" ${i}.*
    cd ../../../
done

cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
for i in {1..4}
do
    nohup bash get_annotation_ldscores_tissue_v2_${i}.sh >get_annotation_ldscores_tissue_v2_${i}.out 2>get_annotation_ldscores_tissue_v2_${i}.err &
done

cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
for i in BD2025OConnel SCZ2022PGC3
do
    bash get_partitioned_h2_tissue_v2.sh /home/lilab/wuyong/data/GWAS/${i}/${i}.sumstats.gz
done


cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
for k in `cat filename.txt`
do
    cd /home/lilab/wuyong/project/scRNA/integrate/downsample500/${k}/LDSC/Bed
    Rscript ../../../get_tissue_pvalue.R
done

#Manage the LDSC data
cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
mkdir LDSC-result; cd LDSC-result
grep 'Sub' ../filename.txt | perl ../Merge_LDSC.pl ../subtype.txt - | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.4f", $F[$i];}}print "\n";}' > LDSC_Subtype.txt
grep -v 'Sub' ../filename.txt | perl ../Merge_LDSC.pl ../type.txt - | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.4f", $F[$i];}}print "\n";}' > LDSC_Type.txt
cut -f 1,2 LDSC_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.4f", $F[$i];}}print "\n";}' > All_Disorder_Subtype.txt
cut -f 1,2 LDSC_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.4f", $F[$i];}}print "\n";}' > All_Disorder_Type.txt


#=====================================================================================================
#let's do EWCE
cd /home/lilab/wuyong/project/scRNA/integrate/downsample500
mkdir EWCE-result; cd EWCE-result

gzip -dc SCHEMA_gene_results.tsv.gz | cut -f 1,16 | sort -g -k2 | awk '$2!="NA"' | perl -e 'my %id; open IN, $ARGV[0] or die; while(<IN>){chomp; my @a=split(/\t/, $_); $id{$a[1]}=$a[4];}close IN; open IN, $ARGV[1] or die; while(<IN>){chomp; my @a=split(/\t/, $_); if(exists $id{$a[0]}){print "$a[0]\t$id{$a[0]}\t$a[1]\n";}else{print "$a[0]\t$a[0]\t$a[1]\n";}}close IN;' /home/lilab/reference/hg38_gencode/annotation_hg38_1.epi - | sed '1d' | awk 'BEGIN{print "ensembl\tsymbol\tp";}{print $0;}' >SCHEMA_Exome.txt
gzip -dc BipEx_gene_results.tsv.gz | awk -vFS="\t" '$2=="Bipolar Disorder"' | cut -f 1,15 | sort -g -k2 | awk '$2!="NA"' | perl -e 'my %id; open IN, $ARGV[0] or die; while(<IN>){chomp; my @a=split(/\t/, $_); $id{$a[1]}=$a[4];}close IN; open IN, $ARGV[1] or die; while(<IN>){chomp; my @a=split(/\t/, $_); if(exists $id{$a[0]}){print "$a[0]\t$id{$a[0]}\t$a[1]\n";}else{print "$a[0]\t$a[0]\t$a[1]\n";}}close IN;' /home/lilab/reference/hg38_gencode/annotation_hg38_1.epi - | awk 'BEGIN{print "ensembl\tsymbol\tp";}{print $0;}' >BipEx_Exome.txt

nohup Rscript EWCE_downsample500.R >EWCE_downsample500.out 2>EWCE_downsample500.err &

sed -i 's/\./_/g' ../subtype.txt
sed -i 's/\./_/g' ../type.txt
sed -i -e 's/"//g' -e 's/,/\t/g' *csv
for dis in BD_Exome SCZ_Exome
do
    awk '{print "'${dis}'""_downsample500_lvl2_"$0".csv"}' filename.txt | perl ../Merge_EWCE.pl ../subtype.txt - ${dis} lvl2
    awk '{print "'${dis}'""_downsample500_lvl1_"$0".csv"}' filename.txt | perl ../Merge_EWCE.pl ../type.txt - ${dis} lvl1
done

for i in p q fc sd
do
    cat *EWCE_lvl1_${i}.txt | awk 'NR==1 || $1!="CellType"' >EWCE_Type_${i}.txt
    cat *EWCE_lvl2_${i}.txt | awk 'NR==1 || $1!="CellType"' >EWCE_Subtype_${i}.txt
    cut -f 1,2 EWCE_Subtype_${i}.txt | perl ../Get_same_age.pl - >All_Disorder_Subtype_${i}.txt
    cut -f 1,2 EWCE_Type_${i}.txt | perl ../Get_same_age.pl - >All_Disorder_Type_${i}.txt
done

sed -i 's/_/./g' ../subtype.txt
sed -i 's/_/./g' ../type.txt

#Run SCENIC
bash 06_SCENIC.sh

#Figures
Rscrit 07_Figure-LDSC-MAGMA.R
Rscrit 08_Figure-EWCE.R
Rscript 09_DEAndClusterGVis.R
