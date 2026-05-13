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
Rscript 04_Downsample.R

#get the top 10% cell type specificity genes
cd /home/lilab/wuyong/project/scRNA/MP-R1
mkdir Mean
cp ./downsampleMin/lineage_SumMean.csv Mean/Mean.csv
cd Mean
Rscript ../get_top10_gene.R Mean

cd /home/lilab/wuyong/project/scRNA/MP-R1
mkdir Mean_SubCellType
cp ./downsampleMin/sub_lineage_SumMean.csv Mean_SubCellType/Mean_SubCellType.csv
cd Mean_SubCellType
Rscript ../get_top10_gene.R Mean_SubCellType

#=====================================================================================================
#let's do LDSC
cd /home/lilab/wuyong/project/scRNA/MP-R1
bash get_annotation_ldscores_tissue_v2.sh

cd /home/lilab/wuyong/project/scRNA/MP-R1
for i in BD2019Stahl BD2021PGC3 BD2025OConnel SCZ2014PGC2 SCZ2019Lam SCZ2022PGC3
do
    bash get_partitioned_h2_tissue_v2.sh /home/lilab/wuyong/data/GWAS/${i}/${i}.sumstats.gz
done

cd /home/lilab/wuyong/project/scRNA/MP-R1
for k in Mean Mean_SubCellType
do
    cd /home/lilab/wuyong/project/scRNA/MP-R1/${k}/LDSC/Bed
    Rscript ../../../get_tissue_pvalue.R
done

cd /home/lilab/wuyong/project/scRNA/MP-R1
cd LDSC-result
echo 'Mean' | perl ../Merge_LDSC.pl ../subtype.txt - | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > LDSC_Subtype.txt
echo 'Mean_SubCellType' | perl ../Merge_LDSC.pl ../type.txt - | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > LDSC_Type.txt
cut -f 1,2 LDSC_Type.txt | sed 's/\-/\t/' | sed '1d' | perl -ane 'BEGIN{print "cell_type\tpvalue\tTrait\tmethod\n"}$F[1]=~s/\./-/g;print "$F[1]\t$F[2]\t$F[0]\tLDSC\n";' >LDSC-type.txt
cut -f 1,2 LDSC_Subtype.txt | sed 's/\-/\t/' | sed '1d' | perl -ane 'BEGIN{print "cell_type\tpvalue\tTrait\tmethod\n"}$F[1]=~s/\./-/g;print "$F[1]\t$F[2]\t$F[0]\tLDSC\n";' >LDSC-subtype.txt

#=====================================================================================================
#let's do MAGMA
cd /home/lilab/wuyong/project/scRNA/MP-R1
mkdir MAGMA-result; cd MAGMA-result
for k in Mean Mean_SubCellType
do
    for i in BD2019Stahl BD2021PGC3 BD2025OConnel SCZ2014PGC2 SCZ2019Lam SCZ2022PGC3
    do
        magma --gene-results /home/lilab/wuyong/data/GWAS/${i}/${i}.annotated_35kbup_10_down.genes.raw --set-annot /home/lilab/wuyong/project/scRNA/MP-R1/${k}/MAGMA/top10.txt --out ${i}-${k}
    done
done

cd /home/lilab/wuyong/project/scRNA/MP-R1/MAGMA-result
for i in BD2019Stahl BD2021PGC3 BD2025OConnel SCZ2014PGC2 SCZ2019Lam SCZ2022PGC3
do
    awk '{print "'${i}'""-Mean_SubCellType"}' | perl ../Merge_MAGMA.pl ../subtype.txt - | cut -f 1,3- >${i}_SubCellType_MAGMA.txt
    awk '{print "'${i}'""-Mean"}' | perl ../Merge_MAGMA.pl ../type.txt - | cut -f 1,3- >${i}_CellType_MAGMA.txt
done

find *SubCellType_MAGMA.txt | perl -ae '/(.*)_SubCellType_MAGMA.txt/; my $dis=$1."-"; open DIS, $_ or die; while(<DIS>){chomp; if($.==1){s/\Q$dis\E//g; print $_."\n";}else{print $dis.$_."\n";}}close DIS;' | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > MAGMA_Subtype.txt
find *_CellType_MAGMA.txt | perl -ae '/(.*)_CellType_MAGMA.txt/; my $dis=$1."-"; open DIS, $_ or die; while(<DIS>){chomp; if($.==1){s/\Q$dis\E//g; print $_."\n";}else{print $dis.$_."\n";}}close DIS;' | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > MAGMA_Type.txt
cut -f 1,2 MAGMA_Type.txt | sed 's/\-/\t/' | sed '1d' | perl -ane 'BEGIN{print "cell_type\tpvalue\tTrait\tmethod\n"}$F[1]=~s/\./-/g;print "$F[1]\t$F[2]\t$F[0]\tMAGMA\n";' >MAGMA-type.txt
cut -f 1,2 MAGMA_Subtype.txt | sed 's/\-/\t/' | sed '1d' | perl -ane 'BEGIN{print "cell_type\tpvalue\tTrait\tmethod\n"}$F[1]=~s/\./-/g;print "$F[1]\t$F[2]\t$F[0]\tMAGMA\n";' >MAGMA-subtype.txt

#=====================================================================================================
#let's do EWCE
cd /home/lilab/wuyong/project/scRNA/MP-R1
mkdir EWCE-result
Rscript 05_EWCE.R

#=====================================================================================================
#Seismic analysis
cd /home/lilab/wuyong/project/scRNA/MP-R1
mkdir Seismic-result
Rscript 06_Seismic.R /home/lilab/wuyong/project/scRNA/MP-R1/downsampleMin/downsampleMin_lineage.rds lineage None All_lineage
Rscript 06_Seismic.R /home/lilab/wuyong/project/scRNA/MP-R1/downsampleMin/downsampleMin_sublineage.rds sub_lineage None All_sublineage
cd Seismic-result
awk -vOFS="\t" 'BEGIN{print "cell_type\tpvalue\tTrait\tmethod"}{print $1,$2,$4,"Seismic"}' All_lineage_Asso.tsv | sed '2d' >Seismic-type.txt
awk -vOFS="\t" 'BEGIN{print "cell_type\tpvalue\tTrait\tmethod"}{print $1,$2,$4,"Seismic"}' All_sublineage_Asso.tsv | sed '2d' >Seismic-subtype.txt

#=====================================================================================================
#Prepare the plot file
cd /home/lilab/wuyong/project/scRNA/MP-R1
mkdir plot; cd plot
cat ../LDSC-result/LDSC-type.txt ../MAGMA-result/MAGMA-type.txt ../Seismic-result/Seismic-type.txt | awk 'NR==1 || $1!="cell_type"' >All-type.txt
cat ../LDSC-result/LDSC-subtype.txt ../MAGMA-result/MAGMA-subtype.txt ../Seismic-result/Seismic-subtype.txt | awk 'NR==1 || $1!="cell_type"' >All-subtype.txt

#plot
Rscript 01-Figure-DiffGWAS.R
Rscript 02-Figure-AllCortex.R
Rscript 03-Figure-AllCortex-EWCE.R
Rscript 04-Frc-Stage.R
Rscript 05-All-EWCE.R
Rscript 06-Frc-Stage-EWCE.R
Rscript 07-Upset-Enrichment.R
