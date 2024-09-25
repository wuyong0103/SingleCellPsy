#!/bin/bash
#2024-02-04
#written by Yong Wu
#wuyong0103@126.com

#Download the scRNA-seq data
cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science
#wget https://datasets.cellxgene.cziscience.com/e9e38c8d-fd21-475d-acf6-f8994b003d70.rds
for i in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz meta.tsv velmeshev_snRNA_seq.h5ad
do
    wget --no-check-certificate https://cells.ucsc.edu/pre-postnatal-cortex/all/rna/${i}
done

cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/Glia
for i in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz meta.tsv
do
    wget --no-check-certificate https://cells.ucsc.edu/pre-postnatal-cortex/glia/rna/${i}
done

cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/Excitatory
for i in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz meta.tsv
do
    wget --no-check-certificate https://cells.ucsc.edu/pre-postnatal-cortex/ex-neu/rna/${i}
done

cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/Interneuron
for i in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz meta.tsv
do
    wget --no-check-certificate https://cells.ucsc.edu/pre-postnatal-cortex/in/rna/${i}
done

cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/Microglia
for i in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz meta.tsv
do
    wget --no-check-certificate https://cells.ucsc.edu/pre-postnatal-cortex/mg/${i}
done

cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/Pericyte
for i in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz meta.tsv
do
    wget --no-check-certificate https://cells.ucsc.edu/pre-postnatal-cortex/per/${i}
done

cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/Vascular
for i in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz meta.tsv
do
    wget --no-check-certificate https://cells.ucsc.edu/pre-postnatal-cortex/end/${i}
done


#Add subtype information to meta data
perl -e 'my %id; open EX, $ARGV[0] or die; my $fl=<EX>; while(<EX>){chomp; my @a=split(/\t/, $_); $id{$a[0]}="EX_".$a[14];}close EX; open IN, $ARGV[1] or die; my $fl=<IN>; while(<IN>){chomp; my @a=split(/\t/, $_); $id{$a[0]}="In_".$a[14];}close IN; open GL, $ARGV[2] or die; my $fl=<GL>; while(<GL>){chomp; my @a=split(/\t/, $_); $id{$a[0]}=$a[15];}close GL; open MI, $ARGV[3] or die; my $fl=<MI>; while(<MI>){chomp; my @a=split(/\t/, $_); $id{$a[0]}="Micro";}close MI; open PE, $ARGV[4] or die; my $fl=<PE>; while(<PE>){chomp; my @a=split(/\t/, $_); $id{$a[0]}="Peri";}close PE; open VA, $ARGV[5] or die; my $fl=<VA>; while(<VA>){chomp; my @a=split(/\t/, $_); $id{$a[0]}="Vasc";}close VA; open ME, $ARGV[6] or die; my $fl=<ME>; chomp($fl); print "$fl\tsub_lineage\n"; while(<ME>){chomp; my @a=split(/\t/, $_); if(exists $id{$a[0]}){print "$_\t$id{$a[0]}\n";}else{print "$_\tNA\n";}}close ME;' ./Excitatory/meta.tsv Interneuron/meta.tsv Glia/meta.tsv Microglia/meta.tsv Pericyte/meta.tsv Vascular/meta.tsv meta.tsv | cut -f 1,6,7,8,13,15,16 | sed -e 's/0-1 years/years0_1/' -e 's/10-20 years/years10_20/' -e 's/1-2 years/years1_2/' -e 's/2-4 years/years2_4/' -e 's/2nd trimester/trimester2nd/' -e 's/3rd trimester/trimester3rd/' -e 's/4-10 years/years4_10/' -e '1s/age(days)/age_days/' >meta_subtype.tsv

#Calculate the mean expression of each cell type
python3 Calculate_MeanExpression.py /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/
for i in `find ./Mean*csv`; do sed -i -e '1s/\-/_/g' -e '1s/^/Gene/' ${i}; done

#Calculate the cell type specificity of each cell type
for i in `find ./Mean*csv | sed -e 's/\.\///' -e 's/\.csv//'`
do
    mkdir ${i}
    mv ${i}.csv ./${i}/
    cd ${i}
    Rscript ../Calculate_CTD.R ${i}
    cd ../
done

find ./Mean* -type d | sed 's/\.\///' | grep -v '/' >filename.txt

#Prepare the GWAS data, take example for schizophrenia
cd /home/lilab/wuyong/data/GWAS/SCZ2022PGC3
#For LDSC
awk 'NR>74' PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv | awk -vOFS="\t" 'BEGIN{print "SNP\tCHR\tPOS\tA1\tA2\tFreq\tINFO\tBeta\tSE\tP\tN_CAS\tN_CON"}{print $2,$1,$3,$4,$5,($6*$14+$7*$15)/($14+$15),$8,$9,$10,$11,$14,$15}' >SCZ2022PGC3.ldsc
$Munge --sumstats SCZ2022PGC3.ldsc --merge-alleles /home/lilab/software/ldsc-master/reference/w_hm3.snplist --N-cas-col N_CAS --N-con-col N_CON --frq Freq --out SCZ2022PGC3

#For MAGMA
awk -vOFS="\t" '{n=$11+$12;print $1,$10,n}' SCZ2022PGC3.ldsc | sed '1s/0/N/' >SCZ2022PGC3.pval
awk -vOFS="\t" '{n=$11+$12;print $1,$2,$3,$3}' SCZ2022PGC3.ldsc | sed '1d' >SCZ2022PGC3.bed
magma --annotate window=35,10 --snp-loc SCZ2022PGC3.bed --gene-loc /home/lilab/reference/hg19_gencode/NCBI37.3.gene.loc.extendedMHCexcluded --out SCZ2022PGC3.annotated_35kbup_10_down
magma --bfile /home/lilab/reference/1000G_Phase3_ForMiXeR/1000G.EUR.QC --pval SCZ2022PGC3.pval ncol=3 --gene-annot SCZ2022PGC3.annotated_35kbup_10_down.genes.annot --out SCZ2022PGC3.annotated_35kbup_10_down


cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science
#LDSC enrichment
#get_annotation_ldscores_tissue_v2.sh and get_partitioned_h2_tissue_v2.sh were got from https://github.com/jbryois/scRNA_disease/tree/master/Code_Paper/LDSC_pipeline; Please refer to paper Bryois et al., 2020, Nat Genet
bash get_annotation_ldscores_tissue_v2.sh

for i in AD2022Bellenguez ADHD2023Demontis ASD2019Grove BD2021PGC3 EatingDisorder2019Watson Insomina2022Watanabe HoardingSymptoms2022Strom MDD2019Howard OCT2018PGC SCZ2022PGC3 Suicide2023Docherty TouretteSyndrome2019Yu
do
    bash get_partitioned_h2_tissue_v2.sh /home/lilab/wuyong/data/GWAS/${i}/${i}.sumstats.gz
done

for k in `cat filename.txt`
do
    cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/${k}/LDSC/Bed
    Rscript ../../../get_tissue_pvalue.R
done


#MAGMA enrichment
mkdir MAGMA-result
cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/MAGMA-result
for k in `cat ../filename.txt`
do
    for i in AD2022Bellenguez ADHD2023Demontis ASD2019Grove BD2021PGC3 EatingDisorder2019Watson Insomina2022Watanabe HoardingSymptoms2022Strom MDD2019Howard OCT2018PGC SCZ2022PGC3 Suicide2023Docherty TouretteSyndrome2019Yu 
    do
        magma --gene-results /home/lilab/wuyong/data/GWAS/${i}/${i}.annotated_35kbup_10_down.genes.raw --set-annot /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/${k}/MAGMA/top10.txt --out ${i}-${k}
    done
done


cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/LDSC-result
cut -f 15 meta.tsv | sed '1d' | sort -u | awk '$1!="OUT"' >type.txt
cut -f 7 meta_subtype.tsv | sed '1d' | sort -u | awk '$1!="NA"' | sed 's/-/_/g' >subtype.txt

#Manage the MAGMA results
cd ~/project/scRNA/data/Velmeshev2023Science/MAGMA-result
for i in AD2022Bellenguez ADHD2023Demontis ASD2019Grove BD2021PGC3 EatingDisorder2019Watson Insomina2022Watanabe HoardingSymptoms2022Strom MDD2019Howard OCT2018PGC SCZ2022PGC3 Suicide2023Docherty TouretteSyndrome2019Yu Anxiety2021Forstner Intelligence2018Savage
do
    awk '$0~/SubCellType/{print "'${i}'""-"$0}' ../filename.txt | perl ../Merge_MAGMA.pl ../subtype.txt - | cut -f 1,3- >${i}_SubCellType_MAGMA.txt
    awk '$0!~/SubCellType/{print "'${i}'""-"$0}' ../filename.txt | perl ../Merge_MAGMA.pl ../type.txt - | cut -f 1,3- >${i}_CellType_MAGMA.txt
done
find *SubCellType_MAGMA.txt | perl -ae '/(.*)_SubCellType_MAGMA.txt/; my $dis=$1."-"; open DIS, $_ or die; while(<DIS>){chomp; if($.==1){s/\Q$dis\E//g; print $_."\n";}else{print $dis.$_."\n";}}close DIS;' | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > MAGMA_Subtype.txt
find *_CellType_MAGMA.txt | perl -ae '/(.*)_CellType_MAGMA.txt/; my $dis=$1."-"; open DIS, $_ or die; while(<DIS>){chomp; if($.==1){s/\Q$dis\E//g; print $_."\n";}else{print $dis.$_."\n";}}close DIS;' | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > MAGMA_Type.txt
cut -f 1,2 MAGMA_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Subtype.txt
cut -f 1,2 MAGMA_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Type.txt
cut -f 1,11 MAGMA_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_Male//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Type_Male.txt
cut -f 1,20 MAGMA_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_Female//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Type_Female.txt
cut -f 1,11 MAGMA_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType_Male//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Subtype_Male.txt
cut -f 1,20 MAGMA_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType_Female//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Subtype_Female.txt


#Manage the LDSC data
cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/LDSC-result
grep 'Sub' ../filename.txt | perl ../Merge_LDSC.pl ../subtype.txt - | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > LDSC_Subtype.txt
grep -v 'Sub' ../filename.txt | perl ../Merge_LDSC.pl ../type.txt - | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > LDSC_Type.txt
cut -f 1,2 LDSC_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Subtype.txt
cut -f 1,2 LDSC_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Type.txt
cut -f 1,11 LDSC_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_Male//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Type_Male.txt
cut -f 1,20 LDSC_Type.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_Female//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Type_Female.txt
cut -f 1,11 LDSC_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType_Male//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}'> All_Disorder_Subtype_Male.txt
cut -f 1,20 LDSC_Subtype.txt | perl ../Get_same_age.pl - | sed '1s/2[0-9]\{3\}[A-Za-z0-9]*_Mean_SubCellType_Female//g' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > All_Disorder_Subtype_Female.txt


#EWCE analysis
cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/EWCE
nohup bash Rscript EWCE_ctd.R >EWCE_ctd.out 2>EWCE_ctd.err &
gzip -dc ../features.tsv.gz >scRNA_All.txt
gzip -dc SCHEMA_gene_results.tsv.gz | cut -f 1,16 | sort -g -k2 | awk '$2!="NA"' | perl -e 'my %id; open IN, $ARGV[0] or die; while(<IN>){chomp; my @a=split(/\t/, $_); $id{$a[1]}=$a[4];}close IN; open IN, $ARGV[1] or die; while(<IN>){chomp; my @a=split(/\t/, $_); if(exists $id{$a[0]}){print "$a[0]\t$id{$a[0]}\t$a[1]\n";}else{print "$a[0]\t$a[0]\t$a[1]\n";}}close IN;' /home/lilab/reference/hg38_gencode/annotation_hg38_1.epi - | sed '1d' | awk 'BEGIN{print "ensembl\tsymbol\tp";}{print $0;}' >SCHEMA_Exome.txt
gzip -dc BipEx_gene_results.tsv.gz | awk -vFS="\t" '$2=="Bipolar Disorder"' | cut -f 1,15 | sort -g -k2 | awk '$2!="NA"' | perl -e 'my %id; open IN, $ARGV[0] or die; while(<IN>){chomp; my @a=split(/\t/, $_); $id{$a[1]}=$a[4];}close IN; open IN, $ARGV[1] or die; while(<IN>){chomp; my @a=split(/\t/, $_); if(exists $id{$a[0]}){print "$a[0]\t$id{$a[0]}\t$a[1]\n";}else{print "$a[0]\t$a[0]\t$a[1]\n";}}close IN;' /home/lilab/reference/hg38_gencode/annotation_hg38_1.epi - | awk 'BEGIN{print "ensembl\tsymbol\tp";}{print $0;}' >BipEx_Exome.txt
gzip -dc ASC_gene_results.tsv.gz | cut -f 1,15 | sort -g -k2 | perl -e 'my %id; open IN, $ARGV[0] or die; while(<IN>){chomp; my @a=split(/\t/, $_); $id{$a[1]}=$a[4];}close IN; open IN, $ARGV[1] or die; while(<IN>){chomp; my @a=split(/\t/, $_); if(exists $id{$a[0]}){print "$a[0]\t$id{$a[0]}\t$a[1]\n";}else{print "$a[0]\t$a[0]\t$a[1]\n";}}close IN;' /home/lilab/reference/hg38_gencode/annotation_hg38_1.epi - | sed '1d' | awk 'BEGIN{print "ensembl\tsymbol\tp";}{print $0;}' > ASC_Exome.txt

nohup Rscript EWCE_enrichment.R >EWCE_enrichment.out 2>EWCE_enrichment.err &

sed -i -e 's/"//g' -e 's/,/\t/g' ./result/*csv
for dis in ASD_Diff ASD_Exome ASD_SFARI BD_Diff BD_Exome MDD_Diff SCZ_Diff SCZ_Exome
do
    awk '{print "result/""'${dis}'""_"$0"_lvl2.csv"}' filename.txt | perl ../Merge_EWCE.pl ../subtype.txt - ${dis} lvl2
    awk '{print "result/""'${dis}'""_"$0"_lvl1.csv"}' filename.txt | perl ../Merge_EWCE.pl ../type.txt - ${dis} lvl1
done

cat *EWCE_lvl1_q.txt | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] == 0){print "\t<1E-4";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > EWCE_Type_q.txt
cat *EWCE_lvl2_q.txt | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] == 0){print "\t<1E-4";}elsif($F[$i]<0.001){printf "\t%.2e", $F[$i];}else{printf "\t%.3f", $F[$i];}}print "\n";}' > EWCE_Subtype_q.txt
cat *EWCE_lvl1_sd.txt | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}else{printf "\t%.3f", $F[$i];}}print "\n";}' > EWCE_Type_sd.txt
cat *EWCE_lvl2_sd.txt | awk 'NR==1 || $1!="CellType"' | perl -ane 'if($.==1){print;}else{print $F[0]; for(my $i=1; $i<@F; $i++){if($F[$i] eq "NA"){print "\tNA";}else{printf "\t%.3f", $F[$i];}}print "\n";}' > EWCE_Subtype_sd.txt

for i in p q fc sd
do
    cut -f 1,2 EWCE_Subtype_${i}.txt | perl ../Get_same_age.pl - >All_Disorder_Subtype_${i}.txt
    cut -f 1,2 EWCE_Type_${i}.txt | perl ../Get_same_age.pl - >All_Disorder_Type_${i}.txt
    cut -f 1,3 EWCE_Subtype_${i}.txt | perl ../Get_same_age.pl - | sed '1s/_male//g' >All_Disorder_Subtype_Male_${i}.txt
    cut -f 1,4 EWCE_Subtype_${i}.txt | perl ../Get_same_age.pl - | sed '1s/_female//g' >All_Disorder_Subtype_Female_${i}.txt
    cut -f 1,3 EWCE_Type_${i}.txt | perl ../Get_same_age.pl - | sed '1s/_male//g' >All_Disorder_Type_Male_${i}.txt
    cut -f 1,4 EWCE_Type_${i}.txt | perl ../Get_same_age.pl - | sed '1s/_female//g' >All_Disorder_Type_Female_${i}.txt
done



#2024-09-24
cd /home/lilab/wuyong/project/scRNA/data/Velmeshev2023Science/EWCE
nohup Rscript EWCE_enrichment_Drug.R >EWCE_enrichment_Drug.out 2>EWCE_enrichment_Drug.err &

sed -i -e 's/"//g' -e 's/,/\t/g' ./result/*Drug*
for dis in BD_Drug MDD_Drug SCZ_Drug
do
    awk '{print "result/""'${dis}'""_"$0"_lvl2.csv"}' filename.txt | perl ../Merge_EWCE.pl ../subtype.txt - ${dis} lvl2;
    awk '{print "result/""'${dis}'""_"$0"_lvl1.csv"}' filename.txt | perl ../Merge_EWCE.pl ../type.txt - ${dis} lvl1;
done
