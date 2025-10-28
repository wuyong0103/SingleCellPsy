#!/bin/bash
#2025-08-24
#Written by Yong Wu
#wuyong0103@126.com

cd /home/lilab/liming/wuyong/project/scRNA/integrate/downsample500

for i in Ex-L56NP In-SST
do
   mkdir ${i}
   cd ${i}
   Rscript ../Seurat2Exp.R -c ${i}
   python ../Exp2loom.py --celltype ${i}

    conda init
    conda activate pyscenic
    pyscenic grn --num_workers 10 \
             --output ${i}_grn.tsv \
             --method grnboost2 ${i}.loom /home/lilab/reference/SCENIC/allTFs_hg38.txt
    pyscenic ctx ${i}_grn.tsv /home/lilab/reference/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
             --annotations_fname /home/lilab/reference/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
             --expression_mtx_fname ${i}.loom \
             --mode "dask_multiprocessing" \
             --output ${i}_ctx.csv \
             --num_workers 10 \
             --mask_dropouts

    pyscenic aucell ${i}.loom ${i}_ctx.csv \
             --output ${i}_aucell.loom \
             --num_workers 10

    Rscript ../calcRSS_by_scenic.R --input_loom ${i}_aucell.loom --input_meta ../downsample500_sublineage.txt --stage stage --celltype ${i}

    Rscript ../plot.R --aucell ${i}_aucell.loom --rss ${i}_rss.rds --seurat ${i}_seurat.rds --annotation stage --grn ${i}_grn.tsv

    cd ../
done
