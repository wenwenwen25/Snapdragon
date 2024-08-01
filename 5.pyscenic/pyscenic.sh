pyscenic grn --num_workers 10 \
--sparse \
--method grnboost2 \
--output sce.adj.csv \
sce.loom \
mm_mgi_tfs.txt

pyscenic grn --num_workers 20 \
--sparse \
--method grnboost2 \
--output sce_TD.adj.csv \
TD_sce.loom \
mm_mgi_tfs.txt

pyscenic grn --num_workers 10 \
--sparse \
--method grnboost2 \
--output sce_TV.adj.csv \
TV_sce.loom \
mm_mgi_tfs.txt

pyscenic grn --num_workers 20 \
--sparse \
--method grnboost2 \
--output sce_WD.adj.csv \
WD_sce.loom \
mm_mgi_tfs.txt

pyscenic grn --num_workers 20 \
--sparse \
--method grnboost2 \
--output sce_WV.adj.csv \
WV_sce.loom \
mm_mgi_tfs.txt

pyscenic ctx --num_workers 10 \
--output WD_sce.regulons.csv \
--expression_mtx_fname WD_sce.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname ~/hh028/Tet2_DSS_scRNAseq/BM/240204/pyscenic/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
sce_WD.adj.csv \
~/hh028/Tet2_DSS_scRNAseq/BM/240204/pyscenic/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather


pyscenic aucell --num_workers 3 \
--output sce_SCENIC.loom \
WD_sce.loom \
WD_sce.regulons.csv


pyscenic grn --num_workers 10 \
--sparse \
--method grnboost2 \
--output LM_sce.adj.csv \
LM_sce.loom \
hs_hgnc_tfs.txt
