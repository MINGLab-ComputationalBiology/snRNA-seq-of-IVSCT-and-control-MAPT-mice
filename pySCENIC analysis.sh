#!/bin/bash

WORKDIR="snSCENIC"
LOOM_FILE="$WORKDIR/Inh_pyscenic.loom"
TF_FILE="$WORKDIR/allTFs_mm.txt"
RANKINGS_DB="$WORKDIR/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOTATIONS_DB="$WORKDIR/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

ADJ_FILE="$WORKDIR/adj.inh.csv"
REGULONS_FILE="$WORKDIR/regulons.inh.csv"
AUCELL_LOOM="$WORKDIR/aucell_output.inh.loom"

pyscenic grn \
    --num_workers 16 \
    --output $ADJ_FILE \
    --method grnboost2 \
    $LOOM_FILE \
    $TF_FILE

pyscenic ctx \
    $ADJ_FILE \
    $RANKINGS_DB \
    --annotations_fname $ANNOTATIONS_DB \
    --expression_mtx_fname $LOOM_FILE \
    --output $REGULONS_FILE \
    --num_workers 16

pyscenic aucell \
    $LOOM_FILE \
    $REGULONS_FILE \
    --output $AUCELL_LOOM \
    --num_workers 16
