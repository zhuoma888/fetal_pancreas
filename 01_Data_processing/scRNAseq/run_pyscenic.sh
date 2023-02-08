#！ /bin/bash
## pyscenic shell script

#load database
dir=/home/marzon/pyscenic
tfs=$dir/allTFs_hg38.txt
feather1=$dir/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
feather2=$dir/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
tbl=$dir/motifs-v9-nr.hgnc.tbl

#input loom file
input_loom=$dir/EP.loom
ls $tfs $feather1 $feather1 $tbl

pyscenic grn \
	--num_workers 10 \
	--output adj.sample.tsv \
	--method grnboost2 \
	$input_loom $tfs

pyscenic ctx \
	adj.sample.tsv $feather1 $feather2 \
	--annotations_fname $tbl \
	--expression_mtx_fname $input_loom \
	--mode "dask_multiprocessing" \
	--output reg.csv \
	--num_workers 10 \
	--mask_dropouts

pyscenic aucell \
	$input_loom \
	reg.csv \
	--output out_SCENIC_EP.loom \
	--num_workers 10
