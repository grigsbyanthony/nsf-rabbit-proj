> [!note]
> Performed 2025-7-21
> Data processing and documentation was performed on 2025-7-21 following the completion of the master metadata file. **This pipeline uses QIIME2 v2025.4**.

```
. 
└── documentation/ 
	└── data_processing.md
```

## Importing and demultiplexing

```
 qiime tools import \
   --type EMPPairedEndSequences \
   --input-path sequences \
   --output-path multiplexed_data.qza

 qiime demux emp-paired \
   --i-seqs multiplexed_data.qza \
   --m-barcodes-file demux_s64.tsv \
   --m-barcodes-column BarcodeSequence \
   --p-rev-comp-mapping-barcodes \
   --p-rev-comp-barcodes \
   --o-per-sample-sequences demux1.qza \
   --o-error-correction-details demux1_details.qza
   
qiime demux summarize \
  --i-data demux1.qza \
  --o-visualization demux1.qzv
```

```
 qiime tools import \
   --type EMPPairedEndSequences \
   --input-path sequences \
   --output-path multiplexed_data.qza

 qiime demux emp-paired \
   --i-seqs multiplexed_data.qza \
   --m-barcodes-file demux_s65.tsv \
   --m-barcodes-column BarcodeSequence \
   --p-rev-comp-mapping-barcodes \
   --p-rev-comp-barcodes \
   --o-per-sample-sequences demux2.qza \
   --o-error-correction-details demux2_details.qza
   
qiime demux summarize \
  --i-data demux2.qza \
  --o-visualization demux2.qzv
```

## Denoising

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux1.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 215 \
  --p-trunc-len-r 218 \
  --o-table table-s64.qza \
  --o-representative-sequences rep-seqs-s64.qza \
  --o-denoising-stats denoising-stats-s64.qza
```

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux2.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 215 \
  --p-trunc-len-r 218 \
  --o-table table-s65.qza \
  --o-representative-sequences rep-seqs-s65.qza \
  --o-denoising-stats denoising-stats-s65.qza
```

## Denoising stats summary

> [!warning]
> Concerning metadata.tsv
> `rab38_Ac_Swab_Tr1` and `rab55_Fr_Swab_Tr2` were both duplicated IDs. `rab38_Ac_Swab_Tr1` had identical entires, while `rab55_Fr_Swab_Tr2` had entries with different `Sample Barcodes`. One was run on sr64 and the other on sr65, so the the former was omitted from the metadata.tsv.

```
qiime feature-table summarize \
  --i-table table-s64.qza \
  --o-visualization table-s64.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-s64.qza \
  --o-visualization rep-seqs-s64.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats-s64.qza \
  --o-visualization denoising-stats-s64.qzv
```

```
qiime feature-table summarize \
  --i-table table-s65.qza \
  --o-visualization table-s65.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-s65.qza \
  --o-visualization rep-seqs-s65.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats-s65.qza \
  --o-visualization denoising-stats-s65.qzv
```

## Merging feature tables and representative sequences

> [!warning]
> The following samples were filtered from `sr64` since they were redundant
> ```
> cat > samples-to-exclude.txt << EOF
> sample-id
> rab38_Ac_Swab_Tr1
> rab55_Fr_Swab_Tr2
> EOF
> ```
> 
> ```
> qiime feature-table filter-samples \
>   --i-table table-s64.qza \
>   --m-metadata-file samples-to-exclude.txt \
>   --p-exclude-ids \
>   --o-filtered-table filtered-table-s64.qza
> ```

```
qiime feature-table merge \
  --i-tables filtered-table-s64.qza \
  --i-tables table-s65.qza \
  --o-merged-table merged_table.qza

qiime feature-table merge-seqs \
  --i-data rep-seqs-s64.qza \
  --i-data rep-seqs-s65.qza \
  --o-merged-data merged_rep-seqs.qza
  
qiime feature-table summarize \
  --i-table merged_table.qza \
  --o-visualization merged_table.qzv \
  --m-sample-metadata-file metadata.tsv
```

## Taxonomy w/ greengenes2

```
qiime feature-classifier classify-sklearn \
  --i-classifier greengenes2-2024.09-515-806-nb-classifier.qza \
  --i-reads merged_rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

## Filtering

```
qiime taxa filter-table \
  --i-table merged_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,sp004296775 \
  --o-filtered-table merged-table-nmnc.qza

# Filtering out archaea
qiime taxa filter-table \
  --i-table merged-table-nmnc.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude archaea \
  --o-filtered-table merged-table-nmncna.qza

# Filtering out unassigned
qiime taxa filter-table \
    --i-table merged-table-nmncna.qza \
    --i-taxonomy taxonomy.qza \
    --p-exclude "Unassigned" \
    --o-filtered-table merged-table-filtered.qza

qiime feature-table filter-features \
    --i-table merged-table-filtered.qza \
    --p-min-frequency 10 \
    --p-min-samples 2 \
    --o-filtered-table merged-table-filtered-all.qza
```

## Taxonomy visualization 

```
qiime taxa barplot \
  --i-table merged-table-filtered-all.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

## Phylogenetic tree creation w/ `fasttree`

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```
