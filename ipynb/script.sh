#!/bin/sh

source /home/chang/miniconda3/etc/profile.d/conda.sh
conda activate SComatic

SCOMATIC=$(pwd)
sample=$1
REF=~/refdata-gex-GRCh38-2020-A/fasta/genome.fa
output_dir=$SCOMATIC/scomatic_out
META=$SCOMATIC/cell_clus.tsv

mkdir -p $output_dir

output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1


python ~/SComatic/scripts/SplitBam/SplitBamCellTypes.py --bam $SCOMATIC/possorted_genome_bam.bam \
        --meta $META \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1


output_dir2=$output_dir/Step2_BaseCellCounts
mkdir -p $output_dir2

for bam in $(ls -d $output_dir1/*bam);do
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp
  # Command line to submit to cluster
  python ~/SComatic/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 16
  rm -rf $temp
done


output_dir3=$output_dir/Step3_BaseCellCountsMerged
mkdir -p $output_dir3

python ~/SComatic/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv

# Step 4.1
output_dir4=$output_dir/Step4_VariantCalling
mkdir -p $output_dir4

python ~/SComatic/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF


# Step 4.2
editing=~/SComatic/RNAediting/AllEditingSites.hg38.txt
PON=~/SComatic/PoNs/PoN.scRNAseq.hg38.tsv

python ~/SComatic/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${sample}.calling.step1.tsv \
          --outfile ${output_dir4}/${sample} \
          --editing $editing \
          --pon $PON

bedtools intersect -header -a ${output_dir4}/${sample}.calling.step2.tsv -b ~/SComatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/${sample}.calling.step2.pass.tsv

STEP4_2_pass=${output_dir4}/${sample}.calling.step2.pass.tsv

output_dir7=$output_dir/SingleCellAlleles
mkdir -p $output_dir7

for bam in $(ls -d $output_dir1/*bam);do  
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
    temp=$output_dir7/temp_${cell_type}
    mkdir -p $temp
    python ~/SComatic/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $bam  \
        --infile ${STEP4_2_pass}   \
        --nprocs 16   \
        --meta $META   \
        --outfile ${output_dir7}/${cell_type}.single_cell_genotype.tsv  \
        --tmp_dir $temp  \
        --ref $REF
    rm -rf $temp
done
