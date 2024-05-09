#!/bin/bash
#SBATCH -c 6
#SBATCH -t 0-4:00
#SBATCH -p priority
#SBATCH --mem=32G
#SBATCH -o brh040_%j.out
#SBATCH -e brh040_%j.err
#SBATCH --mail-type=END,FAIL

## load modules and virtual environment
module load gcc/9.2.0
module load python/3.9.14
source ~/pythonvirtualenv/venv_velocity/bin/activate

#download velocity index
kb ref -d linnarsson -i /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/index.idx -g /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/t2g.txt -c1 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/spliced_t2c.txt -c2 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/unspliced_t2c.txt

#generate velocity count matrices - for 1st donor (HT2), BroadRun230120 and BroadRun230330
kb count --verbose -i /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/index.idx -g /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/t2g.txt -x 10xv3 -m 30G --workflow lamanno --h5ad -c1 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/spliced_t2c.txt -c2 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/unspliced_t2c.txt -o kb_count_out_ht2 --overwrite \
/n/groups/cbdm_lab/scRNAseq/BroadRun230120/230120_SL-NVM_0989_BHN3TGDRX2/BroadRun230120/GEX_Fastqs_DC7/outs/fastq_path/HN3TGDRX2/GEX_S1_L002_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230120/230120_SL-NVM_0989_BHN3TGDRX2/BroadRun230120/GEX_Fastqs_DC7/outs/fastq_path/HN3TGDRX2/GEX_S1_L002_R2_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD1/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD1/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R2_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD1/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD1/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R2_001.fastq.gz

#generate velocity count matrices - for 2nd donor (HT3), BroadRun230120 and BroadRun 230330
kb count --verbose -i /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/index.idx -g /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/t2g.txt -x 10xv3 -m 30G --workflow lamanno --h5ad -c1 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/spliced_t2c.txt -c2 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/unspliced_t2c.txt -o kb_count_out_ht3 --overwrite \
/n/groups/cbdm_lab/scRNAseq/BroadRun230120/230120_SL-NVM_0989_BHN3TGDRX2/BroadRun230120/GEX_Fastqs_DC8/outs/fastq_path/HN3TGDRX2/GEX_S1_L002_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230120/230120_SL-NVM_0989_BHN3TGDRX2/BroadRun230120/GEX_Fastqs_DC8/outs/fastq_path/HN3TGDRX2/GEX_S1_L002_R2_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD2/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD2/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R2_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD2/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD2/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R2_001.fastq.gz

#generate velocity count matrices - for 3rd donor (HT4), BroadRun230330
kb count --verbose -i /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/index.idx -g /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/t2g.txt -x 10xv3 -m 30G --workflow lamanno --h5ad -c1 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/spliced_t2c.txt -c2 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/unspliced_t2c.txt -o kb_count_out_ht4 --overwrite \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD3/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD3/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R2_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD3/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD3/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R2_001.fastq.gz

#generate velocity count matrices - for 4th donor (HT5), BroadRun230330
kb count --verbose -i /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/index.idx -g /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/t2g.txt -x 10xv3 -m 30G --workflow lamanno --h5ad -c1 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/spliced_t2c.txt -c2 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/unspliced_t2c.txt -o kb_count_out_ht5 --overwrite \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD4/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD4/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R2_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD4/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD4/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R2_001.fastq.gz

#generate velocity count matrices - for 5th donor (HT6), BroadRun230330
kb count --verbose -i /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/index.idx -g /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/t2g.txt -x 10xv3 -m 30G --workflow lamanno --h5ad -c1 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/spliced_t2c.txt -c2 /n/groups/cbdm_lab/brh040/analysis/rna_velocity/kb_ref/unspliced_t2c.txt -o kb_count_out_ht6 --overwrite \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD5/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD5/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L001_R2_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD5/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R1_001.fastq.gz \
/n/groups/cbdm_lab/scRNAseq/BroadRun230330/230329_SL-NVV_0705_AH2FFCDRX3/BD/BD5/fastqs/outs/fastq_path/H2FFCDRX3/GEX_S1_L002_R2_001.fastq.gz
