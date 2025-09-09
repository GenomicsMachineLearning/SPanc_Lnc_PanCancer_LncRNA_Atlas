#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=250GB
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --job-name=rnasamba
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --account=a_nguyen_quan

#module load miniconda3/4.12.0
source ~/.bashrc
#source activate hlpiensemble
source activate rnasamba

#rnasamba classify all_cuTARs.tsv /QRISdata/Q4386/polyA_FreshFrozen/uTAR_fasta/all_cuTARs.fa partial_length_weights.hdf5

#for i in /QRISdata/Q4386/polyA_FreshFrozen/uTAR_fasta/*_uTARs.fa; do rnasamba classify $i.rnasamba.txt $i partial_length_weights.hdf5 ; done

#rnasamba classify all_cuTARs_rnasamba.txt ../all_cuTARs_refined.fa partial_length_weights.hdf5

rnasamba classify all_cuTARs_v5_rnasamba.txt ../all_cuTARs_refined_v5.fa partial_length_weights.hdf5


### CPAT
#extract fasta for coding potential 
bedtools getfasta -fi /nfs_master/prakrithi/test/ST/refdata-gex-GRCh38-2020-A/fasta/genome.fa -bed all_cuTARs.bed -s > all_cuTARs.fa

# CPAT - Coding potential
nohup ~/miniconda3/envs/cpat_new/bin/cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData -g all_cuTARs.fa -o all_cuTARs_cpat 
