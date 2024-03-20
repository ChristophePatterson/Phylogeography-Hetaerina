#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=10G
#SBATCH --gres=tmp:10G
#SBATCH -t 72:00:00
#SBATCH --output=slurm-%x.%j.out

/nobackup/tmjj24/apps/G-PhoCS/bin/G-PhoCS /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/Hetaerina_americana_ddRAD_americana_dg/G-Phocs/model_runs/G-Phocs-a1-b20-americana-Iso-run0.config -n 24 > /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/Hetaerina_americana_ddRAD_americana_dg/G-Phocs/model_runs/G-Phocs-a1-b20-americana-Iso-run0.log
