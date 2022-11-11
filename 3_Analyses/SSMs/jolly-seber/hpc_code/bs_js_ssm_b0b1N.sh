#!/bin/bash
#SBATCH --job-name=js_ssm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=qbiol
#SBATCH --mem=8GB                                                         # memory (MB)
#SBATCH --time=0-01:00                                                    # time (D-HH:MM)
#SBATCH -o /home/s4532001/7008_guppies/output/oe_files/js_ssm_b0b1N.%j.%a.out    # STDOUT
#SBATCH -e /home/s4532001/7008_guppies/output/oe_files/js_ssm_b0b1N.%j.%a.err    # STDERR
#SBATCH --array=1-4

echo "Start time: "; date

module load applications/R/4.1.3

Rscript --vanilla /home/s4532001/7008_guppies/scripts/js_ssm_b0b1N.R $SLURM_ARRAY_TASK_ID

echo "End time: "; date