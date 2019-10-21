#!/bin/sh

#  om_run_swjn_analysis_12.sh
#SBATCH --job-name=om_run_swjn_analysis_12
#SBATCH -t 1:00:00
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --array=1
#SBATCH --mem-per-cpu 5000
#SBATCH --exclude node017,node018
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ehoseini@mit.edu
#SBATCH --output=om_run_swjn_analysis_12_%j.out
#SBATCH --error=om_run_swjn_analysis_12_%j.err

module add mit/matlab/2018b
matlab -nodisplay -r "addpath(genpath('/home/ehoseini/MyCodes/ecog-sentence/')); \
swjn_om_analysis_12_align_anatomical_using_cpd;\
quit;"
