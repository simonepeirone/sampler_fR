#!/bin/env bash
##comment out lines by adding at least two `#' at the beginning
#SBATCH --job-name=fR_m2
#SBATCH --account=silvestri
#SBATCH --partition=computation
#SBATCH --output=/home/peirone/%x.out
#SBATCH --error=/home/peirone/%x.err
#SBATCH --time=100-00:00:00
#SBATCH --mem=10000
#SBATCH --ntasks=8
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=FAIL
# send mail to this address
#SBATCH --mail-user=peirone@lorentz.leidenuniv.nl

cd /marisdata/peirone/Stability_fR/

mpirun -np 8 ./eftcosmomc parameter_files/fR_mass2.ini > chains/fR_mass2.log
