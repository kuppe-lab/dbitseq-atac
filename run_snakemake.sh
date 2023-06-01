#!/usr/local_rwth/bin/zsh
#SBATCH --account=rwth1209
#SBATCH --job-name=dbitseq
#SBATCH --chdir=/work/rwth1209/software/DbiT-seq-pipeline/Data_preprocessing
#SBATCH --ntasks=1 --cpus-per-task=20
#SBATCH --mem=64g
#SBATCH --time=12:00:00
#SBATCH --output=dbitseq%J.log

SLURM_ARGS="-p {cluster.partition} -J {cluster.job-name} -n {cluster.ntasks} -c {cluster.cpus-per-task} \
--mem={cluster.mem} -t {cluster.time} \
-o {cluster.output} -e {cluster.error}"

# Insert this after any #SLURM commands
#export CONDA_ROOT=$HOME/mambaforge
#. $CONDA_ROOT/etc/profile.d/conda.sh
#export PATH="$CONDA_ROOT/bin:$PATH"
# but naturally before using any python scripts
export CONDA_ROOT=/work/rwth1209/software/conda
. $CONDA_ROOT/etc/profile.d/conda.sh
export PATH="$CONDA_ROOT/bin:$PATH"


source activate /work/rwth1209/enviroments/dbitseq

snakemake -j 20 --cluster-config cluster.json --cluster "sbatch $SLURM_ARGS" 

