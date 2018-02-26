#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J _p31_
#SBATCH --clusters=mpp2
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=28
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=00:30:00

# modules:
source /etc/profile.d/modules.sh
module load gcc/4.8 

# libraries:
export MY_BASE=$HOME/GliomaAdvance/lib
export LD_LIBRARY_PATH=$MY_BASE/myVTK/lib/vtk-5.4/:$MY_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH

# Threads:
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# brain simulation set up
program=brain
vtk=1
dumpfreq=50
adaptive=1
verbose=1
IC=1
anatomy=HGG_UQ
pID=31

echo "In the directory: $PWD"
echo "Running program on $SLURM_NODES nodes, with $SLURM_CPUS_ON_NODE cores on node, each with $SLURM_CPUS_PER_TASK cores."

./$program -nthreads $SLURM_CPUS_ON_NODE -anatomy $anatomy -IC $IC -vtk $vtk -dumpfreq $dumpfreq -adaptive $adaptive -verbose $verbose -pID $pID 

