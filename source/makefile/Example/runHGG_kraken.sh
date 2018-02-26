# Threads:
export OMP_NUM_THREADS=8

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

