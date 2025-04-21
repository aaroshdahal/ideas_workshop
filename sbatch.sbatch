#!/bin/bash
#SBATCH -JFenics_Test                         # Job name
#SBATCH --account=gts-akumar355-paid            # charge account
#SBATCH -N1 --ntasks-per-node=5                # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=7G                        # Memory per core
#SBATCH -t00:15:00                               # Duration of the job (Ex: 1 hour)
#SBATCH -qinferno                                # QOS Name
##SBATCH -oReport-%j.out                         # Combined output and error messages file
##SBATCH --mail-type=BEGIN,END,FAIL             # Mail preferences
##SBATCH --mail-user=adahal8@gatech.edu       # E-mail address for notifications

cd $SLURM_SUBMIT_DIR                            # Change to working directory created in $HOME

# Load modules
module purge
module load fenics/2018-12-05
 

# Set paths
export PATH=/storage/coda1/p-akumar355/0/shared/hydra/bin:$PATH
export JOBID=`echo $SLURM_JOB_ID | cut -d"." -f1`
export LD_LIBRARY_PATH=/storage/home/hcodaman1/pace-apps-AMD/spack/packages/linux-rhel7-x86_64/gcc-4.8.5/intel-parallel-studio-cluster.2019.5-jnu223pfh6sdf4dhonmiz3iqxtqfacjp/vtune_amplifier_2019.6.0.602217/vpp/collector/hwloc/lib:$LD_LIBRARY_PATH

# Show loaded modules
module list
which mpirun

# Run Fenics MPI Code

cd /storage/home/hcoda1/9/adahal8/p-akumar355-0/glacier/ideas_workshop/
mpirun -np 5 apptainer exec $FENICS_SIF python3 surfing_fenics.py >& fenics.${JOBID}.oe
