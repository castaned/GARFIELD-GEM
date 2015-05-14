module load intel/compiler/13.0.1/64bit
module load intel/mkl/11.0u1/64bit
module load intel/mpi/4.1/64bit

export GARFIELD_HOME=/panfs/vol/GEM/RESTORED/SW/garfield_mpi
export HEED_DATABASE=$GARFIELD_HOME/Heed/heed++/database/
export ROOTSYS=/panfs/vol/GEM/RESTORED/SW/ROOT/root
source $ROOTSYS/bin/thisroot.sh 
