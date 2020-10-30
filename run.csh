#!/bin/csh -f



#SBATCH --partition=compute 
#SBATCH --ntasks=64

setenv GRIDDESC /work/MOD3DATA/SE53BENCH/GRIDDESC
setenv INIT_CONC_1 /work/MOD3DATA/SE53BENCH/icbc/ICON_20160630_bench.nc
setenv CTM_CONC_1 "/work/MOD3DEV/fsidi/trash/ioapi_pio/CCTM_CONC.nc -v"
setenv IOAPI_LOG_WRITE F     #> turn on excess WRITE3 logging [ options: T | F ]


rm -rf CTM_* CCTM_CONC.nc 

cd /work/MOD3DEV/fsidi/trash/ioapi_pio
mpirun -np 2 ./test1.exe

if ($? == 0) then
 echo "     --->> Normal Completion of program driver"
 echo "     Success!"
else
 echo "     --->> ERROR: Program was abort with exit status: "
 echo $?
endif 
