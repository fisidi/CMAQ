program cmaq_driver

use pe_util

implicit none
integer :: i

call pe_init()

call init_logdev()

call pe_decomp(npes) 

if( mype .eq. 0) then

io_pe = .true.

end if

call cio()

call initscen()

do i = 1, nvars3d
  call pwrite3( CTM_CONC_1, vname3d(i), sdate3d, stime3d, cgrid(:,:,:,i) )
end do

call mpi_barrier(MPI_COMM_WORLD, ierr)
if (ierr /= 0) then
  call PM3EXIT('intscen', 0, 0,'MPI_barrier ERROR', xstat1)
end if



call pm3exit('driver',0,0,'Success!',0)



end program cmaq_driver
