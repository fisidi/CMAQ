module pe_util

use mpi
use m3utilio
implicit none

!*******************************!
!         Domain Decomp.        !     
integer, parameter :: nc=2, nr=1!
!*******************************!

integer :: mype_ncols, mype_nrows, mype_nlays 
integer :: mype_sc, mype_ec, mype_sr, mype_er 
integer :: ierr, mype, npes
integer :: logdev = -1
logical :: io_pe
real(8) :: palpha, pbeta, pgama, xcent, ycent, xorig, yorig, xcell, ycell
integer :: ncols, nrows, nthik, nlays, gdtyp
integer,allocatable :: ncols_pe( : ), nrows_pe( : )
integer,allocatable :: ncols_pe_ndx(:, : ), nrows_pe_ndx(:, : )

real, allocatable :: cgrid(:,:,:,:)

character(40), parameter :: CTM_CONC_1 = 'CTM_CONC_1' 

contains 

subroutine pe_init()
implicit none

call MPI_INIT( ierr )

if (ierr /= 0) then
  write(*, '(5x,A)') 'MPI_INIT ERROR'
  call exit(ierr)
end if

call MPI_COMM_SIZE( MPI_COMM_WORLD, npes, ierr )
if (ierr /= 0) then
  write(*, '(5x,A)') 'MPI_COM_SIZE ERROR'
  call exit(ierr)
end if

call MPI_COMM_RANK( MPI_COMM_WORLD, mype, ierr )
if (ierr /= 0) then
  write(*, '(5x,A)') 'MPI_COM_SIZE ERROR'
  call exit(ierr)
end if

if ( mype == 0) then
  io_pe = .true. 
else 
  io_pe = .false.
end if

end subroutine pe_init


subroutine pm3exit(caller, jdate, jtime, msgtxt, exitstat) 

implicit none
!...........   ARGUMENTS and their descriptions:

CHARACTER(*),intent(in) ::  CALLER       ! Name of the caller.
INTEGER     ,intent(in) ::  JDATE        ! Model date for the error.
INTEGER     ,intent(in) ::  JTIME        ! Model time for the error.
CHARACTER(*),intent(in) ::  MSGTXT       ! Error message.
INTEGER     ,intent(in) ::  EXITSTAT     ! Exit status for program.


!...........   LOCAL VARIABLES

INTEGER      LENSTR       ! String length of CALLER.
INTEGER      IERROR       ! Error from MPI_ABORT routine.
CHARACTER(7) PE_STR       ! String suffix to go with processor ID.
CHARACTER(16) CALL16      ! First 16 characters of CALLER.
CHARACTER(26) PCALLER     ! New caller string with PE information.
CHARACTER(3) CMYPE

WRITE (PE_STR,'(A7)') ' on PE '
WRITE(CMYPE, '(I3.3)') MYPE
!.......  Construct new CALLER string.
LENSTR = MIN( 16, LEN_TRIM(CALLER) )
CALL16 = CALLER( 1: LENSTR )
PCALLER = CALL16(1:LENSTR)//PE_STR//trim(CMYPE)






CALL MPI_FINALIZE( IERROR )

if (ierror /= 0) then
  write(*, '(5x,A)') 'MPI_FINALIZE ERROR'
end if

call m3exit(pcaller, jdate, jtime,msgtxt,exitstat) 


end subroutine pm3exit


subroutine init_logdev

implicit none
CHARACTER(  8 ) :: PRESTR  = 'CTM_LOG_'
CHARACTER( 200 ) :: IOLOGEQ
CHARACTER( 90 ) :: CMAQ_HEADER( 21 )
CHARACTER(3) CMYPE
integer :: I, IOST, NIOAPI_HEAD, IHEAD
integer,parameter :: NHEAD=200
CHARACTER( 90 ) :: IOAPI_HEADER( NHEAD )

CMAQ_HEADER(1:21) = (/  & 
   '================================================================================', &
   '|                                                                              |', &
   '|               The Community Multiscale Air Quality (CMAQ) Model              |', &
   '|                                   Version 5.3.2                              |', &
   '|                                                                              |', &
   '|                          Built and Maintained by the                         |', &
   '|                        Office of Research and Development                    |', &
   '|                   United States Environmental Protection Agency              |', &
   '|                                                                              |', &
   '|                            https://www.epa.gov/cmaq                          |', &
   '|                                                                              |', &
   '|       Source Code:   https://www.github.com/USEPA/cmaq/tree/master           |', &
   '|       Documentation: https://www.github.com/USEPA/cmaq/tree/master/DOCS      |', &
   '|                                                                              |', &
   '|         The CMAQ Model is tested and released with cooperation from          |', &
   '|         the Community Modeling and Analysis System (CMAS) Center via         |', &
   '|         contract support. CMAS is managed by the Institute for the           |', &
   '|         Environment, University of North Carolina at Chapel Hill.            |', &
   '|         CMAS URL: (https://www.cmascenter.org)                               |', &
   '|                                                                              |', &
   '================================================================================'  &
      /)                                                                               

WRITE(CMYPE, '(I3.3)') MYPE
IOLOGEQ = PRESTR // CMYPE
IF ( .NOT. SETENVVAR ( 'LOGFILE', IOLOGEQ ) ) THEN
   WRITE( *,* ) '*** Could not set environment variable for ' // IOLOGEQ
   CALL PM3EXIT( 'SETUP_LOGDEV', 0, 0, '', 2 )
END IF

IF (MYPE == 0) then
DO I = 1,21

WRITE(6, '(A)') CMAQ_HEADER(I)

END DO

END IF
! Redirect Standard Output
IF ( MYPE .NE. 0 ) OPEN( UNIT = 6, FILE = "/dev/null", STATUS = "OLD" )

LOGDEV = INIT3() 

! Put Standard Output Back to File Unit 6 By Closing File Unit 6
IF ( MYPE .NE. 0 ) CLOSE( 6 ) 


REWIND( LOGDEV )
DO IHEAD = 1,NHEAD
    READ( LOGDEV, '(A)', IOSTAT=IOST ), IOAPI_HEADER( IHEAD )
    ! Check for End of File
    IF ( IOST .LT. 0 ) THEN
       NIOAPI_HEAD = IHEAD - 1
       EXIT
    END IF
END DO
REWIND( LOGDEV )

DO I = 1,21

WRITE(LOGDEV, '(A)') CMAQ_HEADER(I)

END DO

WRITE(LOGDEV, *)   

DO IHEAD = 1,NIOAPI_HEAD
    WRITE( LOGDEV, '(A)' ) IOAPI_HEADER( IHEAD )
END DO

end subroutine init_logdev

subroutine pe_decomp(nprocs)

implicit none
CHARACTER(16) :: GRID_NAME = '2016_12SE1'       ! grid name selected from GRIDDESC
CHARACTER( 16 ) :: COORD_SYS_NAME
CHARACTER(32) :: XMSG
integer, intent(in) :: nprocs
integer :: i, ndx
character( 49 ) :: colrow = '  PE    #Cols    Col_Range     #Rows    Row_Range'
character( 49 ) :: title



allocate( ncols_pe(nprocs), nrows_pe(nprocs), ncols_pe_ndx(3,nprocs), nrows_pe_ndx(3,nprocs), stat=ierr) 

if (ierr .ne. 0 ) then
  xmsg = ' Failure to allocate ncols_pe'
  call PM3EXIT('pe_decomp', 0, 0, xmsg, xstat1)
end if
 


IF ( .NOT. DSCGRID ( GRID_NAME, &
                     COORD_SYS_NAME, GDTYP, & 
                     PALPHA, PBETA, PGAMA, &
                     XCENT, YCENT, &
                     XORIG, YORIG, XCELL, YCELL, &
                     NCOLS, NROWS, NTHIK ) ) THEN
   XMSG = 'Failure retrieving horizontal grid parameters'
   CALL PM3EXIT ( 'HGRD_INIT', 0, 0, XMSG, XSTAT1 )
   
END IF

ncols_pe = ncols/nc
nrows_pe = nrows/nr 

do i=1,ncols-ncols_pe(1)*nc
  ncols_pe(i) = ncols_pe(i) + 1
end do

do i=1,nrows-nrows_pe(1)*nr
  nrows_pe(i) = nrows_pe(i) + 1
end do

do i=1, nprocs

ndx = mod(i,nc)
if (ndx .eq. 0) ndx=nc

ncols_pe_ndx(1,i) = ncols_pe(ndx)

if ( ndx .eq. 1) then

ncols_pe_ndx(2,i) = 1
ncols_pe_ndx(3,i) = ncols_pe_ndx(1,i) 

else 

ncols_pe_ndx(2,i) = ncols_pe_ndx(3,i-1)+1
ncols_pe_ndx(3,i) = ncols_pe_ndx(3,i-1)+ncols_pe_ndx(1,i)

end if


ndx = (i-1)/nc + 1 

nrows_pe_ndx(1,i) = nrows_pe(ndx)

if ( i .le. nc) then

nrows_pe_ndx(2,i) = 1
nrows_pe_ndx(3,i) = nrows_pe_ndx(1,i) 

else 

nrows_pe_ndx(2,i) = nrows_pe_ndx(3,i-nc)+1
nrows_pe_ndx(3,i) = nrows_pe_ndx(3,i-nc)+nrows_pe_ndx(1,i)

end if






end do 


mype_ncols = ncols_pe_ndx(1,mype+1)
mype_sc = ncols_pe_ndx(2,mype+1) 
mype_ec = ncols_pe_ndx(3,mype+1) 

mype_nrows = nrows_pe_ndx(1,mype+1)
mype_sr = nrows_pe_ndx(2, mype+1)
mype_er = nrows_pe_ndx(3,mype+1)

mype_nlays = 35

title = colrow


if (mype .eq. 0) then
 write( *,* )
 write( *,* ) '         -=-  MPP Processor-to-Subdomain Map  -=-'
 write( *,'(A,I3)' ) '                 Number of Processors = ',nprocs
 write( *,* ) '   ____________________________________________________'
 write( *,* ) '   |                                                  |'
 write( *,* ) '   |' // title // ' |'
 write( *,* ) '   |__________________________________________________|'
 write( *,* ) '   |                                                  |'
 do i = 1, nprocs
    write( *,1003 ) i-1, ncols_pe_ndx(1,i), ncols_pe_ndx(2,i), ncols_pe_ndx(3,i), &
                         nrows_pe_ndx(1,i), nrows_pe_ndx(2,i), nrows_pe_ndx(3,i)
 end do
 write( *,* ) '   |__________________________________________________|'
 write( *,* )

end if



1003  format('    |', i3, 5x, i4, 3x, i4, ':', i4, &
                         7x, i4, 3x, i4, ':', i4, '   |')

end subroutine pe_decomp



subroutine cio()

implicit none
integer :: ierr, i
character(128) :: xmsg
CHARACTER( 40 ), parameter :: ICFILE       = 'INIT_CONC_1'

if (.not. OPEN3(ICFILE, fsread3, 'cio')) then
  XMSG = 'Could not get description of file  '// ICFILE
  CALL PM3EXIT ( 'cio', 0, 0, XMSG, XSTAT1 )
end if

IF ( .NOT. DESC3( ICFILE ) ) THEN
   XMSG = 'Could not get description of file  '// ICFILE 
   CALL PM3EXIT ( 'cio', 0, 0, XMSG, XSTAT1 )
END IF


allocate(cgrid(mype_ncols,mype_nrows,mype_nlays,nvars3d), stat=ierr) 

if (ierr .ne. 0 ) then
   XMSG = 'Could not get allocate CGRID array'
   CALL PM3EXIT ( 'cio', 0, 0, XMSG, XSTAT1 )
end if


do i=1,nvars3d

IF ( .NOT. XTRACT3( ICFILE, vname3d(i), &
                    1, mype_nlays, mype_sr,mype_er, mype_sc, mype_ec, &
                    sdate3d, stime3d, cgrid(:,:,:,i) ) ) THEN
   XMSG = 'Could not extract ' // ICFILE // ' file'
   CALL PM3EXIT ( 'cio', 0, 0, XMSG, XSTAT1 )
END IF

end do 


end subroutine cio

subroutine initscen()

implicit none

character(128) :: MSG
LOGICAL, EXTERNAL :: FLUSH3

IF ( .NOT. OPEN3( CTM_CONC_1, FSRDWR3, 'initscen' ) ) THEN
   BACKSPACE( LOGDEV )
   MSG = 'Could not open ' // TRIM( CTM_CONC_1 ) // ' for update - ' &
        // 'try to open new'
   CALL M3MESG( MSG )
   ! Open the file new
   call mpi_barrier(MPI_COMM_WORLD, ierr)
   if (ierr /= 0) then
     call PM3EXIT('intscen', 0, 0,'MPI_barrier ERROR', xstat1)
   end if

IF ( IO_PE ) THEN   ! open new

   IF ( .NOT. OPEN3( CTM_CONC_1, FSNEW3, 'initscen' ) ) THEN
      MSG = 'Could not open ' // CTM_CONC_1
      CALL PM3EXIT( 'initscen', 0, 0, MSG, XSTAT1 )
   END IF

   IF ( .NOT. FLUSH3 ( CTM_CONC_1 ) ) THEN
      MSG = 'Could not sync to disk ' // CTM_CONC_1
      CALL PM3EXIT( 'initscen', 0,0, MSG, XSTAT1 )
   END IF

END IF

END IF


end subroutine initscen


SUBROUTINE LOG_HEADING( FUNIT, CHEAD_IN )

!  Formats and writes a user-supplied heading to a specific log file.
!  This approach is intended to standardize the log files that are
!  created by CMAQ. The length of the input array is set at 80 because
!  we would like to try limiting lines to 80 characters and a heading
!  should probably just be one line.
!.........................................................................

IMPLICIT NONE

INTEGER, INTENT( IN )           :: FUNIT
CHARACTER( * ), INTENT( IN )    :: CHEAD_IN
CHARACTER( len=: ), ALLOCATABLE :: CHEAD
CHARACTER( 20 )                 :: FMT
CHARACTER( 20 )                 :: FMT2
INTEGER                         :: LDASH
INTEGER :: LOG_MAJOR_TAB   = 5   ! Left tab for all text including headings

! Capitalize the heading
CHEAD = CHEAD_IN
CALL UPCASE( CHEAD )

! Write the heading to the log file
WRITE( FUNIT, * )
WRITE( FMT, '("(", I0, "x,A,A,A)")' ) LOG_MAJOR_TAB
WRITE( FMT2,'("(", I0, "x,A,)")' ) LOG_MAJOR_TAB

LDASH = 2*8 + LEN_TRIM( CHEAD )
WRITE( FUNIT, FMT2 ), REPEAT( '=', LDASH )
WRITE( FUNIT, FMT ), &
'|>---   ',TRIM( CHEAD ),'   ---<|'
WRITE( FUNIT, FMT2 ), REPEAT( '=', LDASH )

END SUBROUTINE LOG_HEADING

subroutine pwrite3(fname, vname, jdate, jtime, buffer)
implicit none

character(*), intent(in) :: fname, vname
integer, intent(in) :: jdate, jtime
real, intent(in) :: buffer(:,:,:)


integer :: wsize, rsize, ierror, i,j,k,m, loc
INTEGER, SAVE :: RECVBUF_SIZE = 0
INTEGER, SAVE :: WRITBUF_SIZE = 0
real, allocatable :: WRITBUF(:,:,:)
real, allocatable :: RECVBUF( : )
character(128) :: MSG
INTEGER        MSGSIZE       ! Message size of subgrid to receive
INTEGER        WHO           ! For identifying sending processor
INTEGER        C0            ! First column in global grid
INTEGER        R0            ! First row in global grid
INTEGER        CX
INTEGER        RX
INTEGER        NC            ! Number of columns in local grid
INTEGER        NR            ! Number of rows in local grid
INTEGER        STATUS( MPI_STATUS_SIZE )   ! MPI status code
INTEGER, PARAMETER :: TAG1 = 901 ! MPI message tag for processor ID
INTEGER, PARAMETER :: TAG2 = 902 ! MPI message tag for data array.




   wsize = ncols3d*nrows3d*nlays3d
   rsize = mype_ncols*mype_nrows*nlays3d

    IF ( ALLOCATED ( WRITBUF ) ) DEALLOCATE ( WRITBUF )
    ALLOCATE ( WRITBUF  ( NCOLS3D, NROWS3D, NLAYS3D ), STAT = IERROR )
    IF ( IERROR .NE. 0 ) THEN
       MSG = 'Failure allocating WRITBUF '
       CALL PM3EXIT( 'PWRGRDD', JDATE, JTIME, MSG, 1 )
    END IF
    WRITBUF = 0.0


    IF ( ALLOCATED ( RECVBUF ) ) DEALLOCATE ( RECVBUF )
    ALLOCATE ( RECVBUF  ( RSIZE ), STAT = IERROR )
    IF ( IERROR .NE. 0 ) THEN
       MSG = 'Failure allocating RECVBUF '
       CALL PM3EXIT( 'PWRGRDD', JDATE, JTIME, MSG, 1 )
    END IF
    RECVBUF = 0.0

if ( io_pe ) then 

  WRITBUF(mype_sc:mype_ec,mype_sr:mype_er,:) = buffer

  do i = 1,npes-1
   CALL MPI_RECV( WHO, 1, MPI_INTEGER, MPI_ANY_SOURCE, TAG1, MPI_COMM_WORLD, STATUS, IERROR )


   IF ( IERROR .NE. 0 ) THEN
      MSG = 'MPI error receiving processor id WHO.'
      CALL PM3EXIT( 'PWRGRDD', JDATE, JTIME, MSG, XSTAT1 )
   END IF

   NC = ncols_pe_ndx(1,WHO+1)
   NR = nrows_pe_ndx(1,WHO+1)
   C0 = ncols_pe_ndx(2,WHO+1)
   CX = ncols_pe_ndx(3,WHO+1)
   R0 = nrows_pe_ndx(2,WHO+1)
   RX = nrows_pe_ndx(3,WHO+1)

   MSGSIZE = NC * NR * NLAYS3D
   CALL MPI_RECV( RECVBUF, MSGSIZE, MPI_REAL, WHO, TAG2, MPI_COMM_WORLD, STATUS, IERROR )


   IF ( IERROR .NE. 0 ) THEN
      MSG = 'MPI error receiving data array RECVBUF.'
      CALL PM3EXIT( 'PWRGRDD', JDATE, JTIME, MSG, XSTAT1 )
   END IF

   loc = 0
   do k = 1, nlays3d
     do j = R0,RX
       do m = C0, CX
         loc = loc+1  
         WRITBUF(m,j,k) = RECVBUF(loc)
        end do
      end do
    end do

  end do

  IF ( .NOT. WRITE3( fname, vname, jdate, jtime, WRITBUF ) ) THEN
     MSG = 'Could not write ' &
     // TRIM( VNAME ) // &
     ' to file '// TRIM( FNAME )
     CALL PM3EXIT( 'PWRGRDD', JDATE, JTIME, MSG, XSTAT1)
  END IF

else 

  WHO = mype
  MSGSIZE = rsize 
  CALL MPI_SEND( WHO, 1, MPI_INTEGER, IO_PE, TAG1, MPI_COMM_WORLD, IERROR )

  IF ( IERROR .NE. 0 ) THEN
       MSG = 'MPI error sending processor id WHO.'
       CALL PM3EXIT( 'PWRGRDD', JDATE, JTIME, MSG, XSTAT1 )
  END IF
  
  CALL MPI_SEND( BUFFER, MSGSIZE, MPI_REAL, IO_PE, TAG2, MPI_COMM_WORLD, IERROR )

  IF ( IERROR .NE. 0 ) THEN
     MSG = 'MPI error sending data array BUFFER.'
     CALL PM3EXIT( 'PWRGRDD', JDATE, JTIME, MSG, XSTAT1 )
  END IF


end if


end subroutine pwrite3

end module pe_util
