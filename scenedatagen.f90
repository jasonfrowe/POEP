program scenedatagenerator
!generate file that is used by the scene simulator
use precision
implicit none
integer :: nunit,filestatus,npt,i,nmax,nf
real(double) :: mindate
real(double), allocatable, dimension(:,:) :: data
character(80) :: modelfile,cline

!! Code Parameters !!
modelfile="MOST_model_20180626.dat"
nmax=400000 !maximum number of data points 
nf=16 !number of fields in MOST data

nunit=10 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

!skip first two lines of file (header)
read(nunit,*) cline
read(nunit,*) cline

!allocate data array to read in MOST file
allocate(data(nf,nmax))

npt=1   !initalize counter for number of data points
do !we do a loop.  If there are memory errors, we can get more this way
	read(nunit,*,iostat=filestatus) (data(i,npt),i=1,nf)	
	if(filestatus == 0) then
		npt=npt+1
	elseif(filestatus == -1) then
      	exit  !successively break from data read loop.
    else
      	write(0,*) "File Error!! Line:",npt
      	write(0,900) "iostat: ",filestatus
      	900 format(A8,I3)
      	stop
    endif
enddo
close(nunit) !close file.
npt=npt-1
write(0,*) "Number of model points: ",npt  !report number of data points read.


mindate=minval(data(1,1:npt))

!writing out datafile
write(6,*) "Date xsig  xsig  xysig"
do i=1,10
	write(6,500) data(1,i)-mindate,data(7,i),data(8,i),data(9,i)
enddo
500 format(F13.8,1X,2(F5.3,1X),F6.3,1X)

end program scenedatagenerator