subroutine readmodel(nunit,filestatus,nstarmax,nmax,npt,data,stars)
use precision
use startype
implicit none
!input vars
integer :: nunit,filestatus,nstarmax,nmax,npt
real(double), dimension(:,:) :: data
type(starpos), dimension(:,:) :: stars
!local vars
integer :: i



npt=1   !initalize counter for number of data points
do !we do a loop.  If there are memory errors, we can get more this way
	read(nunit,*,iostat=filestatus) (data(i,npt),i=1,4),(stars(i,npt)%xcoo,stars(i,npt)%ycoo,& 
	  stars(i,npt)%flux,i=1,nstarmax)
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


return
end subroutine readmodel