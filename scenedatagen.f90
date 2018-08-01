program scenedatagenerator
!generate file that is used by the scene simulator
use precision
use startype
implicit none
integer :: nunit,filestatus,npt,i,j,k,nmax,nf,xout,yout,nsub,nstars,seed, &
 ii,jj,kk
integer, dimension(3) :: now
real(double) :: mindate,bpix,tavg,sigscale,ran2,dumr
real(double), allocatable, dimension(:,:) :: data,pixels
character(80) :: modelfile,cline
type(starpos), allocatable, dimension(:,:) :: stars,stars2

interface
	subroutine genstarpos(npt,nf,nstars,data,stars,stars2,xout,yout,nsub,seed)
		use precision
		use startype
		implicit none
		integer :: npt,nf,xout,yout,nsub,nstars,seed
		real(double), dimension(:,:) :: data
		type(starpos), dimension(:,:) :: stars,stars2
	end subroutine genstarpos
    subroutine displayfits(nxmax,nymax,parray,bpix,tavg,sigscale)
     	use precision
        implicit none
        integer, intent(inout) :: nxmax,nymax
        real(double), dimension(:,:), intent(inout) :: parray
        real(double), intent(inout) :: bpix,tavg
        real(double), intent(in) :: sigscale
    end subroutine displayfits
end interface

!! Code Parameters !!
modelfile="MOST_model_20180626.dat"
nmax=400000 !maximum number of data points 
nf=16 !number of fields in MOST data
!image dimensions
xout=1024  !dimensions for output image.
yout=1024
nsub=16 !number of subrasters for star sims
nstars=nsub*nsub !number of stars to generate.

!Initialization of random number
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)

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


!allocate space to contain properties of simulated stars
!stars - primary target
!stars2 - secondary target
allocate(stars(nstars,npt),stars2(nstars,npt))

!generate positions 
call genstarpos(npt,nf,nstars,data,stars,stars2,xout,yout,nsub,seed)

!generate fluxes
call genstarflux(npt,nstars,stars,stars2)
!stars(:,:)%flux=1.0d0 !give a constant for the moment.  
!stars2(:,:)%flux=0.2d0

!place stars on pixel map

allocate(pixels(xout,yout))
!display fits file
call pgopen('?')
!call pgopen('/xserve')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgask(.false.)
do kk=1,npt
	pixels=0.0
	do i=1,nstars
		do j=-2,2
			do k=-2,2
				ii=int(stars(i,kk)%xcoo)+j
				jj=int(stars(i,kk)%ycoo)+k
				pixels(ii,jj)=pixels(ii,jj)+stars(i,kk)%flux
				ii=int(stars2(i,kk)%xcoo)+j
				jj=int(stars2(i,kk)%ycoo)+k
				pixels(ii,jj)=pixels(ii,jj)+stars2(i,kk)%flux
			enddo
		enddo
	enddo

	bpix=99.9e30
	tavg=0.0d0
	sigscale=3.0
	call displayfits(xout,yout,pixels,bpix,tavg,sigscale)
enddo
call pgclos()

!calculate min date so that the time-series starts at zero.
mindate=minval(data(1,1:npt))

!writing out datafile
write(6,*) "Date xsig  xsig  xysig"
do i=1,10
	write(6,500) data(1,i)-mindate,data(7,i),data(8,i),data(9,i)
enddo
500 format(F13.8,1X,2(F5.3,1X),F6.3,1X)

end program scenedatagenerator