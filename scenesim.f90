program scenesimulator
!Code for generating FITS files given a PSF and point positions
use precision
use startype
implicit none
integer :: xout,yout,nunit,filestatus,nstarmax,nmax,npt,xmax,ymax, &
 nstars,i,j,noversample
real(double) :: flux,xcoo,ycoo
real(double), allocatable, dimension(:,:) :: data,pixels
character(80) :: prefix,fileout,fileout_c,fileout_m,modelfile
type(starpos), allocatable, dimension(:,:) :: stars

interface
	subroutine readmodel(nunit,filestatus,nstarmax,nmax,npt,data,stars)
		use precision
		use startype
		implicit none
		integer :: nunit,filestatus,nstarmax,nmax,npt
		real(double), dimension(:,:) :: data
		type(starpos), dimension(:,:) :: stars
	end subroutine readmodel
end interface

!Image Dimensions
xout=1024
yout=1024
!Image generation
noversample=4 !oversampling when stamping PSFs and convolving 
!input stars
nstarmax=472 !maximum number of stars
nmax=100000  !maximum number of data points in time-series 

!Name for output
prefix="sgen"
fileout=trim(prefix)//".fits"
fileout_m=trim(prefix)//"_m.txt"
fileout_c=trim(prefix)//"_c.fits"

!input file
modelfile="POEP_sim_20180802.dat"

nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

allocate(data(4,nmax)) !store time and PSF
allocate(stars(nstars,nmax))
nstars=nstarmax
call readmodel(nunit,filestatus,nstarmax,nmax,npt,data,stars)
close(nunit)


!lets fill in pixel values
xmax=xout*noversample !size of over sampled detector array.
ymax=yout*noversample
!array to hold detector array values
allocate(pixels(xmax,ymax))

do i=1,npt  !time-series
	do j=1,nstars !stars for a single time.
		flux=stars(j,i)%flux
		xcoo=stars(j,i)%xcoo
		ycoo=stars(j,i)%ycoo
		write(0,*) xcoo,ycoo,flux
		call addflux2pix(xcoo,ycoo,xmax,ymax,pixels,flux)
	enddo


	read(5,*)
enddo


end program scenesimulator