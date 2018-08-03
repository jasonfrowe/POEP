program scenesimulator
!Code for generating FITS files given a PSF and point positions
use precision
use startype
implicit none
integer :: xout,yout,nunit,filestatus,nstarmax,nmax,npt,xmax,ymax, &
 nstars,i,j,noversample,iargc,nKs,ii,jj,seed,now(3)
real(double) :: flux,xcoo,ycoo,bpix,tavg,sigscale,sx,sy,sxy,ran2,dumr
real(double), allocatable, dimension(:,:) :: data,pixels,Kernel,cpixels,&
 opixels
character(80) :: prefix,fileout,fileout_c,fileout_m,modelfile,cnum
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
	subroutine displayfits(nxmax,nymax,parray,bpix,tavg,sigscale,seed)
     	use precision
        implicit none
        integer, intent(inout) :: nxmax,nymax,seed
        real(double), dimension(:,:), intent(inout) :: parray
        real(double), intent(inout) :: bpix,tavg
        real(double), intent(in) :: sigscale
    end subroutine displayfits
    subroutine genKernel(sx,sy,sxy,noversample,nKs,Kernel)
		use precision
		implicit none
		integer :: nKs,noversample
		real(double) :: sx,sy,sxy
		real(double), dimension(:,:) :: Kernel
	end subroutine genKernel
	subroutine convolveft(nKs,Kernel,xmax,ymax,pixels,cpixels)
		use precision
		implicit none
		integer :: nKs,xmax,ymax
		real(double), dimension(:,:) :: Kernel,pixels,cpixels
	end subroutine convolveft
	subroutine writefits(nxmax,nymax,parray,fileout,time)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax
      real(double), intent(inout) :: time
      real(double), dimension(:,:), intent(inout) :: parray
      character(80) :: fileout
   end subroutine
end interface

!Image Dimensions
xout=1024
yout=1024
!Image generation
noversample=4 !oversampling when stamping PSFs and convolving 
nKs=32*noversample !natural size of Kernels times oversampling
!input stars
nstarmax=472 !maximum number of stars
nmax=100000  !maximum number of data points in time-series 

!Name for output
prefix="sgen"

!Initialization of random number
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)

if(iargc().lt.1)then
   write(0,*) "Usage: scenesim <POEP_sim>"
   stop
else
	call getarg(1,modelfile)
endif

nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

allocate(data(4,nmax)) !store time and PSF
allocate(stars(nstarmax,nmax))
nstars=nstarmax
call readmodel(nunit,filestatus,nstarmax,nmax,npt,data,stars)
close(nunit)


!lets fill in pixel values
xmax=xout*noversample !size of over sampled detector array.
ymax=yout*noversample
!array to hold detector array values
allocate(pixels(xmax,ymax),cpixels(xmax,ymax),opixels(xout,yout))
!array to hold Kernel for convolution
allocate(Kernel(nKs,nKs))


bpix=99.9e30 !pixel with values above bpix are considered bad
sigscale=10.0 !clipping for display purposes, sig=0 displays full scale.
call pgopen('?')
call PGPAP (16.0 ,1.0) !use a square 8" across
call pgask(.false.)
do i=1,npt  !time-series
	pixels=0.0d0
	do j=1,nstars !stars for a single time.
		flux=stars(j,i)%flux*32768.0d0
		xcoo=stars(j,i)%xcoo*noversample
		ycoo=stars(j,i)%ycoo*noversample
		!write(0,*) xcoo,ycoo,flux
		call addflux2pix(xcoo,ycoo,xmax,ymax,pixels,flux)
	enddo
	tavg=data(1,i)
	!call displayfits(xmax,ymax,pixels,bpix,tavg,sigscale,seed)

	sx=data(2,i) !PSF parameters
	sy=data(3,i)
	sxy=data(4,i)
	call genKernel(sx,sy,sxy,noversample,nKs,Kernel) !generate Kernel
	!call displayfits(nKs,nKs,Kernel,bpix,tavg,sigscale,seed)
	!read(5,*)

	call convolveft(nKs,Kernel,xmax,ymax,pixels,cpixels)
	
	opixels=0.0d0 !reinitialize the array
	!dnossq=noversample*noversample
	do ii=noversample,xmax,noversample  !resample (bin) the array.
   		do jj=noversample,ymax,noversample
     		opixels(ii/noversample,jj/noversample)=                        &
          	  Sum(cpixels(ii-noversample+1:ii,jj-noversample+1:jj))
   		enddo
	enddo
	opixels=opixels
	call displayfits(xout,yout,opixels,bpix,tavg,sigscale,seed)

	write(0,*) "max, min: ",minval(opixels),maxval(opixels)
	!write(0,*) "Sum: ",Sum(opixels(960:1024,960:1024))
	

	write(cnum,'(I6.6)') i
	fileout=trim(prefix)//trim(cnum)//".fits"
!	fileout_m=trim(prefix)//"_m.txt"
!	fileout_c=trim(prefix)//"_c.fits"
	write(0,*) fileout
	call writefits(xout,yout,opixels,fileout,tavg)

	read(5,*)

enddo
call pgclos()

end program scenesimulator