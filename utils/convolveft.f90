subroutine convolveft(nKs,Kernel,xmax,ymax,pixels,cpixels)
use precision
!iso_c_binding is for FFTW3 interface
use, intrinsic :: iso_c_binding
implicit none
!add in the FFTW3 modules
include 'fftw3.f03'

!import vars
integer :: nKs,xmax,ymax
real(double), dimension(:,:) :: Kernel,pixels, cpixels

!FFTW3 vars
type(C_PTR) :: planA,planB,planC
integer ( kind = 4 ) :: nh
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: AC,BC,CC

!local vars
integer :: XF,YF
real(double), allocatable, dimension(:,:) :: A,B,C

cpixels=0.0d0! initalize convolved image to zero.

XF=nKs+xmax !size of zero padded arrays for convolution
YF=nKs+ymax
allocate(A(XF,YF),B(XF,YF),C(XF,YF)) !arrays to apply FFT
nh=(XF/2)+1 !size for complex array
allocate(AC(nh,YF),BC(nh,YF),CC(nh,YF)) !allocate complex arrays for FT

!only need to FFT pixels once..
A=0.0d0 !initalize to zero
A(1:xmax,1:ymax)=pixels(1:xmax,1:ymax) !assign Image to zero-padded A
planA=fftw_plan_dft_r2c_2d(YF,XF,A,AC,FFTW_ESTIMATE)
call fftw_execute_dft_r2c(planA,A,AC)

call fftw_destroy_plan(planA)

!We can precompute plans..
planB=fftw_plan_dft_r2c_2d(YF,XF,B,BC,FFTW_ESTIMATE)
planC=fftw_plan_dft_c2r_2d(YF,XF,CC,C,FFTW_ESTIMATE)

B=0.0d0
C=0.0d0
!Kernel should be centered around 0,0 on B
B(1:nKs/2      ,1:nKs/2)      =Kernel(nKs/2:nKs,nKs/2:nKs)
B(XF-nKs/2+1:XF,YF-nKs/2+1:YF)=Kernel(1:nKs/2  ,1:nKs/2)
B(XF-nKs/2+1:XF,1:nKs/2)      =Kernel(1:nKs/2  ,nKs/2:nKs)
B(1:nKs/2,      YF-nKs/2+1:YF)=Kernel(nKs/2:nKs,1:nKs/2)


call fftw_execute_dft_r2c(planB,B,BC)

CC=AC*BC
call fftw_execute_dft_c2r(planC,CC,C)
cpixels=C(1:xmax,1:ymax)/dble(xmax*ymax)

call fftw_destroy_plan(planB)
call fftw_destroy_plan(planC)

return
end

