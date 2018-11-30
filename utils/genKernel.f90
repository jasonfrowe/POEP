subroutine genKernel(sx,sy,sxy,noversample,nKs,Kernel)
use precision
implicit none
!import vars
integer :: nKs,noversample
real(double) :: sx,sy,sxy
real(double), dimension(:,:) :: Kernel
!local vars
integer i,j
real(double) :: xc,yc,di,dj,xd,yd,expt1,expt2,expt3

xc=dble(nKs)/2+0.5 !center of PSF
yc=dble(nKs)/2+0.5 

do i=1,nKs
	di=dble(i)
	xd=(di-xc)/noversample !distance from centre of PSF
	do j=1,nKs
		dj=dble(j)
		yd=(dj-yc)/noversample
		expt1=-((xd/sx)**2.0)
		expt2=-((yd/sy)**2.0)
		expt3=2.0*sxy*(xd*yd)/(sx*sy)
		Kernel(i,j)=exp(expt1+expt2+expt3)
        !write(0,*) i,j,Kernel(i,j)
	enddo
	!read(5,*)
enddo
Kernel=Kernel/Sum(Kernel) !normalize

return
end