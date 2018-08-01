subroutine genstarpos(npt,nf,nstars,data,stars,stars2,xout,yout,nsub,seed)
use precision
use startype
implicit none
!import vars
integer :: npt,nf,xout,yout,nsub,nstars,seed
real(double), dimension(:,:) :: data
type(starpos), dimension(:,:) :: stars,stars2
!local vars
integer :: i,j,k
real(double) :: xcoo_med,ycoo_med,median,dx,dy,xcoo_p,ycoo_p,r,dr,dnsub,th, &
 Pi,tPi,ran2
real(double), allocatable, dimension(:) :: xcoo,ycoo,xcoo2,ycoo2


!Constants
Pi=acos(-1.d0)!define Pi and 2*Pi
tPi=2.0d0*Pi

!estimate median of x,y centroids
xcoo_med=median(npt,data(5,1:npt))
ycoo_med=median(npt,data(6,1:npt))

!spacing for primary stars
dx=dble(xout)/dble(nsub)
dy=dble(yout)/dble(nsub)

!positional info for secondary
r=min(dx,dy)/2.0*0.96

!pre-calculate doubles
dnsub=dble(nsub)

allocate(xcoo(nstars),ycoo(nstars),xcoo2(nstars),ycoo2(nstars))
k=0
do i=1,nsub
	xcoo_p=dx*(i-1)+dx/2.0d0
	dr=r*dble(i-1)/dnsub
	do j=1,nsub
		k=k+1
		ycoo_p=dy*(j-1)+dy/2.0d0
		xcoo(k)=xcoo_p
		ycoo(k)=ycoo_p
		th=tPi*ran2(seed)
		xcoo2(k)=xcoo_p+dr*cos(th)
		ycoo2(k)=ycoo_p+dr*sin(th)
	enddo
enddo

do j=1,npt
	do i=1,nstars
		stars(i,j)%xcoo = xcoo(i)+data(5,j)-xcoo_med
		stars(i,j)%ycoo = ycoo(i)+data(6,j)-ycoo_med

		stars2(i,j)%xcoo = xcoo2(i)+data(5,j)-xcoo_med
		stars2(i,j)%ycoo = ycoo2(i)+data(6,j)-ycoo_med
	enddo
enddo

return
end