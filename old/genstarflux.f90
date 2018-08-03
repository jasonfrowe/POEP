subroutine genstarflux(npt,nstars,stars,stars2,xout,yout,nsub)
use precision
use startype
!import vars
integer :: npt,nstars,xout,yout,nsub
type(starpos), dimension(:,:) :: stars,stars2
!local vars
integer :: k
real(double) :: fratio,dmag
real(double), allocatable, dimension(:) :: fluxes,fluxes2

allocate(fluxes(nsub*nsub),fluxes2(nsub*nsub))

k=0
do i=1,nsub
	do j=1,nsub,4
		dmag=dble(j-1)*1.2d0
		fratio=10.0**(dmag/-2.5d0)

		k=k+1
		fluxes(k)=1.0d0
		fluxes2(k)=fratio*fluxes(k)
		if(i.eq.1) fluxes2(k)=0.0d0
		write(0,*) i,j,fluxes(k),fluxes2(k)

		k=k+1
		fluxes2(k)=1.0d0
		fluxes(k)=fratio*fluxes2(k)
		if(i.eq.1) fluxes(k)=0.0d0
		write(0,*) i,j,fluxes(k),fluxes2(k)

		k=k+1
		fluxes(k)=1.0d0
		fluxes2(k)=fratio*fluxes(k)
		if(i.eq.1) fluxes2(k)=0.0d0
		write(0,*) i,j,fluxes(k),fluxes2(k)

		k=k+1
		fluxes2(k)=1.0d0
		fluxes(k)=fratio*fluxes2(k)
		if(i.eq.1) fluxes(k)=0.0d0
		write(0,*) i,j,fluxes(k),fluxes2(k)

	enddo
 enddo

do j=1,npt
	do i=1,nstars
		stars(i,j)%flux = fluxes(i)
		stars2(i,j)%flux = fluxes2(i)
	enddo
enddo


return
end