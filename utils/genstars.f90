subroutine genstars(npt,nf,nstars,data,stars,xout,yout,nsub,seed)
use precision
use startype
implicit none
!import vars
integer :: npt,nf,nstars,xout,yout,nsub,seed
real(double), dimension(:,:) :: data
type(starpos), dimension(:,:) :: stars
!local vars
integer :: nstarmax,nst,ntd,i,j,k,np,ix,iy,nsep,ndmag,nrast,l,ntdur,ii
real(double) :: dx,dy,fratio,th,pi,tpi,ran2,xcoo_med,ycoo_med,median,RedRs, &
 radialvel,albedo,tides
real(double), allocatable, dimension(:) :: rhostar,Ag,rstar,per,fluxes, &
 xcoo,ycoo,rpl,sep,dmag,tdur,ld1,ld2,mstar,mpl
!transit model parameters
integer :: nplanet,nfit,nplanetmax,nmax
integer, allocatable, dimension(:) :: dtype,ntt
real(double), allocatable, dimension(:) :: sol,time,itime,tmodel
real(double), allocatable, dimension(:,:) :: tobs,omc

!Constants
Pi=acos(-1.d0)!define Pi and 2*Pi
tPi=2.0d0*Pi
RedRs = 0.009158 !radius of Sun / Radius of Earth 

!initialize variables for transit model.
nmax=npt !maximum datapoints we will through at transitmodel
nplanet=1 !only consider single planetary systems
nplanetmax=nplanet
allocate(ntt(nplanet),tobs(nplanet,npt),omc(nplanet,npt))
ntt=0 !no TTVs
tobs=0.0
omc=0.0
nfit=108 !max parameters that transitmodel can use.
allocate(sol(nfit))
allocate(dtype(npt),time(npt),itime(npt),tmodel(npt))
dtype=0 !mark all data as photometry 
time(1:npt)=data(1,1:npt) !copy time stamps into 1D array 
itime=median(npt,data(15,1:npt))/86400.0d0 !integration time for transit model.

!estimate median of x,y centroids
xcoo_med=median(npt,data(5,1:npt))
ycoo_med=median(npt,data(6,1:npt))

nstarmax=size(stars%flux,1)  !get size of stars storage to pervent overflow
nstars=0 !counting the number of stars
nrast=0  !counting raster position

!spacing for primary stars
dx=dble(xout)/dble(nsub)
dy=dble(yout)/dble(nsub)

allocate(fluxes(nstarmax),xcoo(nstarmax),ycoo(nstarmax))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!generate single stars with transits!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nst=3
allocate(rhostar(nst),Ag(nst),rstar(nst),mstar(nst),ld1(nst),ld2(nst))
rhostar(1)=1.4
rhostar(2)=1.4
rhostar(3)=50.0
Ag(1)=0.05 !geometric albedo will also set secondary eclipse depth (no thermal)
Ag(2)=0.40
Ag(3)=0.00
rstar(1)=1.0   !stellar radius (Rsun)
rstar(2)=1.0
rstar(3)=0.1
mstar(1)=1.0   !stellar mass (Msun)
mstar(2)=1.0
mstar(3)=0.04
ld1(1)=0.3946  !quadratic limb-darkening
ld2(1)=0.2671
ld1(2)=0.3946
ld2(2)=0.2671
ld1(3)=0.3360
ld2(3)=0.2147
sol(4)=0.0d0 !using quad limb-darkening, so these are set to zero.
sol(5)=0.0d0 
sol(6)=0.0d0 !no dilution needed
sol(7)=0.0d0 !no RVs
sol(8)=0.0d0 !set photometric zero point to zero
sol(13)=0.0d0 !eccentricity
sol(14)=0.0d0
sol(15)=0.0d0 !RV amplitude (not needed)
sol(16)=0.0d0 !secondary eclipse, default is zero. 
sol(17)=0.0d0 !ellipitcal amp
sol(18)=0.0d0 !albedo amp
ntd=4
allocate(rpl(ntd),mpl(ntd))
rpl(1)=10.0 !planet radius (Rearth)
rpl(2)= 4.0
rpl(3)= 1.6
rpl(4)= 1.0
mpl(1)=317.8 !mass of planet (Mearth)
mpl(2)=20.0
mpl(3)=8.0
mpl(4)=1.0
np=3
allocate(per(np))
per(1)=3.0
per(2)=6.0
per(3)=9.0
do i=1,nst !loop over stellar types
	sol(1)=rhostar(i) !set mean stellar density
	sol(2)=ld1(i)
	sol(3)=ld2(i)
	do j=1,ntd !loop over transit-depths
		sol(12)=rpl(j)*RedRs/rstar(i) !Rp/Rs 
		do k=1,np !loop over periods
			nstars=nstars+1
			nrast=nrast+1
			if(nstars.le.nstarmax)then
				!call transit model
				sol(10)=per(k)+1.5d0*ran2(seed) !set orbital period, with extra jitter
				sol(9)=sol(10)*ran2(seed) !set T0 as random phase from 0 to 1
				sol(11)=0.8*ran2(seed) !set impact parameter [0-0.8]
				sol(15)=radialvel(mstar(i),mpl(j),sol(10),sol(11),sol(1))
				sol(18)=albedo(Ag,Rpl(j),sol(10),sol(1),rstar(i))
				sol(16)=sol(18) !secondary is equal to phase curve (no thermal)
				sol(17)=tides(mstar(i),rstar(i),mpl(j),sol(1),sol(10),sol(11))
				!write(0,*) "Tides: ",sol(17)
				call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime, &
 				  ntt,tobs,omc,tmodel,dtype)
				do ii=1,npt
					stars(nstars,ii)%flux=tmodel(ii)
					!write(0,*) tmodel(ii)
					!read(5,*)
				enddo
				ix=(nrast-1)/nsub+1
				iy=nrast-nsub*(ix-1)
				xcoo(nstars)=dx*(ix-1)+dx/2.0d0+ran2(seed) !jitter to move off center
				ycoo(nstars)=dy*(iy-1)+dy/2.0d0+ran2(seed)
				!write(0,*) nstars,xcoo(nstars),ycoo(nstars)
			else 
				write(0,*) "Error, nstars > nstarmax"
				stop
			endif
			!read(5,*)
		enddo
	enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!generate foreground blends (bright star has signal)!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nsep=3
allocate(sep(nsep)) !separation in pixels 
sep(1)=0.5
sep(2)=2.0
sep(3)=5.0
ndmag=4
allocate(dmag(ndmag))
dmag(1)=0.0
dmag(2)=3.0
dmag(3)=9.0
dmag(4)=13.0
!Transit model parameters - using Sun as FG star
sol(1)=rhostar(1) !set mean stellar density
sol(2)=ld1(1)
sol(3)=ld2(1)
sol(4)=0.0d0 !using quad limb-darkening, so these are set to zero.
sol(5)=0.0d0 
sol(6)=0.0d0 !no dilution needed
sol(7)=0.0d0 !no RVs
sol(8)=0.0d0 !set photometric zero point to zero
sol(13)=0.0d0 !eccentricity
sol(14)=0.0d0
sol(15)=0.0d0 !RV amplitude (not needed)
sol(16)=0.0d0 !secondary eclipse, default is zero. 
sol(17)=0.0d0 !ellipitcal amp
sol(18)=0.0d0 !albedo amp
do i=1,nsep !loop over separations
	do j=1,ndmag
		fratio=10.0**(dmag(j)/-2.5d0)
		do k=1,ntd !loop over transit depths
			sol(12)=rpl(k)*RedRs/rstar(1) !Rp/Rs 
			do l=1,np !loop over periods
				nrast=nrast+1
				ix=(nrast-1)/nsub+1
				iy=nrast-nsub*(ix-1)
				
				sol(10)=per(l)+1.5d0*ran2(seed) !set orbital period, with extra jitter
				sol(9)=sol(10)*ran2(seed) !set T0 as random phase from 0 to 1
				sol(11)=0.8*ran2(seed) !set impact parameter [0-0.8]
				sol(15)=radialvel(mstar(1),mpl(k),sol(10),sol(11),sol(1))
				sol(18)=albedo(Ag,Rpl(i),sol(10),sol(1),rstar(i))
				sol(16)=sol(18) !secondary is equal to phase curve (no thermal)
				sol(17)=tides(mstar(i),rstar(i),mpl(k),sol(1),sol(10),sol(11))
				call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime, &
 				  ntt,tobs,omc,tmodel,dtype)

				nstars=nstars+1
				if(nstars.le.nstarmax)then
					xcoo(nstars)=dx*(ix-1)+dx/2.0d0+ran2(seed) !jitter to move off center
					ycoo(nstars)=dy*(iy-1)+dy/2.0d0+ran2(seed)
					do ii=1,npt
						stars(nstars,ii)%flux=tmodel(ii)
						!write(0,*) tmodel(ii)
						!read(5,*)
					enddo
				else
					write(0,*) "Error, nstars > nstarmax"
					stop
				endif

				!write(0,*) nstars,xcoo(nstars),ycoo(nstars)

				nstars=nstars+1
				if(nstars.le.nstarmax)then
					th=tPi*ran2(seed)
					xcoo(nstars)=xcoo(nstars-1)+sep(i)*cos(th)
					ycoo(nstars)=ycoo(nstars-1)+sep(i)*sin(th)
					do ii=1,npt
						stars(nstars,ii)%flux=fratio
						!write(0,*) tmodel(ii)
						!read(5,*)
					enddo
				else
					write(0,*) "Error, nstars > nstarmax"
					stop
				endif

				!write(0,*) nstars,xcoo(nstars),ycoo(nstars)

			enddo
		enddo
	enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generate Background Blends!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ntdur=2
allocate(tdur(ntdur))
tdur(1)=1.0 !transit-duration in hours
tdur(2)=3.0 !transit-duration in hours
!generate false positives
do i=1,nsep !loop over separations
	do j=1,ntd !loop over transit depths
		sol(12)=(0.5+ran2(seed))*RedRs/rstar(3)  !1 to 1.5 Rearth around a BD
		fratio=(rpl(1)*RedRs/rstar(1))**2.0d0 / (sol(12))**2.0d0
		write(0,*) "fratio: ",fratio
		read(5,*)
		do k=1,ntdur !loop over transit durations
			do l=1,np !loop over periods
				nrast=nrast+1
				ix=(nrast-1)/nsub+1
				iy=nrast-nsub*(ix-1)
				
				nstars=nstars+1
				if(nstars.le.nstarmax)then
					xcoo(nstars)=dx*(ix-1)+dx/2.0d0+ran2(seed) !jitter to move off center
					ycoo(nstars)=dy*(iy-1)+dy/2.0d0+ran2(seed)
					fluxes(nstars)=1.0d0
				else
					write(0,*) "Error, nstars > nstarmax"
					stop
				endif

				nstars=nstars+1
				if(nstars.le.nstarmax)then
					th=tPi*ran2(seed)
					xcoo(nstars)=xcoo(nstars-1)+sep(i)*cos(th)
					ycoo(nstars)=ycoo(nstars-1)+sep(i)*sin(th)
					fluxes(nstars)=1.0d0*fratio
				else
					write(0,*) "Error, nstars > nstarmax"
					stop
				endif
			enddo
		enddo
	enddo
enddo


!add in contant stars
j=nrast+1
write(0,*) "ns ",j,nstarmax
do i=j,nsub*nsub
	nrast=nrast+1
	ix=(nrast-1)/nsub+1
	iy=nrast-nsub*(ix-1)
				
	nstars=nstars+1
	if(nstars.le.nstarmax)then
		xcoo(nstars)=dx*(ix-1)+dx/2.0d0+ran2(seed) !jitter to move off center
		ycoo(nstars)=dy*(iy-1)+dy/2.0d0+ran2(seed)
		fluxes(nstars)=1.0d0
	else
		write(0,*) "Error, nstars > nstarmax"
		stop
	endif
	!write(0,*) nstars,xcoo(nstars),ycoo(nstars)
enddo

!generate full Time-series using pointing model.
do j=1,npt
	do i=1,nstars
		stars(i,j)%xcoo = xcoo(i)+data(5,j)-xcoo_med
		stars(i,j)%ycoo = ycoo(i)+data(6,j)-ycoo_med
		!stars(i,j)%flux = fluxes(i)
	enddo
enddo


return
end subroutine genstars