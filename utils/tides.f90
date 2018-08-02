function tides(mstar,rstar,mpl,rhostar,Per,b)
use precision
implicit none
!import vars
real(double) :: tides,mstar,mpl,rhostar,Per,rstar,b
!local vars
integer :: l
real(double) :: Pi,tPi,pid2,G,Msun,Mearth,M1,M2,eps,adrs,Rearth,Rsun, &
 incl,ra(3),ad(3),lambda(3),f(3),P(3),dJ(3),d,Psec,cos2incl,sin2incl, &
 R1,asemi


!constants
Pi=acos(-1.d0)!define Pi and 2*Pi
tPi=2.0d0*Pi 
pid2=Pi/2.0d0
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
Msun=1.9891d30 !mass of Sun (kg)
Mearth=5.9736e24 !mass of Earth (kg)
Rearth=6.378136d6 !radius of Earth (m)
Rsun=6.95508d8    !radius of Sun (m)

tides=0.0

M1=Msun*mstar
M2=Mearth*mpl
R1=Rsun*rstar

Psec=Per*86400.0 !convert period in days to period in seconds

!get a/R*
adrs=1000.0*rhostar*G*(Psec)**2/(3.0d0*Pi)
adrs=adrs**(1.0d0/3.0d0)

eps=M2/M1*(1.0d0/adrs)**3.0d0
d=adrs*rstar*Rsun
asemi=d !assume circular orbit.
incl=acos(1.0d0/adrs) !90 deg is seen edge on for RV
cos2incl=cos(incl)*cos(incl)
sin2incl=sin(incl)*sin(incl)

do l=2,3
	ra(l)=(1.0d0/adrs)**(l-2)
    ad(l)=(asemi/d)**(l+1)
    lambda(l)=dble(l)+2.0d0
enddo


f(2)=-1.3d1*(1.0d0+lambda(2)/4.0d0)/1.0d1
f(3)=-5.0d0*(1.0d0+lambda(3)/1.0d1)/8.0d0
P(2)=2.5d-1*(-(3.0d0*cos2incl-1.0d0)+3.0d0*sin2incl)
P(3)=1.25d-1*sin(incl)*(-3.0d0*(5.0d0*cos2incl-1.0d0)+5.0d0*sin2incl)

do l=2,3
	dJ(l)=ra(l)*ad(l)*f(l)*P(l)
    tides=tides+dJ(l)
enddo

tides=abs(tides*eps)*1.0d6 !convert to ppm

return
end