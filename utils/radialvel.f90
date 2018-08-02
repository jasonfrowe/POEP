function radialvel(mstar,mpl,Per,b,rhostar)
use precision
implicit none
!import vars
real(double) :: radialvel,mstar,mpl,Per,b,rhostar
!local vars
real(double) :: Pi,tPi,pid2,G,adrs,Msun,Mearth,Psec,vr3,incl


!constants
Pi=acos(-1.d0)!define Pi and 2*Pi
tPi=2.0d0*Pi 
pid2=Pi/2.0d0
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
Msun=1.9891d30 !mass of Sun (kg)
Mearth=5.9736e24 !mass of Earth (kg)

Psec=Per*86400.0 !convert period in days to period in seconds

!calculate inclination
!first get a/R*
adrs=1000.0*rhostar*G*(Psec)**2/(3.0d0*Pi)
adrs=adrs**(1.0d0/3.0d0)
!then get inclination
incl=acos(b/adrs) !90 deg is seen edge on for RV

vr3=tpi*G/Psec*(mpl*Mearth)**3.0d0/(mpl*Mearth+mstar*Msun)**2.0 * sin(incl)**3.0

radialvel=vr3**(1.0d0/3.0d0)

return
end


