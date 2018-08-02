function albedo(Ag,Rpl,Per,rhostar,rstar)
!returns Fp/F* given albedo and orbtial properties
use precision
implicit none
!import vars
real(double) :: Ag,Rpl,Per,rhostar,albedo,rstar
!local vars
real(double) :: Pi,tPi,pid2,G,adrs,Psec,fpdfs,Rearth,Rsun

!constants
Pi=acos(-1.d0)!define Pi and 2*Pi
tPi=2.0d0*Pi 
pid2=Pi/2.0d0
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
Rearth=6.378136d6 !radius of Earth (m)
Rsun=6.95508d8    !radius of Sun (m)

Psec=Per*86400.0 !convert period in days to period in seconds

!first get a/R*
adrs=1000.0*rhostar*G*(Psec)**2/(3.0d0*Pi)
adrs=adrs**(1.0d0/3.0d0)

fpdfs = Ag*(Rpl*Rearth/(adrs*rstar*Rsun))**2.0d0

albedo=fpdfs*1.0d6  !convert to ppm

return
end