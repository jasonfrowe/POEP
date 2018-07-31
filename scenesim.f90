program scenesimulator
!Code for generating FITS files given a PSF and point positions
use precision
implicit none
integer :: xout,yout,nunit,filestatus
character(80) :: prefix,fileout,fileout_c,fileout_m,modelfile

!Image Dimensions
xout=1024
yout=1024

!Name for output
prefix="sgen"
fileout=trim(prefix)//".fits"
fileout_m=trim(prefix)//"_m.txt"
fileout_c=trim(prefix)//"_c.fits"

!input file
modelfile="scenemodel_20180731.dat"

nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

call readmodel(nunit,filestatus)


end program scenesimulator