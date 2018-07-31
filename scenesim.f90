program scenesimulator
!Code for generating FITS files given a PSF and point positions
use precision
implicit none
integer :: xout, yout
character(80) :: prefix,fileout,fileout_c,fileout_m,filename

!Image Dimensions
xout=1024
yout=1024

!Name for output
prefix="sgen"
fileout=trim(prefix)//".fits"
fileout_m=trim(prefix)//"_m.txt"
fileout_c=trim(prefix)//"_c.fits"

!input file
filename="scenemodel_20180731.dat"

call readmodel(filename)


end program scenesimulator