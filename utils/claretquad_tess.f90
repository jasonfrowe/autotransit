program TESS_Limb_Darkening
use precision
implicit none
integer iargc,iostatus,npriors,ntype,i
real(double) :: teff,logg,feh
real(double), allocatable, dimension(:) :: priors
character(80) :: cline

interface
  function calcldprior(npriors,teff,logg,feh,ntype)
    use precision
    implicit none
    integer :: npriors,ntype
    real(double) :: teff, logg, feh
    real(double), dimension(:), pointer :: calcldprior
  end function calcldprior
end interface

!parameters 
npriors=6
!Limb-darkening type
ntype=0 !0=Claret u1,u2 ; 1=Kipping q1,q2

iostatus=0
if (iargc().lt.3) then !check command line arguments
  iostatus=1
else
  call getarg(1,cline)
    read(cline,*,iostat=iostatus) teff  !read in Teff
    if(iostatus>0) write(0,*) "*** Invalid Teff entry ***"
    call getarg(2,cline)
    read(cline,*,iostat=iostatus) logg  !read in log(g)
    if(iostatus>0) write(0,*) "*** Invalid log(g) entry ***"
    call getarg(3,cline)
    read(cline,*,iostat=iostatus) feh !read in [Fe/H]
    if(iostatus>0) write(0,*) "*** Invalid [Fe/H] entry ***"
endif

if (iargc().ge.4) then
     call getarg(4,cline)
    read(cline,*,iostat=iostatus) ntype  !read in Teff
    if(iostatus>0) write(0,*) "*** Invalid ntype entry ***"
endif

if (iostatus>0) then
  write(0,*) "Usage: claretquad Teff log(g) [Fe/H] <ntype>"
    write(0,*) "  Teff is in K"
    write(0,*) "  log(g) is cgs"
    write(0,*) "  [Fe/H] is [m/H]"
    write(0,*) "  <ntype> : 1=quadratic (u1,u2), 2=Kipping (q1,q2)"
    stop
endif

allocate(priors(npriors))
priors=0.0
! Priors(1) = u1 or q1
! Priors(2) = u2 or q2
! Priors(3) = (1)+1 sigma
! Priors(4) = (1)-1 sigma
! Priors(5) = (2)+1 sigma
! Priors(6) = (2)-1 sigma

priors=calcldprior(npriors,teff,logg,feh,ntype)

write(6,500) (priors(i),i=1,npriors)
500 format(6(F10.6,1X))

end program TESS_Limb_Darkening
