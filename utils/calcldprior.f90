function calcldprior(npriors,teff,logg,feh,ntype)
use precision
implicit none
!import vars
integer :: npriors,ntype
real(double) :: teff, logg, feh
real(double), dimension(:), pointer :: calcldprior
!local vars
type :: ldpars
	real(double), dimension(4) :: pars
	real(double), dimension(2) :: c
end type ldpars
integer, parameter :: nteff=79, nlogg=11, nfeh=19
integer :: i,j,k,filestatus,nunit,teffidx,loggidx,fehidx,nt,nl,nf,tidx(2),lidx(2), &
 fidx(2),n,ii,jj,kk,nld
real(double) :: teffs(nteff),loggs(nlogg),fehs(nfeh),cmin(2),diff,diffmin,c1(2), &
 c2(2),c3(2),c4(2),c5(2),c6(2),c7(2),c8(2),ld1,ld2,yout1,yout2,q1,q2
real(double), allocatable, dimension(:) :: pars1,cin1,x,y1,y2
type(ldpars), dimension(nteff,nlogg,nfeh) :: ldparsin
real(double), dimension(:), allocatable, target :: ldprior
character :: dumc
character(80) :: workdir,filename
data teffs/3500.,3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500., &
     5750.,6000.,6250.,6500.,6750.,7000.,7250.,7500.,7750.,8000., &
     8250.,8500.,8750.,9000.,9250.,9500.,9750.,10000.,10250.,10500., &
     10750.,11000.,11250.,11500.,11750.,12000.,12250.,12500.,12750., &
     13000.,14000.,15000.,16000.,17000.,18000.,19000.,20000.,21000., &
     22000.,23000.,24000.,25000.,26000.,27000.,28000.,29000.,30000., &
     31000.,32000.,33000.,34000.,35000.,36000.,37000.,37500.,38000., &
     39000.,40000.,41000.,42000.,42500.,43000.,44000.,45000.,46000., &
     47000.,47500.,48000.,49000.,50000./
data loggs/0.,0.5,1,1.5,2.,2.5,3,3.5,4.,4.5,5./
data fehs/-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.3, &
     -0.2,-0.1,0.0,0.1,0.2,0.3,0.5,1.0/

allocate(ldprior(npriors))
ldprior=0.0d0

!Open up limb-darkening table
!workdir="/data/Kepler/tables/"
workdir="./tables/"
filename="Claret-limbquad-tess.txt"  !name of Claret tables

nunit=10 
!First see if file in current working directory.
open(unit=nunit,file=trim(workdir)//trim(filename),iostat=filestatus,status='old')
!If not, then try workdir
if(filestatus>0) open(unit=nunit,file=filename,iostat=filestatus,status='old')
!if we still can't find the find, then we spit out an erro
if(filestatus>0)then !trap missing file errors
   write(0,'(A12,A80)') "Cannot open ",filename
   write(0,'(A13,A80)') " tried . and ",workdir
   stop
endif

!read in Table
allocate(pars1(4),cin1(2))
do i=1,nteff
	do j=1,nlogg
		do k=1,nfeh
			ldparsin(i,j,k)%c(1)=-10.0 !give bad values as default.
      ldparsin(i,j,k)%c(2)=-10.0 
		enddo
	enddo
enddo

do i=1,22
	read(nunit,*) dumc !skip header
enddo

i=0
do
   	read(nunit,*,iostat=filestatus) (pars1(j),j=1,4),(cin1(j),j=1,2)
	if(filestatus == 0) then

      	!find index to fill in limb-darkening results.
      	call locate(teffs,nteff,pars1(2),nt)
      	call locate(loggs,nlogg,pars1(1),nl)
      	call locate(fehs,nfeh,pars1(3),nf)

        ldparsin(nt,nl,nf)%pars(1:4)=pars1(1:4)
        ldparsin(nt,nl,nf)%c(1:2)=cin1(1:2)

   	elseif(filestatus == -1) then
      	exit  !successively break from data read loop.
   	else
      	write(0,*) "File Error!! Line:",i+1
      	write(0,900) "iostat: ",filestatus
      	900 format(A8,I3)
      	stop
   	endif
enddo
close(nunit) !close file

allocate(x(nteff),y1(nteff),y2(nteff))

do i=1,nteff
  do j=1,nlogg
    do k=1,nfeh
      if(ldparsin(i,j,k)%c(1) .le. -9.0)then  
        
        ld1=0.0d0
        ld2=0.0d0
        nld=0
 
        n=0
        do ii=1,nteff
          if(ldparsin(ii,j,k)%c(1) .gt. -9.0)then
            n=n+1
            x(n)=teffs(ii)
            y1(n)=ldparsin(ii,j,k)%c(1)
            y2(n)=ldparsin(ii,j,k)%c(2)
          endif
        enddo
        if (n.gt.1) then
          call lininterp(x,y1,n,teffs(i),yout1)
          call lininterp(x,y2,n,teffs(i),yout2)
          ld1=ld1+yout1
          ld2=ld2+yout2
          nld=nld+1
          !write(0,*) "new: ",n,yout1,yout2
        endif

        n=0
        do jj=1,nlogg
          if(ldparsin(i,jj,k)%c(1) .gt.-9.0)then
            n=n+1
            x(n)=loggs(jj)
            y1(n)=ldparsin(i,jj,k)%c(1)
            y2(n)=ldparsin(i,jj,k)%c(2)
          endif
        enddo
        if (n.gt.1) then
          call lininterp(x,y1,n,loggs(j),yout1)
          call lininterp(x,y2,n,loggs(j),yout2)
          ld1=ld1+yout1
          ld2=ld2+yout2
          nld=nld+1
          !write(0,*) "new: ",n,yout1,yout2
        endif

        n=0
        do kk=1,nfeh
          if(ldparsin(i,j,kk)%c(1) .gt.-9.0)then
            n=n+1
            x(n)=fehs(kk)
            y1(n)=ldparsin(i,j,kk)%c(1)
            y2(n)=ldparsin(i,j,kk)%c(2)
            !write(0,*) "yo3",n,x(n),y1(n),y2(n)
          endif
        enddo
        if (n.gt.1) then
          call lininterp(x,y1,n,fehs(k),yout1)
          call lininterp(x,y2,n,fehs(k),yout2)
          ld1=ld1+yout1
          ld2=ld2+yout2
          nld=nld+1
          !write(0,*) "new: ",n,yout1,yout2
        endif

        if (nld.gt.0) then
          ldparsin(i,j,k)%c(1)=ld1/dble(nld)
          ldparsin(i,j,k)%c(2)=ld2/dble(nld)
          !write(0,*) "new: ",ldparsin(i,j,k)%c(1),ldparsin(i,j,k)%c(2)
        else
          write(0,*) teffs(i),loggs(j),fehs(k)
          write(0,*) ldparsin(i,j,k)%c(1),ldparsin(i,j,k)%c(2)
          write(0,*) "No values can be found"
        endif

        !read(5,*)
      endif
    enddo
  enddo
enddo

!find index around stellar parameters
call locate(teffs,nteff,teff,teffidx)
call locate(loggs,nlogg,logg,loggidx)
call locate(fehs,nfeh,feh,fehidx)

!calculate min/max bracketing indices.
!taking edges into account
if(teffidx.ge.nteff)then
	tidx(1)=nteff-1
	tidx(2)=nteff
else
	tidx(1)=teffidx
	tidx(2)=teffidx+1
endif
if(loggidx.ge.nlogg)then
	lidx(1)=nlogg-1
	lidx(2)=nlogg
else
	lidx(1)=loggidx
	lidx(2)=loggidx+1
endif
if(fehidx.ge.nfeh)then
	fidx(1)=nfeh-1
	fidx(2)=nfeh
else
	fidx(1)=fehidx
	fidx(2)=fehidx+1
endif

c1(1:2)=ldparsin(tidx(1),lidx(1),fidx(1))%c(1:2)
c2(1:2)=ldparsin(tidx(1),lidx(1),fidx(2))%c(1:2)
c3(1:2)=ldparsin(tidx(2),lidx(1),fidx(2))%c(1:2)
c4(1:2)=ldparsin(tidx(2),lidx(1),fidx(1))%c(1:2)
c5(1:2)=ldparsin(tidx(1),lidx(2),fidx(1))%c(1:2)
c6(1:2)=ldparsin(tidx(1),lidx(2),fidx(2))%c(1:2)
c7(1:2)=ldparsin(tidx(2),lidx(2),fidx(2))%c(1:2)
c8(1:2)=ldparsin(tidx(2),lidx(2),fidx(1))%c(1:2)

do k=1,2
	call trilinear(c1(k),c2(k),c3(k),c4(k),c5(k),c6(k),c7(k),c8(k),cmin(k), &
		logg,teff,feh,loggs(lidx(1)),loggs(lidx(2)),teffs(tidx(1)), &
		teffs(tidx(2)),fehs(fidx(1)),fehs(fidx(2)))
enddo

!Convert to q1,q2 if necessary.  
if (ntype.eq.2) then
        q1=(cmin(1)+cmin(2))**2.0
        q2=cmin(1)/(2.0*(cmin(1)+cmin(2)))
        ldprior(1)=q1
        ldprior(2)=q2
else
        ldprior(1)=cmin(1)
        ldprior(2)=cmin(2)
endif



calcldprior => ldprior


return
end function calcldprior
