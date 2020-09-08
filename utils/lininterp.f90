subroutine lininterp(x,y,n,xin,yout)
!linear interpolation.  X is expected to be in order.
use precision
implicit none
integer :: n,i
real(double) :: small
real(double) :: x(n),y(n),xin,yout
 
small=1.5d-8 !to catch identifical values

if(n.eq.1)then
  yout=y(1)
  return
endif

if(xin.le.x(1))then
  i=1
  if(x(i+1)-x(i).lt.small)then
  	yout=y(i)
  else
  	yout=y(i)+(xin-x(i))/(x(i+1)-x(i))*(y(i+1)-y(i))
  endif
elseif(xin.ge.x(n))then
  i=n-1
  if(x(i+1)-x(i).lt.small)then
  	yout=y(i)
  else
  	yout=y(i)+(xin-x(i))/(x(i+1)-x(i))*(y(i+1)-y(i))
  endif
else
  do i=1,n-1
    if((xin.ge.x(i)).and.(xin.lt.x(i+1)))then
    	if(x(i+1)-x(i).lt.small)then
    		yout=y(i) !avoid division by zero.
    	else
      		yout=y(i)+(xin-x(i))/(x(i+1)-x(i))*(y(i+1)-y(i))
      	endif
    endif
  enddo
endif

return
end
