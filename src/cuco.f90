!                                            ********************************
!                                            *         module cuco
!                                            *
!                                            *        R. J. Purser
!                                            *
!                                            *        NOAA/NCEP/EMC
!                                            * 
!                                            *           2017
!                                            *
!                                            *      jim.purser@noaa.gov
!                                            *
!                                            ********************************
!
! Map transformation routines for the conformal cubic grid
! Add to the library, libcuco.a
!
! The index convention for the six square panels (or "tiles") is such
! that,
! when the cube is developed (i.e., unfolded), the panels can be aligned
! with their xm(1) and xm(2) coordinate directions to the right and
! upwards
! in the following arrangement:
! 
! +-------+
! |       |
! |   5   | <-- North pole at xc=(0,0,1)
! |       |
! +-------+-------+-------+-------+
! |       |       |       |       |
! |   1   |   2   |   3   |   4   | <-- Equator, xc(3)=0
! |       |       |       |       |
! +-------+-------+-------+-------+
! |       |
! |   6   | <-- South pole at xc=(0,0,-1)
! |       |
! +-------+
!     |
!     xc(2)=0, the prime meridian.  (Center of panel 1 is at xc=(1,0,0))
! In each panel, the map coordinates xm(1) and xm(2) run between -1 and
! +1
! (the origin being the center of the square. The earth-centered
! cartesians
! are scaled to be unit vectors. The derivatives of xc wrt xm form the
! 3*2 Jacobian matrix dxcdxm and, going the other way, the 2*3 matrix
! dxmdxc of the inverse mapping is the generalized inverse of dxcdxm.
!
! If one wishes to go from xm to xc, but with an xm that might be
! slightly
! over the edge of the map panel ipan, then one should first run xmtoxm
! to find the actual map panel, jpan, that contains the location, and
! use
! the returned new coordinates for the effective xm. The Jacobian of the
! needed affine transformation involved is also returned, which allows
! proper adjustments to be applied (by the chain rule) to xcd, if
! needed.
!
! DIRECT DEPENDENCIES:
! Modules: pmat, pmat4, pcoco, pietc, pkind
! Libraries: pcoco, pmat
!=============================================================================
module cuco
!=============================================================================
! Routines related to the convenient computation of conformal cubic
! mapping
!=============================================================================
use pkind, only: spi,sp,dp,dpc
use pietc, only: T,F,u0,u1,u2
use pcoco, only: set_cuco,cuztoc
use pmat4, only: ztoc
implicit none
private
public:: mcuco,co1,co2, inicuco,xmtoxc,xctoxm,xmtoxm
integer(spi),parameter   :: mcuco=50
real(dp),dimension(mcuco):: co1,co2
real(dp),dimension(3,3,6):: rotp0
logical                  :: lnotinicuco=T ! true when NOT initialized yet
!=============================================================================
interface inicuco;  module procedure inicuco;                    end interface
interface xmtoxc
   module procedure xmtoxc,xmtoxcd,xmtoxc_s,xmtoxcd_s
                                                                end interface
interface xctoxm;  module procedure xctoxm,xctoxmd,xctoxm_s,xctoxmd_s
                                                                end interface
interface xmtoxm;  module procedure xmtoxm;                     end interface
interface reorient;module procedure reorient;                   end interface 
contains

!=============================================================================
subroutine inicuco!  [inicuco]
!=============================================================================
! Initialize conformal cubic expansion coefficients and the rotations,
! rotp0
! for each of the six panel coordinate frames.
!=============================================================================
use pcoco, only: set_cuco
implicit none
!=============================================================================
call set_cuco(mcuco,co1,co2)

rotp0(:,1,1)=(/u0,u1,u0/);rotp0(:,2,1)=(/u0,u0,u1/);rotp0(:,3,1)=(/u1,u0,u0/)
rotp0(:,1,2)=(/-u1,u0,u0/);rotp0(:,2,2)=(/u0,u0,u1/);rotp0(:,3,2)=(/u0,u1,u0/)
rotp0(:,1,3)=(/u0,-u1,u0/);rotp0(:,2,3)=(/u0,u0,u1/);rotp0(:,3,3)=(/-u1,u0,u0/)
rotp0(:,1,4)=(/u1,u0,u0/);rotp0(:,2,4)=(/u0,u0,u1/);rotp0(:,3,4)=(/u0,-u1,u0/)
rotp0(:,1,5)=(/u0,u1,u0/);rotp0(:,2,5)=(/-u1,u0,u0/);rotp0(:,3,5)=(/u0,u0,u1/)
rotp0(:,1,6)=(/u0,u1,u0/);rotp0(:,2,6)=(/u1,u0,u0/);rotp0(:,3,6)=(/u0,u0,-u1/)
lnotinicuco=F
end subroutine inicuco

!=============================================================================
subroutine xmtoxc(xm,xc,ipan)!  [xmtoxc]
!=============================================================================
! From panel index, ipan, and map coordinates, xm, return earth-centered
! cartesian vector, xc.
!=============================================================================
implicit none
real(dp),dimension(2),intent(in ):: xm
real(dp),dimension(3),intent(out):: xc
integer(spi),         intent(in ):: ipan
!-----------------------------------------------------------------------------
complex(dpc)           :: z
!=============================================================================
if(lnotinicuco)call inicuco
z=cmplx(xm(1)+u1,xm(2)+u1,dpc)/u2
call cuztoc(mcuco,co1,co2,z,xc)
xc=matmul(rotp0(:,:,ipan),xc)
end subroutine xmtoxc

!=============================================================================
subroutine xmtoxcd(xm,xc,dxcdxm,ipan)!  [xmtoxc]
!=============================================================================
implicit none
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: dxcdxm
integer(spi),           intent(in ):: ipan
!-----------------------------------------------------------------------------
complex(dpc),dimension(3):: xcd
complex(dpc)             :: z
integer(spi)             :: i
!=============================================================================
if(lnotinicuco)call inicuco
z=cmplx(xm(1)+u1,xm(2)+u1,dpc)/u2
call cuztoc(mcuco,co1,co2,z,xc,xcd)
do i=1,3; dxcdxm(i,1)=real(xcd(i)); dxcdxm(i,2)=-aimag(xcd(i)); enddo
dxcdxm=matmul(rotp0(:,:,ipan),dxcdxm)/u2
xc=matmul(rotp0(:,:,ipan),xc)
end subroutine xmtoxcd
!==============================================================================
subroutine xmtoxc_s(xm_s,xc_s,ipan)!  [xmtoxc]
!==============================================================================
! Single precision wrapper of xmtoxc
!==============================================================================
implicit none
real(sp),dimension(2),intent(in ):: xm_s
real(sp),dimension(3),intent(out):: xc_s
integer(spi),         intent(in ):: ipan
!------------------------------------------------------------------------------
real(dp),dimension(2):: xm
real(dp),dimension(3):: xc
!==============================================================================
xm=xm_s
call xmtoxc(xm,xc,ipan)
xc_s=xc
end subroutine xmtoxc_s
!==============================================================================
subroutine xmtoxcd_s(xm_s,xc_s,dxcdxm_s,ipan)!  [xmtoxc]
!==============================================================================
! Single precision wrapper of xmtoxcd
!==============================================================================
implicit none
real(sp),dimension(2),  intent(in ):: xm_s
real(sp),dimension(3),  intent(out):: xc_s
real(sp),dimension(3,2),intent(out):: dxcdxm_s
integer(spi),           intent(in ):: ipan
!------------------------------------------------------------------------------
real(dp),dimension(2)  :: xm
real(dp),dimension(3)  :: xc
real(dp),dimension(3,2):: dxcdxm
!==============================================================================
xm=xm_s
call xmtoxc(xm,xc,dxcdxm,ipan)
xc_s=xc
dxcdxm_s=dxcdxm
end subroutine xmtoxcd_s

!==============================================================================
subroutine xctoxm(xc, xm,ipan)!  [xctoxm]
!==============================================================================
use pmat, only: inv
implicit none
real(dp),dimension(3),intent(in ):: xc
real(dp),dimension(2),intent(out):: xm
integer(spi),         intent(out):: ipan
!------------------------------------------------------------------------------
integer(spi),parameter   :: nit=20
real(dp),    parameter   :: epss=1.e-20_dp
real(dp),dimension(3,2)  :: dxcdxm
real(dp),dimension(2,2)  :: mat2
real(dp),dimension(2,3)  :: axx
real(dp),dimension(3,6)  :: cens
real(dp),dimension(6)    :: p6
real(dp),dimension(3)    :: xxc,xct,rxc
integer(spi),dimension(1):: ii
integer(spi)             :: it
!=============================================================================
if(lnotinicuco)call inicuco
cens=rotp0(:,3,:)
xxc=xc/sqrt(dot_product(xc,xc))
p6=matmul(xxc,cens)
ii=maxloc(p6)
ipan=ii(1)
xct=matmul(transpose(rotp0(:,:,ipan)),xxc)
xm=xct(1:2)/xct(3)
do it=1,nit
   call xmtoxc(xm,xct,dxcdxm,ipan)
   rxc=xct-xxc
   axx=transpose(dxcdxm)
   mat2=matmul(axx,dxcdxm)
   call inv(mat2)
   axx=matmul(mat2,axx)
   xm=xm-matmul(axx,rxc)
   if(dot_product(rxc,rxc)<epss)exit
enddo
if(it>nit)print'("Warning; In xctoxm, Newton iterations not converging")'
end subroutine xctoxm
   
!==============================================================================
subroutine xctoxmd(xc, xm,dxmdxc,ipan)!  [xctoxm]
!==============================================================================
use pmat, only: inv
implicit none
real(dp),dimension(3),  intent(in ):: xc
real(dp),dimension(2),  intent(out):: xm
real(dp),dimension(2,3),intent(out):: dxmdxc
integer(spi),           intent(out):: ipan
!------------------------------------------------------------------------------
integer(spi),parameter   :: nit=20
real(dp),    parameter   :: epss=1.e-20_dp
real(dp),dimension(3,2)  :: dxcdxm
real(dp),dimension(2,2)  :: mat2
real(dp),dimension(2,3)  :: axx
real(dp),dimension(3,6)  :: cens
real(dp),dimension(6)    :: p6
real(dp),dimension(3)    :: xxc,xct,rxc
integer(spi),dimension(1):: ii
integer(spi)             :: it
!=============================================================================
if(lnotinicuco)call inicuco
cens=rotp0(:,3,:)
xxc=xc/sqrt(dot_product(xc,xc))
p6=matmul(xxc,cens)
ii=maxloc(p6)
ipan=ii(1)
xct=matmul(transpose(rotp0(:,:,ipan)),xxc)
xm=xct(1:2)/xct(3)
do it=1,nit
   call xmtoxc(xm,xct,dxcdxm,ipan)
   rxc=xct-xxc
   axx=transpose(dxcdxm)
   mat2=matmul(axx,dxcdxm)
   call inv(mat2)
   dxmdxc=matmul(mat2,axx)
   xm=xm-matmul(dxmdxc,rxc)
   if(dot_product(rxc,rxc)<epss)exit
enddo
if(it>nit)print'("Warning; In xctoxmd, Newton iterations not converging")'
end subroutine xctoxmd
!==============================================================================
subroutine xctoxm_s(xc_s,xm_s,ipan)!  [xctoxm]
!==============================================================================
implicit none
real(sp),dimension(3),intent(in ):: xc_s
real(sp),dimension(2),intent(out):: xm_s
integer(spi),         intent(out):: ipan
!------------------------------------------------------------------------------
real(dp),dimension(3):: xc
real(dp),dimension(2):: xm
!==============================================================================
xc=xc_s
call xctoxm(xc,xm,ipan)
xm_s=xm
end subroutine xctoxm_s
!==============================================================================
subroutine xctoxmd_s(xc_s,xm_s,dxmdxc_s,ipan)!  [xctoxm]
!==============================================================================
implicit none
real(sp),dimension(3),  intent(in ):: xc_s
real(sp),dimension(2),  intent(out):: xm_s
real(sp),dimension(2,3),intent(out):: dxmdxc_s
integer(spi),           intent(out):: ipan
!------------------------------------------------------------------------------
real(dp),dimension(3)  :: xc
real(dp),dimension(2)  :: xm
real(dp),dimension(2,3):: dxmdxc
!==============================================================================
xc=xc_s
call xctoxm(xc,xm,dxmdxc,ipan)
xm_s=xm
dxmdxc_s=dxmdxc
end subroutine xctoxmd_s

!==============================================================================
subroutine reorient(ipan,xm,rrot)!  [reorient]
!==============================================================================
! If xm lies inside the domain |xm(1)|<1 and |xm(2)<1, leave it
! unchanged
! and return the identity rotation in rrot. If xm lies just outside this
! map domain, switch ipan to the neighboring panel, return xm as the
! map coordinates relative to the origin at the center of this new
! panel,
! and return rrot as the relative rotation, d(xm_new)/d(xm_old).
!==============================================================================
use pkind, only: spi,dp
use pietc, only: u0,u1,mu1,u2,mu2
implicit none
integer(spi),           intent(inout):: ipan
real(dp),dimension(2),  intent(inout):: xm
real(dp),dimension(2,2),intent(  out):: rrot
!------------------------------------------------------------------------------
real(dp),dimension(2,2,0:3)   :: m
real(dp),dimension(2,0:3)     :: v
integer(spi), dimension(0:3,6):: jrot,pnei
integer(spi)                  :: j
data v/mu2,u0, u0,mu2, u2,u0, u0,u2/
data m/u1,u0,u0,u1,  u0,mu1,u1,u0,  mu1,u0,u0,mu1,  u0,u1,mu1,u0/
data jrot/0,0,0,0,  0,3,0,1,  0,2,0,2,  0,1,0,3,  1,2,3,0, 3,0,1,2/
data pnei/2,5,4,6,  3,5,1,6,  4,5,2,6,  1,5,3,6,  2,3,4,1,  2,1,4,3/
!==============================================================================
rrot=m(:,:,0)
if(xm(1)>u1)then
   j=jrot(0,ipan); ipan=pnei(0,ipan)
   rrot=matmul(rrot,m(:,:,j))
   xm=matmul(m(:,:,j),xm+v(:,0))
endif
if(xm(1)<mu1)then
   j=jrot(2,ipan); ipan=pnei(2,ipan)
   rrot=matmul(rrot,m(:,:,j))
   xm=matmul(m(:,:,j),xm+v(:,2))
endif
if(xm(2)>u1)then
   j=jrot(1,ipan); ipan=pnei(1,ipan)
   rrot=matmul(rrot,m(:,:,j))
   xm=matmul(m(:,:,j),xm+v(:,1))
endif
if(xm(2)<mu1)then
   j=jrot(3,ipan); ipan=pnei(3,ipan)
   rrot=matmul(rrot,m(:,:,j))
   xm=matmul(m(:,:,j),xm+v(:,3))
endif
end subroutine reorient

!==============================================================================
subroutine xmtoxm(ipan,xm1,jpan,xm2,dxm2dxm1)!  [xmtoxm]
!==============================================================================
! Given a point whose map coordinates on the cube are given by xm1
! relative
! to the frame centered at the middle of panel ipan, return the same in
! jpan, xm2 if |xm1(1)|<1 and |xm1(2)|<1, but if one of these
! coordinates
! exceeds its bounds, return jpan as the neighboring panel in that
! direction,
! xm2 the map coordinates relative to the frame of jpan, and, in all
! cases,
! return the rotation matrix dxm2dxm1 as the local jacobian matrix of
! the
! affine transformation, d(xm2)/d(xm1).
!==============================================================================
use pkind, only: spi,dp
implicit none
integer(spi),           intent(in ):: ipan
real(dp),dimension(2),  intent(in ):: xm1
integer(spi),           intent(out):: jpan
real(dp),dimension(2),  intent(out):: xm2
real(dp),dimension(2,2),intent(out):: dxm2dxm1
!==============================================================================
if(ipan<1 .or. ipan>6)stop 'In xmtoxm; ipan is an invalid panel index'
jpan=ipan
xm2=xm1
call reorient(jpan,xm2,dxm2dxm1)
end subroutine xmtoxm

end module cuco


