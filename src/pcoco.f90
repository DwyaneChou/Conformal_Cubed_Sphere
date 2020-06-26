!
!                                *********************************************
!                                *  R. J. Purser, NOAA/NCEP/EMC
!                                2014   *
!                                *  jim.purser@noaa.gov
!                                *
!                                *
!                                *
!                                *********************************************
!
! Utility routines for conformal coefficients of highly symmetrical
! polyhedrons mapped to the sphere. Each geometry is expanded about its
! vertex
! (series co1) and about its face-midpoint (series co2) (even though,
! except 
! in the case of the hexagonal dihedron, the first series is formally
! sufficient). In each series, the map is scaled so as to make the
! radius of
! convergence unity.
!
! Codes for the 7 griddable-regular polyhedra are as follows:
! Code;  Polyhedron;    Face symmetry index, F;  Vertex symmetry index,
! V;  V0
! dh: Triangular dihedron        3                           2
! 6
! th: Tetrahedron                3                           3
! 6
! oh: Octahedron                 3                           4
! 6
! ih: Icosahedron                3                           5
! 6
! sh: Square dihedron            4                           2
! 4
! cu: Cube                       4                           3
! 4
! hh: Hexagonal dihedron         6                           2
! 3
!
! Each series expansion for these geometries assumes a scaling in which
! the distance between adjacent vertices is one unit (excep in the
! special
! case of the triangular dihedron), and an orientation
! such that, with an origin at one vertex, the neighboring vertex in
! complex
! map units is to be found at z=c1=(1,0). The stereographic mapping of
! this
! neigbor, with magnitude referred to as "stereographic distance", will
! be
! denoted w1=(w1,0). The face in the positive half-plane with edge
! (0,c1)
! has center at complex ca, and stereographic distance wa; the negative
! half plane counterpart is the face at center cb=conjg(ca). 
! For the 2-series expansions (which we want eventually to be standard)
! the first series, CCco1, is the expansion about 0 of w**V in terms of
! z**V0.
! The second series, coefficients, CCco2, relates to the expansion of
! the
! steroegraphic map w' about wa, but in terms of powers of the map
! coordinate
! z' about ca scaled and rotated to make vertex z=0 the point z'=(1,0).
! In
! this case the function expanded is (w')**F and is in powers of
! (z')**F.
! For both series, the map scaling is designed to ensure that the radius
! of convergence is unity. The arc or "circles" method for
! boot-strapping
! these series uses arcs of radius (in |z| or in |z'|) of r0, slightly
! less
! than 1.
! Note: since extracting w involves applying fractional powers, it is
! crucial
! that the correct alias is selected.
! The construction process requires constantly transforming from one
! stereo-
! graphic framework to another, so we have adopted a standard way of
! doing
! this. If wa is the standard stereographic distance of the center of
! another
! stereographic projection that provides a mapped point w' in the map
! rotated to put z=0 on the positive axis [i.e., at w'=(wa,0)], then the
! mobius transformation, w = (wa-w')/(wa*w'+c1), gives the standard
! sterographic image in the framework oriented to also put the other
! map center on the positive axis, at w=(wa,0). Thus, ANY stereographic
! conversion can be done when this canonical mobius transformation is
! both
! preceded and followed by simple rotations. This is how we always do
! it.
!
! The case of the triangular dihedron is special because the radius of
! convergence of the first series expansion (co1) is restricted, not by
! the distance to the nearest vertex, but by the distance to the
! opposite
! edge-midpoint, which maps to stereographic infinity. As a result of
! this
! rescaling, and the implied smaller test arc for CO1, the alternative
! "sect1=3" series expanding for points on this arc is the CO2 expansion
! cantered on ca, not the point ze usd for the other triangular-face
! geometries.
!=============================================================================
! 
! DIRECT DEPENDENCIES:
! Modules:   pmat4, pfft1, pietc, pkind
! Libraries: pfft, pmat
!=============================================================================
module pcoco
!=============================================================================
use pkind, only: spi,sp,dp,dpc
implicit none
private
public:: set_dhco,set_thco,set_ohco,set_ihco,set_shco,set_cuco,set_hhco, &
         dhztoc,dhztoc0,thztoc,thztoc0,ohztoc,ohztoc0,ihztoc,ihztoc0, &
         shztoc,shztoc0,cuztoc,cuztoc0,hhztoc,hhztoc0
interface set_dhco; module procedure set_dhco2;           end interface
interface set_thco; module procedure set_thco2;           end interface
interface set_ohco; module procedure set_ohco1,set_ohco2; end interface
interface set_ihco; module procedure set_ihco2;           end interface
interface set_shco; module procedure set_shco2;           end interface
interface set_cuco; module procedure set_cuco2;           end interface
interface set_hhco; module procedure set_hhco2;           end interface
interface dhztoc;   module procedure dhztoc,  dhztoc_d;   end interface
interface dhztoc0;  module procedure dhztoc0, dhztoc0_d;  end interface
interface thztoc;   module procedure thztoc,  thztoc_d;   end interface
interface thztoc0;  module procedure thztoc0, thztoc0_d;  end interface
interface ohztoc;   module procedure ohztoc,  ohztoc_d;   end interface
interface ohztoc0;  module procedure ohztoc0, ohztoc0_d;  end interface
interface ihztoc;   module procedure ihztoc,  ihztoc_d;   end interface
interface ihztoc0;  module procedure ihztoc0, ihztoc0_d;  end interface
interface shztoc;   module procedure shztoc,  shztoc_d;   end interface
interface shztoc0;  module procedure shztoc0, shztoc0_d;  end interface
interface cuztoc;   module procedure cuztoc,  cuztoc_d;   end interface
interface cuztoc0;  module procedure cuztoc0, cuztoc0_d;  end interface
interface hhztoc;   module procedure hhztoc,  hhztoc_d;   end interface
interface hhztoc0;  module procedure hhztoc0, hhztoc0_d;  end interface

contains

!==============================================================================
subroutine set_ohco1(m,co)!  [set_ohco]
!==============================================================================
! Get the conformal octahedron expansion coefficients.
! m is the number of Taylor expansion coefficients
! co(1:m) is the array of real coefficients.
! If z is the complex map coordinate in the vertex-centered framework
! in which a neighboring vertex is at z=1, and Real(z)<= 1/2 and 
! |arg(z)| <= 30-degrees, then the complex stereographic mapping, w, of
! this
! z for the map centered at the vertex, z=0, oriented with z=1 mapping
! to w=1,
! is given by S(z**6)**(1/4), where S is the analytic function whose
! Taylor
! series coefficients are given by co.
!
! Note: this is an earlier style of derivation for a single expansion
! series.
! It is intended to be retired once all routines have been converted to
! the
! alternative version, set_ohco2, for the two expansion series.
!==============================================================================
use pietc,only: u0,o2,u1,r2,r3,pi,c1,ci,z060,z210,z225,z315
use pfft1 ,only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- no. of test pts in a 30-degree test arc
                        n2=n*2,&! <- no. of pts that complete a minimal cycle
                        nm=n-1,n2m=n2-1,nit=40
real(dp),parameter           :: pio6=pi/6,pio3=pi/3,dalpha=pio6/n, &! 
                                r0=.96_dp,&  ! radius of arc in map space
                                r06=r0**6,r03=r0**3,or12=o2/r3,r3p=u1+r3, &
                                r2or3p=r2/r3p
complex(dpc),parameter       :: ca=(o2,or12),wa=(r2or3p,u0)
complex(dpc),dimension(0:n)  :: zarg
complex(dpc),dimension(0:n2m):: ww
complex(dpc)                 :: w,z0,z1
integer(spi)                 :: i,it
integer(spi),dimension(0:n)  :: sect
!=============================================================================
do i=0,n
   z0=r0*exp(ci*i*dalpha); z1=c1-z0
   zarg(i)=z1**6
   if(aimag(z1)>=-o2*abs(z1))then;      sect(i)=1
   else;                                sect(i)=2
   endif
   if(i==0.or.i==n)zarg(i)=cmplx(real(zarg(i)),u0,dpc)
enddo

co=u0 ! <- initialial guess can be anything finite -- e.g., zero
do it=1,nit
! Recompute the w of each test point of the 30-degree arc:
   do i=0,n
      call tay(m,co,zarg(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
      w=sqrt(sqrt(w))
      select case(sect(i))
      case(1)
         w=(c1-w)/(c1+w)
      case(2)
         w=(ci-w)/(ci+w)
      end select
      ww(i)=w**4
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r06,co)
enddo
end subroutine set_ohco1

!==============================================================================
subroutine set_dhco2(m,co1,co2)!  [set_dhco]
!==============================================================================
! Triangular dihedron, {3,2}
! In this geometry, the nearest singularity to a given vertex in map
! space
! is the midpoint of the opposite edge, which must map to infinity! This
! is a situation unique to this triangular dihedron and necessitates an
! alternative scaling for the first series of coefficients to put this
! edge
! midpoint a unit distance away. Orientation is not changed though.
!=============================================================================
use pietc, only: u0,o2,u1,pi,s60,r3,c1,ci,                    &
                 p300=>z180,p030=>z090,p330=>z270,p270=>z090, &
                 z030,z060,z120,z210
use pfft1, only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co1,co2
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- number of test points in each test arc
                    n2=n*2,&! <- number of points that complete a minimal cycle
                    nm=n-1,n2m=n2-1,nit=400,& ! <- (Larger nit needed)
                    f=3,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter           :: of=u1/f,ov=u1/v,ov0=u1/v0, &
                                dalpha=pi*ov0/n,dbeta=pi*of/n, &
                                r0=.9,r0v0=r0**v0,r0f=r0**f, &
                                re=u1/r3,o2re=o2*re,mo2re=-o2re
complex(dpc),parameter       :: ca=(o2,o2re),cb=(o2,mo2re),&
                                w1=(r3,u0),wa=(u1,u0)!=c1 ?
complex(dpc),dimension(0:n)  :: zarg1,zarg2
complex(dpc),dimension(0:n2m):: ww
complex(dpc)                 :: w,z0,z1,za,zb
integer(spi)                 :: i,it
integer(spi),dimension(0:n)  :: sect1,sect2
real(sp),dimension(9)        :: tco1,tco2
data tco1/&! (dh co1 seed)
 0.9619036372021E+01_sp,-0.1762397347166E+02_sp, 0.2577038145238E+02_sp, &
-0.3390780117146E+02_sp, 0.4204570080061E+02_sp,-0.5018357716593E+02_sp, &
 0.5832145459619E+02_sp,-0.6645933197943E+02_sp, 0.7459720936461E+02_sp/
data tco2/&! (dh co2 seed)
 0.5513701576720E+01_sp,-0.1520045253854E+02_sp, 0.2993241397445E+02_sp, &
-0.4951151943770E+02_sp, 0.7399776178162E+02_sp,-0.1033787776542E+03_sp, &
 0.1376570774942E+03_sp,-0.1768322023174E+03_sp, 0.2209042299788E+03_sp/
!=============================================================================
co1=0; co2=0; i=min(9,m); co1(1:i)=tco1(1:i); co2(1:i)=tco2(1:i);
do i=0,n
   z0=r0*s60*exp(ci*i*dalpha); z1=(c1-z0)/s60; za=z030*r3*(z0-ca)
   if(abs(z1)<=abs(za))then; zarg1(i)=z1**v0
      if(aimag(z1)>=-o2*abs(z1))then;      sect1(i)=1
      else;                                  sect1(i)=2
      endif
   else;                     zarg1(i)=za**f; sect1(i)=3
   endif
   if(i==0.or.i==n)zarg1(i)=cmplx(real(zarg1(i)),u0,dpc)
enddo
do i=0,n
   z0=ca+z210*re*r0*exp(ci*i*dbeta); zb=z210*r3*(z0-cb); z0=z0/s60
   if(abs(z0)<=abs(zb))then; zarg2(i)=z0**v0
      if(aimag(z0)>=-o2*abs(z0))then;      sect2(i)=1
      else;                                sect2(i)=2
      endif
   else;                     zarg2(i)=zb**f; sect2(i)=3
   endif
   if(i==0.or.i==n)zarg2(i)=cmplx(real(zarg2(i)),u0,dpc)
enddo

do it=1,nit
   do i=0,n
      select case(sect1(i))
      case(1)
         call tay(m,co1,zarg1(i),w)
         if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=(w1-w)/(w1*w+c1)
      case(2)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=p300*w; w=(w1-w)/(w1*w+c1)
      case(3)
         call tay(m,co2,zarg1(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of
         w=z120*w; w=(wa-w)/(wa*w+c1); w=p030*w
      end select
      ww(i)=w**v
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0v0,co1)
   do i=0,n
      select case(sect2(i))
      case(1)
         call tay(m,co1,zarg2(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=p330*w; w=(wa-w)/(wa*w+c1)
      case(2)
         call tay(m,co1,zarg2(i),w)
         w=w**ov
         w=p270*w; w=(wa-w)/(wa*w+c1)
      case(3)
         call tay(m,co2,zarg2(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z060*w; w=c1/w; w=z060*w
      end select
      ww(i)=w**f
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0f,co2)
enddo
end subroutine set_dhco2

!==============================================================================
subroutine set_thco2(m,co1,co2)!  [set_thco]
!==============================================================================
! Tetrahedron, {3,3}
!==============================================================================
use pietc, only: u0,o2,u1,pi,r2,r3,r5,or2,phi,c1,ci, &
                 p300=>z240,p030=>z060,p330=>z300,p270=>z180,  z060,z210
use pfft1 , only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co1,co2
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- number of test points in each test arc
                    n2=n*2,&! <- number of points that complete a minimal cycle
                    nm=n-1,n2m=n2-1,nit=40,&
                    f=3,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter           :: of=u1/f,ov=u1/v,ov0=u1/v0,     &
                                dalpha=pi*ov0/n,dbeta=pi*of/n, &
                                r0=.98,r0v0=r0**v0,r0f=r0**f,  &
                                re=u1/r3,o2re=o2*re,mo2re=-o2re
complex(dpc),parameter       :: ca=(o2,o2re),cb=(o2,mo2re),ce=(u1,re), &
                                w1=(r2,u0),wa=(or2,u0),wg=(r2,u0)
complex(dpc),dimension(0:n)  :: zarg1,zarg2
complex(dpc),dimension(0:n2m):: ww
complex(dpc)                 :: w,z0,z1,zb,ze
integer(spi)                 :: i,it
integer(spi),dimension(0:n)  :: sect1,sect2
real(sp),dimension(9)        :: tco1,tco2
data tco1/&! (th co1 seed)
 0.9068914462853E+01_sp,-0.1246203115328E+02_sp, 0.1097735631252E+02_sp, &
-0.7921574561043E+01_sp, 0.5100275350975E+01_sp,-0.3049657461188E+01_sp, &
 0.1731293360460E+01_sp,-0.9458002097187E+00_sp, 0.5016118969401E+00_sp/
data tco2/&! (th co2 seed)
 0.9746939435828E+00_sp,-0.1007657162539E+01_sp, 0.5456708012914E+00_sp, &
-0.2108345462604E+00_sp, 0.6603031327278E-01_sp,-0.1791356423097E-01_sp, &
 0.4377374424631E-02_sp,-0.9878289452680E-03_sp, 0.2093859642477E-03_sp/
!=============================================================================
co1=0; co2=0; i=min(9,m); co1(1:i)=tco1(1:i); co2(1:i)=tco2(1:i);
do i=0,n
   z0=r0*exp(ci*i*dalpha); z1=c1-z0; ze=ci*r3*(z0-ce)
   if(abs(z1)<=abs(ze))then; zarg1(i)=z1**v0
      if(aimag(z1)>=-o2*abs(z1))then;      sect1(i)=1
      else;                                sect1(i)=2
      endif
   else;                   zarg1(i)=ze**f; sect1(i)=3
   endif
   if(i==0.or.i==n)zarg1(i)=cmplx(real(zarg1(i)),u0,dpc)
enddo
do i=0,n
   z0=ca+z210*re*r0*exp(ci*i*dbeta); zb=z210*r3*(z0-cb)
   if(abs(z0)<=abs(zb))then; zarg2(i)=z0**v0
      if(aimag(z0)>=-o2*abs(z0))then;      sect2(i)=1
      else;                                sect2(i)=2
      endif
   else;                   zarg2(i)=zb**f; sect2(i)=3
   endif
   if(i==0.or.i==n)zarg2(i)=cmplx(real(zarg2(i)),u0,dpc)
enddo

do it=1,nit
   do i=0,n
      select case(sect1(i))
      case(1)
         call tay(m,co1,zarg1(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=(w1-w)/(w1*w+c1)
      case(2)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=p300*w; w=(w1-w)/(w1*w+c1)
      case(3)
         call tay(m,co2,zarg1(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z060*w; w=c1/w; w=p030*w
      end select
      ww(i)=w**v
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0v0,co1)
   do i=0,n
      select case(sect2(i))
      case(1)
         call tay(m,co1,zarg2(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=p330*w; w=(wa-w)/(wa*w+c1)
      case(2)
         call tay(m,co1,zarg2(i),w)
         w=w**ov
         w=p270*w; w=(wa-w)/(wa*w+c1)
      case(3)
         call tay(m,co2,zarg2(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z060*w; w=(wg-w)/(wg*w+c1); w=z060*w
      end select
      ww(i)=w**f
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0f,co2)
enddo
end subroutine set_thco2

!==============================================================================
subroutine set_ohco2(m,co1,co2)!  [set_ohco]
!==============================================================================
! Octahedron, {3,4}
!==============================================================================
use pietc, only: u0,o2,u1,pi,r2,r3,or2,c1,ci, &
                 p300=>z270,p030=>z045,p330=>z315,p270=>z225, z060,z210
use pfft1 , only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co1,co2
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- number of test points in each test arc
                    n2=n*2,&! <- number of points that complete a minimal cycle
                    nm=n-1,n2m=n2-1,nit=40,&
                    f=3,v=4,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter           :: of=u1/f,ov=u1/v,ov0=u1/v0, &
                                dalpha=pi*ov0/n,dbeta=pi*of/n, &
                                r0=.96,r0v0=r0**v0,r0f=r0**f, &
                                re=u1/r3,o2re=o2*re,mo2re=-o2re, &
                                r3p=r3+u1,r2or3p=r2/r3p,&
                                r3por2=r3p/r2
complex(dpc),parameter       :: ca=(o2,o2re),cb=(o2,mo2re),ce=(u1,re), &
                                w1=c1,wa=(r2or3p,u0),we=(r3por2,u0),wg=(or2,u0)
complex(dpc),dimension(0:n)  :: zarg1,zarg2
complex(dpc),dimension(0:n2m):: ww
complex(dpc)                 :: w,z0,z1,zb,ze
integer(spi)                 :: i,it
integer(spi),dimension(0:n)  :: sect1,sect2
real(sp),dimension(9)        :: tco1,tco2
data tco1/&! (oh co1 seed)
 0.1900056567314E+01_sp,-0.1031489988284E+01_sp, 0.1830665854289E+00_sp, &
-0.3476830258531E-01_sp,-0.2360311571181E-02_sp,-0.3406635872892E-02_sp, &
-0.2143436458638E-02_sp,-0.1545115435207E-02_sp,-0.1147317716686E-02_sp/
data tco2/&! (oh co2 seed)
 0.2436734858957E+00_sp,-0.1049642877645E+00_sp, 0.4392236563425E-02_sp, &
-0.1135194453049E-02_sp,-0.9274554954815E-03_sp,-0.5266686615558E-03_sp, &
-0.3582486325748E-03_sp,-0.2544778558661E-03_sp,-0.1881826612429E-03_sp/
!=============================================================================
! These coefficients, ohco1, look the same as hhco2 !
co1=0; co2=0; i=min(9,m); co1(1:i)=tco1(1:i); co2(1:i)=tco2(1:i);
do i=0,n
   z0=r0*exp(ci*i*dalpha); z1=c1-z0; ze=ci*r3*(z0-ce)
   if(abs(z1)<=abs(ze))then; zarg1(i)=z1**v0
      if(aimag(z1)>=-o2*abs(z1))then;      sect1(i)=1
      else;                                sect1(i)=2
      endif
   else;                   zarg1(i)=ze**f; sect1(i)=3
   endif
   if(i==0.or.i==n)zarg1(i)=cmplx(real(zarg1(i)),u0,dpc)
enddo
do i=0,n
   z0=ca+z210*re*r0*exp(ci*i*dbeta); zb=z210*r3*(z0-cb)
   if(abs(z0/zb)<u1)then;   zarg2(i)=z0**v0
      if(aimag(z0)>=-o2*abs(z0))then;      sect2(i)=1
      else;                                sect2(i)=2
      endif
   else;                   zarg2(i)=zb**f; sect2(i)=3
   endif
   if(i==0.or.i==n)zarg2(i)=cmplx(real(zarg2(i)),u0,dpc)
enddo

do it=1,nit
   do i=0,n
      select case(sect1(i))
      case(1)
         call tay(m,co1,zarg1(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=(w1-w)/(w1*w+c1)
      case(2)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=p300*w; w=(w1-w)/(w1*w+c1)
      case(3)
         call tay(m,co2,zarg1(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z060*w; w=(we-w)/(we*w+c1); w=p030*w
      end select
      ww(i)=w**v
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0v0,co1)
   do i=0,n
      select case(sect2(i))
      case(1)
         call tay(m,co1,zarg2(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=p330*w; w=(wa-w)/(wa*w+c1)
      case(2)
         call tay(m,co1,zarg2(i),w)
         w=w**ov
         w=p270*w; w=(wa-w)/(wa*w+c1)
      case(3)
         call tay(m,co2,zarg2(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z060*w; w=(wg-w)/(wg*w+c1); w=z060*w
      end select
      ww(i)=w**f
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0f,co2)
enddo
end subroutine set_ohco2

!==============================================================================
subroutine set_ihco2(m,co1,co2)!  [set_ihco]
!==============================================================================
! Icosahedron, {3,5}
!==============================================================================
use pietc, only: u0,o2,u1,u2,u3,pi,r2,r3,r5,or2,phi,c1,ci, &
                 p300=>z288,p030=>z036,p330=>z324,p270=>z252, z060,z210
use pfft1 , only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co1,co2
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- number of test points in each test arc
                    n2=n*2,&! <- number of points that complete a minimal cycle
                    nm=n-1,n2m=n2-1,nit=40,&
                    f=3,v=5,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter           :: of=u1/f,ov=u1/v,ov0=u1/v0,       &
                                dalpha=pi*ov0/n,dbeta=pi*of/n,   &
                                r0=.96,r0v0=r0**v0,r0f=r0**f,    &
                                rwe=(3+r5)/(r3*sqrt(5+2*r5)+u1), &
                                rwa=4/(r2*r3*sqrt(5+r5)+3+r5),   &
                                re=u1/r3,o2re=o2*re,mo2re=-o2re,ophi=u1/phi, &
                                r5p3=r5+u3,u2or5p3=u2/r5p3
complex(dpc),parameter       :: ca=(o2,o2re),cb=(o2,mo2re), &
                                ce=(u1,re),w1=(ophi,u0),wa=(rwa,u0), &
                                we=(rwe,u0),wg=(u2or5p3,u0)
complex(dpc),dimension(0:n)  :: zarg1,zarg2
complex(dpc),dimension(0:n2m):: ww
complex(dpc)                 :: w,z0,z1,zb,ze
integer(spi)                 :: i,it
integer(spi),dimension(0:n)  :: sect1,sect2
real(sp),dimension(9)        :: tco1,tco2
data tco1/&! (ih co1 seed)
 0.1187535354571E+00_sp,-0.2216091771703E-01_sp,-0.2673752260185E-02_sp, &
-0.1155838847443E-02_sp,-0.6346935826460E-03_sp,-0.3969422965766E-03_sp, &
-0.2701473161594E-03_sp,-0.1949235874619E-03_sp,-0.1468112252273E-03_sp/
data tco2/&! (ih co2 seed)
 0.4931604613269E-01_sp,-0.7477633936271E-02_sp,-0.1416525829555E-02_sp, &
-0.5340151024262E-03_sp,-0.2925461691417E-03_sp,-0.1817895933006E-03_sp, &
-0.1233701250510E-03_sp,-0.8884961066928E-04_sp,-0.6683174099670E-04_sp/
!=============================================================================
co1=0; co2=0; i=min(9,m); co1(1:i)=tco1(1:i); co2(1:i)=tco2(1:i);
do i=0,n
   z0=r0*exp(ci*i*dalpha); z1=c1-z0; ze=ci*r3*(z0-ce)
   if(abs(z1)<=abs(ze))then; zarg1(i)=z1**v0
      if(aimag(z1)>=-o2*abs(z1))then;      sect1(i)=1
      else;                                sect1(i)=2
      endif
   else;                   zarg1(i)=ze**f; sect1(i)=3
   endif
   if(i==0.or.i==n)zarg1(i)=cmplx(real(zarg1(i)),u0,dpc)
enddo
do i=0,n
   z0=ca+z210*re*r0*exp(ci*i*dbeta); zb=z210*r3*(z0-cb)
   if(abs(z0/zb)<u1)then;   zarg2(i)=z0**v0
      if(aimag(z0)>=-o2*abs(z0))then;      sect2(i)=1
      else;                                sect2(i)=2
      endif
   else;                   zarg2(i)=zb**f; sect2(i)=3
   endif
   if(i==0.or.i==n)zarg2(i)=cmplx(real(zarg2(i)),u0,dpc)
enddo

do it=1,nit
   do i=0,n
      select case(sect1(i))
      case(1)
         call tay(m,co1,zarg1(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=(w1-w)/(w1*w+c1)
      case(2)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=p300*w; w=(w1-w)/(w1*w+c1)
      case(3)
         call tay(m,co2,zarg1(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z060*w; w=(we-w)/(we*w+c1); w=p030*w
      end select
      ww(i)=w**v
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0v0,co1)
   do i=0,n
      select case(sect2(i))
      case(1)
         call tay(m,co1,zarg2(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
         w=w**ov
         w=p330*w; w=(wa-w)/(wa*w+c1)
      case(2)
         call tay(m,co1,zarg2(i),w)
         w=w**ov
         w=p270*w; w=(wa-w)/(wa*w+c1)
      case(3)
         call tay(m,co2,zarg2(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z060*w; w=(wg-w)/(wg*w+c1); w=z060*w
      end select
      ww(i)=w**f
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0f,co2)
enddo
end subroutine set_ihco2

!==============================================================================
subroutine set_shco2(m,co1,co2)!  [set_shco]
!==============================================================================
! Square dihedron, {4,2}
!==============================================================================
use pietc, only: u0,o2,u1,mu1, r2,or2,r3,pi,c1,ci,           &
                 p270=>z180,p045=>z090,p315=>z270, z045,z225
use pfft1 , only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co1,co2
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- number of test points in each test arc
                    n2=n*2,&! <- number of points that complete a minimal cycle
                    nm=n-1,n2m=n2-1,nit=40, &
                    f=4,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter           :: of=u1/f,ov=u1/v,ov0=u1/v0,         &
                                dalpha=pi*ov0/n,dbeta=pi*of/n,     &
                                r0=.96,r0v0=r0**v0,r0f=r0**f,mo2=-o2
complex(dpc),parameter       :: ca=(o2,o2), cb=(o2,mo2), &
                                wa=c1,w1=c1
complex(dpc),dimension(0:n)  :: zarg1,zarg2
complex(dpc),dimension(0:n2m):: ww 
complex(dpc)                 :: z0,z1,za,zb,w
integer(spi),dimension(0:n)  :: sect1,sect2
integer(spi)                 :: i,it
real(dp),        dimension(9):: tco1,tco2
data tco1/&! (sh co1 seed)
 0.2954261252019E+01_sp,-0.3491063818073E+01_sp, 0.2406486732081E+01_sp, &
-0.1281252794503E+01_sp, 0.5854228999512E+00_sp,-0.2417224704995E+00_sp, &
 0.9286688073166E-01_sp,-0.3380499078942E-01_sp, 0.1180075671599E-01_sp/
data tco2/&! (sh co2 seed)
 0.2954261252019E+01_sp,-0.3491063818073E+01_sp, 0.2406486732081E+01_sp, &
-0.1281252794503E+01_sp, 0.5854228999512E+00_sp,-0.2417224704995E+00_sp, &
 0.9286688073166E-01_sp,-0.3380499078941E-01_sp, 0.1180075671599E-01_sp/
!=============================================================================
co1=0; co2=0; i=min(9,m); co1(1:i)=tco1(1:i); co2(1:i)=tco2(1:i);
do i=0,n
   z0=r0*exp(ci*i*dalpha); z1=c1-z0; za=z045*r2*(z0-ca)
   if(abs(z1)<=abs(za))then; zarg1(i)=z1**v0
      if(aimag(z1)>=-real(z1))then;           sect1(i)=1
      else;                                   sect1(i)=2
      endif
   else;                     zarg1(i)=za**f      
      if(aimag(z0)<o2)then;                   sect1(i)=3
      else;                                   sect1(i)=4
      endif
   endif
   if(i==0.or.i==n)zarg1(i)=cmplx(real(zarg1(i)),u0,dpc)
enddo
do i=0,n
   z0=ca+z225*or2*r0*exp(ci*i*dbeta)
   zb=z225*r2*(z0-cb)
   if(abs(z0/zb)<u1)then;    zarg2(i)=z0**v0; sect2(i)=1
   else;                     zarg2(i)=zb**f;  sect2(i)=2
   endif
   if(i==0.or.i==n)zarg2(i)=cmplx(real(zarg2(i)),u0,dpc)
enddo

do it=1,nit
! Recompute the w of each test point of each 45-degree arc:
   do i=0,n
      select case(sect1(i))
      case(1)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=(w1-w)/(w1*w+c1)
      case(2)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=p270*w;   w=(w1-w)/(w1*w+c1)
      case(3)
         call tay(m,co2,zarg1(i),w)
         w=w**of
         w=ci*w; w=(wa-w)/(wa*w+c1); w=p045*w
      case(4)
         call tay(m,co2,zarg1(i),w)
         w=w**of
         w=-w;   w=(wa-w)/(wa*w+c1); w=p045*w
      end select
      ww(i)=w**v
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0v0,co1)
   do i=0,n
      select case(sect2(i))
      case(1)
         call tay(m,co1,zarg2(i),w)
         w=w**ov
         w=p315*w;  w=(wa-w)/(wa*w+c1)
      case(2)
         call tay(m,co2,zarg2(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=z045*w;  w=c1/w;  w=z045*w
      end select
      ww(i)=w**f
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0f,co2)
enddo
end subroutine set_shco2

!==============================================================================
subroutine set_cuco2(m,co1,co2)!  [set_cuco]
!==============================================================================
! Cube, {4,3}
!=============================================================================
use pietc, only: u0,o2,u1, r2,or2,r3,pi,c1,ci,                &
                 p270=>z240,p045=>z060,p315=>z300, z045,z225
use pfft1 , only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co1,co2
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- number of test points in each test arc
                    n2=n*2,&! <- number of points that complete a minimal cycle
                    nm=n-1,n2m=n2-1,nit=40, &
                    f=4,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter           :: of=u1/f,ov=u1/v,ov0=u1/v0, &
                                dalpha=pi*ov0/n,dbeta=pi*of/n, &
                                r0=.96,r0v0=r0**v0,r0f=r0**f, &
                                 r3p=r3+u1,r2or3p=r2/r3p,mo2=-o2
complex(dpc),parameter       :: ca=(o2,o2),cb=(o2,mo2), wa=(r2or3p,u0), &
                                w1=(or2,u0)
complex(dpc),dimension(0:n)  :: zarg1,zarg2
complex(dpc),dimension(0:n2m):: ww 
complex(dpc)                 :: z0,z1,za,zb,w
integer(spi),dimension(0:n)  :: sect1,sect2
integer(spi)                 :: i,it
real(sp),    dimension(9)    :: tco1,tco2
data tco1/&! (cu co1 seed)
 0.5222445411749E+00_sp,-0.1349990960563E+00_sp,-0.1970373552294E-01_sp, &
-0.3167426867673E-02_sp,-0.2797723788949E-02_sp,-0.1720480734504E-02_sp, &
-0.1164080730235E-02_sp,-0.8325527861928E-03_sp,-0.6217962129431E-03_sp/
data tco2/&! (cu co2 seed)
 0.1094170834081E+00_sp,-0.3352187479631E-01_sp,-0.9955635663016E-03_sp, &
-0.9915868353376E-03_sp,-0.5372057577482E-03_sp,-0.3367620231221E-03_sp, &
-0.2279393461075E-03_sp,-0.1632484731336E-03_sp,-0.1219559186293E-03_sp/
!=============================================================================
co1=0; co2=0; i=min(9,m); co1(1:i)=tco1(1:i); co2(1:i)=tco2(1:i);
do i=0,n
   z0=r0*exp(ci*i*dalpha); z1=c1-z0; za=z045*r2*(z0-ca)
   if(abs(z1)<=abs(za))then; zarg1(i)=z1**v0
      if(aimag(z1)>=-real(z1))then;           sect1(i)=1
      else;                                   sect1(i)=2
      endif
   else;                     zarg1(i)=za**f
      if(aimag(z0)<o2)then;                 sect1(i)=3
      else;                                 sect1(i)=4
      endif
   endif
   if(i==0.or.i==n)zarg1(i)=cmplx(real(zarg1(i)),u0,dpc)
enddo
do i=0,n
   z0=ca+z225*or2*r0*exp(ci*i*dbeta); zb=z225*r2*(z0-cb)
   if(abs(z0)<=abs(zb))then; zarg2(i)=z0**v0; sect2(i)=1
   else;                     zarg2(i)=zb**f;  sect2(i)=2
   endif
   if(i==0.or.i==n)zarg2(i)=cmplx(real(zarg2(i)),u0,dpc)
enddo

do it=1,nit
! Recompute the w of each test point of each 45-degree arc:
   do i=0,n
      select case(sect1(i))
      case(1)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=(w1-w)/(w1*w+c1)
      case(2)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=p270*w; w=(w1-w)/(w1*w+c1)
      case(3)
         call tay(m,co2,zarg1(i),w)
         w=w**of
         w=ci*w; w=(wa-w)/(wa*w+c1); w=p045*w
      case(4)
         call tay(m,co2,zarg1(i),w)
         w=w**of
         w=-w;   w=(wa-w)/(wa*w+c1); w=p045*w
      end select
      ww(i)=w**v
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0v0,co1)
   do i=0,n
      select case(sect2(i))
      case(1)
         call tay(m,co1,zarg2(i),w)
         w=w**ov
         w=p315*w; w=(wa-w)/(wa*w+c1)
      case(2)
         call tay(m,co2,zarg2(i),w); if(i==n)w=cmplx(real(w),u0,dpc)
         w=w**of;                    if(i==n)w=conjg(w)
         w=w*z045; w=(c1-w)/(w+c1); w=w*z045
      end select
      ww(i)=w**f
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0f,co2)
enddo
end subroutine set_cuco2

!==============================================================================
subroutine set_hhco2(m,co1,co2)!  [set_hhco]
!==============================================================================
! Hexagonal dihedron, {6,2}
!==============================================================================
use pietc, only: u0,o2,u1,pi,or3,or2,c1,ci, &
                 z060,z120,z240,z270,z300
use pfft1 ,only: tay,tay_reset
implicit none
integer(spi),         intent(in ):: m
real(dp),dimension(m),intent(out):: co1,co2
!-----------------------------------------------------------------------------
integer(spi),parameter:: n=120,& ! <- number of test points in each test arc
                    n2=n*2,&! <- number of points that complete a minimal cycle
                    nm=n-1,n2m=n2-1,nit=40, &
                    f=6,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter           :: of=u1/f,ov=u1/v,ov0=u1/v0,     &
                                dalpha=pi*ov0/n,dbeta=pi*of/n, &
                                r0=.98_dp,r0v0=r0**v0,r0f=r0**f,  &
                                re=or3
complex(dpc),parameter       :: ca=z060,cb=z300,w1=(re,u0)
complex(dpc),dimension(0:n)  :: zarg1,zarg2
complex(dpc),dimension(0:n2m):: ww
complex(dpc)                 :: w,z0,z1,za
integer(spi)                 :: i,it
integer(spi),dimension(0:n)  :: sect1
real(sp),dimension(9)        :: tco1,tco2
data tco1/&! (hh co1 seed)
 0.6126335085233E+00_sp,-0.3127665131381E+00_sp, 0.5565850487536E-01_sp, &
-0.1704913256054E-01_sp, 0.6242795394712E-03_sp,-0.1666855938577E-02_sp, &
-0.7516285486059E-03_sp,-0.5972808404004E-03_sp,-0.4307685305614E-03_sp/
data tco2/&! (hh co2 seed)
 0.1900056567314E+01_sp,-0.1031489988284E+01_sp, 0.1830665854290E+00_sp, &
-0.3476830258538E-01_sp,-0.2360311571118E-02_sp,-0.3406635872948E-02_sp, &
-0.2143436458586E-02_sp,-0.1545115435256E-02_sp,-0.1147317716638E-02_sp/
!=============================================================================
co1=u0; co2=u0; i=min(9,m); co1(1:i)=tco1(1:i); co2(1:i)=tco2(1:i);
do i=0,n
   z0=r0*exp(ci*i*dalpha); z1=c1-z0; za=z0-ca
   if(aimag(z0)<=o2*abs(z0))then;   zarg1(i)=z1**v0
      if(real(z1)>=o2*abs(z1))then; sect1(i)=1
      else;                         sect1(i)=2
      endif
   else;                            zarg1(i)=za**f
      if(real(z0)>o2)then;          sect1(i)=3
      else;                         sect1(i)=4
      endif
   endif
   if(i==0.or.i==n)zarg1(i)=cmplx(real(zarg1(i)),u0,dpc)
enddo
do i=0,n
   z0=ca+z240*r0*exp(ci*i*dbeta)
                                      zarg2(i)=z0**v0
   if(i==0)zarg2(i)=cmplx(real(zarg2(i)),u0,dpc)
enddo
do it=1,nit
   do i=0,n
      select case(sect1(i))
      case(1)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=(w1-w)/(w1*w+c1)
      case(2)
         call tay(m,co1,zarg1(i),w)
         w=w**ov
         w=-w;  w=(w1-w)/(w1*w+c1)
      case(3)
         call tay(m,co2,zarg1(i),w)
         w=w**of
         w=z060*w; w=(c1-w)/(w+c1); w=ci*w
      case(4)
         call tay(m,co2,zarg1(i),w)
         w=w**of
         w=(c1-w)/(w+c1); w=ci*w
      end select
      ww(i)=w**v
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0v0,co1)
   do i=0,n
      call tay(m,co1,zarg2(i),w); if(i==0)w=cmplx(real(w),u0,dpc)
      w=w**ov
      w=-ci*w; w=(c1-w)/(w+c1)
      ww(i)=w**f
      if(i==0.or.i==n)then; ww(   i)=cmplx(real(ww(i)),u0,dpc)
      else;                 ww(n2-i)=conjg(ww(i))
      endif
   enddo
   call tay_reset(n2,m,ww,r0f,co2)
enddo
end subroutine set_hhco2

!-----

!=============================================================================
subroutine dhztoc(m,co1,co2,z,xc)!  [dhztoc]
!=============================================================================
! Conformal triangular dihedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pkind, only: dp,dpc
use pietc, only: FF=>F,u0,o2,u1,r3,or3,s60,pi2,c0,c1,ci, &
                 p030=>z090,p330=>z270, z030,z090,z150,z210,z270,z300,z330
use pfft1 ,only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re
complex(dpc),parameter:: ca=(o2,o2re),wa=(u1,u0)
complex(dpc)          :: z0,za,zp,ww,w
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); za=z150*za
case(1); za=z030*za
case(2); za=z270*za
end select
z0=(ca+z210*za)/s60
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w)
      endif
   else; w=c0
   endif
else       ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)z0 =z300*z0      ! <-perform the rotation
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w=conjg(w)
      endif
   else; w=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w
   else;       w =p330*w
   endif
endif
select case(k)
case(0); w =z210*w
case(1); w =z330*w
case(2); w =z090*w
end select
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine dhztoc
!=============================================================================
subroutine dhztoc_d(m,co1,co2,z,xc,xcd)!  [dhztoc]
!=============================================================================
! Conformal triangular dihedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pkind, only: dp,dpc
use pietc, only: FF=>F,u0,o2,u1,r3,or3,s60,pi2,c0,c1,ci, &
                 p030=>z090,p330=>z270,z030,z090,z150,z210,z270,z300,z330
use pfft1 ,only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re
complex(dpc),parameter:: ca=(o2,o2re)
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); zad=z150; za=z150*za
case(1); zad=z030; za=z030*za
case(2); zad=z270; za=z270*za
end select
z0=(ca+z210*za)/s60; z0d=z210*zad/s60
zad=r3*zad;          za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
else       ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)then; z0d=z300*z0d; z0 =z300*z0      ! <-perform the rotation
   endif
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w=conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w=c0; wd=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w; wd=p030*wd
   else;       w =p330*w; wd=p330*wd
   endif
   wd=-wd*(c1*2)/(w+c1)**2; w=(c1-w)/(w+c1)
endif
select case(k)
case(0); w =z210*w; wd=z210*wd
case(1); w =z330*w; wd=z330*wd
case(2); w =z090*w; wd=z090*wd
end select
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine dhztoc_d
!=============================================================================
subroutine dhztoc0(m,co1,co2,z,xc)!  [dhztoc0]
!=============================================================================
! Conformal triangular dihedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r3,or3,s60,pi2,c0,c1,ci, &
                 p030=>z090,p330=>z270,         z150
use pfft1 ,only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re)
complex(dpc)          :: z0,za,zp,ww,w
logical               :: ka,kn1
!=============================================================================
if(real(z)>o2+eps)        stop 'x>1/2; Outside fundamental region'
if(aimag(z)>o2*abs(z)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z)<-eps)         stop 'y<0;  Outside fundamental region'
za=z-ca
za=z150*za
z0=z/s60
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); endif
   else; w=c0
   endif
! Transform from face-centered frame to vertex-centered frame:
   w=(c1-w)/(w+c1)
   w=p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); endif
   else; w =c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine dhztoc0
!=============================================================================
subroutine dhztoc0_d(m,co1,co2,z,xc,xcd)!  [dhztoc0]
!=============================================================================
! Conformal triangular dihedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r3,or3,s60,pi2,c0,c1,ci, &
                 p030=>z090,p330=>z270,         z150
use pfft1 ,only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re)
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
logical               :: ka,kn1
!=============================================================================
if(real(z)>o2+eps)        stop 'x>1/2; Outside fundamental region'
if(aimag(z)>o2*abs(z)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z)<-eps)         stop 'y<0;  Outside fundamental region'
za=z-ca
zad=z150; za=z150*za
z0=z/s60; z0d=c1/s60
zad=r3*zad; za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else
      w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
! Transform from face-centered frame to vertex-centered frame:
   wd=-wd*2/(w+c1)**2; w=(c1-w)/(w+c1)
   wd=p030*wd; w =p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine dhztoc0_d

!=============================================================================
subroutine thztoc(m,co1,co2,z,xc)!  [thztoc]
!=============================================================================
! Conformal tetrahedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z060,p330=>z300, z030,z090,z150,z210,z270,z300,z330
use pfft1 ,only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re
complex(dpc),parameter:: ca=(o2,o2re),wa=(or2,u0)
complex(dpc)          :: z0,za,zp,ww,w
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); za=z150*za
case(1); za=z030*za
case(2); za=z270*za
end select
z0=ca+z210*za
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w)
      endif
   else; w=c0
   endif
else       ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)z0 =z300*z0      ! <-perform the rotation
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w=conjg(w)
      endif
   else; w=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w
   else;       w =p330*w
   endif
   w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z210*w
case(1); w =z330*w
case(2); w =z090*w
end select
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine thztoc
!=============================================================================
subroutine thztoc_d(m,co1,co2,z,xc,xcd)!  [thztoc]
!=============================================================================
! Conformal tetrahedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z060,p330=>z300, z030,z090,z150,z210,z270,z300,z330
use pfft1 ,only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re
complex(dpc),parameter:: ca=(o2,o2re),wa=(or2,u0)
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); zad=z150; za=z150*za
case(1); zad=z030; za=z030*za
case(2); zad=z270; za=z270*za
end select
z0=ca+z210*za; z0d=z210*zad
zad=r3*zad;    za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
else       ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)then
      z0d=z300*z0d; z0 =z300*z0      ! <-perform the rotation
   endif
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w=conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w=c0; wd=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w; wd=p030*wd
   else;       w =p330*w; wd=p330*wd
   endif
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2; w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z210*w; wd=z210*wd
case(1); w =z330*w; wd=z330*wd
case(2); w =z090*w; wd=z090*wd
end select
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine thztoc_d
!=============================================================================
subroutine thztoc0(m,co1,co2,z0,xc)!  [thztoc0]
!=============================================================================
! Conformal tetrahedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z060,p330=>z300,    z150
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re),wa=(or2,u0)
complex(dpc)          :: za,zp,ww,w
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(aimag(z0)>o2*abs(z0)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca; za=z150*za
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w=c0
   endif
! Transform from face-centered frame to vertex-centered frame:
   w=(wa-w)/(wa*w+c1)
   w =p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w)
      endif
   else; w =c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine thztoc0
!=============================================================================
subroutine thztoc0_d(m,co1,co2,z0,xc,xcd)!  [thztoc0]
!=============================================================================
! Conformal tetrahedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z060,p330=>z300,    z150
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re),wa=(or2,u0)
complex(dpc)          :: za,zp,ww,zad,z0d,w,wd
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(aimag(z0)>o2*abs(z0)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca; zad=z150; za=z150*za
z0d=c1
zad=r3*zad; za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
! Transform from face-centered frame to vertex-centered frame:
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2; w=(wa-w)/(wa*w+c1)
   wd=p030*wd;                     w =p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine thztoc0_d

!=============================================================================
subroutine ohztoc(m,co1,co2,z,xc)!  [ohztoc]
!=============================================================================
! Conformal octahedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z045,p330=>z315, z030,z090,z150,z210,z270,z300,z330
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=4,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         r3p=r3+u1,r2or3p=r2/r3p
complex(dpc),parameter:: ca=(o2,o2re),wa=(r2or3p,u0)
complex(dpc)          :: z0,za,zp,ww,w
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); za=z150*za
case(1); za=z030*za
case(2); za=z270*za
end select
z0=ca+z210*za
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
else                      ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)z0 =z300*z0      ! <-perform the rotation
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w =c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w
   else;       w =p330*w
   endif
   w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z210*w
case(1); w =z330*w
case(2); w =z090*w
end select
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine ohztoc
!=============================================================================
subroutine ohztoc_d(m,co1,co2,z,xc,xcd)!  [ohztoc]
!=============================================================================
! Conformal octahedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z045,p330=>z315, z030,z090,z150,z210,z270,z300,z330
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=4,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         r3p=r3+u1,r2or3p=r2/r3p
complex(dpc),parameter:: ca=(o2,o2re),wa=(r2or3p,u0)
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); zad=z150; za=z150*za
case(1); zad=z030; za=z030*za
case(2); zad=z270; za=z270*za
end select
z0=ca+z210*za; z0d=z210*zad
zad=r3*zad;    za =r3*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
else                      ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)then; z0d=z300*z0d; z0 =z300*z0      ! <-perform the rotation
   endif
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov;  wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w; wd=p030*wd
   else;       w =p330*w; wd=p330*wd
   endif
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2
   w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z210*w; wd=z210*wd
case(1); w =z330*w; wd=z330*wd
case(2); w =z090*w; wd=z090*wd
end select
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine ohztoc_d
!=============================================================================
subroutine ohztoc0(m,co1,co2,z0,xc)!  [ohztoc0]
!=============================================================================
! Conformal octahedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z045,p330=>z315,   z150
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=4,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         r3p=r3+u1,r2or3p=r2/r3p,eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re), wa=(r2or3p,u0)
complex(dpc)          :: za,zp,ww,w
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(aimag(z0)>o2*abs(z0)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca
za=z150*za
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w=c0
   endif
! Transform from face-centered frame to vertex-centered frame:
   w=(wa-w)/(wa*w+c1)
   w =p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w =c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine ohztoc0
!=============================================================================
subroutine ohztoc0_d(m,co1,co2,z0,xc,xcd)!  [ohztoc0]
!=============================================================================
! Conformal octahedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,or3,pi2,c0,c1,ci, &
                 p030=>z045,p330=>z315,   z150
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=4,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,o2re=o2*re, &
                         r3p=r3+u1,r2or3p=r2/r3p,eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re),wa=(r2or3p,u0)
complex(dpc)          :: za,zp,ww,zad,z0d,w,wd
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(aimag(z0)>o2*abs(z0)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca
zad=z150; za=z150*za
z0d=c1
zad=r3*zad; za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
! Transform from face-centered frame to vertex-centered frame:
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2; w=(wa-w)/(wa*w+c1)
   wd=p030*wd; w =p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine ohztoc0_d

!=============================================================================
subroutine ihztoc(m,co1,co2,z,xc)!  [ihztoc]
!=============================================================================
! Conformal icosahedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,u3,u4,u5,r2,r3,r5,or3,pi2,c0,c1,ci, &
                 p030=>z036,p330=>z324, z030,z090,z150,z210,z270,z300,z330
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=5,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,     &
                         rwa=u4/(r2*r3*sqrt(u5+r5)+u3+r5),re=or3,o2re=o2*re
complex(dpc),parameter:: ca=(o2,o2re),wa=(rwa,u0)
complex(dpc)          :: z0,za,zp,ww,w
integer               :: k
logical(spi)          :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); za=z150*za
case(1); za=z030*za
case(2); za=z270*za
end select
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
else                      ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)z0 =z300*z0      ! <-perform the rotation
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w
   else;       w =p330*w
   endif
   w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z210*w
case(1); w =z330*w
case(2); w =z090*w
end select
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine ihztoc
!=============================================================================
subroutine ihztoc_d(m,co1,co2,z,xc,xcd)!  [ihztoc]
!=============================================================================
! Conformal icosahedron mapping to the unit sphere.
! Map from the unit triangular map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,u3,u4,u5,r2,r3,r5,or3,pi2,c0,c1,ci, &
                 p030=>z036,p330=>z324, z030,z090,z150,z210,z270,z300,z330
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=5,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,     &
                         rwa=u4/(r2*r3*sqrt(u5+r5)+u3+r5),re=or3,o2re=o2*re
complex(dpc),parameter:: ca=(o2,o2re),wa=(rwa,u0)
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); zad=z150; za=z150*za
case(1); zad=z030; za=z030*za
case(2); zad=z270; za=z270*za
end select
z0=ca+z210*za; z0d=z210*zad
zad=r3*zad;     za =r3*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
else                      ! Expand about the vertex
   kn=aimag(z0)>o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)then; z0d=z300*z0d; z0 =z300*z0      ! <-perform the rotation
   endif
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w=c0; wd=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p030*w; wd=p030*wd
   else;       w =p330*w; wd=p330*wd
   endif
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2;    w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z210*w; wd=z210*wd
case(1); w =z330*w; wd=z330*wd
case(2); w =z090*w; wd=z090*wd
end select
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine ihztoc_d
!=============================================================================
subroutine ihztoc0(m,co1,co2,z0,xc)!  [ihztoc0]
!=============================================================================
! Conformal icosahedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,u3,u4,u5,r2,r3,r5,or3,pi2,c0,c1,ci, &
                 p030=>z036,p330=>z324,   z150
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=5,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,     &
                         rwa=u4/(r2*r3*sqrt(u5+r5)+u3+r5),re=or3,o2re=o2*re, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re),wa=(rwa,u0)
complex(dpc)          :: za,zp,ww,w
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(aimag(z0)>o2*abs(z0)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca
za=z150*za
za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn1)w=conjg(w)
   else; w=c0
   endif
! Transform from face-centered frame to vertex-centered frame:
   w=(wa-w)/(wa*w+c1)
   w =p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn1)w=conjg(w)
   else; w=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine ihztoc0
!=============================================================================
subroutine ihztoc0_d(m,co1,co2,z0,xc,xcd)!  [ihztoc0]
!=============================================================================
! Conformal icosahedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,u3,u4,u5,r2,r3,r5,or3,pi2,c0,c1,ci, &
                 p030=>z036,p330=>z324,   z150
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=3,v=5,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,     &
                         rwa=u4/(r2*r3*sqrt(u5+r5)+u3+r5),re=or3,o2re=o2*re, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2re),wa=(rwa,u0)
complex(dpc)          :: za,zp,ww,zad,z0d,w,wd
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(aimag(z0)>o2*abs(z0)+eps)stop 'y>|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca
zad=z150;   za=z150*za
z0d=c1
zad=r3*zad; za =r3*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w=conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
! Transform from face-centered frame to vertex-centered frame:
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2; w=(wa-w)/(wa*w+c1)
   wd=p030*wd; w =p030*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w=conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w=c0; wd=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine ihztoc0_d

!=============================================================================
subroutine shztoc(m,co1,co2,z,xc)!  [shztoc]
!=============================================================================
! Conformal square-dihedron mapping to the unit sphere.
! Map from the unit square map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,c0,c1,ci, &
                 p045=>z090,p315=>z270,  z045,z135,z225,z315
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2)
complex(dpc)          :: z0,za,zp,ww,w
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0; if(real (za)>u0)k=1; if(aimag(za)>u0)k=k+2
select case(k)
case(1); za=-ci*za
case(2); za= ci*za
case(3); za=   -za
end select
z0=ca+za
za=z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative octant
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
else                    ! Expand about the vertex
   kn=aimag(z0)>real(z0)! negative octant of rotated z0 frame
   if(kn)z0 =-ci*z0     ! <-perform the rotation
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p045*w
   else;       w =p315*w
   endif
   w=(c1-w)/(w+c1)
endif
select case(k)
case(0); w =z225*w
case(1); w =z315*w
case(2); w =z135*w
case(3); w =z045*w
end select
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine shztoc
!=============================================================================
subroutine shztoc_d(m,co1,co2,z,xc,xcd)!  [shztoc]
!=============================================================================
! Conformal square-dihedron mapping to the unit sphere.
! Map from the unit square map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,c0,c1,ci, &
                 p045=>z090,p315=>z270,  z045,z135,z225,z315
use pfft1,  only: tay
use pmat4,  only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2)
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0; if(real (za)>u0)k=1; if(aimag(za)>u0)k=k+2
select case(k)
case(0); zad= c1
case(1); zad=-ci; za=-ci*za
case(2); zad= ci; za= ci*za
case(3); zad=-c1; za=   -za
end select
z0=ca+za;        z0d=zad
zad=z135*r2*zad; za =z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative octant
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
else                    ! Expand about the vertex
   kn=aimag(z0)>real(z0)! negative octant of rotated z0 frame
   if(kn)then; z0d=-ci*z0d; z0 =-ci*z0      ! <-perform the rotation
   endif
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w=c0; wd=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p045*w; wd=p045*wd
   else;       w =p315*w; wd=p315*wd
   endif
   wd=-wd*(c1*2)/(w+c1)**2; w=(c1-w)/(w+c1)
endif
select case(k)
case(0); w =z225*w; wd=z225*wd
case(1); w =z315*w; wd=z315*wd
case(2); w =z135*w; wd=z135*wd
case(3); w =z045*w; wd=z045*wd
end select
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine shztoc_d
!=============================================================================
subroutine shztoc0(m,co1,co2,z0,xc)!  [shztoc0]
!=============================================================================
! Conformal square-dihedron mapping to the unit sphere.
! Map from the square's fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,c0,c1,ci, p045=>z090, z135
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2)
complex(dpc)          :: za,zp,ww,w
logical               :: ka,kn1
!=============================================================================
if(aimag(z0)>real(z0)+eps)stop 'y>x; Outside fundamental region'
if(real(z0)>o2+eps)       stop 'x>1/2; Outside fundamental region'
if(aimag(z0)<-eps)        stop 'y<0; Outside fundamental region'
za=z0-ca
za =z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w =ww**of
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w=c0
   endif
! Transform to vertex-centered frame from face-centered frame
   w=(c1-w)/(w+c1)
   w =p045*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w =c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine shztoc0
!=============================================================================
subroutine shztoc0_d(m,co1,co2,z0,xc,xcd)!  [shztoc0]
!=============================================================================
! Conformal square-dihedron mapping to the unit sphere.
! Map from the square's fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,c0,c1,ci, p045=>z090, z135
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2)
complex(dpc)          :: za,zp,ww,zad,z0d,w,wd
logical               :: ka,kn1
!=============================================================================
if(aimag(z0)>real(z0)+eps)stop 'y>x; Outside fundamental region'
if(real(z0)>o2+eps)       stop 'x>1/2; Outside fundamental region'
if(aimag(z0)<-eps)        stop 'y<0; Outside fundamental region'
za=z0-ca
zad=c1
z0d=zad
zad=z135*r2*zad; za =z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w =ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
! Transform to vertex-centered frame from face-centered frame
   wd=-wd*(c1*2)/(w+c1)**2; w=(c1-w)/(w+c1)
   wd=p045*wd;                     w =p045*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine shztoc0_d

!=============================================================================
subroutine cuztoc(m,co1,co2,z,xc)!  [cuztoc]
!=============================================================================
! Conformal cubic mapping to the unit sphere.
! Map from the unit square map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,c0,c1,ci, &
                 p045=>z060,p315=>z300,  z045,z135,z225,z315
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,eps=1.e-10_dp,& 
                         war=r2/(u1+r3)
complex(dpc),parameter:: ca=(o2,o2),wa=(war,u0)
complex(dpc)          :: z0,za,zp,ww,w
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0; if(real (za)>u0)k=1; if(aimag(za)>u0)k=k+2
select case(k)
case(0)
case(1); za=-ci*za
case(2); za= ci*za
case(3); za=   -za
end select
z0=ca+za
za =z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative octant
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
else                    ! Expand about the vertex
   kn=aimag(z0)>real(z0)! negative sector of rotated z0 frame
   if(kn)z0 =-ci*z0     ! <-perform the rotation
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p045*w
   else;       w =p315*w
   endif
   w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z225*w
case(1); w =z315*w
case(2); w =z135*w
case(3); w =z045*w
end select
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine cuztoc
!=============================================================================
subroutine cuztoc_d(m,co1,co2,z,xc,xcd)!  [cuztoc]
!=============================================================================
! Conformal cubic mapping to the unit sphere.
! Map from the unit square map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,c0,c1,ci, &
                 p045=>z060,p315=>z300,  z045,z135,z225,z315
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,r3p=r3+u1,r2or3p=r2/r3p, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2),wa=(r2or3p,u0)
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0; if(real (za)>u0)k=1; if(aimag(za)>u0)k=k+2
select case(k)
case(0); zad= c1
case(1); zad=-ci; za=-ci*za
case(2); zad= ci; za= ci*za
case(3); zad=-c1; za=   -za
end select
z0=ca+za;        z0d=zad
zad=z135*r2*zad; za =z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative octant
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
else                    ! Expand about the vertex
   kn=aimag(z0)>real(z0)! negative sector of rotated z0 frame
   if(kn)then; z0d=-ci*z0d; z0 =-ci*z0      ! <-perform the rotation
   endif
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p045*w; wd=p045*wd
   else;       w =p315*w; wd=p315*wd
   endif
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2; w=(wa-w)/(wa*w+c1)
endif
select case(k)
case(0); w =z225*w; wd=z225*wd
case(1); w =z315*w; wd=z315*wd
case(2); w =z135*w; wd=z135*wd
case(3); w =z045*w; wd=z045*wd
end select
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine cuztoc_d
!=============================================================================
subroutine cuztoc0(m,co1,co2,z0,xc)!  [cuztoc0]
!=============================================================================
! Conformal cubic mapping to the unit sphere.
! Map from the cube's fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,c0,c1,ci, p045=>z060,z135
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,r3p=r3+u1,r2or3p=r2/r3p,&
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2),wa=(r2or3p,u0)
complex(dpc)          :: za,zp,ww,w
logical               :: ka,kn1
!=============================================================================
if(aimag(z0)>real(z0)+eps)stop 'y>x; Outside fundamental region'
if(real(z0)>o2+eps)       stop 'x>1/2; Outside fundamental region'
if(aimag(z0)<-eps)        stop 'y<0; Outside fundamental region'
za=z0-ca
za =z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w =ww**of
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w=c0
   endif
! Transform to vertex-centered frame from face-centered frame
   w=(wa-w)/(wa*w+c1)
   w =p045*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w =c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine cuztoc0
!=============================================================================
subroutine cuztoc0_d(m,co1,co2,z0,xc,xcd)!  [cuztoc0]
!=============================================================================
! Conformal cubic mapping to the unit sphere.
! Map from the cube's fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,r2,r3,c0,c1,ci, p045=>z060,z135
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=4,v=3,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,r3p=r3+u1,r2or3p=r2/r3p, &
                         eps=1.e-10_dp
complex(dpc),parameter:: ca=(o2,o2),wa=(r2or3p,u0)
complex(dpc)          :: za,zp,ww,zad,z0d,w,wd
logical               :: ka,kn1
!=============================================================================
if(aimag(z0)>real(z0)+eps)stop 'y>x; Outside fundamental region'
if(real(z0)>o2+eps)       stop 'x>1/2; Outside fundamental region'
if(aimag(z0)<-eps)        stop 'y<0; Outside fundamental region'
za=z0-ca;        zad=c1
z0d=zad
zad=z135*r2*zad; za =z135*r2*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w =ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
! Transform to vertex-centered frame from face-centered frame
   wd=-wd*(c1+wa**2)/(wa*w+c1)**2; w=(wa-w)/(wa*w+c1)
   wd=p045*wd; w =p045*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine cuztoc0_d

!=============================================================================
subroutine hhztoc(m,co1,co2,z,xc)!  [hhztoc]
!=============================================================================
! Conformal hexagonal dihedron mapping to the unit sphere.
! Map from the unit hexagonal map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or3,pi2,c0,c1,ci, &
                 p060=>z090,p300=>z270,      z060,z120,z240,z300
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=6,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3
complex(dpc),parameter:: ca=z060
complex(dpc)          :: z0,za,zp,ww,w
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); za=z120*za
case(1); za=z060*za
case(2)
case(3); za=z300*za
case(4); za=z240*za
case(5); za=-za
end select
z0=ca+z240*za
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w=c0
   endif
else                     ! Expand about the vertex
   kn=real(z0)<o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)z0 =z240*z0     ! <-perform the rotation
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)w =conjg(w)
   else; w =c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p060*w
   else;       w =p300*w
   endif
   w=(c1-w)/(w+c1)
endif
select case(k)
case(0); w =z240*w
case(1); w =z300*w
case(2)
case(3); w =z060*w
case(4); w =z120*w
case(5); w =    -w
end select
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine hhztoc
!=============================================================================
subroutine hhztoc_d(m,co1,co2,z,xc,xcd)!  [hhztoc]
!=============================================================================
! Conformal hexagonal dihedron mapping to the unit sphere.
! Map from the unit hexagonal map face to cartesians oriented
! with respect to the face center.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or3,pi2,c0,c1,ci, &
                 p060=>z090,p300=>z270,      z060,z120,z240,z300
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=6,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3
complex(dpc),parameter:: ca=z060
complex(dpc)          :: z0,za,zp,ww,zad,z0d,w,wd
integer(spi)          :: k
logical               :: ka,kn,kn1
!=============================================================================
za=z-ca
k=0;
if(abs(za)/=u0)k=modulo(floor(f*atan2(real(za),-aimag(za))/pi2)+1,f)
select case(k)
case(0); zad=z120; za=z120*za
case(1); zad=z060; za=z060*za
case(2); zad=c1  
case(3); zad=z300; za=z300*za
case(4); zad=z240; za=z240*za
case(5); zad=-c1 ; za=-za
end select
z0=ca+z240*za; z0d=z240*zad
ka=abs(za)<abs(z0)
if(ka)then         ! Expand about the face center:
   kn=aimag(za)<u0 ! negative sector
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
else                     ! Expand about the vertex
   kn=real(z0)<o2*abs(z0)! negative sector of rotated z0 frame
   if(kn)then; z0d=z240*z0d; z0 =z240*z0      ! <-perform the rotation
   endif
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov; wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn.neqv.kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
! Convert from vertex-centered frame to face-centered frame:
   if(kn)then; w =p060*w; wd=p060*wd
   else;       w =p300*w; wd=p300*wd
   endif
   wd=-wd*(c1*2)/(w+c1)**2; w=(c1-w)/(w+c1)
endif
select case(k)
case(0); w =z240*w; wd=z240*wd
case(1); w =z300*w; wd=z300*wd
case(2)
case(3); w =z060*w; wd=z060*wd
case(4); w =z120*w; wd=z120*wd
case(5); w =    -w; wd=    -wd
end select
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine hhztoc_d
!=============================================================================
subroutine hhztoc0(m,co1,co2,z0,xc)!  [hhztoc0]
!=============================================================================
! Conformal hexagonal dihedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or3,c0,c1,ci, &
                 p060=>z090,p300=>z270,     z060,z120,z240
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=6,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,eps=1.e-10_dp
complex(dpc),parameter:: ca=z060
complex(dpc)          :: za,zp,ww,w
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(real(z0)<o2*abs(z0)-eps) stop 'x<|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca
za=z120*za
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww)
   if(abs(ww)>u0)then
      w=ww**of
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w=c0
   endif
   w=(c1-w)/(w+c1)
   w =p060*w
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww)
   if(abs(ww)>u0)then
      w=ww**ov
      kn1=aimag(w)<u0
      if(kn1)w =conjg(w)
   else; w =c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc)
end subroutine hhztoc0
!=============================================================================
subroutine hhztoc0_d(m,co1,co2,z0,xc,xcd)!  [hhztoc0]
!=============================================================================
! Conformal hexagonal dihedron mapping to the unit sphere.
! Map from the fundamental region to cartesians oriented
! with respect to the vertex.
!=============================================================================
use pietc, only: FF=>F,u0,o2,u1,or3,c0,c1,ci, &
                 p060=>z090,p300=>z270,     z060,z120,z240
use pfft1, only: tay
use pmat4, only: ztoc
implicit none
integer(spi),             intent(in ):: m
real(dp),dimension(m),    intent(in ):: co1,co2
complex(dpc),             intent(in ):: z0
real(dp),dimension(3),    intent(out):: xc
complex(dpc),dimension(3),intent(out):: xcd
!-----------------------------------------------------------------------------
integer(spi),parameter:: f=6,v=2,v0=(2*f)/(f-2)! <- Face and vertex Symmetries
real(dp),parameter    :: of=u1/f,ov=u1/v,ov0=u1/v0,re=or3,eps=1.e-10_dp
complex(dpc),parameter:: ca=z060
complex(dpc)          :: za,zp,ww,zad,z0d,w,wd
logical               :: ka,kn1
!=============================================================================
if(real(z0)>o2+eps)         stop 'x>1/2; Outside fundamental region'
if(real(z0)<o2*abs(z0)-eps) stop 'x<|z|/2; Outside fundamental region'
if(aimag(z0)<-eps)          stop 'y<0;  Outside fundamental region'
za=z0-ca
zad=z120; za=z120*za
z0d=z240*zad
ka=abs(za)<abs(z0)
if(ka)then           ! Expand about the face center:
   zp=za**f
   call tay(m,co2,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**of; wd=of*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*f*zp*zad/za
   else; w=c0; wd=cmplx(co2(1)**of,u0,dpc)*zad
   endif
   wd=-wd*(c1*2)/(w+c1)**2; w=(c1-w)/(w+c1)
   w =p060*w; wd=p060*wd
else       ! Expand about the vertex
   zp=z0**v0
   call tay(m,co1,zp,ww,wd)
   if(abs(ww)>u0)then
      w=ww**ov
      wd=ov*w*wd/ww
      kn1=aimag(w)<u0
      if(kn1)then; w =conjg(w); wd=conjg(wd)
      endif
      wd=wd*v0*zp*z0d/z0
   else; w =c0; wd=c0
   endif
endif
! Then convert to cartesians
call ztoc(w,FF,xc,xcd)
xcd=xcd*wd
end subroutine hhztoc0_d


end module pcoco

