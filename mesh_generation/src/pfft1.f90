!
!
!                                                         ******************
!                                                         *   pfft1.f90
!                                                         *
!                                                         *  R. J.
!                                                         Purser  *
!                                                         *
!                                                         NOAA/NCEP/EMC
!                                                         *
!                                                         *  22st Jun
!                                                         2012 *
!                                                         ******************
! MODULE PFFT1: Evaluate Taylor and Laurent series and their
! derivatives,
!             and recover the coefficients of such series from complex
!             Fourier expansions of the functions around circles that
!             lie
!             within the central expansion's disc of convergence.
!
! MODULE:  pfft1.mod                     LIBRARY: libpfft.a
! DEPENDENCIES:
! Libraries:           
! Modules:   pfft, pietc, pkind
! 
! For Taylor series, include the option of evaluating terms in a
! periodic
! pattern (skipping intermedaite evaluations)
!=============================================================================
module pfft1
!=============================================================================
use pkind, only: spi,dp,dpc
implicit none
private
public:: tay,tay0,lau,tay_reset,tay0_reset,lau_reset,lauarr,lausuf, &
         cinvrt,conv1

interface tay
   module procedure tay, tayd, taydd, tays, taysd, taysdd,&
                   ctay,ctayd,ctaydd,ctays,ctaysd,ctaysdd;       end interface

interface tay0
   module procedure tay0, tay0d, tay0dd, tay0s, tay0sd, tay0sdd, &
                   ctay0,ctay0d,ctay0dd,ctay0s,ctay0sd,ctay0sdd; end interface

interface lau
   module procedure lau,laud,laudd, clau,lau2,claud,lau2d,claudd,lau2dd
                                                                 end interface

interface tay_reset
   module procedure tay_reset,tay_resetx,ctay_reset,ctay_resetx; end interface

interface tay0_reset
   module procedure tay0_reset,tay0_resetx, ctay0_reset,ctay0_resetx   
                                                                 end interface

interface lau_reset
   module procedure lau_reset,lau_resetx, &
        clau_reset,lau2_reset,clau_resetx,lau2_resetx;           end interface

interface lauarr
   module procedure lauarr,lau2arr,clauarr;                      end interface

interface lausuf
   module procedure lausuf;                                      end interface

interface cinvrt; module procedure dcinvrt;                      end interface
interface conv1;  module procedure dconv1;                       end interface


contains

!=============================================================================
subroutine tay(n,aco,z,w)!  [tay]
!=============================================================================
integer(spi),         intent(IN ):: n
real(dp),dimension(n),intent(IN ):: aco
complex(dpc),         intent(IN ):: z
complex(dpc),         intent(OUT):: w
!=============================================================================
call lau(1,n,aco,z,w)
end subroutine tay

!=============================================================================
subroutine tayd(n,aco,z,w,wd)!  [tay]
!=============================================================================
integer(spi),         intent(IN ):: n
real(dp),dimension(n),intent(IN ):: aco
complex(dpc),         intent(IN ):: z
complex(dpc),         intent(OUT):: w,wd
!=============================================================================
call lau(1,n,aco,z,w,wd)
end subroutine tayd

!=============================================================================
subroutine taydd(n,aco,z,w,wd,wdd)!  [tay]
!=============================================================================
integer(spi),         intent(IN ):: n
real(dp),dimension(n),intent(IN ):: aco
complex(dpc),         intent(IN ):: z
complex(dpc),         intent(OUT):: w,wd,wdd
!=============================================================================
call lau(1,n,aco,z,w,wd,wdd)
end subroutine taydd

!=============================================================================
subroutine tays(n,s,aco,z,w)!  [tay]
!=============================================================================
! "skipped" version, evaluating every s terms, 1, s+1, etc.
integer(spi),         intent(IN ):: n,s
real(dp),dimension(n),intent(IN ):: aco
complex(dpc),         intent(IN ):: z
complex(dpc),         intent(OUT):: w
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs
!=============================================================================
zs=z**s
w=0
m=1+s*((n-1)/s)
do i=m,1,-s
   w=aco(i)+w*zs
enddo
w=w*z
end subroutine tays

!=============================================================================
subroutine taysd(n,s,aco,z,w,wd)!  [tay]
!=============================================================================
integer(spi),         intent(IN ):: n,s
real(dp),dimension(n),intent(IN ):: aco
complex(dpc),         intent(IN ):: z
complex(dpc),         intent(OUT):: w,wd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
else
   zsd=z**(s-1)
   zs=z*zsd
   zsd=s*zsd
endif
w=0
wd=0
m=1+s*((n-1)/s)
do i=m,1,-s
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
wd=w+wd*z
w=w*z
end subroutine taysd

!=============================================================================
subroutine taysdd(n,s,aco,z,w,wd,wdd)!  [tay]
!=============================================================================
integer(spi),         intent(IN ):: n,s
real(dp),dimension(n),intent(IN ):: aco
complex(dpc),         intent(IN ):: z
complex(dpc),         intent(OUT):: w,wd,wdd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd,zsdd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
   zsdd=0
elseif(s==2)then
   zs=z*z
   zsd=2*z
   zsdd=2
else
   zsdd=z**(s-2)
   zsd=z*zsdd
   zs=z*zsd
   zsd=s*zsd
   zsdd=s*(s-1)*zsdd
endif
w=0
wd=0
wdd=0
m=1+s*((n-1)/s)
do i=m,1,-s
   wdd=w*zsdd+2*wd*zsd+wdd*zs
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
wdd=2*wd+wdd*z
wd=w+wd*z
w=w*z
end subroutine taysdd

!=============================================================================
subroutine ctay(n,aco,z,w)!  [tay]
!=============================================================================
integer(spi),             intent(IN ):: n
complex(dpc),dimension(n),intent(IN ):: aco
complex(dpc),             intent(IN ):: z
complex(dpc),             intent(OUT):: w
!=============================================================================
call lau(1,n,aco,z,w)
end subroutine ctay

!=============================================================================
subroutine ctayd(n,aco,z,w,wd)!  [tay]
!=============================================================================
integer(spi),                  intent(IN ):: n
complex(dpc),dimension(n),intent(IN ):: aco
complex(dpc),             intent(IN ):: z
complex(dpc),             intent(OUT):: w,wd
!=============================================================================
call lau(1,n,aco,z,w,wd)
end subroutine ctayd

!=============================================================================
subroutine ctaydd(n,aco,z,w,wd,wdd)!  [tay]
!=============================================================================
integer(spi),             intent(IN ):: n
complex(dpc),dimension(n),intent(IN ):: aco
complex(dpc),             intent(IN ):: z
complex(dpc),             intent(OUT):: w,wd,wdd
!=============================================================================
call lau(1,n,aco,z,w,wd,wdd)
end subroutine ctaydd

!=============================================================================
subroutine ctays(n,s,aco,z,w)!  [tay]
!=============================================================================
integer(spi),             intent(IN ):: n,s
complex(dpc),dimension(n),intent(IN ):: aco
complex(dpc),             intent(IN ):: z
complex(dpc),             intent(OUT):: w
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs
!=============================================================================
zs=z**s
w=0
m=1+s*((n-1)/s)
do i=m,1,-s
   w=aco(i)+w*zs
enddo
w=w*z
end subroutine ctays

!=============================================================================
subroutine ctaysd(n,s,aco,z,w,wd)!  [tay]
!=============================================================================
integer(spi),             intent(IN ):: n,s
complex(dpc),dimension(n),intent(IN ):: aco
complex(dpc),             intent(IN ):: z
complex(dpc),             intent(OUT):: w,wd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
else
   zsd=z**(s-1)
   zs=z*zsd
   zsd=s*zsd
endif
w=0
wd=0
m=1+s*((n-1)/s)
do i=m,1,-s
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
wd=w+wd*z
w=w*z
end subroutine ctaysd

!=============================================================================
subroutine ctaysdd(n,s,aco,z,w,wd,wdd)!  [tay]
!=============================================================================
integer(spi),             intent(IN ):: n,s
complex(dpc),dimension(n),intent(IN ):: aco
complex(dpc),             intent(IN ):: z
complex(dpc),             intent(OUT):: w,wd,wdd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd,zsdd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
   zsdd=0
elseif(s==2)then
   zs=z*z
   zsd=2*z
   zsdd=2
else
   zsdd=z**(s-2)
   zsd=z*zsdd
   zs=z*zsd
   zsd=s*zsd
   zsdd=s*(s-1)*zsdd
endif
w=0
wd=0
wdd=0
m=1+s*((n-1)/s)
do i=m,1,-s
   wdd=w*zsdd+2*wd*zsd+wdd*zs
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
wdd=2*wd+wdd*z
wd=w+wd*z
w=w*z
end subroutine ctaysdd

!=============================================================================
subroutine tay0(n,aco,z,w)!  [tay0]
!=============================================================================
integer(spi),           intent(IN ):: n
real(dp),dimension(0:n),intent(IN ):: aco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w
!=============================================================================
print'('' In tay0; deprecated routine: change name to lau and add 0 arg'')'
call lau(0,n,aco,z,w)
end subroutine tay0
!=============================================================================
subroutine tay0d(n,aco,z,w,wd)!  [tay0]
!=============================================================================
integer(spi),           intent(IN ):: n
real(dp),dimension(0:n),intent(IN ):: aco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w,wd
!=============================================================================
print'('' In tay0; deprecated routine: change name to lau and add 0 arg'')'
call lau(0,n,aco,z,w,wd)
end subroutine tay0d

!=============================================================================
subroutine tay0dd(n,aco,z,w,wd,wdd)!  [tay0]
!=============================================================================
integer(spi),           intent(IN ):: n
real(dp),dimension(0:n),intent(IN ):: aco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w,wd,wdd
!=============================================================================
print'('' In tay0; deprecated routine: change name to lau and add 0 arg'')'
call lau(0,n,aco,z,w,wd,wdd)
end subroutine tay0dd

!=============================================================================
subroutine tay0s(n,s,aco,z,w)!  [tay0]
!=============================================================================
integer(spi),           intent(IN ):: n,s
real(dp),dimension(0:n),intent(IN ):: aco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs
!=============================================================================
zs=z**s
w=0
m=s*(n/s)
do i=m,0,-s
   w=aco(i)+w*zs
enddo
end subroutine tay0s

!=============================================================================
subroutine tay0sd(n,s,aco,z,w,wd)!  [tay0]
!=============================================================================
integer(spi),           intent(IN ):: n,s
real(dp),dimension(0:n),intent(IN ):: aco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w,wd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
else
   zsd=z**(s-1)
   zs=z*zsd
   zsd=s*zsd
endif
w=0
wd=0
m=s*(n/s)
do i=m,0,-s
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
end subroutine tay0sd

!=============================================================================
subroutine tay0sdd(n,s,aco,z,w,wd,wdd)!  [tay0]
!=============================================================================
integer(spi),           intent(IN ):: n,s
real(dp),dimension(0:n),intent(IN ):: aco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w,wd,wdd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd,zsdd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
   zsdd=0
elseif(s==2)then
   zs=z*z
   zsd=2*z
   zsdd=2
else
   zsdd=z**(s-2)
   zsd=z*zsdd
   zs=z*zsd
   zsd=s*zsd
   zsdd=s*(s-1)*zsdd
endif
w=0
wd=0
wdd=0
m=s*(n/s)
do i=m,0,-s
   wdd=w*zsdd+2*wd*zsd+wdd*zs
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
end subroutine tay0sdd

!=============================================================================
subroutine ctay0(n,aco,z,w)!  [tay0]
!=============================================================================
integer(spi),               intent(IN ):: n
complex(dpc),dimension(0:n),intent(IN ):: aco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w
!=============================================================================
print '(" In tay0; deprecated routine: change name to lau and add 0 arg")'
call lau(0,n,aco,z,w)
end subroutine ctay0

!=============================================================================
subroutine ctay0d(n,aco,z,w,wd)!  [tay0]
!=============================================================================
integer(spi),               intent(IN ):: n
complex(dpc),dimension(0:n),intent(IN ):: aco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w,wd
!=============================================================================
print '(" In tay0; deprecated routine: change name to lau and add 0 arg")'
call lau(0,n,aco,z,w,wd)
end subroutine ctay0d

!=============================================================================
subroutine ctay0dd(n,aco,z,w,wd,wdd)!  [tay0]
!=============================================================================
integer(spi),               intent(IN ):: n
complex(dpc),dimension(0:n),intent(IN ):: aco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w,wd,wdd
!=============================================================================
print '(" In tay0; deprecated routine: change name to lau and add 0 arg")'
call lau(0,n,aco,z,w,wd,wdd)
end subroutine ctay0dd

!=============================================================================
subroutine ctay0s(n,s,aco,z,w)!  [tay0]
!=============================================================================
integer(spi),               intent(IN ):: n,s
complex(dpc),dimension(0:n),intent(IN ):: aco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs
!=============================================================================
zs=z**s
w=0
m=s*(n/s)
do i=m,0,-s
   w=aco(i)+w*zs
enddo
end subroutine ctay0s

!=============================================================================
subroutine ctay0sd(n,s,aco,z,w,wd)!  [tay0]
!=============================================================================
integer(spi),               intent(IN ):: n,s
complex(dpc),dimension(0:n),intent(IN ):: aco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w,wd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
else
   zsd=z**(s-1)
   zs=z*zsd
   zsd=s*zsd
endif
w=0
wd=0
m=s*(n/s)
do i=m,0,-s
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
end subroutine ctay0sd

!=============================================================================
subroutine ctay0sdd(n,s,aco,z,w,wd,wdd)!  [tay0]
!=============================================================================
integer(spi),               intent(IN ):: n,s
complex(dpc),dimension(0:n),intent(IN ):: aco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w,wd,wdd
!-----------------------------------------------------------------------------
integer(spi):: i,m
complex(dpc):: zs,zsd,zsdd
!=============================================================================
if(s==1)then
   zs=z
   zsd=1
   zsdd=0
elseif(s==2)then
   zs=z*z
   zsd=2*z
   zsdd=2
else
   zsdd=z**(s-2)
   zsd=z*zsdd
   zs=z*zsd
   zsd=s*zsd
   zsdd=s*(s-1)*zsdd
endif
w=0
wd=0
wdd=0
m=s*(n/s)
do i=m,0,-s
   wdd=w*zsdd+2*wd*zsd+wdd*zs
   wd=w*zsd+wd*zs
   w=aco(i)+w*zs
enddo
end subroutine ctay0sdd

! Routines for evaluating Laurent (double-sided) power series
!=============================================================================
subroutine lau(m,n,fco,z,w)!  [lau]
!=============================================================================
! Evaluate the Laurent series with powers from -m to n, at complex z.
!=============================================================================
use pietc, only: u1
integer(spi),           intent(IN ):: m,n
real(dp),dimension(m:n),intent(IN ):: fco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zi
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=fco(0)
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      w=w+fco(-k)*zp
   enddo
endif
zp=u1
do k=1,n
   zp=zp*z
   w=w+fco(k)*zp
enddo
end subroutine lau
!=============================================================================
subroutine laud(m,n,fco,z,w,wd)!  [lau]
!=============================================================================
! Like lau, but evaluate also the derivative
!=============================================================================
use pietc, only: u1
integer(spi),           intent(IN ):: m,n
real(dp),dimension(m:n),intent(IN ):: fco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w,wd
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zp1,zi,tk,tk1
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=fco(0)
wd=0
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      tk=fco(-k)*zp
      w=w+tk
      wd=wd-k*tk
   enddo
   wd=wd*zi
endif
zp=u1
do k=1,n
   zp1=zp
   zp=zp*z
   tk=fco(k)*zp
   tk1=k*fco(k)*zp1
   w=w+tk
   wd=wd+tk1
enddo
end subroutine laud
!=============================================================================
subroutine laudd(m,n,fco,z,w,wd,wdd)!  [lau]
!=============================================================================
! Like lau, but also evaluate two derivatives
!=============================================================================
use pietc, only: u1
integer(spi),           intent(IN ):: m,n
real(dp),dimension(m:n),intent(IN ):: fco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w,wd,wdd
!-----------------------------------------------------------------------------
integer(spi)       :: k
complex(dpc)      :: zp,zp1,zp2,zi,tk,tk1,tk2
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=fco(0)
wd=0
wdd=0
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      tk=fco(-k)*zp
      w=w+tk
      tk=-k*tk
      wd=wd+tk
      wdd=wdd-(k+1)*tk
   enddo
   wd=wd*zi
   wdd=wdd*zi*zi
endif
zp=u1
zp1=0
do k=1,n
   zp2=zp1
   zp1=zp
   zp=zp*z
   tk=fco(k)*zp
   tk1=k*fco(k)*zp1
   tk2=(k-1)*k*fco(k)*zp2
   w=w+tk
   wd=wd+tk1
   wdd=wdd+tk2
enddo
end subroutine laudd

!=============================================================================
subroutine clau(m,n,fco,z,w)!  [lau]
!=============================================================================
! Evaluate the Laurent series with powers from -m to n, at complex z.
!=============================================================================
use pietc, only: u1=>c1
integer(spi),               intent(IN ):: m,n
complex(dpc),dimension(m:n),intent(IN ):: fco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zi
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m>=0)w=fco(0)
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      w=w+fco(-k)*zp
   enddo
endif
zp=u1
do k=1,n
   zp=zp*z
   w=w+fco(k)*zp
enddo
end subroutine clau
!=============================================================================
subroutine lau2(m,n,rco,qco,z,w)!  [lau]
!=============================================================================
! Like clau, except coefficients stored as pair of real arrays
!=============================================================================
use pietc, only: u1=>c1
integer(spi),           intent(IN ):: m,n
real(dp),dimension(m:n),intent(IN ):: rco,qco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zi
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=cmplx(rco(0),qco(0),dpc)
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      w=w+cmplx(rco(-k),qco(-k),dpc)*zp
   enddo
endif
zp=u1
do k=1,n
   zp=zp*z
   w=w+cmplx(rco(k),qco(k),dpc)*zp
enddo
end subroutine lau2 
!=============================================================================
subroutine claud(m,n,fco,z,w,wd)!  [lau]
!=============================================================================
! Like lau, but evaluate also the derivative
!=============================================================================
use pietc, only: u1=>c1
integer(spi),               intent(IN ):: m,n
complex(dpc),dimension(m:n),intent(IN ):: fco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w,wd
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zp1,zi,tk,tk1
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=fco(0)
wd=0
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      tk=fco(-k)*zp
      w=w+tk
      wd=wd-k*tk
   enddo
   wd=wd*zi
endif
zp=u1
do k=1,n
   zp1=zp
   zp=zp*z
   tk=fco(k)*zp
   tk1=k*fco(k)*zp1
   w=w+tk
   wd=wd+tk1
enddo
end subroutine claud
!=============================================================================
subroutine lau2d(m,n,rco,qco,z,w,wd)!  [lau]
!=============================================================================
! Like lau, but evaluate also the derivative
!=============================================================================
use pietc, only: u1=>c1
integer(spi),           intent(IN ):: m,n
real(dp),dimension(m:n),intent(IN ):: rco,qco
complex(dpc),           intent(IN ):: z
complex(dpc),           intent(OUT):: w,wd
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zp1,zi,tk,tk1,fco
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=cmplx(rco(0),qco(0),dpc)
wd=0
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      fco=cmplx(rco(-k),qco(-k),dpc)
      tk=fco*zp
      w=w+tk
      wd=wd-k*tk
   enddo
   wd=wd*zi
endif
zp=u1
do k=1,n
   zp1=zp
   zp=zp*z
   fco=cmplx(rco(k),qco(k),dpc)
   tk=fco*zp
   tk1=k*fco*zp1
   w=w+tk
   wd=wd+tk1
enddo
end subroutine lau2d
!=============================================================================
subroutine claudd(m,n,fco,z,w,wd,wdd)!  [lau]
!=============================================================================
! Like lau, but also evaluate two derivatives
!=============================================================================
use pietc, only: u1=>c1
integer(spi),               intent(IN ):: m,n
complex(dpc),dimension(m:n),intent(IN ):: fco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w,wd,wdd
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zp1,zp2,zi,tk,tk1,tk2
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=fco(0)
wd=0
wdd=0
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      tk=fco(-k)*zp
      w=w+tk
      tk=-k*tk
      wd=wd+tk
      wdd=wdd-(k+1)*tk
   enddo
   wd=wd*zi
   wdd=wdd*zi*zi
endif
zp=u1
zp1=0
do k=1,n
   zp2=zp1
   zp1=zp
   zp=zp*z
   tk=fco(k)*zp
   tk1=k*fco(k)*zp1
   tk2=(k-1)*k*fco(k)*zp2
   w=w+tk
   wd=wd+tk1
   wdd=wdd+tk2
enddo
end subroutine claudd
!=============================================================================
subroutine lau2dd(m,n,rco,qco,z,w,wd,wdd)!  [lau]
!=============================================================================
! Like lau, but also evaluate two derivatives
!=============================================================================
use pietc, only: u1=>c1
integer(spi),               intent(IN ):: m,n
real(dp),dimension(m:n),    intent(IN ):: rco,qco
complex(dpc),               intent(IN ):: z
complex(dpc),               intent(OUT):: w,wd,wdd
!-----------------------------------------------------------------------------
integer(spi)      :: k
complex(dpc)      :: zp,zp1,zp2,zi,tk,tk1,tk2,fco
!=============================================================================
if(m>1)stop 'In lau; lower index bound of coefficient array must be <= 1'
w=0; if(m<=0)w=cmplx(rco(0),qco(0),dpc)
wd=0
wdd=0
if(m<0)then
   zp=u1
   zi=u1/z
   do k=1,-m
      zp=zp*zi
      fco=cmplx(rco(-k),qco(-k),dpc)
      tk=fco*zp
      w=w+tk
      tk=-k*tk
      wd=wd+tk
      wdd=wdd-(k+1)*tk
   enddo
   wd=wd*zi
   wdd=wdd*zi*zi
endif
zp=u1
zp1=0
do k=1,n
   zp2=zp1
   zp1=zp
   zp=zp*z
   fco=cmplx(rco(k),qco(k),dpc)
   tk=fco*zp
   tk1=k*fco*zp1
   tk2=(k-1)*k*fco*zp2
   w=w+tk
   wd=wd+tk1
   wdd=wdd+tk2
enddo
end subroutine lau2dd

!============================================================================
subroutine tay_reset(nf,nco,zzf,r,co)!  [tay_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed real, from the results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r.
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
real(dp),dimension(1:nco),     intent(  OUT):: co
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp
integer(spi)      :: i
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=1,nco
   rp=rp*ri
   co(i)=real(zzf(i))*rp
enddo
end subroutine tay_reset

!============================================================================
subroutine tay_resetx(nf,nco,zzf,r,co,norm)!  [tay_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed real, from the results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
real(dp),dimension(1:nco),     intent(  OUT):: co
real(dp),                      intent(  OUT):: norm
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp,s,eqpos,ezer,eqneg,erneg
integer(spi)      :: i
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=1,nco
   rp=rp*ri
   co(i)=real(zzf(i))*rp
enddo
s=maxval(abs(real(zzf(1:nco))))
eqpos=maxval(abs(aimag(zzf(0:nco))))
ezer=        abs( real(zzf(0)))
erneg=maxval(abs( real(zzf(nf-nco:nf-1))))
eqneg=maxval(abs(aimag(zzf(nf-nco:nf-1))))
norm=max(eqpos,ezer,erneg,eqneg)/s
end subroutine tay_resetx

!============================================================================
subroutine tay0_reset(nf,nco,zzf,r,co)!  [tay0_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed real, from the results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r 
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
real(dp),dimension(0:nco),     intent(  OUT):: co
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp
integer(spi)      :: i
!============================================================================
print'(" In tay0_reset; deprecated routine; use lau_reset and add 0 arg")'
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=0,nco
   co(i)=real(zzf(i))*rp
   rp=rp*ri
enddo
end subroutine tay0_reset

!============================================================================
subroutine tay0_resetx(nf,nco,zzf,r,co,norm)!  [tay0_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed real, from the results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
real(dp),dimension(0:nco),     intent(  OUT):: co
real(dp),                      intent(  OUT):: norm
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp,s,eqpos,erneg,eqneg
integer(spi)      :: i
!============================================================================
print'(" In tay0_reset; deprecated routine; use lau_reset and add 0 arg")'
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=0,nco
   co(i)=real(zzf(i))*rp
   rp=rp*ri
enddo
s=maxval(abs(real(zzf(0:nco))))
eqpos=maxval(abs(aimag(zzf(0:nco))))
erneg=maxval(abs( real(zzf(nf-nco:nf-1))))
eqneg=maxval(abs(aimag(zzf(nf-nco:nf-1))))
norm=max(eqpos,erneg,eqneg)/s
end subroutine tay0_resetx

!============================================================================
subroutine lau_reset(nf,lco,nco,zzf,r,co)!  [lau_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed real, from the results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,lco,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
real(dp),dimension(lco:nco),   intent(  OUT):: co
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp
integer(spi)      :: i,j
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=r**(-lco)
do i=lco,nco
   j=mod(i+nf,nf)
   co(i)=real(zzf(j))*rp
   rp=rp*ri
enddo
end subroutine lau_reset

!============================================================================
subroutine lau_resetx(nf,lco,nco,zzf,r,co,norm)!  [lau_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed real, from the results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,lco,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
real(dp),dimension(lco:nco),   intent(  OUT):: co
real(dp),                      intent(  OUT):: norm
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp,s,eqpos,erneg,eqneg
integer(spi)      :: i,j,plco
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=r**(-lco)
do i=lco,nco
   j=mod(i+nf,nf)
   co(i)=real(zzf(j))*rp
   rp=rp*ri
enddo
plco=min(lco,0)
s=max(maxval(abs( real(zzf(0:nco)))),maxval(abs(real(zzf(nf+lco:nf-1)))))
eqpos=maxval(abs(aimag(zzf(0:nco))))
erneg=maxval(abs( real(zzf(nf-nco:nf+plco-1))))
eqneg=maxval(abs(aimag(zzf(nf-nco:nf-1))))
norm=max(eqpos,eqneg,erneg)/s
end subroutine lau_resetx

!============================================================================
subroutine ctay_reset(nf,nco,zzf,r,co)!  [tay_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed complex, from the
! results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
complex(dpc),dimension(1:nco), intent(  OUT):: co
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp
integer(spi)      :: i
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=1,nco
   rp=rp*ri
   co(i)=zzf(i)*rp
enddo
end subroutine ctay_reset

!============================================================================
subroutine ctay_resetx(nf,nco,zzf,r,co,norm)!  [tay_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed complex, from the
! results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
complex(dpc),dimension(1:nco), intent(  OUT):: co
real(dp),                      intent(  OUT):: norm
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp,s,ezer,erneg,eqneg
integer(spi)      :: i
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=1,nco
   rp=rp*ri
   co(i)=zzf(i)*rp
enddo
s=max(maxval(abs(real(zzf(1:nco)))),maxval(abs(aimag(zzf(1:nco)))))
ezer=max(abs(real(zzf(0))),abs(aimag(zzf(0))))
erneg=maxval(abs( real(zzf(nf-nco:nf-1))))
eqneg=maxval(abs(aimag(zzf(nf-nco:nf-1))))
norm=max(ezer,erneg,eqneg)/s
end subroutine ctay_resetx

!============================================================================
subroutine ctay0_reset(nf,nco,zzf,r,co)!  [tay0_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed complex, from the
! results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
complex(dpc),dimension(0:nco), intent(  OUT):: co
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp
integer(spi)      :: i
!============================================================================
print'(" In tay0_reset; deprecated routine; use lau_reset and add 0 arg")'
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=0,nco
   co(i)=zzf(i)*rp
   rp=rp*ri
enddo
end subroutine ctay0_reset

!============================================================================
subroutine ctay0_resetx(nf,nco,zzf,r,co,norm)!  [tay0_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed complex, from the
! results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
complex(dpc),dimension(0:nco), intent(  OUT):: co
real(dp),                      intent(  OUT):: norm
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp,s,erneg,eqneg
integer(spi)      :: i
!============================================================================
print'(" In tay0_reset; deprecated routine; use lau_reset and add 0 arg")'
call dfft(nf,zzf)
ri=u1/r
rp=u1
do i=0,nco
   co(i)=zzf(i)*rp
   rp=rp*ri
enddo
s=max(maxval(abs(real(zzf(0:nco)))),maxval(abs(aimag(zzf(0:nco)))))
erneg=maxval(abs( real(zzf(nf-nco:nf-1))))
eqneg=maxval(abs(aimag(zzf(nf-nco:nf-1))))
norm=max(erneg,eqneg)/s
end subroutine ctay0_resetx

!============================================================================
subroutine clau_reset(nf,lco,nco,zzf,r,co)!  [lau_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed complex, from the
! results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                   intent(IN   ):: nf,lco,nco
complex(dpc),dimension(0:nf-1), intent(INOUT):: zzf
real(dp),                       intent(IN   ):: r
complex(dpc),dimension(lco:nco),intent(  OUT):: co
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp
integer(spi)      :: i,j
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=r**(-lco)
do i=lco,nco
   j=mod(i+nf,nf)
   co(i)=zzf(j)*rp
   rp=rp*ri
enddo
end subroutine clau_reset
!============================================================================
subroutine lau2_reset(nf,lco,nco,zzf,r,rco,qco)!  [lau_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed real, from the results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                   intent(IN   ):: nf,lco,nco
complex(dpc),dimension(0:nf-1), intent(INOUT):: zzf
real(dp),                       intent(IN   ):: r
real(dp),dimension(lco:nco),    intent(  OUT):: rco,qco
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp
complex(dpc)      :: co
integer(spi)      :: i,j
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=r**(-lco)
do i=lco,nco
   j=mod(i+nf,nf)
   co=zzf(j)*rp
   rco(i)=real(co)
   qco(i)=aimag(co)
   rp=rp*ri
enddo
end subroutine lau2_reset
!============================================================================
subroutine clau_resetx(nf,lco,nco,zzf,r,co,norm)!  [lau_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed complex, from the
! results
! of a complex fourier transformation of cycle nf. The nominal radius of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                   intent(IN   ):: nf,lco,nco
complex(dpc),dimension(0:nf-1), intent(INOUT):: zzf
real(dp),                       intent(IN   ):: r
complex(dpc),dimension(lco:nco),intent(  OUT):: co
real(dp),                       intent(  OUT):: norm
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp,s,erneg,eqneg
integer(spi)      :: i,j,plco
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=r**(-lco)
do i=lco,nco
   j=mod(i+nf,nf)
   co(i)=zzf(j)*rp
   rp=rp*ri
enddo
plco=min(lco,0)
s=max(maxval(abs( real(zzf(0:nco)))),maxval(abs( real(zzf(nf+lco:nf-1)))),&
      maxval(abs(aimag(zzf(0:nco)))),maxval(abs(aimag(zzf(nf+lco:nf-1)))))
erneg=maxval(abs( real(zzf(nf-nco:nf+plco-1))))
eqneg=maxval(abs(aimag(zzf(nf-nco:nf+plco-1))))
norm=max(erneg,eqneg)/s
end subroutine clau_resetx 
!============================================================================
subroutine lau2_resetx(nf,lco,nco,zzf,r,rco,qco,norm)!  [lau_reset]
!============================================================================
! Reset the Taylor series coefficients, assumed complex, from the
! results
! of a complex fourier transformation of cycle nf, but store real and
! imag
! components separately in real arrays, rco and qco. The nominal radius
! of
! the circuit is r
!============================================================================
use pietc, only: u1
use pfft,  only: dfft
integer(spi),                  intent(IN   ):: nf,lco,nco
complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf
real(dp),                      intent(IN   ):: r
real(dp),dimension(lco:nco),   intent(  OUT):: rco,qco
real(dp),                      intent(  OUT):: norm
!-----------------------------------------------------------------------------
real(dp)          :: ri,rp,s,erneg,eqneg
complex(dpc)      :: co
integer(spi)      :: i,j,plco
!============================================================================
call dfft(nf,zzf)
ri=u1/r
rp=r**(-lco)
do i=lco,nco
   j=mod(i+nf,nf)
   co=zzf(j)*rp
   rco(i)=real(co)
   qco(i)=aimag(co)
   rp=rp*ri
enddo
plco=min(lco,0)
s=max(maxval(abs( real(zzf(0:nco)))),maxval(abs( real(zzf(nf+lco:nf-1)))),&
      maxval(abs(aimag(zzf(0:nco)))),maxval(abs(aimag(zzf(nf+lco:nf-1)))))
erneg=maxval(abs( real(zzf(nf-nco:nf+plco-1))))
eqneg=maxval(abs(aimag(zzf(nf-nco:nf+plco-1))))
norm=max(erneg,eqneg)/s
end subroutine lau2_resetx

!============================================================================
subroutine lauarr(lco,nco,eps,aco,arr)!  [lauarr]
!============================================================================
! From a given array ACO of Laurent expansion coefficients, create
! another
! array, ARR, of squared-magnitudes rr such that, at each index i, when
! the
! squared-abs-argument is no greater than arr(i), the terms
! contributed to the sum of the expansion terms at index j>i is 
! necessarily smaller than EPS-times
! the largest individual term formed taken only from ACO(lco:i).
! In other words, if eps is suitably small (such as 1.e-17) and the
! magnitude 
! of the argument is known and its square absolute value is smaller than 
! ARR(i) for some i, then
! we need only use aco(lco:i), not necessarily the whole of aco, in
! order to 
! get a numerical approximation to the function expanded with an error
! compatible with (or smaller than) the round-off error we could expect
! to
! have obtained by using the series in its entirety.
!============================================================================
integer(spi),               intent(IN ):: lco,nco
real(dp),                   intent(IN ):: eps
real(dp),dimension(lco:nco),intent(IN ):: aco
real(dp),dimension(lco:nco),intent(OUT):: arr
!-----------------------------------------------------------------------------
integer(spi),parameter     :: nit=8
real(dp),parameter         :: ten=10
real(dp),dimension(lco:nco):: absaco,laco,terms,lterms
real(dp)                   :: leps,leps4,nudge,lratio, &
                              lsumterms,maxlterms, &
                              la,lb,lc
integer(spi)                  :: i,j,it
!============================================================================
leps=log10(eps); leps4=leps*4; nudge=ten**leps4
absaco=abs(aco)
laco=leps4
do i=lco,nco; if(absaco(i)>nudge)laco(i)=log10(absaco(i));enddo
lc=leps4
do i=lco,nco-1
   la=lc
   lb=1
   do it=1,nit
      lc=(la+lb)/2
      do j=lco,nco
         lterms(j)=j*lc+laco(j)
         terms(j)=ten**lterms(j)
      enddo
      lsumterms=log10(sum(terms(i+1:nco)))
      maxlterms=maxval(lterms(lco:i))
      lratio=lsumterms-leps-maxlterms
      if(lratio>0)then; lb=lc
      else            ; la=lc; endif
   enddo
   arr(i)=ten**lb
enddo
arr(nco)=arr(nco-1)
arr=arr**2
end subroutine lauarr

!============================================================================
subroutine clauarr(lco,nco,eps,aco,arr)!  [lauarr]
!============================================================================
! From a given array ACO of Laurent expansion coefficients, create
! another
! array, ARR, of squared-magnitudes rr such that, at each index i, when
! the
! squared-abs-argument is no greater than arr(i), the terms
! contributed to the sum of the expansion terms at index j>i is 
! necessarily smaller than EPS-times
! the largest individual term formed taken only from ACO(lco:i).
! In other words, if eps is suitably small (such as 1.e-17) and the
! magnitude 
! of the argument is known and its square absolute value is smaller than 
! ARR(i) for some i, then
! we need only use aco(lco:i), not necessarily the whole of aco, in
! order to 
! get a numerical approximation to the function expanded with an error
! compatible with (or smaller than) the round-off error we could expect
! to
! have obtained by using the series in its entirety.
!============================================================================
integer(spi),                   intent(IN ):: lco,nco
real(dp),                       intent(IN ):: eps
complex(dpc),dimension(lco:nco),intent(IN ):: aco
real(dp),dimension(lco:nco),    intent(OUT):: arr
!-----------------------------------------------------------------------------
integer(spi),parameter     :: nit=8
real(dp),parameter         :: ten=10
real(dp),dimension(lco:nco):: absaco,laco,terms,lterms
real(dp)                   :: leps,leps4,nudge,lratio, &
                              lsumterms,maxlterms, &
                              la,lb,lc
integer(spi)               :: i,j,it
!============================================================================
leps=log10(eps); leps4=leps*4; nudge=ten**leps4
absaco=abs(aco)
laco=leps4
do i=lco,nco; if(absaco(i)>nudge)laco(i)=log10(absaco(i));enddo
lc=leps4
do i=lco,nco-1
   la=lc
   lb=1
   do it=1,nit
      lc=(la+lb)/2
      do j=lco,nco
         lterms(j)=j*lc+laco(j)
         terms(j)=ten**lterms(j)
      enddo
      lsumterms=log10(sum(terms(i+1:nco)))
      maxlterms=maxval(lterms(lco:i))
      lratio=lsumterms-leps-maxlterms
      if(lratio>0)then; lb=lc
      else            ; la=lc; endif
   enddo
   arr(i)=ten**lb
enddo
arr(nco)=arr(nco-1)
arr=arr**2
end subroutine clauarr

!============================================================================
subroutine lau2arr(lco,nco,eps,rco,qco,arr)!  [lauarr]
!============================================================================
! From a given array ACO of Laurent expansion coefficients, create
! another
! array, ARR, of squared-magnitudes rr such that, at each index i, when
! the
! squared-abs-argument is no greater than arr(i), the terms
! contributed to the sum of the expansion terms at index j>i is 
! necessarily smaller than EPS-times
! the largest individual term formed taken only from ACO(lco:i).
! In other words, if eps is suitably small (such as 1.e-17) and the
! magnitude 
! of the argument is known and its square absolute value is smaller than 
! ARR(i) for some i, then
! we need only use aco(lco:i), not necessarily the whole of aco, in
! order to 
! get a numerical approximation to the function expanded with an error
! compatible with (or smaller than) the round-off error we could expect
! to
! have obtained by using the series in its entirety.
!============================================================================
integer(spi),                        intent(IN ):: lco,nco
real(dp),                       intent(IN ):: eps
real(dp),dimension(lco:nco),    intent(IN ):: rco,qco
real(dp),dimension(lco:nco),    intent(OUT):: arr
!-----------------------------------------------------------------------------
integer(spi),parameter     :: nit=8
real(dp),parameter         :: ten=10
real(dp),dimension(lco:nco):: absaco,laco,terms,lterms
real(dp)                   :: leps,leps4,nudge,lratio, &
                              lsumterms,maxlterms, &
                              la,lb,lc
integer(spi)               :: i,j,it
!============================================================================
leps=log10(eps); leps4=leps*4; nudge=ten**leps4
absaco=abs(cmplx(rco,qco,dpc))
laco=leps4
do i=lco,nco; if(absaco(i)>nudge)laco(i)=log10(absaco(i));enddo
lc=leps4
do i=lco,nco-1
   la=lc
   lb=1
   do it=1,nit
      lc=(la+lb)/2
      do j=lco,nco
         lterms(j)=j*lc+laco(j)
         terms(j)=ten**lterms(j)
      enddo
      lsumterms=log10(sum(terms(i+1:nco)))
      maxlterms=maxval(lterms(lco:i))
      lratio=lsumterms-leps-maxlterms
      if(lratio>0)then; lb=lc
      else            ; la=lc; endif
   enddo
   arr(i)=ten**lb
enddo
arr(nco)=arr(nco-1)
arr=arr**2
end subroutine lau2arr

!=============================================================================
subroutine lausuf(lco,nco,arr,z,lb)!  [lausuf]
!=============================================================================
! find how many terms are sufficient when the complex argument of a
! Laurent
! expansion is z. The critical squared-magnitudes |z|^2 at each
! putative cut-off index lb is given in the table ARR which is
! precalculated
! in subroutine lauarr from the expansion coefficients that define the
! expansion. 
!=============================================================================
integer(spi),               intent(IN ):: lco,nco
real(dp),dimension(lco:nco),intent(IN ):: arr
complex(dpc),               intent(IN ):: z
integer(spi),               intent(OUT):: lb
!-----------------------------------------------------------------------------
integer(spi),parameter:: nit=15
integer(spi)          :: la,lc,it
real(dp)              :: ax
!=============================================================================
ax=real(z)**2+aimag(z)**2
la=lco
lb=nco
do it=1,nit
   lc=(la+lb)/2
   if(lb-la<=4)return
   if(arr(lc)<=ax)then;      la=lc
   else;                     lb=lc
   endif
enddo
end subroutine lausuf

!=============================================================================
SUBROUTINE dcinvrt(m,w,z)!  [cinvrt]
!=============================================================================
!   R.J.PURSER, NATIONAL METEOROLOGICAL CENTER, WASHINGTON D.C.  1994
!                   SUBROUTINE DCINVRT
!  compute the Taylor series coefficients z for the functional-inverse
!  of
!  the function whose taylor series coefficients are w, for the case
!  where
!  the constant coefficients of both z and w are 0
!
! --> m:    number of taylor series coefficients computed
! --> w:    taylor coefficients of original function (starting with
! linear)
! <-- z:    taylor coefficients of inverse function
!=============================================================================
INTEGER(SPI),         INTENT(IN ):: m
REAL(DP),DIMENSION(m),INTENT(IN ):: w
REAL(DP),DIMENSION(m),INTENT(OUT):: z
!-----------------------------------------------------------------------------
REAL(DP),DIMENSION(m)            :: work1,work2,work3
REAL(DP)                         :: w1,zj
INTEGER(SPI)                     :: i,j
!=============================================================================
w1=w(1)
DO i=1,m
   work1(i)=w(i)
   work2(i)=0
ENDDO
work2(1)=1
DO j=1,m
   zj=work2(j)/work1(j)
   z(j)=zj
   DO i=j,m
      work2(i)=work2(i)-zj*work1(i)
   ENDDO
   CALL conv1(m,work1,w,work3)
   DO i=1,m
      work1(i)=work3(i)
   ENDDO
ENDDO
END SUBROUTINE dcinvrt

!=============================================================================
SUBROUTINE dconv1(m,a,b,c)!  [conv1]
!=============================================================================
!   R.J.PURSER, NATIONAL METEOROLOGICAL CENTER, WASHINGTON D.C.  1994
!                   SUBROUTINE CONV1
!   convolve double-precision series a with b to form c, up to m terms
!   starting with element 1
!
! --> m:   number of elements of a, b, c
! --> a,b: inputs (convolution factors)
! <-- c:   output (convolution product)
!==========================================================================
INTEGER(SPI),         INTENT(IN ):: m
REAL(DP),DIMENSION(m),INTENT(IN ):: a,b
REAL(DP),DIMENSION(m),INTENT(OUT):: c
!-----------------------------------------------------------------------------
INTEGER(SPI)                     :: i,j,k
!=============================================================================
DO k=1,m
   c(k)=0
ENDDO
DO i=1,m
   DO j=1,m-i
      k=i+j
      c(k)=c(k)+a(i)*b(j)
   ENDDO
ENDDO
END SUBROUTINE dconv1

end module pfft1


