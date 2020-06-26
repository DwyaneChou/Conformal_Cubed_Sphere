!
!                                                         ******************
!                                                         *    pfft.f90
!                                                         *
!                                                         *  R. J.
!                                                         Purser  *
!                                                         *
!                                                         NOAA/NCEP/EMC
!                                                         *
!                                                         *      1994
!                                                         *
!                                                         ******************
! Latest revisions: April 2003
!                   April 2011 (include radix-7, 11, 13)
!                   Sept  2012 (replace "entry" statements by
!                   subroutines) 
!
! NOAA/NCEP/Environmental Modeling Center.
! jim.purser@noaa.gov
! Suite of Fast Fourier Transform routines based on factors 2, 3, 5, 7,
! 11 ("B") and 13 ("D")..
! Data assumed from 0 to n-1 and transforms performed in place.
! For complex transforms of real data, real components of wavenumber k
! are found at positions k [0,..,n/2], imag components of wavenumber -k
! are found at positions n-k [n/2+1,..n-1]. 
!
! DEPENDENCIES
! Modules: pkind, pietc
!=============================================================================
module fft23457bd
!=============================================================================
use pkind, only: spi
implicit none
integer(spi) ln2,ln3,ln4,ln5,ln7,lnb,lnd,ln,n,nm,nh
end module fft23457bd

!=============================================================================
module pfft
!=============================================================================
use pkind, only: spi,sp,dp
implicit none
private
public :: fftco,fftcop,cfft,dfft,csft,dsft,dsfe,hsfe, &
          rfft,hfft,fftcnv,fconv

interface snfftln;    module procedure infftln;                  end interface
interface sget2357bd; module procedure iget2357bd;               end interface
interface get23457bd; module procedure iget23457bd;              end interface
interface rumble;     module procedure srumble,drumble;          end interface
interface fftco;      module procedure sfftco, dfftco;           end interface
interface fftcop;     module procedure sfftcop, dfftcop;         end interface
interface cfft
   module procedure scfft,  dcfft, scfftp, dcfftp, zcfft;        end interface
interface dfft
   module procedure sdfft,  ddfft, sdfftp, ddfftp, zdfft;        end interface
interface gfft;     module procedure sgfft,  dgfft;              end interface
interface csft;     module procedure scsft,  dcsft;              end interface
interface dsft;     module procedure sdsft,  ddsft;              end interface
interface dsfe;     module procedure sdsfe,  ddsfe;              end interface
interface hsfe;     module procedure shsfe,  dhsfe;              end interface
interface rfft     
   module procedure srfft,  drfft, srfftp, drfftp;               end interface
interface hfft
   module procedure shfft,  dhfft, shfftp, dhfftp;               end interface
interface fftcnv
   module procedure sfftcnv,dfftcnv,sfftcnvp, dfftcnvp;          end interface
interface fconv;    module procedure sfconv, dfconv;             end interface

contains

!=============================================================================
subroutine infftln(n1,m,nfftln)!  [snfftln]
!=============================================================================
! Divide out as many factors, m, of given n1 as is possible, returning
! this number of factors as nfftln, and returning n1 as the quotient.
!=============================================================================
integer(spi),intent(INOUT):: n1
integer(spi),intent(IN   ):: m
integer(spi),intent(  OUT):: nfftln
!-----------------------------------------------------------------------------
integer(spi)              :: i,n2
!=============================================================================
nfftln=0
do i=1,30
   n2=n1/m
   if(n2*m /= n1)return
   n1=n2
   nfftln=i
enddo
end subroutine infftln 

!=============================================================================
subroutine iget2357bd(nin,get2357bd)!  [sget2357bd]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!
!   Factorize NIN in terms of 2's, 3's, 4's, 5's, 7's, 11's and 13's
!   and 
!   verify that these include the only prime factors. Set LN2, LN3, LN4,
!   LN5, 
!   LN7, LNB, LND respectively to the number of powers of: 2, 3, 4, 5,
!   7, 
!   11, 13; set LN to be the sum, 
!   LN2+LN3+LN4*2+LN5+LN7+LNB+LND, and initialize n, nm, nh, 
!   all in fft23457bd.
!
! --> NIN      Number of data along the line of Fourier transformation
!=============================================================================
use fft23457bd
logical,     intent(OUT):: get2357bd
integer(spi),intent(IN ):: nin
!-----------------------------------------------------------------------------
integer(spi)            ::  k
!=============================================================================
k=nin
call snfftln(k,13,lnd)
call snfftln(k,11,lnb)
call snfftln(k,7,ln7)
call snfftln(k,5,ln5)
call snfftln(k,4,ln4)
call snfftln(k,3,ln3)
call snfftln(k,2,ln2)
ln=ln2+ln4*2+ln3+ln5+ln7+lnb+lnd
n=(2**ln2)*(3**ln3)*(4**ln4)*(5**ln5)*(7**ln7)*(11**lnb)*(13**lnd)
nm=n-1; nh=n/2
if(n /= nin)stop 'prime factors of fft period are not only 2, 3, 5, 7, 11, 13'
get2357bd=(k==1)
end subroutine iget2357bd

!=============================================================================
subroutine iget23457bd(pl2,pl3,pl4,pl5,pl7,plb,pld)!  [get23457bd]
!=============================================================================
use fft23457bd
integer(spi), intent(OUT):: pl2,pl3,pl4,pl5,pl7,plb,pld
pl2=ln2; pl3=ln3; pl4=ln4; pl5=ln5; pl7=ln7; plb=lnb; pld=lnd
end subroutine iget23457bd

!=============================================================================
subroutine srumble(jumble,tumble)!  [rumble] 
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!
!   Initialize coefficient arrays JUMBLE and TUMBLE for use with the
!   fast-Fourier-transform routines when the number N of data has the
!   prime-
!   factorization 2**LN2*3**LN3*5**LN5*7**LN7*11**LNB*13**LND
!
! <-- JUMBLE:  Permutation of N data encoded as a sequence of
! transpositions
! <-- TUMBLE:  Trigonometric coefficients for use in FFT. The first half
! are
!              the cosines, the second half the sines, of uniformly
!              increasing
!              relevant angles.
!=============================================================================
use pietc_s, only: pi2
use fft23457bd
integer(spi), dimension(0:*),intent(OUT):: jumble
real(sp),     dimension(0:*),intent(OUT):: tumble
!-----------------------------------------------------------------------------
integer(spi), parameter    :: ml=30
integer(spi)               :: i,j,l,id,is,ir,kd
integer(spi), dimension(ml):: nd,md
real(sp)                   :: ang,pi2on
!=============================================================================
pi2on=pi2/n
do i=0,nh-1
   ang=pi2on*i
   tumble(i)   =cos(ang); tumble(i+nh)=sin(ang)
enddo
id=1;  is=0
do i=1,lnd;       is=is+1; md(is)=id; id=id*13; enddo
do i=1,lnb;       is=is+1; md(is)=id; id=id*11; enddo
do i=1,ln7;       is=is+1; md(is)=id; id=id*7;  enddo
do i=1,ln5;       is=is+1; md(is)=id; id=id*5;  enddo
do i=1,ln3;       is=is+1; md(is)=id; id=id*3;  enddo
do i=1,ln2+ln4*2; is=is+1; md(is)=id; id=id*2;  enddo
id=1
do i=1,ln2+ln4*2; nd(is)=id; id=id*2;  is=is-1; enddo
do i=1,ln3;       nd(is)=id; id=id*3;  is=is-1; enddo
do i=1,ln5;       nd(is)=id; id=id*5;  is=is-1; enddo
do i=1,ln7;       nd(is)=id; id=id*7;  is=is-1; enddo
do i=1,lnb;       nd(is)=id; id=id*11; is=is-1; enddo
do i=1,lnd;       nd(is)=id; id=id*13; is=is-1; enddo
jumble(0)=n
do i=1,nm
   ir=i; j=0
   do l=1,ln; kd=ir/nd(l); ir=ir-kd*nd(l); j=j+kd*md(l); enddo
   jumble(i)=j
enddo
do i=1,nm
   j=jumble(i)
   do
      if(j >= i)exit
      j=jumble(j)
   enddo
   jumble(i)=j
enddo
end subroutine srumble 

!=============================================================================
subroutine drumble(jumble,tumble)!  [rumble]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!=============================================================================
use pietc, only: pi2
use fft23457bd
integer(spi), dimension(0:*),intent(OUT):: jumble
real(dp),     dimension(0:*),intent(OUT):: tumble
!-----------------------------------------------------------------------------
integer(spi), parameter    :: ml=30
integer(spi)               :: i,j,l,id,is,ir,kd
integer(spi), dimension(ml):: nd,md
real(dp)                   :: ang,pi2on
!=============================================================================
pi2on=pi2/n
do i=0,nh-1
   ang=pi2on*i
   tumble(i)   =cos(ang); tumble(i+nh)=sin(ang)
enddo
id=1;  is=0
do i=1,lnd;       is=is+1; md(is)=id; id=id*13; enddo
do i=1,lnb;       is=is+1; md(is)=id; id=id*11; enddo
do i=1,ln7;       is=is+1; md(is)=id; id=id*7;  enddo
do i=1,ln5;       is=is+1; md(is)=id; id=id*5;  enddo
do i=1,ln3;       is=is+1; md(is)=id; id=id*3;  enddo
do i=1,ln2+ln4*2; is=is+1; md(is)=id; id=id*2;  enddo
id=1
do i=1,ln2+ln4*2; nd(is)=id; id=id*2;  is=is-1; enddo
do i=1,ln3;       nd(is)=id; id=id*3;  is=is-1; enddo
do i=1,ln5;       nd(is)=id; id=id*5;  is=is-1; enddo
do i=1,ln7;       nd(is)=id; id=id*7;  is=is-1; enddo
do i=1,lnb;       nd(is)=id; id=id*11; is=is-1; enddo
do i=1,lnd;       nd(is)=id; id=id*13; is=is-1; enddo
jumble(0)=n
do i=1,nm
   ir=i; j=0
   do l=1,ln; kd=ir/nd(l); ir=ir-kd*nd(l); j=j+kd*md(l); enddo
   jumble(i)=j
enddo
do i=1,nm
   j=jumble(i)
   do
      if(j >= i)exit
      j=jumble(j)
   enddo
   jumble(i)=j
enddo

end subroutine drumble

!=============================================================================
subroutine sfftco(n,j,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd)!  [fftco] 
!=============================================================================
!   R. J. Purser, NCEP, Washington D.C. 1999. jim.purser@noaa.gov
!                   SUBROUTINE SFFTCO
!   Provide the FFT integer coefficients, j, and real coefficients, w,
! corresponding to a period, n, which factors into powers of 2, 3, 5, 7,
! 11, 13
! only. For short period transforms, retain (in js, ws) up to two
! previously 
! used sets of these coefficients to avoid the cost of recalculating
! them
! in repeated applications. Then the first element of the relevant
! column of js is always the period, n, making it easy to recognize
! whether or not the new requirements for j and w are met in previously
! recorded pairs, js and ws. 
!   For long period transforms, it is worth supplying an integer array,
!   js,
! and a real array, ws, both of size n, to store the transform
! parameters
! between applications so that they do not have to be recomputed each
! time.
! This is found to save significant time in large scale repeated
! applications.
! In that case, the relevant version of this routine is SFFTCOP.
!
! --> n:      period of data
! --> j:      array of n permutation indices to unscramble the fft
! output.
! --> w:      array of n real trigonometric coefficients for the fft.
!=============================================================================
use pietc_s, only: u0
integer(spi), parameter :: bsize=2048,bsizem=bsize-1,bsize2=bsize*2
integer(spi),                 intent(IN) :: n
integer(spi),dimension(0:n-1),intent(OUT):: j
real(sp),    dimension(0:n-1),intent(OUT):: w
integer(spi),                 intent(OUT):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
!-----------------------------------------------------------------------------
logical                           :: get2357bd
integer(spi)                      :: ult,nm
integer(spi),dimension(0:bsizem,2):: js
real(sp),    dimension(0:bsizem,2):: ws
data ult/1/ ! Column index of js and ws for latest fft coefficients used.
data js/bsize2*0/,ws/bsize2*u0/
!=============================================================================
nm=n-1
call sget2357bd(n,get2357bd)
if(.not. get2357bd)stop 'prime factors are not only 2, 3, 5, 7, 11, 13'
call get23457bd(ln2,ln3,ln4,ln5,ln7,lnb,lnd)
if(n <= bsize)then
   if(n == js(0,ult))then
      j(0:nm)=js(0:nm,ult); w(0:nm)=ws(0:nm,ult) ! copy existing values
      return
   endif
   ult=3-ult                             ! Try the alternative location
   if(n == js(0,ult))then
      j(0:nm)=js(0:nm,ult); w(0:nm)=ws(0:nm,ult) ! copy existing values
      return
   endif
   call rumble(j,w)
   js(0:nm,ult)=j(0:nm); ws(0:nm,ult)=w(0:nm)
else
   call rumble(j,w) ! Too big to record in js and ws
endif
end subroutine sfftco

!=============================================================================
subroutine dfftco(n,j,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd)!  [fftco]
!=============================================================================
!   R. J. Purser, NCEP, Washington D.C. 1999. jim.purser@noaa.gov
!                   SUBROUTINE DFFTCO: double precision version of fftco
!=============================================================================
use pietc, only: u0
integer(spi), parameter:: bsize=2048,bsizem=bsize-1,bsize2=bsize*2
integer(spi),                 intent(IN) :: n
integer(spi),dimension(0:n-1),intent(OUT):: j
real(dp),    dimension(0:n-1),intent(OUT):: w
integer(spi),                 intent(OUT):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
!-----------------------------------------------------------------------------
logical                            :: get2357bd
integer(spi)                       :: ult,nm
integer(spi), dimension(0:bsizem,2):: js
real(dp),dimension(0:bsizem,2)     :: ws
data ult/1/ ! Column index of js and ws for latest fft coefficients used.
data js/bsize2*0/,ws/bsize2*u0/
!=============================================================================
nm=n-1
call sget2357bd(n,get2357bd)
if(.not. get2357bd)stop 'prime factors are not only 2, 3, 5, 7, 11, 13'
call get23457bd(ln2,ln3,ln4,ln5,ln7,lnb,lnd)
if(n <= bsize)then
   if(n == js(0,ult))then
      j(0:nm)=js(0:nm,ult); w(0:nm)=ws(0:nm,ult) ! copy existing values
      return
   endif
   ult=3-ult                             ! Try the alternative location
   if(n == js(0,ult))then
      j(0:nm)=js(0:nm,ult); w(0:nm)=ws(0:nm,ult) ! copy existing values
      return
   endif
   call rumble(j,w)
   js(0:nm,ult)=j(0:nm); ws(0:nm,ult)=w(0:nm)
else
   call rumble(j,w) ! Too big to record in js and ws
endif
end subroutine dfftco

!=============================================================================
subroutine sfftcop(n,js,ws,ln2,ln3,ln4,ln5,ln7,lnb,lnd)!  [fftcop]
!=============================================================================
integer(spi),                   intent(IN   ):: n
integer(spi),dimension(0:n-1),  intent(INOUT):: js
real(sp),    dimension(0:n-1),  intent(INOUT):: ws
integer(spi),                   intent(  OUT):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
!-----------------------------------------------------------------------------
logical:: get2357bd
!=============================================================================
call sget2357bd(n,get2357bd)
if(.not. get2357bd)stop 'prime factors are not only 2, 3, 5, 7, 11, 13'
call get23457bd(ln2,ln3,ln4,ln5,ln7,lnb,lnd)
if(n /= js(0))call rumble(js,ws)
end subroutine sfftcop
!=============================================================================
subroutine dfftcop(n,js,ws,ln2,ln3,ln4,ln5,ln7,lnb,lnd)!  [fftcop]
!=============================================================================
!   R. J. Purser, NCEP, Washington D.C. 1999. jim.purser@noaa.gov
!  double precision version of fftcop with the saved parameter arrays
!  provided through the parameter list. This allows the user to provide
!  these parameter arrays of the appropriate size for the duration of
!  the
!  fft tasks and does not require permanent storage to be set aside for
!  them.
!  This is advantageous when the FFT problem size might be very large,
!  or
!  where many different-sized FFT tasks are alternated repetitively
!  (which
!  would necessitate redundant recalculations of these parameters is
!  only
!  a fixed size pair of work arrays were reserved).
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(dp),    dimension(0:n-1),intent(INOUT):: ws
integer(spi),                 intent(  OUT):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
!-----------------------------------------------------------------------------
logical:: get2357bd
!=============================================================================
call sget2357bd(n,get2357bd)
if(.not. get2357bd)stop 'prime factors are not only 2, 3, 5, 7, 11, 13'
call get23457bd(ln2,ln3,ln4,ln5,ln7,lnb,lnd)
if(n /= js(0))call rumble(js,ws)
end subroutine dfftcop

!=============================================================================
subroutine zcfft(n,z)!  [cfft]
!=============================================================================
! Like cfft, but with r and q combined into complex z.
!=============================================================================
use pkind,only: dpc
integer(spi),                 intent(IN   ):: n
complex(dpc),dimension(0:n-1),intent(INOUT):: z
!-----------------------------------------------------------------------------
real(dp),dimension(0:n-1):: r,q
!=============================================================================
r=real(z); q=aimag(z)
call cfft(n,r,q)
z=cmplx(r,q,dpc)
end subroutine zcfft

!=============================================================================
subroutine zdfft(n,z)!  [dfft]
!=============================================================================
! Like dfft, but with r and q combined into complex z.
!=============================================================================
use pkind,only: dp,dpc
integer(spi),                 intent(IN   ):: n
complex(dpc),dimension(0:n-1),intent(INOUT):: z
!-----------------------------------------------------------------------------
real(dp),dimension(0:n-1):: r,q
!=============================================================================
r=real(z); q=aimag(z)
call dfft(n,r,q)
z=cmplx(r,q,dpc)
end subroutine zdfft

!=============================================================================
subroutine sdfft(n,r,q)!  [dfft]
!=============================================================================
use pietc, only: u1
integer(spi),             intent(IN   ):: n
real(sp),dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi):: nm
real(sp)    :: rfac
!=============================================================================
rfac=u1/n; nm=n-1
!  for fourier synthesis, scale, and reverse the order of wavenumbers:
r(0)=r(0)*rfac; r(1:nm)=r(nm:1:-1)*rfac
q(0)=q(0)*rfac; q(1:nm)=q(nm:1:-1)*rfac
call scfft(n,r,q)
end subroutine sdfft

!=============================================================================
subroutine scfft(n,r,q)!  [cfft]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
! (Revised for f90, 1999; radix-7, 11, 13 included, 2011)
!=============================================================================
integer(spi),             intent(IN   ):: n
real(sp),dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi)                  :: ln2,ln3,ln4,ln5,ln7,lnb,lnd
integer(spi),dimension(0:n-1) :: jumble
real(sp),    dimension(0:n-1) :: w
!=============================================================================
call fftco(n,jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd)
call  gfft(n,jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd,r,q)
end subroutine scfft

!=============================================================================
subroutine ddfft(n,r,q)!  [dfft]
!=============================================================================
use pietc, only: u1
integer(spi),             intent(IN   ):: n
real(dp),dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi):: nm
real(dp)    :: rfac
!=============================================================================
rfac=u1/n; nm=n-1
!  for fourier synthesis, scale, and reverse the order of wavenumbers:
r(0)=r(0)*rfac; r(1:nm)=r(nm:1:-1)*rfac
q(0)=q(0)*rfac; q(1:nm)=q(nm:1:-1)*rfac
call dcfft(n,r,q)
end subroutine ddfft

!=============================================================================
subroutine dcfft(n,r,q)!  [cfft]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
! (Revised for f90, 1999; radix-7, 11, 13 included, 2011)
!                   SUBROUTINES DCFFT, DDFFT
! Double precision versions of cfft, dfft
!=============================================================================
integer(spi),             intent(IN   ):: n
real(dp),dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi)                 :: ln2,ln3,ln4,ln5,ln7,lnb,lnd
integer(spi),dimension(0:n-1):: jumble
real(dp),    dimension(0:n-1):: w
!=============================================================================
call fftco(n,jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd)
call  gfft(n,jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd,r,q)
end subroutine dcfft

!=============================================================================
subroutine sdfftp(n,js,ws,r,q)!  [dfft]
!=============================================================================
use pietc_s, only: u1
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(sp),    dimension(0:n-1),intent(INOUT):: ws
real(sp),    dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi):: nm
real(sp)    :: rfac
!=============================================================================
rfac=u1/n; nm=n-1
!  for fourier synthesis, scale, and reverse the order of wavenumbers:
r(0)=r(0)*rfac; r(1:nm)=r(nm:1:-1)*rfac
q(0)=q(0)*rfac; q(1:nm)=q(nm:1:-1)*rfac
call scfftp(n,js,ws,r,q)
end subroutine sdfftp
!=============================================================================
subroutine ddfftp(n,js,ws,r,q)!  [dfft]
!=============================================================================
use pietc, only: u1
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(dp),    dimension(0:n-1),intent(INOUT):: ws
real(dp),    dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi):: nm
real(dp)    :: rfac
!=============================================================================
rfac=u1/n; nm=n-1
!  for fourier synthesis, scale, and reverse the order of wavenumbers:
r(0)=r(0)*rfac; r(1:nm)=r(nm:1:-1)*rfac
q(0)=q(0)*rfac; q(1:nm)=q(nm:1:-1)*rfac
call dcfftp(n,js,ws,r,q)
end subroutine ddfftp

!=============================================================================
subroutine scfftp(n,js,ws,r,q)!  [cfft]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
! (Revised for f90, 1999; radix-7, 11, 13, included and user-provided
! parameter
!  arrays an option, 2011)
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(sp),    dimension(0:n-1),intent(INOUT):: ws
real(sp),    dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
!=============================================================================
call fftcop(n,js,ws,ln2,ln3,ln4,ln5,ln7,lnb,lnd)
call   gfft(n,js,ws,ln2,ln3,ln4,ln5,ln7,lnb,lnd,r,q)
end subroutine scfftp
!=============================================================================
subroutine dcfftp(n,js,ws,r,q)!  [cfft]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
! (Revised for f90, 1999; radix-7, 11, 13, included and user-provided
! parameter
!  arrays an option, 2011)
! Double precision versions of cfft, dfft
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(dp),    dimension(0:n-1),intent(INOUT):: ws
real(dp),    dimension(0:n-1),intent(INOUT):: r,q
!------------------------------------------------------------------------------
integer(spi):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
!=============================================================================
call fftcop(n,js,ws,ln2,ln3,ln4,ln5,ln7,lnb,lnd)
call   gfft(n,js,ws,ln2,ln3,ln4,ln5,ln7,lnb,lnd,r,q)
end subroutine dcfftp

!=============================================================================
subroutine sgfft(n,jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd,r,q)!  [gfft]
!=============================================================================
use pietc_s, only: u1,u3,s30,ms30,s60
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(IN   ):: jumble
real(sp),    dimension(0:n-1),intent(IN   ):: w
integer(spi),                 intent(IN   ):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
real(sp),    dimension(0:n-1),intent(INOUT):: r,q
!-----------------------------------------------------------------------------
integer(spi):: nm,nh,i,j,l,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc, &
               ma,ma2,ma3,ma4,ma5,ma7,mab,mad,mb,mb2, &
               jmb,jmb2,jmb3,jmb4,jmb5,jmb6, &
               nze
real(sp)    :: r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,ra,rb,rc, &
               q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,qa,qb,qc,  &
               rf1,rf2,rf3,rf4,rf5,rf6,rf7,rf8,rf9,rfa,rfb,rfc,&
               qf1,qf2,qf3,qf4,qf5,qf6,qf7,qf8,qf9,qfa,qfb,qfc, &
               rep,rze,rzc,ret,rec,rth,rtc,rio,ric,rka,rkc,rla,rlc, &
               qep,qze,qet,qth,qio,qka,qla, &
               t1,t2,t3,t4,t5,t6
!=============================================================================
nm=n-1
nh=n/2  
! PERMUTE THE DATA:
do i=1,nm
   j=jumble(i)
   if(j > i)then
      t1=r(i); r(i)=r(j); r(j)=t1
      t1=q(i); q(i)=q(j); q(j)=t1
   endif
enddo

!  TRANSFORM THE DATA:
ma=1; mb=n
!  RADIX 4
do l=1,ln4
   mb=mb/4; ma4=ma*4; mb2=mb*2
   do j=0,ma-1
      jmb=j*mb; jmb2=j*mb2
      rf1=w(jmb);          qf1=w(nh+jmb)        ! f**1
      rf2=w(jmb2);         qf2=w(nh+jmb2)       ! f**2
      rf3=rf1*rf2-qf1*qf2; qf3=rf1*qf2+qf1*rf2  ! f**3
      do i=0,nm,ma4
         k0=i+j;    k1=k0+ma;   k2=k1+ma;   k3=k2+ma
         r0=r(k0); r1=r(k1);  r2=r(k2);  r3=r(k3)
         q0=q(k0); q1=q(k1);  q2=q(k2);  q3=q(k3)
         t1    =r3*rf3-q3*qf3 ! q13
         q3    =r3*qf3+q3*rf3 ! r13
         r3    =r2*rf1-q2*qf1 ! r12
         r2    =r2*qf1+q2*rf1 ! q12
         q2    =q3    -r2     ! r23
         r2    =q3    +r2     ! q22
         q3    =r3    +t1     ! r22
         t1    =r3    -t1     ! q23
         r3    =r1*rf2-q1*qf2 ! r11
         q1    =r1*qf2+q1*rf2 ! q11
         r1    =r0    -r3     ! r21
         r0    =r0    +r3     ! r20
         r(k3) =r1    -q2     ! r43
         r(k1) =r1    +q2     ! r41
         q2    =q0    +q1     ! q20
         q1    =q0    -q1     ! q21
         q(k0) =q2    +r2     ! q40
         q(k2) =q2    -r2     ! q42
         r(k2) =r0    -q3     ! r42
         r(k0) =r0    +q3     ! r40
         q(k3) =q1    -t1     ! q43
         q(k1) =q1    +t1     ! q41
      enddo
   enddo
   ma=ma4
enddo
if(ln2==1)then
!  RADIX 2
   mb=mb/2; ma2=ma*2
   do j=0,ma-1
      jmb=j*mb
      do i=0,nm,ma2
         k0=j+i;        k1=k0+ma
         rf1=w(jmb);    qf1=w(nh+jmb) ! f**1
         r0=r(k0);     q0=q(k0)
         r1=r(k1);     q1=q(k1)
         t1   =r1*qf1+q1*rf1 ! q11
         q1   =r1*rf1-q1*qf1 ! r11
         r(k1)=r0    -q1     ! r41
         r(k0)=r0    +q1     ! r40
         q(k1)=q0    -t1     ! q41
         q(k0)=q0    +t1     ! q40
      enddo
   enddo
   ma=ma2
endif
!  RADIX 3
rep=ms30;  rec=u3*s30;  qep=s60 ! <- epsilon
do l=1,ln3
   mb=mb/3; ma3=ma*3
   do j=0,ma-1
      jmb=j*mb
      rf1=w(jmb);          qf1=w(nh+jmb) ! f**1
      rf2=rf1*rf1-qf1*qf1; qf2=2*rf1*qf1 ! f**2
      do i=0,nm,ma3
         k0=i+j;     k1=k0+ma;    k2=k1+ma
         r1=r(k1);  q1=q(k1)
         r2=r(k2);  q2=q(k2)
         t1    = r2*qf2+q2 *rf2 ! r12
         q2    = r2*rf2-q2 *qf2 ! q12
         r2    = r1*qf1+q1 *rf1 ! q11
         r1    = r1*rf1-q1 *qf1 ! r11
         q1    = r2    +t1      ! q21
         r2    =(r2    -t1)*qep ! r22
         t1    = r1    +q2      ! r21
         r1    =(r1    -q2)*qep ! q22
         r(k0) = r(k0)+t1       ! r40
         q(k0) = q(k0)+q1       ! q40
         t1    = r(k0)-t1 *rec  ! r21
         q1    = q(k0)-q1 *rec  ! q21
         q(k2) = q1    -r1      ! q42
         q(k1) = q1    +r1      ! q41
         r(k1) = t1    -r2      ! r41
         r(k2) = t1    +r2      ! r42
      enddo
   enddo
   ma=ma3
enddo

if(ln5 > 0)then
!  RADIX 5
   nze=n/5;  rze=w(nze); qze=w(nh+nze); rzc=u1-rze ! <- zeta
   ret=rze*rze-qze*qze;  qet=2*rze*qze; rec=u1-ret ! <- eta
   do l=1,ln5
      mb=mb/5; ma5=ma*5
      do j=0,ma-1
         jmb=j*mb;             jmb2=jmb*2
         rf1=w(jmb);           qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);          qf2=w(nh+jmb2)      ! f**2
         rf3=rf1*rf2-qf1*qf2;  qf3=rf1*qf2+qf1*rf2 ! f**3
         rf4=rf2*rf2-qf2*qf2;  qf4=2*rf2*qf2       ! f**4
         do i=0,nm,ma5
            k0=i+j;   k1=k0+ma;  k2=k1+ma;  k3=k2+ma;  k4=k3+ma
            r1=r(k1); r2=r(k2); r3=r(k3); r4=r(k4)
            q1=q(k1); q2=q(k2); q3=q(k3); q4=q(k4)
            t1    =r1*qf1+q1*rf1       ! q11
            r1    =r1*rf1-q1*qf1       ! r11
            q1    =r4*rf4-q4*qf4       ! r14
            r4    =r4*qf4+q4*rf4       ! q14
            q4    =r1    -q1           ! q24
            r1    =r1    +q1           ! r21
            q1    =t1    +r4           ! q21
            r4    =t1    -r4           ! q24
            t1    =r3*rf3-q3*qf3       ! r13
            r3    =r3*qf3+q3*rf3       ! q13
            q3    =r2*qf2+q2*rf2       ! q12
            r2    =r2*rf2-q2*qf2       ! r12
            q2    =q3    +r3           ! q22
            r3    =q3    -r3           ! q23
            q3    =r2    -t1           ! r23
            r2    =r2    +t1           ! r22
            r(k0) =r(k0)+r1    +r2     ! r40
            q(k0) =q(k0)+q1    +q2     ! q40
            t1    =r4*qze+r3*qet       ! q34
            r3    =r4*qet-r3*qze       ! q33
            r4    =r(k0)-r1*rec-r2*rzc ! r32
            r1    =r(k0)-r1*rzc-r2*rec ! r31
            r(k2) =r4    -r3           ! r42
            r(k3) =r4    +r3           ! r43
            r(k4) =r1    +t1           ! r44
            r(k1) =r1    -t1           ! r41
            t1    =q(k0)-q1*rzc-q2*rec ! q31
            q2    =q(k0)-q1*rec-q2*rzc ! q32
            q1    =q4*qet-q3*qze       ! r33
            q4    =q4*qze+q3*qet       ! r34
            q(k3) =q2    -q1           ! q43
            q(k2) =q2    +q1           ! q42
            q(k1) =t1    +q4           ! q41
            q(k4) =t1    -q4           ! q44
         enddo
      enddo
      ma=ma5
   enddo
endif

if(ln7 > 0)then
!  RADIX 7
   nze=n/7
   rze=w(nze);          qze=w(nh+nze);       rzc=u1-rze ! <- zeta 
   ret=rze*rze-qze*qze; qet=2*rze*qze;       rec=u1-ret ! <- eta 
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=u1-rth ! <- theta 
   do L=1,Ln7
      mb=mb/7; ma7=ma*7
      do j=0,ma-1
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3
         rf1=w(jmb);          qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);         qf2=w(nh+jmb2)      ! f**2
         rf3=w(jmb3);         qf3=w(nh+jmb3)      ! f**3
         rf4=rf2*rf2-qf2*qf2; qf4=2*rf2*qf2       ! f**4
         rf5=rf2*rf3-qf2*qf3; qf5=rf2*qf3+qf2*rf3 ! f**5
         rf6=rf3*rf3-qf3*qf3; qf6=2*rf3*qf3       ! f**6
         do i=0,nm,ma7
            k0=i+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            r0=r(k0)
            r1=r(k1);r2=r(k2);r3=r(k3);r4=r(k4);r5=r(k5);r6=r(k6)
            q0=q(k0)
            q1=q(k1);q2=q(k2);q3=q(k3);q4=q(k4);q5=q(k5);q6=q(k6)
            t1=r1*qf1+q1*rf1           ! q11
            r1=r1*rf1-q1*qf1           ! r11
            q1=r6*rf6-q6*qf6           ! r16
            r6=r6*qf6+q6*rf6           ! q16
            q6=r1-q1                   ! r26
            r1=r1+q1                   ! r21
            q1=t1+r6                   ! q21
            r6=t1-r6                   ! q26
            t2=r2*qf2+q2*rf2           ! q12
            r2=r2*rf2-q2*qf2           ! r12
            q2=r5*rf5-q5*qf5           ! r15
            r5=r5*qf5+q5*rf5           ! q15
            q5=r2-q2                   ! r25
            r2=r2+q2                   ! r22
            q2=t2+r5                   ! q22
            r5=t2-r5                   ! q25
            t3=r3*qf3+q3*rf3           ! q13
            r3=r3*rf3-q3*qf3           ! r13
            q3=r4*rf4-q4*qf4           ! r14
            r4=r4*qf4+q4*rf4           ! q14
            q4=r3-q3                   ! r24
            r3=r3+q3                   ! r23
            q3=t3+r4                   ! q23
            r4=t3-r4                   ! q24
            r0=r0+r1+r2+r3             ! r40
            q0=q0+q1+q2+q3             ! q40
            t1=   r6*qze+r5*qet+r4*qth ! q36
            t2=   r6*qet-r5*qth-r4*qze ! q35
            t3=   r6*qth-r5*qze+r4*qet ! q34
            r6=r0-r1*rzc-r2*rec-r3*rtc ! r31
            r5=r0-r1*rec-r2*rtc-r3*rzc ! r32
            r4=r0-r1*rtc-r2*rzc-r3*rec ! r33
            r(k0)=r0                   ! r40
            r(k1)=r6-t1                ! r41
            r(k6)=r6+t1                ! r46
            r(k2)=r5-t2                ! r42
            r(k5)=r5+t2                ! r45
            r(k3)=r4-t3                ! r43
            r(k4)=r4+t3                ! r44
            t1=   q6*qze+q5*qet+q4*qth ! r36
            t2=   q6*qet-q5*qth-q4*qze ! r35
            t3=   q6*qth-q5*qze+q4*qet ! r34
            q6=q0-q1*rzc-q2*rec-q3*rtc ! q31
            q5=q0-q1*rec-q2*rtc-q3*rzc ! q32
            q4=q0-q1*rtc-q2*rzc-q3*rec ! q33
            q(k0)=q0                   ! q40
            q(k1)=q6+t1                ! q41
            q(k6)=q6-t1                ! q46
            q(k2)=q5+t2                ! q42
            q(k5)=q5-t2                ! q45
            q(k3)=q4+t3                ! q43
            q(k4)=q4-t3                ! q44
         enddo
      enddo
      ma=ma7
   enddo
endif
if(lnb > 0)then
!  RADIX 11
   nze=n/11
   rze=w(nze);          qze=w(nh+nze);       rzc=u1-rze ! <- zeta
   ret=rze*rze-qze*qze; qet=2*rze*qze;       rec=u1-ret ! <- eta 
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=u1-rth ! <- theta
   rio=ret*ret-qet*qet; qio=2*ret*qet;       ric=u1-rio ! <- iota
   rka=ret*rth-qet*qth; qka=ret*qth+qet*rth; rkc=u1-rka ! <- kappa
   do L=1,Lnb
      mb=mb/11; mab=ma*11
      do j=0,ma-1
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3; jmb4=jmb*4; jmb5=jmb*5
         rf1=w(jmb);          qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);         qf2=w(nh+jmb2)      ! f**2
         rf3=w(jmb3);         qf3=w(nh+jmb3)      ! f**3
         rf4=w(jmb4);         qf4=w(nh+jmb4)      ! f**4
         rf5=w(jmb5);         qf5=w(nh+jmb5)      ! f**5
         rf6=rf3*rf3-qf3*qf3; qf6=2*rf3*qf3       ! f**6
         rf7=rf3*rf4-qf3*qf4; qf7=rf3*qf4+qf3*rf4 ! f**7
         rf8=rf4*rf4-qf4*qf4; qf8=2*rf4*qf4       ! f**8
         rf9=rf4*rf5-qf4*qf5; qf9=rf4*qf5+qf4*rf5 ! f**9
         rfa=rf5*rf5-qf5*qf5; qfa=2*rf5*qf5       ! f**10
         do i=0,nm,mab
            k0=i+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            k7=k6+ma; k8=k7+ma; k9=k8+ma; ka=k9+ma
            r0=r(k0)
            r1=r(k1);r2=r(k2);r3=r(k3);r4=r(k4);r5=r(k5);r6=r(k6)
            r7=r(k7);r8=r(k8);r9=r(k9);ra=r(ka)
            q0=q(k0)
            q1=q(k1);q2=q(k2);q3=q(k3);q4=q(k4);q5=q(k5);q6=q(k6)
            q7=q(k7);q8=q(k8);q9=q(k9);qa=q(ka)
            t1=r1*qf1+q1*rf1           ! q11
            r1=r1*rf1-q1*qf1           ! r11
            q1=ra*rfa-qa*qfa           ! r1a
            ra=ra*qfa+qa*rfa           ! q1a
            qa=r1-q1                   ! r2a
            r1=r1+q1                   ! r21
            q1=t1+ra                   ! q21
            ra=t1-ra                   ! q2a
            t2=r2*qf2+q2*rf2           ! q12
            r2=r2*rf2-q2*qf2           ! r12
            q2=r9*rf9-q9*qf9           ! r19
            r9=r9*qf9+q9*rf9           ! q19
            q9=r2-q2                   ! r29
            r2=r2+q2                   ! r22
            q2=t2+r9                   ! q22
            r9=t2-r9                   ! q29
            t3=r3*qf3+q3*rf3           ! q13
            r3=r3*rf3-q3*qf3           ! r13
            q3=r8*rf8-q8*qf8           ! r18
            r8=r8*qf8+q8*rf8           ! q18
            q8=r3-q3                   ! r28
            r3=r3+q3                   ! r23
            q3=t3+r8                   ! q23
            r8=t3-r8                   ! q28
            t4=r4*qf4+q4*rf4           ! q14
            r4=r4*rf4-q4*qf4           ! r14
            q4=r7*rf7-q7*qf7           ! r17
            r7=r7*qf7+q7*rf7           ! q17
            q7=r4-q4                   ! r27
            r4=r4+q4                   ! r24
            q4=t4+r7                   ! q24
            r7=t4-r7                   ! q27
            t5=r5*qf5+q5*rf5           ! q15
            r5=r5*rf5-q5*qf5           ! r15
            q5=r6*rf6-q6*qf6           ! r16
            r6=r6*qf6+q6*rf6           ! q16
            q6=r5-q5                   ! r26
            r5=r5+q5                   ! r25
            q5=t5+r6                   ! q25
            r6=t5-r6                   ! q26
            r0=r0+r1+r2+r3+r4+r5       ! r40
            q0=q0+q1+q2+q3+q4+q5       ! q40
            t1=   ra*qze+r9*qet+r8*qth+r7*qio+r6*qka ! q3a
            t2=   ra*qet+r9*qio-r8*qka-r7*qth-r6*qze ! q39
            t3=   ra*qth-r9*qka-r8*qet+r7*qze+r6*qio ! q38
            t4=   ra*qio-r9*qth+r8*qze+r7*qka-r6*qet ! q37
            t5=   ra*qka-r9*qze+r8*qio-r7*qet+r6*qth ! q36
            ra=r0-r1*rzc-r2*rec-r3*rtc-r4*ric-r5*rkc ! r31
            r9=r0-r1*rec-r2*ric-r3*rkc-r4*rtc-r5*rzc ! r32
            r8=r0-r1*rtc-r2*rkc-r3*rec-r4*rzc-r5*ric ! r33
            r7=r0-r1*ric-r2*rtc-r3*rzc-r4*rkc-r5*rec ! r34
            r6=r0-r1*rkc-r2*rzc-r3*ric-r4*rec-r5*rtc ! r35
            r(k0)=r0                   ! r40
            r(k1)=ra-t1                ! r41
            r(ka)=ra+t1                ! r4a
            r(k2)=r9-t2                ! r42
            r(k9)=r9+t2                ! r49
            r(k3)=r8-t3                ! r43
            r(k8)=r8+t3                ! r48
            r(k4)=r7-t4                ! r44
            r(k7)=r7+t4                ! r47
            r(k5)=r6-t5                ! r45
            r(k6)=r6+t5                ! r46
            t1=   qa*qze+q9*qet+q8*qth+q7*qio+q6*qka ! r3a
            t2=   qa*qet+q9*qio-q8*qka-q7*qth-q6*qze ! r39
            t3=   qa*qth-q9*qka-q8*qet+q7*qze+q6*qio ! r38
            t4=   qa*qio-q9*qth+q8*qze+q7*qka-q6*qet ! r37
            t5=   qa*qka-q9*qze+q8*qio-q7*qet+q6*qth ! r36
            qa=q0-q1*rzc-q2*rec-q3*rtc-q4*ric-q5*rkc ! q31
            q9=q0-q1*rec-q2*ric-q3*rkc-q4*rtc-q5*rzc ! q32
            q8=q0-q1*rtc-q2*rkc-q3*rec-q4*rzc-q5*ric ! q33
            q7=q0-q1*ric-q2*rtc-q3*rzc-q4*rkc-q5*rec ! q34
            q6=q0-q1*rkc-q2*rzc-q3*ric-q4*rec-q5*rtc ! q35
            q(k0)=q0                   ! q40
            q(k1)=qa+t1                ! q41
            q(ka)=qa-t1                ! q4a
            q(k2)=q9+t2                ! q42
            q(k9)=q9-t2                ! q49
            q(k3)=q8+t3                ! q43
            q(k8)=q8-t3                ! q48
            q(k4)=q7+t4                ! q44
            q(k7)=q7-t4                ! q47
            q(k5)=q6+t5                ! q45
            q(k6)=q6-t5                ! q46
         enddo
      enddo
      ma=mab
   enddo
endif
if(lnd > 0)then
!  RADIX 13
   nze=n/13
   rze=w(nze);          qze=w(nh+nze);       rzc=u1-rze ! <- zeta
   ret=rze*rze-qze*qze; qet=2*rze*qze;       rec=u1-ret ! <- eta 
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=u1-rth ! <- theta
   rio=ret*ret-qet*qet; qio=2*ret*qet;       ric=u1-rio ! <- iota
   rka=ret*rth-qet*qth; qka=ret*qth+qet*rth; rkc=u1-rka ! <- kappa
   rla=rth*rth-qth*qth; qla=2*rth*qth;       rlc=u1-rla ! <- lambda
   do L=1,Lnd
      mb=mb/13; mad=ma*13
      do j=0,ma-1
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3; jmb4=jmb*4; jmb5=jmb*5;
jmb6=jmb*6
         rf1=w(jmb);          qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);         qf2=w(nh+jmb2)      ! f**2
         rf3=w(jmb3);         qf3=w(nh+jmb3)      ! f**3
         rf4=w(jmb4);         qf4=w(nh+jmb4)      ! f**4
         rf5=w(jmb5);         qf5=w(nh+jmb5)      ! f**5
         rf6=w(jmb6);         qf6=w(nh+jmb6)      ! f**6
         rf7=rf3*rf4-qf3*qf4; qf7=rf3*qf4+qf3*rf4 ! f**7
         rf8=rf4*rf4-qf4*qf4; qf8=2*rf4*qf4       ! f**8
         rf9=rf4*rf5-qf4*qf5; qf9=rf4*qf5+qf4*rf5 ! f**9
         rfa=rf5*rf5-qf5*qf5; qfa=2*rf5*qf5       ! f**10
         rfb=rf5*rf6-qf5*qf6; qfb=rf5*qf6+qf5*rf6 ! f**11
         rfc=rf6*rf6-qf6*qf6; qfc=2*rf6*qf6       ! f**12
         do i=0,nm,mad
            k0=i+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            k7=k6+ma; k8=k7+ma; k9=k8+ma; ka=k9+ma; kb=ka+ma; kc=kb+ma
            r0=r(k0)
            r1=r(k1);r2=r(k2);r3=r(k3);r4=r(k4);r5=r(k5);r6=r(k6)
            r7=r(k7);r8=r(k8);r9=r(k9);ra=r(ka);rb=r(kb);rc=r(kc)
            q0=q(k0)
            q1=q(k1);q2=q(k2);q3=q(k3);q4=q(k4);q5=q(k5);q6=q(k6)
            q7=q(k7);q8=q(k8);q9=q(k9);qa=q(ka);qb=q(kb);qc=q(kc)
            t1=r1*qf1+q1*rf1           ! q11
            r1=r1*rf1-q1*qf1           ! r11
            q1=rc*rfc-qc*qfc           ! r1c
            rc=rc*qfc+qc*rfc           ! q1c
            qc=r1-q1                   ! r2c
            r1=r1+q1                   ! r21
            q1=t1+rc                   ! q21
            rc=t1-rc                   ! q2c
            t2=r2*qf2+q2*rf2           ! q12
            r2=r2*rf2-q2*qf2           ! r12
            q2=rb*rfb-qb*qfb           ! r1b
            rb=rb*qfb+qb*rfb           ! q1b
            qb=r2-q2                   ! r2b
            r2=r2+q2                   ! r22
            q2=t2+rb                   ! q22
            rb=t2-rb                   ! q2b
            t3=r3*qf3+q3*rf3           ! q13
            r3=r3*rf3-q3*qf3           ! r13
            q3=ra*rfa-qa*qfa           ! r1a
            ra=ra*qfa+qa*rfa           ! q1a
            qa=r3-q3                   ! r2a
            r3=r3+q3                   ! r23
            q3=t3+ra                   ! q23
            ra=t3-ra                   ! q2a
            t4=r4*qf4+q4*rf4           ! q14
            r4=r4*rf4-q4*qf4           ! r14
            q4=r9*rf9-q9*qf9           ! r19
            r9=r9*qf9+q9*rf9           ! q19
            q9=r4-q4                   ! r29
            r4=r4+q4                   ! r24
            q4=t4+r9                   ! q24
            r9=t4-r9                   ! q29
            t5=r5*qf5+q5*rf5           ! q15
            r5=r5*rf5-q5*qf5           ! r15
            q5=r8*rf8-q8*qf8           ! r18
            r8=r8*qf8+q8*rf8           ! q18
            q8=r5-q5                   ! r28
            r5=r5+q5                   ! r25
            q5=t5+r8                   ! q25
            r8=t5-r8                   ! q28
            t6=r6*qf6+q6*rf6           ! q16
            r6=r6*rf6-q6*qf6           ! r16
            q6=r7*rf7-q7*qf7           ! r17
            r7=r7*qf7+q7*rf7           ! q17
            q7=r6-q6                   ! r27
            r6=r6+q6                   ! r26
            q6=t6+r7                   ! q26
            r7=t6-r7                   ! q27
            r0=r0+r1+r2+r3+r4+r5+r6    ! r40
            q0=q0+q1+q2+q3+q4+q5+q6    ! q40
            t1=   rc*qze+rb*qet+ra*qth+r9*qio+r8*qka+r7*qla ! q3c
            t2=   rc*qet+rb*qio+ra*qla-r9*qka-r8*qth-r7*qze ! q3b
            t3=   rc*qth+rb*qla-ra*qio-r9*qze+r8*qet+r7*qka ! q3a
            t4=   rc*qio-rb*qka-ra*qze+r9*qth-r8*qla-r7*qet ! q39
            t5=   rc*qka-rb*qth+ra*qet-r9*qla-r8*qze+r7*qio ! q38
            t6=   rc*qla-rb*qze+ra*qka-r9*qet+r8*qio-r7*qth ! q37
            rc=r0-r1*rzc-r2*rec-r3*rtc-r4*ric-r5*rkc-r6*rlc ! r31
            rb=r0-r1*rec-r2*ric-r3*rlc-r4*rkc-r5*rtc-r6*rzc ! r32
            ra=r0-r1*rtc-r2*rlc-r3*ric-r4*rzc-r5*rec-r6*rkc ! r33
            r9=r0-r1*ric-r2*rkc-r3*rzc-r4*rtc-r5*rlc-r6*rec ! r34
            r8=r0-r1*rkc-r2*rtc-r3*rec-r4*rlc-r5*rzc-r6*ric ! r35
            r7=r0-r1*rlc-r2*rzc-r3*rkc-r4*rec-r5*ric-r6*rtc ! r36
            r(k0)=r0                   ! r40
            r(k1)=rc-t1                ! r41
            r(kc)=rc+t1                ! r4c
            r(k2)=rb-t2                ! r42
            r(kb)=rb+t2                ! r4b
            r(k3)=ra-t3                ! r43
            r(ka)=ra+t3                ! r4a
            r(k4)=r9-t4                ! r44
            r(k9)=r9+t4                ! r49
            r(k5)=r8-t5                ! r45
            r(k8)=r8+t5                ! r48
            r(k6)=r7-t6                ! r46
            r(k7)=r7+t6                ! r47
            t1=   qc*qze+qb*qet+qa*qth+q9*qio+q8*qka+q7*qla ! r3c
            t2=   qc*qet+qb*qio+qa*qla-q9*qka-q8*qth-q7*qze ! r3b
            t3=   qc*qth+qb*qla-qa*qio-q9*qze+q8*qet+q7*qka ! r3a
            t4=   qc*qio-qb*qka-qa*qze+q9*qth-q8*qla-q7*qet ! r39
            t5=   qc*qka-qb*qth+qa*qet-q9*qla-q8*qze+q7*qio ! r38
            t6=   qc*qla-qb*qze+qa*qka-q9*qet+q8*qio-q7*qth ! r37
            qc=q0-q1*rzc-q2*rec-q3*rtc-q4*ric-q5*rkc-q6*rlc ! q31
            qb=q0-q1*rec-q2*ric-q3*rlc-q4*rkc-q5*rtc-q6*rzc ! q32
            qa=q0-q1*rtc-q2*rlc-q3*ric-q4*rzc-q5*rec-q6*rkc ! q33
            q9=q0-q1*ric-q2*rkc-q3*rzc-q4*rtc-q5*rlc-q6*rec ! q34
            q8=q0-q1*rkc-q2*rtc-q3*rec-q4*rlc-q5*rzc-q6*ric ! q35
            q7=q0-q1*rlc-q2*rzc-q3*rkc-q4*rec-q5*ric-q6*rtc ! q36
            q(k0)=q0                   ! q40
            q(k1)=qc+t1                ! q41
            q(kc)=qc-t1                ! q4c
            q(k2)=qb+t2                ! q42
            q(kb)=qb-t2                ! q4b
            q(k3)=qa+t3                ! q43
            q(ka)=qa-t3                ! q4a
            q(k4)=q9+t4                ! q44
            q(k9)=q9-t4                ! q49
            q(k5)=q8+t5                ! q45
            q(k8)=q8-t5                ! q48
            q(k6)=q7+t6                ! q46
            q(k7)=q7-t6                ! q47
         enddo
      enddo
      ma=mad
   enddo
endif
end subroutine sgfft

!=============================================================================
subroutine dgfft(n,jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd,r,q)!  [gfft]
!=============================================================================
use pietc, only: u1,u2,u3,s30,ms30,s60
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(IN   ):: jumble
real(dp),    dimension(0:n-1),intent(IN   ):: w
integer(spi),                 intent(IN   ):: ln2,ln3,ln4,ln5,ln7,lnb,lnd
real(dp),    dimension(0:n-1),intent(INOUT):: r,q
!-----------------------------------------------------------------------------
integer(spi) :: nm,nh,i,j,l,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc, &
                ma,ma2,ma3,ma4,ma5,ma7,mab,mad,mb,mb2, &
                jmb,jmb2,jmb3,jmb4,jmb5,jmb6, &
                nze
real(dp)     :: r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,ra,rb,rc, &
                q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,qa,qb,qc,  &
                rf1,rf2,rf3,rf4,rf5,rf6,rf7,rf8,rf9,rfa,rfb,rfc,&
                qf1,qf2,qf3,qf4,qf5,qf6,qf7,qf8,qf9,qfa,qfb,qfc, &
                rep,rze,rzc,ret,rec,rth,rtc,rio,ric,rka,rkc,rla,rlc, &
                qep,qze,qet,qth,qio,qka,qla, &
                t1,t2,t3,t4,t5,t6
!=============================================================================
nm=n-1
nh=n/2  
! PERMUTE THE DATA:
do i=1,nm
   j=jumble(i)
   if(j > i)then
      t1=r(i); r(i)=r(j); r(j)=t1
      t1=q(i); q(i)=q(j); q(j)=t1
   endif
enddo

!  TRANSFORM THE DATA:
ma=1; mb=n
!  RADIX 4
do l=1,ln4
   mb=mb/4; ma4=ma*4; mb2=mb*2
   do j=0,ma-1
      jmb=j*mb; jmb2=j*mb2
      rf1=w(jmb);          qf1=w(nh+jmb)        ! f**1
      rf2=w(jmb2);         qf2=w(nh+jmb2)       ! f**2
      rf3=rf1*rf2-qf1*qf2; qf3=rf1*qf2+qf1*rf2  ! f**3
      do i=0,nm,ma4
         k0=i+j;    k1=k0+ma;   k2=k1+ma;   k3=k2+ma
         r0=r(k0); r1=r(k1);  r2=r(k2);  r3=r(k3)
         q0=q(k0); q1=q(k1);  q2=q(k2);  q3=q(k3)
         t1=r3*rf3-q3*qf3 ! q13
         q3=r3*qf3+q3*rf3 ! r13
         r3=r2*rf1-q2*qf1 ! r12
         r2=r2*qf1+q2*rf1 ! q12
         q2=q3    -r2     ! r23
         r2=q3    +r2     ! q22
         q3=r3    +t1     ! r22
         t1=r3    -t1     ! q23
         r3=r1*rf2-q1*qf2 ! r11
         q1=r1*qf2+q1*rf2 ! q11
         r1=r0    -r3     ! r21
         r0=r0    +r3     ! r20
         r(k3)=r1    -q2     ! r43
         r(k1)=r1    +q2     ! r41
         q2=q0    +q1     ! q20
         q1=q0    -q1     ! q21
         q(k0)=q2    +r2     ! q40
         q(k2)=q2    -r2     ! q42
         r(k2)=r0    -q3     ! r42
         r(k0)=r0    +q3     ! r40
         q(k3)=q1    -t1     ! q43
         q(k1)=q1    +t1     ! q41
      enddo
   enddo
   ma=ma4
enddo
if(ln2==1)then
!  RADIX 2
   mb=mb/2; ma2=ma*2
   do j=0,ma-1
      jmb=j*mb
      do i=0,nm,ma2
         k0=j+i;        k1=k0+ma
         rf1=w(jmb);    qf1=w(nh+jmb) ! f**1
         r0=r(k0);     q0=q(k0)
         r1=r(k1);     q1=q(k1)
         t1=r1*qf1+q1*rf1 ! q11
         q1=r1*rf1-q1*qf1 ! r11
         r(k1)=r0    -q1     ! r41
         r(k0)=r0    +q1     ! r40
         q(k1)=q0    -t1     ! q41
         q(k0)=q0    +t1     ! q40
      enddo
   enddo
   ma=ma2
endif
!  RADIX 3
rep=ms30;  rec=u3*s30;  qep=s60 ! <- epsilon
do l=1,ln3
   mb=mb/3; ma3=ma*3
   do j=0,ma-1
      jmb=j*mb
      rf1=w(jmb);          qf1=w(nh+jmb) ! f**1
      rf2=rf1*rf1-qf1*qf1; qf2=2*rf1*qf1 ! f**2
      do i=0,nm,ma3
         k0=i+j;     k1=k0+ma;    k2=k1+ma
         r1=r(k1);  q1=q(k1)
         r2=r(k2);  q2=q(k2)
         t1= r2*qf2+q2 *rf2 ! r12
         q2= r2*rf2-q2 *qf2 ! q12
         r2= r1*qf1+q1 *rf1 ! q11
         r1= r1*rf1-q1 *qf1 ! r11
         q1= r2    +t1      ! q21
         r2=(r2    -t1)*qep ! r22
         t1= r1    +q2      ! r21
         r1=(r1    -q2)*qep ! q22
         r(k0)= r(k0)+t1    ! r40
         q(k0)= q(k0)+q1    ! q40
         t1= r(k0)-t1 *rec  ! r21
         q1= q(k0)-q1 *rec  ! q21
         q(k2)= q1    -r1   ! q42
         q(k1)= q1    +r1   ! q41
         r(k1)= t1    -r2   ! r41
         r(k2)= t1    +r2   ! r42
      enddo
   enddo
   ma=ma3
enddo

if(ln5 > 0)then
!  RADIX 5
   nze=n/5;  rze=w(nze); qze=w(nh+nze); rzc=u1-rze ! <- zeta
   ret=rze*rze-qze*qze;  qet=2*rze*qze; rec=u1-ret ! <- eta
   do l=1,ln5
      mb=mb/5; ma5=ma*5
      do j=0,ma-1
         jmb=j*mb;             jmb2=jmb*2
         rf1=w(jmb);           qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);          qf2=w(nh+jmb2)      ! f**2
         rf3=rf1*rf2-qf1*qf2;  qf3=rf1*qf2+qf1*rf2 ! f**3
         rf4=rf2*rf2-qf2*qf2;  qf4=u2*rf2*qf2      ! f**4
         do i=0,nm,ma5
            k0=i+j;   k1=k0+ma;  k2=k1+ma;  k3=k2+ma;  k4=k3+ma
            r1=r(k1); r2=r(k2); r3=r(k3); r4=r(k4)
            q1=q(k1); q2=q(k2); q3=q(k3); q4=q(k4)
            t1=r1*qf1+q1*rf1       ! q11
            r1=r1*rf1-q1*qf1       ! r11
            q1=r4*rf4-q4*qf4       ! r14
            r4=r4*qf4+q4*rf4       ! q14
            q4=r1    -q1           ! q24
            r1=r1    +q1           ! r21
            q1=t1    +r4           ! q21
            r4=t1    -r4           ! q24
            t1=r3*rf3-q3*qf3       ! r13
            r3=r3*qf3+q3*rf3       ! q13
            q3=r2*qf2+q2*rf2       ! q12
            r2=r2*rf2-q2*qf2       ! r12
            q2=q3    +r3           ! q22
            r3=q3    -r3           ! q23
            q3=r2    -t1           ! r23
            r2=r2    +t1           ! r22
            r(k0)=r(k0)+r1+r2      ! r40
            q(k0)=q(k0)+q1+q2      ! q40
            t1=      r4*qze+r3*qet ! q34
            r3=      r4*qet-r3*qze ! q33
            r4=r(k0)-r1*rec-r2*rzc ! r32
            r1=r(k0)-r1*rzc-r2*rec ! r31
            r(k2)=r4   -r3         ! r42
            r(k3)=r4   +r3         ! r43
            r(k4)=r1   +t1         ! r44
            r(k1)=r1   -t1         ! r41
            t1=q(k0)-q1*rzc-q2*rec ! q31
            q2=q(k0)-q1*rec-q2*rzc ! q32
            q1=      q4*qet-q3*qze ! r33
            q4=      q4*qze+q3*qet ! r34
            q(k3)=q2   -q1         ! q43
            q(k2)=q2   +q1         ! q42
            q(k1)=t1   +q4         ! q41
            q(k4)=t1   -q4         ! q44
         enddo
      enddo
      ma=ma5
   enddo
endif

if(ln7 > 0)then
!  RADIX 7
   nze=n/7
   rze=w(nze);          qze=w(nh+nze);       rzc=u1-rze ! <- zeta 
   ret=rze*rze-qze*qze; qet=u2*rze*qze;      rec=u1-ret ! <- eta 
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=u1-rth ! <- theta 
   do L=1,Ln7
      mb=mb/7; ma7=ma*7
      do j=0,ma-1
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3
         rf1=w(jmb);          qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);         qf2=w(nh+jmb2)      ! f**2
         rf3=w(jmb3);         qf3=w(nh+jmb3)      ! f**3
         rf4=rf2*rf2-qf2*qf2; qf4=u2*rf2*qf2      ! f**4
         rf5=rf2*rf3-qf2*qf3; qf5=rf2*qf3+qf2*rf3 ! f**5
         rf6=rf3*rf3-qf3*qf3; qf6=u2*rf3*qf3      ! f**6
         do i=0,nm,ma7
            k0=i+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            r0=r(k0)
            r1=r(k1);r2=r(k2);r3=r(k3);r4=r(k4);r5=r(k5);r6=r(k6)
            q0=q(k0)
            q1=q(k1);q2=q(k2);q3=q(k3);q4=q(k4);q5=q(k5);q6=q(k6)
            t1=r1*qf1+q1*rf1           ! q11
            r1=r1*rf1-q1*qf1           ! r11
            q1=r6*rf6-q6*qf6           ! r16
            r6=r6*qf6+q6*rf6           ! q16
            q6=r1-q1                   ! r26
            r1=r1+q1                   ! r21
            q1=t1+r6                   ! q21
            r6=t1-r6                   ! q26
            t2=r2*qf2+q2*rf2           ! q12
            r2=r2*rf2-q2*qf2           ! r12
            q2=r5*rf5-q5*qf5           ! r15
            r5=r5*qf5+q5*rf5           ! q15
            q5=r2-q2                   ! r25
            r2=r2+q2                   ! r22
            q2=t2+r5                   ! q22
            r5=t2-r5                   ! q25
            t3=r3*qf3+q3*rf3           ! q13
            r3=r3*rf3-q3*qf3           ! r13
            q3=r4*rf4-q4*qf4           ! r14
            r4=r4*qf4+q4*rf4           ! q14
            q4=r3-q3                   ! r24
            r3=r3+q3                   ! r23
            q3=t3+r4                   ! q23
            r4=t3-r4                   ! q24
            r0=r0+r1+r2+r3             ! r40
            q0=q0+q1+q2+q3             ! q40
            t1=   r6*qze+r5*qet+r4*qth ! q36
            t2=   r6*qet-r5*qth-r4*qze ! q35
            t3=   r6*qth-r5*qze+r4*qet ! q34
            r6=r0-r1*rzc-r2*rec-r3*rtc ! r31
            r5=r0-r1*rec-r2*rtc-r3*rzc ! r32
            r4=r0-r1*rtc-r2*rzc-r3*rec ! r33
            r(k0)=r0                   ! r40
            r(k1)=r6-t1                ! r41
            r(k6)=r6+t1                ! r46
            r(k2)=r5-t2                ! r42
            r(k5)=r5+t2                ! r45
            r(k3)=r4-t3                ! r43
            r(k4)=r4+t3                ! r44
            t1=   q6*qze+q5*qet+q4*qth ! r36
            t2=   q6*qet-q5*qth-q4*qze ! r35
            t3=   q6*qth-q5*qze+q4*qet ! r34
            q6=q0-q1*rzc-q2*rec-q3*rtc ! q31
            q5=q0-q1*rec-q2*rtc-q3*rzc ! q32
            q4=q0-q1*rtc-q2*rzc-q3*rec ! q33
            q(k0)=q0                   ! q40
            q(k1)=q6+t1                ! q41
            q(k6)=q6-t1                ! q46
            q(k2)=q5+t2                ! q42
            q(k5)=q5-t2                ! q45
            q(k3)=q4+t3                ! q43
            q(k4)=q4-t3                ! q44
         enddo
      enddo
      ma=ma7
   enddo
endif
if(lnb > 0)then
!  RADIX 11
   nze=n/11
   rze=w(nze);          qze=w(nh+nze);       rzc=u1-rze ! <- zeta
   ret=rze*rze-qze*qze; qet=u2*rze*qze;      rec=u1-ret ! <- eta 
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=u1-rth ! <- theta
   rio=ret*ret-qet*qet; qio=u2*ret*qet;      ric=u1-rio ! <- iota
   rka=ret*rth-qet*qth; qka=ret*qth+qet*rth; rkc=u1-rka ! <- kappa
   do L=1,Lnb
      mb=mb/11; mab=ma*11
      do j=0,ma-1
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3; jmb4=jmb*4; jmb5=jmb*5
         rf1=w(jmb);          qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);         qf2=w(nh+jmb2)      ! f**2
         rf3=w(jmb3);         qf3=w(nh+jmb3)      ! f**3
         rf4=w(jmb4);         qf4=w(nh+jmb4)      ! f**4
         rf5=w(jmb5);         qf5=w(nh+jmb5)      ! f**5
         rf6=rf3*rf3-qf3*qf3; qf6=u2*rf3*qf3      ! f**6
         rf7=rf3*rf4-qf3*qf4; qf7=rf3*qf4+qf3*rf4 ! f**7
         rf8=rf4*rf4-qf4*qf4; qf8=u2*rf4*qf4      ! f**8
         rf9=rf4*rf5-qf4*qf5; qf9=rf4*qf5+qf4*rf5 ! f**9
         rfa=rf5*rf5-qf5*qf5; qfa=u2*rf5*qf5      ! f**10
         do i=0,nm,mab
            k0=i+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            k7=k6+ma; k8=k7+ma; k9=k8+ma; ka=k9+ma
            r0=r(k0)
            r1=r(k1);r2=r(k2);r3=r(k3);r4=r(k4);r5=r(k5);r6=r(k6)
            r7=r(k7);r8=r(k8);r9=r(k9);ra=r(ka)
            q0=q(k0)
            q1=q(k1);q2=q(k2);q3=q(k3);q4=q(k4);q5=q(k5);q6=q(k6)
            q7=q(k7);q8=q(k8);q9=q(k9);qa=q(ka)
            t1=r1*qf1+q1*rf1           ! q11
            r1=r1*rf1-q1*qf1           ! r11
            q1=ra*rfa-qa*qfa           ! r1a
            ra=ra*qfa+qa*rfa           ! q1a
            qa=r1-q1                   ! r2a
            r1=r1+q1                   ! r21
            q1=t1+ra                   ! q21
            ra=t1-ra                   ! q2a
            t2=r2*qf2+q2*rf2           ! q12
            r2=r2*rf2-q2*qf2           ! r12
            q2=r9*rf9-q9*qf9           ! r19
            r9=r9*qf9+q9*rf9           ! q19
            q9=r2-q2                   ! r29
            r2=r2+q2                   ! r22
            q2=t2+r9                   ! q22
            r9=t2-r9                   ! q29
            t3=r3*qf3+q3*rf3           ! q13
            r3=r3*rf3-q3*qf3           ! r13
            q3=r8*rf8-q8*qf8           ! r18
            r8=r8*qf8+q8*rf8           ! q18
            q8=r3-q3                   ! r28
            r3=r3+q3                   ! r23
            q3=t3+r8                   ! q23
            r8=t3-r8                   ! q28
            t4=r4*qf4+q4*rf4           ! q14
            r4=r4*rf4-q4*qf4           ! r14
            q4=r7*rf7-q7*qf7           ! r17
            r7=r7*qf7+q7*rf7           ! q17
            q7=r4-q4                   ! r27
            r4=r4+q4                   ! r24
            q4=t4+r7                   ! q24
            r7=t4-r7                   ! q27
            t5=r5*qf5+q5*rf5           ! q15
            r5=r5*rf5-q5*qf5           ! r15
            q5=r6*rf6-q6*qf6           ! r16
            r6=r6*qf6+q6*rf6           ! q16
            q6=r5-q5                   ! r26
            r5=r5+q5                   ! r25
            q5=t5+r6                   ! q25
            r6=t5-r6                   ! q26
            r0=r0+r1+r2+r3+r4+r5       ! r40
            q0=q0+q1+q2+q3+q4+q5       ! q40
            t1=   ra*qze+r9*qet+r8*qth+r7*qio+r6*qka ! q3a
            t2=   ra*qet+r9*qio-r8*qka-r7*qth-r6*qze ! q39
            t3=   ra*qth-r9*qka-r8*qet+r7*qze+r6*qio ! q38
            t4=   ra*qio-r9*qth+r8*qze+r7*qka-r6*qet ! q37
            t5=   ra*qka-r9*qze+r8*qio-r7*qet+r6*qth ! q36
            ra=r0-r1*rzc-r2*rec-r3*rtc-r4*ric-r5*rkc ! r31
            r9=r0-r1*rec-r2*ric-r3*rkc-r4*rtc-r5*rzc ! r32
            r8=r0-r1*rtc-r2*rkc-r3*rec-r4*rzc-r5*ric ! r33
            r7=r0-r1*ric-r2*rtc-r3*rzc-r4*rkc-r5*rec ! r34
            r6=r0-r1*rkc-r2*rzc-r3*ric-r4*rec-r5*rtc ! r35
            r(k0)=r0                   ! r40
            r(k1)=ra-t1                ! r41
            r(ka)=ra+t1                ! r4a
            r(k2)=r9-t2                ! r42
            r(k9)=r9+t2                ! r49
            r(k3)=r8-t3                ! r43
            r(k8)=r8+t3                ! r48
            r(k4)=r7-t4                ! r44
            r(k7)=r7+t4                ! r47
            r(k5)=r6-t5                ! r45
            r(k6)=r6+t5                ! r46
            t1=   qa*qze+q9*qet+q8*qth+q7*qio+q6*qka ! r3a
            t2=   qa*qet+q9*qio-q8*qka-q7*qth-q6*qze ! r39
            t3=   qa*qth-q9*qka-q8*qet+q7*qze+q6*qio ! r38
            t4=   qa*qio-q9*qth+q8*qze+q7*qka-q6*qet ! r37
            t5=   qa*qka-q9*qze+q8*qio-q7*qet+q6*qth ! r36
            qa=q0-q1*rzc-q2*rec-q3*rtc-q4*ric-q5*rkc ! q31
            q9=q0-q1*rec-q2*ric-q3*rkc-q4*rtc-q5*rzc ! q32
            q8=q0-q1*rtc-q2*rkc-q3*rec-q4*rzc-q5*ric ! q33
            q7=q0-q1*ric-q2*rtc-q3*rzc-q4*rkc-q5*rec ! q34
            q6=q0-q1*rkc-q2*rzc-q3*ric-q4*rec-q5*rtc ! q35
            q(k0)=q0                   ! q40
            q(k1)=qa+t1                ! q41
            q(ka)=qa-t1                ! q4a
            q(k2)=q9+t2                ! q42
            q(k9)=q9-t2                ! q49
            q(k3)=q8+t3                ! q43
            q(k8)=q8-t3                ! q48
            q(k4)=q7+t4                ! q44
            q(k7)=q7-t4                ! q47
            q(k5)=q6+t5                ! q45
            q(k6)=q6-t5                ! q46
         enddo
      enddo
      ma=mab
   enddo
endif
if(lnd > 0)then
!  RADIX 13
   nze=n/13
   rze=w(nze);          qze=w(nh+nze);       rzc=u1-rze ! <- zeta
   ret=rze*rze-qze*qze; qet=u2*rze*qze;      rec=u1-ret ! <- eta 
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=u1-rth ! <- theta
   rio=ret*ret-qet*qet; qio=u2*ret*qet;      ric=u1-rio ! <- iota
   rka=ret*rth-qet*qth; qka=ret*qth+qet*rth; rkc=u1-rka ! <- kappa
   rla=rth*rth-qth*qth; qla=u2*rth*qth;      rlc=u1-rla ! <- lambda
   do L=1,Lnd
      mb=mb/13; mad=ma*13
      do j=0,ma-1
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3; jmb4=jmb*4; jmb5=jmb*5;
jmb6=jmb*6
         rf1=w(jmb);          qf1=w(nh+jmb)       ! f**1
         rf2=w(jmb2);         qf2=w(nh+jmb2)      ! f**2
         rf3=w(jmb3);         qf3=w(nh+jmb3)      ! f**3
         rf4=w(jmb4);         qf4=w(nh+jmb4)      ! f**4
         rf5=w(jmb5);         qf5=w(nh+jmb5)      ! f**5
         rf6=w(jmb6);         qf6=w(nh+jmb6)      ! f**6
         rf7=rf3*rf4-qf3*qf4; qf7=rf3*qf4+qf3*rf4 ! f**7
         rf8=rf4*rf4-qf4*qf4; qf8=u2*rf4*qf4      ! f**8
         rf9=rf4*rf5-qf4*qf5; qf9=rf4*qf5+qf4*rf5 ! f**9
         rfa=rf5*rf5-qf5*qf5; qfa=u2*rf5*qf5      ! f**10
         rfb=rf5*rf6-qf5*qf6; qfb=rf5*qf6+qf5*rf6 ! f**11
         rfc=rf6*rf6-qf6*qf6; qfc=u2*rf6*qf6      ! f**12
         do i=0,nm,mad
            k0=i+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            k7=k6+ma; k8=k7+ma; k9=k8+ma; ka=k9+ma; kb=ka+ma; kc=kb+ma
            r0=r(k0)
            r1=r(k1);r2=r(k2);r3=r(k3);r4=r(k4);r5=r(k5);r6=r(k6)
            r7=r(k7);r8=r(k8);r9=r(k9);ra=r(ka);rb=r(kb);rc=r(kc)
            q0=q(k0)
            q1=q(k1);q2=q(k2);q3=q(k3);q4=q(k4);q5=q(k5);q6=q(k6)
            q7=q(k7);q8=q(k8);q9=q(k9);qa=q(ka);qb=q(kb);qc=q(kc)
            t1=r1*qf1+q1*rf1           ! q11
            r1=r1*rf1-q1*qf1           ! r11
            q1=rc*rfc-qc*qfc           ! r1c
            rc=rc*qfc+qc*rfc           ! q1c
            qc=r1-q1                   ! r2c
            r1=r1+q1                   ! r21
            q1=t1+rc                   ! q21
            rc=t1-rc                   ! q2c
            t2=r2*qf2+q2*rf2           ! q12
            r2=r2*rf2-q2*qf2           ! r12
            q2=rb*rfb-qb*qfb           ! r1b
            rb=rb*qfb+qb*rfb           ! q1b
            qb=r2-q2                   ! r2b
            r2=r2+q2                   ! r22
            q2=t2+rb                   ! q22
            rb=t2-rb                   ! q2b
            t3=r3*qf3+q3*rf3           ! q13
            r3=r3*rf3-q3*qf3           ! r13
            q3=ra*rfa-qa*qfa           ! r1a
            ra=ra*qfa+qa*rfa           ! q1a
            qa=r3-q3                   ! r2a
            r3=r3+q3                   ! r23
            q3=t3+ra                   ! q23
            ra=t3-ra                   ! q2a
            t4=r4*qf4+q4*rf4           ! q14
            r4=r4*rf4-q4*qf4           ! r14
            q4=r9*rf9-q9*qf9           ! r19
            r9=r9*qf9+q9*rf9           ! q19
            q9=r4-q4                   ! r29
            r4=r4+q4                   ! r24
            q4=t4+r9                   ! q24
            r9=t4-r9                   ! q29
            t5=r5*qf5+q5*rf5           ! q15
            r5=r5*rf5-q5*qf5           ! r15
            q5=r8*rf8-q8*qf8           ! r18
            r8=r8*qf8+q8*rf8           ! q18
            q8=r5-q5                   ! r28
            r5=r5+q5                   ! r25
            q5=t5+r8                   ! q25
            r8=t5-r8                   ! q28
            t6=r6*qf6+q6*rf6           ! q16
            r6=r6*rf6-q6*qf6           ! r16
            q6=r7*rf7-q7*qf7           ! r17
            r7=r7*qf7+q7*rf7           ! q17
            q7=r6-q6                   ! r27
            r6=r6+q6                   ! r26
            q6=t6+r7                   ! q26
            r7=t6-r7                   ! q27
            r0=r0+r1+r2+r3+r4+r5+r6    ! r40
            q0=q0+q1+q2+q3+q4+q5+q6    ! q40
            t1=   rc*qze+rb*qet+ra*qth+r9*qio+r8*qka+r7*qla ! q3c
            t2=   rc*qet+rb*qio+ra*qla-r9*qka-r8*qth-r7*qze ! q3b
            t3=   rc*qth+rb*qla-ra*qio-r9*qze+r8*qet+r7*qka ! q3a
            t4=   rc*qio-rb*qka-ra*qze+r9*qth-r8*qla-r7*qet ! q39
            t5=   rc*qka-rb*qth+ra*qet-r9*qla-r8*qze+r7*qio ! q38
            t6=   rc*qla-rb*qze+ra*qka-r9*qet+r8*qio-r7*qth ! q37
            rc=r0-r1*rzc-r2*rec-r3*rtc-r4*ric-r5*rkc-r6*rlc ! r31
            rb=r0-r1*rec-r2*ric-r3*rlc-r4*rkc-r5*rtc-r6*rzc ! r32
            ra=r0-r1*rtc-r2*rlc-r3*ric-r4*rzc-r5*rec-r6*rkc ! r33
            r9=r0-r1*ric-r2*rkc-r3*rzc-r4*rtc-r5*rlc-r6*rec ! r34
            r8=r0-r1*rkc-r2*rtc-r3*rec-r4*rlc-r5*rzc-r6*ric ! r35
            r7=r0-r1*rlc-r2*rzc-r3*rkc-r4*rec-r5*ric-r6*rtc ! r36
            r(k0)=r0                   ! r40
            r(k1)=rc-t1                ! r41
            r(kc)=rc+t1                ! r4c
            r(k2)=rb-t2                ! r42
            r(kb)=rb+t2                ! r4b
            r(k3)=ra-t3                ! r43
            r(ka)=ra+t3                ! r4a
            r(k4)=r9-t4                ! r44
            r(k9)=r9+t4                ! r49
            r(k5)=r8-t5                ! r45
            r(k8)=r8+t5                ! r48
            r(k6)=r7-t6                ! r46
            r(k7)=r7+t6                ! r47
            t1=   qc*qze+qb*qet+qa*qth+q9*qio+q8*qka+q7*qla ! r3c
            t2=   qc*qet+qb*qio+qa*qla-q9*qka-q8*qth-q7*qze ! r3b
            t3=   qc*qth+qb*qla-qa*qio-q9*qze+q8*qet+q7*qka ! r3a
            t4=   qc*qio-qb*qka-qa*qze+q9*qth-q8*qla-q7*qet ! r39
            t5=   qc*qka-qb*qth+qa*qet-q9*qla-q8*qze+q7*qio ! r38
            t6=   qc*qla-qb*qze+qa*qka-q9*qet+q8*qio-q7*qth ! r37
            qc=q0-q1*rzc-q2*rec-q3*rtc-q4*ric-q5*rkc-q6*rlc ! q31
            qb=q0-q1*rec-q2*ric-q3*rlc-q4*rkc-q5*rtc-q6*rzc ! q32
            qa=q0-q1*rtc-q2*rlc-q3*ric-q4*rzc-q5*rec-q6*rkc ! q33
            q9=q0-q1*ric-q2*rkc-q3*rzc-q4*rtc-q5*rlc-q6*rec ! q34
            q8=q0-q1*rkc-q2*rtc-q3*rec-q4*rlc-q5*rzc-q6*ric ! q35
            q7=q0-q1*rlc-q2*rzc-q3*rkc-q4*rec-q5*ric-q6*rtc ! q36
            q(k0)=q0                   ! q40
            q(k1)=qc+t1                ! q41
            q(kc)=qc-t1                ! q4c
            q(k2)=qb+t2                ! q42
            q(kb)=qb-t2                ! q4b
            q(k3)=qa+t3                ! q43
            q(ka)=qa-t3                ! q4a
            q(k4)=q9+t4                ! q44
            q(k9)=q9-t4                ! q49
            q(k5)=q8+t5                ! q45
            q(k8)=q8-t5                ! q48
            q(k6)=q7+t6                ! q46
            q(k7)=q7-t6                ! q47
         enddo
      enddo
      ma=mad
   enddo
endif
end subroutine dgfft

!=============================================================================
subroutine scsft(n,rb,qb)!  [csft]
!=============================================================================
!           SUBROUTINE SCSFT
! Complex slow Fourier transform. For inverse, see SDSFT.
! --> n:     Period
! <-> rb,qb: Real and imaginary parts of data and transform.
!=============================================================================
use pietc_s, only: u0,u1,pi2
integer(spi),             intent(IN   ):: n
real(sp),dimension(0:n-1),intent(INOUT):: rb,qb
!-----------------------------------------------------------------------------
integer(spi)             :: nm,i,k,ik
real(sp)                 :: dang,ang
real(sp),dimension(0:n-1):: ra,qa,wc,ws
!=============================================================================
nm=n-1
dang=pi2/n
do k=0,nm; ang=k*dang; wc(k)=cos(ang); ws(k)=sin(ang); enddo
ra=u0; qa=u0
do k=0,nm
   do i=0,nm
      ik=mod(i*k,n)
      ra(k)=ra(k)+wc(ik)*rb(i)-ws(ik)*qb(i)
      qa(k)=qa(k)+wc(ik)*qb(i)+ws(ik)*rb(i)
   enddo
enddo
rb=ra; qb=qa
end subroutine scsft 

!=============================================================================
subroutine dcsft(n,rb,qb)!  [csft]
!=============================================================================
!           SUBROUTINE DCSFT: Double precision version of csft
!=============================================================================
use pietc, only: u0,pi2
integer(spi),             intent(IN   ):: n
real(dp),dimension(0:n-1),intent(INOUT):: rb,qb
!-----------------------------------------------------------------------------
integer(spi)             :: nm,i,k,ik
real(dp)                 :: dang,ang
real(dp),dimension(0:n-1):: ra,qa,wc,ws
!=============================================================================
nm=n-1
dang=pi2/n
do k=0,nm; ang=k*dang; wc(k)=cos(ang); ws(k)=sin(ang); enddo
ra=u0; qa=u0
do k=0,nm
   do i=0,nm
      ik=mod(i*k,n)
      ra(k)=ra(k)+wc(ik)*rb(i)-ws(ik)*qb(i)
      qa(k)=qa(k)+wc(ik)*qb(i)+ws(ik)*rb(i)
   enddo
enddo
rb=ra; qb=qa
end subroutine dcsft 

!=============================================================================
subroutine sdsft(n,rb,qb)!  [dsft]
!=============================================================================
!           SUBROUTINE SDSFT
! Complex slow inverse Fourier transform. For direct transform, see
! CSFT.
! --> n:     Period
! <-> rb,qb: Real and imaginary parts of transform and synthesis.
!=============================================================================
use pietc_s, only: u0,u1,pi2
integer(spi),             intent(IN   ):: n
real(sp),dimension(0:n-1),intent(INOUT):: rb,qb
!-----------------------------------------------------------------------------
integer(spi)             :: nm,i,k,ik
real(sp)                 :: dang,ang,rfac
real(sp),dimension(0:n-1):: ra,qa,wc,ws
!=============================================================================
nm=n-1
rfac=u1/n
dang=pi2/n
do k=0,nm; ang=k*dang; wc(k)=cos(ang); ws(k)=sin(ang); enddo
ra=u0; qa=u0
do k=0,nm
   do i=0,nm
      ik=mod(i*k,n)
      ra(k)=ra(k)+wc(ik)*rb(i)+ws(ik)*qb(i)
      qa(k)=qa(k)+wc(ik)*qb(i)-ws(ik)*rb(i)
   enddo
enddo
rb=ra*rfac; qb=qa*rfac
end subroutine sdsft 

!=============================================================================
subroutine ddsft(n,rb,qb)!  [dsft]
!=============================================================================
!           SUBROUTINE DDSFT: Double precision version of dsft
!=============================================================================
use pietc, only: u1,pi2
integer(spi),             intent(IN   ):: n
real(dp),dimension(0:n-1),intent(INOUT):: rb,qb
!-----------------------------------------------------------------------------
integer(spi)             :: nm,i,k,ik
real(dp)                 :: dang,ang,rfac
real(dp),dimension(0:n-1):: ra,qa,wc,ws
!=============================================================================
nm=n-1
rfac=u1/n
dang=pi2/n
do k=0,nm; ang=k*dang; wc(k)=cos(ang); ws(k)=sin(ang); enddo
ra=0; qa=0
do k=0,nm
   do i=0,nm
      ik=mod(i*k,n)
      ra(k)=ra(k)+wc(ik)*rb(i)+ws(ik)*qb(i)
      qa(k)=qa(k)+wc(ik)*qb(i)-ws(ik)*rb(i)
   enddo
enddo
rb=ra*rfac; qb=qa*rfac
end subroutine ddsft 

!=============================================================================
subroutine sdsfe(n,rb,qb,x,y,r,q,dr,dq)!  [dsfe]
!=============================================================================
! Single precision complex "slow" Fourier evaluation.
! Evaluate the (complex) Fourier transform at a single point, together
! with
! its complex derivative there. Assume the period of the series is 2*pi.
!
! --> n     : number of imaginary Fourier coefficients of a complex
! function
! --> rb,qb : real and imaginary Fourier coefficients
! --> x,y   : real and imaginary components of complex representation of
! point.
! <-- r,q   : real and imaginary components of the function at the
! point.
! <-- dr,dq : real and imaginary components of derivative of function
! there.
!=============================================================================
use pietc_s, only: u0,u1
integer(spi),             intent(IN ):: n
real(sp),dimension(0:n-1),intent(IN ):: rb,qb
real(sp),                 intent(IN ):: x,y
real(sp),                 intent(OUT):: r,q,dr,dq
!-----------------------------------------------------------------------------
real(sp)    :: rfac,angx,angy,cangx,sangx,eangy,rbk,qbk,rbl,qbl,rterm,qterm
integer(spi):: k,nm,nh,nhm,L
!============================================================================
nm=n-1;  nh=(n+1)/2;  nhm=nh-1
r=u0; q=u0; dr=u0; dq=u0
rfac=u1/n
do k=0,nhm
   L=mod(n-k,n); rbk=rb(k); qbk=qb(k); rbl=rb(L); qbl=qb(L)
   angx=k*x; cangx=cos(angx); sangx=sin(angx); angy=k*y
   if(.not.(rbk==u0 .and. qbk==u0))then
      eangy=exp(angy)
      rterm=(rbk*cangx+qbk*sangx)*eangy;
qterm=(-rbk*sangx+qbk*cangx)*eangy
      r =r +rterm;    q =q +qterm;  if(k==0)cycle
      dr=dr+qterm*k;  dq=dq-rterm*k
   endif
   if(rbl==u0 .and. qbl==u0)cycle
   eangy=exp(-angy)
   rterm=(rbl*cangx-qbl*sangx)*eangy; qterm=( rbl*sangx+qbl*cangx)*eangy
   r=r  +rterm;    q =q +qterm
   dr=dr-qterm*k;  dq=dq+rterm*k
enddo
r=r*rfac; q=q*rfac; dr=dr*rfac; dq=dq*rfac
end subroutine sdsfe

!=============================================================================
subroutine ddsfe(n,rb,qb,x,y,r,q,dr,dq)!  [dsfe]
!=============================================================================
! Double precision complex "slow" Fourier evaluation.
! Evaluate the (complex) Fourier transform at a single point, together
! with
! its complex derivative there. Assume the period of the series is 2*pi.
!
! --> n     : number of imaginary Fourier coefficients of a complex
! function
! --> rb,qb : real and imaginary Fourier coefficients
! --> x,y   : real and imaginary components of complex representation of
! point.
! <-- r,q   : real and imaginary components of the function at the
! point.
! <-- dr,dq : real and imaginary components of derivative of function
! there.
!=============================================================================
use pietc, only: u0,u1
integer(spi),             intent(IN ):: n
real(dp),dimension(0:n-1),intent(IN ):: rb,qb
real(dp),                 intent(IN ):: x,y
real(dp),                 intent(OUT):: r,q,dr,dq
!-----------------------------------------------------------------------------
real(dp)    :: rfac,angx,angy,cangx,sangx,eangy,rbk,qbk,rbl,qbl,rterm,qterm
integer(spi):: k,nm,nh,nhm,L
!============================================================================
nm=n-1;  nh=(n+1)/2;  nhm=nh-1
r=u0; q=u0; dr=u0; dq=u0
rfac=u1/n
do k=0,nhm
   L=mod(n-k,n); rbk=rb(k); qbk=qb(k); rbl=rb(L); qbl=qb(L)
   angx=k*x; cangx=cos(angx); sangx=sin(angx); angy=k*y
   if(.not.(rbk==u0 .and. qbk==u0))then
      eangy=exp(angy)
      rterm=(rbk*cangx+qbk*sangx)*eangy;
qterm=(-rbk*sangx+qbk*cangx)*eangy
      r =r +rterm;    q =q +qterm;  if(k==0)cycle
      dr=dr+qterm*k;  dq=dq-rterm*k
   endif
   if(rbl==u0 .and. qbl==u0)cycle
   eangy=exp(-angy)
   rterm=(rbl*cangx-qbl*sangx)*eangy; qterm=( rbl*sangx+qbl*cangx)*eangy
   r=r  +rterm;    q =q +qterm
   dr=dr-qterm*k;  dq=dq+rterm*k
enddo
r=r*rfac; q=q*rfac; dr=dr*rfac; dq=dq*rfac
end subroutine ddsfe

!============================================================================
subroutine shsfe(n,b,x,r,dr)!  [hsfe]
!============================================================================
! Slow Hermitian Fourier synthesis evaluation and derivative
!============================================================================
use pietc_s, only: u2
integer(spi),             intent(IN ):: n
real(sp),dimension(0:n-1),intent(IN ):: b
real(sp),                 intent(IN ):: x
real(sp),                 intent(OUT):: r,dr
!----------------------------------------------------------------------------
integer(spi):: nm,nh,nhm,k,l
real(sp)    :: ang,cang,sang,bk,bl,rfac
!============================================================================
nm=n-1;  nh=(n+1)/2;  nhm=nh-1
rfac=u2/n
r=(b(0)+b(nh)*cos(x*nh))/u2
dr=0
do k=1,nhm
   l=n-k
   bk=b(k); bl=b(l)
   ang=x*k; cang=cos(ang); sang=sin(ang)
   r=r+bk*cang-bl*sang
   dr=dr-(bk*sang+bl*cang)*k
enddo
r=r*rfac; dr=dr*rfac
end subroutine shsfe
!
!============================================================================
subroutine dhsfe(n,b,x,r,dr)!  [hsfe]
!============================================================================
use pietc, only: u2
integer(spi),             intent(IN ):: n
real(dp),dimension(0:n-1),intent(IN ):: b
real(dp),                 intent(IN ):: x
real(dp),                 intent(OUT):: r,dr
!----------------------------------------------------------------------------
integer(spi):: nm,nh,nhm,k,l
real(dp)    :: ang,cang,sang,bk,bl,rfac
!============================================================================
nm=n-1;  nh=(n+1)/2;  nhm=nh-1
rfac=u2/n
r=(b(0)+b(nh)*cos(x*nh))/u2
dr=0
do k=1,nhm
   l=n-k
   bk=b(k); bl=b(l)
   ang=x*k; cang=cos(ang); sang=sin(ang)
   r=r+bk*cang-bl*sang
   dr=dr-(bk*sang+bl*cang)*k
enddo
r=r*rfac; dr=dr*rfac
end subroutine dhsfe

!=============================================================================
subroutine srfft(n,r)!  [rfft] 
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!   (Revised for f90, 1999)
!                   SUBROUTINES SRFFT
!   Fourier analyze (RFFT) a line of real data (for inverse, see HFFT)
!   The real coefficients of the transforms are stored in direct order
!   of wavenumber in elements 0 through n/2. The imaginary coefficients
!   of the transform are stored in reverse order of wavenumber in
!   the remaining elements. In this way, the structure of the transform 
!   of real data is made to conform more immediately to that of the 
!   complex data transform (provided in cfft).
!
! --> N       Number of data of series (product of 2s,3s,5s,7s,11s,13s)
! <-> R       Data and transform
!=============================================================================
integer(spi),             intent(IN   ):: n
real(sp),dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nhp
real(sp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nhp=n/2+1
t=0
call cfft(n,r,t)
r(nhp:nm)=t(nhp:nm)
end subroutine srfft 

!=============================================================================
subroutine drfft(n,r)!  [rfft] 
!=============================================================================
integer(spi),             intent(IN   ):: n
real(dp),dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nhp
real(dp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nhp=n/2+1
t=0
call cfft(n,r,t)
r(nhp:nm)=t(nhp:nm)
end subroutine drfft 

!=============================================================================
subroutine srfftp(n,js,ws,r)!  [rfft] 
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(sp),    dimension(0:n-1),intent(INOUT):: ws
real(sp),    dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nhp
real(sp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nhp=n/2+1
t=0
call cfft(n,js,ws,r,t)
r(nhp:nm)=t(nhp:nm)
end subroutine srfftp
!=============================================================================
subroutine drfftp(n,js,ws,r)!  [rfft] 
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(dp),    dimension(0:n-1),intent(INOUT):: ws
real(dp),    dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nhp
real(dp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nhp=n/2+1
t=0
call cfft(n,js,ws,r,t)
r(nhp:nm)=t(nhp:nm)
end subroutine drfftp 

!=============================================================================
subroutine shfft(n,r)!  [hfft]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!   (Revised for f90, 1999)
!                   SUBROUTINES SHFFT
!   Fourier synthesize (HFFT) a line of real data
!
! --> N       Number of data of series (product of 2s,3s,5s,7s,11s,13s)
! <-> R       Data and transform
!=============================================================================
integer(spi),             intent(IN   ):: n
real(sp),dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nh,nmh,nhp
real(sp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nh=n/2; nmh=nm/2; nhp=nh+1
t(0)=0; t(nh)=0
t(nhp:nm)=r(nhp:nm); t(1:nmh)=-t(nm:nhp:-1); r(nhp:nm)=r(nmh:1:-1)
call dfft(n,r,t)
end subroutine shfft 
!=============================================================================
subroutine dhfft(n,r)!  [hfft] 
!=============================================================================
integer(spi),             intent(IN   ):: n
real(dp),dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nh,nmh,nhp
real(dp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nh=n/2; nmh=nm/2; nhp=nh+1
t(0)=0; t(nh)=0
t(nhp:nm)=r(nhp:nm); t(1:nmh)=-t(nm:nhp:-1); r(nhp:nm)=r(nmh:1:-1)
call dfft(n,r,t)
end subroutine dhfft 

!=============================================================================
subroutine shfftp(n,js,ws,r)!  [hfft] 
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(sp),    dimension(0:n-1),intent(INOUT):: ws
real(sp),    dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nh,nmh,nhp
real(sp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nh=n/2; nmh=nm/2; nhp=nh+1
t(0)=0; t(nh)=0
t(nhp:nm)=r(nhp:nm); t(1:nmh)=-t(nm:nhp:-1); r(nhp:nm)=r(nmh:1:-1)
call dfft(n,js,ws,r,t)
end subroutine shfftp
!=============================================================================
subroutine dhfftp(n,js,ws,r)!  [hfft] 
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(dp),    dimension(0:n-1),intent(INOUT):: ws
real(dp),    dimension(0:n-1),intent(INOUT):: r
!-----------------------------------------------------------------------------
integer(spi)             :: nm,nh,nmh,nhp
real(dp),dimension(0:n-1):: t
!=============================================================================
nm=n-1; nh=n/2; nmh=nm/2; nhp=nh+1
t(0)=0; t(nh)=0
t(nhp:nm)=r(nhp:nm); t(1:nmh)=-t(nm:nhp:-1); r(nhp:nm)=r(nmh:1:-1)
call dfft(n,js,ws,r,t)
end subroutine dhfftp

!=============================================================================
subroutine sfftcnv(n,a,b,d)!  [fftcnv]
!=============================================================================
integer(spi),                  intent(IN ):: n
real(sp),     dimension(0:n-1),intent(IN ):: a,b
real(sp),     dimension(0:n-1),intent(OUT):: d
real(sp),     dimension(0:n-1)            :: ws
integer(spi), dimension(0:n-1)            :: js
!=============================================================================
js(0)=0
call fftcnv(n,js,ws,a,b,d)
end subroutine sfftcnv
!=============================================================================
subroutine dfftcnv(n,a,b,d)!  [fftcnv]
!=============================================================================
integer(spi),                 intent(IN ):: n
real(dp),    dimension(0:n-1),intent(IN ):: a,b
real(dp),    dimension(0:n-1),intent(OUT):: d
real(dp),    dimension(0:n-1)            :: ws
integer(spi),dimension(0:n-1)            :: js
!=============================================================================
js(0)=0
call fftcnv(n,js,ws,a,b,d)
end subroutine dfftcnv
!=============================================================================
subroutine sfftcnvp(n,js,ws,a,b,d)!  [fftcnv]
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(sp),    dimension(0:n-1),intent(INOUT):: ws
real(sp),    dimension(0:n-1),intent(IN   ):: a,b
real(sp),    dimension(0:n-1),intent(  OUT):: d
real(sp),    dimension(0:n-1)              :: at,bt
!=============================================================================
at=a
bt=b
call rfft(n,js,ws,at)
call rfft(n,js,ws,bt)
call fconv(n,at,bt,d)
call hfft(n,js,ws,d)
end subroutine sfftcnvp
!=============================================================================
subroutine dfftcnvp(n,js,ws,a,b,d)!  [fftcnv]
!=============================================================================
integer(spi),                 intent(IN   ):: n
integer(spi),dimension(0:n-1),intent(INOUT):: js
real(dp),    dimension(0:n-1),intent(INOUT):: ws
real(dp),    dimension(0:n-1),intent(IN   ):: a,b
real(dp),    dimension(0:n-1),intent(  OUT):: d
real(dp),    dimension(0:n-1)              :: at,bt
!=============================================================================
at=a
bt=b
call rfft(n,js,ws,at)
call rfft(n,js,ws,bt)
call fconv(n,at,bt,d)
call hfft(n,js,ws,d)
end subroutine dfftcnvp

!=============================================================================
subroutine sfconv(n,a,b,c)!  [fconv]
!=============================================================================
! Assuming Hermitian representations, A and B, of period n, form their
! product, C, that represents the transform of the convolution of the a
! and b
! whose transforms are aforementioned A and B.
!=============================================================================
integer(spi),             intent(IN ):: n
real(sp),dimension(0:n-1),intent(IN ):: a,b
real(sp),dimension(0:n-1),intent(OUT):: c
!-----------------------------------------------------------------------------
integer(spi)                              :: i,ic,nmnh,nh
!=============================================================================
nh=n/2
nmnh=n-nh      ! =nh for n even; =nh+1 for n odd
do i=0,nh,nmnh ! <- trick index stride and limit allows for no 2-wave at n odd
   c(i)=a(i)*b(i)
enddo
do i=1,nmnh-1; ic=n-i
   c(i )=a(i)*b(i )-a(ic)*b(ic)
   c(ic)=a(i)*b(ic)+a(ic)*b(i )
enddo
end subroutine sfconv
!=============================================================================
subroutine dfconv(n,a,b,c)!  [fconv]
!=============================================================================
integer(spi),             intent(IN ):: n
real(dp),dimension(0:n-1),intent(IN ):: a,b
real(dp),dimension(0:n-1),intent(OUT):: c
!-----------------------------------------------------------------------------
integer(spi)                              :: i,ic,nmnh,nh
!=============================================================================
nh=n/2
nmnh=n-nh      ! =nh for n even; =nh+1 for n odd
do i=0,nh,nmnh ! <- trick index stride and limit allows for no 2-wave at n odd
   c(i)=a(i)*b(i)
enddo
do i=1,nmnh-1; ic=n-i
   c(i )=a(i)*b(i )-a(ic)*b(ic)
   c(ic)=a(i)*b(ic)+a(ic)*b(i )
enddo
end subroutine dfconv

end module pfft


