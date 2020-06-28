    ! Conformal Cubed Sphere mesh generator
    program CCS
      use constants_mod
      use precis_m
      use jim_utils
      use nfft_m
      implicit none
      
      integer i,j,k
      integer kit
      
      complex(kind=rx) :: cirs2, z, w, cang, cdang1, cdang2
      complex(kind=rx) :: cir, ciq, ciap
      real(kind=rx), dimension(40) :: a1, a2, b1, b2
      integer :: n1
      real(kind=rx) :: a, b, cx, cy, ssr, sc, sci, zac, st, sti
      real(kind=rx), dimension(3,3) :: rotm
      
      real(kind=rx) :: orc, sw1, sw2, rs1, rs2, rs3, rs4, pio4, u, v, rs3p, rs4p
      
      integer, parameter :: nt   = 50
      integer, parameter :: ntay = 30
      real   , parameter :: r    = 0.95
      integer, parameter :: n    = 128
      integer, parameter :: nh   = n/2
      
       real(kind=rx), dimension(  n) :: work, jumble
       real(kind=rx), dimension(0:n)  :: ra, qa
      
      complex :: ci = cmplx(0.0, 1.0, kind=rx)

      n1 = ntay
      
      !  First guess for first few taylor coefficients (based on a=1/3)
      !  needed to begin the boot-strap procedure:
      a1(1:n1) = 0.0
      a1(1) =  1.05060_rx
      a1(2) = -0.172629_rx
      a1(3) =  0.125716_rx
      a1(4) =  0.00015657_rx
      a1(5) = -0.0017866_rx
      a1(6) = -0.0031168_rx
      
      a2(1:n1) = 0.0
      a2(1) =  0.688411_rx
      a2(2) = -0.180924_rx
      a2(3) = -0.186593_rx
      a2(4) =  0.024811_rx
      a2(5) = -0.044888_rx
      a2(6) = -0.049190_rx
      b     = -0.147179_rx
      
      a      = 1./3.
      
      st     = 4.0/3.0
      sti    = 3.0/4.0
      ciap   = ci*a + 1.0
      cir    = (-ci)**0.75_rx
      ciq    = ci**0.25_rx
      orc    = 0.5_rx            ! relaxation factor
      sw1    = sqrt(1.0+a*a)
      sw2    = min ( a*2, sqrt(2.0_rx)*(1.0-a) )  ! Really just min?????
      ssr    = (sw1/sw2)**2
      rs1    = 0.95_rx*sw1
      rs2    = 0.95_rx*sw2
      rs3    = rs1**4
      rs4    = rs2**st
      
      pio4   = atan(1.0_rx)
      cdang1 = ci*pio4/nh
      cdang2 = -3.0*cdang1
      cirs2  = ci*rs2
      
      jumble(1) = 0.0
      
      if ( 3*ntay > n ) then
        print*, "too small in inoct, choose a larger power of 2"
        stop
      end if
      
      do kit = 1, nt
        do i = 0, nh
           cang = cdang1*i
           z = rs1*exp(cang)
           call toct ( z, w )
           w = w**4
           u = real(w)
           v = aimag(w)
           ra(i)   = u
           qa(i)   = -v
           ra(n-i) = u
           qa(n-i) = v
        enddo
        qa(0)  = 0.0
        qa(nh) = 0.0
        call cfft ( ra, qa, n, 1.0_rx, work, jumble )
        rs3p = 1.0
        do i = 1,n1
           rs3p = rs3p*rs3
           a1(i) = a1(i) + orc*(ra(i)/rs3p-a1(i))
        enddo
        
        do i = 0,nh
           cang = cdang2*i
           z = ciap-cirs2*exp(cang)
           call toct(z,w)
           w = ci*(w-1.0)/(w+1.0)
           u = real(w)
           v = aimag(w)
           ra(i) = u
           qa(i) = v
           ra(n-i) = u
           qa(n-i) = -v
        enddo
        qa(0) = 0.0
        qa(nh) = 0.0
        call cfft ( ra, qa, n, 1.0_rx, work, jumble )
        b = b + orc*(ra(0)-b)
        rs4p = 1.0
        do i=1,n1
           rs4p = rs4p*rs4
           a2(i) = a2(i) +orc*(ra(i)/rs4p-a2(i))
        enddo
      enddo
      
      do i = 1,n1
        print*,i,a1(i),a2(i)
      enddo
    contains
    
      !--------------------------------------------------------------------------
      !   r.j.purser, national meteorological center, washington d.c.  1995
      !                   subroutine  toct
      !   transform from complex-z in the standard unit-octagon to complex-w in the
      !   unit-circle
      !----------------------------------------------------------------------------
      subroutine toct ( z, w )
         complex(kind=rx), intent(in)  :: z
         complex(kind=rx), intent(out) :: w
         complex(kind=rx) :: zt
         logical :: kx, ky, kxy, kx1, kxy1, ks
         real(kind=rx) :: x, y, t, dd1, dd2
         
         x = real(z)
         y = aimag(z)
         kx = x < 0.0
         if ( kx ) then
            x = -x
         end if
         ky = y < 0.0
         if ( ky ) then
            y = -y
         end if
         kxy  =  y > x
         if ( kxy ) then
            t = x
            x = y
            y = t
         endif
         kx1 = x > 1.0
         kxy1 = x+y > 1.0+a
         if ( kx1 ) then
            x = 2.0 - x
         else if ( kxy1 ) then
            t = x
            x = 1.0 + a - y
            y = 1.0 + a - t
         endif
         dd1 = x**2 + y**2
         dd2 = (1.0-x)**2 + (a-y)**2
         
         zt = cmplx(x,y,kind=rx)
         ks = dd1 < ssr*dd2
         if ( ks ) then
            zt = zt**4
            call tay(zt,a1,n1,w)
            if ( w /= (0.0,0.0) ) then
               w = ciq*(-w*ci)**0.25_rx
            end if
         else
      !     jlm: note that DEC alpha has trouble with (0.,0.)**4/3
            zt = a+ci*zt-ci
            if ( zt /= cmplx(0.0,0.0,kind=rx) ) then
               zt = zt**st   !  jlm fix for DEC computers
            end if
            call tay(zt,a2,n1,w)
            w = w+b
            w = (w+ci)/(ci-w)
         endif
      
         if ( kx1 .or. kxy1 ) then
            w = w/abs(w)**2
         end if
      
         x = real(w)
         y = aimag(w)
         if ( kxy ) then
            t = x
            x = y
            y = t
         endif
         if ( ky ) then
            y = -y
         end if
         if ( kx ) then
            x = -x
         end if
         w = cmplx(x,y,kind=rx)
      end subroutine toct
    
    end program CCS
