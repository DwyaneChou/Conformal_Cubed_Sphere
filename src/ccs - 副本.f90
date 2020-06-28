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
      
      integer, parameter :: nt   = 2
      integer, parameter :: ntay = 30
      real   , parameter :: r    = 0.95
      integer, parameter :: n    = 128
      integer, parameter :: nh   = n/2
      
       real(kind=rx), dimension(  n) :: work,jumble
       real(kind=rx), dimension(0:n) :: ra, qa
       
       real(kind=rx), dimension(ntay) :: Ak
       real(kind=rx), dimension(ntay) :: Ak0
       real CA
      
      complex :: ci = cmplx(0.0, 1.0, kind=rx)
      
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
      rs1    = r*sw1
      rs2    = r*sw2
      rs3    = rs1**4
      rs4    = rs2**st
      
      if ( 3*ntay > n ) then
        print*, "too small in inoct, choose a larger power of 2"
        stop
      end if
      
      work   = 0
      jumble = 0.0
      
      cdang1 = ci*pi/2/nh
      
      Ak(1) = 1.47713062600964
      CA   = sqrt(3.) * gamma(1./3.) * Ak(1)**(1./3.) / (256.**(1./3.)*pi)
      
      do k = 2,ntay
        Ak(k) = -CA * k**(-7./3.)
      enddo
      
      Ak0 = Ak
      
      do kit = 1, nt
        do i = 0, nh
          cang = - 0.25 * pi * ci + cdang1*i
          !z = rs1*exp(cang)
          z = r*exp(cang)
          z = 1. - z
          call tay(z**4,Ak,ntay,w)
          w = w**(1./3.)
          w = ( 1. - w ) / ( 1. + w/2. )
          w = w**3.
                  
          u       = real(w)
          v       = aimag(w)
          ra(i)   = u
          qa(i)   = -v
          ra(n-i) = u
          qa(n-i) = v
        enddo
        qa(0)  = 0.0
        qa(nh) = 0.0
        
        call cfft ( ra, qa, n, 1.0_rx, work, jumble )

        rs3p = 1.0
        do k = 1,ntay
          !rs3p = rs3p*rs3
          rs3p = rs3p*r
          Ak(k) = Ak(k) + orc*(ra(k)/rs3p-Ak(k))
          !Ak(k) = Ak(k) + 0.5 * ( ra(k) / r**(4*k) - Ak(k) )
          !Ak(k) = ra(k) / r**(4*k)
        enddo
      enddo
      
      do k = 1,ntay
        print*,k,Ak(k),Ak0(k)
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
        logical :: kx, ky, kxy, kx1, kxy1
        real(kind=rx) :: x, y, t
        
        x = real(z)
        y = aimag(z)
        kx = x > 0.5_rx
        ky = y > 0.5_rx
        kxy = y > x
        
        if ( kx ) then
           x = 1.0 - x
        endif
        if ( ky ) then
           y = 1.0 - y
        endif
        if ( y > x ) then
           t = x
           x = y
           y = t
        endif
        
        zt = cmplx(x,y,kind=rx)
        zt = zt**4
        call tay(zt,Ak,ntay,w)
        w = w**(1./3.)
        w = ( 1. - w ) / ( 1. + w/2. )
        w = w**3.
        
        x = real(w)
        y = aimag(w)
        
        if ( kx ) then
           x = -x
        end if
        if ( ky ) then
           y = -y
        end if
        if ( kxy ) then
           t = x
           x = y
           y = t
        endif
        
        w = cmplx(x,y,kind=rx)
        
      end subroutine toct
    
    end program CCS
