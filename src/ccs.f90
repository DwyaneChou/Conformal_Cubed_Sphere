    ! Conformal Cubed Sphere mesh generator
    program CCS
      use constants_mod
      use precis_m
      use jim_utils
      use nfft_m
      implicit none
      
      integer i,j,k
      integer kit
      real    u,v
      complex w,z
      
      integer, parameter :: nt   = 20
      integer, parameter :: ntay = 30
      real   , parameter :: r    = 0.95
      integer, parameter :: n    = 360
      integer, parameter :: nh   = n/2
      
       real(kind=rx), dimension(1200) :: work
       real(kind=rx), dimension(0:n)  :: ra, qa
       real(kind=rx), dimension(n)    :: jumble
       
       real(kind=rx), dimension(ntay) :: A
       real CA
      
      complex :: ci = cmplx(0.0, 1.0, kind=rx)
      complex :: cdang1
      complex :: cang
      
      if ( 3*ntay > n ) then
        print*, "too small in inoct, choose a larger power of 2"
        stop
      end if
      
      cdang1 = ci*pi/2/nh
      
      a(1) = 1.47713062600964
      CA   = sqrt(3.) * gamma(1./3.) * a(1)**(1./3.) / (256.**(1./3.)*pi)
      
      do k = 2,ntay
        a(k) = -CA * k**(-7./3.)
      enddo
      
      do kit = 1, nt
        do i = 0, nh
          cang = cdang1*i - 0.25 * pi * ci
          !print*,i,dble(cang/ci*180./pi)
          z = r*exp(cang)
          z = 1. - z
          !call toct ( z, w )
          call tay(z**4,a,ntay,w)
          w       = w**(1./3.)
          w       = ( 1. - w ) / ( 1. + w/2. )
          w       = w**3.
                  
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

        do k = 1,ntay
          !a(k) = a(k) + 0.5 * ( ra(k) / r**(4*k) - a(k) )
          a(k) = ra(k) / r**(4*k)
        enddo
      enddo
      
      do k = 1,ntay
        print*,k,a(k)
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
        call tay(zt,a,ntay,w)
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
