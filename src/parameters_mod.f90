module parameters_mod
  use constants_mod
  implicit none
  
  ! Namelist parameters
  ! Domain
  integer :: ntang  ! ntang=1 for tang. vectors by finite diffs, ntang=2 for map factors by finite diffs too
  real    :: schmidt! schmidt coef 1 for uniform cubed sphere
  real    :: dx     !  grid-spacing in the x-direction
  real    :: dy     !  grid-spacing in the y-direction
  
  ! Index parameter
  integer :: ids      ! The starting index in the x-direction (Physical domain)
  integer :: ide      ! The ending index in the x-direction  (Physical domain)
  integer :: jds      ! The starting index in the y-direction  (Physical domain)
  integer :: jde      ! The ending index in the y-direction  (Physical domain)
  
  integer :: ifs      ! The starting index of patch(face)
  integer :: ife      ! The ending index of patch(face)
                      
  integer :: Nx       ! Element numbers in the x-direction
  integer :: Ny       ! Element numbers in the y-direction
  
  integer, parameter :: Nf = 6           ! Number of cube faces
  
  real, parameter :: x_min = -45.   !  start location of x-direction
  real, parameter :: x_max =  45.   !  end location of x-direction
  real, parameter :: y_min = -45.   !  start location of y-direction
  real, parameter :: y_max =  45.   !  end location of y-direction
  
  namelist /domain/ dx,dy,schmidt,ntang
  
  contains
  
  subroutine readNamelist
    
    open(1, file = 'namelist.input',status='old')
    read(1, nml  = domain       )
    close(1)
    
  end subroutine readNamelist
  
  subroutine initParameters
    
    ! Setting default values
    dx = 2.
    dy = 2.
    
    ! Read namelist
    call readNamelist
    
    ! Check if dx and dy are avaliable to achieve the integral element number
    if( (x_max - x_min)/dx - int((x_max - x_min)/dx) /= 0. )then
      stop '90 divide dx must be integer, choose another dx'
    end if
    
    if( (y_max - y_min)/dy - int((y_max - y_min)/dy) /= 0. )then
        stop '90 divide dy must be integer, choose another dy'
    end if
    
    ! Calculate element numbers on x/y direction
    Nx = int((x_max - x_min)/dx)
    Ny = int((y_max - y_min)/dy)
    nx = 2 * nx + 1
    ny = 2 * ny + 1
    
    ! Calculate starting and ending index for physical domain
    ids  = 1
    ide  = nx
    jds  = 1
    jde  = ny
    
    ! Setting the starting patch index and ending patch index
    ifs = 1
    ife = Nf
    
    !! Convert Degree to map coordinate
    !dx = dx * R2M
    !dy = dy * R2M
    
  end subroutine initParameters
  
end module parameters_mod
    