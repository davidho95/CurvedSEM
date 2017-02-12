module TwoDimWaveModule

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI = 3.141592653589793

  integer, dimension(2), parameter :: NUM_SPEC_EL = (/2, 2/)
  integer, parameter :: NUM_GLL = 3
  integer, dimension(2), parameter :: NUM_GLOBAL_POINTS = NUM_SPEC_EL * (NUM_GLL - 1) + (/1, 1/)
  integer :: total_num_glob
  contains

  subroutine Setup_mesh_2d()
    real(dp), dimension(0:NUM_SPEC_EL(1), 0:NUM_SPEC_EL(2), 2) :: anchor_points
    real(dp), dimension(NUM_GLL) :: gll_points, gll_weights
    real(dp), dimension(NUM_GLL, NUM_GLL) :: hprime

    integer, parameter :: TOTAL_NUM_SPEC = NUM_SPEC_EL(1) * NUM_SPEC_EL(2)
    integer, parameter :: TOTAL_NUM_NODES = NUM_GLL * NUM_GLL * TOTAL_NUM_SPEC
    real(dp), dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2) :: nodes
    integer, dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC) :: i_bool
    real(dp), dimension(TOTAL_NUM_NODES) :: temp_points_x, temp_points_y
    integer, dimension(TOTAL_NUM_NODES) :: locval
    logical, dimension(TOTAL_NUM_NODES) :: ifseg
    real(dp), dimension(:, :), allocatable :: global_points 

    integer i_gll, j_gll, i_spec, i_spec_x, i_spec_y, i_eoff, i_loc, i_glob

    ! get the GLL points and weights
    call zwgljd(gll_points,gll_weights,NUM_GLL,0.0_dp,0.0_dp)
    if(mod(NUM_GLL,2) /= 0) gll_points((NUM_GLL-1)/2+1) = 0.0_dp
    
    ! get the derivatives of the Lagrange polynomials at 
    ! the GLL points; recall that  hprime(i,j)=h'_{j}(gll_points_{i}) 
    do j_gll=1,NUM_GLL
       do i_gll=1,NUM_GLL
          hprime(i_gll,j_gll) = lagrange_deriv_GLL(j_gll-1,i_gll-1,gll_points,NUM_GLL)
       end do
    end do

    ! Setup anchor points
    do i_spec_x = 0, NUM_SPEC_EL(1)
      do i_spec_y = 0, NUM_SPEC_EL(2)
        anchor_points(i_spec_x, i_spec_y, 1) = 2 * PI &
         * dble(i_spec_x) / dble(NUM_SPEC_EL(1))
        anchor_points(i_spec_x, i_spec_y, 2) = 2 * PI &
         * dble(i_spec_y) / dble(NUM_SPEC_EL(2)) 
      enddo
    enddo

    ! Setup global points
    i_spec = 0
    do i_spec_x = 1, NUM_SPEC_EL(1)
      do i_spec_y = 1, NUM_SPEC_EL(2)
        i_spec = i_spec + 1
        do i_gll = 1, NUM_GLL
           do j_gll = 1, NUM_GLL
              nodes(i_gll,j_gll,i_spec, 1) =&
                anchor_points(i_spec_x-1, i_spec_y-1, 1) +&
                (gll_points(i_gll)+1.0_dp)*0.5_dp*&
                (anchor_points(i_spec_x, i_spec_y, 1) -&
                anchor_points(i_spec_x - 1, i_spec_y - 1, 1))
              nodes(i_gll,j_gll,i_spec, 2) =&
                anchor_points(i_spec_x-1, i_spec_y-1, 2) +&
                (gll_points(i_gll)+1.0_dp)*0.5_dp*&
                (anchor_points(i_spec_x, i_spec_y, 2) -&
                anchor_points(i_spec_x - 1, i_spec_y - 1, 2))
           enddo
        enddo
      enddo
    enddo

    do i_spec = 1, TOTAL_NUM_SPEC
       i_eoff = NUM_GLL*NUM_GLL*(i_spec-1)
       i_loc = 0
       do j_gll = 1,NUM_GLL
          do i_gll = 1,NUM_GLL
             i_loc = i_loc + 1
             temp_points_x(i_loc+i_eoff) = nodes(i_gll,j_gll,i_spec, 1)
             temp_points_y(i_loc+i_eoff) = nodes(i_gll, j_gll, i_spec, 2)
          end do
       end do
    end do

    locval = 0.
    ifseg = .FALSE.

    call get_global(NUM_GLL,TOTAL_NUM_SPEC,temp_points_x,temp_points_y,&
      i_bool,locval,ifseg,total_num_glob,TOTAL_NUM_NODES)

    allocate(global_points(total_num_glob, 2))

    do i_spec = 1, TOTAL_NUM_SPEC
      do j_gll = 1, NUM_GLL
        do i_gll = 1, NUM_GLL
          i_glob = i_bool(i_gll, j_gll, i_spec)
          global_points(i_glob, :) = nodes(i_gll, j_gll, i_spec, :)
        enddo
      enddo
    enddo

    do i_glob = 1, total_num_glob
      print *, global_points(i_glob, 1), global_points(i_glob, 2)
    enddo

  end subroutine Setup_mesh_2d

end module TwoDimWaveModule