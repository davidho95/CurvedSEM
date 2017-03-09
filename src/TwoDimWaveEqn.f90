
!===============================================================================
!
! TwoDimWaveEqn.f90
!
! Solves the wave equation on a (2 + 1) dimensional pseudo-Riemannian manifold
!
! David Ho (2017)
! 
! Adapted from Specfem1d v1.0, Komatisch & Tromp (2007)
!
!===============================================================================
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!===============================================================================

!===============================================================================
! Main program
!===============================================================================
program TwoDimWaveEqn

  use MeshModule
  use WaveModule

  implicit none

  type(mesh) :: my_mesh
  real(dp), allocatable :: coordinates(:,:)
  real(dp) delta_t
  integer, parameter :: NUM_TIME_STEPS = 5000
  integer i_step, ispecx, ispecy, ispec, inode, jnode, iglob
  integer :: output_step = 50

  character(len=43) :: output_path = "/home/davidho/WaveEqCurvedSpacetime/output/"
  character(len=50) snapshot_file

  call initialise_mesh(my_mesh, (/25, 25/), 4, 2d0, 1d0)

  call initial_conditions(my_mesh)

  allocate(coordinates(my_mesh%total_num_glob, 3))

  delta_t = estimate_timestep(my_mesh)

  do i_step = 1, NUM_TIME_STEPS
    call increment_system(my_mesh, delta_t)
    if (mod(i_step, output_step) == 0) then
      write(snapshot_file, "('snapshot',i5.5)") i_step
      open(unit=10, file=output_path//snapshot_file, action="write", status="unknown")

      coordinates = displ_on_torus(my_mesh, 0.5d0)

      do ispecx = 1,my_mesh%num_spec_el(1)
        do inode = 1,my_mesh%num_gll
          do ispecy = 1,my_mesh%num_spec_el(2)
            ispec = (ispecx-1)*my_mesh%num_spec_el(2) + ispecy
            do jnode = 1,my_mesh%num_gll
              iglob = my_mesh%i_bool(inode,jnode,ispec)
              write(10,*) coordinates(iglob,1),coordinates(iglob,2),&
                coordinates(iglob,3),my_mesh%vel(iglob),my_mesh%accel(iglob)
            end do
          end do
          write(10,*)
        end do
      end do
      close(10)
    endif
  enddo
end program TwoDimWaveEqn