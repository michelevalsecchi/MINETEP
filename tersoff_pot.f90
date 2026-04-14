! This module calculates the potential energy of silicon atoms in a unit cell
! using a Tersoff potential given the atomic coordinates and unit cell geometry.
!
! The functional form of the Tersoff potential and the
! empirical parameters used can be found at: J. Tersoff, Phys. Rev. B 38, 9902 (1988).


module potential

contains

!-------------------------------------------------------------------------------
!Defining the repulsive and attractive terms of the potential

  subroutine f_a(r, f) !Attractive term of potential
    implicit none
    real(kind=8), intent(in) :: r
    real(kind=8) :: B=471.18, l2=1.7322
    real(kind=8), intent(out) :: f

    f = -1.d0*B*exp(-l2*r)

  end subroutine f_a


  subroutine f_r(r, f)  !Repulsive term of potential
    implicit none
    real(kind=8), intent(in) :: r
    real(kind=8) :: A=1830.8, l1=2.4799
    real(kind=8), intent(out) :: f

    f = A*exp(-l1*r)

  end subroutine f_r


  !Defining the cut-off component of the potential & angular component of the coef-
  !-cients
  subroutine f_c(r, f) !Cut-off potential
    implicit none
    real(kind=8), intent(in) :: r
    real(kind=8) :: Rc=2.85, D=0.15
    real(kind=8), intent(out) :: f

    if (r < Rc-D) then
     f = 1
    else if (r > Rc+D) then
     f = 0
    else
     f = 0.5 - 0.5*(sin(0.5*acos(-1.d0)*(r-Rc)/D))
    end if

  end subroutine f_c

  subroutine g(cost, f) !Angular component of coefficients
    implicit none
    real(kind=8), intent(in) :: cost
    real(kind=8) :: c=1.0039*1e5, d=16.218, h=-0.59826
    real(kind=8), intent(out) :: f

    f = 1 + (c/d)**2 - c**2/(d**2 + (h - cost)**2)

  end subroutine g

!-------------------------------------------------------------------------------
!Tersoff potential function
subroutine u_pot(var, x0, Nc, V_tot)  !Full potential
  implicit none
  !Input variables
  integer, intent(in) :: Nc  !Number of atoms in the unit cell
  real(kind=8), intent(in) :: var(3*Nc+3)
  real(kind=8), intent(in), dimension(3) :: x0
  real(kind=8) :: beta=1.0999*1e-6, n=0.78734, l3=1.7322  !Tersoff potential parameters

  !Dummy variables
  real(kind=8) :: upos(Nc, 3)  !Motif atom positions (Nc, 3)
  real(kind=8), dimension(3) :: a, b, c  !Bravais lattice vectors
  real(kind=8) :: npos(3*Nc), dscal(27*Nc)  !Motif atom positions (3*Nc)
  real(kind=8), dimension(6) :: lpar  !Lattice parameters (a,b,c,θ,φ,ψ)
  integer, allocatable, dimension(:) :: neigh, near_neigh  !Array of neighbour indices in dscal and dvec
  integer, dimension(3) :: x, vec
  integer :: count, i1, j1, k1, l1, m1,n1,p1, s     !Loops
  integer, dimension(27, 3) :: p          !Permutations of [0, 1, -1]
  real(kind=8) :: pos(27*Nc, 3), dvec(27*Nc, 3) !Positions of all atoms within 27 unit cells of origin
  logical :: bool(27*Nc)    !Boolean
  real(kind=8) :: f_cut, g_res, f_rep, f_att
  real(kind=8) :: zeta, b_ij, cosjk        !Bond order coefficients

  !Output variables
  real(kind=8), intent(out) :: V_tot

  !-----------------------------------------------------------------------------

  !Extracting the atomic positions (npos) and lattice parameters (lpar)
  npos = [x0, var(:size(var)-6)]
  lpar = var(size(var)-6+1:)

  V_tot = 0.d0    !Initialise the potential energy
  upos = reshape(npos, shape(upos), order=[2, 1])  !Re-arrange input into N x 3 array

  !Write the Bravais lattice vectors in a cartesian coordinate system
  a = lpar(1)*[1.d0, 0.d0, 0.d0]
  b = lpar(2)*[cos(lpar(4)), sin(lpar(4)), 0.d0]
  c = lpar(3)*[cos(lpar(5)), sin(lpar(5))*cos(lpar(6)), sin(lpar(5))*sin(lpar(6))]

  !Generate the permutations with repetition of [0,1,-1] to obtain all
  !the unit cells that are neighbours of the first one
  x = [0, 1, -1]
  count = 0     !Initialise counter

  do i1=1,size(x)
    do j1=1,size(x)
      do k1=1,size(x)
        p(count + k1,:) = [x(i1), x(j1), x(k1)]
      end do
    count = count + 3
    end do
  end do

  !Find the positions of the atoms in the first 27 unit cells. The resulting
  !array is still an array of triples
  count = 0

  do i1 = 1,27
    vec = p(i1, :)
    do j1 = 1,Nc
      pos(j1 + count, :) = upos(j1, :) + vec(1)*a + vec(2)*b + vec(3)*c
    end do
  count = count + Nc
  end do

  !Compute the total potential energy of all the atoms in the first unit cell.
  !First iterate over the i atoms, then over the j neighbours of i and finally over
  !the k neighbours of i that don't include j.
  do i1 = 1,Nc
    do n1=1,3
      dvec(:,n1) = pos(:,n1) - upos(i1,n1)
    end do
    dscal = norm2(dvec, dim=2)

    !Creating array of booleans corresponding to where 0 < dscal < 3
    s = size(pack(dscal, dscal < 3 .and. dscal /= 0))

    if(.not. allocated(neigh)) allocate(neigh(s))
    if(.not. allocated(near_neigh)) allocate(near_neigh(s-1))
    bool = dscal < 3 .and. dscal /= 0

    !Generating array with indices of neighbours of i1 where 0 < dscal < 3
    count = 1
    do p1=1,size(bool)
      if (bool(p1) .eqv. .true.) then
        neigh(count) = p1
        count = count + 1
      end if
    end do

    !Performing summation over neighbours
    do j1 = 1,size(neigh)

      l1 = neigh(j1)
      zeta = 0              !Initialise zeta
      near_neigh = pack(neigh, neigh /= l1) !Delete jth neighbour of the ith atom from the list of atoms over which we iterate to compute zeta(ij)
      do k1 = 1,size(near_neigh)

        m1 = near_neigh(k1)
        cosjk = dot_product(dvec(l1, :), dvec(m1, :))/(dscal(l1)*dscal(m1))

        !Call all necessary functions for computation of zeta
        call f_c(dscal(m1), f_cut)
        call g(cosjk, g_res)
        !Compute zeta
        zeta = zeta + f_cut*g_res*exp((l3*(dscal(l1) - dscal(m1)))**3)
      end do

      b_ij = (1 + (beta*zeta)**n)**(-1/(2*n))

      !Compute V_tot accounting for double-counting
      !Call f_c, f_r, f_a
      call f_c(dscal(l1), f_cut)
      call f_r(dscal(l1), f_rep)
      call f_a(dscal(l1), f_att)

      V_tot = V_tot + f_cut*(f_rep + b_ij*f_att)
    end do
    if(allocated(neigh)) deallocate(neigh)
    if(allocated(near_neigh)) deallocate(near_neigh)
  end do

  V_tot = V_tot/2 !Accounting for double-counting
end subroutine u_pot

end module potential
