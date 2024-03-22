program mol2_amarrage
    use chargeur_covalence
    use lecture_mol2
    use affiche_topologie

    implicit none

    ! Variables
    character(len=128) :: filename, line

    integer :: i, j, k, ok, end, population, max_iter
    real :: merge_rate, mutate_rate
    type(Covalence) :: cov
    type(CovTable) :: table
    type(Molecule) :: ligand, site
    type(AtomXYZ) :: atom_xyz
    type(Topology) :: topo
    type(AtomXYZ), allocatable, dimension(:, :) :: solutions

    ! Mandatory arguments
    if (iargc() < 3) then
        call getarg(0, filename)
        print *, 'Error: Unvalid arguments'
        print '(2x, a, 1x, a)', trim(filename), &
            '<Cov_radii> <ligand_mol2> <target_mol2> ' // &
            '<population:100> <merge_rate:0.95> <mutate_rate:0.05>' // &
            '<max_iter:500>'
        stop 10
    end if

    ! Optional arguments
    population = 100
    merge_rate = 0.95
    mutate_rate = 0.05
    max_iter = 500

    if (iargc() >= 4) then
        call getarg(4, line)
        read(line, *) population
    endif

    if (iargc() >= 5) then
        call getarg(5, line)
        read(line, *), merge_rate
    endif

    if (iargc() >= 6) then
        call getarg(6, line)
        read(line, *), mutate_rate
    endif

    if (iargc() >= 7) then
        call getarg(7, line)
        read(line, *), max_iter
    endif

    ! Main routine

    ! Loading covalence table
    call getarg(1, filename)
    filename = trim(filename)
    table = charge_covalence(filename)

    ! Loading ligand file
    call getarg(2, filename)
    filename = trim(filename)
    ligand = lecture_fichier_mol2(filename)

    ! Loading target file
    call getarg(3, filename)
    filename = trim(filename)
    site = lecture_fichier_mol2(filename)

    ! Genetic algorithm

    ! Arg summary
    print '(a, a, a, i4, a)', "[mol2_amarrage] Ligand molecule : ", &
        trim(ligand%mol_name), " (", ligand%num_atoms, " atoms)"
    print '(a, a, a, i4, a)', "[mol2_amarrage] Site molecule   : ", &
        trim(site%mol_name), " (", site%num_atoms, " atoms)"
    print '(a, i4)', "[mol2_amarrage] population  = ", population
    print '(a, f3.2)', "[mol2_amarrage] merge_rate  = ", merge_rate
    print '(a, f3.2)', "[mol2_amarrage] mutate_rate = ", mutate_rate
    print '(a, i4)', "[mol2_amarrage] max_iter    = ", max_iter

    ! Initialization
    allocate(solutions(population, ligand%num_atoms))

    ! Find ligand topology
    topo = calcul_topologie(table, ligand, 0.10, 0.35, 0.05)
    ! do k = 1, topo%num_bonds
    !     print ('(i6, 1x, i4, 1x, i4, 1x, i1)'), k, &
    !         topo%bonds(1, k), topo%bonds(2, k), topo%bonds(3, k)
    ! end do
    ! TODO: What must we do with it ???

    ! Generate initial solutions

    print '(a)', "[mol2_amarrage] Generating initial solutions..."
    call RANDOM_SEED()  ! Init random
    do i = 1, population
        do j = 1, ligand%num_atoms
            call generate_atom(table, site, solutions(i, j))
        end do
    end do

    ! Wrap-up
    deallocate(solutions)

contains

    subroutine sample(scale, random)
        real, intent(in) :: scale
        real, intent(out) :: random

        call random_number(random)
        random = (random * (2 * scale + 1)) - scale
    end subroutine sample

    subroutine generate_atom(table, site, atom)
        ! Function signature
        type(Molecule), intent(in) :: site
        type(CovTable), intent(in) :: table
        type(AtomXYZ), intent(inout) :: atom

        integer :: i
        real :: distance, delta_simple, delta_double, delta_triple, tol
        type(Covalence) :: radii_a, radii_b

10      do while (.true.)
            ! exit = break

            ! Sample atom position
            call sample(10.0, atom%x)
            call sample(10.0, atom%y)
            call sample(10.0, atom%z)

            ! Check generated atom conformity
            tol = 0.10 ! Tolerance for checking bond between generated atom and molecule
            radii_a = table%get_cov_radii( atom%atom )
            do i = 1, site%num_atoms
                ! Compute distance in picometers
                distance = atom%compute_dist( site%atoms(i) )

                ! Get covalence radii for site atom
                radii_b = table%get_cov_radii( site%atoms(i)%atom )

                delta_simple = abs(radii_a%simple + radii_b%simple - distance) &
                    / (radii_a%simple + radii_b%simple)
                delta_double = abs(radii_a%double + radii_b%double - distance) &
                    / (radii_a%double + radii_b%double)
                delta_triple = abs(radii_a%triple + radii_b%triple - distance) &
                    / (radii_a%triple + radii_b%triple)

                ! Check if generated atom collide with target module
                if ( &
                    delta_simple <= tol .or. &
                    delta_double <= tol .or. &
                    delta_triple <= tol &
                    ) then
                    go to 10
                end if
            end do

            exit
        end do

    end subroutine generate_atom

    ! Calculates the angle in radians between two vectors in 3D space
    real function compute_angle(atom1, atom2, atom3)
        type(AtomXYZ), intent(in) :: atom1, atom2, atom3

        real :: ux, uy, uz, vx, vy, vz, norm_u, norm_v, dot_uv

        ! Calculate vector u (atom2 - atom1)
        ux = atom2%x - atom1%x
        uy = atom2%y - atom1%y
        uz = atom2%z - atom1%z

        ! Calculate vector v (atom3 - atom1)
        vx = atom3%x - atom1%x
        vy = atom3%y - atom1%y
        vz = atom3%z - atom1%z

        ! Calculate norms of vectors u and v
        norm_u = sqrt(ux * ux + uy * uy + uz * uz)
        norm_v = sqrt(vx * vx + vy * vy + vz * vz)

        ! Calculate dot product of u and v
        dot_uv = ux * vx + uy * vy + uz * vz

        ! Calculate angle using arccosine formula
        compute_angle = acos(dot_uv / (norm_u * norm_v))
    end function compute_angle

    logical function check_hydrogen_bond(atom_solution, atom_site)
        type(AtomXYZ), intent(in) :: atom_solution, atom_site

        ! TODO
    end function check_hydrogen_bond

    subroutine eval_solution(ligand, ligand_num_atoms, site, count)
        type(Molecule) :: site
        type(AtomXYZ), dimension(:) :: ligand
        integer :: count, ligand_num_atoms
        integer :: i, j

        count = 0
        do i = 1, ligand_num_atoms
            do j = 1, site%num_atoms
                if (check_hydrogen_bond(ligand(i), site%atoms(j))) then
                    count = count + 1
                end if
            end do
        end do

    end subroutine eval_solution

end program mol2_amarrage
