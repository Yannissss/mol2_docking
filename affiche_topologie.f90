module affiche_topologie
    use chargeur_covalence
    use lecture_mol2
    use omp_lib

    implicit none

    integer, parameter :: BOND_SIMPLE = 1
    integer, parameter :: BOND_DOUBLE = 2
    integer, parameter :: BOND_TRIPLE = 3

    type, public :: Topology
        ! (atom_a_num, atom_b_num, bond_type)
        integer, allocatable, dimension(:, :) :: bonds
        integer :: num_bonds
    end type Topology

contains

    type(Topology) function calcul_topologie(table, mol, tol_min, tol_max, tol_step)
        ! Function signature
        type(CovTable), intent(in) :: table
        type(Molecule), intent(in) :: mol
        real, intent(in) :: tol_min, tol_max, tol_step

        ! Variables
        integer :: i, j, k
        type(AtomXYZ) :: atom_a, atom_b
        type(Covalence) :: radii_a, radii_b
        real :: distance, tol, delta_simple, delta_double, delta_triple
        logical :: valid_topology

        ! OpenMP vars
        double precision :: start_time, end_time, elapsed ! Time utils
        integer :: num_threads

        ! Auxiliary matrices
        ! Delta matrix
        real, allocatable, dimension(:, :, :) :: deltas
        ! Type of chemical bound between atoms i & j / 0 = NO BOUND
        integer, allocatable, dimension(:, :) :: atoms_bond


        ! Allocate auxiliary tables
        allocate(deltas(3, mol%num_atoms, mol%num_atoms))
        allocate(atoms_bond(mol%num_atoms, mol%num_atoms))

        ! Main loop
        print '(a, a)', "[affiche_topologie] Solving topology of ", &
            trim(mol%mol_name)
        print '(a, f3.2, a, f3.2, a, f3.2)', "[affiche_topologie] w/ tol_min = ", &
            tol_min, ", tol_max = ", tol_max, ", tol_step = ", tol_step

        !$omp parallel shared(num_threads)

        !$omp master
        num_threads = omp_get_num_threads()
        print '(a, i3, a)', "[affiche_topologie] and using ", &
            num_threads, " thread(s)"
        !$omp end master

        !$omp end parallel

        ! Start stopwatch
        start_time = omp_get_wtime()

        ! First we compute the deltas matrix that is independant of the tolerance used

        !$omp parallel &
        !$omp shared(mol, table, deltas) &
        !$omp private(i, j, atom_a, atom_b, radii_a, radii_b, distance, delta_simple, delta_double, delta_triple)
        !$omp do schedule(guided) collapse(1)
        do i=1, mol%num_atoms
            do j=i + 1, mol%num_atoms
                ! We check if there is a bond between atom A & B
                atom_a = mol%atoms(i) ! Read atom positions
                atom_b = mol%atoms(j)

                radii_a = table%get_cov_radii(atom_a%atom) ! Read elements
                radii_b = table%get_cov_radii(atom_b%atom) ! radii

                ! Compute distance in picometers
                distance = atom_a%compute_dist(atom_b)

                ! Small maths to compute relative distance deltas between two atoms
                ! d: distance, t: tol, r: radius_a + radius_b
                ! d in [r(1 - t); r(1 + t)]
                ! <=> r * (1 - t) <= d <= r * (1 + t)
                ! <=>  -t * r <= d - r <= t * r
                ! <=>  -t <= (d - r) / r <= t
                ! <=> |d - r| / r <= t

                ! Compute delta for differents bond types
                delta_simple = abs(radii_a%simple + radii_b%simple - distance) &
                    / (radii_a%simple + radii_b%simple)
                delta_double = abs(radii_a%double + radii_b%double - distance) &
                    / (radii_a%double + radii_b%double)
                delta_triple = abs(radii_a%triple + radii_b%triple - distance) &
                    / (radii_a%triple + radii_b%triple)

                ! Store deltas in matrix
                deltas(1, i, j) = delta_simple
                deltas(2, i, j) = delta_double
                deltas(3, i, j) = delta_triple
            end do
        end do
        !$omp end do
        !$omp end parallel

        ! Next we try to construct a compliant topology by increasing our
        ! delta tolerance step by stem until atoms are bound
        atoms_bond = 0
        valid_topology = .false.
        tol = tol_min
        k = 0 ! Number of bonds found
        do while ( (.not. valid_topology) .and. (tol <= tol_max) )
            ! Topology enhancement loop

            !$omp parallel &
            !$omp shared(mol, deltas, atoms_bond) &
            !$omp private(i, j) &
            !$omp reduction(+ : k)
            !$omp do schedule(guided) collapse(1)
            do i=1, mol%num_atoms
                do j=i + 1, mol%num_atoms
                    ! Check if there is already a bond between these 2 atoms
                    if (atoms_bond(i, j) == 0) then
                        ! We always check for triple, double and then simple
                        ! because delta_triple < delta_double < delta_simple
                        if (deltas(BOND_TRIPLE, i, j) <= tol) then
                            atoms_bond(i, j) = BOND_TRIPLE
                            k = k + 1
                        else if (deltas(BOND_DOUBLE, i, j) <= tol) then
                            atoms_bond(i, j) = BOND_DOUBLE
                            k = k + 1
                        else if (deltas(BOND_SIMPLE, i, j) <= tol) then
                            atoms_bond(i, j) = BOND_SIMPLE
                            k = k + 1
                        end if
                    end if
                end do
            end do
            !$omp end do
            !$omp end parallel

            ! Check topology compliance
            ! We just check is there is any atom no bound to anything
            valid_topology = .true.

            !$omp parallel &
            !$omp shared(atoms_bond) &
            !$omp private(i) &
            !$omp reduction(.and. : valid_topology)
            !$omp do schedule(static)
            do i = 1, mol%num_atoms
                ! Check if atom i has no bond
                if (all(atoms_bond(i, :) == 0) &
                    .and. all(atoms_bond(:, i) == 0)) then

                    ! and continue to search for a solution
                    valid_topology = .false.
                end if
            end do
            !$omp end do
            !$omp end parallel

            ! We increase the tolerance if no valid topology was found
            if (.not. valid_topology) then
                tol = tol + tol_step
            end if
        end do

        ! Stop stopwatch
        end_time = omp_get_wtime()

        ! Check if we found a valid topology for molecule
        if (valid_topology) then
            ! Build final topology
            calcul_topologie%num_bonds = k
            allocate(calcul_topologie%bonds(3, k))

            k = 1
            do i = 1, mol%num_atoms
                do j = (i + 1), mol%num_atoms
                    if (atoms_bond(i, j) /= 0) then
                        calcul_topologie%bonds(1, k) = i ! atom_a
                        calcul_topologie%bonds(2, k) = j ! atom_b
                        calcul_topologie%bonds(3, k) = atoms_bond(i, j) ! bond type
                        k = k + 1
                    end if
                end do
            end do

            ! Post execution stats
            print '(a, i2, a)', "[affiche_topologie] Done with max tolerance of ", &
                int(100.0 * tol), "%"
            print '(a, i4, a)', "[affiche_topologie] Found ", &
                calcul_topologie%num_bonds, ' bond(s)'

        else ! Otherwise we could not find a valid topology with our tolerance
            print '(a, a)', &
                "[affiche_topologie] Could not found a valid topology for molecule ", &
                trim(mol%mol_name)
            print '(a, i2, a)', "[affiche_topologie] with max tolerance of ", &
                int(100.0 * tol_max), "%"

            ! This means that no valid topology was found
            calcul_topologie%num_bonds = 0
        end if

        ! Display stopwatch
        elapsed = end_time - start_time

        write (*, '(a)', advance='no'), "[affiche_topologie] Solver execution done in "
        if (elapsed >= 1.0d-0) then
            print '(f5.1, a)', elapsed * 1.0d+0, " s"
        else if (elapsed >= 1.0d-3) then
            print '(f5.1, a)', elapsed * 1.0d+3, " ms"
        else if (elapsed >= 1.0d-6) then
            print '(f5.1, a)', elapsed * 1.0d+6, " µs"
        else if (elapsed >= 1.0d-9) then
            print '(f5.1, a)', elapsed * 1.0d+9, " ns"
        end if

        ! Clean-up
        deallocate(deltas)
        deallocate(atoms_bond)

    end function calcul_topologie

end module affiche_topologie
