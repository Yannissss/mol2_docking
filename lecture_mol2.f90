module lecture_mol2
    implicit none

    type, public :: AtomXYZ
        character(len=2) :: atom
        real :: x, y, z
    contains
        procedure compute_dist
    end type AtomXYZ

    type, public :: Molecule
        character(len=128) :: mol_name
        integer :: num_atoms
        type(AtomXYZ), dimension(:), allocatable :: atoms
    end type Molecule

contains

    ! Compute L2 norm between two atoms in picometers
    real function compute_dist(self, other)
        class(AtomXYZ), intent(in) :: self, other

        real :: dx, dy, dz

        ! Convert Angstrom to picometers (1 A° = 100 pm)
        dx = 100.0 * (self%x - other%x)
        dy = 100.0 * (self%y - other%y)
        dz = 100.0 * (self%z - other%z)

        ! Return distance
        compute_dist = sqrt(dx * dx + dy * dy + dz * dz)

    end function compute_dist

    type(Molecule) function lecture_fichier_mol2(filename)
        ! Signature
        character(len=*), intent(in) :: filename

        ! Variables
        integer :: i, j, num_args, fd, end, ok
        character(len=128) :: line

        character(len=128) :: mol_name
        integer :: num_atoms, atom_num
        character(len=12) :: atom_name, atom_type
        real :: x, y, z
        type(AtomXYZ) :: atom_xyz
        type(AtomXYZ), dimension(:), allocatable :: atoms

        ! Open file
        fd = 101
        print '(a,a)', "[lecture_mol2] File to read  = ", trim(filename)
        open(unit = fd, file = filename, iostat = ok, status = 'old')
        if ( ok /= 0 ) then
            print '(a,4x,a)', "Error during opening", filename
            stop 20
        end if

        ! First section read : get number of atoms
        lecture_fichier_mol2%num_atoms = 0
        do
            read (fd, '(a)', iostat = end), line
            if(end /= 0) then
                print '(a,4x,a)', "Error while reading", filename
                exit
            else
                if ( trim(line) == "@<TRIPOS>MOLECULE" ) then
                    ! Get molecule name
                    read (fd, '(a)', iostat = end), mol_name
                    ! Get atoms numbers
                    read (fd, '(i5)', iostat = end), num_atoms
                    ! Exit to section 2
                    exit
                end if
            end if
        end do

        ! Print information
        print '(a, a)', "[lecture_mol2] Molecule name : ", trim(mol_name)
        print '(a, i5)', "[lecture_mol2] Num of atoms  : ", num_atoms

        ! Read atom data

        ! Skip to atoms data
        do
            read (fd, '(a)', iostat = end), line
            if(end/=0)then
                print '(a, a)', "Could not find atom data of: ", filename
                stop 1
            else
                if ( trim(line) == "@<TRIPOS>ATOM" ) then
                    exit
                end if
            end if
        end do

        ! Read atom data
        allocate(atoms(num_atoms))

        do j = 1, num_atoms
            read (fd, '(i7, 1x, a8, 3f10.4, 1x, a8)', iostat = end), atom_num, atom_name, x, y, z, atom_type

            ! Normalize chemical element
            i = index(atom_type, '.')
            if ( i > 0 ) then
                i = i - 1
                atom_type = atom_type(:i)
            end if

            ! Store atom data
            atom_xyz%atom = trim(atom_type)
            atom_xyz%x = x
            atom_xyz%y = y
            atom_xyz%z = z

            atoms(j) = atom_xyz

            if (end/=0) then
                print '(a, a)', "Invalid mol2 format of ", filename
                stop 1
            end if
        end do

        ! Close file descriptor
        close(fd)

        ! Returns
        lecture_fichier_mol2%mol_name = mol_name
        lecture_fichier_mol2%num_atoms = num_atoms
        lecture_fichier_mol2%atoms = atoms

    end function lecture_fichier_mol2

end module lecture_mol2
