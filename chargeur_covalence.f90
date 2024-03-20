module chargeur_covalence
    implicit none

    type, public :: Covalence
        character(len=2) :: atom
        integer :: atom_num, simple, double, triple
    end type Covalence

    type, public :: CovTable
        integer :: num_cov
        type(Covalence), dimension(120) :: data
    contains
        procedure get_cov_radii
    end type CovTable

contains

    type(Covalence) function get_cov_radii(self, elem)
        class(CovTable), intent(in) :: self
        character(len=*), intent(in) :: elem
        integer :: i

        do i = 1, self%num_cov
            if ( self%data(i)%atom == elem ) then
                get_cov_radii%atom     = elem
                get_cov_radii%atom_num = i
                get_cov_radii%simple   = self%data(i)%simple
                get_cov_radii%double   = self%data(i)%double
                get_cov_radii%triple   = self%data(i)%triple
                exit
            end if
        end do

    end function get_cov_radii

    type(CovTable) function charge_covalence(fileName)
        ! Signature & variables
        character(len=128) :: fileName, line
        integer :: covfile, ok, end, i, j, k
        logical :: test_comment

        character(len=2) :: cov_atom
        integer :: cov_atom_num, cov_simple, cov_double, cov_triple

        type(CovTable) :: table
        integer :: num_cov

        ! Open file
        covfile = 10
        print '(a,a)', "[chargeur_covalence] File to read = ", trim(fileName)
        open(unit = covfile, file = fileName, iostat = ok, status = 'old')
        if ( ok /= 0 ) then
            print '(a,4x,a)', "Error during opening", fileName
            stop 20
        end if

        ! skip first section
        do i=1,10
            read (covfile, '(a)', iostat = end), line
            if(end/=0)then
                exit
            end if
        end do

        ! read data
        i = 1
        do
            read (covfile, '(i3, a2, i3, i3, i3)', iostat = end), &
                cov_atom_num, cov_atom, cov_simple, cov_double, cov_triple
            if (end /= 0) then
                num_cov = i - 1 ! Store num of covalence read
                exit
            else
                table%data(i)%atom_num = cov_atom_num
                table%data(i)%atom     = cov_atom
                table%data(i)%simple   = cov_simple
                table%data(i)%double   = cov_double
                table%data(i)%triple   = cov_triple
                i = i + 1
            end if
        end do

        table%num_cov = num_cov
        charge_covalence = table

    end function charge_covalence

end module chargeur_covalence
