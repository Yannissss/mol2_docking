program mol2_amarrage
    use chargeur_covalence
    use lecture_mol2

    implicit none

    ! Variables
    character(len=128) :: filename, line

    integer :: i, j, k, ok, end, population, max_iter
    real :: merge_rate, mutate_rate
    type(Covalence) :: cov
    type(CovTable) :: table
    type(Molecule) :: ligand, site
    type(AtomXYZ) :: atom_xyz

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
    print '(a, a, a, i4, a)', "[mol2_amarrage] Ligand molecule : ", ligand%mol_name, " / ", ligand%num_atoms, " atoms"
    print '(a, a, a, i4, a)', "[mol2_amarrage] Site molecule : ", site%mol_name, " / ", site%num_atoms, " atoms"
    print '(a, i4)', "[mol2_amarrage] population  = ", population
    print '(a, f3.2)', "[mol2_amarrage] merge_rate  = ", merge_rate
    print '(a, f3.2)', "[mol2_amarrage] mutate_rate = ", mutate_rate
    print '(a, i4)', "[mol2_amarrage] max_iter    = ", max_iter

    ! Wrap-up

end program mol2_amarrage
