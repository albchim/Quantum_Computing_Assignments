MODULE matrices

    TYPE cmatrix
        INTEGER, dimension(2) :: dim
        DOUBLE COMPLEX, dimension(:,:), ALLOCATABLE :: elem, adj
        DOUBLE COMPLEX :: trace
    END TYPE cmatrix

    INTERFACE OPERATOR (.tr.)
        MODULE PROCEDURE matrix_trace
    END INTERFACE

    INTERFACE OPERATOR (.adj.)
        MODULE PROCEDURE matrix_adjoint
    END INTERFACE

    contains

    subroutine deall_cmatrix(M)
        type(cmatrix) :: M
        deallocate(M%elem)
        deallocate(M%adj)
    end subroutine

    ! Trace operator function
    FUNCTION matrix_trace(M) result(trace)
        type(cmatrix), intent(in) :: M
        double complex :: trace
        integer :: ii
        trace=0.
        ! Simple square matrix trace
        if (M%dim(1)==M%dim(2)) then
            do ii=1, M%dim(1)
                trace = trace + M%elem(ii,ii)
            end do
            return
        else
            print*,"Warning: Matrix is not square! Cannot compute trace"
            trace=0.
            return
        end if
    END FUNCTION

    ! Adjoint operator function
    FUNCTION matrix_adjoint(M) result(adjoint)
        type(cmatrix), intent(in) :: M
        double complex, dimension(:,:), allocatable :: adjoint
        integer :: ii
        
        allocate(adjoint(M%dim(1), M%dim(2)))
        do ii=1,M%dim(1)
            do jj=ii+1,M%dim(2)
                adjoint(jj,ii)=dconjg(M%elem(ii,jj))
                adjoint(ii,jj)=dconjg(M%elem(jj,ii))
            end do
            adjoint(ii,ii)=dconjg(M%elem(ii,ii))
        end do

    END FUNCTION

    subroutine check_hermitianity(H)
        ! Checks whether the matrix is hermitian
        type(cmatrix) :: H
        print*, "Checking hermitianity...."
        if (all((real(H%elem) - real(H%adj))<1e-5) .and. all((aimag(H%elem) - aimag(H%adj))<1e-5)) then
            print*, "  ---> Matrix is Hermitian"
        else 
            print*, "WARNING: ---> The matrix is not hermitian!"
            stop
        end if
    end subroutine

    subroutine write_matrix(M, string, file_number)
        character(40) :: string
        integer :: file_number
        integer :: ii, jj
        integer, dimension(2) :: dim
        type(cmatrix), intent(in) :: M
        
        write(file_number,*) string
        write(file_number, *) ""
        if (M%dim(1)/=M%dim(2)) then
            write(file_number,*) "MATRIX IS NOT SQUARE!"
        end if
        
        write(file_number,*) "[dim], ", M%dim
        write(file_number,*) "[trace], ", M%trace
        write(file_number,"(A8)", advance='no') "[elem], "
        do ii=1, M%dim(2)
            write(file_number,"(A17, i3, A2)",advance='no') "column[", ii, "], "
        end do
        write(file_number,*) "" ! allows new line
        do ii=1, M%dim(1)
            if (sum(imag(M%elem(ii,:)))==0.) then
                write(file_number,"(A4, i3, A2)",advance='no') "row[", ii, "],"
                do jj=1, M%dim(2)
                    write(file_number,"(F20.5, A)", advance='no') real(M%elem(ii,jj)),","
                    if (jj==M%dim(2)) then
                        ! stops line when row is finished
                        write(file_number,*) ""
                    end if
                end do
            else
                write(file_number,"(A4, i3, A2)",advance='no') "row[", ii, "],"
                do jj=1, M%dim(2)
                    write(file_number,'(F15.3,SP,F5.3,A)', advance='no') real(M%elem(ii,jj)), aimag(M%elem(ii,jj)), "i,"
                    if (jj==M%dim(2)) then
                        ! stops line when row is finished
                        write(file_number,*) ""
                    end if
                end do
            end if
        end do
        write(file_number, *) ""
        
    END subroutine

END MODULE matrices

MODULE quantum_states
    use matrices

    implicit none

    TYPE qstate
        integer :: n_subs ! number of subsystems (N)
        integer :: dim ! subsystem dimension (d)
        logical :: sep ! boolean for separable state
        double complex, dimension(:), allocatable :: coeff ! values stored N*d
        integer :: len ! lenght of the statevector
    END TYPE qstate

    contains

    subroutine deall_qstate(state)
        type(qstate) :: state
        deallocate(state%coeff)
    end subroutine

    function init_qstate(N, d, sep, debug) result(psi)
        ! initializes the quantum state filling it with random numbers
        ! it differs in case the state is separable or not

        type(qstate) :: psi
        integer, intent(in) :: N
        integer, intent(in) :: d
        logical, intent(in) :: sep
        logical, optional :: debug
        !dummy variables
        integer :: istat
        integer :: ii, jj
        double precision :: Re, Im, sum

        psi%n_subs = N
        psi%dim = d
        psi%sep = sep

        ! SEPARABLE CASE
        if (psi%sep .eqv. .TRUE.) then
            psi%len = N * d

            allocate(psi%coeff(psi%len), stat=istat)
            if (istat/=0) then
                print*, "Error in allocating the coefficients"
            end if

            ! init coefficients
            do ii=1, psi%n_subs
                sum = 0.
                do jj=1, psi%dim
                    call random_number(Re)
                    call random_number(Im)
                    ! actual element
                    psi%coeff(jj + (ii-1)*psi%dim) = cmplx(Re, Im)
                    ! calc normalization
                    sum = sum + abs(cmplx(Re, Im))**2
                end do
                ! normalize each subsystem to 1/N
                do jj=1, psi%dim
                    psi%coeff(jj + (ii-1)*psi%dim) = psi%coeff(jj + (ii-1)*psi%dim) / (sqrt(psi%n_subs*sum))
                end do
                
                if (present(debug) .and. debug .eqv. .TRUE.) then
                    print*, ""
                    print*, "--------------DEBUG CELL------------------"
                    print*, "Normalization of subsystem", ii
                    sum = 0.
                    do jj=1, psi%dim
                        sum = sum + abs(psi%coeff(jj + (ii-1)*psi%dim))**2
                    end do
                    print*, sum
                    print*, "------------------------------------------"
                    print*, ""
                end if
            end do
            
            if (present(debug) .and. debug .eqv. .TRUE.) then
                print*, ""
                print*, "--------------DEBUG CELL------------------"
                print*, "Normalization separable state:"
                sum = 0.
                do jj=1, psi%dim*psi%n_subs
                    sum = sum + abs(psi%coeff(jj))**2
                end do
                print*, sum
                print*, "------------------------------------------"
                print*, ""
            end if

        ! NON-SEPARABLE CASE
        else
            psi%len = d ** N
            allocate(psi%coeff(psi%len), stat=istat)
            if (istat/=0) then
                print*, "Error in allocating the coefficients"
            end if

            ! init coefficients
            do ii=1, psi%dim ** psi%n_subs
                call random_number(Re)
                call random_number(Im)
                ! actual element
                psi%coeff(ii) = cmplx(Re, Im)
            end do

            ! normalization coeff
            sum = 0.
            do ii=1, psi%dim ** psi%n_subs
                sum = sum + abs(psi%coeff(ii))**2
            end do

            psi%coeff = psi%coeff/sqrt(sum)

            if (present(debug) .and. debug .eqv. .TRUE.) then
                print*, ""
                print*, "--------------DEBUG CELL------------------"
                print*, "Normalization non-separable state:"
                sum = 0.
                do jj=1, psi%dim**psi%n_subs
                    sum = sum + abs(psi%coeff(jj))**2
                end do
                print*, sum
                print*, "------------------------------------------"
                print*, ""
            end if
        end if
    end function


    subroutine sep_to_gen(psi_out, debug)
        ! Writes the N*D states vector in D basis

        type(qstate), intent(inout) :: psi_out
        logical, optional :: debug
        double precision :: sum
        integer :: ii, jj, kk, qq, istat
        integer, dimension(:), allocatable :: coeff
        double complex, dimension(:), allocatable :: temp

        if (psi_out%len==psi_out%dim**psi_out%n_subs .and. psi_out%sep .eqv. .TRUE.) then
            print*, "The state seems to be already general..."
        end if

        ! change the lenght of the vector
        psi_out%len = psi_out%dim**psi_out%n_subs
        psi_out%sep = .FALSE.

        ! set temp vector to keep old values
        allocate(temp(psi_out%len), stat=istat)
        if (istat/=0) then
            print*, "Error in allocating the vector"
        end if
        temp = psi_out%coeff

        ! coefficients of the representation
        allocate(coeff(psi_out%n_subs), stat=istat) ! coefficients for each subsystem
        if (istat/=0) then
            print*, "Error in allocating the vector"
        end if

        ! create new allocation for the full vector
        deallocate(psi_out%coeff)
        allocate(psi_out%coeff(psi_out%dim**psi_out%n_subs), stat=istat)
        if (istat/=0) then
            print*, "Error in allocating the vector"
        end if
        
        
        sum = 0. ! needed to normalize
        psi_out%coeff=cmplx(1., 0.)
        do ii=1, psi_out%len !cycle over all the elements of the full vector
            kk = ii - 1 !due to fortran indexing
            do jj=1, psi_out%n_subs !cycle over the subsystems
                qq = psi_out%n_subs - jj
                coeff(qq+1) = int((kk-mod(kk,psi_out%dim**qq))/psi_out%dim**qq)+1 !+1 to readjust the indexing
                kk = kk - (coeff(qq+1)-1)*psi_out%dim**qq
            end do
            do jj=1, psi_out%n_subs
                psi_out%coeff(ii) = psi_out%coeff(ii) * temp(coeff(jj)+(jj-1)*psi_out%dim)
            end do
            sum = sum + abs(psi_out%coeff(ii))**2
        end do

        ! normalization
        psi_out%coeff = psi_out%coeff/sqrt(sum)

        deallocate(coeff)
        deallocate(temp)

        if (present(debug) .and. debug .eqv. .TRUE.) then
            print*, ""
            print*, "--------------DEBUG CELL------------------"
            print*, "Normalization enlarged state:"
            sum = 0.
            do jj=1, psi_out%len
                sum = sum + abs(psi_out%coeff(jj))**2
            end do
            print*, sum
            print*, "------------------------------------------"
            print*, ""
        end if

    end subroutine


    function init_density(psi, debug) result(density)
        ! outputs the density matrix from a given state vector
        ! it differs in case the state is separable or not

        type(qstate) :: psi
        type(cmatrix) :: density
        logical, optional :: debug
        integer :: ii, jj, istat

        if (psi%sep .eqv. .TRUE.)then
            call sep_to_gen(psi, debug)
        end if
            
        allocate(density%elem(psi%len, psi%len))
        if (istat/=0) then
            print*, "Error in allocating the matrix"
        end if
        do ii=1, psi%len
            do jj=1, psi%len
                density%elem(ii,jj) = psi%coeff(ii)*dconjg(psi%coeff(jj))
            end do
        end do

        density%dim = SHAPE(density%elem)
        density%trace = .tr.density
        density%adj = .adj.density

        if (present(debug) .and. debug .eqv. .TRUE.) then
            print*, ""
            print*, "--------------DEBUG CELL------------------"
            print*, "Density matrix Trace:"
            print*, density%trace
            call check_hermitianity(density)
            print*, "------------------------------------------"
            print*, ""
        end if
            
    end function


    function reduce_density_bipartite(density, sys_to_keep, debug) result(rdensity)
        ! computes the partial trace on a bipartite system density matrix

        type(cmatrix) :: density ! density matrix
        type(cmatrix) :: rdensity ! reduced density matrix
        integer :: sys_to_keep ! integer corresponding to the system one wants to keep (1 or 2)
        integer :: dim ! dimension of the subsystems (supposed to be equal)
        integer, dimension(2) :: temp
        logical, optional :: debug

        integer :: ii, jj, kk

        temp = SHAPE(density%elem)

        dim = INT(temp(1)/2)
        allocate(rdensity%elem(dim, dim))
        rdensity%elem = 0.
        
        if (sys_to_keep==1) then
            do ii=1, dim
                do jj=1, dim
                    do kk=1, dim
                        rdensity%elem(ii,jj) = rdensity%elem(ii,jj) + density%elem((kk-1)*dim + ii,(kk-1)*dim + jj)
                    end do
                end do
            end do

        elseif (sys_to_keep==2) then
            do ii=1, dim
                do jj=1, dim
                    do kk=1, dim
                        rdensity%elem(ii,jj) = rdensity%elem(ii,jj) + density%elem((ii-1)*dim + kk,(jj-1)*dim + kk)
                    end do
                end do
            end do
        
        else
            print*, "Wrong sys_to_keep input....exiting"
            stop
        end if

        rdensity%dim = (/dim,dim/)
        rdensity%trace = .tr.density
        rdensity%adj = .adj.rdensity

        if (present(debug) .and. debug .eqv. .TRUE.) then
            print*, ""
            print*, "--------------DEBUG CELL------------------"
            print*, "Reduced density matrix Trace:"
            print*, rdensity%trace
            call check_hermitianity(rdensity)
            print*, "------------------------------------------"
            print*, ""
        end if

    end function



END MODULE quantum_states

PROGRAM main

    use quantum_states

    implicit none

    type(qstate) :: psi
    type(cmatrix) :: density, rdensity
    integer :: N, d, jj
    logical :: sep, debug

    N = 2
    d = 2
    debug = .FALSE.

    print*, "Number of subsystems:", N
    print*, "Dimension of subsystem hilbert spaces:", d
    print*, ""

    ! Non separable case
    psi = init_qstate(N, d, .FALSE., debug)
    print*, "General State vector"
    print*, (psi%coeff(jj), jj=1, psi%len)
    print*, ""
    density = init_density(psi, debug)
    rdensity = reduce_density_bipartite(density, 2, debug)

    open(unit=100, file="density_gen.txt", status="REPLACE")
    call write_matrix(density, "Density_matrix_General_Pure", 100)
    call write_matrix(rdensity, "Reduced_density_matrix_General_Pure", 100)
    close(100)

    call deall_qstate(psi)
    call deall_cmatrix(density)
    call deall_cmatrix(rdensity)

    ! Separable case
    psi = init_qstate(N, d, .TRUE., debug)
    print*, "Separable State vector"
    print*, (psi%coeff(jj), jj=1, psi%len)
    print*, ""
    density = init_density(psi, debug)
    rdensity = reduce_density_bipartite(density, 2, debug)

    open(unit=200, file="density_sep.txt", status="REPLACE")
    call write_matrix(density, "Density_matrix_Separable_Pure", 200)
    call write_matrix(rdensity, "Reduced_density_matrix_Separable_Pure", 200)
    close(200)

    call deall_qstate(psi)
    call deall_cmatrix(density)
    call deall_cmatrix(rdensity)

END PROGRAM



