MODULE matrices

    TYPE cmatrix
        INTEGER, dimension(2) :: dim
        DOUBLE COMPLEX, dimension(:,:), ALLOCATABLE :: elem, adj
        DOUBLE COMPLEX :: trace
        double precision, dimension(:), allocatable :: eigenvals
        double complex, dimension(:,:), allocatable :: eigenvects
    END TYPE cmatrix

    INTERFACE OPERATOR (.tr.)
        MODULE PROCEDURE matrix_trace
    END INTERFACE

    INTERFACE OPERATOR (.adj.)
        MODULE PROCEDURE matrix_adjoint
    END INTERFACE

    contains

    subroutine allocate_cmatrix(matrix, N, M, opt)
        type(cmatrix) :: matrix
        integer :: N,M, istat
        logical, optional :: opt

        matrix%dim = (/N,M/)
        allocate(matrix%elem(N,M), stat= istat)
        if (istat/=0) then
            print*, "Error in allocating the matrix elements"
        end if

        if (present(opt) .and. opt .eqv. .TRUE.)then
            allocate(matrix%adj(M,N), stat= istat)
            if (istat/=0) then
                print*, "Error in allocating the matrix elements"
            end if
        end if
    end subroutine


    subroutine deallocate_cmatrix(M, opt)
        type(cmatrix) :: M
        logical, optional :: opt
        deallocate(M%elem, stat= istat)
        if (istat/=0) then
            print*, "Error in deallocating the matrix elements"
        end if
        if (present(opt) .and. opt .eqv. .TRUE.)then
            deallocate(M%adj, stat= istat)
            if (istat/=0) then
                print*, "Error in deallocating the matrix elements"
            end if
            deallocate(M%eigenvals, stat= istat)
            if (istat/=0) then
                print*, "Error in deallocating the matrix elements"
            end if
            deallocate(M%eigenvects, stat= istat)
            if (istat/=0) then
                print*, "Error in deallocating the matrix elements"
            end if
        end if
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
        double complex, dimension(M%dim(2), M%dim(1)) :: adjoint
        integer :: ii
        
        do ii=1,M%dim(1)
            do jj=1,M%dim(2)
                adjoint(jj,ii)=dconjg(M%elem(ii,jj))
            end do
        end do
    END FUNCTION

    subroutine check_hermitianity(H)
        ! Checks whether the matrix is hermitian
        type(cmatrix) :: H
        print*, "Checking hermitianity...."
        if (all((real(H%elem) - real(H%adj))<1e-10) .and. all((aimag(H%elem) - aimag(H%adj))<1e-10)) then
            print*, "  ---> Matrix is Hermitian"
        else 
            print*, "WARNING: ---> The matrix is not hermitian!"
            stop
        end if
    end subroutine


    subroutine compute_eigenvals(H, eig, eigve)
        ! Computes the eigenvalues of the matrix
        type(cmatrix) :: H
        ! used to check convergence of the output of zheev
        integer :: info, lwork
        double complex, dimension(:), allocatable :: work
        double precision, dimension(:), allocatable :: rwork
        double complex, dimension(:,:), allocatable, intent(out) :: eigve
        double precision, dimension(:), allocatable, intent(out) :: eig
        
        ! Checks for square matrix
        if (H%dim(1)==H%dim(2)) then
            
            allocate(eigve(H%dim(1),H%dim(1)))
            allocate(eig(H%dim(1)))
            eigve=H%elem
            eig=0. ! initialize to zero

            ! compute optimal lwork
            lwork=-1
	        allocate(work(1))
	        allocate(rwork(max(1, 3*H%dim(1)-2)))
            call zheev('V', 'U', H%dim(1), eigve, H%dim(1), eig, work, lwork, rwork, info)
            ! optimal lwork
            lwork = int(real(work(1)))
	        deallocate(work)
            deallocate(rwork)
            
            ! re initialize the outputs
            eigve=H%elem
            eig=0.
	        ! actual diag
	        allocate(work(max(1,lwork)))
	        allocate(rwork(max(1, 3*H%dim(1)-2)))
	        call zheev('V', 'U', H%dim(1), eigve, H%dim(1), eig, work, lwork, rwork, info)

	        deallocate(work)
            deallocate(rwork)

            ! Checks the output status using the INFO variable
            if (info==0) then
                !print*, "  ---> Diagonalization successfully converged!"
            elseif (info<0) then
                print*, "WARNING: ---> Value of the", -info,"th element is illegal"
                stop
            else
                print*, "WARNING: ---> The diagonalizing procedure did not converge and stopped at element:", info
                stop
            end if
        else
            print*, "WARNING: ---> Matrix is not square! Exiting..."
            stop
        end if
    end subroutine

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

MODULE ising
    use matrices

    implicit none

    contains

    function tensor_product(A, B)result(M)
        type(cmatrix) :: A, B, M
        integer :: ii, jj, kk, qq

        call allocate_cmatrix(M, A%dim(1)*B%dim(1), A%dim(2)*B%dim(2))
        do ii=1, B%dim(1)
            do jj=1, B%dim(2)
                do kk=1, A%dim(1)
                    do qq=1, A%dim(2)
                        M%elem(kk+(ii-1)*A%dim(1), qq+(jj-1)*A%dim(2)) = B%elem(ii,jj)*A%elem(kk,qq)
                    end do
                end do
            end do
        end do   
    end function


    function field_interaction_op(mat, index, N)result(self)
        type(cmatrix) :: mat, self, idn, temp
        integer :: index, ii, N

        ! Generate identity matrix
        call allocate_cmatrix(idn, mat%dim(1), mat%dim(2))
        idn%elem=0.
        do ii=1, mat%dim(1)
            idn%elem(ii,ii) = cmplx(1., 0.)
        end do

        call allocate_cmatrix(self, mat%dim(1), mat%dim(2))

        if (index==1) then
            self%elem = mat%elem
        else
            self%elem = idn%elem
        end if
            
        ! Perform tensor product
        do ii=2, N
            if (ii==index) then
                temp = tensor_product(self, mat)
            else
                temp = tensor_product(self, idn)
            end if

            call deallocate_cmatrix(self)
            call allocate_cmatrix(self, mat%dim(1)**ii, mat%dim(2)**ii)

            self%elem = temp%elem
            call deallocate_cmatrix(temp)
        end do
    end function


    function pair_interaction_op(mat, index1, index2, N)result(pair)
        type(cmatrix) :: mat, pair, idn, temp
        integer :: index1, index2, ii, N

        ! Generate identity matrix
        call allocate_cmatrix(idn, mat%dim(1), mat%dim(2))
        idn%elem=0.
        do ii=1, mat%dim(1)
            idn%elem(ii,ii) = cmplx(1., 0.)
        end do

        call allocate_cmatrix(pair, mat%dim(1), mat%dim(2))

        if (index1==1 .or. index2==1) then
            pair%elem = mat%elem
        else
            pair%elem = idn%elem
        end if
            
        ! Perform tensor product
        do ii=2, N
            if (ii==index1 .or. ii==index2) then
                temp = tensor_product(pair, mat)
            else
                temp = tensor_product(pair, idn)
            end if

            call deallocate_cmatrix(pair)
            call allocate_cmatrix(pair, mat%dim(1)**ii, mat%dim(2)**ii)

            pair%elem = temp%elem
            call deallocate_cmatrix(temp)
        end do
    end function


    subroutine real_RG(H, N, subsystems, lambda, pauli_x, pauli_z)
        type(cmatrix), intent(inout) :: H
        type(cmatrix), intent(in) :: pauli_x, pauli_z
        double precision, intent(in) :: lambda
        integer, intent(in) :: N, subsystems
        type(cmatrix) :: subint1, subint2, tempdouble, idnn, temp
        type(cmatrix) :: P, tempP
        integer :: ii, dd

        ! REAL-RG ALGORITHM

        ! Pair interaction part between subsystems
        subint1 = field_interaction_op(pauli_z, N, N)
        subint2 = field_interaction_op(pauli_z, 1, N)

        ! Generate identity matrix (dim N)
        call allocate_cmatrix(idnn, 2**N, 2**N)
        idnn%elem=0.
        do ii=1, 2**N
            idnn%elem(ii,ii) = cmplx(1., 0.)
        end do

        call allocate_cmatrix(P, 2**(2*N), 2**N, .TRUE.)
        call allocate_cmatrix(tempP, 2**(2*N), 2**N)

        do dd=1, subsystems
            if (MOD(dd, 20)==0) then
                print*, "Simulated particles number N=", N, "**", dd
            end if
            call allocate_cmatrix(tempdouble, 2**(2*N), 2**(2*N), .TRUE.)


            ! Field interaction
            ! first part
            temp = tensor_product(H, idnn)
            tempdouble%elem = temp%elem
            call deallocate_cmatrix(temp)
            !second part
            temp = tensor_product(idnn, H)
            tempdouble%elem = tempdouble%elem + temp%elem
            call deallocate_cmatrix(temp)

            ! Pair interaction part
            temp = tensor_product(subint1, subint2)
            tempdouble%elem = tempdouble%elem + lambda*temp%elem
            call deallocate_cmatrix(temp)

            ! Compute eigenvectors
            call compute_eigenvals(tempdouble, tempdouble%eigenvals, tempdouble%eigenvects)

            ! Projector
            P%elem = tempdouble%eigenvects(:, 1:2**N)
            P%adj = .adj.P

            ! Project hamiltonian
            tempP%elem = matmul(tempdouble%elem, P%elem)
            H%elem = matmul(P%adj, tempP%elem)

            ! Project interaction subsystem terms
            temp = tensor_product(idnn, subint1)
            tempP%elem = matmul(temp%elem, P%elem)
            subint1%elem = matmul(P%adj, tempP%elem)
            call deallocate_cmatrix(temp)

            temp = tensor_product(subint2, idnn)
            tempP%elem = matmul(temp%elem, P%elem)
            subint1%elem = matmul(P%adj, tempP%elem)
            call deallocate_cmatrix(temp)

            call deallocate_cmatrix(tempdouble, .TRUE.)

        end do

        call deallocate_cmatrix(subint1)
        call deallocate_cmatrix(subint2)
        call deallocate_cmatrix(idnn)
        call deallocate_cmatrix(P)
        deallocate(P%adj)
        call deallocate_cmatrix(tempP)

    end subroutine

END MODULE ising

PROGRAM main

    use ising

    implicit none

    type(cmatrix) :: pauli_x, pauli_z, H, tempH
    double precision :: lambda
    integer :: N, subsystems, ii, jj, istat, l_resolution

    ! Inputs!!!!!!!!
    print*, "Enter the number of qbits for each subsystem (INTEGER NUMBER): "
    read(*,*) N
    print*, "Enter the number of desired subsystems to add (INTEGER NUMBER):"
    read(*,*) subsystems

    lambda = 0
    l_resolution = 50

    call allocate_cmatrix(pauli_x, 2, 2)
    call allocate_cmatrix(pauli_z, 2, 2)

    ! Initialize Pauli matrices
    pauli_x%elem(1,1)=(0.,0.)
	pauli_x%elem(1,2)=(1.,0.)
	pauli_x%elem(2,1)=(1.,0.)
	pauli_x%elem(2,2)=(0.,0.)

	pauli_z%elem(1,1)=(1.,0.)
	pauli_z%elem(1,2)=(0.,0.)
	pauli_z%elem(2,1)=(0.,0.)
    pauli_z%elem(2,2)=(-1.,0.)

    open(unit=50, file="l_eigenvalues.txt", status="REPLACE")
    
    do jj=1, l_resolution+1

        print*, "Iteration with lambda=", lambda

        ! Hamiltonian construction
        call allocate_cmatrix(H, 2**N, 2**N, .TRUE.)
        H%elem = 0.
        
        do ii=1, N
            tempH = field_interaction_op(pauli_z, ii, N)

            H%elem = H%elem + tempH%elem
            call deallocate_cmatrix(tempH)
        end do

        do ii=1, N-1
            tempH = pair_interaction_op(pauli_x, ii, ii+1, N)

            H%elem = H%elem + lambda*tempH%elem
            call deallocate_cmatrix(tempH)
        end do

        ! REAL-RG ALGORITHM
        call real_RG(H, N, subsystems, lambda, pauli_x, pauli_z)

        call compute_eigenvals(H, H%eigenvals, H%eigenvects)

        ! write to text file
        write(50,*, iostat=istat) lambda, (H%eigenvals(ii), ii=1, 1)
        if (istat /= 0) then
            print*, "WARNING: ---> Error in wrinting, shutting down..." 
            stop
        end if

        ! update lambda
        lambda = lambda + 5./(FLOAT(l_resolution))
        ! deallocate the old hamiltonian
        call deallocate_cmatrix(H, .TRUE.)

    end do

    call deallocate_cmatrix(pauli_x)
    call deallocate_cmatrix(pauli_z)



END PROGRAM