MODULE matrices
    implicit none

    !contants
    double precision :: pi=3.14159265359

    ! Type definition
    TYPE cmatrix
        INTEGER, dimension(2) :: dim
        DOUBLE COMPLEX, dimension(:,:), ALLOCATABLE :: elem, adj
        DOUBLE COMPLEX :: trace
        DOUBLE COMPLEX :: det
        ! Optional type object to store things tidy
        double precision, dimension(:), allocatable :: eigenvals
    END TYPE cmatrix 

    INTERFACE OPERATOR (.tr.)
        MODULE PROCEDURE matrix_trace
    END INTERFACE

    INTERFACE OPERATOR (.adj.)
        MODULE PROCEDURE matrix_adjoint
    END INTERFACE

    INTERFACE OPERATOR (.det.)
        MODULE PROCEDURE matrix_det
    END INTERFACE

    contains

    ! Initialization function
    FUNCTION init(temp) result(M)

        type(cmatrix) :: M
        integer, dimension(2) :: dim
        double complex, dimension(:,:), allocatable :: temp
        dim=SHAPE(temp)
        allocate(M%elem(dim(1),dim(2)))
        M%elem=temp
        M%dim=dim
        M%trace=.tr.M
        !M%det=.det.M  ! Avoid to speed up the computation
        M%adj=.adj.M
        if (M%dim(1)==M%dim(2)) then
            allocate(M%eigenvals(M%dim(1)))
            M%eigenvals = 0.
        end if
        !deallocate(temp)
        return
    END FUNCTION

    FUNCTION conjugate(matrix) result(conj)
        integer, dimension(2) :: dim
        integer :: ii, jj
        double complex, dimension(:,:), allocatable :: matrix
        double complex, dimension(:,:), allocatable :: conj
        dim=SHAPE(matrix)
        allocate(conj(dim(1),dim(2)))
        do ii=1,dim(1)
            do jj=1,dim(2)
                conj(ii,jj)=cmplx(real(matrix(ii,jj)), -imag(matrix(ii,jj)))
            end do
        end do
        return
    END FUNCTION

    FUNCTION transpose(matrix) result(transp)
        double complex, dimension(:,:), allocatable :: matrix
        double complex, dimension(:,:), allocatable :: transp
        integer :: ii, jj
        integer, dimension(2) :: dim
        dim=SHAPE(matrix)
        ! Square matrix case
        if (dim(1)==dim(2)) then 
            allocate(transp(dim(1),dim(2)))
            do ii=1,dim(2)
                do jj=ii+1,dim(1)
                    transp(jj,ii)=matrix(ii,jj)
                    transp(ii,jj)=matrix(jj,ii)
                end do
                transp(ii,ii)=matrix(ii,ii)
            end do
            return
        ! Non-square matrix case
        else
            allocate(transp(dim(2), dim(1)))
            ! select the bigger dimention to reduce loop repetitions
            if (dim(1) > dim(2)) then
                do ii=1,dim(2)
                    do jj=ii,dim(1)
                        transp(ii,jj)=matrix(jj,ii)
                        transp(jj,ii)=matrix(ii,jj)
                    end do
                end do
            else !if (dim(2) > dim(1))) then
                do ii=1,dim(1)
                    do jj=ii,dim(2)
                        transp(ii,jj)=matrix(jj,ii)
                        transp(jj,ii)=matrix(ii,jj)
                    end do
                end do
            end if
            return
        end if
    END FUNCTION
        
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
        ! Square matrix case
        if (M%dim(1)==M%dim(2)) then
            allocate(adjoint(M%dim(1), M%dim(2)))
            adjoint=transpose(conjugate(M%elem))
            return
        ! Non square case
        else
            allocate(adjoint(M%dim(2), M%dim(1))) ! overrides allocation in main porgram
            adjoint=transpose(conjugate(M%elem))
            return
        end if
    END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!OPTIONAL FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine LUdecon(matrix, L, U, debug)
        ! No need to allocate L,U when calling the subroutine
        integer, optional :: debug
        logical :: d
        integer, dimension(2) :: dim
        double complex, dimension(:,:), allocatable :: matrix
        double complex, dimension(:,:), allocatable :: temp, test ! test
        double complex, dimension(:,:), allocatable :: L, U
        double complex :: coeff
        integer :: ii, jj, kk, n

        dim=SHAPE(matrix)
        if (dim(1)==dim(2)) then
            allocate(test(dim(1),dim(2))) ! test
            allocate(temp(dim(1),dim(2)))
            allocate(L(dim(1),dim(2)))
            allocate(U(dim(1),dim(2)))

            temp=matrix
            L=0.0
            U=0.0
            n=dim(1)
            
            ! L matrix (lower triangular matrix) calculation
            do kk=1, n-1
                do ii=kk,n ! Wrong index original code
                   coeff=matrix(ii,kk)/matrix(kk,kk)
                   L(ii,kk) = coeff
                   do jj=kk+1,n
                      temp(ii,jj) = temp(ii,jj)-coeff*temp(kk,jj)
                   end do
                end do
            end do
            ! Set diagonal to 1 to avoid precision problems
            do ii=1,n
                L(ii,ii) = 1.0
            end do

            ! U matrix is the upper triangular part of A
            do jj=1,n
                do ii=1,jj
                    U(ii,jj) = temp(ii,jj)
                end do
            end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !DEBUG verify the decomposition
            if (present(debug)) then
                print*, "Checking the decomposition element by element"
                test=matmul(L,U) !inv used as a temp object to store the LU product
                do ii=1, n
                    do jj=1, n
                        if (real(test(ii,jj))-real(temp(ii,jj))<1e-10) then
                            print*,"Ok", ii, jj
                        else
                            print*,"Wrong!", ii, jj
                        end if
                    end do
                end do
            end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            deallocate(test)
            deallocate(temp)
        else
            print*, "Warning: Matrix is not square! Cannot LU decompose"
        end if
    END subroutine

    FUNCTION matrix_det(M) result(det)
        integer :: debug ! debug
        integer :: ii
        type(cmatrix), intent(in) :: M
        double complex, dimension(:,:), allocatable :: L
        double complex, dimension(:,:), allocatable :: U
        double complex :: det

        if (M%dim(1)==M%dim(2)) then
            ! No need to allocate L,U when calling the subroutine
            call LUdecon(M%elem, L, U)!, debug)  ! Uncomment debug to test LUdecomposition
            do ii=1, M%dim(1)
                det = det + U(ii,ii)
            end do
            deallocate(L)
            deallocate(U)
            return
        else
            print*, "Matrix is not square!"
            return
        end if
    END function
!!!!!!!!!!!!!!!!!!!!!!!!!OPTIONAL FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    FUNCTION random_hermitian(N, gauss) result(H)
        integer :: N
        ! Generates non symmetric hermitian random matrix
        double complex :: H(N,N)
        ! Temp values for random entries
        double precision :: a, b, c, d
        ! Optional logical for gaussian distribution
        logical, optional :: gauss
        integer :: ii, jj

        ! Box-Muller for gaussian random numbers
        if (present(gauss) .and. gauss.eqv..TRUE.) then
            do ii=1, N
                call random_number(a)
                ! Sets real random diagonal entries
                H(ii,ii) = cmplx(a,0)
                do jj=ii+1, N
                    call random_number(c)
                    call random_number(d)
                    a = sqrt(-2*log(c))*cos(2*pi*d)
                    b = sqrt(-2*log(c))*sin(2*pi*d)
                    H(ii,jj)=cmplx(a, b)
                    H(jj,ii)=cmplx(a, -b)
                end do
            end do

        ! Uniformly distributed
        else
            do ii=1, N
                call random_number(a)
                ! Sets real random diagonal entries
                H(ii,ii) = cmplx(a,0)
                do jj=ii+1, N
                    call random_number(a)
                    call random_number(b)
                    H(ii,jj)=cmplx(a, b)
                    H(jj,ii)=cmplx(a, -b)
                end do
            end do
        end if
        return
    END function


    subroutine check_hermitianity(H)
        ! Checks whether the matrix is hermitian
        type(cmatrix) :: H
        print*, "Checking hermitianity...."
        if (all((real(H%elem) - real(H%adj))<1e-15) .and. all((aimag(H%elem) - aimag(H%adj))<1e-15)) then
            print*, "  ---> Matrix is Hermitian"
        else 
            print*, "WARNING: ---> The matrix is not hermitian!"
            stop
        end if
    END SUBROUTINE


    function random_real_diag(N, gauss) result(D)

        integer :: N
        logical, optional :: gauss
        double complex :: D(N,N)
        double precision :: a, b
        integer :: ii
        D=0.
        if (present(gauss) .and. gauss.eqv..TRUE.) then
            do ii=1, N
                call random_number(a)
                call random_number(b)
                a = sqrt(-2*log(a))*cos(2*pi*b)
                D(ii,ii)=cmplx(a,0.)
            end do
        else
            do ii=1, N
                call random_number(a)
                D(ii,ii)=cmplx(a,0.)
            end do
        end if
        return
    end function


    function compute_eigenvals(H) result(eig)
        ! Computes the eigenvalues of the matrix
        type(cmatrix) :: H
        ! used to check convergence of the output of zheev
        integer :: info, lwork
        double complex, dimension(:), allocatable :: work
        double precision, dimension(:), allocatable :: rwork
        double complex, dimension(:,:), allocatable :: temp_mat
        double precision, dimension(:), allocatable :: eig
        
        ! Checks for square matrix
        if (H%dim(1)==H%dim(2)) then
            
            ! allocates temporary matrix copy and eigenvalues vector
            allocate(temp_mat(H%dim(1),H%dim(2)))
            allocate(eig(H%dim(1)))
            temp_mat=H%elem
            eig=0. ! initialize to zero
            ! needed inputs for zheev subroutine
            lwork= 2*(H%dim(1))
            allocate(work(lwork))
            allocate(rwork(lwork))

            ! Computes the eigenvalues only ('N'), 'U' overwrites the upper triangular matrix in place
            ! of the input matrix.
            call zheev('N', 'U', H%dim(1), temp_mat, H%dim(1), eig, work, lwork, rwork, info)

            ! Checks the output status using the INFO variable
            if (info==0) then
                print*, "  ---> Diagonalization successfully converged!"
            elseif (info<0) then
                print*, "WARNING: ---> Value of the", -info,"th element is illegal"
                stop
            else
                print*, "WARNING: ---> The diagonalizing procedure did not converge and stopped at element:", info
            end if

            deallocate(temp_mat)
            deallocate(work)
            deallocate(rwork)

        else
            print*, "WARNING: ---> Matrix is not square! Exiting..."
            stop
        end if
    end function


    function compute_spacings(eig, N) result(norm_spacings)

        double precision, dimension(:), allocatable :: eig
        double precision, dimension(:), allocatable :: norm_spacings
        double precision :: mean_spacing
        integer :: ii, N

        allocate(norm_spacings(N-1))
        mean_spacing = 0.

        do ii=1, (N - 1)
            norm_spacings(ii) = eig(ii+1) - eig(ii) !computes space (NOT normalized)
            mean_spacing = mean_spacing + norm_spacings(ii) !starts computing mean spacing
        end do
        mean_spacing = mean_spacing /  (N - 1)
        norm_spacings = norm_spacings / mean_spacing
    
    end function


END MODULE matrices


PROGRAM MAIN

    use matrices

    type(cmatrix) :: H, D
    integer :: N
    double complex, dimension(:,:), allocatable :: temp
    double precision, dimension(:), allocatable :: space_herm, space_diag
    double precision :: mean_spacing

    integer :: ii
    character(20) :: out_h, out_d
    !Checks i/o status
    integer :: istat

    ! stores file names
    out_h = 'out_h.txt'
    out_d = 'out_d.txt'

    print*, "Enter integer for dimension of the square matrix: "
    read(*,*) N
    print*, "Matrix dimension is NxN with N =", N
    !print*, ""
    !print*, "The spacings will be saved at:"
    !print*, out_h, out_d
    !print*, ""

    open(unit=100, file=out_h, action='WRITE', access="sequential", position='append')
    open(unit=200, file=out_d, action='WRITE', access="sequential", position='append')

    !print*, "*HERMITIAN CASE-------------------------------------------"

    !!!!!!!!!RANDOM HERMITIAN!!!!!!!!!
    temp = random_hermitian(N, gauss=.TRUE.)
    H=init(temp)
    call check_hermitianity(H)

    ! compute eigenvalues
    H%eigenvals = compute_eigenvals(H)
    ! compute spacings
    space_herm = compute_spacings(H%eigenvals, N)
    ! writes to file
    do ii=1, N-1
        write(100, "(F15.7)", iostat=istat) space_herm(ii)
        if (istat /= 0) then
            print*, "WARNING: ---> Error in wrinting, shutting down..." 
            stop
        end if
    end do

    !print*, "*REAL DIAGONAL CASE---------------------------------------"

    !!!!!!!!!RANDOM REAL DIAGONAL!!!!!
    temp = random_real_diag(N, gauss=.TRUE.)
    D=init(temp)

    ! compute eigenvalues
    D%eigenvals = compute_eigenvals(D)
    ! compute spacings
    space_diag = compute_spacings(D%eigenvals, N)
    ! writes to file
    do ii=1, N-1
        write(200, "(F15.7)", iostat=istat) space_diag(ii)
        if (istat /= 0) then
            print*, "WARNING: ---> Error in wrinting, shutting down..." 
            stop
        end if
    end do

    deallocate(H%elem)
    deallocate(D%elem)
    deallocate(H%eigenvals)
    deallocate(D%eigenvals)
    deallocate(space_diag)
    deallocate(space_herm)

END PROGRAM
