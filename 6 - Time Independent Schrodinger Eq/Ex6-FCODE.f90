MODULE matrices
    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!OLD CODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        double complex, dimension(:,:), allocatable :: eigenvect
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
            allocate(M%eigenvect(M%dim(1),M%dim(2)))
            M%eigenvect = 0.
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
!!!!!!!!!!!!!!!!!!!!!!!!!OLD CODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    
    subroutine onedim_ho(x_min, x_max, N, M, spacing)
        ! LEAVES 1/h**2 out in order to avoid precision problems
        ! to be multiplied by the eigenvectors
        ! derivative estimation of order h**4
        integer, intent(in) :: N
        double precision, intent(in) :: x_min, x_max
        double complex, dimension(:,:), allocatable, intent(inout) :: M
        double precision, intent(inout) :: spacing
        integer :: ii, io

        allocate(M(N,N))
        spacing=(x_max - x_min)/N
        ! filling first-Last row
        M(1,1) = cmplx(2 + (spacing**2)*(x_min**2), 0)
        M(1,2) = cmplx((-1), 0)
        M(N,N) = cmplx((2 + (spacing**2)*((x_min + N*spacing)**2)), 0)
        M(N,N-1) = cmplx((-1), 0)
        ! filling rows
        do ii=2, N-1
            M(ii,ii) = cmplx((2 + (spacing**2)*((x_min + (ii-1)*spacing)**2)), 0)
            M(ii,ii-1) = cmplx((-1), 0)
            M(ii,ii+1) = cmplx((-1), 0)
        end do
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
                print*, "  ---> Diagonalization successfully converged!"
            elseif (info<0) then
                print*, "WARNING: ---> Value of the", -info,"th element is illegal"
                stop
            else
                print*, "WARNING: ---> The diagonalizing procedure did not converge and stopped at element:", info
            end if
        else
            print*, "WARNING: ---> Matrix is not square! Exiting..."
            stop
        end if
    end subroutine

END MODULE matrices


PROGRAM MAIN

    use matrices

    implicit none

    type(cmatrix) :: M
    integer :: N, n_eigenvects
    double complex, dimension(:,:), allocatable :: temp
    double precision :: x_min, x_max, spacing
    double complex, dimension(:,:), allocatable :: temp_eigvect
    double precision, dimension(:), allocatable :: temp_eig

    integer :: ii, jj
    !Checks i/o status
    integer :: istat

    print*, "Enter the number of points contained in the mesh (INTEGER NUMBER): "
    read(*,*) N

    x_min = -5
    x_max = 5

    call onedim_ho(x_min, x_max, N, temp, spacing)
    print*, "Spacing: ", spacing
    M=init(temp)

    ! compute eigenvalues
    call compute_eigenvals(M, temp_eig, temp_eigvect)

    M%eigenvals = temp_eig/spacing**2
    M%eigenvect = temp_eigvect/sqrt(spacing)

    ! writes to file
    open(unit=100, file='eigenvalues.txt')!, status="REPLACE")
    do ii=1, N
        write(100, "(F15.7)", iostat=istat) M%eigenvals(ii)
        if (istat /= 0) then
            print*, "WARNING: ---> Error in wrinting, shutting down..." 
            stop
        end if
    end do
    close(100)

    n_eigenvects = 4
    open(unit=200, file="eigenvectors.txt")!, status="REPLACE")
    do ii=1, N
        write(200,*) x_min+(ii-1)*spacing, ((real(M%eigenvect(ii, jj))), jj = 1, n_eigenvects)
    end do
    close(200)

    deallocate(M%elem)
    deallocate(M%eigenvals)
    deallocate(M%eigenvect)
    deallocate(temp)

END PROGRAM