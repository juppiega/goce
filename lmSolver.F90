#include "fintrf.h"

module lmSolverModule
implicit none

! Interface for function returning a vector of (weighted) residuals.
interface
    function costFunction (x)
        implicit none
        real(kind = 8), intent(in) :: x(:)
        real(kind = 8), allocatable :: costFunction(:)
    end function
end interface

contains

subroutine lmSolve(FUN, X0, tolX, tolFun, tolOpt, lambda0, maxFuncEvals, maxIter, solution, &
                   funVec, firstOrderOptAtSolution, JacobianAtSolution, JTWJ, exitFlag)
    implicit none
    ! Input arguments
    procedure(costFunction) :: FUN ! Cost function. Returns a vector of (weighted) residuals.
    real(kind = 8), intent(in) :: X0(:)
    real(kind = 8), optional :: tolX, tolFun, tolOpt, lambda0 ! X-tolerance (1E-6), function tolerance (1E-6), first-order opt. tolerance (1E-14)
                                                              ! (all tolerances are optional. TolX and TolFun are relative, TolOpt is absolute), initial damping factor (1E-2) (optional).
    integer, optional :: maxFuncEvals, maxIter ! Maximum number of function evaluations (1000 * numVars), maximum iterations (100).

    ! Output arguments
    real(kind = 8), allocatable, optional, intent(out) :: solution(:), funVec(:), JacobianAtSolution(:,:),JTWJ(:,:) ! Final X, final function value, Jacobian at solution.
    real(kind = 8), optional, intent(out) :: firstOrderOptAtSolution ! Infinity norm of gradient at solution.
    integer, optional, intent(out) :: exitFlag ! Exit statue. -2 = Exceeded maxFuncEvals, -1 = Exceeded maxIter,
                                     ! 1 = First order optimality smaller than tolOpt, 2 = norm(step) < tolX * norm(X), 3 = changeSumSq < tolFun * oldSumSq.

    ! Local variables
    logical :: successfulStep, done, containsInfs
    real(kind = 8), allocatable :: trialX(:), trialCostFunc(:), gradient(:), A(:), JTJ(:), step(:), costFuncVec(:)
    real(kind = 8), allocatable :: X(:), Jacobian(:,:), JTJ_diag(:)
    integer, allocatable :: AdiagElem(:)
    real(kind = 8) :: lambda, firstOrderOpt
    integer :: k, iter, numFuncEvals, flag
    integer(kind = 8) :: ierr, one = 1, numVars, numResid
    integer(kind = 4) :: mexPrintf, mexStat, mexEvalString, mexCallMATLAB
    character(len = 400) :: line, numChar
    character(len = 1) :: lower
    mwPointer mxCreateString

    lower = 'L'

    write(line, '(a)') 'iter  Fcount    sumSq     firstOrdOpt     lambda     norm(step)'
    mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(line), 'disp')

    ! Set default arguments.
    if (.not. present(tolX)) tolX = 1E-6
    if (.not. present(tolFun)) tolFun = 1E-6
    if (.not. present(tolOpt)) tolOpt = 1E-14
    if (.not. present(lambda0)) lambda0 = 1E-2
    if (.not. present(maxIter)) maxIter = 100

    ! Set initial values
    successfulStep = .true.
    done = .false.
    X = X0
    lambda = lambda0
    numVars = size(X0)
    iter = 0
    numFuncEvals = 0
    allocate(step(numVars))
    step = 0
    flag = 0

    if (.not. present(maxFuncEvals)) maxFuncEvals = 1000 * numVars

    ! Evaluate initial cost function
    costFuncVec = FUN(X)
    if (any(.not. isfinite(costFuncVec))) stop 'Initial point undefined. Levenberg-Marquardt cannot continue!'
    numResid = size(costFuncVec)
    numFuncEvals = 1

    ! Allocate Jacobian
    allocate(Jacobian(numResid, numVars))

    ! Allocate gradient
    allocate(gradient(numVars))
    allocate(JTJ_diag(numVars))

    ! Allocate lower triangular A
    allocate(A((numVars*(numVars+1))/2))
    allocate(JTJ(size(A)))
    allocate(AdiagElem(numVars))

    ! Compute, which elements of the lower diagonal vector are on the main diagonal.
    do k = 1,numVars
        AdiagElem(k) = k + (2*numVars-k)*(k-1)/2
    end do

    ! Compute Jacobian at initial point.
    call computeJacobian(FUN, X, tolX, numVars, numResid, Jacobian)
    numFuncEvals = numFuncEvals + 2*numVars
    call computeJTJ(Jacobian, numVars, JTJ) ! JTJ = Jacobian^T * Jacobian (lower triangular).
    !$omp parallel do
    do k = 1, numVars
        gradient(k) = dot_product(Jacobian(:,k), costFuncVec)
    end do
    !$omp end parallel do

    firstOrderOpt = maxval(abs(gradient))

        ! Main solver loop.
    do
        ! Print progress.
        write (line, '(I3,3x,I5,2x,ES11.5,2x,ES11.5,2x,ES11.5,2x,ES11.5)') iter, numFuncEvals, sum(costFuncVec**2), &
                             firstOrderOpt, lambda, norm2(step)
        mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(line), 'disp')

        iter = iter + 1

        ! Construct coeficient matrix.
        A = JTJ
        JTJ_diag = JTJ(AdiagElem)
        ! Add lambda*JTJ to A's diagonal elements.
        A(AdiagElem) = JTJ(AdiagElem) * (lambda + 1.0)

        step = gradient ! Input gradient to DPOTRS and output the solution.

        ! Solve the system.
        mexStat = mexCallMATLAB(0, 0, 1, mxCreateString('Performing Cholesky'), 'disp')
        call DPPTRF(lower, numVars, A, ierr) ! Perform Cholesky factorization. A = L*L^T on output.
        if (ierr /= 0) then
            line = 'Cholesky failed! JTJ_diag('
            write(numChar,*) ierr
            line = trim(trim(line)//adjustl(trim(numChar)))//') = '
            write(numChar,*) JTJ_diag(ierr)
            line = trim(line)//adjustl(trim(numChar))
            mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(line), 'disp')
            do k = 1, numVars
                if (JTJ_diag(k) < 1E-14) then
                    write(numChar,*) k
                    line = trim(adjustl(trim(line))//adjustl(trim(numChar)))//','
                end if
            end do
            mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(line), 'disp')
            exit
        end if
        !mexStat = mexCallMATLAB(0, 0, 1, mxCreateString('Solving the system'), 'disp')
        call DPPTRS(lower, numVars, one, A, step, numVars, ierr) ! Solve system using the above factorization.
        !write(line,*) ierr
        !mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(line), 'disp')

        trialX = X - step ! Possible next point.
        where (abs(trialX) < 1E-9)
            trialX = sign(dble(1E-9), trialX)
        end where

        trialCostFunc = FUN(trialX) ! Evaluate function at trial point.
        numFuncEvals = numFuncEvals + 1
        containsInfs = any(.not. isfinite(costFuncVec)) ! Check, whether trial point contains Infs.

        ! Check whether trial point decreased the sum of squares.
        if (sum(trialCostFunc**2) < sum(costFuncVec**2) .and. .not. containsInfs) then
            X = trialX

            call computeJacobian(FUN, X, tolX, numVars, numResid, Jacobian)
            numFuncEvals = numFuncEvals + 2*numVars
            call computeJTJ(Jacobian, numVars, JTJ) ! JTJ = Jacobian^T * Jacobian (lower triangular).
            !$omp parallel do
            do k = 1, numVars
                gradient(k) = dot_product(Jacobian(:,k), trialCostFunc)
            end do
            !$omp end parallel do
            firstOrderOpt = maxval(abs(gradient))

            if (successfulStep) lambda = max(lambda/10, dble(1E-10)) ! Decrease lambda.
            successfulStep = .true.
            flag = converged(costFuncVec, trialCostFunc, tolFun, norm2(X), norm2(step), tolX, &
                      firstOrderOpt, tolOpt, lambda, iter, maxIter, numFuncEvals, maxFuncEvals)
            costFuncVec = trialCostFunc
            if (flag /= 0) exit
        else ! Trial point faired worse than the current point.
            successfulStep = .false.
            lambda = min(lambda*10, dble(1E10)) ! Increase lambda.
            if (norm2(step) < tolX * (sqrt(epsilon(tolX)) + norm2(X))) then
                flag = 2
            else if (iter > maxIter) then
                flag = -1
            else if (numFuncEvals > maxFuncEvals) then
                flag = -2
            end if
            !print *, maxval(gradient), step, X, sum(trialCostFunc**2)
            step = 0
            if (flag /= 0) exit
        end if

    end do ! Main loop.

    ! Print progress.
    write (line, '(I3,3x,I5,2x,ES11.5,2x,ES11.5,2x,ES11.5,2x,ES11.5,a)') iter, numFuncEvals, sum(costFuncVec**2), &
                   firstOrderOpt, lambda, norm2(step), new_line('a')
    mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(line), 'disp')

    if (all(x == X0)) then
        write(line, '(a)') 'Stayed at the initial point. Solver did not converge to desired accuracy.'
    elseif (flag == 1) then
        write(line, '(a)') 'First-order optimality below TolOpt. Solver converged.'
    elseif (flag == 2) then
        write(line, '(a)') 'Chenge in X below TolX. Solver converged.'
    elseif (flag == 3) then
        write(line, '(a)') 'Change in sum of squares below TolFun. Solver converged.'
    elseif (flag == -1) then
        write(line, '(a)') 'Maximum iterations exceeded. Solver did not converge to desired accuracy.'
    elseif (flag == -2) then
        write(line, '(a)') 'Maximum function evaluations exceeded. Solver did not converge to desired accuracy.'
    end if
    mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(line), 'disp')
    mexStat = mexCallMATLAB(0, 0, 1, mxCreateString(''), 'disp')

    mexStat = mexCallMATLAB(0, 0, 1, mxCreateString('Before Jacobian'), 'disp')
    if (present(funVec)) funVec = costFuncVec
    if (present(solution)) solution = X
    if (present(JTWJ)) then
        JTWJ = computeCovarianceInverse(costFuncVec, JTJ, size(Jacobian, 2))
    end if
    if (present(firstOrderOptAtSolution)) firstOrderOptAtSolution = firstOrderOpt
    if (present(JacobianAtSolution)) JacobianAtSolution = Jacobian
    if (present(exitFlag)) exitFlag = flag
    mexStat = mexCallMATLAB(0, 0, 1, mxCreateString('Passed Jacobian'), 'disp')

end subroutine

! Subroutine to compute Jacobian using central finite difference.
subroutine computeJacobian(FUN, X, tolX, numVars, numResid, Jacobian)
    use omp_lib
    implicit none
    procedure(costFunction) :: FUN
    real(kind = 8), intent(in) :: X(:), tolX
    integer(kind = 8), intent(in) :: numVars, numResid
    real(kind = 8), intent(inout) :: Jacobian(:,:)
    real(kind = 8) :: dx(numVars), xForw(numVars), xBackw(numVars)
    real(kind = 8), allocatable :: deriv(:)
    integer :: i
    logical, allocatable :: infInd(:)
    mwPointer :: mxCreateString
    integer(kind = 4) :: mexStat, mexCallMATLAB
    
    mexStat = mexCallMATLAB(0, 0, 1, mxCreateString('Computing Jacobian'), 'disp')

    allocate(deriv(numResid), infInd(numResid))

    dx = max(0.25*TolX*abs(X), dble(1E-10))

    !$omp parallel do private(xForw, xBackw, deriv, infInd)
    do i = 1, numVars
        xForw = X; xBackw = X
        xForw(i) = X(i) + dx(i)
        xBackw(i) = X(i) - dx(i)

        deriv = (FUN(xForw) - FUN(xBackw)) / (2.0 * dx(i))

        infInd = .not. isfinite(deriv)
        if(any(infInd)) then
            block
            integer, allocatable :: indVec(:)
            real(kind = 8) :: vecMean
            indVec = find(.not. infInd)
            vecMean = sum(deriv(indVec)) / size(indVec)
            indVec = find(infInd)
            deriv(indVec) = vecMean
            end block
        end if

        Jacobian(:,i) = deriv
    end do
    !$omp end parallel do
end subroutine

! Subroutine to compute lower diagonal coefficient matrix
subroutine computeJTJ(Jacobian, numVars,JTJ)
    use omp_lib
    implicit none
    real(kind = 8), intent(in) :: Jacobian(:,:)
    integer(kind = 8), intent(in) :: numVars
    real(kind = 8), intent(inout) :: JTJ(:)
    integer :: i, j, ind

    !$omp parallel do private(j, ind)
    do i = 1, numVars
        do j = 1, i
            ind = i + (2*numVars-j)*(j-1)/2
            JTJ(ind) = dot_product(Jacobian(:,i), Jacobian(:,j))
        end do
    end do
    !$omp end parallel do
end subroutine

! Function to check convergence.
integer function converged(costFuncVec, trialCostFunc, tolFun, norm2X, norm2step, tolX, &
                       firstOrderOpt, tolOpt, lambda, iter, maxIter, numFuncEvals, maxFuncEvals)
    implicit none
    real(kind = 8), intent(in) :: costFuncVec(:), trialCostFunc(:), tolFun, norm2X, norm2step, tolX, &
                                  firstOrderOpt, tolOpt, lambda
    integer, intent(in) :: iter, maxIter, numFuncEvals, maxFuncEvals
    real(kind = 8) :: trialSumSq, sumSq

    sumSq = sum(costFuncVec**2)
    trialSumSq = sum(trialCostFunc**2)

    if (firstOrderOpt < tolOpt) then
        converged = 1
    elseif (1.0/lambda > tolX .and. norm2step < tolX * (sqrt(epsilon(tolX)) + norm2X)) then
        converged = 2
    else if (1.0/lambda > tolFun .and. abs(sumSq - trialSumSq) < tolFun * sumSq) then
        converged = 3
    else if (iter >= maxIter) then
        converged = -1
    else if (numFuncEvals >= maxFuncEvals) then
        converged = -2
    else
        converged = 0
    end if

end function

! Function, which returns a logical vector containing true, if the value is finite (not NaN and not +/-infinity).
function isfinite(vec)
    implicit none
    real(kind = 8), intent(in) :: vec(:)
    logical :: isfinite(size(vec))
    real(kind = 8) :: maxReal8
    maxReal8 = huge(vec(1))

    isfinite = (abs(vec) >= 0 .and. abs(vec) <= maxReal8)
end function

! Find true indices.
function find(logicVec) result(trueIndices)
    implicit none
    logical, intent(in) :: logicVec(:)
    integer :: trueCount, i, k
    integer, allocatable :: trueIndices(:)

    trueCount = 0
    do i = 1, size(logicVec)
        if (logicVec(i) .eqv. .true.) trueCount = trueCount + 1
    end do

    allocate(trueIndices(trueCount))

    k = 1
    do i = 1, size(logicVec)
        if (logicVec(i) .eqv. .true.) then
            trueIndices(k) = i
            k = k + 1
        end if
    end do
end function

function computeCovarianceInverse(costFuncVec, JTJ, numParams) result(JTWJ)
    implicit none
    real(kind = 8), intent(inout) :: costFuncVec(:), JTJ(:)
    integer, intent(in) :: numParams
    real(kind = 8), allocatable :: JTWJ(:,:)
    real(kind = 8) :: fitVariance
    integer(kind = 8) :: numData, i, j
    
    numData = size(costFuncVec)
    numParams = size(Jacobian,2)

    allocate(JTWJ(numParams, numParams))

    fitVariance = sum(costFuncVec**2) / (numData - numParams + 1)
    JTJ = JTJ / fitVariance

    do i = 1, numParams
        do j = 1, i
            ind = i + (2*numParams-j)*(j-1)/2
            JTWJ(i,j) = JTJ(ind)
        end do
    end do
    
    do j = 2, numParams
        do i = 1, j-1
            JTWJ(i,j) = JTWJ(j,i)
        end do
    end do

end function

end module lmSolverModule
