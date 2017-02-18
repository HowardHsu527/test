program LearnMat
	implicit none
	double precision, dimension(:,:), allocatable :: A, L, U
	integer :: i, N, M

	print *, "Enter number of grid points N:"
	read *, N
	print *, "Enter number of grid points M:"
	read *, M

	!Assemble A matrix
	allocate(A(N,N))
	allocate(L(N,N))
	allocate(U(N,N))

	A = 0.0d0 !all elements zero
	do i = 1, N
		A(i,i) = 1.0d0
		if (i/=N) then
			A(i,i+1) = 2.0d0
			A(i+1,i) = 2.0d0
		end if
	end do

	do i = 1, M
		if (i/=M) then
			A(i,i+M-1) = 5
			A(i+M-1,i) = 5
		end if
	end do

	!print *, A
	
	!!Generate Lower and Upper triangle Matrix
	L = 0 ! all elements zeros
	U = 0
	do i = 1, N
		L(i,:i) = A(i,:i)
		U(i,i+1:) = A(i,i+1:)
	end do

	print *, L
	!print *, u
end program LearnMat