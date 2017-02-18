!! HW1-- ECE697NA -  E. Polizzi
!! This program proposes to solve the Schrodinger eigenvalue equation in 1D
!!                 -alpha.Psi''+U(x).Psi(x)=E.Psi(x)  along [0,L]
!! with alpha=(hbar**2)/(2.me) and Dirichlet boundary conditions:  Psi(0)=0;  Psi(L)=0
!!
!! This simple problem can be solved analytically if U(x)=0, so we can get a reference solution
!! This is the particle in a box problem (infinite qauntum well) 
!! 
!! We will use the physics SI units: L=1d-9 (1 nm box), hbar= 1.0545718d-34 , me=9.10938356d-31
!! The energy we will be obtained in Joule, we will consider dividing E by q=1.602d-19 when postprocessing to get in eV
!!
!! We propose to use the Finite Difference Method (FDM) to discretize the problem
!! considering a segment of [0,L] and 1 to N points of discretization (uniform grid)
!! We propose to solve the resulting linear system both directly and iteratively,
!!            look at accuracy, timings, etc. 
!!
!!
!! The code can be compiled as follows: "ifort -o hw1 hw1.f90 -mkl" using intel Fortran compiler and MKL
!!                               or   : "gfortran -o hw1 hw1.f90 ....." using Gnu Fortran compiler, requires link to lapack and blas
!!
!! You can run the code as follows: "./hw1"
!!
!! To do: complete the call to the LAPACK routine for solving the eigenvalue

program hw1

!!!!!!! variable declaration
  implicit none
  double precision,parameter :: pi=3.1415926535897932d0,hbar= 1.0545718d-34,me=9.10938356d-31,q=1.602d-19
  double precision :: L,alpha,h,t 
  integer :: i,j,N 
  double precision,dimension(:),allocatable :: x,E,Ea,r,normr
  double precision,dimension(:,:),allocatable :: psia,psi,A
  integer :: t1,t2,tim
  !! for LAPACK
  integer :: info,lwork
  double precision,dimension(:),allocatable :: work
!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Constant !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  L=1d-9
  alpha=hbar**2/(2.0d0*me)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! READ the values of N  !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"Enter number of grid points N:"
  read *,N


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Setup uniform grid !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! grid step
  h=L/(N-1)

  !! find the coordinates x of each node
  allocate(x(1:N))
  do i=1,N !all grid points
     x(i)=(i-1)*h
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Obtain Analytical Solution!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(psia(1:N,1:N)) ! eigenvectors
  do j=1,N ! loop over the eigenvalues (N is the max number)
     do i=1,N !all grid points   (at each grid point)   !!!!!!!
        psia(i,j)=sqrt(2.0d0/L)*sin(j*pi*x(i)/L)
     enddo
  enddo


  allocate(Ea(1:N)) ! eigenvalues
  do j=1,N
     Ea(j)=((j*pi*hbar)**2)/(2.0d0*me*L**2)
  enddo

  print *,'Analytical Soltuion: '
  do j=1,min(10,N)
    print *,j,Ea(j)/q
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Construct Ax=b !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! memory allocations
  allocate(A(N,N)) ! matrix (dense storage)


  !! Assemble Matrix A (Laplacian)
  A=0.0d0  ! all elements zero
  t=(hbar**2)/(2.0d0*me*h**2)
  do i=1,N 
     A(i,i)=2.0d0*t ! diagonal elements
     if (i/=N) then
        A(i,i+1)=-1.0d0*t ! superdiagonal elements
        A(i+1,i)=-1.0d0*t ! subdiagonal elements
     endif
  enddo


  !! Include Boundary conditions at the edges
  A(1,1:N)=0.0d0
  A(1:N,1)=0.0d0
  A(1,1)=1.0d0 ! diagonal element

  A(N,1:N)=0.0d0
  A(1:N,N)=0.0d0
  A(N,N)=1.0d0 ! diagonal element

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Solve the eigenvalue problem A.psi=E.psi !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(psi(N,N)) ! eigenvector calculated
  allocate(E(1:N))   ! eigenvalue calculated

!!!!  using a Direct solver (and Lapack routine DSYEV)!!!
  lwork=3*N-1 ! for lapack
  allocate(work(1:lwork)) ! for lapack
  !! we make a copy of A  (since it is overwritten calling LAPACK by the eigenvectors)
  psi=A ! copy of A 

  call system_clock(t1,tim)
  print *,'start LAPACK procedure' 


!!!!!!!!!!>>>>>>>>>>>>>>>> to complete (input: matrix A (copy into psi), output: psi contains the eigenvector and E the eigenvalues
  
  !character :: JOBZ, UPLO
  !integer :: LDA
  !JOBZ = 'V'
  !UPLO = 'U'
  !LDA = N
  call dsyev('V', 'U', N, A, N, E, work, lwork, info)
  lwork = work(1)
  deallocate(work)
  if (info==0) then
     print *,'LAPACK success'
  else
     print *,'LAPACK error',info
  end if


  call system_clock(t2,tim)

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Print Total time for solving Ax=b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  print *,'Total time',(t2-t1)*1.0d0/tim


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Check relative error of the first ten eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print *,'Relative Error on eigenvalues:'
  do j=1,min(10,N)
     print *,j,E(j)/q,abs(E(j)-Ea(j))/abs(Ea(j))
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Check Residual of first 10 eigenpairs !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'Eigenpairs residuals:'
  allocate(r(N))
  !! perform r=A*psi_i-E(i)*psi_i
  ! Remark: We can use an intrinsic Fortran function matmul
  !         however this is not going to provide the best performances as
  !         we will discuss later in class  
  !r=matmul(A,psi)-e*psi
  do j=1,min(10,N) ! loop over all eigenvalues
     r=psi(:,j)
     call DGEMM('N','N', N, 1, N, 1.0d0, A, N, psi(1,j), N, -1.0d0*E(j), r, N ) !! BLAS function
     ! calculate the norm of r  (we use Fortran intrinsic function)
     print *,j,maxval(abs(r))/maxval(abs(E(1)*psi(:,1)))!we take the norm relative (for example)
  enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Save results !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! save the first N-1 energie solutions (in eV) 
  open(10,file='E.out',status='replace')
  do j=1,N-2
     write(10,*) j,Ea(j)/q,E(j)/q
  end do
  close(10)



  !! save the two first eigenvectors (probability density)
  open(10,file='Psi12.out',status='replace')
  do i=1,N
     write(10,*) x(i),abs(psi(i,1))**2,abs(psi(i,2))**2
  end do
  close(10)




end program hw1
