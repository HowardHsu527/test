!! HW1-- ECE697NA -  E. Polizzi
!! This program proposes to solve the Poisson equation in 2D
!!                 -Lapacian(V)=f(x,y)
!! with Dirichlet boundary conditions:  V=0 at the boundaries.
!! We will also be using  f(x,y)=(alpha^2+beta^2)*sin(alpha*x)*sin(beta*y), 
!!                 where  alpha=2.0d0*pi/Lx;  beta=2.0d0*pi/Ly                                    
!! And Lx=Ly=L=10.0d0
!!
!! This simple problem can be solved analytically, so we can get a reference solution
!! 
!! We propose to use the Finite Difference Method (FDM) to discretize the problem
!! considering a rectangular domain defined by [0,Lx] and [0,Ly] with respectively Nx and Ny points
!!                                                                 of discretization (uniform grid)
!! We propose to solve the resulting linear system both directly and iteratively,
!!            look at accuracy, timings, etc. 
!!
!!
!! The code can be compiled as follows: "ifort -o hw1_2d hw1_2d.f90 -mkl" using intel Fortran compiler and MKL
!!
!! You can run the code as follows:  "./hw1_2d direct" or "./hw1_2d iter"


program hw1_2d

!!!!!!! variable declaration
  implicit none
  integer :: i,j,N,it,itmax,k,Nxy,e
  integer :: t1,t2,tim
  double precision :: h,L,alpha,beta,pi,normr,normb
  double precision, dimension(:,:), allocatable :: A,cA
  double precision,dimension(:),allocatable :: b,v,x,r,xi,yi,xold
  character(len=100) :: name
  !! for LAPACK
  integer :: info
  integer,dimension(:),allocatable :: IPIV
!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! read command line argument name="iter" or name="direct"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,name)

  if ((trim(name)/="iter").and.(trim(name)/="direct")) then
     print *,"command line argument missing"
     stop
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Constant !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  L=10.0d0
  pi=3.1415926535897932d0
  alpha=2.0d0*pi/L
  beta=alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! READ the values of N  !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"Enter number of grid points along x or y (Nx=Ny):"
  read *,Nxy

  N=(Nxy)**2 ! total number of unknowns

  if (N>6400) then
   print *,'Nx too large'
   stop
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Setup uniform grid!!! !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! grid step
  h=L/(Nxy-1)

  !! find the coordinates x and y of each nodes
  allocate(xi(N))
  allocate(yi(N))
  k=0 
  do i=1,Nxy
     do j=1,Nxy
        k=(i-1)*Nxy+j ! global numbering of the node
        xi(k)=(i-1)*h
        yi(k)=(j-1)*h
     enddo
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Obtain Analytical Solution!!!!! 
!!!!!!   (at each grid point)   !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(v(1:N))

!<<<<<<<<<<<<<<<<<<<<< To complete (Analytical)
!  

!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Construct Ax=b !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! memory allocations
  allocate(A(N,N)) ! matrix (dense storage)
  allocate(b(N))   ! right-hand-side
  allocate(x(N))   ! numerical solution


  !! Assemble Matrix A (Laplacian)
  A=0.0d0  ! all elements zero


!!!<<<<<<<<<<<<<<<<<< To complete A matrix (Hint: you can loop over blocks first - similarly to "include BC" below)
!
  do i = 1, N
    A(i,i) = 4
    if (i/=N) then
      A(i,i+1) = -1
      A(i+1,i) = -1
    end if
  end do

  do i = 1, N
    if (i/=(N-Nxy)) then
      A(i,i+Nxy) = -1
      A(i+Nxy,i) = -1
    end if
  end do
!
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !! Assemble right-hand side RHS 
  b=0.0d0
  do i=1,N
     b(i)=(alpha**2+beta**2)*sin(alpha*xi(i))*sin(beta*yi(i))
  enddo


  !! Include Boundary conditions at the edges
  do i=1,Nxy
     do j=1,Nxy
        k=(i-1)*Nxy+j ! global numbering of the node
        if ((i==1).or.(i==Nxy).or.(j==1).or.(j==Nxy)) then
           A(k,1:N)=0.0d0
           A(1:N,k)=0.0d0
           A(k,k)=1.0d0 ! diagonal element
           b(k)=0.0d0 ! B.C. for RHS
        end if
     end do
  end do

  print *, A

end program hw1_2d