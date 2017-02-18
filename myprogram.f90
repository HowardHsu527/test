program nameofprogram
implicit none

integer :: i, j, k
real :: x, y, z

x = 3.61
y = cos(x)
z = x + y

i = 3
j = i**2
k = i - j

open(10, file = 'mydata.dat')
write(10,*) i, j, k, x, y, z
close(10)
end program nameofprogram