import numpy as np
from scipy.linalg import null_space

#x_c holds the x co-ordinate of the masses
#y_c holds the y co-ordinate of the masses

N = int(input("Enter the number of masses:"))
n = N*3

print("Enter the masses in vector form:")
m = np.array(list(map(float,input().split())))

grid = np.zeros(shape=(n,n+1))

b = np.zeros(n)
spring_c = np.zeros(n)
print("Enter the x co-ordinates of the masses:")
x_c= np.array(list((map(float,input().split()))))

#print(x_c.shape)
#print(x_c.size)
#print("The x co-ordinates:",x_c)

print("Enter the y co-ordinates of the masses:")
y_c = np.array(list(map(float,input().split())))
#print("The y co-ordinates:",y_c)

print("Enter the S displacement vector:")
s_c = np.array(list(map(float,input().split())))
#print("The S displacements:",s_c)

print("Enter the total length of the rod:")
L = float(input())

i = 0
j = 0
k = 0
g = 9.8

while i<n:
    grid[i][j]=-1*x_c[k]
    l = 0
    while l <= k:
        grid[i][j+2] += s_c[l]-x_c[l]
        l += 1
    if j<n-3:
        grid[i][j+3] = x_c[k+1]
    j = j+3
    k = k+1
    i += 3



j = 0
k = 0
i = 1
while i< n:
    grid[i][j] = y_c[k]
    l = 0
    while l <=k :
        grid[i][j+2] += y_c[l]
        l += 1
    if j<n-3:
        grid[i][j+3] = -y_c[k+1]
    j = j+3
    k = k+1
    i += 3



j = 0
k =0
i = 2
while i<n:
    grid[i][j+1] = -s_c[k]
    l = 0
    while l <= k:
        grid[i][j+2] += -s_c[l]+x_c[l]
        l += 1
    if j<n-3:
        grid[i][j+4] = s_c[k+1]
    j = j+3
    k = k+1
    i += 3

t = 0
for i in range(len(s_c)):
    t += s_c[i]

if t>=L:
    print("S values exceeded the total length of the rod")
    exit()
    
grid[-1][-1]= L - t


i = 1
j = 0
while i<n:
    b[i] = g*m[j]
    i += 3
    j += 1

print("A: ")
print(grid,"\n")

ns = null_space(grid)
ns = ns*np.sign(ns[0,0])
print("Null space of A:",ns,"\n")

print("b:",b,"\n")

        
spring_c = np.linalg.lstsq(grid,b)

error = spring_c[3]
spring_c = spring_c[0]

print("Least square error:")
print(error,"\n")
print("Spring stiffness vector for the 3n springs:")
s = spring_c[:-1]
s = s.reshape(N,3)
print(s)

print("Spring stiffness of (n+1)th beta spring:")
print(spring_c[-1])
