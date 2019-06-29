import numpy as np
from linear_elephant import *
from show_trunk import show_trunk
import sys

prompt = 'Enter the number of masses:'
N = int(input(prompt))
prompt = 'Enter the mass values:'
mass = np.array(list(map(float,input(prompt).split())))
prompt = 'Enter the values of alpha:'
a = np.array(list(map(float,input(prompt).split())))
prompt = 'Enter the values of beta:'
b = np.array(list(map(float,input(prompt).split())))
prompt = 'Enter the values of gamma:'
c = np.array(list(map(float,input(prompt).split())))
prompt = 'Enter the total length of the rod:'
L = int(input(prompt))

if len(sys.argv)>1:
    if int(sys.argv[1])==1:
        u_a = list(map(float,input("Unstretched length of alpha's:").split()))
        u_b = list(map(float,input("Unstretched length of beta's:").split()))
        u_c = list(map(float,input("Unstretched length of gamma's:").split()))


    if int(sys.argv[1])==2:
        un_a = float(input("Enter unstretched length of all alpha's"))
        un_b = float(input("Enter unstretched length of all beta's"))
        un_c = float(input("Enter unstretched length of all beta's"))


    X = stretched(N,a,b,c,mass,L,u_a,u_b,u_c)        

else:
    X = linear_elephant(N,a,b,c,mass,L)
    #X = find_xys_given_abg(N, a, b, c,mass,L)

print("Mass Co-ordinates")
s = X.reshape(N,3)
print(s)
    
show_trunk(X,len(X),L)
