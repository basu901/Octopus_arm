import numpy as np
import show_trunk
import scipy.linalg as linalg
import copy,math

#Used for zero un-stretched length
def linear_elephant(N,a1,c1,a2,c2,b,mass1,mass2,L):
    n=N*2*2+N
    grid=np.zeros(shape=(n,n))
    val = np.zeros(n)

    #populating x1
    k=0
    m=0 #maintains row number
    while(m<n):
        l=0 
        j=0 
        while(j<n):
            if m==j:
                grid[m][j]=-(a1[k]+c1[k]+c2[k])
                grid[m][j+2]=c2[k]
                grid[m][j+4]=c1[k]
                if j+5<n:
                    grid[m][j+5]=a1[k+1]                
                break
            else:
                grid[m][j]=-(c1[k]+c2[k])
                grid[m][j+2]=c2[k]
                grid[m][j+4]=c1[k]
            j=j+5
        k=k+1
        m=m+5
                    
    #populating y1
    m=1
    k=0
    while m<n:
        j=1
        while j<n:
            if m==j:
                grid[m][j]=a1[k]+c1[k]+c2[k]
                grid[m][j+2]=-c2[k]
                if j+5<n:
                    grid[m][j+5]=-a1[k+1]
                break
            else:
                count = 0
                while(count<k):
                    grid[m][j]=c1[k]+c2[k]
                    grid[m][j+2]=-c2[k]
                    count = count+1
                    j=j+5
            
        m = m+5
        k=k+1


    #populating x2
    m = 2
    k = 0
    while m<n:
        j=0
        while j<n:
            if m==j:
                grid[m][j]=-(c2[k]+a2[k])
                if j+5<n:
                    grid[m][j+5]=a2[k+1]
                break
            else:
                count = 0
                while(count<k):
                    grid[m][j]=c2[k]
                    grid[m][j+2]=-c2[k]
                    count = count+1
                    j = j+5
                grid[m][j]=c2[k]
                j = j+2

        m = m+5
        k = k+1

    #populating y2
    m = 3
    k = 0
    while m<n:
        j=1
        while j<n:
            if m==j:
                grid[m][j]=c2[k]+a2[k]
                if j+5<n:
                    grid[m][j+5]=-a2[k+1]
                break
            else:
                count = 0
                while(count<k):
                    grid[m][j] = -c2[k]
                    grid[m][j+2]=c2[k]
                    count = count+1
                    j = j+5
                grid[m][j]=-c2[k]
                j = j+2

        m = m+5
        k = k+1

    #populating s  

    m=4
    k =0
    
    while(m<n):
        j=0
        while(j<n):
            if m==j:
                grid[m][j]=-(b[k]+c1[k])
                if j+5<n:
                    grid[m][j+5]=b[k+1]
                break
            else:
                count = 0
                while(count<k):
                    grid[m][j]=c1[k]
                    grid[m][j+4]=-c1[k]
                    j = j+5
                    count = count+1
                grid[m][j]=c1[k]
                j = j+4
                
        k=k+1
        m=m+5


    #Modifying last row of grid
    
    j=4
    while j<n:
        grid[-1][j]=grid[-1][j]-b[-1]
        j=j+5


    #Populating constant terms    
    m1 = 1 #Index of mass 1
    m2 = 3 #Index of mass 2
    k1 = 0
    k2 = 0
    while m1<n and m2<n:
        val[m1] = mass1[k1]*9.8
        val[m2] = mass2[k2]*9.8
        k1 = k1+1
        k2 = k2+1
        m1 = m1+5
        m2 = m2+5
    
    val[-1]=-1*b[-1]*L

    #print("A:")
    #print(grid)

    #print("b:")
    #print(val)

    LU = linalg.lu_factor(grid)
    
    x = linalg.lu_solve(LU,val)

    return x


#Used for non-zero unstretched length
def stretched(N,a,b,c,mass,L,u_a,u_b,u_c):
    n=N*3
    grid=np.zeros(shape=(n,n))
    val = np.zeros(n)
    k=0
    m=0
    while(m<n):
        l=0 #Used to change sign of gamma
        j=0 #Keeping track of the column index
        while(j<n):
            if m==j:
                grid[m][j]=-(a[k]+c[k])
                grid[m][j+2]=c[k]
                if j+3<n:
                    grid[m][j+3]=a[k+1]                
                break
            else:
                grid[m][j]=c[k]*((-1)**(l+1))
                if l==0:
                    l=1
                else:
                    j=j-1
                    l=0
            j=j+2
        k=k+1
        m=m+3
                    
    m=1
    k=0
    while m<n:
        j=1
        while j<n:
            if m==j:
                grid[m][j]=a[k]+c[k]
                if j+3<n:
                    grid[m][j+3]=a[k+1]*-1
                break
            else:
                grid[m][j]=c[k]
                j=j+2
            j=j+1
        m = m+3
        k=k+1

    #print("On adding y terms")  
    #print(grid)

    k=0
    m=2
    while(m<n):
        l=0
        j=0
        while(j<n):
            if m==j:
                grid[m][j]=-(b[k]+c[k])
                if j+3<n:
                    grid[m][j+3]=b[k+1]
                break
            else:
                grid[m][j]=c[k]*((-1)**l)
                if l==0:
                    l=1
                else:
                    j=j-1
                    l=0
                    #if j==m:
                    #    continue
            j=j+2
        k=k+1
        m=m+3

    m=n-1
    j=2
    while j<n:
        grid[m][j]=grid[m][j]-b[N]
        j=j+3
        
    #print('On adding s terms')  
    #print(grid)

    #z = len(u_a)
    #grid=np.zeros(shape=(n,n))
    d=np.zeros(len(u_a))
    #print(u_a)
    #print(u_c)
    #print(len(d))
    i = 0
    while i<len(u_a):
        for j in range(0,i):
            d[i] += u_a[j]-u_c[j]
        #print(i)
        d[i] += (u_a[i]-u_c[i])
        i = i+1
        
    m = 0
    k = 0
    while m<n:
        if k+1<len(u_a):
            val[m]= a[k+1]*u_a[k+1]-(a[k]*u_a[k])-c[k]*d[k]
        else:
            val[m] = -(a[k]*u_a[k])-c[k]*d[k]
        m += 3
        k = k+1


    e = np.zeros(len(u_b))

    for i in range(0,len(u_b)):
        for j in range(0,i):
            e[i] += u_b[j]
        e[i] += u_b[i]

    
    m = 1
    k = 0
    while m<n:
        if k+1<len(u_a):
            val[m] = u_b[k]*a[k]+u_b[k+1]*a[k+1]+c[k]*e[k]
        else:
            val[m] = u_b[k]*a[k]+c[k]*e[k]
        m += 3
        k = k+1


    f = np.zeros(len(u_c))

    for i in range(0,len(u_c)):
        for j in range(0,i):
            f[i] += u_c[j]
        f[i] += u_c[i]

    m =2
    k= 0
    while m<n:
        if k+1<len(u_c):
            val[m] = u_c[k+1]*b[k+1]-u_c[k]*b[k]+c[k]*f[k]
        else:
            val[m] = -u_c[k]*b[k]+c[k]*f[k]
        m += 3
        k = k+1
    
    m=1
    k=0
    while m<n:
        val[m]=mass[k]*9.8
        k=k+1
        m=m+3
    
    val[n-1]=-1*b[N]*L

    print("A:")
    print(grid)

    print("b:")
    print(val)
    
    g=np.linalg.solve(grid,val)

    return g


def check_positions(X,N):
    x = copy.deepcopy(X)
    x = x.reshape(N,5)
    x = x[:,2:4].T

    a = x[0]<0
    b = x[1]<0
    if np.dot(a,b)==True:
        return 4

    a = x[0]<0
    b = x[1]>0
    if np.dot(a,b)==True:
        return 2

    a = x[0]>0
    b = x[1]<0
    if np.dot(a,b)==True:
        return 3

    a = x[0]>0
    b = x[1]>0
    if np.dot(a,b)==True:
        return 1

    return 0



def store_xy(X,n):
    i=2
    j=3
    k=4
    x = list()
    y=list()
    s=list()
    while i<n:
        x.append(X[i])
        y.append(X[j])
        s.append(X[k])
        i = i+5
        j = j+5
        k = k+5

    for i in range(1,len(x)):
        x[i] = x[i]+x[i-1]
        y[i] = y[i]+y[i-1]
        s[i] = s[i]+s[i-1]

    return x[-1],y[-1]

L =20
N = 6
mass1 = np.ones(N)
mass2 = np.ones(N)

x1_val = list()
y1_val = list()

x2_val = list()
y2_val = list()

x3_val = list()
y3_val = list()

x4_val = list()
y4_val = list()

loops = 200



#0-199:a1 = [0.1,5] a2=[0.1,20] c1 = [0.1,5] c2 = [0.1,20]


for i in range(loops):
    a1 = np.random.uniform(low=0,high = 1, size = N)
    a1 = 0.1*a1+(1-a1)*5
    #a1= np.array([0.5,0.5,0.5])
    a2 = np.random.uniform(low=0,high = 1, size = N)
    a2 = 0.1*a2 + (1-a2)*5
    #a2 = np.array([0.5,0.5,0.5])
    b= np.random.uniform(low =0, high =1, size = N+1)
    b = 0.1*b + (1-b)*20
    #b = np.array([10,10,10,0.1])
    c1 = np.random.uniform(low = 0, high=1, size = N)
    c1 = 5*c1 + (1-c1)*20
    #c1 = np.array([5,5,5])
    c2 = np.random.uniform(low=0,high = 1, size = N)
    c2 = 0.1*c2 + (1-c2)*5
    #c2 = np.array([0.5,0.5,0.5])
    X = linear_elephant(N,a1,c1,a2,c2,b,mass1,mass2,L)

    show_trunk.show_two_layer(X,len(X),L,i)

    color = check_positions(X,N)
    if color == 1:
        x_val,y_val = store_xy(X,len(X))
        x1_val.append(x_val)
        y1_val.append(y_val)
    if color == 2:
        x_val,y_val = store_xy(X,len(X))
        x2_val.append(x_val)
        y2_val.append(y_val)
    if color == 3:
        x_val,y_val = store_xy(X,len(X))
        x3_val.append(x_val)
        y3_val.append(y_val)
    if color == 4:
        x_val,y_val = store_xy(X,len(X))
        x4_val.append(x_val)
        y4_val.append(y_val)
        
    del X,a1,a2,b,c1,c2
#print('X:',X)

#show_trunk.show_end_point(x1_val,y1_val,x2_val,y2_val,x3_val,y3_val,x4_val,y4_val,L)

