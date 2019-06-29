import numpy as np
import show_trunk
import copy,math
import scipy.linalg as linalg

#Used for zero un-stretched length
def linear_elephant(N,a,b,c,mass,L):
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
    
    m=1
    k=0
    while m<n:
        val[m]=mass[k]*9.8
        k=k+1
        m=m+3
    
    val[n-1]=-1*b[N]*L

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

    #print("A:")
    #print(grid)

    #print("b:")
    #print(val)

    LU = linalg.lu_factor(grid)
    
    x = linalg.lu_solve(LU,val)

    return x

def check_positions(X,N):
    x = copy.deepcopy(X)
    x = x.reshape(N,3)
    x = x[:,0:2].T

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



    
    

def classifier(X,N):
    x = copy.deepcopy(X)
    x = x.reshape(N,3)
    x = x[:,0:2]
    class_num = 0
    for i in range(0,N-1):
        print (math.atan2(x[i][1],x[i][0])-math.atan2(x[i+1][1],x[i+1][0]))
        if abs(math.atan2(x[i][1],x[i][0])-math.atan2(x[i+1][1],x[i+1][0]))>=(math.pi/4):
            class_num += 1

    del x
    print('\n')
    return class_num

def store_xy(X,n):
    i=0
    j=1
    k=2
    x = list()
    y=list()
    s=list()
    while i<n:
        x.append(X[i])
        y.append(X[j])
        s.append(X[k])
        i = i+3
        j = j+3
        k = k+3

    x_last = x[-1]
    y_last = y[-1]

    for i in range(1,len(x)):
        x[i] = x[i]+x[i-1]
        y[i] = y[i]+y[i-1]
        s[i] = s[i]+s[i-1]

    return x[-1],y[-1],x_last,y_last
    

loops = 10000

#number of masses
N = 10
mass = np.ones(N)

L =100

u_a = np.array([0.1,0.1,0.5,0.6,1.0,0.3,0.8,0.5,0.2,1.0])
u_b = np.array([0.5,0.1,0.5,0.1,1.0,0.9,0.45,0.5,0.1,0.2,0.1])
u_c = np.array([2,1,0.5,0.2,1,1.5,0.8,0.2,0.1,0.5])

x1_val = list()
y1_val = list()

x2_val = list()
y2_val = list()

x3_val = list()
y3_val = list()

x4_val = list()
y4_val = list()

x_info = list()
y_info = list()


for i in range(loops):
    a = np.random.uniform(low=0.1,high = 1, size = N)
    a = a*0.1+(1-a)*2
    b= np.random.uniform(low =0.1, high =1, size =N+1)
    b = b*0.1+(1-b)*2
    c = np.random.uniform(low = 0.1, high=1, size =N)
    c = c*0.1+(1-c)*2
    #X = linear_elephant(N,a,b,c,mass,L)
    X = stretched(N,a,b,c,mass,L,u_a,u_b,u_c)

    #show_trunk.show_trunk_unstretched(X,len(X),L,i)
  
    color = check_positions(X,N)
    if color == 1:
        x_val,y_val,x_last,y_last = store_xy(X,len(X))
        x1_val.append(x_val)
        y1_val.append(y_val)
    if color == 2:
        x_val,y_val,x_last,y_last = store_xy(X,len(X))
        x2_val.append(x_val)
        y2_val.append(y_val)
    if color == 3:
        x_val,y_val,x_last,y_last = store_xy(X,len(X))
        x3_val.append(x_val)
        y3_val.append(y_val)
    if color == 4:
        x_val,y_val,x_last,y_last = store_xy(X,len(X))
        x4_val.append(x_val)
        y4_val.append(y_val)
        
    
    x_total,y_total,x_last,y_last = store_xy(X,len(X))

    #print(x_total,y_total,x_last,y_last)

    x_info.append(x_last)
    x_info.append(x_total)

    y_info.append(y_last)
    y_info.append(y_total)

    del X,a,b,c


show_trunk.show_grid_end_effector_orientation(x_info,y_info,L,loops)
show_trunk.show_end_point(x1_val,y1_val,x2_val,y2_val,x3_val,y3_val,x4_val,y4_val,L)
