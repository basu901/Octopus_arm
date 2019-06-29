import numpy as np


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

    print("A:")
    print(grid)

    print("b:")
    print(val)
    
    g=np.linalg.solve(grid,val)

    return g


def find_xys_given_abg(n, alpha_arr, beta_arr, gamma_arr,m_arr,L):
        A = np.zeros((3*n, 3*n))
        b = np.zeros((3*n, 1))
        
        # Build A and b
        for k in range(0, n):
            x_k_idx = (k*3)
            y_k_idx = (k*3) + 1
            s_k_idx = (k*3) + 2
            b[(k*3)] = 0
            b[(k*3) + 1] = m_arr[k]*9.8
            if (k != n-1):
                b[(k*3) + 2] = 0
            else:
                b[(k*3) + 2] = -1* beta_arr[k+1] * L          
            # M - X Forces
            A[(k*3), x_k_idx] = -1*(alpha_arr[k] + gamma_arr[k])
            if (k != n-1):
                A[(k*3), (k+1)*3] = alpha_arr[k+1]
            for i in range(0, k):
                A[(k*3), i*3] = -1*gamma_arr[k]
            for i in range(0, k+1):
                A[(k*3), (i*3) + 2] = gamma_arr[k]
                
            # M - Y Forces
            A[(k*3) + 1, y_k_idx] = (alpha_arr[k] + gamma_arr[k])
            if (k != n-1):
                A[(k*3) + 1, ((k+1)*3) + 1] = -1*alpha_arr[k+1]
            for i in range(0, k):
                A[(k*3) + 1, (i*3) + 1] = gamma_arr[k]  
                
            # L -X Forces
            A[(k*3) + 2, s_k_idx] = -1*(beta_arr[k] + gamma_arr[k])
            if (k != n-1):
                A[(k*3) + 2, ((k+1)*3) + 2] = beta_arr[k+1]
            else:
                for i in range(0, k+1):
                    A[(k*3) + 2, (i*3) + 2] += -1*beta_arr[k+1]
            for i in range(0, k):
                A[(k*3) + 2, (i*3) + 2] += -1*gamma_arr[k] 
            for i in range(0, k+1):
                A[(k*3) + 2, (i*3)] = gamma_arr[k]
                
        # Solve system
        #print A
        return np.linalg.solve(A, b)




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


