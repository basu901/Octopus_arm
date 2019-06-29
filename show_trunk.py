import matplotlib.pyplot as plt
import numpy as np
import math as m

def show_trunk_unstretched(X,n,l,pic_index):
    i = 0
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

    for i in range(1,len(x)):
        x[i] = x[i]+x[i-1]
        y[i] = y[i]+y[i-1]
        s[i] = s[i]+s[i-1]


    #x1 , y1 = [0,x[0]],[0,(-1)*y[0]]
    x1 , y1 = [0,x[0]],[0,y[0]]
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(x1,y1,marker = 'o')
    ax.plot(x,y,color = 'b',marker = 'o')


    x2 , y2 = [0,s[0]],[0,0]
    plt.plot(x2,y2,marker = 'x')
    y2 = np.zeros(len(s))
    ax.plot(s,y2,color = 'r',marker = 'x')
    ax.plot([s[-1],l],[0,0],color = 'r',marker = 'x')

    for i in range(len(s)):
        x1,y1 = [x[i],s[i]],[y[i],0]
        ax.plot(x1,y1,color='g')
    
    plt.gca().invert_yaxis()
    #plt.show()
    fig.savefig('two_mass_one_layer\iter_2_500\plot'+str(pic_index)+'.png')
    plt.close(fig)


def show_two_layer(X,n,l,pic_index):
    i = 0
    j = 1
    k = 4
    x1 = list() #storing x pos of 1st line of masses
    y1=list()   #storing y pos of 1st line of masses
    s=list()    #storing s pos of masses
    x2 = list() #storing x pos of 2nd line of masses
    y2 = list() #storing y pos of 2nd line of masses
    while i<n:
        x1.append(X[i])
        y1.append(X[j])
        x2.append(X[i+2])
        y2.append(X[j+2])
        s.append(X[k])
        i = i+5
        j = j+5
        k = k+5

    for i in range(1,len(x1)):
        x1[i] = x1[i]+x1[i-1]
        y1[i] = y1[i]+y1[i-1]
        x2[i] = x2[i]+x2[i-1]
        y2[i] = y2[i]+y2[i-1]
        s[i] = s[i]+s[i-1]


    x , y = [0,x1[0]],[0,y1[0]]
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(x,y,color = 'r',marker = 'o')
    ax.plot(x1,y1,color = 'r',marker = 'o')

    x , y = [0,x2[0]],[0,y2[0]]
    ax.plot(x,y,color = 'b',marker = 'o')
    ax.plot(x2,y2,color = 'b',marker = 'o')

    for i in range(len(x1)):
        x,y = [x2[i],x1[i]],[y2[i],y1[i]]
        ax.plot(x,y,color = 'g')


    x2 , y2 = [0,s[0]],[0,0]
    plt.plot(x2,y2,color = 'k',marker = 'x')
    y2 = np.zeros(len(s))
    ax.plot(s,y2,color = 'k',marker = 'x')

    for i in range(len(s)):
        x,y = [x1[i],s[i]],[y1[i],0]
        ax.plot(x,y,color='g')
        
    ax.plot(l,0,color = 'k',marker = 'x')
    plt.gca().invert_yaxis()
    #plt.show()
    fig.savefig('two_layer_structs_rectified\plot20'+str(pic_index)+'.png')

def show_end_point(x1_val,y1_val,x2_val,y2_val,x3_val,y3_val,x4_val,y4_val,l):
    
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(x1_val,y1_val,color = 'r',marker = 'o')
    ax.scatter(x2_val,y2_val,color = 'b',marker = 'o')
    ax.scatter(x3_val,y3_val,color = 'g',marker = 'o')
    ax.scatter(x4_val,y4_val,color = 'y',marker = 'o')
    ax.plot([0,l],[0,0],color ='k')
    
    #x2 = np.zeros(n)
    #y2 = np.zeros(len(s))
   
    #ax.plot(s,y2,color = 'r',marker = 'x')
    #ax.plot(l,0,color = 'r',marker = 'x')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.show()
    #fig.savefig('classifier_'+str(class_num)+'\plot_1'+str(pic_index)+'.png')


def show_grid_end_effector_orientation(X,Y,L,n):
    x = np.asarray(X)
    #print(x)
    x = x.reshape(n,2) #contains the position of the end effector in the 2nd column and the last x(delta x) in the first column
    y = np.asarray(Y)
    y = y.reshape(n,2) #contains the position of the end effector in the 2nd column and the last x(delta x) in the first column

    x_points = [m.floor(i)for i in x[:,1]]
    y_points = [m.floor(i)for i in y[:,1]]

    x_plot = np.asarray(x_points)
    y_plot = np.asarray(y_points)

    x_min = m.floor(min(x_plot))-1
    x_max = m.ceil(max(x_plot))+1

    y_min = m.floor(min(y_plot))-1
    y_max = m.ceil(max(y_plot))+2

    if x_max<L:
        x_max = L+1

    x_plot = x_plot + 0.5   #x value of the centroid of the grid cells
    y_plot = y_plot + 0.5   #y value of the centroid of the grid cells

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(x_plot,y_plot,color = 'k' ,marker = 'o')

    x_d = np.zeros(n) #contains the delta x values
    x_d = x[:,0]       

    y_d = np.zeros(n)  #contains the delta y values
    y_d = y[:,0]

    x_add = np.zeros(n)
    y_add = np.zeros(n)

    plt.plot([0,L],[0,0],color = 'k')
    
    for i in range(n):
        angle = m.atan2(y_d[i],x_d[i])
        #print(angle)
        x_add = 0.5* m.cos(angle)
        y_add = 0.5* m.sin(angle)
        #print(x_plot)
        #print(str(type(x_plot)))
        x_second = x_plot[i]-x_add
        y_second = y_plot[i]-y_add
        
        plt.plot([x_second,x_plot[i]],[y_second,y_plot[i]])

    plt.gca().invert_yaxis()
    plt.yticks(np.arange(y_min,y_max))
    plt.xticks(np.arange(x_min,x_max))
    plt.grid()
    plt.show()
                        
