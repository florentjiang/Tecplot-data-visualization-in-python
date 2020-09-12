import re
import math
import pdb
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# load flow field data file
def readAirfoil(path):
    if not os.path.isfile(path):
        print('file does not exist')
        return
    file_in = open(path,'r')
    x = []
    y = []
    for line in file_in:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))
    return x, y

# load flow field data file
def loadData(path, z_slice):
    eps = 1E-5
    if not os.path.isfile(path):
        print('file does not exist')
        return
    file_in=open(path,'r')
    solution = []
    for line in file_in:
        if re.match("  [ -]0",line) and len(re.split('\s+|,',line)) >= 3 \
            and float(re.split('\s+|,',line)[3]) >= z_slice - eps \
            and float(re.split('\s+|,',line)[3]) <= z_slice + eps:
            solution.append([float(i) for i in re.split('\s+|,',line)[1:10]])
    df_solution = pd.DataFrame(solution, columns=['x','y','z','r','ru','rv','rw','re','p'])
    df_solution.drop_duplicates(inplace=True) # data points could be duplicated on the internal block boundaries, removal needed 
    df_solution['velocity'] = np.sqrt((df_solution['ru']/df_solution['r'])**2 + (df_solution['rv']/df_solution['r'])**2)
    return df_solution[['x','y','velocity']]

# Calculation of geometric properties of boundary element segments
def geometry(x_list,y_list,seg_list):
    Ns = np.sum(seg_list) # total no. of segments
    Np = Ns+1 # total no. of segment end-points
    
    # total no. of segments at the beginning of each boundary element
    seg_num = [0 for _ in range(seg_list.size)]
    for i in range(1,seg_list.size):
        seg_num[i] = seg_num[i-1] + seg_list[i-1]

    x, y = [np.zeros(Np) for i in range(2)]
    x[0] = x[-1] = x_list[0]; y[0] = y[-1] = y_list[0]
    for i in range(seg_list.size):
        x[seg_num[i]:seg_num[i]+seg_list[i]+1] = np.linspace(x_list[i],x_list[i+1],seg_list[i]+1)
        y[seg_num[i]:seg_num[i]+seg_list[i]+1] = np.linspace(y_list[i],y_list[i+1],seg_list[i]+1)

    return x, y

# Function used to decide if a data point is inside the airfoil
# Use Ray casting algorithm to determine if a point is inside a polygon
def isInsidePolygon(p, poly):
    ans = False
    i = -1
    l = len(poly)
    j = l - 1
    while i < l - 1:
        i += 1
        if ((poly[i][0] <= p[0] and p[0] < poly[j][0]) or (poly[j][0] <= p[0] and p[0] < poly[i][0])):
            if (p[1] < (poly[j][1] - poly[i][1]) * (p[0] - poly[i][0]) / (poly[j][0] - poly[i][0]) + poly[i][1]):
                ans = not ans
        j = i
    return ans

# The data extraction pipeline
def pipeline(alpha, scale, z_slice, Nx, Ny, data_dir, airfoil_path):
    # define rotation matrix
    D2R = np.pi/180
    rot = np.array([np.cos(alpha*D2R), -np.sin(alpha*D2R), np.sin(alpha*D2R), np.cos(alpha*D2R)]).reshape(2,2) # rotation matrix

    # load coordinates of external boundary points
    x_list1 = np.array([-10.,40.,40.,-10.,-10.])
    y_list1 = np.array([-12.,-12.,12,12,-12])
    seg_list1 = np.array([40,20,40,20])
    x1, y1 = geometry(x_list1,y_list1,seg_list1)


    # load coordinates of airfoil (clockwise / last pt = first pt)
    x_list2, y_list2 = readAirfoil(airfoil_path)
    x_list2 = scale*np.array(x_list2)[::-1] # clockwise
    y_list2 = scale*np.array(y_list2)[::-1] # clockwise
    x_list2 = np.dot(np.c_[x_list2,y_list2],rot)[:,0] # rotate airfoil according to alpha
    y_list2 = np.dot(np.c_[x_list2,y_list2],rot)[:,1]
    Ns2 = x_list2.size - 1
    seg_list2 = np.array([1 for _ in range(Ns2)])
    x2, y2 = geometry(x_list2,y_list2,seg_list2)

    # Combine the internal & external boundaries
    x = np.append(x1[:-1],x2[:-1])
    y = np.append(y1[:-1],y2[:-1])
    X = np.linspace(x.min(),x.max(),Nx)
    Y = np.linspace(y.min(),y.max(),Ny)
    X,Y = np.meshgrid(X[1:-1],Y[1:-1])
    X = X.ravel(); Y = Y.ravel()

    # Operation for each time step:
    G = []
    for flow_file in os.listdir(data_dir):
        # load flow field solution data and correct data points inside the airfoil region
        data_path = os.path.join(data_dir,flow_file)
        df_focus = loadData(data_path, z_slice)
        I = []
        for i in range(X.size):
            if isInsidePolygon([X[i],Y[i]], np.c_[x_list2,y_list2]):
                df_focus = df_focus.append({'x': X[i], 'y': Y[i], 'velocity':0.0}, ignore_index=True)
                I.append(i)

        coord = df_focus[['x','y']].values
        df_focus[['x','y']] = np.dot(coord,rot)

        # Interpolate on flow field solution data on visualization data points
        xi = np.linspace(x.min(),x.max(),Nx)
        yi = np.linspace(y.min(),y.max(),Ny)
        zi = griddata((scale*df_focus['x'], scale*df_focus['y']), df_focus['velocity'], (xi[None,:], yi[:,None]), method='linear')
        G.append(zi)
        return G, xi, yi, x2, y2