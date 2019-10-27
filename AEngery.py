# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 15:34:34 2016

@author: Lina492375qW1188 and mnb11338
"""
import numpy as np
import math
from scipy import spatial

ATOM_NUMBER = 92482

def ReadAngle(filename):
    f = open('.\%s' %filename, 'r')
    adata = []
    for line in f:
        if 'ITEM: ATOMS id c_1 c_2 ' in line: 
            for line in f:
                adata.append(tuple([float(cell) for cell in line.split()]))
    f.close()
    adata = sorted(adata, key=lambda adata: adata[0]) # sort=排序 data
    return adata

def ReadUnrolled(filename):
    f = open('.\%s' %filename, 'r')
    unrolled=[]
    for line in f:
        if 'Atoms' in line:
            for i, line in enumerate(f): #enumerate = 照順序列出來
                unrolled.append(tuple([float(cell) for cell in line.split()]))
                if i > ATOM_NUMBER:
                    break
    f.close()

    del unrolled[0]
    return unrolled

def flatten(unrolled, data_S1, data_S, xmin, xmax):
    # unrolled: row1[3] = x coordinate
    #           row1[4] = y coordinate
    # data_S1:  row2[1] = eangle/atom at step 1
    #           row2[2] = ebond/atom at step 1
    # data_S:   row3[1] = eangle/atom at step S
    # data_S:   row3[2] = ebond/atom at step S
    separation = ([tuple([row1[3], row1[4], row2[1], row2[2], row3[1], row3[2]])
                     for (row1, row2, row3) in zip(unrolled, data_S1, data_S) 
                     if row1[3] > xmin and row1[3] < xmax])
   #假設separation=[(1,2,3),(4,5,6),(7,8,9)]，x=[1,4,7] 
    x = [row[0] for row in separation]
    y = [row[1] for row in separation]
    angle = [row[4]-row[2] for row in separation]
    bond = [row[5]-row[3] for row in separation]

    flatten_coord = list(zip(x, y, angle, bond))
    return flatten_coord

def OnlyCreases(angle_emap, bond_emap, threshold, edge, edge_threshold):
    
    L_a = len(angle_emap)
    
    for index, row_a, row_b in zip(range(L_a), angle_emap, bond_emap):
        
        if row_a[2] < threshold:
            angle_emap[index] = (row_a[0], row_a[1], 0.0)
            bond_emap[index] = (row_b[0], row_b[1], 0.0)
            
        elif abs(row_a[0]) > edge and row_a[2] < edge_threshold:
                angle_emap[index] = (row_a[0], row_a[1], 0.0)
                bond_emap[index] = (row_b[0], row_b[1], 0.0)
                
    return angle_emap, bond_emap

def RemoveRegion(angle_emap, bond_emap, xmin, xmax, ymin, ymax):
    
    L_a = len(angle_emap)
    
    for index, row_a, row_b in zip(range(L_a), angle_emap, bond_emap):
        
        if (row_a[0] < xmax and row_a[0] > xmin and 
            row_a[1] < ymax and row_a[1] > ymin):
            angle_emap[index] = (row_a[0], row_a[1], 0.0)
            bond_emap[index] = (row_b[0], row_b[1], 0.0)
            
    return angle_emap, bond_emap

# Data Input
output_S1 = ReadAngle('pe0.R80.5_Kb10000_Ka2000_step0') # 未擠壓
output_S = ReadAngle('pe0.R80.5_Kb10000_Ka2000_step7000') # 已擠壓，每次都要改！
unrolled = ReadUnrolled('project92482')         #84,85 扭過-未扭=扭的結果
#Kb20000 - project 91756 // Kb10000 - project 93958
                       
# Set parameters
Lside = 80 # left side
Rside = -80 # right side
darkness = 1.8
brightness = 0.2
magnification = 1.0 # 幾乎不用改
Nx = 70     #80
Ny = 70     #80

# Seperate angle and bond energy map.
x, y, angle, bond = zip(*flatten(unrolled, output_S1, output_S, Rside, Lside))
angle_tot = list(zip(x, y, angle))
bond_tot = list(zip(x, y, bond))

import PlotCollection
intensity = PlotCollection.Intensity2D(darkness, brightness, magnification)

# Plot total angle and bond energy map.
#intensity.ScipyBinned(angle_tot, Nx, Ny) 
#intensity.ScipyBinned(bond_tot , Nx, Ny) 

#eangle_tot = sum([row[2] for row in flat_angle_tot]) # calculate total angle energy.
#ebond_tot = sum([row[2] for row in flat_bond_tot])   # calculate total bond energy.



# Creases Angle and Bond Energy:
# OnlyCreases would cut those useless part first.
# useless part: including those ones with too high energy or where you don't want.
# Remember that you need this command to create angle_crease and bond_crease in
# order to use all the following functions.

angle_crease, bond_crease = OnlyCreases(angle_tot, bond_tot, 0, 100, 0)

#intensity.ScipyBinned(angle_crease, Nx, Ny) 
#intensity.ScipyBinned(bond_crease, Nx, Ny)

#eangle_crease = sum([row[2] for row in angle_crease])
#ebond_crease = sum([row[2] for row in bond_crease])

# return the data of Energy map
#bs = intensity.EnergyMap(angle_crease, Nx, Ny)


###############################################################################
# separate x, y, angle energy and bond energy from angle_crease and bond_crease.
# 
x, y, a_crease = zip(*angle_crease)
x, y, b_crease = zip(*bond_crease)
ab = [a+b for a, b in zip(a_crease, b_crease)]
ab_crease = list(zip(x, y, ab))
bs = intensity.EnergyMap(ab_crease, Nx, Ny)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

bz = bs.statistic

norm = (cm.colors.Normalize(vmin=np.nanmin(bz)*darkness,
                            vmax=np.nanmax(bz)*brightness))

fig = plt.figure()

(plt.imshow(bz, extent=(np.amin(bs.y_edge), np.amax(bs.y_edge),
                                np.amin(bs.x_edge), np.amax(bs.x_edge)),
                                interpolation='sinc', origin='lower',
                                aspect='auto', cmap=cm.afmhot, norm=norm))
#, origin='lower'
d = (np.nanmax(bz)-np.nanmin(bz))/4
(plt.colorbar(ticks=[np.nanmin(bz), np.nanmin(bz)+d, np.nanmin(bz)+2*d, 
                     np.nanmin(bz)+3*d, np.nanmax(bz)]))

plt.xlabel('y')
plt.ylabel('x')

# FIND_CLOSEST_POINT
def Find(x, arr):
    """Find the point in arr which is closest to value."""
    ref = np.array([x]*len(arr))
    d = [(index, abs(diff)) for index, diff in enumerate(ref-arr)]    
    # sort array d and make the smallest difference be the first one. And then
    # return the index of it.
    i_smallest = sorted(d, key=lambda b:b[1])[0][0]
    return i_smallest

# DEFINE_POLYNOMIAL
def PolyFunc(press_arr):
    """press_arr should contain at least 2 elements. This function will return
       the polynomial function that connect all points in press_arr."""
    polyfunc_list = []
    # Reset polyfunc_list everytime for each call, which is used to avoid 
    # ambiguity and make each call of PolyFunc return list of polynomial 
    # functions that corresponds to input array press_arr.
    
    for i in range(0, len(press_arr)-1):
        x1, y1 = press_arr[i][2], press_arr[i][3]
        x2, y2 = press_arr[i+1][2], press_arr[i+1][3]
               
        if y1 != y2:
            m = (x2-x1)/(y2-y1) # slop
            b = x1 - m*y1       # intercept
            polyfunc_list.append(np.poly1d([m, b]))
        else:
            b = y1
            polyfunc_list.append(np.poly1d([b]))
        
    return polyfunc_list

def pos_on_pline(ends, polyfunc):
    #ends = [(x_i, y_i, xdata_i, ydata_i), 
    #        (x_{i+1}, y_{i+1}, xdata_{i+1}, ydata_{i+1})]
# Pixel:
    #x_begin, y_begin = ends[0][0], ends[0][1]
    #x_end, y_end = ends[1][0], ends[1][1]

    xdata_begin, ydata_begin = ends[0][2], ends[0][3]
    xdata_end, ydata_end = ends[1][2], ends[1][3]

# how many value in bs.x_edge between two point.
    n = abs(Find(ydata_end, bs.x_edge)-Find(ydata_begin, bs.x_edge))
# how many pixel between two points (b.t.p.).
    #n = int(abs(y_end - y_begin)) 
    dx = (xdata_end - xdata_begin)/n # the data's real difference b.t.p.
    dy = (ydata_end - ydata_begin)/n

    pos_list = [] # position list stores all the positions of point on line.
    for i in range(n):
        x_pos = polyfunc(ydata_begin + i*dy)  # xpos is belong to y_edge
        y_pos = ydata_begin + i*dy            # ypos is belong to x_edge
        #z_pos = bz[Find(y_pos, bs.x_edge)][Find(x_pos, bs.y_edge)]
        
        x_real = bs.x_edge[Find(y_pos, bs.x_edge)]
        y_real = bs.y_edge[Find(x_pos, bs.y_edge)]
        # x_real is in x_edge and y_real is in y_edge.
        pos_list.append((x_real, y_real))
        
    return pos_list


### GAUSSIAN_SUM
from scipy.optimize import curve_fit

def Gaussian(x, a, x0, w):
    e0 = -(x-x0)**2/(2*w**2)
    return a*np.exp(e0)

def Area(xc, yc):
    """xc is in x_edge and yc is in y_edge."""
    t = 0.5 # threshold of yc.

    # Y-fit:
    Y = bz[Find(xc, bs.x_edge)]
    dY = np.gradient(np.log10(Y))     # 1st derivative
    ddY = np.gradient(dY)             # 2nd derivative
    k = abs(ddY[Find(yc, bs.y_edge)]) # curvature at point nearest to yc
    w = np.sqrt(1/(2*k))              # the estimate of width
    
    lower = [0, yc-t, 0]
    higher = [np.inf, yc+t, w]
    
    height = bz[Find(xc, bs.x_edge)]
    
    ydiff = bs.y_edge[1]-bs.y_edge[0]
    horiz = np.delete(bs.y_edge, -1) + ydiff
    
    para = curve_fit(Gaussian, horiz, height, bounds=(lower, higher))
    a, x0, w = para[0]

    # X-fit:
    X = bz[:, Find(yc, bs.y_edge)]
    dX = np.gradient(np.log10(X))
    ddX = np.gradient(dX)
    kX = abs(ddX[Find(xc, bs.x_edge)])
    Xw = np.sqrt(1/(2*kX))
    
    Xlower = [a-t, xc-t, 0]
    Xhigher = [a+t, xc+t, Xw]
    
    Xheight = bz[:, Find(yc, bs.y_edge)]
    
    xdiff = bs.x_edge[1]-bs.x_edge[0]
    Xhoriz = np.delete(bs.x_edge, -1) + xdiff
    
    para = curve_fit(Gaussian, Xhoriz, Xheight, bounds=(Xlower, Xhigher))
    a, Xx0, Xw = para[0]

    #return np.sqrt(2*np.pi)*a*w
    return np.pi*Xw*w*a

### INTERACTIVE_PRESS
L = []
area = []
press_arr = []

def on_press(event, polyfunc = []):
    
    press = event.x, event.y, event.xdata, event.ydata
    press_arr.append(press)

    try:
        polyfunc = PolyFunc(press_arr)
        ends = [press_arr[-2], press_arr[-1]]
#press1:[event.x][event.y] ...        = press_arr[-2][0]press_arr[-2][1] ..
#press2:... [event.xdata][event.ydata]= .. press_arr[-1][2]press_arr[-1][3]

# Calculate length: *math.cos(press_arr[-1][3]*math.pi/160)
        Li = (np.sqrt((press_arr[-1][2]*math.cos(press_arr[-1][3]*math.pi/160)-press_arr[-2][2]*math.cos(press_arr[-2][3]*math.pi/160))**2+
                     ((press_arr[-1][3]-press_arr[-2][3])**2)/2.5))
        COS = ((80**2+80**2-Li**2)/(2*80*80))                #餘弦定理
        Lr = (2*math.pi*80)*(math.acos(COS)/(2*math.pi))         #球面長度
        L.append(Lr)

# Extract and draw points:
        pos_list = pos_on_pline(ends, polyfunc[-1])
        pos_x, pos_y = zip(*pos_list)
        plt.plot(pos_y, pos_x, marker = 'p', markersize = 6, color = 'b')
        plt.xlim((np.amin(bs.y_edge), np.amax(bs.y_edge)))
        plt.ylim((np.amin(bs.x_edge), np.amax(bs.x_edge)))
        fig.canvas.draw()

# Calculate area:        
        area.append(sum([Area(xc, yc) for xc, yc in pos_list]))

    except:
        print("Only 1 point in press_arr!")

def on_key(event): # 空白鍵會執行兩次... Don't ask me why!!!

    print("=================================================")
    print("The length of all line segment:", sum(L))
    print("The sum of all line segment:", sum(area))
    print("=================================================")

    f = open('.\data.txt','a')       
    f.write(str(sum(L))+' '+str(sum(area))+'\n')
    f.close
    
    del L[:]         # Reset L
    del area[:]      # Reset area
    del press_arr[:] # Reset press_arr


cidpress = fig.canvas.mpl_connect('button_press_event', on_press)
cidkey = fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()
