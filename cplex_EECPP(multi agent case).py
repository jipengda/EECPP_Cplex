# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 18:50:35 2020

@author: 冀鹏
"""
from docplex.mp.model import Model
import Data
import numpy as np
import math
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
#supplementary functions
#------------------------------------------------------------------------------
def check_obstacle(obstacles, lastNode, newNode, colNumber, rowNumber):
    sideLength = 1
    flag=1 # initialize flag=1
    for obstacle in obstacles:
        x3 = obstacle % colNumber
        y3 = obstacle // colNumber
        r_xmin = x3 - sideLength/2
        r_xmax = x3 + sideLength/2
        r_ymin = y3 - sideLength/2
        r_ymax = y3 + sideLength/2
        if RSIntersection(r_xmin, r_xmax, r_ymin, r_ymax, lastNode, newNode, colNumber, rowNumber) == 1:
            flag = 0 #flag is variable to replace directly return False of True because we need to test all obstacles
        else:
            flag = 1
        if flag == 0:#once flag change to 0, means one obstacle is hit, thus new node can not be added
            return False
    return True# This means all obstacles are not hit, new node can be considered to be added
"""
 * if [x1, x2] and [x3, x4] (x4 maybe smaller than x3) has interrection or has not，if yes: return 1， if no: return 0
"""


def IntervalOverlap(x1, x2, x3, x4):
    t = 0    
    if x3 > x4:
       t = x3
       x3 = x4	   
       x4 = t
    if x3 >= x2 or x4 <= x1:
        return 0
    else:
        return 1
    
"""
 * judge rectangular r and line segment AB has intersection or not，if yes: return 1，if no: return 0
"""
def RSIntersection(r_xmin,r_xmax,r_ymin, r_ymax, nodeA, nodeB, colNumber, rowNumber):
    A_x = nodeA % colNumber
    A_y = nodeA // colNumber
    B_x = nodeB % colNumber
    B_y = nodeB // colNumber
    if (A_y == B_y):# line segement is parallel to x axis//线段平行于x轴
        if A_y <= r_ymax and A_y >= r_ymin:
            return IntervalOverlap(r_xmin, r_xmax, A_x,B_x)
        else:
            return 0

	# Echange point A and point B to let point B has bigger y coordinate//AB两点交换，让B点的y坐标最大

    # Exchange node A and node B, let B's y value is bigger
    t = 0
    if A_y > B_y:
       t = A_y
       A_y = B_y
       B_y = t
       t= A_x
       A_x = B_x
       B_x=t
	
    # In line segment//xianduan AB, to find point C and D
    # Two points secure a line: (x-x1)/(x2-x1)=(y-y1)/(y2-y1)
    k = (B_x - A_x)/(B_y - A_y)
    if A_y < r_ymin:
       D_y = r_ymin
       D_x = k*(D_y - A_y) + A_x
    else:
       D_y=A_y
       D_x=A_x
    if B_y > r_ymax:
       C_y = r_ymax
       C_x = k*(C_y-A_y) + A_x
    else:
       C_y = B_y
       C_x = B_x
    if C_y >= D_y: # y axis has overlap
       return IntervalOverlap(r_xmin, r_xmax,D_x, C_x)
    else:
       return 0

#------------------------------------------------------------------------------
#Data
#------------------------------------------------------------------------------

#add obstalces part
agentNumber = 10
colNumber = 4
rowNumber = 4
Battery_capacity_constraint = 11
nodesNumber = colNumber * rowNumber
departurePoint = 0
distance_lambda = 0.1164
turn_gamma = 0.0173
obstacles=[12,13]
Nodes=[i for i in range(nodesNumber) if i not in obstacles and i!= departurePoint]
NodesAndDeparturePoint = Nodes + [departurePoint]
AllNodes = NodesAndDeparturePoint + obstacles
radians_to_degrees = 180/(math.pi)
distanceValue=0
theta_radians=0
theta_degrees=0
agents = [(a) for a in range(agentNumber)]
#Nodes = range(1, nodesNumber) # Except departurePoint( arrival point is same as departurePoint)
agents_NodesAndDeparturePoint = [(a,n) for a in agents for n in NodesAndDeparturePoint]
edges = [(i,j) for i in NodesAndDeparturePoint for j in NodesAndDeparturePoint ]
agents_edges = [(a,i,j) for a in agents for i,j in edges]
arcs = [(i,j,k) for i in NodesAndDeparturePoint for j in NodesAndDeparturePoint for k in NodesAndDeparturePoint]
agents_arcs = [(a,i,j,k) for a in agents for i,j,k in arcs]

coord_x = Data.create_coord_x(colNumber, rowNumber)
coord_y = Data.create_coord_y(colNumber, rowNumber)


D=Data.create_D(nodesNumber, coord_x, coord_y)
c = {(i,j):0 for i,j in edges}
for i,j in edges:
    distance_cost = D[i][j] * distance_lambda
    c[(i,j)] = distance_cost
    
for o,p in edges:
    View = check_obstacle(obstacles, o, p, colNumber, rowNumber)
    if View == 0:
        c[(o,p)] = math.inf
    else:
        pass
    
seq=[-10,-9,-8,-7,-6,-5,-4,-3-2,-1,0,1,2,3,4,5,6,7,8,9,10]
fixed_turn_gamma=0.0173
turn_factor=0.0001
random.seed(10) 
q={(i,j,k):0 for i,j,k in arcs}
for i,j,k in arcs:
#   turn_gamma = fixed_turn_gamma + random.choice(seq) * turn_factor
    theta_radians=math.pi-np.arccos(round((D[i][j]**2+D[j][k]**2-D[i][k]**2)/(2*D[i][j]*D[j][k]),2))
#    if math.isnan(theta_radians) == True:
#        theta_radians=0
    theta_degrees=theta_radians * radians_to_degrees
    turning_cost = theta_degrees * turn_gamma
    q[(i,j,k)] = turning_cost

# An arc flow model for the basic EECPP
model = Model("EECPP PROBLEM")
x = model.binary_var_dict(agents_edges, name = 'X') # flow variables, 1 if the agent goes directly from node i to node j
I = model.binary_var_dict(agents_arcs, name = "I") # use I in order to linearly replace x (i,j) x x (j,k)
d = model.continuous_var_dict(agents_NodesAndDeparturePoint, name = "D",lb=0, ub=nodesNumber) # d is a dummy variable associated with node i for subtour elimination

Z = model.sum(( c[(i,j)] * x[(a,i,j)] ) for a,i,j in agents_edges) + model.sum((q[(i,j,k)]*I[(a,i,j,k)]) for a,i,j,k in agents_arcs)
real_Z = Z - model.sum( (q[(i, 0, k)] * I[(a, i, 0, k)]) for a, i,k in agents_edges) # when agents comes back to depot
# it does not need to make further turning for the second node!
model.minimize(real_Z) #(1)


for j in Nodes:
    model.add_constraint(model.sum( x[(a,i,j)] for a,i in agents_NodesAndDeparturePoint if i!=j ) == 1, ctname = None) #(2)

model.add_constraint(model.sum( x[(a, i, 0)] for a,i in agents_NodesAndDeparturePoint if i!= departurePoint ) <= agentNumber, ctname = None) #(3)

#once a drone visits a node, it also departs from the same node. Each node can and can only be visited only once.
for a in agents:
    for j in NodesAndDeparturePoint:
        model.add_constraint( model.sum( x[(a,i,j)] for i in NodesAndDeparturePoint if i!=j )==model.sum( x[(a,j,k)] for k in NodesAndDeparturePoint if k!=j ), ctname=None ) #(4)
    
for a in agents:
    for i,j,k in arcs:
        model.add_constraint( I[(a,i,j,k)] >= x[(a,i,j)] + x[(a,j,k)] - 1, ctname = None ) #(5)
        model.add_constraint( I[(a,i,j,k)] <= x[(a,i,j)], ctname = None ) #(6)
        model.add_constraint( I[(a,i,j,k)] <= x[(a,j,k)], ctname = None ) #(7)

#(8)
for a in agents:
    u= model.sum(( c[(i,j)] * x[(a,i,j)] ) for i,j in edges) + model.sum(( q[(i,j,k)]*I[(a,i,j,k)] ) for i,j,k in arcs if j!=0)
    expr = model.sum((q[(i,0,k)] * I[(a,i,0,k)]) for i,k in edges)
    model.add_constraint( u - expr <= Battery_capacity_constraint, ctname = 'agent %d BCC Limit'%a)  

#(9) subtour elimination constraint
for a in agents:
    for i,j in edges:    
        if j!=departurePoint:
            model.add_indicator( x[(a,i,j)],d[a,i] + 1 == d[a,j],name=None)

model.parameters.timelimit=1800
#model.parameters.timelimit=3600

solution = model.solve(log_output=True)

m=coord_x
n=coord_y

draw_edges =[i for i in agents_edges if x[i].solution_value>0.9]
width = 2*colNumber - 2
height = 2*rowNumber - 2
plt.figure(figsize=(width, height)) #(width, height) or (colNumber-1, rowNumber-1)
plt.xlabel("Coordinate X")
plt.ylabel("Coordinate Y")
plt.title("Solution of Cplex for Multi-agent EECPP Problem")

for n in NodesAndDeparturePoint:
    if n!=departurePoint:
        plt.scatter(x=coord_x[n],y=coord_y[n],color='blue')
    else:
        plt.scatter(x=coord_x[n],y=coord_y[n],color='green')
    
for n in range(len(coord_x)):# name of node
    plt.annotate(str(n),xy=(coord_x[n],coord_y[n]),xytext=(coord_x[n]+0.01,coord_y[n]+0.1),color='red')

for a,i,j in draw_edges:
    if a==0:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='blue')
    if a==1:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='orange')
    if a==2:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='green')    
    if a==3:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='red')
    if a==4:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='purple')
    if a==5:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='brown')
    if a==6:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='pink')
    if a==7:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='gray')
    if a==8:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='olive')
    if a==9:
        plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='cyan')

#plt.grid(True)
#plt.savefig("size 4 x 3.png")
#plt.savefig("size 4 x 3.eps")
plt.show()
