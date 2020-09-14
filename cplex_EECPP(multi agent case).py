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
#Data
#------------------------------------------------------------------------------
agentNumber = 10
colNumber = 3
rowNumber = 3
Battery_capacity_constraint = 11
nodesNumber = colNumber * rowNumber
departurePoint = 0
distance_lambda = 0.1164
turn_gamma = 0.0173
turn_gamma = 0.015
radians_to_degrees = 180/(math.pi)
distanceValue=0
theta_radians=0
theta_degrees=0
agents = [(a) for a in range(agentNumber)]
Nodes = range(1, nodesNumber) # Except departurePoint( arrival point is same as departurePoint)
NodesAndDeparturePoint = [n for n in range(nodesNumber)]
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


q={(i,j,k):0 for i,j,k in arcs}
for i,j,k in arcs:
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

model.parameters.timelimit=150
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
