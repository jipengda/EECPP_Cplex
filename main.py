# calculate total_distance and total_degrees for the test map
# optimal_path replace gen_model_label
# based on / refer ReproduceJalil_cplex.py
from docplex.mp.model import Model
import Data
#import random
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches

colNumber = 8
rowNumber = 15
sideLength=1 * 10
nodesNumber = colNumber * rowNumber # size 4 x 5
departurePoint = 0
arrivalPoint = departurePoint

agentNumber = 1
A = agentNumber
control_value = 1000 # 1000>>127.7, big enough to cover all the nodes.
Battery_capacity_constraint = control_value
C = Battery_capacity_constraint
lamda = 0.1164
garma = 0.0173
radians_to_degrees = 180/(math.pi)
distanceValue=0
theta_radians=0
theta_degrees=0
#Nodes = range(1, nodesNumber) # Except departurePoint( arrival point is same as departurePoint)
Nodes = [node for node in range(1, nodesNumber)]
NodesAndDeparturePoint =[node for node in range(0, nodesNumber)]
obstacles = [36,37,
             44,45,
             52,53,
             60,61,
             68,69,
             74,75,76,77,
             82,83,84,
             99,100,
             107,108]
NodesAndDeparturePoint = list(set(NodesAndDeparturePoint) - set(obstacles))
Nodes = list(set(Nodes) - set(obstacles))
AllNodes = [node for node in range(0, nodesNumber)]
edges = [(i,j) for i in NodesAndDeparturePoint for j in NodesAndDeparturePoint]
arcs = [(i,j,k) for i in NodesAndDeparturePoint for j in NodesAndDeparturePoint for k in NodesAndDeparturePoint]


coord_x = Data.create_coord_x(colNumber, rowNumber)
coord_y = Data.create_coord_y(colNumber, rowNumber)
distance={(i,j):0 for i,j in edges}
degrees={(i,j,k):0 for i,j,k in arcs}
c = {(i,j):0 for i,j in edges}
for i,j in edges:
    # when i, j are not adjacent, distance is infinity.
    if abs(coord_x[i] - coord_x[j]) > 10 or abs(coord_y[i] - coord_y[j]) > 10: # wrong, or right judgement of condition?
        distanceValue = math.inf
    else:
        distanceValue = np.hypot(coord_x[i] - coord_x[j],coord_y[i] - coord_y[j])

    # or we do it in the opposite way
#    if abs(coord_x[i] - coord_x[j]) <=1 and abs(coord_y[i] - coord_y[j]) <=1 :
#    	distanceValue = np.hypot(coord_x[i] - coord_x[j], coord_y[i]-coord_y[j])
#    else:
#    	distanceValue = math.inf
    distance[(i,j)]=distanceValue
    distance_cost = lamda * distanceValue 
    c[(i,j)] = distance_cost



q={(i,j,k):0 for i,j,k in arcs}
for i,j,k in arcs:
    theta_radians=math.pi-np.arccos(round((distance[(i,j)]**2+distance[(j,k)]**2-distance[(i,k)]**2)/(2*distance[(i,j)]*distance[(j,k)]),2))
    if math.isnan(theta_radians) is True:
        theta_radians = 0
    else:
        pass
    theta_degrees=theta_radians * radians_to_degrees
    degrees[(i,j,k)]=theta_degrees
    turning_cost = garma * theta_degrees
    q[(i,j,k)] = turning_cost

# An arc flow model for the basic EECPP
model = Model("EECPP PROBLEM")
x = model.binary_var_dict(edges, name = 'X') # flow variables, 1 if the agent goes directly from node i to node j
I = model.binary_var_dict(arcs, name = "I") # use I in order to linearly replace x (i,j) x x (j,k)
d = model.continuous_var_list(AllNodes, name = "D") # d is a dummy variable associated with node i for subtour elimination

Z = model.sum((c[(i,j)] * x[(i,j)]) for i,j in edges) + model.sum((q[(i,j,k)]*I[(i,j,k)]) for i,j,k in arcs)
real_Z = Z - model.sum( (q[(i,0,k)] * I[(i,0,k)]) for i,k in edges) # when agent comes back to depot, it does not need to make further turning
model.minimize(real_Z) #(1)

model.add_constraint(model.sum( x[(0,j)] for j in Nodes) == 1, ctname = None)

for j in Nodes:
    model.add_constraint(model.sum( x[(i,j)] for i in NodesAndDeparturePoint ) == 1, ctname = None) #(2)

#model.add_constraint(model.sum( x[(a, i, 0)] for a,i in agents_NodesAndDeparturePoint ) <= A, ctname = None) #(3)

#once a drone visits a node, it also departs from the same node. Each node can and can only be visited only once.
for j in NodesAndDeparturePoint:
    model.add_constraint(model.sum( x[(i,j)] for i in NodesAndDeparturePoint)==1, ctname = None)
    model.add_constraint(model.sum( x[(j,k)] for k in NodesAndDeparturePoint)==1, ctname = None) #(4) 
    
for i,j,k in arcs:
    model.add_constraint( I[(i,j,k)] >= x[(i,j)] + x[(j,k)] - 1, ctname = None ) #(5)
    model.add_constraint( I[(i,j,k)] <= x[(i,j)], ctname = None ) #(6)
    model.add_constraint( I[(i,j,k)] <= x[(j,k)], ctname = None ) #(7)

#(8)
model.add_constraint( (model.sum(( x[(i,j)] * c[(i,j)] ) for i,j in edges) + model.sum(( q[(i,j,k)]*I[(i,j,k)] ) for i,j,k in arcs) )<= C, ctname = None)  

#(9) subtour elimination constraint
for i,j in edges:    
  if j!=departurePoint:
    model.add_indicator( x[(i,j)],d[i] + 1 == d[j],name=None)

solution = model.solve(log_output=True)

m=coord_x
n=coord_y

s=[]
for n in range(len(coord_x)):
    s_temp = []
    s_temp.append("%.1f"%coord_x[n])
    s_temp.append("%.1f"%coord_y[n])
    s.append(s_temp)
    
draw_edges =[i for i in edges if x[i].solution_value>0.9]

#plt.figure(figsize=(2*colNumber-2, 2*rowNumber-2)) #(width, height) or (colNumber-1, rowNumber-1)
plt.figure()
plt.xlabel("Coordinate X")
plt.ylabel("Coordinate Y")
plt.title("Solution EECPP PROBLEM")

        
for n in NodesAndDeparturePoint:
    if n!=departurePoint:
        plt.scatter(x=coord_x[n],y=coord_y[n],color='blue')
        currentAxis=plt.gca()
        rect=patches.Rectangle( (coord_x[n]-1/2*sideLength,coord_y[n]-1/2*sideLength),sideLength,sideLength,linewidth=1,edgecolor='k',facecolor='none' )
        currentAxis.add_patch(rect)
    else:
        plt.scatter(x=coord_x[n],y=coord_y[n],color='green')
        currentAxis=plt.gca()
        rect=patches.Rectangle( (coord_x[n]-1/2*sideLength,coord_y[n]-1/2*sideLength),sideLength,sideLength,linewidth=1,edgecolor='k',facecolor='none' )
        currentAxis.add_patch(rect)
        
# obstacle part
for obstacle in obstacles:
    # use x0,y0 to be short for coord_x[obstacle],coord_y[obstacle]
    x0=coord_x[obstacle]
    y0=coord_y[obstacle]
    plt.scatter(x=x0,y=y0,color='red')
    #second, pathes.Rectangle
    currentAxis=plt.gca()
    rect=patches.Rectangle( (x0-1/2*sideLength,y0-1/2*sideLength),sideLength,sideLength,linewidth=1,edgecolor='k',facecolor='grey' )
    currentAxis.add_patch(rect)        
# obstacle part end

for i,j in draw_edges:
    plt.plot([coord_x[i],coord_x[j]],[coord_y[i],coord_y[j]],color='grey')

    
for n in range(len(coord_x)):# name of node
    plt.annotate(str(n),xy=(coord_x[n],coord_y[n]),xytext=(coord_x[n]+0.01,coord_y[n]+0.1),color='red')

#plt.grid(True)
plt.savefig("size 5 x 4.png")
#plt.savefig("size 5 x 4.eps")
plt.show()


draw_edges =[i for i in edges if x_values[i].solution_value>0.9]
optimal_path=[0]    
sad=0
iteration = len(draw_edges)
x=x_values
for i in range(iteration):
	for happy in NodesAndDeparturePoint:
    	if x[(sad, happy)].solution_value > 0.9:
        	optimal_path.append(happy)
            sad = happy
            break



optimal_length = len(optimal_path)
total_distance = 0
for index in range(0, optimal_length-1):
	total_distance = total_distance + distance[(optimal_path[index], optimal_path[index+1])]




total_degrees = 0
for index in range(0, optimal_length-2):
	total_degrees = total_degrees + degrees[(optimal_path[index]),(optimal_path[index+1]),(optimal_path[index+2])]