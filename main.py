# calculate total_distance and total_degrees for the test map
# optimal_path replace gen_model_label


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