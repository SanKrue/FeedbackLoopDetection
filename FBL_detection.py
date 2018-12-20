import networkx as nx
import numpy as np
import pandas as pd

def get_all_loops(jacobian):
	df = pd.DataFrame(columns=["loop","length","sign"]) #generate an empty pandas data frame
	g = nx.DiGraph(np.transpose(jacobian)) #get a directed graph from the transposed jacobian matrix
	cycles = nx.simple_cycles(g) #get all cycles out of the graph
	weights = nx.get_edge_attributes(g, 'weight') #store the weights of the edges
	for cycle in cycles: #for each cycle in the list cycles
		cycle.append(cycle[0]) #append the first entry
		w = np.prod([weights[(cycle[i-1], cycle[i])] for i in range(1, len(cycle))]) #calculate the product of the weights to get 																						  the loop sign
		if w<0: #if the weight is less than 0, store -1 for the sign
			sign = -1
		else: #if the weight is greater than 0, store 1 for the sign
			sign = 1 
		length = len(cycle)-1 #store length -1 in the list 
		df.loc[len(df)] = [cycle, length, sign] #put all lists in the data frame
	return df #return data frame


