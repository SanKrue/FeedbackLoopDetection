library(deSolve)
library(rootSolve)
library(igraph)

get_all_loops <- function(jacobian, num_loops){ #function for loop detection in an ODE, argument: jacobian matrix of an ODE, and 													 maximal number of loops
	if(missing(num_loops)) { #if the user do not use the value for maximal number of loops, set the value to 1 Million
		num_loops <- 1000000
    }
	#get all loops and the sign of the loops
	tmp_path <- c() #variable to store the current path
	loop <- c() #variable to store the current loop
	all_loops <- list() #list to store all loops
	loop_signs <- c() #vector to store if the loops are positive or negative
	loop_length <- c() #vector to store the length of all loops
	species_done <- c(1)
	d <- dim(jacobian)[1]
	for (i in 1:d){ #iterate over the dimension of the jacobian
		if (jacobian[i,i] < 0){ #if the entry smaller than zero, there is a negative self loop
			loop_signs <- c(loop_signs,-1)
			loop_length <- c(loop_length, 1)
			loop <- c(i,i)
			all_loops <- c(all_loops, list(loop))
		}
		if (jacobian[i,i] > 0){ #if the entry greater than zero, there is a positive self loop
			loop_signs <- c(loop_signs,1)
			loop_length <- c(loop_length, 1)
			loop <- c(i,i)
			all_loops <- c(all_loops, list(loop))
		}
		else{ #if the entry is zero, there isn't a self loop
			next		
		}
	}

	for (i in 2:d){ #initialize a backward loop from the last row/column of the jacobian matrix to 2 (the last 										matrix is then a 2x2 matrix)
		if (length(all_loops) <= num_loops){
			tmp_jacobian <- jacobian[1:i,1:i] #to keep the original jacobian, store the current jacobian in a temporary matrix
			t_tmp_jacobian <- t(tmp_jacobian)
			jac_graph <- graph_from_adjacency_matrix(t_tmp_jacobian,mode="directed",diag=TRUE, weighted=TRUE,  add.colnames = NA) 																		#build the directed graph from the current jacobian matrix
			tmp_path <- all_simple_paths(jac_graph,i) #search all simples paths from the current node to all other nodes
			if (length(tmp_path)!=0){ #only interate over the loop if the list is not empty			
				for (j in length(tmp_path):1){ #all possible paths from the current node are stored in a list called tmp_path, 													iterate over all paths in the list, starting from the last path (this is the 													longest path in the list)
					first_node <- as_ids(tmp_path[[j]][1]) #the first node of the loop is the first entry (of the vector) in the 														        current list entry (as_ids is a function from igraph to convert a 																vertex or edge sequence to an ordinary vector)


					last_node <- as_ids(tmp_path[[j]][length(tmp_path[[j]])]) #the last node is the last entry of the vector in 																			   the current list entry
					if (t_tmp_jacobian[last_node,first_node] != 0){ #if there is an entry at the position of the first and last 																	 node in the current jacobian which is not 0, then there is an 																		 edge between these nodes (e.g. an entry != 0 at the position 																		 4,1 in the jacobian means that there is an edge between node 1 																	 and node 4); if this is the case, then there is a loop
						loop <- c(last_node, as_ids(tmp_path[[j]])) #store the current loop in the variable (e.g. last node = 1, 																         tmp_path[[j]] = 4 3 2 1 --> loop = 1 4 3 2 1)
						all_loops <- c(all_loops, list(loop)) #store the loop in the list of all paths
						loop_length <- c(loop_length, (length(loop)-1)) #store the loop length (is is loop length-1, because the 																			edges should be counted and not the nodes in the loop 																			vector)
						edge_vec <- rep(loop, each=2)[-1] #repeat every entry of the vector loop except the first entry
						edge_vec <- edge_vec[-length(edge_vec)] #delete last entry 	
						if (prod(E(jac_graph)$weight[get.edge.ids(jac_graph, edge_vec)])>0){ #build the product of all weights of 																								  the edges which are part of the 																								  current loop and proof if this 																								  product is greater than 0 
							loop_signs <- c(loop_signs,1) #1 if the loop is positive
						}
						else{
							loop_signs <- c(loop_signs,-1) #-1 if the loop is negative
						}
					}
				}
			}
		species_done <- c(species_done,i)
		}
		#if the length of the loop list is greater than the value for number of loops the user gets a warning
		else{
			warning(paste('More loops found than specified. Please adjust the amount of loops.' , sprintf('Loops for the first %s of %d species found.', species_done[length(species_done)], d)))
		}
	}
result <- data.frame(length=loop_length,sign=loop_signs) #store the each loop length and sign of the loop in a dataframe
result$loop <- all_loops #store all loops as new column "loop" in the dataframe
result <- result[c("loop","length","sign")] #reorder the data frame
return(result) #return the data frame
}

get_edges <- function(loop_list,from,to){ #function for edge detection in a list of loops, arguments: list of loops (loops), 		                                       starting node (from), end node (to)
	bool <- c() #vector to store booleans
	for (k in 1:length(loop_list)){
		tmp_loop <- unlist(loop_list[k]) #unlist each loop
		edge_vec <- rep(tmp_loop, each=2)[-1] #repeat every entry of the vector loop except the first entry
		edge_vec <- edge_vec[-length(edge_vec)] #delete last entry
		edge_df <- as.data.frame(edge_vec) #store the vector as data frame
		edges <- cbind(edge_df[seq(1,dim(edge_df)[1]-1,by=2),],edge_df[seq(2,dim(edge_df)[1],by=2),])#store the edges in   																									a matrix, "from" is the first 																									column, "to" the second column 
		if (any(edges==from)){ #try if the matrix "edges" contain the value "from"
			bool <- c(bool,edges[which(edges[,1]==from),2]==to) #store in the vector "bool" a "TRUE" if the value "from" is in the 																	 first column of the matrix "edges" and if the value "to" is in 																 the same row and in the second column and "FALSE" if this is not  																	 the case
		}
		else{ #if the value "from" is not in the matrix store a "FALSE" in the vector "bool"
			bool <- c(bool,FALSE)
		}
	}
return(bool) #return the boolean vector
}


sort_loops <- function(loop_list){ #function to sort each loop within a list, so that every loop starts with the smallest node
	sorted_loops <- list() #list for sorted loops
	for (i in 1:length(loop_list)){ 
		sorted_loop <- c() #vector for each loop
		ind_min <- which.min(unlist(loop_list[i])) #store the index of the smallest node
		sorted_loop <- c(sorted_loop,unlist(loop_list[i])[ind_min:length(unlist(loop_list[i]))]) #store the values from the index 																									  of the smallest node to the end 																									  of the loop
		if (length(sorted_loop) < length(unlist(loop_list[i]))){ #if the sorted_loop is smaller than the whole loop, store the 																	  values from the second entry to the index of the smallest node 																	  minus 1
			sorted_loop <- c(sorted_loop, unlist(loop_list[i])[2:(ind_min-1)], unlist(loop_list[i])[ind_min])
		}
		sorted_loops <- c(sorted_loops, list(sorted_loop)) #store the sorted loop in the list of all sorted loops
	}
return(sorted_loops)
}
