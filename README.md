**Implementation 1:**

Files: 

•	approxDP: This is the main function, and the program is run by executing the command:

**python3 approxDP.py 5 10000 graph_sample_10k**

in which 5 is the value of k (number of segments), and 10000 is the limit of file lines that we want to be read.

We can also run the command:

**python3 approxDP.py 5 10000000 main_graph**

to execute on the main graph which has slightly less than 10 million lines. But with this command change the parameter **Span** (generator_pars['wholeSpan']) manually in approxDP.py to generate a random graph with 10 million lines.

In this file, first a dictionary based on the input file is used such that its key will be the unique timestamps and edges are assigned to these timestamps. Then we have these procedures that are executed respectively:

1)	Class initialization: The ApproxDP class is initialized with parameters eps (approximation constant value), k (the number of desired subgraphs), and timestamps (a dictionary where keys are timestamps and values are lists of edges).
2)	Running DP: The run_DP method is the core dynamic programming algorithm. It uses an incremental densest subgraph algorithm (IncDensest) to efficiently compute the density of subgraphs at each timestamp. An instance of IncDensest is created for each timestamp, and the densest subgraph is updated by adding edges for that timestamp. The best_density obtained from the densest subgraph at each timestamp is used in the dynamic programming process. The algorithm considers different segments and timestamps, storing the results in a dynamic programming table (DP).
3)	SPRS: The __SPRS method (Semi-Partition Removal Scheme) is a subroutine used in the dynamic programming process. It prunes unnecessary subproblems based on the computed densities.
4)	The get_sol_intervals method retrieves the starting indices of the selected subgraphs from the dynamic programming table.
5)	The get_sol_graphs method constructs the densest subgraphs using charikarHeap method based on the selected intervals returned from get_sol_intervals function and computes their densities.

•	incremental_densest: The IncDensest class is an implementation of an incremental algorithm  for finding densest subgraphs in a dynamic graph. Its procedures are as follows:

1)	Class initialization: The class is initialized with parameters eps (approximation constant value), G (a NetworkX graph representing the current state of the evolving graph), node2set (a dictionary mapping nodes to their corresponding set (cluster) index), delta (a dictionary to keep track of changes in degrees for each node and set), best_density and best_iter (used to store the best density and the corresponding iteration during the process).
2)	Adding nodes and edges: this is done by __add_new_node method and add_edge method that adds an edge between two nodes, updating the graph and triggering a potential reevaluation of the densest subgraph. Then returns the current best density.
3)	Densest subgraph computation: The __add method is used to update the state of the graph and check whether the densest subgraph needs to be recomputed. It adjusts degrees and sets, potentially leading to rebuilding the densest subgraph. Then the __find_densest method iteratively refines the densest subgraph's density until convergence. This method uses the __find method that employs a Fibonacci heap data structure to efficiently find the densest subgraph. It iteratively removes nodes with the lowest degrees, updating the degrees of remaining nodes. Then __get_degrees method updates the delta dictionary to reflect the current degrees of nodes with respect to their sets.
•	utils: Some utility functions in which readdata_dict_limit function is utilized to read data from a file containing timestamped edges, and it stores the information in a dictionary.

•	charikar: In this file, the charikarHeap function is used to obtain the densest subgraphs for each selected interval. The function iteratively removes nodes with the lowest degrees from the graph S using the Fibonacci heap. This is done until the graph is empty. At each iteration, the average degree (avg_degree) is computed, and if it is greater than the current best average degree (best_avg), it updates best_avg and best_iter to store the information about the iteration with the densest subgraph so far. The removed nodes are stored in the order list. After the iterations, a new copy of the input graph G is created (S) and nodes up to the best iteration (best_iter) are removed to obtain the final densest subgraph. The function returns the densest subgraph (S) and the corresponding best average degree (best_avg).

•	fibonacci_heap_mod: An official implementation of a priority queue backed by a Fibonacci heap.

Now, after identifying dense subgraphs and specifying their intervals and densities (some prints and plots are done), then we need to compare the results with random graphs to verify the usefulness of the approach and the significance of the results. The random graphs should be generated by keeping some metrics same as the original graph.

•	generate: To do the above-mentioned task, we used the function generate which works based on the following concepts. 
1.	Background Graph (G):
•	The background graph G is generated using the generateRG function. This function creates a random graph with a specified number of nodes (n), edge probability (edges_p), and node labels (nodes). The fast_gnp_random_graph function from NetworkX is used for this purpose. This graph serves as the overall structure or skeleton onto which additional structures will be added.
•	The edges of the background graph G are then shuffled, and a time series (TS) is constructed by associating each edge with a timestamp. This process is implemented in the plantBackground_uniform function.
•	The idea of having a background graph is to provide a baseline structure to the synthetic network. This structure may represent some underlying connectivity pattern that is common across the entire network.
2.	Planted Graphs (C):
•	The code then generates community structures, called "planted graphs" or "C," using the generateRG function again. Each community structure C is associated with a specific time interval.
•	The plant_uniform function is used to "plant" these community structures into the background graph by associating edges from the community structures with specific timestamps within their respective intervals.
•	The motivation for adding community structures is to simulate the presence of groups or communities within the network. Each community structure is generated separately, and the edges are distributed across time to create a time-evolving network.

Note that you can tune the parameter **Noise** (generator_pars['noise']) manually in approxDP.py.

**Implementation 2:**
This program is run by typing the command:

**python3 approxWeighted.py 5 10000 graph_sample_10k**

The overall methodology is the same as the previous implementation except the following:

•	find_densest function: This function computes the densest subgraph in a given graph by solving a minimum cut problem. It sets up a directed graph H where each edge has a capacity. The minimum cut is found using the nx.minimum_cut function from NetworkX, and the density of the resulting cut is used to update the lower and upper bounds for the densest subgraph. The algorithm iteratively refines the bounds until they converge within a small epsilon.

•	IncDensest Class: This class now has a G attribute representing the evolving graph. The best_density attribute is updated by calling the find_densest function each time a new edge is added through the add_edge method. The add_edge method checks if the nodes and edge are already present in the graph and adds them if not. Then, it updates the best_density by recomputing the densest subgraph using the find_densest function.

