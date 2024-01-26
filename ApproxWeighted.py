import numpy as np
import networkx as nx
import incremental_densest_weighted as opid
import os.path
import utils
import time
import matplotlib.pyplot as plt
import sys
import generate_weighted
import evaluation
from charikar import charikarHeap_weighted

    
class ApproxDP:
    def __init__(self, eps, k, timestamps):
        self.k = k
        self.eps = eps
        self.timestamps = timestamps
        self.sorted_TS = sorted(timestamps.keys())
        self.m = len(self.sorted_TS)
        self.DP = np.empty((self.k, self.m))
        self.C = np.empty((self.k, self.m))
        self.AD = {}
    
    def run_DP(self):
        ds = opid.IncDensest()
        print("Segment =  0")
        for i in range(0, self.m):
            for (n1,n2,weight) in self.timestamps[self.sorted_TS[i]]:
                best_density =  ds.add_edge(n1, n2, weight)
            self.DP[0, i] = best_density
            self.C[0, i] = 0
            
        for l in range(1, self.k):
            self.AD = {}
            print("Segment = ",l)
            for i in range(0, self.m):
                self.AD[i] = opid.IncDensest()
                best_candidate = -1
                best_index = -1
                for idx in self.AD: 
                    for (n1,n2,weight) in self.timestamps[self.sorted_TS[i]]:
                        best_density =  self.AD[idx].add_edge(n1, n2, weight)
                    candidate = self.DP[l-1, idx-1] + best_density if idx > 0 else best_density
                    if candidate >= best_candidate:
                        best_candidate = candidate
                        best_index = idx                
                if i > 0 and self.DP[l, i-1] > self.DP[l-1, i] and self.DP[l, i-1] > best_candidate:
                    best_candidate = self.DP[l, i-1]
                    best_index = self.C[l, i-1]
                elif ((i > 0 and self.DP[l-1, i] > self.DP[l, i-1]) or i == 0) and self.DP[l-1, i] > best_candidate:
                    best_candidate = self.DP[l, i-1]
                    best_index = -1
                    
                self.DP[l, i] = best_candidate
                self.C[l, i] = best_index
               
                self.__SPRS(best_candidate, l)        

    def __SPRS(self, best_candidate, l):
        delta = best_candidate * self.eps / (self.k + l * self.eps)
        sorted_A = sorted(self.AD.keys()) 
        j = 0
        while j <= len(sorted_A) - 3:
            
            if self.DP[l-1, sorted_A[j+2]] - self.DP[l-1, sorted_A[j]] <= delta:
                del self.AD[sorted_A[j+1]]
                del sorted_A[j+1]
            else:
                j += 1               
        return
        
        
    def get_sol_intervals(self):
        starts = []
        end = self.m - 1
        for l in range(self.k-1,-1, -1):
            choice = int(self.C[l, end])
            if choice > -1:
                starts.append(choice)
                end = choice -1
            if end < 0:
                return starts[::-1]
        return starts[::-1]
    
        
    def get_sol_graphs(self):
        starts = self.get_sol_intervals()
        graphs = [-1]*len(starts)
        densities = [-1]*len(starts)
        intervals = [-1]*len(starts)
        for s in range(len(starts)):
            G = nx.Graph()
            end = len(self.sorted_TS)-1 if s == len(starts)-1 else starts[s+1]-1
            intervals[s] = (starts[s], end)
            for j in range(starts[s], end+1):
                for (n1, n2, weight) in self.timestamps[self.sorted_TS[j]]:
                    G.add_edge(n1, n2 , weight=weight)
            G = nx.Graph(G)
            G.remove_edges_from(nx.selfloop_edges(G))
            S, best_avg = charikarHeap_weighted(G)
            graphs[s] = S
            densities[s] = best_avg
        return graphs, densities, intervals

        
if __name__ == "__main__":

    assert len(sys.argv) == 4, "Usage: python3 approxWeighted.py 5 10000 graph_sample_10k" 
    #python3 approxWeighted.py 5 10000000 main_graph
    
    # INPUT READING
    k = sys.argv[1]  # 5, 100, 10
    assert k.isdigit(), "k must be an integer"
    k = int(k)
    limit = sys.argv[2] # limit in lines of code
    assert limit.isdigit(), "limit must be an integer"
    limit = int(limit)
    dataset = sys.argv[3] # main_graph / graph_sample_10k

    dp_eps = 0.1
    
    alg_pars = {}
    alg_pars['dp_eps'] = dp_eps
    
    filename = os.path.join('data', dataset+'.txt')
    TS, unique_original_timestamps = utils.readdata_dict_limit_weight(filename, limit)
 
    t1 = time.time()
    ADP = ApproxDP(alg_pars['dp_eps'], k, TS)
    ADP.run_DP()
    t2 = time.time()
    print ('Running time:', t2 - t1)
    
    print("Number of timestamps = ",len(unique_original_timestamps))
    
    graphs, densities, intervals = ADP.get_sol_graphs()
    print ('Densities:', densities)
    # Print number of edges and nodes for each graph
    print('Nodes and edges in graphs:')
    for i in graphs:
        print(f'{i}')
        
    print("Time intervals:")
    for interval in intervals:
            a,b = interval
            print(unique_original_timestamps[a],"-",unique_original_timestamps[b])
    max_index = densities.index(max(densities))
    s,e = intervals[max_index]
    print("Maximum Density =", max(densities),", in the interval =",unique_original_timestamps[s],"-",unique_original_timestamps[e])

    for n in range(len(graphs)):
        # Plot the graph using matplotlib
        nx.draw(graphs[n], with_labels=True, font_weight='bold')
        # Show the plot for the current graph without waiting
        plt.figure()
        plt.show(block=False)
        plt.pause(0.1)
        # Show the plot
    plt.show()
    ####################################################################
    #GENERATING RANDOM GRAPH
    print("****Generating Random Graph****")
    n =  np.average([graphs[i].number_of_nodes() for i in range(len(graphs))])
    d = np.average(densities)

    generator_pars = {}
    generator_pars['k'] = k #10 #2 #5
    generator_pars['noise'] = 1  #0.5 # by lowering this the densities will be more accurate but by increasing to 100 it will be bad
    generator_pars['innerdegree'] = round(d)
    generator_pars['nodesInCom'] = round(n) 
    generator_pars['backgoundN'] = round(n) #1000 #100
    generator_pars['wholeSpan'] = 10000 #10000 #1000
        
    alg_pars = {}
    alg_pars['dp_eps'] = 0.1

    
    TS_random, backNoise, innerNoise, generated_C_random = generate_weighted.generate(generator_pars)
    t1 = time.time()
    ADP = ApproxDP(alg_pars['dp_eps'], generator_pars['k'], TS_random)
    ADP.run_DP()
    t2 = time.time()
    print ('Running time:', t2 - t1)
    
    graphs_random, densities_random, intervals_random = ADP.get_sol_graphs()
    #print ('intervals: ', intervals_random)
    print('Nodes and edges in graphs:')
    for i in graphs_random:
        print(f'{i}')
    print ('densities: ', densities_random)
    print ('total density:', sum(densities_random))
    
    evaluation.calculate_z_scores(densities,densities_random)
    

    