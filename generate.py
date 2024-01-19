import networkx as nx
import numpy as np
import random
from copy import deepcopy
    
def plant_uniform(interval, C):

    edges = list(C.edges())
    random.shuffle(edges)
    TS = []
    t = int(interval[0])
    for i in range(int(interval[1]) + 1 - int(interval[0])):
        if i == len(edges):
            return TS
        TS.append((t, edges[i]))
        t += 1   
    return TS
    
def plantBackground_uniform(span, G):
    edges = list(G.edges())
    
    TS = []
    for i in range(span):
        if i%len(edges) == 0:
            random.shuffle(edges)
        TS.append((i, edges[i%len(edges)]))    
    return TS
    
        
def generateRG(n, edges_p, nodes):
    G = nx.fast_gnp_random_graph(n, edges_p, seed=None, directed=False)
    gnodes = list(G.nodes())
    relabel = {gnodes[i]: nodes[i] for i in range(n)}    
    G = nx.relabel_nodes(G, relabel) 
    return G
    
def generate(generator_pars):
    
    k = generator_pars['k']
    #B = generator_pars['B']  #length of each interval
    noise = generator_pars['noise']  #background average degree
    innerdegree = generator_pars['innerdegree']  #plant average degree
    nodesInCom = generator_pars['nodesInCom'] #plant number of nodes
    backgoundN = generator_pars['backgoundN'] #background nodes
    wholeSpan = generator_pars['wholeSpan']
    
    com_number = k
    
    nCom = nodesInCom
    desired_avgCom = innerdegree

    n = backgoundN
    desired_avg_degree = noise   
        
    start = 0
    end= 1
    span = wholeSpan

    TS = []
    edges_p = float(desired_avg_degree)/(n-1)
    G = generateRG(n, edges_p, range(n))
    TS = plantBackground_uniform(span, G)

    innerNoise = []
    allocated_interval_length = wholeSpan/k
    generated_C = []
    for i in range(com_number):
        start = allocated_interval_length * i
        end = (allocated_interval_length * (i+1)) -1
        #interval = (allocated_interval_length*(i), allocated_interval_length*(i)+B-1)
        interval = (start, end)
        edges_pCom = float(desired_avgCom)/(nCom - 1)
        C = generateRG(nCom, edges_pCom, range(i*nCom,(i+1)*nCom))
        
        generated_C.append((deepcopy(C), interval))
        TS += plant_uniform(interval, C)  
        t = 2.0*C.number_of_edges()/C.number_of_nodes()
        innerNoise.append(t)
            
    TS.sort()
    TS_out = {}
    for i in TS:
        if i[0] not in TS_out:
            TS_out[i[0]] = []
        TS_out[i[0]].append(i[1])
        
    backNoise = 2.0*G.number_of_edges()/G.number_of_nodes()
    innerNoiseout = np.mean(innerNoise)
    
    return TS_out, backNoise, innerNoiseout, generated_C
