import networkx as nx
import random


def get_solution(G, g=None):
    if not g:
         g = find_densest(G)

    m, n = 1. * G.number_of_edges(), 1. * G.number_of_nodes()
    degrees = G.degree()
    nodes = G.nodes()

    H = G.to_directed()
    for u, v in H.edges():
        H[u][v]['cap'] = G[u][v]['weight'] #if 'weight' in G[u][v] else 1.0

    for v in nodes:
        # Compute node weight based on the sum of weights of incident edges
        incident_edge_weights = [G[u][v]['weight'] for u in G.neighbors(v) ]
        weight_v = sum(incident_edge_weights) 

        H.add_edge('s', v, cap= m + weight_v)
        weight_t =  m + 2. * g - weight_v * degrees[v]
        H.add_edge(v, 't', cap=weight_t)

    cut_value, (S, T) = nx.minimum_cut(H, 's', 't', capacity='cap')
    S.remove('s')
    D = G.subgraph(S)
    return D

def find_densest(G, g=None):
    ub = 0.5 * (G.number_of_nodes() - 1.)
    lb = 0.
    if not g:
        g = (lb + ub) / 2.
    m, n = 1. * G.number_of_edges(), 1. * G.number_of_nodes()
    degrees = G.degree()
    nodes = G.nodes()

    H = G.to_directed()
    for u, v in H.edges():
        H[u][v]['cap'] = G[u][v]['weight'] #if 'weight' in G[u][v] else 1.0

    for v in nodes:
        # Compute node weight based on the sum of weights of incident edges
        incident_edge_weights = [G[u][v]['weight'] for u in G.neighbors(v) ]
        weight_v = sum(incident_edge_weights)
        
        H.add_edge('s', v, cap= m + weight_v)
        weight_t =  m + 2. * g - weight_v * degrees[v]
        H.add_edge(v, 't', cap=weight_t)

    eps = 1e-8

    while ub - lb > eps:
        cut_value, (S, T) = nx.minimum_cut(H, 's', 't', capacity='cap')
        if len(S) - 1 > 0:
            lb = g
        else:
            ub = g
        g = (lb + ub) / 2.
        for v in nodes:
            # Compute node weight based on the sum of weights of incident edges
            incident_edge_weights = [G[u][v]['weight'] for u in G.neighbors(v)]
            weight_v = sum(incident_edge_weights) 
            
            H[v]['t']['cap'] =  m + 2. * g - weight_v * degrees[v]
    return lb

if __name__ == "__main__":
    G = nx.erdos_renyi_graph(100, 0.7)
    
    for u, v in G.edges():
        G[u][v]['weight'] = random.uniform(0.1, 1.0)  # Assign random weights in the range [0.1, 1.0]
    
    
    print("Original Graph Density:", 1. * G.number_of_edges() / G.number_of_nodes())
    
    densest_density = find_densest(G)
    D = get_solution(G, densest_density)
    
    print("Densest Subgraph Density:", 1. * D.number_of_edges() / D.number_of_nodes())