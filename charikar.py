import networkx as nx
import fibonacci_heap_mod
import copy


def charikarHeap_weighted(G):
    import copy
    import fibonacci_heap_mod

    E = G.number_of_edges()
    N = G.number_of_nodes()
    fib_heap = fibonacci_heap_mod.Fibonacci_heap()
    entries = {}
    order = []
    S = copy.deepcopy(G)
    
    for node in G.nodes():
        sum_weights = -sum([G[node][neighbor]['weight'] for neighbor in G.neighbors(node)])
        entries[node] = fib_heap.enqueue(node, sum_weights)
        
    best_avg = 0.0    
    iter = 0
    
    while fib_heap:
        avg_weight = (2.0 * E) / N
        if best_avg <= avg_weight:
            best_avg = avg_weight
            best_iter = iter
            
        min_weight_obj = fib_heap.dequeue_min()
        min_weight_node = min_weight_obj.get_value()
        order.append(min_weight_node)
        
        for n in S.neighbors(min_weight_node):            
            fib_heap.decrease_key(entries[n], 1)
        
        S.remove_node(min_weight_node)
        E -= sum([G[min_weight_node][neighbor]['weight'] for neighbor in G.neighbors(min_weight_node)])
        N -= 1
        iter += 1
        
    S = copy.deepcopy(G)       
    for i in range(best_iter):
        S.remove_node(order[i])
    return S, best_avg




def charikarHeap(G):
         
    E = G.number_of_edges()
    N = G.number_of_nodes()
    fib_heap = fibonacci_heap_mod.Fibonacci_heap()
    entries = {}
    order = []
    S = copy.deepcopy(G)
    
    for node, deg in G.degree():
        entries[node] = fib_heap.enqueue(node, deg)
        
    best_avg = 0.0    
    iter = 0
    
    while fib_heap:
        avg_degree = (2.0 * E)/N
        if best_avg <= avg_degree:
            best_avg = avg_degree
            best_iter = iter
            
        min_deg_obj = fib_heap.dequeue_min()
        min_deg_node = min_deg_obj.get_value()
        order.append(min_deg_node)
        for n in S.neighbors(min_deg_node):            
            fib_heap.decrease_key(entries[n], 1)
        
        
        S.remove_node(min_deg_node)
        E -= min_deg_obj.get_priority()
        N -= 1
        iter += 1
        
        
    S = copy.deepcopy(G)       
    for i in range(best_iter):
        S.remove_node(order[i])
    return S, best_avg
    
def charikar_densest(G):
         
    E = G.number_of_edges()
    N = G.number_of_nodes()
    fib_heap = fibonacci_heap_mod.Fibonacci_heap()
    entries = {}
    order = []
    S = copy.deepcopy(G)
    
    for node, deg in G.degree():
        entries[node] = fib_heap.enqueue(node, deg)
        
    best_avg = 0.0    
    iter = 0
    
    while fib_heap:
        avg_degree = (1.0 * E)/N
        if best_avg <= avg_degree:
            best_avg = avg_degree
            best_iter = iter
            
        min_deg_obj = fib_heap.dequeue_min()
        min_deg_node = min_deg_obj.get_value()
        order.append(min_deg_node)
        for n in S.neighbors(min_deg_node):            
            fib_heap.decrease_key(entries[n], 1)
        
        
        S.remove_node(min_deg_node)
        E -= min_deg_obj.get_priority()
        N -= 1
        iter += 1
        
    return best_avg
    
    
def charikarDicts(G):
 
    S = copy.deepcopy(G)
    
    E = G.number_of_edges()
    N = G.number_of_nodes()
    
    nodes = {}
    best_avg = 0.0    
    iter = 0
    order = []
    
    for node, deg in G.degree():
        nodes[node] = S[node]
        
    while nodes.keys():
        avg_degree = (2.0 * E)/N
        
        if best_avg <= avg_degree:
            best_avg = avg_degree
            best_iter = iter
            
        min_deg = N
        for n, neigh in nodes.iteritems():
            if min_deg > len(neigh):
                min_deg = len(neigh)
                min_deg_node = n
        
        order.append(min_deg_node)
            
        for neigh in nodes[min_deg_node].keys():
            del nodes[neigh][min_deg_node] 
            
        del nodes[min_deg_node]
        E -= min_deg
        N -= 1
        iter += 1
    
    S = copy.deepcopy(G)
    for i in range(best_iter):
        S.remove_node(order[i])
    return S, best_avg
    