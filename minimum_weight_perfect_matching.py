import networkx as nx
import copy
"""
Minimum weight perfect matching of a graph 'graph'
"""
def mwpm(graph):
    G = copy.deepcopy(graph)
    edge_weight = []
    for u,v in G.edges:
        edge_weight.append(G.edges[u,v]['weight'])
    max_weight = max(edge_weight, default = 9999999)
    for u,v in G.edges:
        G.edges[u,v]['weight'] = max_weight +1- G.edges[u,v]['weight']
    result =  nx.algorithms.max_weight_matching(G) 
    return list(result)

