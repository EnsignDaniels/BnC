import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import add_2path_cuts
from parameters import *
import numpy as np


SOURCE = 's'
TARGET = 't'
MAX_CAPACITY = 1000


global_auxGraphsGDDL = {}



def createClusterGraph(clusters,  wrapped_u):
    auxgraph = nx.DiGraph()
    auxgraph.add_nodes_from(clusters)
        
    for e in wrapped_u.keys():
        auxgraph.add_edge(*e, capacity=wrapped_u[e])

    return auxgraph
    

def create_2path_Constraint(auxgraph,p,q,r,wrapped_y, depot_cluster):
    
    def make_a_path(auxgraph,source, dest, nodes_to_out, thresh):
        auxG = auxgraph.copy()
        auxG.remove_nodes_from(nodes_to_out)
        
        cut_val, part = nx.minimum_cut(auxG,source,dest)
        if cut_val < thresh:
            S,barS = part
            return tuple(S), tuple(barS)
        else:
            return None  
    
    thresh = wrapped_y[p,q] + wrapped_y[q,r] - 1
    p1 = make_a_path(auxgraph, depot_cluster, p, [q,r], thresh)
    p2 = make_a_path(auxgraph, p, q, [depot_cluster,r], thresh)
    p3 = make_a_path(auxgraph, q, r, [depot_cluster,p], thresh)
    p4 = make_a_path(auxgraph, r, depot_cluster, [p,q], thresh)
    
    return [(*pp, p,q,r) for pp in [p1, p2, p3, p4] if pp != None]




def generate_2path_cuts(model, depot_cluster=1):

    clusters = model._clusters
    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]
    aux_G = createClusterGraph(clusters, wrapped_u)
    
    for p,q,r in permutations(cluslist,3):
        #print(f'2-path: {p} {q} {r}')
        if wrapped_y[p,q] + wrapped_y[q,r] > 1:    
            paths  = create_2path_Constraint(aux_G, p, q, r, wrapped_y, depot_cluster)
            for cut in paths:
                cuts.add(cut)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = add_2path_cuts(model,cuts)
    return new_model, number_of_cuts  

















