import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import add_3vGDDL_cuts
from parameters import *
import numpy as np


SOURCE = 's'
TARGET = 't'
MAX_CAPACITY = 1000


global_auxGraphsGDDL = {}



def createClusterGraph(clusters,  wrapped_u):
    auxgraph = nx.DiGraph()
    auxgraph.add_nodes_from(clusters)
    auxgraph.add_nodes_from([SOURCE,TARGET])
        
    for e in wrapped_u.keys():
        auxgraph.add_edge(*e, capacity=wrapped_u[e])

    return auxgraph
    

def create_3vGDDL_Constraint(auxgraph,p,q,r,s, wrapped_y, depot_cluster):
    
    def make_a_path(auxgraph,connect_with_source, connect_with_target, node_to_out, thresh):
        auxG = auxgraph.copy()
        auxG.remove_node(node_to_out)
        for nd in connect_with_source:
            auxG.add_edge(SOURCE, nd, capacity=MAX_CAPACITY)
        for nd in connect_with_target:
            auxG.add_edge(nd, TARGET, capacity=MAX_CAPACITY)
        
        cut_val, part = nx.minimum_cut(auxG,SOURCE,TARGET)
        if cut_val < thresh:
            S,barS = part
            S.remove(SOURCE)
            barS.remove(TARGET)
            return tuple(S), tuple(barS)
        else:
            return None  
    
    thresh = wrapped_y[p,q] + wrapped_y[q,r] + wrapped_y[r,s] - 1
    
    p1 = make_a_path(auxgraph, [p,r], [q,s], depot_cluster, thresh)
    p2 = make_a_path(auxgraph, [p,s], [q,depot_cluster], r, thresh)
    p3 = make_a_path(auxgraph, [depot_cluster,r], [p,s], q, thresh)
        
    
    return [(*pp, p,q,r,s) for pp in [p1, p2, p3] if pp != None]




def generate_3vGDDL_cuts(model, depot_cluster=1):
    clusters = model._clusters

    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]
    aux_G = createClusterGraph(clusters, wrapped_u)
    for p,q,r,s in permutations(cluslist,4):
        if wrapped_y[p,q] + wrapped_y[q,r] + wrapped_y[r,s] > 1:
            paths  = create_3vGDDL_Constraint(aux_G, p, q, r, s, wrapped_y, depot_cluster)
            for cut in paths:
                cuts.add(cut)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = add_3vGDDL_cuts(model,cuts)
    return new_model, number_of_cuts  

















