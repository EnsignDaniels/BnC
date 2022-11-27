import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import add_4vGDDL_cuts
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
    

def create_4vGDDL_Constraint(auxgraph,p,q,k,r,s, wrapped_y, depot_cluster):
    
    def make_a_path(auxgraph,connect_with_source, connect_with_target, nodes_to_out, thresh):
        auxG = auxgraph.copy()
        auxG.remove_nodes_from(nodes_to_out)
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
            
            
    thresh = wrapped_y[p,q] + wrapped_y[q,k] + wrapped_y[k,r] + wrapped_y[r,s] - 2
    
    p1 = make_a_path(auxgraph, [p,r], [q,s], [depot_cluster,k], thresh)
    #p2 = make_a_path(auxgraph, [depot_cluster,q,r], [p,k,s], [], thresh)
        
    
    return [(*pp,p,q,k,r,s) for pp in [p1] if pp != None]




def generate_4vGDDL_cuts(model, depot_cluster=1):
    clusters = model._clusters

    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]
    aux_G = createClusterGraph(clusters, wrapped_u)
    for p,q,k,r,s in permutations(cluslist,5):
        if wrapped_y[p,q] + wrapped_y[q,k] + wrapped_y[k,r] + wrapped_y[r,s] > 2:
            paths  = create_4vGDDL_Constraint(aux_G, p,q,k,r, s, wrapped_y, depot_cluster)
            for cut in paths:
                cuts.add(cut)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = add_4vGDDL_cuts(model,cuts)
    return new_model, number_of_cuts  

















