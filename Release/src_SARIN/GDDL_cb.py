import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import addGDDL_cuts
from parameters import *
import numpy as np


SOURCE = 's'
TARGET = 't'
MAX_CAPACITY = 1000


global_auxGraphsGDDL = {}



def createClusterGraphGDDL(p, q, r, clusters,  wrapped_u, depot_cluster):
    auxgraph = None
    auxgraph = nx.DiGraph()
    auxgraph.add_nodes_from(clusters)
    auxgraph.add_node(SOURCE)
    auxgraph.add_node(TARGET)
    
    auxgraph.add_edge(SOURCE,depot_cluster,capacity=MAX_CAPACITY)
    auxgraph.add_edge(SOURCE,r,capacity=MAX_CAPACITY)
    auxgraph.add_edge(p,TARGET,capacity=MAX_CAPACITY)
    auxgraph.add_edge(q,TARGET,capacity=MAX_CAPACITY)
    
    for e in wrapped_u.keys():
        auxgraph.add_edge(*e, capacity=wrapped_u[e])

    return auxgraph
    

def createConstraintGDDL(auxgraph,p,q,r,wrapped_y):
    cut_val, part = nx.minimum_cut(auxgraph,SOURCE,TARGET)
    if cut_val < wrapped_y[r,q] + wrapped_y[p,r]:
        S,barS = part
        S.remove(SOURCE)
        barS.remove(TARGET)
        # cut = tuple(sorted(list(U - Pi_U))), tuple(sorted(list(V - Pi_U)))
        return tuple(S), tuple(barS), p, q, r #wrapped_y[r,q] + wrapped_y[p,r]
    else:
        return None          



def generateGDDL_cuts(model, depot_cluster=1):
    clusters = model._clusters
    tree_closure = model._tree_closure
    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]
    for p,q,r in permutations(cluslist,3):
        if wrapped_y[r,q] + wrapped_y[p,r] > 0:
            aux_G = createClusterGraphGDDL(p, q, r, clusters, wrapped_u, depot_cluster)
            cut = createConstraintGDDL(aux_G, p, q, r, wrapped_y)
        if cut:
            cuts.add(cut)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = addGDDL_cuts(model,cuts)
    return new_model, number_of_cuts  

















