import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import addSimple_cuts
from parameters import *
import numpy as np


SOURCE = 's'
TARGET = 't'
MAX_CAPACITY = 1000




def createClusterGraphSimple(clusters,  wrapped_u):
    auxgraph = None
    auxgraph = nx.DiGraph()
    auxgraph.add_nodes_from(clusters)
    
    for e in wrapped_u.keys():
        auxgraph.add_edge(*e, capacity=wrapped_u[e])

    return auxgraph
    

def createConstraintSimple(auxgraph,p,q,wrapped_y, depot_cluster):
    
    def min_cut(graph, source, dest, to_remove, thresh, p, q):
        graph.remove_node(to_remove)
        cut_val, part = nx.minimum_cut(graph, source, dest)
        if cut_val < thresh:
            S,barS = part
            # cut = tuple(sorted(list(U - Pi_U))), tuple(sorted(list(V - Pi_U)))
            return (tuple(S), tuple(barS), p, q)
        else:
            return None

    DepP = auxgraph.copy()
    res1 = min_cut(DepP,depot_cluster,p,q, wrapped_y[p,q], p, q)
    
    PQ = auxgraph.copy()
    res2 = min_cut(PQ,p,q,depot_cluster,wrapped_y[p,q], p, q)
    
    QDep = auxgraph.copy()
    res3 = min_cut(QDep,q,depot_cluster,p,wrapped_y[p,q], p, q)
    
    result = [res for res in [res1, res2, res3] if res != None]
    return result



def generateSimple_cuts(model, depot_cluster=1):
    clusters = model._clusters
    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]

    aux_G = createClusterGraphSimple(clusters, wrapped_u)

    for p,q in permutations(cluslist,2):
        if wrapped_y[p,q] > 0:
            result = createConstraintSimple(aux_G, p, q, wrapped_y, depot_cluster)
            for cut in result:
                cuts.add(cut)

    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = addSimple_cuts(model,cuts)
    return new_model, number_of_cuts  

















