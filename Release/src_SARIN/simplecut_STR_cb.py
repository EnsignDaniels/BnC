import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import addSimple_cuts
from parameters import *
import numpy as np

from sampler_v2 import get_a_sample, update_rates


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
    

def createConstraintSimple(auxgraph,p,q,wrapped_y, depot_cluster, tree_closure):
    
    def min_cut(graph, source, dest, nodes_to_remove, thresh, p, q):
        graph.remove_nodes_from(nodes_to_remove)
        cut_val, part = nx.minimum_cut(graph, source, dest)
        if cut_val < thresh:
            S,barS = part
            # cut = tuple(sorted(list(U - Pi_U))), tuple(sorted(list(V - Pi_U)))
            return (tuple(S), tuple(barS), p, q)
        else:
            return None

    DepP = auxgraph.copy()
    nodes_to_remove = set(tree_closure.successors(p)) | set(tree_closure.successors(q)) | {q}
    #print(f'nodes_to_remove = {nodes_to_remove}')
    res1 = min_cut(DepP,depot_cluster,p,nodes_to_remove, wrapped_y[p,q], p, q)
    #print(f'res1 = {res1}')
    
    PQ = auxgraph.copy()
    nodes_to_remove = set(tree_closure.predecessors(p)) | set(tree_closure.successors(q))
    #print(f'nodes_to_remove = {nodes_to_remove}')
    res2 = min_cut(PQ,p,q,nodes_to_remove, wrapped_y[p,q], p, q)
    #print(f'res2 = {res2}')
    
    QDep = auxgraph.copy()
    nodes_to_remove = set(tree_closure.predecessors(q)) | set(tree_closure.predecessors(p)) | {p}
    nodes_to_remove.remove(depot_cluster)
    #print(f'nodes_to_remove = {nodes_to_remove}')
    res3 = min_cut(QDep,q,depot_cluster,nodes_to_remove,wrapped_y[p,q], p, q)
    #print(f'res3 = {res3}')
        
    result = [res for res in [res1, res2, res3] if res != None]
    return result



def generateSimple_cuts(model, depot_cluster=1):
    clusters = model._clusters
    tree_closure = model._tree_closure
    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    epoch = model._epoch
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]

    if NEED_SIMPLE_CUTS_SAMPLING:
        rates = model._Simple_rates
        combs = list(permutations(cluslist,2))
        iter = get_a_sample(epoch, cluslist, combs, rates, FRACTION_SIMPLE_CUTS)
    else:
        iter = permutations(cluslist,2) 

    aux_G = createClusterGraphSimple(clusters, wrapped_u)
    counter = 0

    for p,q in iter:
        #print(p,q)
        counter += 1
        if wrapped_y[p,q] > 0:
            result = createConstraintSimple(aux_G, p, q, wrapped_y, depot_cluster, tree_closure)
            #print(f'result={result}')
            if result and NEED_SIMPLE_CUTS_SAMPLING:
                #model._Simple_rates[(p,q)] += RATES_SIMPLE_CUTS_STEP
                update_rates(model._Simple_rates, (p,q), RATES_SIMPLE_CUTS_STEP)
            for cut in result:
                cuts.add(cut)

    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = addSimple_cuts(model,cuts)
    return new_model, number_of_cuts, counter*3  
    #return cuts, number_of_cuts

















