import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import addGDDL_cuts
from parameters import *
import numpy as np

from sampler_v2 import get_a_sample, update_rates


SOURCE = 's'
TARGET = 't'
MAX_CAPACITY = 1000


#global_auxGraphsGDDL = {}



def createClusterGraphGDDL(p, q, r, clusters,  wrapped_u, depot_cluster, tree_closure):
    auxgraph = None
    auxgraph = nx.DiGraph()
    p_suc = set(tree_closure.successors(p))
    q_suc = set(tree_closure.successors(q))
    r_pre = set(tree_closure.predecessors(r))
    r_suc = set(tree_closure.successors(r))
    
    nodes_to_out = (p_suc & q_suc) | (r_pre & p_suc) | (q_suc & r_suc)  - {depot_cluster,p,q,r} 
    auxgraph.add_nodes_from(clusters)
    auxgraph.remove_nodes_from(nodes_to_out)
    
    auxgraph.add_node(SOURCE)
    auxgraph.add_node(TARGET)
    
    auxgraph.add_edge(SOURCE,depot_cluster,capacity=MAX_CAPACITY)
    auxgraph.add_edge(SOURCE,r,capacity=MAX_CAPACITY)
    auxgraph.add_edge(p,TARGET,capacity=MAX_CAPACITY)
    auxgraph.add_edge(q,TARGET,capacity=MAX_CAPACITY)
    
    for e in wrapped_u.keys():
        auxgraph.add_edge(*e, capacity=wrapped_u[e])

    return auxgraph
    

def createConstraintGDDL(auxgraph,p,q,r,depot_cluster,wrapped_y):
    cut_val, part = nx.minimum_cut(auxgraph,SOURCE,TARGET)
    if cut_val < wrapped_y[r,q] + wrapped_y[p,r]:
        S,barS = part
        S.remove(SOURCE)
        barS.remove(TARGET)
        
        if r in S and depot_cluster in S and p in barS and q in barS:
            return tuple(S), tuple(barS), p, q, r #wrapped_y[r,q] + wrapped_y[p,r]
    return None          



def generateGDDL_cuts(model, depot_cluster=1):
    clusters = model._clusters
    tree_closure = model._tree_closure
    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    epoch = model._epoch
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]
    if NEED_GDDL_CUTS_SAMPLING:
        rates = model._GDDL_rates
        combs = list(permutations(cluslist,3))
        iter = get_a_sample(epoch,cluslist,combs,rates, FRACTION_GDDL_CUTS)
    else:
        iter = permutations(cluslist,3) 

    counter = 0
    iter = [(p,q,r) for (p,q,r) in iter if wrapped_y[r,q] + wrapped_y[p,r] > 0]
    #print(f'{len(iter)} combinations found')

    for p,q,r in iter:
        counter += 1
        if wrapped_y[r,q] + wrapped_y[p,r] > 0:
            #print(f'{(p,q,r)}:  wrapped_y[{r},{q}] + wrapped_y[{p},{r}] = {wrapped_y[r,q] + wrapped_y[p,r]}')
            aux_G = createClusterGraphGDDL(p, q, r, clusters, wrapped_u, depot_cluster,tree_closure)
            cut = createConstraintGDDL(aux_G, p, q, r,depot_cluster,wrapped_y)
            if cut:
                cuts.add(cut)
                if NEED_GDDL_CUTS_SAMPLING:
                    #model._GDDL_rates[(p,q,r)] += RATES_GDDL_CUTS_STEP
                    update_rates(model._GDDL_rates, (p,q,r), RATES_GDDL_CUTS_STEP)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = addGDDL_cuts(model,cuts)
    return new_model, number_of_cuts, counter
    #return cuts, number_of_cuts  

















