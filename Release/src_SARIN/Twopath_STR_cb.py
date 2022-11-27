import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import add_2path_cuts
from parameters import *
import numpy as np

from sampler_v2 import get_a_sample, update_rates


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
    

def create_2path_Constraint(auxgraph,p,q,r,wrapped_y, depot_cluster, tree_closure):
    
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
    
    p_suc = set(tree_closure.successors(p))
    q_suc = set(tree_closure.successors(q))
    r_suc = set(tree_closure.successors(r))
    
    p_pre = set(tree_closure.predecessors(p))
    q_pre = set(tree_closure.predecessors(q))
    r_pre = set(tree_closure.predecessors(r))
    
    nodes_to_out = {nd for nd in p_suc | q_suc | r_suc | {q,r} if not nd in [depot_cluster,p]}
    #print(f'nodes_to_out: {nodes_to_out}')
    p1 = make_a_path(auxgraph, depot_cluster, p, nodes_to_out, thresh)
    #p1 = make_a_path(auxgraph, depot_cluster, p, [q,r], thresh)
    
    nodes_to_out = {nd for nd in p_pre | q_suc | r_suc | {depot_cluster, r} if not nd in [p,q]}
    #print(f'nodes_to_out: {nodes_to_out}')
    p2 = make_a_path(auxgraph, p, q, nodes_to_out, thresh)
    #p2 = make_a_path(auxgraph, p, q, [depot_cluster,r], thresh)
    
    nodes_to_out = {nd for nd in q_pre | r_suc | p_pre | {depot_cluster, p} if not nd in [q,r]}
    #print(f'nodes_to_out: {nodes_to_out}')
    p3 = make_a_path(auxgraph, q, r, nodes_to_out, thresh)
    #p3 = make_a_path(auxgraph, q, r, [depot_cluster,p], thresh)
    
    nodes_to_out = {nd for nd in r_pre | p_pre | q_pre | {p,q} if not nd in [depot_cluster, r]}
    #print(f'nodes_to_out: {nodes_to_out}')
    p4 = make_a_path(auxgraph, r, depot_cluster,nodes_to_out, thresh)
    #p4 = make_a_path(auxgraph, r, depot_cluster, [p,q], thresh)
    
    return [(*pp, p,q,r) for pp in [p1, p2, p3, p4] if pp != None]




def generate_2path_cuts(model, depot_cluster=1):

    clusters = model._clusters
    tree_closure = model._tree_closure
    
    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)
    epoch = model._epoch
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]

    if NEED_2PATH_CUTS_SAMPLING:
        rates = model._2path_rates
        combs = list(permutations(cluslist,3))
        iter = get_a_sample(epoch,cluslist,combs,rates, FRACTION_2PATH_CUTS)
    else:
        iter = permutations(cluslist,3) 


    aux_G = createClusterGraph(clusters, wrapped_u)
    counter = 0
    
    for p,q,r in iter:
        #print(f'2-path: {p} {q} {r}')
        counter += 1
        if wrapped_y[p,q] + wrapped_y[q,r] > 1 + EPSILON:    
            paths  = create_2path_Constraint(aux_G, p, q, r, wrapped_y, depot_cluster, tree_closure)
            for cut in paths:
                cuts.add(cut)
            if paths and NEED_2PATH_CUTS_SAMPLING:
                #model._2path_rates[(p,q,r)] += RATES_2PATH_CUTS_STEP
                update_rates(model._2path_rates, (p,q,r), RATES_2PATH_CUTS_STEP)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = add_2path_cuts(model,cuts)
    return new_model, number_of_cuts, counter * 4 
    #return cuts, number_of_cuts 

















