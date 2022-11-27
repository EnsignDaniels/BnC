import networkx as nx
from itertools import combinations, permutations
from gurobiHelper import add_3vGDDL_cuts
from parameters import *
import numpy as np

from sampler_v2 import get_a_sample, update_rates, get_sample_from_a_huge_population


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
    
def test(model,S,barS,p,q,r,s):
    u = model._test_u
    y = model._test_y
    
    if sum(u[p,q] for p in S for q in barS) < y[p,q] + y[q,r] + y[r,s] - 1:
        print('ERROR - Constraint violated')
        print(f'S: {S}, BarS: {barS}')


def create_3vGDDL_Constraint(auxgraph,p,q,r,s, wrapped_y, depot_cluster, tree_closure,model):
    
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
            if set(connect_with_source).issubset(S) and set(connect_with_target).issubset(barS):
                #test(model,S,barS,p,q,r,s)
                return tuple(S), tuple(barS)
        return None  
    
    thresh = wrapped_y[p,q] + wrapped_y[q,r] + wrapped_y[r,s] - 1
    
    p_pre = set(tree_closure.predecessors(p))
    q_pre = set(tree_closure.predecessors(q))
    r_pre = set(tree_closure.predecessors(r))
    s_pre = set(tree_closure.predecessors(s))
    
    p_suc = set(tree_closure.successors(p))
    q_suc = set(tree_closure.successors(q))
    r_suc = set(tree_closure.successors(r))
    s_suc = set(tree_closure.successors(s))
    
    nodes_to_out = {depot_cluster} | ((q_pre | r_pre | s_suc) & (p_pre | q_suc)) - {p,q,r,s} 
    #p1 = make_a_path(auxgraph, [p,r], [q,s], [depot_cluster], thresh)
    p1 = make_a_path(auxgraph, [p,r], [q,s], nodes_to_out, thresh)
    
    nodes_to_out = {r} | ((r_pre | s_pre) & (p_pre | q_suc | r_suc)) - {depot_cluster, p, q, s}
    #p2 = make_a_path(auxgraph, [p,s], [q,depot_cluster], [r], thresh)
    p2 = make_a_path(auxgraph, [p,s], [q,depot_cluster], nodes_to_out, thresh)
    
    nodes_to_out = (p_suc | q_suc) & (q_pre | r_pre | s_suc) | {q} - {depot_cluster, p, r, s}
    #p3 = make_a_path(auxgraph, [depot_cluster,r], [p,s], [q], thresh)
    p3 = make_a_path(auxgraph, [depot_cluster,r], [p,s], nodes_to_out, thresh)
        
    
    #return [(*pp, p,q,r,s) for pp in [p1, p2, p3] if pp != None]
    return [(*pp, p,q,r,s) for pp in [p1,p2,p3] if pp != None]



def generate_3vGDDL_cuts(model, depot_cluster=1):
    clusters = model._clusters
    tree_closure = model._tree_closure

    wrapped_u = model.cbGetNodeRel(model._u)
    wrapped_y = model.cbGetNodeRel(model._y)

    epoch = model._epoch
    
    cuts = set()
    cluslist = [c for c in clusters if c != depot_cluster]

    if NEED_3vGDDL_CUTS_SAMPLING:
        #rates = model._3vGDDL_rates
        #combs = list(permutations(cluslist,4))
        #iter = get_a_sample(cluslist,combs,rates, FRACTION_3vGDDL_CUTS)
        iter = get_sample_from_a_huge_population(epoch,cluslist,4,FRACTION_3vGDDL_CUTS)
    else:
        iter = permutations(cluslist,4)

    aux_G = createClusterGraph(clusters, wrapped_u)
    counter = 0
    for p,q,r,s in iter:
        counter += 1
        if wrapped_y[p,q] + wrapped_y[q,r] + wrapped_y[r,s] > 1:
            paths  = create_3vGDDL_Constraint(aux_G, p, q, r, s, wrapped_y, depot_cluster, tree_closure,model)
            for cut in paths:
                cuts.add(cut)
            #if paths and NEED_3vGDDL_CUTS_SAMPLING:
                #model._3vGDDL_rates[(p,q,r,s)] += RATES_3vGDDL_CUTS_STEP
                #update_rates(model._3vGDDL_rates, (p,q,r,s), RATES_3vGDDL_CUTS_STEP)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = add_3vGDDL_cuts(model,cuts)
    return new_model, number_of_cuts, counter*3   

















