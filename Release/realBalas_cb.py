import networkx as nx
from itertools import combinations
from gurobiHelper import addRealBalas_cuts
from parameters import *


pi_node_dict = {}
sigma_node_dict = {}
pi_sigma_node_dict = {}

global_cuts = set()

def getPiAuxGraphNodes(cluster, clusters, tree_closure, depot_cluster):
    global pi_node_dict
    if not cluster in pi_node_dict:
        pi_node_dict[cluster] = [cl for cl in clusters if not cl in tree_closure.predecessors(cluster)] + [depot_cluster]
    return pi_node_dict[cluster]

def getSigmaAuxGraphNodes(cluster, clusters, tree_closure):
    global sigma_node_dict
    if not cluster in sigma_node_dict:
        sigma_node_dict[cluster] = [cl for cl in clusters if not cl in tree_closure.successors(cluster)]
    return sigma_node_dict[cluster]

def getPiSigmaAuxGraphNodes(cluster1, cluster2, clusters, tree_closure):
    global pi_sigma_node_dict
    if not (cluster1, cluster2) in pi_sigma_node_dict:
        pi_sigma_node_dict[cluster1, cluster2] = [cl for cl in clusters if not cl in tree_closure.predecessors(cluster1) and not cl in tree_closure.successors(cluster2)]
    return pi_sigma_node_dict[cluster1, cluster2]


def createClusterGraphPi(cluster, clusters, tree_closure, wrapped_u, depot_cluster):
    
    cl_star = [cl for cl in clusters if not cl in tree_closure.predecessors(cluster)] + [depot_cluster]
    auxgraph = nx.DiGraph()
    auxgraph.add_nodes_from(cl_star)
    for e in combinations(auxgraph.nodes, 2):
        if e in wrapped_u:
            auxgraph.add_edge(*e, capacity=wrapped_u[e])
        e = e[::-1]
        if e in wrapped_u:
            auxgraph.add_edge(*e, capacity=wrapped_u[e])
    return auxgraph



def createClusterGraphSigma(cluster, clusters, tree_closure, wrapped_u):
    
    cl_star = [cl for cl in clusters if not cl in tree_closure.successors(cluster)]
    auxgraph = nx.DiGraph()
    auxgraph.add_nodes_from(cl_star)
    for e in combinations(auxgraph.nodes, 2):
        if e in wrapped_u:
            auxgraph.add_edge(*e, capacity=wrapped_u[e])
        e = e[::-1]
        if e in wrapped_u:
            auxgraph.add_edge(*e, capacity=wrapped_u[e])
    return auxgraph



def createClusterGraphPiSigma(cluster1, cluster2, clusters, tree_closure, wrapped_u):
    
    cl_star = [cl for cl in clusters if not cl in tree_closure.predecessors(cluster1) and not cl in tree_closure.successors(cluster2)]
    auxgraph = nx.DiGraph()
    auxgraph.add_nodes_from(cl_star)
    for e in combinations(auxgraph.nodes, 2):
        if e in wrapped_u:
            auxgraph.add_edge(*e, capacity=wrapped_u[e])
        e = e[::-1]
        if e in wrapped_u:
            auxgraph.add_edge(*e, capacity=wrapped_u[e])
    return auxgraph
    
    

def createConstraintPi(auxgraph,cluster, tree_closure, ampl, depot_cluster):
    cut_val, part = nx.minimum_cut(auxgraph,cluster,depot_cluster)
    if cut_val < 1 - ampl * EPSILON:
        U,V = part
        Pi_U = set([cl for cl in U if set(tree_closure.successors(cl)).intersection(U)])
        
        cut = tuple(sorted(list(U - Pi_U))), tuple(sorted(list(V - Pi_U)))
        return tuple(U - Pi_U), tuple(V - Pi_U)
    else:
        return None          
    
def createConstraintSigma(auxgraph,cluster, tree_closure, ampl, depot_cluster):
    cut_val, part = nx.minimum_cut(auxgraph,depot_cluster,cluster)
    if cut_val < 1 - ampl * EPSILON:
        U,V = part
        Sigma_V = set([cl for cl in V if set(tree_closure.predecessors(cl)).intersection(V)])
        return tuple(U - Sigma_V), tuple(V - Sigma_V)
    else:
        return None         


def createConstraintPiSigma(auxgraph, cluster1, cluster2, ampl):
    cut_val, part = nx.minimum_cut(auxgraph, cluster1, cluster2)
    if cut_val < 1 - ampl * EPSILON:
        U, V = part
        return tuple(U), tuple(V)
    else:
        return None



        


def generatePi_cuts(model,ampl=1, depot_cluster=1):
    ############################
    clusters = model._clusters
    tree_closure = model._tree_closure
    wrapped_u = model.cbGetNodeRel(model._u)

    
    cuts = set()
    for clus in [cl for cl in clusters if cl != depot_cluster]:
        aux_G = createClusterGraphPi(clus, clusters, tree_closure, wrapped_u, depot_cluster)
        cut = createConstraintPi(aux_G, clus, tree_closure, ampl, depot_cluster)
        if cut:
            cuts.add(cut)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = addRealBalas_cuts(model,cuts)
    return new_model, number_of_cuts



def generateSigma_cuts(model,ampl=1,depot_cluster=1):
    ##########################
    clusters = model._clusters
    tree_closure = model._tree_closure
    wrapped_u = model.cbGetNodeRel(model._u)

    
    cuts = set()
    for clus in [cl for cl in clusters if cl != depot_cluster]:
        aux_G = createClusterGraphSigma(clus, clusters, tree_closure, wrapped_u)
        cut = createConstraintSigma(aux_G, clus, tree_closure, ampl, depot_cluster)
        if cut:
            cuts.add(cut)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = addRealBalas_cuts(model,cuts)
    return new_model, number_of_cuts


def generatePiSigma_cuts(model,ampl=1, depot_cluster=1):
    ######################################
   
    
    clusters = model._clusters
    tree = model._tree
    tree_closure = model._tree_closure
    

    wrapped_u = model.cbGetNodeRel(model._u)
    cuts = set()
    for cluster1, cluster2 in tree.edges:
        if cluster1 == depot_cluster:
            continue
        aux_G = createClusterGraphPiSigma(cluster1, cluster2, clusters, tree_closure, wrapped_u)
        cut = createConstraintPiSigma(aux_G, cluster1, cluster2,ampl)
        if cut:
            cuts.add(cut)
    number_of_cuts = len(cuts)
    new_model = model
    if cuts:
        new_model = addRealBalas_cuts(model,cuts)
    return new_model, number_of_cuts 














