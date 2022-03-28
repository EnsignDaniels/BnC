import numpy as np
import networkx as nx
import gurobipy as gp
from gurobipy import GRB, quicksum
from itertools import permutations, combinations


def createExactModel(model_name, G, clusters, tree,  mips = None, first_cluster = 1):
    with gp.Env(empty = True) as env:
        env.setParam('LogToConsole', 1)
        env.start()
        
        model=gp.Model(model_name, env=env)
        
        n = len(G)
        n_list = list(G)
        n_dict = {n_list[idx]: idx for idx in range(n)}
        m = len(clusters)
        
        tree_closure = nx.transitive_closure_dag(tree)
        edgedict = {t[:2]: t[2] for t in G.edges(data='weight')}
 

        clust_pairs = {(c1,c2):None for c1,c2 in permutations(clusters,2)}
        x_ij, cost = gp.multidict(edgedict)
        z_i,_ = gp.multidict(n_dict)
        u_pq,_ = gp.multidict(clust_pairs)
        y_pq,_ = gp.multidict(clust_pairs)
 
        ### VARIABLES
        x = model.addVars(x_ij, vtype=GRB.BINARY, name='x')
        z = model.addVars(z_i, vtype=GRB.BINARY, name='z')
        y = model.addVars(y_pq, name='y')
        u = model.addVars(u_pq, name='u')

               
        ### OBJECTIVE
        objective = model.setObjective(x.prod(cost), GRB.MINIMIZE)

        ### CONSTRAINTS
        ### Constraint 2
        for c in clusters:
            model.addConstr(sum(z[i] for i in clusters[c]) == 1, f'(2_{c})')

        ### Constraint 3
        for i in G:
            model.addConstr(sum(x[i,j] for j in G.successors(i)) == z[i], f'(3_{i})')

        ### Constraint 4
        for i in G:
            model.addConstr(sum(x[j,i] for j in G.predecessors(i)) == z[i], f'(4_{i})')

        ### Constraint 5 (1+2)
        for p in clusters:
            model.addConstr(sum(u[p,q] for q in clusters if q != p) == 1, f'(5_1_{p})')
            model.addConstr(sum(u[q,p] for q in clusters if q != p) == 1, f'(5_2_{p})')

        ### Constraint 6

        for p,q in permutations(clusters,2):
            model.addConstr(u[p,q] == sum(x[i,j] for i in clusters[p] for j in clusters[q] if G.has_edge(i,j)), f'(6_{p}_{q})')

        ### Constraint 7

        clust_keys = list(clusters.keys())
        for p_idx in range(m):
            p = clust_keys[p_idx] 
            if p == first_cluster:
                continue
            for q_idx in range(p_idx+1,m):
                q = clust_keys[q_idx]
                if q == first_cluster:
                    continue
                for r_idx in range(p_idx+1,m):
                    r = clust_keys[r_idx]
                    if r_idx == q_idx or r == first_cluster:
                        continue
                    model.addConstr(y[p,q] + y[q,r] + y[r,p] <= 2, f'(7_{p}_{q}_{r})')

        ### Constraint 8

        for p, q in permutations(clusters, 2):
            if first_cluster in (p,q):
                continue
            model.addConstr(u[p,q] - y[p,q] <= 0, f'(8_{p}_{q})')

        ### Constraint 9

        for p, q in combinations(clusters, 2):
            if first_cluster in (p,q):
                continue
            model.addConstr(y[p,q] + y[q,p] == 1, f'(9_{p}_{q})')


        ### Constraint 10

        for p, q in permutations(clusters, 2):
            if first_cluster in (p,q):
                continue
            if tree_closure.has_edge(p,q):
                model.addConstr(y[p,q] == 1, f'(10_{p}_{q})')

        ### MIP start
        if mips:
            for key in x:
                x[key].start = 0
            for key in z:
                z[key].start = 0

            for key, val in mips['x'].items():
                x[key].start = val
            for key, val in mips['z'].items():
                z[key].start = val

        ### incoporate variables into the model
        model._x = x
        model._z = z
        model._u = u
        model._y = y


    return model, x, y, u, z 



def optimizeModel(model, time_limit, threads, callback = None):
    model.setParam(GRB.Param.TimeLimit,time_limit)
    model.setParam(GRB.Param.Threads, threads)
    model.setParam(GRB.Param.MIPFocus, 2)

    if callback:
        model.setParam(GRB.Param.Cuts, 0)
        model.setParam(GRB.Param.PreCrush, 1)
        model.optimize(callback)
    else:
        model.optimize()
    
    assert model.status == GRB.OPTIMAL, f'Optimum value has not been found, status code:{model.status}'
    return model.status

def addGSEC_cuts(model,cut_counter):
    u = model._u
    for U,V in cut_counter:
        model.cbCut(sum(u[p,q] for p in U for q in V) >= 1)
    return model

def addBalas_cuts(model,cut_counter, Pi_or_Sigma):
    G = model._G
    x = model._x
    z = model._z
    for U,V, node in cut_counter:
        model.cbCut(sum(x[i,j] for i in U for j in V if G.has_edge(i,j)) - z[node] >= 0)
    return model


def addBalasPS_cuts(model,cuts):
    G = model._G
    x = model._x
    z = model._z
    for u,v in cuts:
        I,J = cuts[(u,v)]
        model.cbCut(sum(x[i,j] for i in I for j in J if G.has_edge(i,j)) - z[u] - z[v] >= -1)
    return model


def addBalasPS_cluster_cuts(model,cuts):
    G = model._G
    x = model._x
    for u,v in cuts:
        I,J = cuts[(u,v)]
        model.cbCut(sum(x[i,j] for i in I for j in J if G.has_edge(i,j)) >= 1)
    return model

def addBalasPSfull_cluster_cuts(model,cuts):
    u = model._u
    for c1,c2 in cuts:
        P,Q = cuts[(c1,c2)]
        model.cbCut(sum(u[p,q] for p in P for q in Q if (p,q) in u) >= 1)
    return model

def addRealBalas_cuts(model, cuts):
    u = model._u
    for P,Q in cuts:
        model.cbCut(sum(u[p,q] for p in P for q in Q) >= 1)
    return model






def main(model_name, G, clusters, tree, tree_closure):
    mips = {}
    mips['x'], mips['z'] = getMIPS(model_name)
    model, x, y, u, z = createExactModel(f'{model_name}-exact', G, clusters, tree, mips)


    model.write(f'models/{model_name}.lp')

    status = optimizeModel(model)
    print(status)



if __name__ == '__main__':
    from fromPCGLNS import getInstance
    from preprocess import preprocessInstance
    from heuristicHelper import getMIPS
    

    model_name='ESC12'
    G, clusters, tree = getInstance(f'PCGLNS/PCGLNS_PCGTSP/{model_name}.pcglns')
    G, clusters, tree, tree_closure = preprocessInstance(G,clusters,tree)


    main(model_name, G, clusters, tree, tree_closure)

