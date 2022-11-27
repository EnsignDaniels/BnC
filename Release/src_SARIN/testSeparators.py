import sys
import networkx as nx
from itertools import permutations, combinations
from math import log10

from fromPCGLNS import getInstance
from preprocess import preprocessInstance
from heuristicHelper import parseResult

from realBalas_cb import generatePi_cuts, generateSigma_cuts, generatePiSigma_cuts
from GBCP_STR_cb import generate_GBCP_cuts
from Gouveia_cuts_cb import init_Gouveia_25_29_rates, generate_25_cuts, generate_26_cuts, generate_27_cuts, generate_28_cuts, generate_29_cuts
from GDDL_STR_cb import generateGDDL_cuts
from simplecut_STR_cb import generateSimple_cuts
from Twopath_STR_cb import generate_2path_cuts

model_name = 'ft53.4'
path_name = 'input/'
PCGLNS_EXT='.pcglns'
DEPOT_CLUSTER = 1


def testPi_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generatePi_cuts(model)
	if number_of_cuts == 0:
		print('Balas Pi    cuts validation test:\t passed')
	else:
		print(f'Balas Pi    cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')

def testSigma_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generateSigma_cuts(model)
	if number_of_cuts == 0:
		print('Balas Sigma cuts validation test:\t passed')
	else:
		print(f'Balas Sigma cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')

def testPiSigma_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generatePiSigma_cuts(model)
	if number_of_cuts == 0:
		print('Balas PiSi  cuts validation test:\t passed')
	else:
		print(f'Balas PiSi  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')

def testGBCP_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generate_GBCP_cuts(model)
	if number_of_cuts == 0:
		print('Balas GBCP  cuts validation test:\t passed')
	else:
		print(f'Balas GBCP  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')

def testGouveia25_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generate_25_cuts(model)
	if number_of_cuts == 0:
		print('Gouveia 25  cuts validation test:\t passed')
	else:
		print(f'Gouveia 25  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')
		u = model._u
		for p,q,r,s in cuts:
			print(f'u[{p},{q}]+u[{q},{p}]+u[{r},{s}]+u[{s},{r}] = {u[(p,q)]+u[(q,p)]+u[(r,s)]+u[(s,r)]}')

def testGouveia26_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generate_26_cuts(model)
	if number_of_cuts == 0:
		print('Gouveia 26  cuts validation test:\t passed')
	else:
		print(f'Gouveia 26  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')

def testGouveia27_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generate_27_cuts(model)
	if number_of_cuts == 0:
		print('Gouveia 27  cuts validation test:\t passed')
	else:
		print(f'Gouveia 27  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')

def testGouveia28_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generate_28_cuts(model)
	if number_of_cuts == 0:
		print('Gouveia 28  cuts validation test:\t passed')
	else:
		print(f'Gouveia 28  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')

def testGouveia29_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generate_29_cuts(model)
	if number_of_cuts == 0:
		print('Gouveia 29  cuts validation test:\t passed')
	else:
		print(f'Gouveia 29  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')	

def testGDDL_STR_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generateGDDL_cuts(model)
	if number_of_cuts == 0:
		print('Strn  GDDL  cuts validation test:\t passed')
	else:
		print(f'Strn  GDDL  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')	

def testSimple_STR_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generateSimple_cuts(model)
	if number_of_cuts == 0:
		print('Strn Simple cuts validation test:\t passed')
	else:
		print(f'Strn Simple cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')	

def test2path_STR_cuts(model, depot_cluster = DEPOT_CLUSTER):
	cuts, number_of_cuts = generate_2path_cuts(model)
	if number_of_cuts == 0:
		print('Strn 2path  cuts validation test:\t passed')
	else:
		print(f'Strn 2path  cuts validation test:\t failed, number of invalid ineqns: {number_of_cuts}')	



class Model:
	_clusters = None
	_tree = None
	_tree_closure = None
	_u = None
	_y = None
	_PiRates = None
	_SigmaRates = None
	_PiSigmaRates = None
	_GBCP_rates = None
	_epoch = None
	_rates25 = None
	_rates26 = None
	_rates27 = None
	_rates28 = None
	_rates29 = None
	_GDDL_rates = None
	_Simple_rates = None
	_2path_rates = None


	def cbGetNodeRel(self,val):
		return val


def main():
	print('Branch-and-Cut separation tester')
	print(f'model name: {model_name}')

	G, clusters, tree = getInstance(f'{path_name}{model_name}{PCGLNS_EXT}')
	G, clusters, tree, tree_closure = preprocessInstance(G,clusters,tree)

	#print(f'clusters: {clusters}')
	#print(f'order_DAG: {tree.edges}')

	n, m, tour, obj = parseResult(model_name)
	print(f'number of nodes: {n}')
	print(f'number of clusters: {m}')
	print(f'\ntour: {tour}')
	print(f'objective: {obj}')

	cluster_tour = []
	for v in tour:
		for c in clusters:
			if v in clusters[c]:
				cluster_tour.append(c)
				break
	print(f'\ncluster_tour: {cluster_tour}')

	u = {(p,q): 0 for (p,q) in permutations(clusters,2)}
	y = {(p,q): 0 for (p,q) in permutations(clusters,2)}

	c_pred = 1
	idx_cur = 1
	for c_cur in cluster_tour[1:]:
		u[(c_pred,c_cur)] = 1
		for c2 in cluster_tour[idx_cur:-1]:
			y[(c_pred,c2)] = 1
		idx_cur += 1
		c_pred = c_cur

	model = Model()
	model._clusters = clusters
	model._tree = tree
	model._tree_closure = tree_closure
	model._u = u
	model._y = y	


	##### init GBCP tests #######	
	T = tree.copy()
	T.remove_node(DEPOT_CLUSTER)
	
	roots =  [v for v,d in T.in_degree()  if d == 0]
	leaves = [v for v,d in T.out_degree() if d == 0]
	ap_rates = {(r,l):1 for r in roots for l in leaves if nx.has_path(T,r,l)}
			
	model._GBCP_rates = ap_rates
	##############################

	###### init Gouveia tests ####
	init_Gouveia_25_29_rates(model)
	##############################

	
	#print(f'\nu reconstructed\n:{[(p,q) for p in clusters for q in clusters if u[(p,q)] == 1 ]}')
	#print(f'\ny reconstructed\n:{[(p,q) for p in clusters for q in clusters if y[(p,q)] == 1 ]}')


	print('Test initialized successfully')
	testPi_cuts(model)
	testSigma_cuts(model)
	testPiSigma_cuts(model)
	testGBCP_cuts(model)
	testGouveia25_cuts(model)
	testGouveia26_cuts(model)
	testGouveia27_cuts(model)
	testGouveia28_cuts(model)
	testGouveia29_cuts(model)
	#testSimple_STR_cuts(model)
	#testGDDL_STR_cuts(model)
	#test2path_STR_cuts(model)


	# print(f'u:\n{u}')



if __name__ == '__main__':
	try:
		for arg in sys.argv:
			if '=' in arg:
				parts = arg.split('=')
				if parts[0] == '--input' or parts[0] == '-i':
					model_name = parts[1]
				if parts[0] == '--path' or parts[0] == '-p':
					path_name = f'{parts[1]}/'
	except:
		print('SYNTAX: python testSeparators.py -i=<model_name> -p=<path_name> ')

	main()