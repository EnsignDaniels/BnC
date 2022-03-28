import sys
import networkx as nx
from itertools import permutations, combinations
from math import log10

from fromPCGLNS import getInstance
from preprocess import preprocessInstance
from heuristicHelper import getMIPS, runHeuristic

from gurobipy import GRB

from gurobiHelper import createExactModel, optimizeModel

from heuristicHelper import getUB, getMIPS
from realBalas_cb import generatePi_cuts, generateSigma_cuts, generatePiSigma_cuts




model_name = 'ESC25'
path_name = 'PCGLNS/PCGLNS_PCGTSP/'
PCGLNS_EXT='.pcglns'
MODELS_PATH='models/'


DEPOT_CLUSTER = 1
use_callbacks = True
time_limit = 36000
threads = 3
	

def cutGeneratorCallback(model, where):
    if where == GRB.Callback.MIPNODE:
        nodecnt = int(model.cbGet(GRB.callback.MIPNODE_NODCNT))
 
        status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
        if status == GRB.OPTIMAL:
        	amplificator = 1 / log10(nodecnt + 10)
        	
        	number_of_BalasPi_cuts, number_of_BalasSigma_cuts = 0, 0
        	_, number_of_BalasPi_cuts = generatePi_cuts(model,amplificator, DEPOT_CLUSTER)
        	_, number_of_BalasSigma_cuts = generateSigma_cuts(model,amplificator, DEPOT_CLUSTER)
        	_, number_of_BalasPS_cuts = generatePiSigma_cuts(model, amplificator, DEPOT_CLUSTER)



def main(need_callbacks = True):
	print('Branch-and-Cut PCGTSP solver')
	print(f'Task name:\t{model_name}')
	print(f'Time limit:\t{time_limit} sec')
	print(f'Threads:\t{threads}')

	G, clusters, tree = getInstance(f'{path_name}{model_name}{PCGLNS_EXT}')
	G, clusters, tree, tree_closure = preprocessInstance(G,clusters,tree)

	runHeuristic(model_name, path_name)
	mips = {}
	mips['x'], mips['z'] = getMIPS(model_name)

	model, x, y, u, z = createExactModel(f'{model_name}-exact', G, clusters, tree, mips)

############################	
	model._G = G
	model._clusters = clusters
	model._tree = tree
	model._tree_closure = tree_closure
############################
	
	print('Model is ready\nPCGLNS MIP start incorporated')
	model.write(f'{MODELS_PATH}{model_name}_exact.lp')

	try:
		if need_callbacks:
			status = optimizeModel(model, time_limit, threads, cutGeneratorCallback)
		else:
			status = optimizeModel(model, time_limit, threads,)
	except AssertionError as msg:
		print(msg)




if __name__ == '__main__':
	try:
		for arg in sys.argv:
			if arg == '--without_callbacks' or arg == '-c':
				use_callbacks = False
			if '=' in arg:
				parts = arg.split('=')
				if parts[0] == '--input' or parts[0] == '-i':
					model_name = parts[1]
				if parts[0] == '--path' or parts[0] == '-p':
					path_name = f'{parts[1]}/'
				if parts[0] == '--time_limit' or parts[0] == '-t':
					time_limit = int(parts[1])
				if parts[0] == '--threads' or parts[0] == '-th':
					threads = int(parts[1])
	except:
		print('SYNTAX: python BnC_cb.py -i=<model_name> [-c or --without_callbacks] [-t or -time_limit=<in secs>] [-th or --threads=<# of threads>] [-p or --path=<path to GLNS soruces>] ')

	main(use_callbacks)


