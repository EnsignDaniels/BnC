import numpy as np
from itertools import permutations

SEED = 1701


def init_Gouveia_cuts_ratings(model, depot_cluster = 1):
	clusters = model._clusters
	clus = [c for c in clusters if c != depot_cluster]
	model._GDDL_rates = {c: 1 for c in permutations(clus,3)}
	model._Simple_rates = {c: 1 for c in permutations(clus,2)}
	model._2Path_rates = {c: 1 for c in permutations(clus,3)}
	model._3vGDDL_rates = {c: 1 for c in permutations(clus,4)}
	model._4vGDDL_rates = {c: 1 for c in permutations(clus,5)}
	

###############
# Function samples a fraction subset from the set of ordered combinations of the given set
# according to the probability measure given by the dictionary rates
#
# rates
#		key: a combination, e.g. (p,q,r)
#		val: rating of the key


def get_a_sample(groundset, combs, rates, fraction):
	np.random.seed(SEED)
	rates_array = np.array([rates[e] for e in combs])
	prob = rates_array / np.sum(rates_array)
	L = len(rates)
	gs_index = list(range(L))

	sample_size = round(L * fraction / 100)

	sample_index=np.random.choice(gs_index, sample_size, replace=False, p=prob)
	return [combs[i] for i in sample_index]




def test_list():
	groundset=[1,2,3,4,5,6,7,8,9]
	comb_len = 3
	combs = [c for c in permutations(groundset,comb_len)]
	L = len(combs)
	rates = {c: 1 for c in combs}
	fraction = 1
	sample = get_a_sample(groundset, combs, rates, fraction)

	print(f'groundset = {groundset}')
	print(f'comb_len = {comb_len}')
	print(f'total number of combinations {L}')
	print(f'fraction = {fraction}')
	print(f'sample = {sample}')


def main():
	test_list()
	# test_compinations()



if __name__ == '__main__':
	main()