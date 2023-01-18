EPSILON = 0.01
GOUVEIA_EPSILON = 0.000001
MIP_FOCUS = 3
PI_RATES_STEP = 13
SIGMA_RATES_STEP = 13
PI_SIGMA_RATES_STEP = 10

NEED_MIP_START = True
PARALLEL_HEURISTIC = True

DEPOT_CLUSTER = 1
PI_SAMPLE_SIZE = 5
SIGMA_SAMPLE_SIZE = 5
PI_SIGMA_SAMPLE_SIZE = 5
SAMPLING_THRESHOLD = 3
CUT_OFF_GAP = 0.0
ENDGAME_GAP = -1

USE_GUROBI_HEURISTICS = False

CANCEL_OBSOLETE_CUTS = False
MAXIMUM_NUMBER_OF_CUTS = -1
ROOT_TIME_LIMIT = 5400000
REASON_ROOT_EVAC = 101
REASON_ENDGAME = 102

CAN_LEAVE_ROOT = True
BEST_BOUND_RATE = 1.02

GUROBI_CUTS = False


ENABLE_BALAS_CUTS = True
ENABLE_GBCP_CUTS = True
ENABLE_GOUVEIA_25_29_CUTS = True
ENABLE_GDDL_CUTS = False
ENABLE_SIMPLE_CUTS = False
ENABLE_2PATH_CUTS = False
ENABLE_3vGDDL_CUTS = False
ENABLE_4vGDDL_CUTS = False


##### SAMPLING OF GOUVEIA CUTS
NEED_GDDL_CUTS_SAMPLING = True
NEED_SIMPLE_CUTS_SAMPLING = True
NEED_2PATH_CUTS_SAMPLING = True
NEED_3vGDDL_CUTS_SAMPLING = True
NEED_4vGDDL_CUTS_SAMPLING = True

NEED_GBCP_CUTS_SAMPLING = False

GBCP_MAX_PATH_COUNT = 2

NEED_GOUVEIA25_CUTS_SAMPLING = True
NEED_GOUVEIA26_CUTS_SAMPLING = True
NEED_GOUVEIA27_CUTS_SAMPLING = True
NEED_GOUVEIA28_CUTS_SAMPLING = True
NEED_GOUVEIA29_CUTS_SAMPLING = True


##### SAMPLE FRACTIONS
FRACTION_GDDL_CUTS   = 0.1
FRACTION_SIMPLE_CUTS = 1
FRACTION_2PATH_CUTS  = 0.1
FRACTION_3vGDDL_CUTS = 0.01
FRACTION_4vGDDL_CUTS = 0.0001

FRACTION_GBCP_CUTS = 2

FRACTION_GOUVEIA25_CUTS = 50
FRACTION_GOUVEIA26_CUTS = 50
FRACTION_GOUVEIA27_CUTS = 50
FRACTION_GOUVEIA28_CUTS = 50
FRACTION_GOUVEIA29_CUTS = 50

##### RATINGS
RATES_GDDL_CUTS_STEP   = 3
RATES_SIMPLE_CUTS_STEP = 1
RATES_2PATH_CUTS_STEP  = 5
RATES_3vGDDL_CUTS_STEP = 10
RATES_4vGDDL_CUTS_STEP = 10

RATES_GBCP_CUTS_STEP = 10

RATES_GOUVEIA25_CUTS_STEP = 1
RATES_GOUVEIA26_CUTS_STEP = 1
RATES_GOUVEIA27_CUTS_STEP = 1
RATES_GOUVEIA28_CUTS_STEP = 1
RATES_GOUVEIA29_CUTS_STEP = 1

MAXIMUM_SIZE_OF_RATINGS_GB = 1


def printparameters():
	print(f'Current MIP focus:\t{MIP_FOCUS}')
	print(f'MIP start:\t{NEED_MIP_START}')
	print(f'Pi rates step:\t{PI_RATES_STEP}')
	print(f'Sigma rates step:\t{SIGMA_RATES_STEP}')
	print(f'Pi+Sigma rates step:\t{PI_SIGMA_RATES_STEP}')
	print(f'Pi sample size:\t{PI_SAMPLE_SIZE}')
	print(f'Sigma sample size:\t{SIGMA_SAMPLE_SIZE}')
	print(f'Pi+Sigma sample size:\t{PI_SIGMA_SAMPLE_SIZE}')
	print(f'Sampling threshold:\t{SAMPLING_THRESHOLD}')
	print(f'Gurobi cuts:\t{GUROBI_CUTS}')
	print(f'Root time limit:\t{ROOT_TIME_LIMIT}')
	print(f'Endgame gap:\t{ENDGAME_GAP}')
	print(f'Parallel heuristic worker:\t{PARALLEL_HEURISTIC}')
	print(f'Leave or not Root node:\t{CAN_LEAVE_ROOT}')
	print(f'Balas cuts:\t{ENABLE_BALAS_CUTS}')
	print(f'GDDL cuts:\t{ENABLE_GDDL_CUTS}')
	print(f'Simple cuts:\t{ENABLE_SIMPLE_CUTS}')
	print(f'Two-path cuts:\t{ENABLE_2PATH_CUTS}')
	print(f'3vGDDL cuts:\t{ENABLE_3vGDDL_CUTS}')
	print(f'4vGDDL cuts:\t{ENABLE_4vGDDL_CUTS}')
	print(f'Best bound rate:\t{BEST_BOUND_RATE}') 
	print('Gouveia cuts sampling')
	print(f'GDDL cuts:\t{NEED_GDDL_CUTS_SAMPLING}')
	if NEED_GDDL_CUTS_SAMPLING:
		print(f'GDDL cuts:\tFraction={FRACTION_GDDL_CUTS}%, STEP={RATES_GDDL_CUTS_STEP}')
	print(f'Simple cuts:\t{NEED_SIMPLE_CUTS_SAMPLING}')
	if NEED_SIMPLE_CUTS_SAMPLING:
		print(f'Simple cuts:\tFraction={FRACTION_SIMPLE_CUTS}%, STEP={RATES_SIMPLE_CUTS_STEP}')
	print(f'2path cuts:\t{NEED_2PATH_CUTS_SAMPLING}')
	if NEED_2PATH_CUTS_SAMPLING:
		print(f'2-Path cuts:\tFraction={FRACTION_2PATH_CUTS}%, STEP={RATES_2PATH_CUTS_STEP}')
	print(f'3vGDDL cuts:\t{NEED_3vGDDL_CUTS_SAMPLING}')
	if NEED_3vGDDL_CUTS_SAMPLING:
		print(f'3vGDDL cuts:\tFraction={FRACTION_3vGDDL_CUTS}%, STEP={RATES_3vGDDL_CUTS_STEP}')
	print(f'4vGDDL cuts:\t{NEED_4vGDDL_CUTS_SAMPLING}')
	if NEED_4vGDDL_CUTS_SAMPLING:
		print(f'4vGDDL cuts:\tFraction={FRACTION_4vGDDL_CUTS}%, STEP={RATES_4vGDDL_CUTS_STEP}')
	if NEED_GDDL_CUTS_SAMPLING or NEED_SIMPLE_CUTS_SAMPLING or NEED_2PATH_CUTS_SAMPLING or NEED_3vGDDL_CUTS_SAMPLING or NEED_4vGDDL_CUTS_SAMPLING:
		print(f'Maximum size of each rating:\t{MAXIMUM_SIZE_OF_RATINGS_GB} GB')

