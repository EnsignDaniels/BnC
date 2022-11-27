EPSILON = 0.01
MIP_FOCUS = 2
PI_RATES_STEP = 1
SIGMA_RATES_STEP = 1
PI_SIGMA_RATES_STEP = 1

DEPOT_CLUSTER = 1
PI_SAMPLE_SIZE = 3
SIGMA_SAMPLE_SIZE = 3
PI_SIGMA_SAMPLE_SIZE = 2
SAMPLING_THRESHOLD = 500
CUT_OFF_GAP = 0.05 

GUROBI_CUTS = False
CANCEL_OBSOLETE_CUTS = False
MAXIMUM_NUMBER_OF_CUTS = 40

ROOT_TIME_LIMIT = 3600
REASON_ROOT_EVAC = 101
REASON_ENDGAME = 102 
ENDGAME_GAP = 0.0

def printparameters():
	print(f'Current MIP focus:\t{MIP_FOCUS}')
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