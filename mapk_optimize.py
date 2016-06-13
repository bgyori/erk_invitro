from run_mapk_mcmc import *
import numpy as np
import scipy.optimize
from pyDOE import *
from scipy.stats.distributions import norm
import sys
import pickle

model_funs = [build_markevich_2step, build_erk_autophos_any,
	      build_erk_autophos_uT, build_erk_autophos_phos,
	      build_erk_activate_mkp]

method_names = ['Nelder-Mead', 'Powell', 'COBYLA','TNC']

def generate_init(model, ns, output_file, lognorm=False):
	#select the parameters
	p_to_fit = [p for p in model.parameters if p.name[0] == 'k']
	nominal_values = np.array([q.value for q in p_to_fit])
	log_nominal_values = np.log10(nominal_values)


	#latin hypercube sampling for picking a starting point
	num_ind = len(p_to_fit)*ns
	ini_val = lhs(len(p_to_fit), samples = num_ind/len(p_to_fit))
	means = log_nominal_values
	stdvs = np.ones(len(p_to_fit))
	if lognorm:
		for ind in range(len(p_to_fit)):
			ini_val[:,ind] = norm(loc=means[ind], scale=stdvs[ind]).ppf(ini_val[:,ind])
	else:
		# First shift unit hypercube to be centered around origin
		# Then scale hypercube along each dimension by 2 stdevs
		# Finally shift hypercube to be centered around nominal values
		ini_val = means + 2 * stdvs * (ini_val - 0.5)
		
	np.save(output_file, ini_val)

def generate_all_inits(ns):
	for mf in model_funs:
		model = mf()
		fname = mf.__name__ + '_init'
		generate_init(model, ns, fname)

def neg_likelihood(*args, **kwargs):
	return -1.0*likelihood(*args, **kwargs)

if __name__ == '__main__':
	if len(sys.argv) < 5:
		print 'Not enough input arguments.'
		sys.exit()
	from_idx = int(sys.argv[1])
	to_idx = int(sys.argv[2])
	if to_idx < from_idx:
		print 'Invalid from-to pair.'
		sys.exit()
	model_id = int(sys.argv[3])
	if model_id >= len(model_funs):
		print 'Invalid model id.'
		sys.exit()
	method_id = int(sys.argv[4])
	if method_id >= len(method_names):
		print 'Invalid method id.'
		sys.exit()
	#read data from file
	data = read_data()

	model = model_funs[model_id]()
	method = method_names[method_id]

	ini_val = np.load(model_funs[model_id].__name__ + '_init.npy')

	results = []
	for i in range(from_idx, to_idx):
		result = scipy.optimize.minimize(neg_likelihood, ini_val[i],
						 args=(model, data),
					 	 method=method)
		results.append(result)
	fname = '%s_%s_%d_%d.pkl' % (method_id, model_id, from_idx, to_idx)
	with open(fname, 'wb') as fh:
		pickle.dump(results, fh)
