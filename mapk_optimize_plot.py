import os
from run_mapk_mcmc import *
import numpy as np
import pickle
from scipy.stats.distributions import norm
from matplotlib import pyplot as plt

#list the algorithms
method_list = ['Nelder-Mead', 'Powell', 'COBYLA', 'TNC']

#list the models
model_list = [build_markevich_2step, build_erk_autophos_any,
	      build_erk_autophos_uT, build_erk_autophos_phos,
	      build_erk_activate_mkp]

#no. of initial values sampled using Latin Hypercube Sampling
num_ini = 1000
from_id = 0
to_id = 1000

#read the data
data = read_data()

obj_func = np.ones((len(method_list),len(model_list),num_ini))
obj_func[:] = np.nan
func_eval = np.ones((len(method_list),len(model_list),num_ini))
func_eval[:] = np.nan
ave_obj_func = np.ones((len(method_list),len(model_list)))
ave_func_eval = np.ones((len(method_list),len(model_list)))

run_ind = np.arange(num_ini)

colors = ["red","blue","green","yellow","cyan","darkviolet","sienna","darkgray","black"]
mark = ["o","^","x","s","D","H","+","."]
plt.ion()

#read the pkl files and store the result
for i in range(len(method_list)):
	for j in range(len(model_list)):
		model = model_list[j]()
		for k in range(0,1000,10):
			fname = "output/%s_%s_%d_%d.pkl" % (i,j,k,k+10)
			if not os.path.exists(fname):
				print fname + ' does not exist.'
				continue
			results = pickle.load(open(fname, "rb"))
			for l in range(10):
				result = results[l]
				if not result.success:
					continue
				obj_func[i][j][k+l] = result.fun
				func_eval[i][j][k+l] = result.nfev
		ave_obj_func[i][j] = np.nanmean(obj_func[i][j])
		ave_func_eval[i][j] = np.nanmean(func_eval[i][j])
		

#for each algorithm plot the results for different models
for a in range(len(method_list)):
	plt.figure()
	plt.suptitle(method_list[a], fontsize=14, fontweight='bold')
	for b in range(len(model_list)):
		plt.plot(run_ind+1, np.sort(obj_func[a][b]), linestyle='',
			marker='o', color = colors[b],
			label = model_list[b].__name__)
		plt.legend(loc='lower right')
	plt.xlabel("Run index")
	plt.ylabel("Objective function value")
	plt.ylim(-1e5,0)
	plt.grid(True)
		
#for each model plot the results for different algorithms
for aa in range(len(model_list)):
	plt.figure()
	plt.suptitle(model_list[aa].__name__, fontsize=14, fontweight='bold')
	for bb in range(len(method_list)):
		plt.plot(run_ind+1, np.sort(obj_func[bb][aa]), linestyle='', 
			marker='o', color = colors[bb],
			label = method_list[bb])
		plt.legend(loc='lower right')
	plt.xlabel("Run index")
	plt.ylabel("Objective function value")
	plt.ylim(-1e5,0)
	plt.grid(True)


#plot heatmap for objective function values
column_labels = model_list.__name__
row_labels = method_list
data = ave_obj_func
fig, ax = plt.subplots()
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

ax.set_xticklabels(row_labels, minor=False)
ax.set_yticklabels(column_labels, minor=False)
plt.show()

"""To avoid confusion"""
"""
#plot heatmap for no. of function evaluations
column_labels = method_list
row_labels = model_list.__name__
data = ave_func_eval
fig, ax = plt.subplots()
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

ax.set_xticklabels(row_labels, minor=False)
ax.set_yticklabels(column_labels, minor=False)
plt.show()

"""
#objective function vs function evaluations
plt.figure()
for aaa in range(len(method_list)):
	for bbb in range(len(model_list)):
		plt.plot(ave_obj_func[aaa][bbb], ave_func_eval[aaa][bbb], linestyle='', 
			marker=mark[aaa], color=colors[bbb],label=model_list[bbb].__name__)
		plt.legend()
		plt.xlabel("Objective function value")
		plt.ylabel("Calls to the objective function")
		plt.grid(True)






	
				

				
