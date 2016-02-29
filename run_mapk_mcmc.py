import sys
import csv
import emcee
import numpy
import matplotlib.pyplot as plt
import mapk_model

from pysb.integrate import Solver
from pysb import *

class MapkExperiment(object):
    def __init__(self):
        self.ts = []
        self.ERK = []
        self.ERKpY = []
        self.ERKpT = []
        self.ERKpTpY = []
        self.ERK_std = []
        self.ERKpY_std = []
        self.ERKpT_std = []
        self.ERKpTpY_std = []
        self.ERKtot = 0
        self.MEKtot = 0
        self.MKPtot = 0
    
    def process_data(self):
        self.ERKpTpY = numpy.array(self.ERKpTpY)
        self.ERKpT = numpy.array(self.ERKpT)
        self.ERKpY = numpy.array(self.ERKpY)
        self.ERK = numpy.array(self.ERK)
     
        sigma_min = 0.05
        self.ERK_std = numpy.array([d if d > sigma_min else 
                                    sigma_min for d in self.ERK_std])
        self.ERKpY_std = numpy.array([d if d > sigma_min else 
                                      sigma_min for d in self.ERKpY_std])
        self.ERKpT_std = numpy.array([d if d > sigma_min else 
                                      sigma_min for d in self.ERKpT_std])
        self.ERKpTpY_std = numpy.array([d if d > sigma_min else 
                                        sigma_min for d in self.ERKpTpY_std])

def getfloat(s):
    if s == '':
        return numpy.nan
    else:
        return float(s)

def read_data():
    nexp = 9
    nt = 23
    ns = 4
    nrep = 9

    initial_conditions = [[930, 0, 0], [930, 0, 0], [930, 0, 93], 
                          [930, 93, 93], [930, 186, 93], [930, 465, 93],
                          [930, 93, 93], [930, 93, 0], [465, 93, 93]]

    data = []

    with open('../Data/combined_data_9.csv', 'r') as fh:
        lines = fh.readlines()
        blocksize = 1 + nt
        for i in range(nexp):
            lines_exp = lines[i*blocksize : (i+1)*blocksize]
            csv_reader = csv.DictReader(lines_exp)
            exp = MapkExperiment()
            for line in csv_reader:
                erk_data = numpy.array([getfloat(line['A%d' % (r + 1)]) for 
                                        r in range(nrep)])
                erk_mean = numpy.nanmean(erk_data)
                # Assumption: all values will be nan so it's enough to test ERK
                if not numpy.isnan(erk_mean):
                    exp.ERK.append(erk_mean)
                    exp.ERK_std.append(numpy.nanstd(erk_data))
                    erkpy_data = numpy.array([getfloat(line['B%d' % (r + 1)]) 
                                              for r in range(nrep)])
                    erkpt_data = numpy.array([getfloat(line['C%d' % (r + 1)]) 
                                              for r in range(nrep)])
                    erkptpy_data = numpy.array([getfloat(line['D%d' % (r + 1)]) 
                                                for r in range(nrep)])
                    erkptpy_data -= erkpt_data
                    erkpt_data -= erkpy_data
                    erkpy_data -= erk_data
                    exp.ERKpY.append(numpy.nanmean(erkpy_data))
                    exp.ERKpY_std.append(numpy.nanstd(erkpy_data))
                    exp.ERKpT.append(numpy.nanmean(erkpt_data))
                    exp.ERKpT_std.append(numpy.nanstd(erkpt_data))
                    exp.ERKpTpY.append(numpy.nanmean(erkptpy_data))
                    exp.ERKpTpY_std.append(numpy.nanstd(erkptpy_data))
                    exp.ts.append(getfloat(line['R%d' % (i + 1)]) * 60.0)
            
            exp.ERKtot = initial_conditions[i][0]
            exp.MEKtot = initial_conditions[i][1]
            exp.MKPtot = initial_conditions[i][2]
            exp.process_data()
            data.append(exp)
    return data

def sim_experiment(model, exp, pd=None):
    if pd is None:
        pd = {}
    pd['ERK_0'] = exp.ERKtot * exp.ERK[0]
    pd['ERKpT_0'] = exp.ERKtot * exp.ERKpT[0]
    pd['ERKpY_0'] = exp.ERKtot * exp.ERKpY[0]
    pd['ERKpTpY_0'] = exp.ERKtot * exp.ERKpTpY[0]
    pd['MEK_0'] = exp.MEKtot
    pd['MKP_0'] = exp.MKPtot

    solver = Solver(model, exp.ts)
    solver.run(pd)
    return solver.yobs

def gauss_lh(x, mu, sigma):
    return -numpy.sum((x-mu)**2 / sigma**2)

def parameter_dict(model, p):
    pd = {}
    pe = [pp for pp in model.parameters if pp.name[0]=='k']
    for pnew, pold in zip(p, pe):
        pd[pold.name] = numpy.power(10.0, pnew)
    return pd

def plot_fit(model, data, pd=None):
    fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True)
    for i, exp in enumerate(data):
        if i == 1:
            continue
        ax = axs[int(numpy.floor(i/3)), i%3]
        yobs = sim_experiment(model, exp, pd)
        erkpt = yobs['ERKpT'] / exp.ERKtot
        erkpy = yobs['ERKpY'] / exp.ERKtot
        erkptpy = yobs['ERKpTpY'] / exp.ERKtot
        
        erk0 = exp.ERKtot * exp.ERK[0]
        erkpt0 = exp.ERKtot * exp.ERKpT[0]
        erkpy0 = exp.ERKtot * exp.ERKpY[0]
        erkptpy0 = exp.ERKtot * exp.ERKpTpY[0]
        mek0 = exp.MEKtot
        mkp0 = exp.MKPtot

        ax.set_title('R%d: [(%d, %d, %d, %d), %d, %d]' % 
            (i+1, erk0, erkpt0, erkpy0, erkptpy0, mek0, mkp0))

        ax.set_xlim([0, 19800])
        ax.set_ylim([0, 1])

        ax.errorbar(exp.ts, exp.ERKpY, yerr=exp.ERKpY_std, 
                    fmt='bo', label='ERKpY')
        ax.errorbar(exp.ts, exp.ERKpT, yerr=exp.ERKpT_std, 
                    fmt='go', label='ERKpT')
        ax.errorbar(exp.ts, exp.ERKpTpY, yerr=exp.ERKpTpY_std, 
                    fmt='ro', label='ERKpTpY')

        ax.plot(exp.ts, erkpy, 'b')
        ax.plot(exp.ts, erkpt, 'g')
        ax.plot(exp.ts, erkptpy, 'r')
    fig.show()

def plot_best(model, data, sampler):
    idx = numpy.argmax(sampler.flatlnprobability)
    p = sampler.flatchain[idx]
    pd = parameter_dict(model, p)
    plot_fit(model, data, pd)

def likelihood(p, model, data, plot=False):
    pd = parameter_dict(model, p)
    lh = 0
    for i, exp in enumerate(data):
        if i == 1:
            continue
        yobs = sim_experiment(model, exp, pd)
        erk = yobs['ERKu'] / exp.ERKtot
        erkpt = yobs['ERKpT'] / exp.ERKtot
        erkpy = yobs['ERKpY'] / exp.ERKtot
        
        lh += gauss_lh(erk[1:], exp.ERK[1:], exp.ERK_std[1:])
        lh += gauss_lh(erkpt[1:], exp.ERKpT[1:], exp.ERKpT_std[1:])
        lh += gauss_lh(erkpy[1:], exp.ERKpY[1:], exp.ERKpY_std[1:])
   
    print lh
    if numpy.isnan(lh):
        return -numpy.inf
    else:
        return lh

def prior(p, model):
    pd = parameter_dict(model, p)
    lp = 0
    for pn, pv in pd.iteritems():
        lp += -(numpy.log10(pv) - numpy.log10(model.parameters[pn].value))**2 / 4.0
    print lp
    return lp
    
def posterior(p, model, data):
    lpri = prior(p, model)
    llh = likelihood(p, model, data)
    lp = lpri + llh
    print lpri, llh
    return lp

def build_markevich_2step():
    Model()
    mapk_model.mapk_monomers()
    mapk_model.mek_phos_erk_2_step()
    mapk_model.mkp_dephos_erk_2_step()
    mapk_model.mapk_initials()
    mapk_model.mapk_observables()
    return model

def build_erk_autophos_any():
    Model()
    mapk_model.mapk_monomers()
    mapk_model.mek_phos_erk_2_step()
    mapk_model.mkp_dephos_erk_2_step()
    mapk_model.erk_dimerize_any()
    mapk_model.erk_autophos()
    mapk_model.mapk_initials()
    mapk_model.mapk_observables()
    return model
    
def build_erk_autophos_uT():
    Model()
    mapk_model.mapk_monomers()
    mapk_model.mek_phos_erk_2_step()
    mapk_model.mkp_dephos_erk_2_step()
    mapk_model.erk_dimerize_uT()
    mapk_model.erk_autophos()
    mapk_model.mapk_initials()
    mapk_model.mapk_observables()
    return model

def build_erk_autophos_phos():
    Model()
    mapk_model.mapk_monomers()
    mapk_model.mek_phos_erk_2_step()
    mapk_model.mkp_dephos_erk_2_step()
    mapk_model.erk_dimerize_uT()
    mapk_model.erk_autophos()
    mapk_model.mapk_initials()
    mapk_model.mapk_observables()
    return model

def build_erk_activate_mkp():
    Model()
    mapk_model.mapk_monomers()
    mapk_model.mek_phos_erk_2_step()
    mapk_model.mkp_dephos_erk_2_step()
    mapk_model.erk_dimerize_uT()
    mapk_model.erk_autophos()
    mapk_model.erk_activate_mkp()
    mapk_model.mapk_initials()
    mapk_model.mapk_observables()
    return model 

def run_one_model(model, data, ns):
    # Vector of nominal parameters
    p = numpy.log10(numpy.array([pp.value for pp in model.parameters 
                                 if pp.name[0]=='k']))
    print posterior(p, model, data)
    
    # Number of temperatures, dimensions and walkers
    ntemps = 20
    ndim = len(p)
    nwalkers = (ndim+1)*2
    
    # Instantiate sampler
    # sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior,
    #                                 args=[model, data])
    
    sampler = emcee.PTSampler(ntemps, nwalkers, ndim, likelihood, prior, 
         threads=1, pool=None, betas=None, a=2.0, Tmax=None, 
         loglargs=[model, data], logpargs=[model], 
         loglkwargs={}, logpkwargs={})
    
    # Random initial parameters for walkers
    p0 = numpy.ones((ntemps, nwalkers, ndim))
    for i in range(ntemps):
        for j in range(nwalkers):
            p0[i, j, :] = p + 1.0*numpy.random.rand(ndim)
    # Run sampler
    sampler.run_mcmc(p0, ns)
    return sampler

if __name__ == '__main__':
    ns = 1000
    if len(sys.argv) > 1:
        model_number = int(sys.argv[1])
        if len(sys.argv) > 2:
            ns = int(sys.argv[2])
    # Read experimental data
    data = read_data()
    # Build model of interest
    if model_number == 1:
        model = build_markevich_2step()
    elif model_number == 2:
        model = build_erk_autophos_phos()
    elif model_number == 3:
        model = build_erk_autophos_uT()
    elif model_number == 4:
        model = build_erk_activate_mkp()

    sampler = run_one_model(model, data, ns)
