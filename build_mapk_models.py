from pysb import Model
from pysb.bng import generate_equations

import itertools
import mapk_model

def build_model(mek_seq_rand, mkp_seq_rand, erk_dimerization, mkp_activation):
    Model()
    mapk_model.mapk_monomers()
    if mek_seq_rand == 'seq':
        mapk_model.mek_phos_erk_seq()
    else:
        mapk_model.mek_phos_erk_random()
    if mkp_seq_rand == 'seq':
        mapk_model.mkp_dephos_erk_seq()
    else:
        mapk_model.mkp_dephos_erk_random()
    if erk_dimerization == 'any':
        mapk_model.erk_dimerize_any()
    elif erk_dimerization == 'uT':
        mapk_model.erk_dimerize_uT()
    elif erk_dimerization == 'phos':
        mapk_model.erk_dimerize_phos()
    if erk_dimerization:
        mapk_model.erk_autophos()
    if mkp_activation:
        mapk_model.erk_activate_mkp()

    mapk_model.mapk_initials()
    mapk_model.mapk_observables()
    model.name = 'mek_%s_mkp_%s_erkauto_%s_mkpact_%s' % \
        (mek_seq_rand, mkp_seq_rand, erk_dimerization, mkp_activation)
    return model

def build_all_models():
    mek_seq_rand = ['seq', 'random']
    mkp_seq_rand = ['seq', 'random']
    erk_dimerization = [False, 'any', 'uT', 'phos']
    mkp_activation = [False, True]

    model_combinations = itertools.product(*(mek_seq_rand, mkp_seq_rand,
                                             erk_dimerization, mkp_activation))

    models = {}
    for mc in model_combinations:
        model = build_model(*mc)
        generate_equations(model)
        models[model.name] = model
    return models
