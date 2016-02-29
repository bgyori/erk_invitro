from pysb import *
from pysb.util import alias_model_components

Model()

def mapk_monomers():
    Monomer('ERK', ['b', 'T', 'Y'], {'T': ['u', 'p'], 'Y': ['u', 'p']})
    Monomer('MEK', ['b'])
    Monomer('MKP', ['b'])
    alias_model_components()

# Questions
# * Is the phosphorylation and dephosphorylation distributive or processive
# ** If distributive then is the same binding site used for both T and Y phosphorylation?
# * Do MEK and MKP compete for the same binding site on ERK?

def mek_phos_erk_2_step():
    """Based on Table S4 of Markevich et al. Assuming distributive random 
    kinase mechanism.
    
    First and second order rate constants are in s^{-1} and nM^{-1}s^{-1}
    respectively.
    """
    Parameter('kf_mek', 4e-2)     # = k_1 + k_5
    Parameter('kr_mek', 1)        # = k_{-1} = k_{-3}
    Parameter('kf_mek_p', 3.2e-2) # = k_3 = k_7
    Parameter('kr_mek_p', 1)      # = k_{-5} = k_{-7}
    Parameter('kc_mek', 1e-2)     # = k_2 = k_6
    Parameter('kc_mek_p', 15)     # = k_4 = k_8
    alias_model_components()

    Rule('MEK_bind_ERK_uu', MEK(b=None) + ERK(b=None, T='u', Y='u') <>\
                            MEK(b=1) % ERK(b=1, T='u', Y='u'), kf_mek, kr_mek)
    Rule('MEK_bind_ERK_pT', MEK(b=None) + ERK(b=None, T='p', Y='u') <>\
                            MEK(b=1) % ERK(b=1, T='p', Y='u'), kf_mek_p, kr_mek_p)
    Rule('MEK_bind_ERK_pY', MEK(b=None) + ERK(b=None, T='u', Y='p') <>\
                            MEK(b=1) % ERK(b=1, T='u', Y='p'), kf_mek_p, kr_mek_p)
    Rule('MEK_bind_ERK_pp', MEK(b=None) + ERK(b=None, T='p', Y='p') <>\
                            MEK(b=1) % ERK(b=1, T='p', Y='p'), kf_mek_p, kr_mek_p)
    Rule('MEK_phos_ERK_T', MEK(b=1) % ERK(b=1, T='u', Y='u') >>\
                            MEK(b=None) + ERK(b=None, T='p', Y='u'), kc_mek)
    Rule('MEK_phos_ERK_Y', MEK(b=1) % ERK(b=1, T='u', Y='u') >>\
                            MEK(b=None) + ERK(b=None, T='u', Y='p'), kc_mek)
    Rule('MEK_phos_ERK_TpY', MEK(b=1) % ERK(b=1, T='u', Y='p') >>\
                            MEK(b=None) + ERK(b=None, T='p', Y='p'), kc_mek_p)
    Rule('MEK_phos_ERK_pTY', MEK(b=1) % ERK(b=1, T='p', Y='u') >>\
                            MEK(b=None) + ERK(b=None, T='p', Y='p'), kc_mek_p)

def mkp_dephos_erk_2_step():
    """Based on and simplified from Table S4 of Markevich et al. 
    Assuming distributive random phosphatase mechanism. It is simplified to 
    2-step mechanisms from the 3-step ones in the original table.
    
    First and second order rate constants are in s^{-1} and nM^{-1}s^{-1}
    respectively.
    """
    Parameter('kf_mkp', 2.9e-3)    # = h_{-6} + h_{-9}
    Parameter('kf_mkp_p', 2e-2)    # = h_{-3} + h_4 = h_7 + h_{-12}
    Parameter('kf_mkp_pp', 4.5e-2) # = h_1 = h_{10}
    Parameter('kr_mkp', 1.4e-1)    # = h_9 which is also used 
                                   #    in the reaction for the original h_6
    Parameter('kr_mkp_p', 1)       # = h_{3} = h_{-4} = h_{-7} = h_{12}
    Parameter('kr_mkp_pp', 1)      # = h_{-1} = h_{-10}
    Parameter('kc_mkp', 5e-1)      # = h_5 which is also used
                                   #    in the reaction for the original h_8
    Parameter('kc_mkp_p', 9.2e-2)  # = h_2 = h_{11}
    alias_model_components()
    
    Rule('MKP_bind_ERK_uu', MKP(b=None) + ERK(b=None, T='u', Y='u') <>\
                            MKP(b=1) % ERK(b=1, T='u', Y='u'), kf_mkp, 
                            kr_mkp)
    Rule('MKP_bind_ERK_pT', MKP(b=None) + ERK(b=None, T='p', Y='u') <>\
                            MKP(b=1) % ERK(b=1, T='p', Y='u'), kf_mkp_p,
                            kr_mkp_p)
    Rule('MKP_bind_ERK_pY', MKP(b=None) + ERK(b=None, T='u', Y='p') <>\
                            MKP(b=1) % ERK(b=1, T='u', Y='p'), kf_mkp_p,
                            kr_mkp_p)
    Rule('MKP_bind_ERK_pp', MKP(b=None) + ERK(b=None, T='p', Y='p') <>\
                            MKP(b=1) % ERK(b=1, T='p', Y='p'), 
                            kf_mkp_pp, kr_mkp_pp)
    
    Rule('MKP_dephos_ERK_ppT', MKP(b=1) % ERK(b=1, T='p', Y='p') >>\
                               MKP(b=None) + ERK(b=None, T='u', Y='p'),
                               kc_mkp_p)
    Rule('MKP_dephos_ERK_ppY', MKP(b=1) % ERK(b=1, T='p', Y='p') >>\
                               MKP(b=None) + ERK(b=None, T='p', Y='u'),
                               kc_mkp_p)
    Rule('MKP_dephos_ERK_pT', MKP(b=1) % ERK(b=1, T='p', Y='u') >>\
                              MKP(b=None) + ERK(b=None, T='u', Y='u'),
                              kc_mkp)
    Rule('MKP_dephos_ERK_pY', MKP(b=1) % ERK(b=1, T='u', Y='p') >>\
                              MKP(b=None) + ERK(b=None, T='u', Y='u'),
                              kc_mkp)

def erk_dimerize_any():
    Parameter('kf_erk', 10)
    Parameter('kr_erk', 1e-5)
    alias_model_components()
    Rule('ERK_bind_ERK', ERK(b=None) + ERK(b=None) <>\
                        ERK(b=1) % ERK(b=1), kf_erk, kr_erk)

def erk_dimerize_uT():
    Parameter('kf_erk', 10)
    Parameter('kr_erk', 1e-5)
    alias_model_components()
    Rule('ERK_bind_ERK', ERK(b=None, T='u') + ERK(b=None, T='u') <>\
                        ERK(b=1, T='u') % ERK(b=1, T='u'), kf_erk, kr_erk)

def erk_dimerize_phos():
    Parameter('kf_erk', 10)
    Parameter('kr_erk', 1e-5)
    alias_model_components()
    # These 3 rules cover all the combinations for any phosphorylated form
    # of ERK to dimerize. Only unphosphorylated for should not be able to
    # dimerize.
    Rule('ERK_bind_ERK', ERK(b=None, Y='p') + ERK(b=None, Y='p') 
                         <> ERK(b=1, Y='p') % ERK(b=1, Y='p'), 
                         kf_erk, kr_erk)
    Rule('ERK_bind_ERK', ERK(b=None, T='p', Y='u') + ERK(b=None, T='p', Y='u') 
                         <> ERK(b=1, T='p', Y='u') % ERK(b=1, T='p', Y='u'), 
                         kf_erk, kr_erk)
    Rule('ERK_bind_ERK', ERK(b=None, Y='p') + ERK(b=None, T='p', Y='u') 
                         <> ERK(b=1, Y='p') % ERK(b=1, T='p', Y='u'), 
                         kf_erk, kr_erk)

def erk_autophos():
    # * How does ERK autophosphorylate?
    # Does an already phosphorylated ERK phosphorylate another one?
    Parameter('kp_erk', 1e-8)
    alias_model_components()
    Rule('ERK_autophos', ERK(b=1, Y='u') % ERK(b=1, Y='p') >>\
                        ERK(b=1, Y='p') % ERK(b=1, Y='p'), kp_erk)

def erk_activate_mkp():
    Monomer('MKP_act', ['b'])
    Parameter('kc_erk_act_mkp', 1)
    Parameter('kc_mkp_act', 1)
    Parameter('kc_mkp_deact', 1)
    alias_model_components()
    
    Rule('ERK_release_act_MKP', MKP(b=1) % ERK(b=1, T='p', Y='p') >>\
                                MKP_act(b=None) + ERK(b=None, T='p', Y='p'),
                                kc_erk_act_mkp)

    Rule('MKP_act_bind_ERK_uu', MKP_act(b=None) + ERK(b=None, T='u', Y='u') <>\
                            MKP_act(b=1) % ERK(b=1, T='u', Y='u'), kf_mkp, 
                            kr_mkp)
    Rule('MKP_act_bind_ERK_pT', MKP_act(b=None) + ERK(b=None, T='p', Y='u') <>\
                            MKP_act(b=1) % ERK(b=1, T='p', Y='u'), kf_mkp_p, 
                            kr_mkp_p)
    Rule('MKP_act_bind_ERK_pY', MKP_act(b=None) + ERK(b=None, T='u', Y='p') <>\
                            MKP_act(b=1) % ERK(b=1, T='u', Y='p'), kf_mkp_p, 
                            kr_mkp_p)
    Rule('MKP_act_bind_ERK_pp', MKP_act(b=None) + ERK(b=None, T='p', Y='p') <>\
                            MKP_act(b=1) % ERK(b=1, T='p', Y='p'), kf_mkp_pp, 
                            kr_mkp_pp)
    
    Rule('MKP_act_dephos_ERK_ppT', MKP_act(b=1) % ERK(b=1, T='p', Y='p') >>\
                                MKP_act(b=None) + ERK(b=None, T='u', Y='p'), 
                                kc_mkp_act)
    Rule('MKP_act_dephos_ERK_ppY', MKP_act(b=1) % ERK(b=1, T='p', Y='p') >>\
                                MKP_act(b=None) + ERK(b=None, T='p', Y='u'), 
                                kc_mkp_act)
    Rule('MKP_act_dephos_ERK_pT', MKP_act(b=1) % ERK(b=1, T='p', Y='u') >>\
                                MKP(b=None) + ERK(b=None, T='u', Y='u'), 
                                kc_mkp_act)
    Rule('MKP_act_dephos_ERK_pY', MKP_act(b=1) % ERK(b=1, T='u', Y='p') >>\
                                MKP(b=None) + ERK(b=None, T='u', Y='u'), 
                                kc_mkp_act)
    
    Rule('MKP_deact', MKP_act(b=None) >> MKP(b=None), kc_mkp_deact)
    Observable('MKP_act_', MKP_act())

def mek_dephos_erk():
    Parameter('kc_mekdp', 1e-5)
    alias_model_components()
    Rule('MEK_dephos_ERKpTpY', MEK(b=1) % ERK(b=1, T='p', Y='p') >>\
                               MEK(b=None) + ERK(b=None, T='u', Y='p'), 
                               kc_mekdp)

def erk_autodephos():
    Parameter('kc_erkdp', 1e-5)
    alias_model_components()
    Rule('ERK_dephos_ERKpTpY', ERK(T='p', Y='p') % ERK(T='p', Y='p') >>\
                               ERK(T='p', Y='p') % ERK(T='u', Y='p'),
                               kc_erkdp)

def mapk_initials():
    Parameter('ERK_0', 1e3)
    Parameter('ERKpT_0', 0)
    Parameter('ERKpY_0', 0)
    Parameter('ERKpTpY_0', 0)
    Parameter('MEK_0', 1e2)
    Parameter('MKP_0', 1e2)
    alias_model_components()

    Initial(ERK(b=None, T='u', Y='u'), ERK_0)
    Initial(ERK(b=None, T='u', Y='p'), ERKpY_0)
    Initial(ERK(b=None, T='p', Y='u'), ERKpT_0)
    Initial(ERK(b=None, T='p', Y='p'), ERKpTpY_0)
    Initial(MEK(b=None), MEK_0)
    Initial(MKP(b=None), MKP_0)
    

def mapk_observables():
    Observable('ERKu', ERK(T='u', Y='u'))
    Observable('ERKpY', ERK(T='u', Y='p'))
    Observable('ERKpT', ERK(T='p', Y='u'))
    Observable('ERKpTpY', ERK(T='p', Y='p'))
    Observable('ERKd', ERK(b=1) % ERK(b=1))
    Observable('ERKdu', ERK(b=1, T='u', Y='u') % ERK(b=1, T='u', Y='u'))
    Observable('ERKdpY', ERK(b=1, T='u', Y='p') % ERK(b=1, T='u', Y='p'))
    Observable('ERKdpT', ERK(b=1, T='p', Y='u') % ERK(b=1, T='p', Y='u'))
    Observable('ERKdpTpY', ERK(b=1, T='p', Y='p') % ERK(b=1, T='p', Y='p'))
    alias_model_components()

if __name__ == '__main__':
    from pysb.integrate import Solver
    import numpy
    import matplotlib.pyplot as plt
   
    mapk_monomers()
    mek_phos_erk_2_step()
    mkp_dephos_erk_2_step()
    erk_autophos()
    mapk_initials()
    mapk_observables()
   
    ts = numpy.linspace(0, 7200, 100)
    solver = Solver(model, ts)

    solver.run()

    plt.plot(ts, solver.yobs['ERKpTpY'])
