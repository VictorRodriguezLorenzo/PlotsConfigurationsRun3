# variables

# 0 = not fold (default)
# 1 = fold underflow bin
# 2 = fold overflow bin
# 3 = fold underflow and overflow bins

# Toggle SR data blinding on/off
blindSR = True

blindCuts = {
    f'ttdm_sr_{category}': 'full'
    for category in cuts['ttdm_sr']['categories']
} if blindSR else {}

# variables
variables = {}

variables['events'] = {
    'name'  : '1',      
    'range' : (1,0,2),  
    'xaxis' : 'events', 
    'fold'  : 3,
    'blind' : blindCuts
}

#variables['nvtx'] = {     
#    'name'  : 'PV_npvsGood',      
#    'range' : (100, 0, 100),  
#    'xaxis' : 'number of vertices', 
#    'fold'  : 3,
#    'blind' : blindCuts
#}
#
#variables['mll'] = {
#    'name': 'mll',    
#    'range' : (60,60,120), 
#    'xaxis' : 'm_{ll} [GeV]',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['mth'] = {
#    'name': 'mth',
#    'range' : (60,0,300),
#    'xaxis' : 'm_{T}^{H} [GeV]',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['mtw1']  = {
#    'name': 'mtw1',
#    'range' : (50, 0,200),
#    'xaxis' : 'm_{T}^{W_{1}} [GeV]',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['mtw2']  = {
#    'name': 'mtw2',
#    'range' : (50, 0,100),
#    'xaxis' : 'm_{T}^{W_{2}} [GeV]',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['ptll']  = {  
#    'name': 'ptll',     
#    'range' : (20, 0,200),   
#    'xaxis' : 'p_{T}^{ll} [GeV]',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['drll']  = {
#    'name': 'drll',
#    'range' : (50, 0,5),
#    'xaxis' : '#Delta R_{ll}',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['dphill']  = {
#    'name': 'abs(dphill)',
#    'range' : (50, 0,3.15),
#    'xaxis' : '#Delta #phi_{ll}',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['detall'] = {
#    'name'  : 'abs(detall)',
#    'range' : (40, 0., 3.15),
#    'xaxis' : '|#Delta#eta_{ll}|',
#    'fold'  : 3,
#    'blind' : blindCuts
#}
#
#variables['pt1']  = { 
#    'name': 'Lepton_pt[0]',     
#    'range' : (20,25,150),
#    'xaxis' : 'p_{T} 1st lep',
#    'fold'  : 3,                         
#    'blind' : blindCuts
#}
#
#variables['pt2']  = {
#    'name': 'Lepton_pt[1]',     
#    'range' : (20,0,125),   
#    'xaxis' : 'p_{T} 2nd lep',
#    'fold'  : 3,                        
#    'blind' : blindCuts
#}
#
#variables['pt3']  = {
#    'name': 'Lepton_pt[2]',
#    'range' : (20,0,100),
#    'xaxis' : 'p_{T} 3rd lep',
#    'fold'  : 3,
#    'blind' : blindCuts
#}
#
#variables['eta1']  = {
#    'name': 'Lepton_eta[0]',     
#    'range' : (40,-3,3),   
#    'xaxis' : '#eta 1st lep',
#    'fold'  : 3,                         
#    'blind' : blindCuts
#}
#
#variables['eta2']  = {
#    'name': 'Lepton_eta[1]',     
#    'range' : (40,-3,3),   
#    'xaxis' : '#eta 2nd lep',
#    'fold'  : 3,                        
#    'blind' : blindCuts
#}
#
#                        
#variables['phi1']  = {
#    'name': 'Lepton_phi[0]',
#    'range' : (20,-3.2,3.2),
#    'xaxis' : '#phi 1st lep',
#    'fold'  : 3,
#    'blind' : blindCuts
#}
#
#variables['phi2']  = {
#    'name': 'Lepton_phi[1]',
#    'range' : (20,-3.2,3.2),
#    'xaxis' : '#phi 2nd lep',
#    'fold'  : 3,
#    'blind' : blindCuts
#}
#
##variables['jetdeepb']  = {
##    'name': 'Alt(Take(Jet_btagDeepFlavB, CleanJet_jetIdx), 0, -99)',
##    'range' : (30,0,1),
##    'xaxis' : 'B tagger 1st jet (DeepB)',
##    'fold' : 0,
##    'blind' : blindCuts
##}
##variables['jetpnetb']  = {
##    'name': 'Alt(Take(Jet_btagPNetB, CleanJet_jetIdx), 0, -99)',
##    'range' : (30,0,1),
##    'xaxis' : 'B tagger 1st jet (PNetB)',
##    'fold' : 0,
##    'blind' : blindCuts
##}
##variables['jetupartb']  = {
##    'name': 'Alt(Take(Jet_btagUParTAK4B, CleanJet_jetIdx), 0, -99)',
##    'range' : (30,0,1),
##    'xaxis' : 'B tagger 1st jet (UParT)',
##    'fold' : 0,
##    'blind' : blindCuts
##}
#
## MET
#variables['trkMet']  = { 
#    'name': 'TkMET_pt',
#    'range' : (20,0,200),
#    'xaxis' : 'trk met [GeV]',
#    'fold' : 3,
#    'blind' : blindCuts
#}
#
#variables['puppimet']  = {
#    'name': 'PuppiMET_pt',
#    'range' : (20,0,200),
#    'xaxis' : 'Puppi MET p_{T} [GeV]',
#    'fold' : 3,
#    'blind' : blindCuts
#}
#
#variables['mpmet']  = {
#    'name': 'mpmet',
#    'range' : (40,0,100),
#    'xaxis' : 'Min. proj. MET p_{T} [GeV]',
#    'fold' : 3,
#    'blind' : blindCuts
#}
#
#variables['njet']  = {
#    'name': 'njets',
#    'range' : (5,1,5),
#    'xaxis' : 'Number of jets',
#    'fold' : 2,
#    'blind' : blindCuts
#}
#
#variables['nbjet']  = {
#    'name': 'nbjets',
#    'range' : (5,1,5),
#    'xaxis' : 'Number of b-jets',
#    'fold' : 2,
#    'blind' : blindCuts
#}
#
#variables['vht_pt'] = {
#    'name': 'vht_pt',
#    'range': (40, 0, 800),
#    'xaxis': 'p_{T}(#ell#ell + jets) [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['jetpt1']  = {
#    #'name': 'Alt(CleanJet_pt, 0, -99) - 9999.9*(CleanJet_pt[0]<30)',
#    'name': 'Alt(CleanJet_pt, 0, -99)',
#    'range' : (40,0,200),
#    'xaxis' : 'p_{T} 1st jet',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['jetpt2']  = {
#    #'name': 'Alt(CleanJet_pt, 1, -99)  - 9999.9*(CleanJet_pt[1]<30)',
#    'name': 'Alt(CleanJet_pt, 1, -99)',
#    'range' : (40,0,200),
#    'xaxis' : 'p_{T} 2nd jet',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['jeteta1']  = {
#    #'name': 'Alt(CleanJet_eta, 0, -99) - 9999.9*(CleanJet_pt[0]<30)',
#    'name': 'Alt(CleanJet_eta, 0, -99)',
#    'range' : (30,-4.7,4.7),
#    'xaxis' : '#eta 1st jet',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['jeteta2']  = {
#    #'name': 'Alt(CleanJet_eta, 1, -99) - 9999.9*(CleanJet_pt[1]<30)',
#    'name': 'Alt(CleanJet_eta, 1, -99)',
#    'range' : (30,-4.7,4.7),
#    'xaxis' : '#eta 2nd jet',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#"""
####### Nº b-jets
#
#btagging_WPs = {
#    "DeepFlavB" : {
#        "loose"    : "0.0479",
#        "medium"   : "0.2435",
#    },
#    "RobustParTAK4B" : {
#        "loose"    : "0.0681",
#        "medium"   : "0.3494",
#    },
#    "PNetB" : {
#        "loose"    : "0.0358",
#        "medium"   : "0.1919",
#    }
#}
#
## Algo / SF name
#btagging_SFs = {
#    "DeepFlavB"      : "deepjet",
#    "RobustParTAK4B" : "RobustParT",
#    "PNetB"          : "partNet",
#}
#
#for bAlgo in btagging_SFs.keys():
#
#    bWP = btagging_WPs[bAlgo]["loose"]
#    
#    variables[f'nbjet_{bAlgo}']  = {
#        'name': f'Sum(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{bAlgo}, CleanJet_jetIdx) > {bWP})',
#        'range' : (5,0,5),
#        'xaxis' : f'Number of b-jets ({bAlgo})',
#        'fold' : 2
#    }
#"""
#
#
################################
#### ttDM specific variables ###
################################
#variables['mT2']  = {
#    'name': 'mT2',
#    'range' : (30,0,200),
#    'xaxis' : 'M_{ll}^{T2} [GeV]',
#    'fold' : 2,
#    'blind' : blindCuts
#}
#
#variables['mt2blbl'] = {
#    'name': 'mt2blbl',
#    'range': (40, 0, 400),
#    'xaxis': 'M_{T2}^{bl,bl} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['pdark']  = {
#    'name': 'doubleNu_producer[8]',
#    'range' : (30,0,600),
#    'xaxis' : 'p_{T}^{dark} [GeV]',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['chel']  = {
#    'name': 'doubleNu_producer[6]',
#    'range' : (30,-1,1),
#    'xaxis' : 'cos(#phi_{ll}^{*})',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['dphi_ttbar']  = {
#    'name': 'doubleNu_producer[7]',
#    'range' : (30,0,3.15),
#    'xaxis' : '|#Delta#phi(t,#bar{t})|',
#    'fold' : 0,
#    'blind' : blindCuts
#}
#
#variables['tt_reco']  = {
#    'name': 'doubleNu_producer[9]',
#    'range' : (2,0,1),
#    'xaxis' : 'tt-reco success?',
#    'fold' : 3,
#    'blind' : blindCuts
#}

### DISCRIMINANTS FOR ttDM DNN ###
for i, phi in enumerate(mPhi):
    variables[f'evaluate_dnn_ttDM_ps_{phi}']  =  { 
            'name': f'evaluate_dnn_ttDM_ps[{i}]',
            'range' : (30,0,1),
            'xaxis': rf'DNN Discriminant ttDM ps m_{{#Phi}} = {phi} GeV',
            'fold' : 3,
            'blind' : blindCuts
            }

for i, phi in enumerate(mPhi):
    variables[f'evaluate_dnn_ttDM_s_{phi}']  =  { 
            'name': f'evaluate_dnn_ttDM_s[{i}]',
            'range' : (30,0,1),
            'xaxis': rf'DNN Discriminant ttDM s m_{{#Phi}} = {phi} GeV',
            'fold' : 3,
            'blind' : blindCuts
            }



### DISCRIMINANTS FOR tWDM DNN ###
for i, phi in enumerate(mPhi):
    variables[f'evaluate_dnn_tWDM_ps_{phi}']  =  { 
            'name': f'evaluate_dnn_tWDM_ps[{i}]',
            'range' : (30,0,1),
            'xaxis': rf'DNN Discriminant tWDM ps m_{{#Phi}} = {phi} GeV',
            'fold' : 3,
            'blind' : blindCuts
            }

for i, phi in enumerate(mPhi):
    variables[f'evaluate_dnn_tWDM_s_{phi}']  =  { 
            'name': f'evaluate_dnn_tWDM_s[{i}]',
            'range' : (30,0,1),
            'xaxis': rf'DNN Discriminant tWDM s m_{{#Phi}} = {phi} GeV',
            'fold' : 3,
            'blind' : blindCuts
            }

#### New topDM variables ###
#
#variables['dphi_met_ll'] = {
#    'name': 'dphi_met_ll',
#    'range': (30, 0, 3.15),
#    'xaxis': '|#Delta#phi(MET,ll)|',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['ht'] = {
#    'name': 'ht',
#    'range': (40, 0, 800),
#    'xaxis': 'H_{T} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['st'] = {
#    'name': 'st',
#    'range': (40, 0, 1200),
#    'xaxis': 'S_{T} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['met_over_sqrt_ht'] = {
#    'name': 'met_over_sqrt_ht',
#    'range': (40, 0, 40),
#    'xaxis': 'p_{T}^{miss}/#sqrt{H_{T}}',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['met_over_st'] = {
#    'name': 'met_over_st',
#    'range': (40, 0, 1),
#    'xaxis': 'p_{T}^{miss}/S_{T}',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['dphi_min_j_met'] = {
#    'name': 'dphi_min_j_met',
#    'range': (30, 0, 3.15),
#    'xaxis': 'min |#Delta#phi(j,MET)|',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['pt_llb'] = {
#    'name': 'pt_llb',
#    'range': (40, 0, 700),
#    'xaxis': 'p_{T}^{llb} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['pt_llbb'] = {
#    'name': 'pt_llbb',
#    'range': (40, 0, 800),
#    'xaxis': 'p_{T}^{llbb} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['dphi_met_llbb'] = {
#    'name': 'dphi_met_llbb',
#    'range': (30, 0, 3.15),
#    'xaxis': '|#Delta#phi(MET,llbb)|',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['m_llbb'] = {
#    'name': 'm_llbb',
#    'range': (40, 0, 1000),
#    'xaxis': 'm_{llbb} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['mT_llbb_met'] = {
#    'name': 'mT_llbb_met',
#    'range': (40, 0, 1200),
#    'xaxis': 'm_{T}(llbb,MET) [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['met_over_pt_llbb'] = {
#    'name': 'met_over_pt_llbb',
#    'range': (40, 0, 5),
#    'xaxis': 'p_{T}^{miss}/p_{T}^{llbb}',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['mt2_bell_l'] = {
#    'name': 'mt2_bell_l',
#    'range': (40, 0, 400),
#    'xaxis': 'M_{T2}^{b#ell,#ell} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['mbl_min'] = {
#    'name': 'mbl_min',
#    'range': (40, 0, 300),
#    'xaxis': 'min m_{b#ell} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['mbl_max'] = {
#    'name': 'mbl_max',
#    'range': (40, 0, 500),
#    'xaxis': 'max m_{b#ell} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['mtbl'] = {
#    'name': 'mtbl',
#    'range': (40, 0, 400),
#    'xaxis': 'm_{b#ell}^{minimax} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['mbb'] = {
#    'name': 'mbb',
#    'range': (40, 0, 500),
#    'xaxis': 'm_{bb} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['drbb'] = {
#    'name': 'drbb',
#    'range': (40, 0, 5),
#    'xaxis': '#DeltaR_{bb}',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['ptbb'] = {
#    'name': 'ptbb',
#    'range': (40, 0, 600),
#    'xaxis': 'p_{T}^{bb} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['max_nonleading_btag'] = {
#    'name': 'max_nonleading_btag',
#    'range': (40, 0, 1),
#    'xaxis': 'max non-leading b tag score',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['pt_b2'] = {
#    'name': 'pt_b2',
#    'range': (40, 0, 300),
#    'xaxis': 'p_{T} 2nd b jet [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['nForwardJet'] = {
#    'name': 'nForwardJet',
#    'range': (5, 0, 5),
#    'xaxis': 'N_{forward jets}',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['leadingForwardJet_pt'] = {
#    'name': 'leadingForwardJet_pt',
#    'range': (40, 0, 300),
#    'xaxis': 'leading forward jet p_{T} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['leadingForwardJet_absEta'] = {
#    'name': 'leadingForwardJet_absEta',
#    'range': (25, 2.5, 5.0),
#    'xaxis': 'leading forward jet |#eta|',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['deta_forwardJet_b'] = {
#    'name': 'deta_forwardJet_b',
#    'range': (40, 0, 7),
#    'xaxis': '|#Delta#eta(j_{fwd},b)|',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['dphi_forwardJet_met'] = {
#    'name': 'dphi_forwardJet_met',
#    'range': (30, 0, 3.15),
#    'xaxis': '|#Delta#phi(j_{fwd},MET)|',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['top1_pt_reco'] = {
#    'name': 'top1_pt_reco',
#    'range': (40, 0, 600),
#    'xaxis': 'reco top 1 p_{T} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['top2_pt_reco'] = {
#    'name': 'top2_pt_reco',
#    'range': (40, 0, 600),
#    'xaxis': 'reco top 2 p_{T} [GeV]',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['pdark_over_met'] = {
#    'name': 'pdark_over_met',
#    'range': (40, 0, 5),
#    'xaxis': 'p_{T}^{dark}/p_{T}^{miss}',
#    'fold': 2,
#    'blind' : blindCuts
#}
#
#variables['angle_ll_llbb_rf'] = {
#    'name': 'angle_ll_llbb_rf',
#    'range': (30, 0, 3.15),
#    'xaxis': '#angle(#ell,#ell) in llbb rest frame',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['dphi_ll_llbb_rf'] = {
#    'name': 'dphi_ll_llbb_rf',
#    'range': (30, 0, 3.15),
#    'xaxis': '#Delta#phi(#ell,#ell) in llbb rest frame',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['cos_l1_llbb_rf'] = {
#    'name': 'cos_l1_llbb_rf',
#    'range': (30, -1, 1),
#    'xaxis': 'cos#theta(#ell_{1}, llbb axis)',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['cos_l2_llbb_rf'] = {
#    'name': 'cos_l2_llbb_rf',
#    'range': (30, -1, 1),
#    'xaxis': 'cos#theta(#ell_{2}, llbb axis)',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['angle_ll_llmet_rf'] = {
#    'name': 'angle_ll_llmet_rf',
#    'range': (30, 0, 3.15),
#    'xaxis': '#angle(#ell,#ell) in ll+MET rest frame',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['dphi_ll_llmet_rf'] = {
#    'name': 'dphi_ll_llmet_rf',
#    'range': (30, 0, 3.15),
#    'xaxis': '#Delta#phi(#ell,#ell) in ll+MET rest frame',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['cos_l1_llmet_rf'] = {
#    'name': 'cos_l1_llmet_rf',
#    'range': (30, -1, 1),
#    'xaxis': 'cos#theta(#ell_{1}, ll+MET axis)',
#    'fold': 0,
#    'blind' : blindCuts
#}
#
#variables['cos_l2_llmet_rf'] = {
#    'name': 'cos_l2_llmet_rf',
#    'range': (30, -1, 1),
#    'xaxis': 'cos#theta(#ell_{2}, ll+MET axis)',
#    'fold': 0,
#    'blind' : blindCuts
#}
