cuts = {}

preselections = '((Lepton_pt[0] > 35 && abs(Lepton_pdgId[0]) == 11) || (Lepton_pt[0] > 30 && abs(Lepton_pdgId[0]) == 13)) \
            && abs(Lepton_eta[0]) < 2.4 \
            && mll > 20 \
            && noJetInHorn \
            && njets >= 2'

# CUTS

cuts["all"] = "Alt(Lepton_pt, 1, 0) < 10"

#######################
#### Signal region ####
#######################

cuts['ttdm_sr']  = {
   'expr': 'sr',
    # Sub-categorization of ttDM SR
   'categories' : {
      '1l_0f_t1' : 'topness[0] < 0 && nForwardjets == 0 && nbjets == 1',
      '1l_1f_t1' : 'topness[0] < 0 && nForwardjets >=1 && nbjets == 1',
      '1l_2b_t1' : 'topness[0] < 0 && nbjets >= 2',
      '1l_0f_t2' : 'topness[0] > 0 && nForwardjets == 0 && nbjets == 1',
      '1l_1f_t2' : 'topness[0] > 0 && nForwardjets >=1 && nbjets == 1',
      '1l_2b_t2' : 'topness[0] > 0 && nbjets >= 2',
   }
}

########################
#### Control regions ###
########################

cuts['ttcr']  = {
   'expr': 'ttcr',  
    # Sub-categorization of tt(2l) CR
   'categories' : {
       'inclusive' : '1',
       }
}

cuts['WjetsCR']  = {
   'expr': 'wjetscr',
    # Sub-categorization of W+jets CR
   'categories' : { 
       'Incl' : '1',
       }
}
