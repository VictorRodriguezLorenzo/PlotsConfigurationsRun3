cuts = {}

preselections = '((abs(Lepton_pdgId[0]) == 11 || abs(Lepton_pdgId[0]) == 13) && (abs(Lepton_pdgId[1]) == 11 || abs(Lepton_pdgId[1]) == 13)) \
            && Lepton_pt[0] > 25 \
            && Lepton_pt[1] > 20 \
            && abs(Lepton_eta[0]) < 2.4 && abs(Lepton_eta[1]) < 2.4 \
            && mll > 20 \
            && noJetInHorn_pT15 \
            && bReq'

# CUTS

#cuts["all"] = "Alt(Lepton_pt, 2, 0) < 10"

#######################
#### Signal region ####
#######################

cuts['ttdm_sr']  = {
   'expr': 'sr',
    # Sub-categorization of SR
   'categories' : {
      '2l_1b' : 'nLepton == 2 && nbjets == 1',
      '2l_2b' : 'nLepton == 2 && nbjets >= 2 && doubleNu_producer[9]',
   }
}

#######################
### Control regions ###
#######################

#cuts['ttvr']  = {
#   'expr': 'ttvr',  
#    # Sub-categorization of VR
#   'categories' : {
#       'inclusive' : '1',
#       }
#}
#
#cuts['dycr']  = {
#   'expr': 'dycr',
#    # Sub-categorization of DY CR
#   'categories' : { 
#       'inclusive' : '1',
#       }
#}
#
#cuts['ttZcr']  = {
#   'expr': 'ttZcr',
#    # Sub-categorization of ttZ CR
#   'categories' : {
#       'inclusive' : '1',
#       }
#}
