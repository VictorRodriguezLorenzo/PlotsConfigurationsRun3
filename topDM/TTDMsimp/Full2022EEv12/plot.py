groupPlot = {}

###########################################
#############  BACKGROUNDS  ###############
###########################################

groupPlot['DY']  = {  
    'nameHR'   : 'DY',
    'isSignal' : 0,
    'color'    : 418, #kGreen+4
    'samples'  : ['DY']
}

groupPlot['ST']  = {  
    'nameHR'   : 'Single Top',
    'isSignal' : 0,
    'color'    : 619, 
    'samples'  : ['ST']
}

groupPlot['TTTo2L2Nu']  = {  
    'nameHR'   : 't#bar{t}(2l)',
    'isSignal' : 0,
    'color'    : 400, 
    'samples'  : ['TTTo2L2Nu']
}

#groupPlot['TTToSemiLeptonic']  = {  
#    'nameHR'   : 't#bar{t}(1l)',
#    'isSignal' : 0,
#    'color'    : 810,   # kOrange + 10
#    'samples'  : ['TTToSemiLeptonic']
#}

groupPlot['VV']  = {  
    'nameHR'   : 'VV',
    'isSignal' : 0,
    'color'    : 857, # kAzure -3
    'samples'  : ['WW', 'WZ', 'ZZ']
}

groupPlot['Fake']  = {
    'nameHR' : 'nonprompt',
    'isSignal' : 0,
    'color': 921,    # kGray + 1
    'colorPlt': "#778899",
    'samples'  : ['Fake']
}

###########################################
###############  SIGNALS  #################
###########################################
'''
mPhi = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600' '700', '800', '1000']

for Zp in mZp:
    groupPlot['DH_' + hs  +  '_'   + DM + '_' + Zp]  = {
            'nameHR' :  'm_{s}=' + hs + ', m_{x}=' + DM + ', m_{Z}=' + Zp ,
            'isSignal' : 2,
            'color': 100, # kRed 
            'samples'  : ['DH_mhs_' + hs + '_mx_' + DM +  '_mZp_' + Zp]
            }

'''
plot = {}

###########################################
#############  BACKGROUNDS  ###############
###########################################

plot['DY']  = {  
    'color'    : 418,    # kGreen+2
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

plot['ST']  = {  
    'color'    : 619,
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

plot['TTTo2L2Nu']  = {  
    'color'    : 400, 
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

#plot['TTToSemiLeptonic']  = {  
#    'color'    : 810,   # kOrange + 10 
#    'isSignal' : 0,
#    'isData'   : 0, 
#    'scale'    : 1.0,
#}

plot['WW']  = {  
    'color'    : 857, # kAzure -3
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

plot['WZ']  = {  
    'color'    : 857, # kAzure -3
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

plot['ZZ']  = {  
    'color'    : 857, # kAzure -3
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

plot['Fake']  = {
    'color': 921,    # kGray + 1
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0
}

###########################################
###############  SIGNALS  #################
###########################################
'''
# tt+DM dilepton scalar
for Zp in mZp:
    plot['TTto2LDMsimpSpin0_s_mchi_{phi}']  = {
            'color': 100, # kRed
            'isSignal' : 2,
            'isData'   : 0,
            'scale'    : 1.0
            }

# tt+DM inclusive scalar
for Zp in mZp:
    plot[f'TTDMsimpSpin0_s_mchi_{phi}']  = {
            'color': 100, # kRed
            'isSignal' : 2,
            'isData'   : 0,
            'scale'    : 1.0
            }

# tt+DM dilepton pseudoscalar
for Zp in mZp:
    plot[f'TTto2LDMsimpSpin0_ps_mchi_{phi}']  = {
            'color': 100, # kRed
            'isSignal' : 2,
            'isData'   : 0,
            'scale'    : 1.0
            }
# tt+DM inclusive pseudoscalar
for Zp in mZp:
    plot[f'TTDMsimpSpin0_ps_mchi_{phi}']  = {
            'color': 100, # kRed
            'isSignal' : 2,
            'isData'   : 0,
            'scale'    : 1.0
            }
'''
# data

plot['DATA']  = { 
    'nameHR'   : 'Data',
    'color'    : 1 ,  
    'isSignal' : 0,
    'isData'   : 1 ,
    'isBlind'  : 0
}


# Legend definition
legend = {}
legend['lumi'] = 'L =  26.67 fb^{-1}'
legend['sqrt'] = '#sqrt{s} = 13.6 TeV'
