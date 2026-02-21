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

groupPlot['TTToSemiLeptonic']  = {  
    'nameHR'   : 't#bar{t}(1l)',
    'isSignal' : 0,
    'color'    : 810,   # kOrange + 10
    'samples'  : ['TTToSemiLeptonic']
}

groupPlot['TTV']  = {  
    'nameHR'   : 'TTV',
    'isSignal' : 0,
    'color'    : 616,   # kMagenta
    'samples'  : ['TTV']
}

groupPlot['TXX']  = {
    'nameHR'   : 't + X',
    'isSignal' : 0,
    'color'    : 870,   # kAzure + 10
    'samples'  : ['TXX']
}

groupPlot['VV']  = {  
    'nameHR'   : 'VV',
    'isSignal' : 0,
    'color'    : 857, # kAzure -3
    'samples'  : ['WW', 'WZ', 'ZZ']
}

groupPlot['VVV']  = {
    'nameHR'   : 'VVV',
    'isSignal' : 0,
    'color'    : 432, # kCyan
    'samples'  : ['WWW', 'WWZ', 'WZZ', 'ZZZ']
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
#mPhi = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600' '700', '800', '1000']
'''
mPhi = ['600']
'''
# tt+DM dilepton scalar
for phi in mPhi:
    groupPlot['TTto2LDMsimpSpin0_s_mphi-{phi}']  = {
            'nameHR' :  'TTto2LDM scalar m_{phi}=' + phi,
            'isSignal' : 2,
            'color': 100, # kRed 
            'samples'  : [f'TTto2LDMsimpSpin0_s_mphi-{phi}']
            }
'''
# tt+DM inclusive scalar
for phi in mPhi:
   groupPlot[f'TTDMsimpSpin0_s_mphi-{phi}']  = {
            'nameHR' :  'TTtoDM scalar m_{#phi}=' + phi + ' (x 10^{3})',
            'isSignal' : 2,
            'color': 632, # kRed 
            'samples'  : [f'TTDMsimpSpin0_s_mphi-{phi}']
            }
'''
# tt+DM dilepton pseudoscalar
for phi in mPhi:
    groupPlot['TTto2LDMsimpSpin0_ps_mphi-{phi}']  = {
            'nameHR' :  'TTto2LDM pseudoscalar m_{phi}=' + phi,
            'isSignal' : 2,
            'color': 100, # kRed 
            'samples'  : [f'TTto2LDMsimpSpin0_ps_mphi-{phi}']
            }

# tt+DM inclusive pseudoscalar
for phi in mPhi:
    groupPlot[f'TTDMsimpSpin0_ps_mphi-{phi}']  = {
            'nameHR' :  'TTtoDM pseudoscalar m_{phi}=' + phi,
            'isSignal' : 2,
            'color': 100, # kRed 
            'samples'  : [f'TTDMsimpSpin0_ps_mphi-{phi}']
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

plot['TTToSemiLeptonic']  = {  
    'color'    : 810,   # kOrange + 10 
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

plot['TTV']  = {  
    'color'    : 616,   # kMagenta
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

plot['TXX']  = {  
    'color'    : 870,   # kMagenta
    'isSignal' : 0,
    'isData'   : 0, 
    'scale'    : 1.0,
}

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

plot['WWW']  = {
    'color'    : 432, # kCyan
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}
plot['WWZ']  = {
    'color'    : 432, # kCyan
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['WZZ']  = {
    'color'    : 432, # kCyan
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['ZZZ']  = {
    'color'    : 432, # kCyan
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
for phi in mPhi:
    plot['TTto2LDMsimpSpin0_s_mphi-{phi}']  = {
            'color': 100, # kRed
            'isSignal' : 2,
            'isData'   : 0,
            'scale'    : 1.0
            }
'''
# tt+DM inclusive scalar
for phi in mPhi:
    plot[f'TTDMsimpSpin0_s_mphi-{phi}']  = {
            'color': 100, # kRed
            'isSignal' : 2,
            'isData'   : 0,
            'scale'    : 1.0
            }
'''
# tt+DM dilepton pseudoscalar
for phi in mPhi:
    plot[f'TTto2LDMsimpSpin0_ps_mphi-{phi}']  = {
            'color': 100, # kRed
            'isSignal' : 2,
            'isData'   : 0,
            'scale'    : 1.0
            }
# tt+DM inclusive pseudoscalar
for phi in mPhi:
    plot[f'TTDMsimpSpin0_ps_mphi-{phi}']  = {
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
legend['lumi'] = 'L =  9.451 fb^{-1}'
legend['sqrt'] = '#sqrt{s} = 13.6 TeV'
