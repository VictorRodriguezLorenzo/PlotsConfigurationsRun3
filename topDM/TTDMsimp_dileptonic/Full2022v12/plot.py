import ROOT

def C(hex_code):
    return ROOT.TColor.GetColor(hex_code)

###########################################
######## CMS / CVD-FRIENDLY COLORS ########
###########################################

cms_cb = {
    # Official 10-color CVD-friendly sequence
    'blue'        : C('#3F90DA'),  # 63, 144, 218
    'orange'      : C('#FFA90E'),  # 255, 169, 14
    'red'         : C('#BD1F01'),  # 189, 31, 1
    'gray'        : C('#94A4A2'),  # 148, 164, 162
    'purple'      : C('#832DB6'),  # 131, 45, 182
    'brown'       : C('#A96B59'),  # 169, 107, 89
    'dark_orange' : C('#E76300'),  # 231, 99, 0
    'tan'         : C('#B9AC70'),  # 185, 172, 112
    'dark_gray'   : C('#717581'),  # 113, 117, 129
    'light_blue'  : C('#92DADD'),  # 146, 218, 221
    'extra_pink'  : C('#C849A9'),  # 200, 73, 169

    'sig_red_1'   : C('#FF0000'),  # pure red
    'sig_red_2'   : C('#990000'),  # dark red
    'sig_red_3'   : C('#FF4500'),  # orange red
    'sig_red_4'   : C('#B00050'),  # wine / magenta-red
    'sig_red_5'   : C('#FF1493'),  # deep pink
    'sig_red_6'   : C('#8B0000'),  # dark red
    'sig_red_7'   : C('#DC143C'),  # crimson

    'black'       : C('#000000'),
}


###########################################
############### GROUP PLOT ################
###########################################

groupPlot = {}

###########################################
#############  BACKGROUNDS  ###############
###########################################

# smaller backgrounds first -> lower in stack and first in legend

groupPlot['VVV']  = {
    'nameHR'   : 'VVV',
    'isSignal' : 0,
    'color'    : cms_cb['light_blue'],
    'samples'  : ['WWW', 'WWZ', 'WZZ', 'ZZZ']
}

groupPlot['VV']  = {
    'nameHR'   : 'VV',
    'isSignal' : 0,
    'color'    : cms_cb['blue'],
    'samples'  : ['WW', 'WZ', 'ZZ']
}

groupPlot['TXX']  = {
    'nameHR'   : 't + X',
    'isSignal' : 0,
    'color'    : cms_cb['purple'],
    'samples'  : ['tHW', 'tHQ']
}

groupPlot['TTToSemiLeptonic']  = {
    'nameHR'   : 't#bar{t}(1l)',
    'isSignal' : 0,
    'color'    : cms_cb['gray'],
    'samples'  : ['TTToSemiLeptonic']
}

groupPlot['ttH'] = {
    'nameHR'   : 'ttH',
    'isSignal' : 0,
    'color'    : cms_cb['dark_orange'],
    'samples'  : ['ttH']
}

groupPlot['ttW'] = {
    'nameHR'   : 'ttW',
    'isSignal' : 0,
    'color'    : cms_cb['brown'],
    'samples'  : ['ttW']
}

groupPlot['ttZ'] = {
    'nameHR'   : 'ttZ',
    'isSignal' : 0,
    'color'    : cms_cb['extra_pink'],
    'samples'  : ['ttZ']
    #'samples'  : ['TTNuNu', 'TTLL_MLL-4to50', 'TTLL_MLL-50', 'TTZ-ZtoQQ']
    #'samples'  : ['TTLL_MLL-4to50', 'TTLL_MLL-50', 'TTZ-ZtoQQ']
    #'samples'  : ['TTNuNu']
}

groupPlot['DY']  = {
    'nameHR'   : 'DY',
    'isSignal' : 0,
    'color'    : cms_cb['tan'],
    'samples'  : ['DY']
}

groupPlot['ST']  = {
    'nameHR'   : 'Single Top',
    'isSignal' : 0,
    'color'    : cms_cb['orange'],
    'samples'  : ['ST']
}

groupPlot['Fake']  = {
    'nameHR'   : 'nonprompt',
    'isSignal' : 0,
    'color'    : cms_cb['dark_gray'],
    'colorPlt' : '#717581',
    'samples'  : ['Fake']
}

groupPlot['TTTo2L2Nu']  = {
    'nameHR'   : 't#bar{t}(2l)',
    'isSignal' : 0,
    'color'    : cms_cb['red'],
    'samples'  : ['TTTo2L2Nu']
}


###########################################
###############  SIGNALS  #################
###########################################

#mPhi = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600', '700', '800', '1000']
mPhi = ['10', '400', '1000']
#mPhi = ['1000']

signal_colors = [
    cms_cb['sig_red_1'],
    cms_cb['sig_red_2'],
    cms_cb['sig_red_3'],
    cms_cb['sig_red_4'],
    cms_cb['sig_red_5'],
    cms_cb['sig_red_6'],
    cms_cb['sig_red_7'],
]

color_by_mphi = {
    phi: color for phi, color in zip(mPhi, signal_colors)
}

'''
# tt+DM dilepton scalar
for phi in mPhi:
    groupPlot[f'TTto2LDMsimpSpin0_s_mphi-{phi}'] = {
        'nameHR': rf'#scale[0.90]{{TT2LDM s, m_{{#phi}}={phi}(x 10^{{2}})}}',
        'isSignal': 2,
        'color': color_by_mphi[phi],
        'samples': [f'TTto2LDMsimpSpin0_s_mphi-{phi}']
    }

'''

# tt+DM dilepton pseudoscalar
for phi in mPhi:
    groupPlot[f'TTto2LDMsimpSpin0_ps_mphi-{phi}'] = {
        'nameHR': rf'#scale[0.90]{{TT2LDM ps, m_{{#phi}}={phi}(x 10^{{2}})}}',
        'isSignal': 2,
        'color': color_by_mphi[phi],
        'samples': [f'TTto2LDMsimpSpin0_ps_mphi-{phi}']
    }

'''
# tW+DM dilepton scalar
for phi in mPhi:
    groupPlot[f'TWto2LDMsimpSpin0_s_mphi-{phi}'] = {
        'nameHR': rf'#scale[0.90]{{TW2LDM s, m_{{#phi}}={phi}(x 10^{{2}})}}',
        'isSignal': 2,
        'color': color_by_mphi[phi],
        'samples': [f'TWto2LDMsimpSpin0_s_mphi-{phi}']
    }

'''

# tW+DM dilepton pseudoscalar
for phi in mPhi:
    groupPlot[f'TWto2LDMsimpSpin0_ps_mphi-{phi}'] = {
        'nameHR': rf'#scale[0.90]{{TW2LDM ps, m_{{#phi}}={phi}(x 10^{{2}})}}',
        'isSignal': 2,
        'color': color_by_mphi[phi],
        'samples': [f'TWto2LDMsimpSpin0_ps_mphi-{phi}']
    }

###########################################
################## PLOT ###################
###########################################

plot = {}

###########################################
#############  BACKGROUNDS  ###############
###########################################

# smaller backgrounds first -> lower in stack

plot['WWW']  = {
    'color'    : cms_cb['light_blue'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['WWZ']  = {
    'color'    : cms_cb['light_blue'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['WZZ']  = {
    'color'    : cms_cb['light_blue'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['ZZZ']  = {
    'color'    : cms_cb['light_blue'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['WW']  = {
    'color'    : cms_cb['blue'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['WZ']  = {
    'color'    : cms_cb['blue'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['ZZ']  = {
    'color'    : cms_cb['blue'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['tHW']  = {
    'color'    : cms_cb['purple'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['tHQ']  = {
    'color'    : cms_cb['purple'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['TTToSemiLeptonic']  = {
    'color'    : cms_cb['gray'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['ttH']  = {
    'color'    : cms_cb['dark_orange'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['ttW']  = {
    'color'    : cms_cb['brown'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['ttZ'] = {
    'color'    : cms_cb['extra_pink'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

#plot['TTNuNu'] = {
#    'color'    : cms_cb['extra_pink'],
#    'isSignal' : 0,
#    'isData'   : 0,
#    'scale'    : 1.0,
#}
#
#plot['TTLL_MLL-4to50'] = {
#    'color'    : cms_cb['extra_pink'],
#    'isSignal' : 0,
#    'isData'   : 0,
#    'scale'    : 1.0,
#}
#
#plot['TTLL_MLL-50'] = {
#    'color'    : cms_cb['extra_pink'],
#    'isSignal' : 0,
#    'isData'   : 0,
#    'scale'    : 1.0,
#}
#
#plot['TTZ-ZtoQQ'] = {
#    'color'    : cms_cb['extra_pink'],
#    'isSignal' : 0,
#    'isData'   : 0,
#    'scale'    : 1.0,
#}

plot['DY']  = {
    'color'    : cms_cb['tan'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['ST']  = {
    'color'    : cms_cb['orange'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['Fake']  = {
    'color'    : cms_cb['dark_gray'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

plot['TTTo2L2Nu']  = {
    'color'    : cms_cb['red'],
    'isSignal' : 0,
    'isData'   : 0,
    'scale'    : 1.0,
}

###########################################
###############  SIGNALS  #################
###########################################
'''
# tt+DM dilepton scalar
for phi in mPhi:
    plot[f'TTto2LDMsimpSpin0_s_mphi-{phi}'] = {
        'color': color_by_mphi[phi],
        'isSignal': 2,
        'isData': 0,
        'scale': 1e2
    }
'''
# tt+DM dilepton pseudoscalar
for phi in mPhi:
    plot[f'TTto2LDMsimpSpin0_ps_mphi-{phi}'] = {
        'color': color_by_mphi[phi],
        'isSignal': 2,
        'isData': 0,
        'scale': 1e2
    }
'''
# tW+DM dilepton scalar
for phi in mPhi:
    plot[f'TWto2LDMsimpSpin0_s_mphi-{phi}'] = {
        'color': color_by_mphi[phi],
        'isSignal': 2,
        'isData': 0,
        'scale': 1e2
    }
'''
# tW+DM dilepton pseudoscalar
for phi in mPhi:
    plot[f'TWto2LDMsimpSpin0_ps_mphi-{phi}'] = {
        'color': color_by_mphi[phi],
        'isSignal': 2,
        'isData': 0,
        'scale': 1e2
    }

###########################################
#################  DATA  ##################
###########################################

plot['DATA']  = {
    'nameHR'   : 'Data',
    'color'    : cms_cb['black'],
    'isSignal' : 0,
    'isData'   : 1,
    'isBlind'  : 0
}

###########################################
################ LEGEND ###################
###########################################

legend = {}
legend['lumi'] = 'L = 7.98 fb^{-1}'
legend['sqrt'] = '#sqrt{s} = 13.6 TeV'
