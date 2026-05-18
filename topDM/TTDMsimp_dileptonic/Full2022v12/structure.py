# structure configuration for datacard

structure = {}

# keys here must match keys in samples.py    

structure['DY']  = {  
                  'isSignal' : 0,
                  'isData'   : 0,
              }

structure['TTTo2L2Nu'] = {   
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['TTToSemiLeptonic'] = {   
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['ST'] = {
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['TTToSemiLeptonic'] = {
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['ttW']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['ttZ']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['ttH']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['TXX']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['WWW']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['WWZ']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['WZZ']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['ZZZ']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['WW']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['WZ']  = { 
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['ZZ']  = {
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

structure['Fake']  = {
                  'isSignal' : 0,
                  'isData'   : 0,
                  }

mPhi = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600', '700', '800', '1000']

# tt+DM dilepton scalar
for phi in mPhi:
    structure[f'TTto2LDMsimpSpin0_s_mphi-{phi}'] = {
        'isSignal': 2,
        'isData': 0
    }

# tt+DM inclusive scalar
for phi in mPhi:
    structure[f'TTDMsimpSpin0_s_mphi-{phi}'] = {
        'isSignal': 2,
        'isData': 0
    }

# tt+DM dilepton pseudoscalar
for phi in mPhi:
    structure[f'TTto2LDMsimpSpin0_ps_mphi-{phi}'] = {
        'isSignal': 2,
        'isData': 0
    }

# tt+DM inclusive pseudoscalar
for phi in mPhi:
    structure[f'TTDMsimpSpin0_ps_mphi-{phi}'] = {
        'isSignal': 2,
        'isData': 0
    }

# data

structure['DATA']  = { 
                  'isSignal' : 0,
                  'isData'   : 1,
              }
