"""
demo all types of response functions
make one instance of each type of response function
also demo loading/saving JSON and H5 files
"""

import numpy as np
import datetime as dt
import rfl

inst_info = {}

#%% top-level info
inst_info['FORMAT_VERSION'] = '1.1.0'
inst_info['L_UNIT'] = 'cm'
inst_info['E_UNIT'] = 'MeV'
inst_info['REFERENCES'] = ['Made-up sensors of all types',
    'Response file prepared by Paul O''Brien, paul.obrien@aero.org',
    'Created %s' % dt.datetime.utcnow().isoformat()]
inst_info['DEAD_TIME_PER_COUNT'] = 0 # unknown
inst_info['DEAD_TYPE'] = 'BLOCKING' # unknown
inst_info['COUNTS_MAX'] = np.inf # unknown
inst_info['SPECIES'] = ['PROT']
inst_info['CROSSCALIB'] = 1 # unknown
inst_info['CROSSCALIB_RMSE'] = np.log(2)/2 # unknown: assume 2 standard deviations is a factor of 2

tmp_chans = {}


# %% [E,TH,PH] - arbitrary geometry
# % "ARB" ETP_TYPE TBL
tmp_chans['ARB'] = {}
tmp_chans['ARB']['RESP_TYPE'] = '[E,TH,PH]' # E response depends on theta, phi
tmp_chans['ARB']['ETP_TYPE'] = 'TBL'
tmp_chans['ARB']['E_GRID'] = np.linspace(1,5,5)
tmp_chans['ARB']['TH_GRID'] = np.linspace(0,15,16) # response out to 15 degrees
tmp_chans['ARB']['PH_GRID'] = np.linspace(0,360,361) # response at all angles
EI,THI,PHI = np.meshgrid(tmp_chans['ARB']['E_GRID'],
                         tmp_chans['ARB']['TH_GRID'],
                         tmp_chans['ARB']['PH_GRID'],indexing='ij')
tmp_chans['ARB']['R'] = 3*np.maximum(0,np.arctan((EI-2)*10))*np.exp(-3*(THI/15)**2)*(1+rfl.cosd(PHI)**2)


# %% [E,TH] - arbitrary geometry, cylindrically symmetric
# % "ARBC" ET_TYPE TBL
tmp_chans['ARBC'] = {}
tmp_chans['ARBC']['RESP_TYPE'] = '[E,TH]' # E response depends on theta, cylindrically symmetric
tmp_chans['ARBC']['ET_TYPE'] = 'TBL'
tmp_chans['ARBC']['E_GRID'] = np.linspace(1,5,5)
tmp_chans['ARBC']['TH_GRID'] = np.linspace(0,15,16) # response out to 15 degrees
[EI,THI] = np.meshgrid(tmp_chans['ARBC']['E_GRID'],tmp_chans['ARBC']['TH_GRID'],indexing='ij')
tmp_chans['ARBC']['R'] = 3*np.maximum(0,np.arctan((EI-2)*10))*np.exp(-3*(THI/15)**2)

# %% [E],[TH,PH] - separable energy/angle response
# % "ESEP" E_TYPE INT, TP_TYPE TBL
tmp_chans['ESEP'] = {}
tmp_chans['ESEP']['RESP_TYPE']= '[E],[TH,PH]' # E independent of theta
tmp_chans['ESEP']['E_TYPE'] = 'INT';
tmp_chans['ESEP']['E0'] = 20;
tmp_chans['ESEP']['EPS'] = 1;
tmp_chans['ESEP']['TP_TYPE'] = 'TBL' #  table of theta,phi response
tmp_chans['ESEP']['TH_GRID'] = np.linspace(0,15,16) # response out to 15 degrees
tmp_chans['ESEP']['PH_GRID'] = np.linspace(0,360,361) # response at all angles
[THI,PHI] = np.meshgrid(tmp_chans['ESEP']['TH_GRID'],tmp_chans['ESEP']['PH_GRID'],indexing='ij')
tmp_chans['ESEP']['A'] = np.exp(-3*(THI/15)**2)*(1+rfl.cosd(PHI)**2)

# % "RECT_TELE" E_TYPE WIDE, TP_TYPE RECT_TELE, BIDIRECTIONAL FALSE
tmp_chans['RECT_TELE'] = {}
tmp_chans['RECT_TELE']['RESP_TYPE'] = '[E],[TH,PH]' # E independent of theta
tmp_chans['RECT_TELE']['E_TYPE'] = 'WIDE'
tmp_chans['RECT_TELE']['E0'] = 20
tmp_chans['RECT_TELE']['E1'] = 22
tmp_chans['RECT_TELE']['EPS'] = 1
tmp_chans['RECT_TELE']['TP_TYPE'] = 'RECT_TELE' # two-element rectangular telescope
tmp_chans['RECT_TELE']['W1'] = 1
tmp_chans['RECT_TELE']['H1'] = 2
tmp_chans['RECT_TELE']['W2'] = 3
tmp_chans['RECT_TELE']['H2'] = 4
tmp_chans['RECT_TELE']['D'] = 5
tmp_chans['RECT_TELE']['BIDIRECTIONAL'] = 'FALSE'

# %% [E],[TH] - separable energy/angle response, cylindrically symmetric
# % "PINHOLE" E_TYPE DIFF, TH_TYPE PINHOLE, BIDIRECTIONAL FALSE
tmp_chans['PINHOLE'] = {}
tmp_chans['PINHOLE']['RESP_TYPE'] = '[E],[TH]' # E independent of theta, cylindrically symmetric
tmp_chans['PINHOLE']['E_TYPE'] = 'WIDE'
tmp_chans['PINHOLE']['E0'] = 20
tmp_chans['PINHOLE']['E1'] = 22
tmp_chans['PINHOLE']['EPS'] = 1
tmp_chans['PINHOLE']['TH_TYPE'] = 'PINHOLE' # pinhole angular response
tmp_chans['PINHOLE']['G'] = 1
tmp_chans['PINHOLE']['BIDIRECTIONAL ']= 'FALSE'

# % "CYL_TELE" E_TYPE WIDE, TH_TYPE CYL_TELE, BIDIRECTIONAL FALSE
tmp_chans['CYL_TELE'] = {}
tmp_chans['CYL_TELE']['RESP_TYPE'] = '[E],[TH]' # E independent of theta, cylindrically symmetric
tmp_chans['CYL_TELE']['E_TYPE'] = 'WIDE'
tmp_chans['CYL_TELE']['E0'] = 20
tmp_chans['CYL_TELE']['E1'] = 22
tmp_chans['CYL_TELE']['EPS'] = 1
tmp_chans['CYL_TELE']['TH_TYPE'] = 'CYL_TELE' # Cylindrical telescope
tmp_chans['CYL_TELE']['R1'] = 1
tmp_chans['CYL_TELE']['R2'] = 2
tmp_chans['CYL_TELE']['D'] = 6
tmp_chans['CYL_TELE']['BIDIRECTIONAL'] = 'FALSE'

# % "DISK" E_TYPE WIDE, TH_TYPE DISK, BIDIRECTIONAL FALSE
tmp_chans['DISK'] = {}
tmp_chans['DISK']['RESP_TYPE'] = '[E],[TH]' # E independent of theta, cylindrically symmetric
tmp_chans['DISK']['E_TYPE'] = 'WIDE'
tmp_chans['DISK']['E0'] = 20
tmp_chans['DISK']['E1'] = 22
tmp_chans['DISK']['EPS'] = 1
tmp_chans['DISK']['TH_TYPE'] = 'DISK' # single-element disk
tmp_chans['DISK']['R1'] = 1
tmp_chans['DISK']['BIDIRECTIONAL'] = 'FALSE'

#% "SLAB" E_TYPE WIDE, TH_TYPE SLAB, BIDIRECTIONAL FALSE
tmp_chans['SLAB'] = {}
tmp_chans['SLAB']['RESP_TYPE'] = '[E],[TH]' # E independent of theta, cylindrically symmetric
tmp_chans['SLAB']['E_TYPE'] = 'WIDE'
tmp_chans['SLAB']['E0'] = 20
tmp_chans['SLAB']['E1'] = 22
tmp_chans['SLAB']['EPS'] = 1
tmp_chans['SLAB']['TH_TYPE'] = 'SLAB' # single-element slab
tmp_chans['SLAB']['W1'] = 1
tmp_chans['SLAB']['H1'] = 2
tmp_chans['SLAB']['BIDIRECTIONAL'] = 'FALSE'

#%% [E] - omnidirectional response
#% "OMNI" E_TYPE TBL, BIDIRECTIONAL FALSE
tmp_chans['OMNI'] = {}
tmp_chans['OMNI']['RESP_TYPE'] = '[E]' # E independent of angle
tmp_chans['OMNI']['E_TYPE'] = 'TBL'
tmp_chans['OMNI']['E_GRID'] = np.linspace(1,5,5)
tmp_chans['OMNI']['EPS'] = np.array([0,0.5,1,1,1])
tmp_chans['OMNI']['G'] = 3
tmp_chans['OMNI']['BIDIRECTIONAL'] = 'FALSE'

#%% make bidirectional versions of all where appropriate

new_chans = {}
for key,val in tmp_chans.items():
    if 'BIDIRECTIONAL' in val:
        key2 = key+'B'
        val2 = {**val,'BIDIRECTIONAL':'TRUE'}
        new_chans[key2] = val2
tmp_chans = {**tmp_chans,**new_chans} # mege

# % set channel names
inst_info['CHANNEL_NAMES'] = sorted(tmp_chans.keys())
for key in inst_info['CHANNEL_NAMES']:
    inst_info[key] = {'PROT':tmp_chans[key]}

jsonname = 'all_types.json'
h5name = 'all_types.h5'

rfl.write_JSON(inst_info,jsonname)
rfl.write_h5(inst_info,h5name)

inst_info = rfl.load_inst_info(inst_info)

# write them again, in case load_inst_info propagated some stuff downward
jsonnamebig = 'all_types_big.json'
h5namebig = 'all_types_big.h5'

rfl.write_JSON(inst_info,jsonnamebig)
rfl.write_h5(inst_info,h5namebig)

sp = 'PROT'

for info in [inst_info,jsonname,jsonnamebig,h5name,h5namebig]:
    if isinstance(info,str):
        print(info)
    else:
        print('object')

    info = rfl.load_inst_info(info)
    
    for v in info['CHANNEL_NAMES']:
        resp = info[v][sp]
        
        if hasattr(resp,'ar'):
            G = resp.ar.hA0
        elif hasattr(resp,'G'):
            G = resp.G
        elif resp.RESP_TYPE == '[E,TH]':
            grids = rfl.broadcast_grids(resp.E_GRID,resp.TH_GRID,0)
            tmp = resp.R(*grids)
            tmp = tmp[:,:,0]*2*np.pi # implicit integral of phi
            tmp = np.trapz(tmp,x=-rfl.cosd(resp.TH_GRID),axis=1)
            G = max(tmp)
        elif resp.RESP_TYPE == '[E,TH,PH]':
            grids = rfl.broadcast_grids(resp.E_GRID,resp.TH_GRID,resp.PH_GRID)
            tmp = resp.R(*grids)
            tmp = np.trapz(tmp,x=np.radians(resp.PH_GRID),axis=2)
            tmp = np.trapz(tmp,x=-rfl.cosd(resp.TH_GRID),axis=1)
            G = max(tmp)
        else:
            raise Exception('Could not figure out how to compute G: %s,%s' % (v,resp.RESP_TYPE))
        
        print('%s.%s.Geff = %g %s^2 sr' % (v,sp,G,resp.L_UNIT))
