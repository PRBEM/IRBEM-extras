"""
demonstrate claculations that lead to RPP result using CREME sample data

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import creme as cr
from degraded_spectra import SRIMTable,SphereLayer
from rpp import ProtonPart, IonPartRPP

# prepare figure path
_TOP_PATH = os.path.dirname(os.path.abspath(__file__))
figpath = os.path.join(_TOP_PATH,'figures')
if not os.path.exists(figpath):
    os.mkdir(figpath)

data_path = os.path.join(_TOP_PATH,'data')

AlTable = SRIMTable('SRIM_Al.csv')
shield = SphereLayer(AlTable,'Al',100,units='mils')
SiTable = SRIMTable('SRIM_Si.csv')

# read CREME files
# flx - differential flux outside spacecraft, #/m^s/sr/s/(MeV/nuc)
# tfx - differential flux inside shielding, #/m^s/sr/s/(MeV/nuc)
# dlt - differential flux w.r.t LET, #/m^s/sr/s (LET in MeV/(g/cm^2))
files = {}
files['flx'] = cr.read_flux(os.path.join(data_path,'q800a90i.FLX'))
files['tfx'] = cr.read_flux(os.path.join(data_path,'q800a90i100.TFX'))
files['dlt'] = cr.read_flux(os.path.join(data_path,'q800a90i100Z2.DLT')) # excludeZ=1

# compute species-specific energies & differential flux
for suffix in ['flx','tfx']:
    f = files[suffix]
    A = np.round(np.array([AlTable.data[AlTable.Z2sym[Z]]['MAI'] for Z in f['Z']])) # MAI mass number
    f['MeV'] = np.outer(f['MeV/nuc'],A) # MeV/nuc -> MeV
    f['dflux'] = f['flux']/A.reshape((1,len(A)))*1e-4*4*np.pi # -> #/cm^2/s/MeV

files['dlt']['dflux'] = files['dlt']['flux']*1e-4*4*np.pi # -> #/cm^s/s/(MeV/(g/cm^2)

for f in files.values(): # ensure positivity
    f['dflux'] = np.maximum(f['dflux'],0) 

syms = [AlTable.Z2sym[Z] for Z in files['flx']['Z']]
degraded = np.full(files['flx']['dflux'].shape,np.nan)

for (i,s) in enumerate(syms):
    degraded[:,i] = shield.degraded_spectrum(files['flx']['MeV'][:,i],files['flx']['dflux'][:,i],fast=True,species=s)

plt.close('all')

#%% plot CREME fluxes
plt.figure()
for s in ['H','He','C','Fe','U']:
    Z = AlTable.data[s]['Z']
    icol = np.where(files['flx']['Z']==Z)[0][0]
    h1 = plt.loglog(files['flx']['MeV/nuc'],files['flx']['flux'][:,icol],'-',label='%s incident' % s)
    color = h1[0].get_color()
    jcol = np.where(files['tfx']['Z']==Z)[0][0]
    plt.loglog(files['tfx']['MeV/nuc'],files['tfx']['flux'][:,jcol],':',label='%s CREME' % s,color=color)
    
plt.xlabel('MeV/nuc')
plt.ylabel('Incident Flux #/m$^2$/sr/s/(MeV/nuc)')
plt.subplots_adjust(right=0.75)
plt.legend(loc=(1.04,0))
plt.savefig(os.path.join(figpath,'creme_fluxes.png'))

#%% plot
plt.figure()
for s in ['H','He','C','Fe','U']:
    Z = AlTable.data[s]['Z']
    icol = np.where(files['flx']['Z']==Z)[0][0]
    h1 = plt.loglog(files['flx']['MeV/nuc'],files['flx']['dflux'][:,icol],'-',label='%s incident' % s)
    color = h1[0].get_color()
    plt.loglog(files['flx']['MeV/nuc'],degraded[:,icol],'--',label='%s CSDA' % s,color=color)
    jcol = np.where(files['tfx']['Z']==Z)[0][0]
    plt.loglog(files['tfx']['MeV/nuc'],files['tfx']['dflux'][:,jcol],':',label='%s CREME' % s,color=color)
    
plt.xlabel('MeV/nuc')
plt.ylabel('Flux #/cm$^2$/s/MeV')
plt.subplots_adjust(right=0.75)
plt.legend(loc=(1.04,0))
plt.savefig(os.path.join(figpath,'csda_fluxes.png'))


#%% now convert energy spectrum to LET spectrum
LETgrid = files['dlt']['LET']

fluxLET = SiTable.LET_spectrum(degraded,files['flx']['MeV'],LETgrid,Z=files['flx']['Z'],fast=False)
fluxLETm = SiTable.LET_spectrum(degraded,files['flx']['MeV'],LETgrid,Z=files['flx']['Z'],fast=True)

plt.figure()
plt.loglog(LETgrid,files['dlt']['dflux'],'k-',label='CREME')
plt.loglog(LETgrid,fluxLET,'b-',label='Quadrature')
plt.loglog(LETgrid,fluxLETm,'g--',label='Trapezoidal')
plt.legend()
plt.xlabel('LET in Si, MeV/(g/cm$^2$)')
plt.ylabel('Flux #/cm$^2$/s/(MeV/(g/cm$^2$))')
plt.savefig(os.path.join(figpath,'let_spectra.png'))
           
# read ion parts into dict of dicts
hups = cr.read_report(os.path.join(data_path,'q800a90i100Z2.HUP'))

for d in hups.values(): # d is for "dict"
    # p = dict_to_part(d) # another way to initialize part from CREME dict
    p = IonPartRPP(d['h'],L0=d['L0'],W=d['W'],S=d['S'],slim=d['slim'],aspect=d['aspect'],volumes=d['volumes'])
    rcreme = d['SEE/s']
    
    rs = p.see_rate(fluxLET,LET=LETgrid,fast=False)
    rf = p.see_rate(fluxLET,LET=LETgrid,fast=True)
    rfc = p.see_rate(files['dlt']['dflux'],LET=LETgrid,fast=True)
    print(str(p))
    print('ion rate: SEE/s: slow %g, fast %g, fast-CREME %g, all-CREME %g' % (rs,rf,rfc,rcreme))
    print('ion rate errors fast vs slow: fast %.4f%%' % (100*(rf/rs-1)))
    print('ion rate errors vs CREME: slow %.1f%%, fast %.1f%%, fast-CREME %.1f%%' 
          % (100*(rs/rcreme-1),100*(rf/rcreme-1),100*(rfc/rcreme-1)))

# read proton parts into dict of dicts
pups = cr.read_report(os.path.join(data_path,'q800a90i100.PUP'))

for d in pups.values(): # d is for "dict"
    # p = dict_to_part(d) # part object
    p = ProtonPart(E0=d['E0'],W=d['W'],S=d['S'],slim=d['slim'],volumes=d['volumes'])
    rcreme = d['SEE/s']
               
    icol = np.where(files['flx']['Z'] == 1)[0][0]
    Egrid = files['flx']['MeV'][:,icol]
    pflux = degraded[:,icol]
    rs = p.see_rate(pflux,energy=Egrid,fast=False)
    rf = p.see_rate(pflux,energy=Egrid,fast=True)
    # use CREME degraded flux
    jcol = np.where(files['tfx']['Z']==1)[0][0]
    rfc = p.see_rate(files['tfx']['dflux'][:,jcol],energy=Egrid,fast=True)
    
    print(str(p))
    print('proton rate: SEE/s: slow %g, fast %g, fast-CREME %g, all-CREME %g' % (rs,rf,rfc,rcreme))
    print('proton rate errors fast vs slow: fast %.4f%%' % (100*(rf/rs-1)))
    print('proton rate errors vs CREME: slow %.1f%%, fast %.1f%%, fast-CREME %.1f%%' 
          % (100*(rs/rcreme-1),100*(rf/rcreme-1),100*(rfc/rcreme-1)))
    
plt.show()
