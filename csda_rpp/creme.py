"""
utilities to read CREME96 data files

Public Functions:
    read_flux - reads CREME96 flux files
    read_report - reads CREME96 PUP and HUP report files
"""

import os
import re
import numpy as np

def read_flux(filename):
    """
    data = read_flux(filename)
    read CREME .flx,.tfx,.dlt,.let files
    for flx and tfx files returns dict with fields:
        'MeV/nuc' - energy, MeV/nuc (NE,)
        Z - atomic number (NZ,)
        flux - differential flux (NE,NZ) #/m^s/sr/s/(MeV/nuc)
    for .let files returns dict with fields:
        'LET' - LET MeV/(g/cm^2), (NE,)
        flux - integral flux (NE,)  #/m^s/sr/s
    for .dlt files returns dict with fields:
        'LET' - LET MeV/(g/cm^2), (NE,)
        flux - differential flux (NE,)  #/m^s/sr/s/(MeV/(g/cm^2))
    """
    if not os.path.exists(filename):
        raise Exception('%s does not exist' % filename)
    isLET = filename.lower().endswith('.let') or filename.lower().endswith('.dlt')
    # headers with and without % until a blank line is found
    # states: 'FIRST' skip first line until '%' encountered
    # '%' reading '%' headers
    # 'BLANKS' - reading blank lines
    # 'DATA' - reading data
    state = 'FIRST'
    Emin = None
    Z = 0
    flux = None
    MeV = None
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            if state == 'FIRST':
                if line.startswith('%'):
                    state = '%'
                continue # skip first line
            elif state == '%':
                if not line.startswith('%'):
                    #  1.0000E-01  1.0000E+05   1002      1     92  UNSHIELDED    1977.000  0 2  200 3                    
                    parts = line.split()
                    Emin = float(parts[0])
                    Emax = float(parts[1])
                    NE = int(parts[2])
                    MeV = np.exp(np.linspace(np.log(Emin),np.log(Emax),NE))
                    Zmin = int(parts[3])
                    Zmax = int(parts[4])
                    NZ = Zmax-Zmin+1
                    if isLET:
                        NZ = 1 # ignore Zmax
                    flux = np.full((NE,NZ),np.nan)
                    state = 'BLANKS'
                continue
            elif state == 'BLANKS':
                if len(line) == 0:
                    continue
                Z = Z+1
                irow = 0 # row in flux
                state = 'DATA'
            if (state != 'DATA'):
                raise Exception('Unexpected line in %s :',(filename,line))
            data = [float(x) for x in line.split()]
            flux[irow:(irow+len(data)),Z-1] = data
            irow += len(data)
            if irow >= NE:
                state = 'BLANKS'
    if np.any(np.isnan(flux)):
        raise Exception('file %s read incomplete' % filename)
    if isLET:
        flux = flux.reshape((NE,))
        return {'LET':MeV,'flux':flux}
    else:
        Z = np.arange(Zmin,1+Zmax) # now Z is an array
        return {'Z':Z,'MeV/nuc':MeV,'flux':flux}
            
def read_report(filename):
    """
    data = read_report(filename)
    read CREME .pup and .hup files
    returns dict of parts. Each part is a dict:
    pup (proton) part has fields:
        '#' - report number
        part - part tag, str
        _key - part or part-# if part is not unique within report
        E0 - Weibull energy threshold, MeV
        W - Weibull width, MeV
        S - Weibull power, dimensionless
        slim - Weibull cross section, cm^-2
        volumes - # of sensitive volumes, int
        SEE/vol/s - SEE rate for single sensitive volume (e.g., one bit)
        SEE/vol/day - same, but per day
        SEE/s - SEE rate for system accounting for multiple sensitive volumes
        SEE/day - same, but per day
    hup (ion) part has
        h - thin dimension of part, um
        aspect - aspect ratio of part
        L0 instead of E0 - Weibull energy threshold, MeV/(g/cm^2)
        W - Weibull width, MeV/(g/cm^2)
        otherwise, same as pup part
    """
    
    """
    PUP:
      REPORT NO.    1:   LM7912                                  
      CROSS-SECTION INPUT   4  WEIBULL FIT: 
          ONSET   =    30.000 MeV
          WIDTH   =    10.000 MeV
          POWER   =     1.000 (dimensionless)
          PLATEAU =    10.000 x 10**-12 cm2/bit
      Number of bits =   1.00000E+00
      Rates:      SEEs/bit/second     /bit/day    /device/second    /device/day
      *****   1     1.10564E-09     9.55271E-05     1.10564E-09     9.55271E-05

    HUP:
      REPORT NO.    1:   lm7912                                  
      RPP Dimensions: X =   44.72136  Y =   44.72136  Z =    1.00000 microns.
      Funnel length =     .00000 microns.
      CROSS-SECTION INPUT   4  WEIBULL FIT: 
          ONSET   =      .500 MeV-cm2/milligram
          WIDTH   =    30.000 MeV-cm2/milligram
          POWER   =     1.000 (dimensionless)
          PLATEAU =  2000.000 square microns/bit
      Number of bits =   1.00000E+02
      Rates:      SEEs/bit/second     /bit/day    /device/second    /device/day
      *****   1     1.90340E-08     1.64454E-03     1.90340E-06     1.64454E-01
    
      REPORT NO.    2:   ATMEGA1280                              
      RPP Dimensions: X =   70.71068  Y =   70.71068  Z =    1.00000 microns.
      Funnel length =     .00000 microns.
      CROSS-SECTION INPUT   4  WEIBULL FIT: 
          ONSET   =    10.000 MeV-cm2/milligram
          WIDTH   =   120.000 MeV-cm2/milligram
          POWER   =     1.500 (dimensionless)
          PLATEAU =  5000.000 square microns/bit
      Number of bits =   1.00000E+02
      Rates:      SEEs/bit/second     /bit/day    /device/second    /device/day
      *****   2     4.78966E-11     4.13826E-06     4.78966E-09     4.13826E-04
 
    """
    if not os.path.exists(filename):
        raise Exception('%s does not exist' % filename)
    isPUP = filename.lower().endswith('.pup')
    if isPUP:
        unit_conversion = 1 # MeV->MeV
        slim_conversion = 1e-12 # PLATEAU given as X *10^-12 cm2/bit
        X0key = 'E0' # energy threshold
    else:
        unit_conversion = 1e3 # (MeV/(mg/cm^2)) -> (MeV/(g/cm^2))
        X0key = 'L0' # LET threshold
        slim_conversion = 1e-4**2 # PLATEAU given in um^2
    # headers with and without % until a blank line is found
    # states: 'FIRST' skip first line until '%' encountered
    # '%' reading '%' headers
    # 'BLANKS' - reading blank lines between reports
    # 'DATA' - reading data
    state = 'FIRST'
    reports = {}
    part = None
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            if state == 'FIRST':
                if line.startswith('%'):
                    state = '%'
                continue # skip first line
            elif state == '%':
                if not line.startswith('%'):
                    state = 'BLANKS'
                continue
            elif state == 'BLANKS':
                part = None
                if len(line) == 0:
                    continue
                state = 'DATA'
            if (state != 'DATA'):
                raise Exception('Unexpected line in %s :',(filename,line))
            if len(line) == 0:
                state = 'BLANKS'
                continue
            found = re.search('.*REPORT\s*NO\.\s*(\d+):\s*(\S+)', line)
            if found:
                assert part is None,'Reports ran into each other:\n%s' % line
                part = {'is_proton_part':isPUP}
                part['#'] = int(found.group(1)) # report number
                part['part'] = found.group(2)
                key = part['part']
                if key in reports:
                    key = key + '-' + str(part['#'])                    
                reports[key] = part
                part['_key'] = key
                continue
            
            found = re.search('.*RPP\s+Dimensions:\s*X\s*=\s*(\S+)\s*Y\s*=\s*(\S+)\s*Z\s*=\s*(\S+)', line)
            if found:
                assert not isPUP,'RPP dimensions given for proton part (PUP report):\n%s' % line
                part['h'] = float(found.group(3)) # Z is thickness
                part['aspect'] = float(found.group(1))/float(found.group(2)) # aspect ratio
                # do not retain X,Y separately as a, b because produces possible inconsistency with slim
                if part['aspect']<1: # ensure aspect ratio is >=1
                    part['aspect'] = 1/part['aspect']
                continue
            found = re.search('.*ONSET\s*=\s*(\S+)', line)
            if found:
                part[X0key] = float(found.group(1))*unit_conversion
            found = re.search('.*WIDTH\s*=\s*(\S+)', line)
            if found:
                part['W'] = float(found.group(1))*unit_conversion
                continue
            found = re.search('.*POWER\s*=\s*(\S+)', line)
            if found:
                part['S'] = float(found.group(1))
                continue
            found = re.search('.*PLATEAU\s*=\s*(\S+)', line)
            if found:
                part['slim'] = float(found.group(1))*slim_conversion # -> cm^2
                continue
            found = re.search('.*Number\s+of\s+bits\s*=\s*(\S+)', line)
            if found:
                part['volumes'] = int(float(found.group(1)))
                continue
            found = re.search('\*\*\*\*\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
            if found:
                assert int(found.group(1)) == part['#'], 'Rates rreport # inconsistent with block header:\n%s' % line
                part['SEE/vol/s'] = float(found.group(2))
                part['SEE/vol/day'] = float(found.group(3))
                part['SEE/s'] = float(found.group(4))
                part['SEE/day'] = float(found.group(5))
                continue
            
    return reports

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    _TOP_PATH = os.path.dirname(os.path.abspath(__file__))    
    creme_path = os.path.join(_TOP_PATH,'data')
    
    pups = read_report(os.path.join(creme_path,'q800a90i100.PUP'))
    hups = read_report(os.path.join(creme_path,'q800a90i100Z2.HUP'))
    print('pups',pups)
    print('hups',hups)
    
    flx = read_flux(os.path.join(creme_path,'q800a90i.FLX'))
    print('flx %dx%d MeV/nuc: %e-%e (%d). Z: %d-%d (%d). NaNs:%g%%' % 
          (*flx['flux'].shape,flx['MeV/nuc'][0],flx['MeV/nuc'][-1],len(flx['MeV/nuc']),
           flx['Z'][0],flx['Z'][-1],len(flx['Z']),100*np.mean(np.isnan(flx['flux']))))
    
        
    tfx = read_flux(os.path.join(creme_path,'q800a90i100.TFX'))
    print('tfx %dx%d MeV/nuc: %e-%e (%d). Z: %d-%d (%d). NaNs:%g%%' % 
          (*tfx['flux'].shape,tfx['MeV/nuc'][0],tfx['MeV/nuc'][-1],len(tfx['MeV/nuc']),
           tfx['Z'][0],tfx['Z'][-1],len(tfx['Z']),100*np.mean(np.isnan(tfx['flux']))))

    let = read_flux(os.path.join(creme_path,'q800a90i100.LET'))
    print('let %dx1 MeV/(g/cm^2): %e-%e (%d). NaNs:%g%%' % 
          (*let['flux'].shape,let['LET'][0],let['LET'][-1],len(let['LET']),
           100*np.mean(np.isnan(let['flux']))))

    dlt = read_flux(os.path.join(creme_path,'q800a90i100.DLT'))
    print('let %dx1 MeV/(g/cm^2): %e-%e (%d). NaNs:%g%%' % 
          (*let['flux'].shape,let['LET'][0],let['LET'][-1],len(let['LET']),
           100*np.mean(np.isnan(let['flux']))))
    
    plt.close('all')
    plt.figure()
    plt.loglog(flx['MeV/nuc'],flx['flux'][:,0],'k-',label='Incident H')
    plt.loglog(flx['MeV/nuc'],flx['flux'][:,1],'b-',label='Incident He')
    plt.loglog(flx['MeV/nuc'],flx['flux'][:,2],'r-',label='Incident Li')
    plt.loglog(tfx['MeV/nuc'],tfx['flux'][:,0],'k.--',label='Transmitted H')
    plt.loglog(tfx['MeV/nuc'],tfx['flux'][:,1],'b.--',label='Transmitted He')
    plt.loglog(tfx['MeV/nuc'],tfx['flux'][:,2],'r.--',label='Transmitted Li')
    plt.xlabel('MeV/nuc')
    plt.ylabel('Flux, #/m$^2$/sr/s/(MeV/nuc)')
    plt.legend()

    plt.figure()
    plt.loglog(dlt['LET'],dlt['flux'])
    plt.xlabel('LET: MeV/(g/cm$^2$)')
    plt.ylabel('Differential Flux, #/m$^2$/sr/s/(MeV/(g/cm$^2$))')               
               
    plt.figure()
    plt.loglog(let['LET'],let['flux'])
    plt.xlabel('LET: MeV/(g/cm$^2$)')
    plt.ylabel('Integral Flux, #/m$^2$/sr/s')

    plt.show()
               