"""
bowtie analysis for Response Function Library

see __main__ block at the end for demo (i.e., run this file)

""" 

import numpy as np
import datetime as dt
import rfl

import warnings
warnings.filterwarnings('ignore','invalid value encountered in')
warnings.filterwarnings('ignore','divide by zero encountered in')


def bowtie(inst_info,type,species,channels=None,exponents=[],Ts=[],
           method='mindG',E0='info',tail=None,plot=False):
    """
    results = bowtie(inst_info,type,species,channels=None,exponents=[],Ts=[],
                     method='mindG',E0='info',tail=None,plot=False)
    perform bowtie analysis
    inst_info  - string giving file name or dict containing inst_info
    type - what type of channel to approximate with the bowtie
      'int' for integral channel
      'diff' for differential channel
      'wide' for wide differential channel
          a wide channel requires an E0 (lower threshold) provided by the
          'E0' keyword
    species - string giving species to study (e.g., 'PROT' or 'ELE')
    channels - list of channel names to do (None for all)
    exponents - power law exponents for differential energy spectrum
      e.g., [3,4,5] for E^-3, E^-4, E^-5
    Ts, energy for exponential spectra, exp(-E/T), must have same units as energy in inst_info
    At least two total exponents or Ts must be given
    results - dict with members given by channels or inst_info.CHANNEL_NAMES
      each channel has E0, E1 and G0
      for diff type, E0 is the effective energy for the channel and
        E1 is E0
        G0 is the geometric-efficiency factor times the energy bandwidth (e.g., cm^2 sr MeV)
      for int type, E0 is the effective energy threshold and
        E1 is infinity
        G0 the geometric-efficiency factor (e.g., cm^2 sr)
      for wide type, E0 is the lower energy threshold
        E1 is the effective upper energy limit for the channel and
        G0 is the geometric-efficiency factor (e.g., cm^2 sr)
      other fields:
      type - copy of input type
      fig - figure object (only present if plotting requested)
    errors on E0 and G0 are given as E0err and G0err
      these errors are the standard deviation of the natural log of the E0, G0
      intersections determined in the bow tie analysis (no sqrt(N) factor)
    options:
     method:
        mindG - Choose E0 or E1 that minimizes spread of G0 
            Recommended method for determining E0 or E1
            E0err or E1err will be 0
        intersection - drafted but not tested
    E0:
        list [...] of same length as CHANNEL_NAMES giving E0 in E_UNIT used
        '50%' - determine E0 as the first time the response exceeds 50% max
        '10%' - determine E0 as the first time the response exceeds 10% max
        info' - from the inst_info structure (e.g., inst_info['Elec1']['ELE']['E0']
      (if E0 is omitted, 'info' is assumed, and then '50%' is used if E0 is
      not in the inst_info structure)
    tail - None or dict - {'E':E_MeV,'fun': tail function handle)
       splice high energy tail: j(E) ~tail(E) for E>=E_MeV
       NOTE: this only works for diff or int bowties
       the tail spectrum is differential
    plot - bool - make diagnostic plot
    """
    
    class PowerLaw(object):
        """
        Power-law spectrum with methods flux and G
        """
        def __init__(self,type,n,tail=None):
            self.type = type # channel type
            self.n = n # power law exponent
            self.tail = tail # tail
        def flux(self,E):
            """flux = flux(E)"""
            n = self.n
            tail = self.tail
            flux = E**(-n)
            if self.tail:
                itail = E>tail['E']
                if any(itail):
                    flux[itail] = tail['fun'](E[itail])/tail['fun'](tail['E'])*tail['E']**(-n)
            return flux
        def G(self,c,E,E0,Emax):
            """
            G = .G(c,E,E0,Emax)
            """
            if self.type == 'diff':
                return self.G_diff(c,E,E0,Emax)
            elif self.type == 'int':
                return self.G_int(c,E,E0,Emax)
            elif self.type == 'wide':
                return self.G_wide(c,E,E0,Emax)
            else:
                raise rfl.KeywordError('Channel type %s unknown' % type)
        def G_diff(self,c,E,E0,Emax):
            """G = G_diff(c,E,E0,Emax)"""
            n = self.n
            tail = self.tail
            G = c*E**n
            if tail is not None:
                iE = E>tail['E']
                if any(iE):
                    G[iE] = c/tail['fun'](E[iE])*tail['E']**n*tail['fun'](tail['E'])
            return G
        def G_int(self,c,E,E0,Emax):
            """G = G_int(c,E,E0,Emax)"""
            n = self.n
            tail = self.tail
            # jdiff = E^-n
            # G = c/Jint
            Elim = Emax
            if tail:
                if tail['E']>Emax:
                    tail = None
                else:
                    Elim = tail['E'] # Elim serves as Emax
            # integrals up to Elim~Emax or tail.E
            if n<1:
                # G = c*(n-1)./(Emax.^(n-1)-E.^(n-1));
                Jint = (Elim**(n-1)-E**(n-1))/(n-1)
            elif n==1: 
                #G = c./(log(Emax)-log(E));
                Jint = np.log(Elim)-np.log(E)
            else: 
                # G = c*E.^(n-1)*(n-1);
                Jint = (E**(1-n)-Elim**(1-n))/(n-1)
            # now do tail integral if needed
            if tail:
                from scipy.integrate import quad
                J0int = Jint
                jtail = tail['fun'](tail['E'])
                j0Etail = tail['E']**-n
                itail = E>tail['E']
                if np.any(~itail):
                    integral,abserr = quad(tail['fun'],tail['E'],Emax)
                    JprimeEtail = j0Etail/jtail*integral
                    Jint[~itail] = J0int[~itail]+JprimeEtail;
                if np.any(itail):
                    for i in np.where(itail)[0]:
                        integral,abserr = quad(tail['fun'],E[i],Emax)
                        Jint[i] = j0Etail/jtail*integral
            
            G = c/Jint;
            return G
        def G_wide(self,c,E,E0,Emax):
            """G = G_wide(c,E,E0,Emax)"""
            n = self.n
            tail = self.tail
            if n==1: # need to include Emax and special caase n
                G = c/np.log(E/E0)
            else: # treat Emax as infinity
                G = c*(n-1)/(E0**(1-n)-E**(1-n)) # positive if E>E0, n>1
            if tail:
                raise NotImplementedError
            return G
            
    class Exponential(PowerLaw):
        """
        Exponential spectrum with methods flux and G
        """
        def __init__(self,type,T,tail=None):
            self.type = type # channel type
            self.T = T # constant in exp(-E/T)
            self.tail = tail # tail
        def flux(self,E):
            """flux = flux(E)"""
            T = self.T
            tail = self.tail
            flux = np.exp(-E/T)
            if tail:
                itail = E>tail['E']
                if any(itail):
                    flux[itail] = tail['fun'](E[itail])/tail['fun'](tail['E'])*np.exp(-tail['E']/T)
            return flux
        def G_diff(self,c,E,E0,Emax):
            """G = G_diff(c,E,E0,Emax)"""
            T = self.T
            tail = self.tail
            G = c*np.exp(E/T)
            if tail:
                iE = E>tail['E']
                if any(iE):
                    G[iE] = c/tail['fun'](E[iE])*np.exp(tail['E']/T)*tail['fun'](tail['E'])
            return G
        def G_int(self,c,E,E0,Emax):
            """G = G_int(c,E,E0,Emax)"""
            T = self.T
            tail = self.tail
            # jdiff = exp(-E/T)
            # G = c/Jint
            # G = c.*exp(E./T)./T;
            Jint = T*np.exp(-E/T);
            # now do tail integral if needed
            if tail:
                from scipy.integrate import quad
                J0int = Jint;
                jtail = tail['fun'](tail['E'])
                j0Etail = np.exp(-tail['E']/T)
                itail = E>tail.E
                if np.any(~itail):
                    integral,abserr = quad(tail['fun'],tail['E'],Emax)
                    JprimeEtail = j0Etail/jtail*integral
                    Jint[~itail] = J0int[~itail]+JprimeEtail
                if np.any(itail):
                    for i in np.where(itail)[0]:
                        integral,abserr = quad(tail['fun'],E[i],Emax)
                        Jint[i] = j0Etail/jtail*integral
            
            G = c/Jint;
            return G
        def G_wide(self,c,E,E0,Emax):
            """G = G_wide(c,E,E0,Emax)"""
            T = self.T
            tail = self.tail
            G = c/T/(np.exp(-E0/T)-np.exp(-E/T))
            if tail:
                raise NotImplementedError
            return G
    
    methods = ['mindg','intersection'] # supported methods
    method = method.lower()
    assert method in methods,'method %s not supported' % method
    
    types = ['diff','wide','int']
    type = type.lower()
    assert type in types,'type %s not supported' % type
    
    if channels is None:
        channels = inst_info['CHANNEL_NAMES']
    
    assert len(exponents)+len(Ts)>=2,'Must specify at least 2 test spectra'
    
    if tail is not None:
        assert method == 'mindg','tail option only supported for method mindG'
        assert type in ['diff','int'],'tail option only supported for type int or diff'
    
    E0mode = E0 # better variable name, E0 will be used as a varaible, not a control switch
    sinfo = {}
    sinfo['PL'] = {'type':'PL','linestyle':'k-','special':False}
    sinfo['EXP'] = {'type':'EXP','linestyle':'k--','special':False}
    
    Nn = len(exponents)
    NT = len(Ts)
    Ns = NT+Nn;
    spectra = []
    for n in exponents:
        spec = {**sinfo['PL'],'param':n}
        spectra.append(spec)
        spec['SpecObj'] = PowerLaw(type,n,tail)
        spec['special'] = (n<=1) or (type=='wide')
    for T in Ts:
        spec = {**sinfo['EXP'],'param':T}
        spectra.append(spec)
        spec['SpecObj'] = Exponential(type,T,tail)
    
    # build intersect function table
    if type == 'int':
        intersect_func = {
            ('PL','PL'):intersect_INT_PL_PL,
            ('PL','EXP'):intersect_INT_PL_EXP,
            ('EXP','PL'):intersect_INT_EXP_PL,
            ('EXP','EXP'):intersect_INT_EXP_EXP
        }
    elif type == 'diff':
        intersect_func = {
            ('PL','PL'):intersect_DIFF_PL_PL,
            ('PL','EXP'):intersect_DIFF_PL_EXP,
            ('EXP','PL'):intersect_DIFF_EXP_PL,
            ('EXP','EXP'):intersect_DIFF_EXP_EXP
        }
    elif type == 'wide':
        intersect_func = {
            ('PL','PL'):intersect_WIDE_PL_PL,
            ('PL','EXP'):intersect_WIDE_PL_EXP,
            ('EXP','PL'):intersect_WIDE_EXP_PL,
            ('EXP','EXP'):intersect_WIDE_EXP_EXP
        }
    
    inst_info = rfl.load_inst_info(inst_info)
    
    if isinstance(E0mode,list):
        assert len(E0mode) == len(channels),'Supplied E0 list must have one entry per channel (%d)' % len(channels)
    
    results = {}
    for ichan,chan in enumerate(channels):
        assert chan in inst_info,'Channel %s not found in inst_info' % chan
        results[chan] = {'G0':np.NaN,'E0':np.NaN,'E1':np.NaN,'type':type}
        assert species in inst_info[chan]['SPECIES'],'Species %s not found in %s' % (species,chan)
        sp = species # shorthand        
        resp = inst_info[chan][sp]
        hE = resp.hE(resp.E_GRID)
        R = hE/rfl.make_deltas(resp.E_GRID)
        if type == 'wide':
            if isinstance(E0mode,list):
                E0 = E0mode[ichan]
            elif isinstance(E0mode,str):
                E0 = E0mode
                if E0 == 'info':
                    E0 = resp.get('E0','50%') # default to 50%
                if E0 == '50%':
                    E0 = resp.E_GRID[R>=np.max(R)/2][0]
                elif E0 == '10%':
                    E0 = resp.E_GRID[R>=np.max(R)/10][0]
                else:
                    raise rfl.KeywordError('Unable to resolve E0 %s for %s' % (E0mode,chan))
            else:
                raise rfl.KeywordError('Unable to resolve E0 %s for %s' % (str(E0mode),chan))
            E_GRID = resp.E_GRID[resp.E_GRID>=E0] # limit E_GRID for E1 search to > E0
        else:
            E0 = np.nan
            E_GRID = resp.E_GRID
        Emax = E_GRID[-1]
        
        if plot:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            results[chan]['fig'] = fig
        
        C = []
        # find reference counts
        for spec in spectra:
            C.append((hE.ravel()*spec['SpecObj'].flux(resp.E_GRID.ravel())).sum())
        # find intersections
        intersect = [];
        if method == 'intersection':
            for ipass in [0,1]:
                # pass 0: do same type
                # pass 1: do hybrid type w/ guess since these can have multiple crossings
                if ipass == 0:
                    Eguess = np.nan
                else:
                    Eguess = np.median(intersect[:,0]) # initial guess from same-type intersections

                for is1 in range(Ns):
                    s1 = spectra[is1]
                    c1 = C[is1]
                    if plot and (ipass==0):
                        plt.figure(fig)
                        plt.loglog(E_GRID,s1['SpecObj'].G(c1,E_GRID,E0,Emax),s1['linestyle'],label='%s: %g' % (spectra[is1]['type'],spectra[is1]['param']))

                    for is2 in range(is1+1,Ns):
                        s2 = spectra[is2]
                        ishybrid = (s1['type']==s2['type'])
                        if (ipass==0) == (ishybrid or s1['special'] or s2['special']): # do this on pass 0 for same types, on pass 1 for hybrid types
                            c2 = C[is2]
                            E = intersect_func[s1['type'],s2['type']](s1,s2,c1,c2,E0,Emax,Eguess)
                            G0 = s1.G(c1,E,E0,Emax)
                            if (E<=0) or (not np.isreal(E)) or (not np.isfinite(E)) or (G0<=0) or (not np.isreal(G0)) or (not np.isfinite(G0)):
                                print('No real, positive, finite fit found for %s/%g - %s/%g' % (s1.type,s1.param,s2.type,s2.param))
                                continue;
                            intersect = [intersect,[E,G0,is1,is2]]
        elif method == 'mindg':
            G0 = np.full((len(E_GRID),Ns),np.nan)
            for js in range(Ns):
                G0[:,js] = spectra[js]['SpecObj'].G(C[js],E_GRID,E0,Emax)
                if plot:
                    plt.figure(fig)
                    plt.loglog(E_GRID,G0[:,js],spectra[js]['linestyle'],label='%s: %g' % (spectra[js]['type'],spectra[js]['param']))
            # find best guess
            stdG = np.std(np.log(G0),axis=1)
            stdG[~np.isfinite(stdG)] = np.inf #ignore nans in argmin
            imin = np.argmin(stdG)
            stdG = stdG[imin]
            intersect = np.full((Ns,4),np.nan)
            intersect[:,0] = E_GRID[imin]
            intersect[:,1] = G0[imin,:]
            intersect[:,2] = range(Ns)
            intersect[:,3] = range(Ns)
        if plot:
            plt.figure(fig);
            plt.loglog(intersect[:,0],intersect[:,1],'bo')

        solution = np.exp(np.median(np.log(intersect[:,0:2]),axis=0))
        sol_error = np.std(np.log(intersect[:,0:2]),axis=0)
        results[chan]['G0'] = solution[1]
        results[chan]['G0err'] = sol_error[1]
        if type == 'int':
            Evar = 'E0';
            results[chan]['E0'] = solution[0]
            results[chan]['E0err'] = sol_error[0]
            results[chan]['E1'] = np.inf
            results[chan]['E1err'] = 0
        elif type == 'diff':
            Evar = 'E0';
            results[chan]['E0'] = solution[0]
            results[chan]['E0err'] = sol_error[0]
            results[chan]['E1'] = results[chan]['E0']
            results[chan]['E1err'] = results[chan]['E0err']
        elif type == 'wide':
            Evar = 'E1';
            results[chan]['E1'] = solution[0]
            results[chan]['E1err'] = sol_error[0]
            results[chan]['E0'] = E0
            results[chan]['E0err'] = 0
        
        if plot:
            plt.figure(fig);
            if any(R.flatten()>0):
                plt.loglog(solution[0],solution[1],'rx',lw=3)
                axmin = np.min(intersect[:,0:2],axis=0)
                axmax = np.max(intersect[:,0:2],axis=0)
                axmin = np.maximum(axmin,E_GRID[0])
                axmax = np.minimum(axmax,E_GRID[-1])
                axmin[0] = np.minimum(axmin[0],E_GRID[0])
                axmax[0] = np.maximum(axmax[0],E_GRID[-1])
                axmin[1] = np.minimum(axmin[1],np.min(R[R>0]))
                axmax[1] = np.maximum(axmax[1],np.max(R))
                axmin = 10**np.floor(-0.5+np.log10(axmin))
                axmax = 10**np.ceil(0.5+np.log10(axmax))
                plt.loglog(resp.E_GRID,np.maximum(R,axmin[1]),'g-',lw=2) # epsdEG vs E
                if not np.isfinite(axmax[1]):
                    print(axmin,axmax)
                plt.axis([axmin[0],axmax[0],axmin[1],axmax[1]])
            else:
                plt.axis([0,1,0,1])
                plt.text(0.5,0.5,'No positive entries in response',ha='center',va='center')
            if type in ['int','wide']:
                Gsym = 'G0'
                G_UNIT = '%s$^2$sr' % resp.L_UNIT
            elif type == 'diff':
                Gsym = 'G0dE'
                G_UNIT = '%s%s$^2$sr' % (resp.E_UNIT,resp.L_UNIT)
            plt.title('%s,%s : %s=%g (%.1f%%), %s=%g (%.1f%%)' % (
                chan,type,Evar,solution[0],sol_error[0]*100,
                Gsym,results[chan]['G0'],results[chan]['G0err']*100));
            plt.xlabel('%s, %s' %(Evar,resp.E_UNIT))
            plt.ylabel('%s, %s'%(Gsym,G_UNIT))
            plt.grid(True)
    return results
# intersect functions all have same syntax:
# input two spectra and two counts
# and E0 which is often ignored
# Eguess is also often ignored
def intersect_INT_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_INT_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess)"""
    n1 = s1['param']
    n2 = s2['param']
    if n1 > n2: # try again with n1 < n2
        return intersect_INT_PL_PL(s2,s1,c2,c1,E0,Emax,Eguess);

    if n1<1: # have to deal with Emax
        # non-analytic cases
        # n2 < 1 cases
        # G0 = c1*(n1-1)/(E^(n1-1)-Emax^(n1-1))
        # G0 = c2*(n2-1)/(E^(n2-1)-Emax^(n2-1))
        # c1*(n1-1)/(E^(n1-1)-Emax(n2-1)) = c2*(n2-1)/(E^(n2-1)-Emax^(n1-1))
        # n2 = 1, n2>1 cases also not analytic
        if not np.isfinite(Eguess):
            Eguess = Emax*0.95
        E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess) # non-analytic
    elif n1 == 1: # n1 = 1, n2 > 1        
        # G0 = c1/(log(Emax)-log(E))
        # G0 = c2*E0^(n2-1)*(n2-1)
        # c1*E0^(n1-1)*(n1-1) = c1/(log(Emax)-log(E))
        # no analytical solution
        if not np.isfinite(Eguess): # guess at n=1.01 or half-way to n2
            n1 = np.minimum(1.01,(n1+n2)/2)
            Eguess = (c1/c2*(n1-1)/(n2-1))**(1/(n2-n1)) # E0
        E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess) # non-analytic
    else:
        # G0 = c1*E0^(n1-1)*(n1-1)
        # G0 = c2*E0^(n2-1)*(n2-1)
        # c1*E0^(n1-1)*(n1-1) = c2*E0^(n2-1)*(n2-1)
        # c1/c2*(n1-1)/(n2-1) = E0^(n2-n1)
        # E0 = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1))
        E = (c1/c2*(n1-1)/(n2-1))**(1/(n2-n1)) # E0
    return E

def intersect_INT_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_INT_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)"""
    # c2/T2*exp(E/T2) = (n1-1)*c1*E^(n1-1)
    return intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)

def intersect_INT_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_INT_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess)"""
    return intersect_INT_PL_EXP(s2,s1,c2,c1,E0,Emax,Eguess)

def intersect_INT_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_INT_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess)"""
    T1 = s1['param']
    T2 = s2['param']
    return (np.log(c2)-np.log(T2)-np.log(c1)+np.log(T1))/(1/T1-1/T2)

def intersect_DIFF_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_DIFF_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess)"""
    n1 = s1['param']
    n2 = s2['param']
    # G0 = c1*E0^n1
    # G0 = c2*E0^n2
    # c1*E0^n1 = c2*E0^n2
    # c1/c2 = E0^(n2-n1)
    # E0 = (c1/c2)^(1/(n2-n1))
    return (c1/c2)**(1/(n2-n1)) # E0

def intersect_DIFF_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_DIFF_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)"""
    # log(c1)-log(c2) = E/T2-n1*log(E)
    return intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)

def intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)"""
    from scipy.optimize import fsolve
    func = lambda logE: np.log(s1['SpecObj'].G(c1,np.exp(logE),E0,Emax))-np.log(s2['SpecObj'].G(c2,np.exp(logE),E0,Emax))
    logE,infodict,ier,mesg = fsolve(func,np.log(Eguess),full_output=True)
    if ier==1:
        E = np.exp(logE[0])
    else:
        E = np.nan
    return E
    
def intersect_DIFF_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_DIFF_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess)"""
    return intersect_DIFF_PL_EXP(s2,s1,c2,c1)

def intersect_DIFF_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_DIFF_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess)"""
    T1 = s1['param']
    T2 = s2['param']
    return (np.log(c1)-np.log(c2))/(1/T2-1/T1)

def intersect_WIDE_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_WIDE_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess)"""
    # G0 = c1*(n1-1)/[E0^(1-n1)-E1^(1-n1)]
    # G0 = c2*(n2-1)/[E0^(1-n2)-E1^(1-n2)]
    # c1*(n1-1)/[E0^(1-n1)-E1^(1-n1)] = c2*(n2-1)/[E0^(1-n2)-E1^(1-n2)]
    # no analytical solution, solve numerically
    return intersect_WIDE(s1,s2,c1,c2,E0,Emax,Eguess)

def intersect_WIDE_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_WIDE_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)"""
    # no analytical solution
    return intersect_WIDE(s2,s1,c2,c1,E0,Emax,Eguess)

def intersect_WIDE(s1,s2,c1,c2,E0,Emax,Eguess): # generic
    """E = intersect_WIDE(s1,s2,c1,c2,E0,Emax,Eguess) generic"""
    if np.isfinite(Eguess):
        logdE = np.log(Eguess-E0)
    else:
        logdE = np.log(E0)-2

    from scipy.optimize import fsolve
    func = lambda logE: np.log(s1['SpecObj'].G(c1,E0+np.exp(logdE),E0,Emax))-np.log(s2['SpecObj'].G(c2,E0+np.exp(logdE),E0,Emax))
    logE,infodict,ier,mesg = fsolve(func,np.log(logdE),full_output=True)

    if ier==1:
        E = E0+np.exp(logdE[0])
    else:
        E = np.nan
    return E

def intersect_WIDE_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_WIDE_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess)"""
    return intersect_WIDE_PL_EXP(s2,s1,c2,c1,E0,Emax,Eguess)

def intersect_WIDE_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess):
    """E = intersect_WIDE_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess)"""
    # c1*(E-E0)/T1/(exp(-E0/T1)-exp(-E/T1)) = c2*(E-E0)/T2/(exp(-E0/T2)-exp(-E/T2));
    return intersect_WIDE(s1,s2,c1,c2,E0,Emax,Eguess)



if __name__ == '__main__':

    # uncomment these when developing/debugging
    # without it, too many figures get created if we re-run it in the same console
    # import matplotlib.pyplot as plt
    # plt.close('all')

    # set up energy grid
    MeV = 10**np.linspace(-1,np.log10(30),1000)
    
    # top-level sensor (instrument) info
    inst_info = {
        'FORMAT_VERSION' : '1.1.0',
        'L_UNIT' : 'cm',
        'E_UNIT' : 'MeV',
        'REFERENCES' : ['Created %s' % dt.datetime.utcnow().isoformat()],
        'DEAD_TIME_PER_COUNT' : 0,
        'DEAD_TYPE' : 'BLOCKING',
        'COUNTS_MAX' : np.inf, # the actual COUNTS_MAX is something else, but I don't know what.
        'SPECIES' : ['ELE'],    
    }
    
    # electron stuff
    inst_info['ELE'] = {
        'RESP_TYPE' : '[E]', # omni with efficiency table
        'E_TYPE' : 'TBL',
        'BIDIRECTIONAL' : 'FALSE',
        'E_GRID' : MeV
    }
    
    # Channel-specific info
    inst_info['CHANNEL_NAMES'] = []
    wideE0s = [] # holds wide channel lower threshold
    wideE1s = [] # holds wide channel upper threshold
    for iE in [1,3,9]:
        for dE in [0.1,0.5,1]:
            chan = 'E%dMEV_d%dkeV' % (iE,dE*1e3)
            inst_info['CHANNEL_NAMES'].append(chan)
            inst_info[chan] = {}
            inst_info[chan]['ELE'] = {
                'G' : 1,
                'E0' : iE,
                'DE' : dE,
                'CROSSCALIB' : 1,
                'CROSSCALIB_RMSE' : np.log(2)/2,
                'EPS' : np.array(np.abs(MeV-iE) < dE/2,dtype=float)
            }
            wideE0s.append(iE-dE/2)
            wideE1s.append(iE+dE/2)
        chan = 'GT%dMEV' % iE
        inst_info['CHANNEL_NAMES'].append(chan)
        inst_info[chan] = {}
        inst_info[chan]['ELE'] = {
            'G' : 1,
            'E0' : iE,
            'DE' : np.inf,
            'CROSSCALIB' : 1,
            'CROSSCALIB_RMSE' : np.log(2)/2,
            'EPS' : np.array(MeV>=iE,dtype=float)
        }
        wideE0s.append(iE)
        wideE1s.append(np.inf)
    
    inst_info = rfl.load_inst_info(inst_info)
    
    
    # diff fits
    results = bowtie(inst_info,'diff','ELE',exponents=[2,3,4,5],plot=True,Ts=[0.1,0.5,1])
    for chan in inst_info['CHANNEL_NAMES']:
        ideal = inst_info[chan]['ELE']
        btie = results[chan]
        print('diff %15s, ideal(bowtie): E0=%.2f(%.2f), dE = %.2f(%.2f)' % (
            chan,ideal.E0,btie['E0'],ideal.DE,btie['G0']/ideal.G))
    
    # wide fits
    results = bowtie(inst_info,'wide','ELE',exponents = list(range(2,6)),E0='50%',plot=True,Ts=[0.1,0.5,1])
    for i,chan in enumerate(inst_info['CHANNEL_NAMES']):
        ideal = inst_info[chan]['ELE']
        btie = results[chan]
        print('wide %15s, ideal(bowtie): E0=%.2f(%.2f), E1=%.2f(%.2f), G0 = %.2f(%.2f)' % (
            chan,wideE0s[i],btie['E0'],wideE1s[i],btie['E1'],ideal.G,btie['G0']))

    # int fits
    results = bowtie(inst_info,'int','ELE',exponents=list(range(2,6)),plot=True,Ts=[0.1,0.5,1])
    for chan in inst_info['CHANNEL_NAMES']:
        ideal = inst_info[chan]['ELE']
        btie = results[chan]
        print('int %15s, ideal(bowtie): E0=%.2f(%.2f), G0 = %.2f(%.2f)' % (        
            chan,ideal.E0,btie['E0'],ideal.G,btie['G0']))

