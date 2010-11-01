#!/usr/lib/env python
# -*- coding: utf-8 -*-

"""Python interface to the IRBEM InvLib library

Calls to ana_spec_inv and omni2uni are currently supported
Angular inversions and spectral inversions are performed using
different classes. It has not yet been implemented, but the plan
is to have a AngInv object as a possible input to a SpecInv
object, such that the fluxes, etc are easy to input for the 
spectral inversion.
Running this wrapper from the command line calls a test suite to 
reproduce the specinv_test and omni2uni_test where the setup and 
function calls are from within Python. The test suite also tests
that the class is working correctly and should be run for both 
Python 2.6+ and 3.X whenever changes are made.

Author: Steve Morley, Los Alamos National Laboratory (smorley@lanl.gov)
Date Created: 23 Sept. 2010

To Do:
Enable read support for PRBEM response files
Generalize setting of _real_params
Add support for calling pc_spec_inv
Improve error trapping
Update unit tests for additional code
Documentation
Add support for calling of *spec_inv_multi (within extant functions)
Angular Inversion code (add wide2uni)
"""

import ctypes, os, sys, csv, math, numbers
import numpy as np

if sys.platform == 'linux2':
    ext = 'so'
elif sys.platform == 'darwin':
    ext = 'dylib' #although invlib doesn't currently compile on mac...
else:
    raise NotImplementedError('Your platform is not currently supported')

floc = os.path.split(locals()['__file__'])
libpath = floc[0]

##following two functions only used in invlib_test.py
def ravel(mylist):
    """Pure Python implementation of numpy ravel method"""
    newlist = []
    for row in mylist:
        # may have nested tuples or lists
        if type(row) in (list, tuple):
            newlist.extend(ravel(row))
        else:
            newlist.append(row)
    return newlist

def transpose(arr):
    """Pure Python implementation of numpy transpose method"""
    return [[r[col] for r in arr] for col in range(len(arr[0]))]
#----------------------


class InvBase(object):
    """
    Base class for INVLIB calls. Not to be used directly.
    
    This base class is here to provide shared functionality for the SpecInv
    and AngInv classes, which do spectral inversion and angular inversion, 
    respectively. For help on using these, please see the help for those 
    classes.
    """
    def __init__(self, *args, **kwargs):
        try:
            self._invlib = ctypes.CDLL(os.path.join(libpath, 'invlib.'+ext))
        except: #case for testing in IRBEM dist
            self._invlib = ctypes.CDLL('../invlib.'+ext)
        if 'verbose' in kwargs:
            if kwargs['verbose'] == True:
                self._verb = 1
            elif kwargs['verbose'] == False:
                self._verb = 0
            else:
                self._verb = kwargs['verbose'] 
        else:
            self._verb = 1
            
        self.H = None
        
        return None
        
    def __repr__(self):
        attrs = [att for att in self.__dict__ if '_' not in att]
        attrs.sort()
        
        clstr = str(self.__class__).split('.')[-1].split('\'')[0]
        
        outstr = clstr + ' instance: ' + '\n'
        for ent in attrs:
            try:
                le = '[' + str(len(getattr(self, ent))) + ']\n'
                outstr = outstr + str(ent) + ' :  ' + le
            except TypeError:
                le = str(type(getattr(self, ent))) + '\n'
                outstr = outstr + str(ent) + ' :  ' + le
                
        return outstr
            
    def _readLANLResp(self, fname, **kwargs):
        """read 'standard' LANL format response files
        
        The MCNP-derived response files contain G(E), that is, the geometric
        factor as a function of energy. More specifically, this is the 
        effective geometric factor (taking into account efficiencies, etc.)
        from Monte Carlo N-Particle simulation (see e.g. Tuszewski et al., 
        Nucl. Instrum. Methods Phys. Res. A, 2002) [cm^2 sr].
        
        To get H in the units required by *_spec_inv [cm^2 sr s keV] the 
        response function must be multiplied by the accumulation interval,
        assuming the counts are per accumulation. To include the energy units
        the energy integral weighting setting in INVLIBs int_params array must 
        be set. Here the trapezoidal integration (int_params[7]=1) should work
        """
        fobj = open(fname, 'r')
        header, data, Earray = [], [], []
        c_symb = ';'
        for line in fobj:
            if c_symb in line: #if prefixed by ';' it's header
                header.append(line.rstrip())
            else:
                if line.rstrip()!='': #skip blank lines
                    dum = line.rstrip().split()
                    data.append(dum[1:])
                    try:
                        Earray.append(float(dum[0])) #get energies from first column
                    except ValueError:
                        pass #skip text input
        
        #move first line to column header list
        col_hdr = data.pop(0)
        
        #set attributes
        if 'dt' in kwargs:
            dt = kwargs['dt']
        else:
            #default to 1-sec for now
            dt = 1
        print('Using accumulation time %g' % dt)
        Htmp = np.array(data, dtype=float).transpose().ravel()*dt
        self.H = Htmp.tolist()
        self.Egrid = Earray
        self.Hcols = col_hdr
        
        if self._verb:
            print('LANL response file read:\n%s' % header[:2])
        try:
            assert self.Eout#==self.Egrid
        except AttributeError:
            print('Energy grid for output not set: Using grid from response function')
            self.Eout = self.Egrid

        return None
        
    def readRespFunc(self, fname=None, std='LANL', dt=1):
        """method for reading PRBEM and LANL format response functions
        
        *_spec_inv requires that response is multiplied by the accumulation
        time, dt. *_spec_inv_multi requires that response is unaltered but 
        takes dt as a parameter and performs the multiplication itself
        
        As *spec_inv_multi is not yet implemented, for now dt is passed as a 
        keyword to the response function read method.
        """
        #perhaps set as single or multi on object instantiation so only 
        #appropriate method is available
        
        if std.upper() not in ['LANL', 'PRBEM']:
            raise ValueError('Unknown response file format specified')
        if fname == None or not os.path.exists(fname):
            raise IOError('Requested file %s does not exist' % (fname))
        if std.upper() == 'LANL':
            #default dt is 1s; 10.24s is spin period of SOPA
            self._readLANLResp(fname, dt=dt)
            #formatting currently hacked into readLANLResp
        elif std.upper() == 'PRBEM':
            try:
                from spacepy import pycdf
                rfun = pycdf.CDF(fname)
                #now get appropriate numbers from CDF to pass to *_spec_inv
                #then set self.H
            except ImportError:
                raise ImportError('''SpacePy is not installed;
                CDF support is unavailable without SpacePy\n
                Please get SpacePy from http://spacepy.lanl.gov''')
        
        return None
        
    
class SpecInv(InvBase):
    """Spectral Inversion class using INVLIB
    
    This class is designed to hold all necessary inputs to *_spec_inv* from
    the INVLIB library. The calls to INVLIB routines are implemented as
    object methods. 
    
    Example use:
    import package (required)
    >>> import pyinvlib as pinv
    instantiate and populate attributes
    >>> dum = pinv.SpecInv(verbose=1)
    >>> dum.counts = [4.15524e+02, 3.70161e+02, 2.42137e+02, 2.12097e+02, \
        1.47379e+02, 1.40524e+02, 9.37500e+01, 5.08064e+01, 1.93548e+01, 3.22581e+00]
    >>> dum.dcounts=[0.3466]*10
    >>> dum.readRespFunc(fname='../../projects/responses/sopa/mcp_output/version_3/electrons/sopa_LANL-97A_t2_HSP_elec.001')
    >>> dum.setParams(fnc=8)
    >>> dum.anaSpecInv()
    plot the output spectrum with errorbars (using matplotlib)
    >>> fig = dum.plot()
    """
    def __init__(self, *args, **kwargs):
        super(SpecInv, self).__init__(self, *args, **kwargs)
        #set some defaults
        self._int_params = [None]*10
        self._real_params = [None]*10
        if 'rme' in kwargs:
            self.rme = kwargs['rme'] 
        else:
            self.rme = 0.511 #defaults to electron with all units in MeV

    def readTestInput(self, verbose=True, func=1+2, minim=0, niter=1000):
        """reads test data for specinv_test.c
        
        To perform the test carried out by specinv_test.c, instantiate a SpecInv object
        and call this method. Then call the anaSpecInv method -- this is all done from
        the testing suite as one of the unit tests as they have the expected output to 
        compare against.
        """
        self._verb = verbose
        try:
            fnamein = os.path.join(libpath,"specinv_test.in1")
            assert os.path.exists(fnamein)
        except AssertionError:
            try: #case for testing in IRBEM dist
                fnamein = os.path.join("..","specinv_test.in1")
                print('Trying %s' % fnamein)
                assert os.path.exists(fnamein)
            except:
                raise IOError('Input test file not found')
        infile = open(fnamein, 'r')
        
        #read file
        deflen = infile.readline().rstrip().split() #1st line, newline removed, split on whitespace
        NC, NE, NEout = int(deflen[0]), int(deflen[1]), int(deflen[2])

        #populate vars
        #remove all trailing \n and store each non-empty line as a string in a list
        datlist = [float(d.rstrip()) for d in infile if len(d)>2]
        self.counts = datlist[:NC]
        self.dcounts = datlist[NC:NC*2]
        self.Egrid = datlist[NC*2:(NC*2)+NE]
        dum = ((NC*2)+NE) #set counter for readability
        self.H = datlist[dum:dum+(NE*NC)]
        dum = dum+(NE*NC)
        self.b = datlist[dum:dum+NC]
        self.Eout = datlist[dum+NC:]
        print('Python Analytic Spectral Inversion test: NC=%d, NE=%d, NEout=%d' % (NC, NE, NEout))
        
        self._int_params[0] = NC
        self._int_params[1] = NE
        self._int_params[2] = NEout
        self._int_params[5] = niter
        self._int_params[6] = self._verb
        self.setParams(fnc=func, Hint=0, minim=minim)
        
        return None
        
    def setParams(self, fittype='ana', fnc=None, Hint=1, minim=None):
        """Generic method for setting parameter inputs to *_spec*
        
        Checks that the object is populated with sufficient data and fills in 
        default values where necessary.
        """
        try:
            assert self.counts
            assert len(self.dcounts)==len(self.counts)
            if not self._int_params[0]:
                self._int_params[0] = len(self.Hcols) #how to calc number of channels? NC
            if not self._int_params[1]:
                self._int_params[1] = len(self.Egrid)  #NE
            #if not self._int_params[2]:
            self._int_params[2] = len(self.Eout)
        except (AssertionError, AttributeError):
            raise TypeError('Counts, Relative Error on Counts and Energy grid must be defined') 
        try:
            assert self.H
        except:
            raise AttributeError('Response functions must be specified')
        try:
            assert len(self.b) == len(self.counts)
        except:
            print('Background not specified, defaulting to zero')
            self.b = [0]*len(self.counts)
        
        #define outputs
        NC, NE = int(self._int_params[0]), int(self._int_params[1])
        NEout = int(self._int_params[2])
        try:
            c, dc = NC*ctypes.c_double, NC*ctypes.c_double
            Egrid = NE*ctypes.c_double
            H = NE*NC*ctypes.c_double
            Eout = NEout*ctypes.c_double
            b = NC*ctypes.c_double
            flux = (NEout * ctypes.c_double)(0)
            dlogflux = (NEout * ctypes.c_double)(0)
            Eout = Eout(*self.Eout)
        except IndexError:
            raise IndexError('Error: Mismatched counts/energy grid')
        NEout = int(len(Eout))
        
        flux = (NEout * ctypes.c_double)(0)
        dlogflux = (NEout * ctypes.c_double)(0)
        
        #check for absence of keywords and set defaults
        if fnc == None:
            self._int_params[3] = 1+2 #defaults to power law and exponential
        else:
            self._int_params[3] = fnc
        if minim == None:
            self._int_params[4] = 0 #default to BFGS minimizer
        else:
            self._int_params[4] = minim 
        if self._int_params[5] == None: self._int_params[5] = 5000
        if self._int_params[6] == None: self._int_params[6] = self._verb
        if self._int_params[7] == None: self._int_params[7] = Hint
        if self._int_params[8] == None: self._int_params[8] = 0
        if self._int_params[9] == None: self._int_params[9] = 0
        
        intp = ctypes.c_long*10
        realp = ctypes.c_double*10
        if fittype.lower() == 'ana':
            intp = intp(*self._int_params)
            realp = realp(*[self.rme, 100, 345, 0, 0, 0, 0, 0, 0, 0])
        elif fittype.lower() == 'pc':
            raise NotImplementedError('Principal Component method not yet implemented')
        
        counts = c(*self.counts)
        dcounts = dc(*self.dcounts)
        Egrid = Egrid(*self.Egrid)
        H = H(*self.H)
        b = b(*self.b)
        
        self._params = (counts, dcounts, Egrid, 
            H, b, intp, realp, None, Eout, flux, dlogflux, None, None)
        return None

    def retCodeAna(self, retval):
        """Function to parse INVLIB return codes"""
        #maybe replace this with defined exceptions
        err_codes = {1: 'Success', 0: 'Unknown Error (BUG)',
            -101: 'Null passed where pointer expected',
            -102: 'One or fewer valid data points',
            -103: 'One or more invalid (NaN or Inf) in input arrays, or unphysical flux/energy input',
            -104: 'No counts in any channel',
            -201: 'No function selected in function bitmap',
            -202: 'Invalid function selected in function bitmap',
            -301: 'Relativistic Maxwellian selected with invalid real_params',
            -302: 'Power law with exponential tail requested with NULL or -ve real_params',
            -401: 'Invalid minimizer selected',
            -402: 'Invalid iterations requested (N<=0)',
            -501: 'Invalid verbose flag or NULL outfile',
            -502: 'User-requested output file could not be opened'}
            
        return err_codes[retval]
    
    def anaSpecInv(self):
        """Analytic spectral inversion
        
        Calls ana_spec_inv from the INVLIB library.
        
        """
        try:
            assert self._params
        except AttributeError:
            warn = 'Warning: setup not completed'
            mssg = 'Attempting to call setParams() with defaults'
            raise AttributeError(warn+'\n'+mssg)
        retval = self._invlib.ana_spec_inv(*self._params)
        print('Analytic Spectral Inversion: %s\n' % self.retCodeAna(retval))
        
        if retval == 1:
            #if successful, save (flux, dlogflux)
            self.flux = list(self._params[-4])
            self.dlogflux = list(self._params[-3])
        
        return retval
    
    def _writeAnaSpec(self, retval, fnameout="specinv_test_py.out1"):
        """Method (unfinished) to write output E, flux, dlogflux"""
        #currently not called
        if retval == 1:
            Eout = list(self._params[-5])
            flux = list(self._params[-4])
            dlogflux = list(self._params[-3])
            
        outlist = []
        for eo, fl, dlf in zip(Eout, flux, dlogflux):
            #add time to 1st position in nested list
            row = ['%6f' % eo, '%6g' % fl, '%6g' % dlf]
            outlist.append(row)
        csv.writer(open(fnameout,'w'), delimiter=',').writerows(outlist)
        print('Output written to file: %s' % fnameout)
        
        return None
        
    def plot(self, **kwargs):
        """Method for quicklook plot of output (if extant in object)
        
        Accepts keyword args:
        title - string containing plot title
        """
        
        try:
            #needs matplotlib for plotting
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('Error: MatPlotLib import failed - please check install')
        
        try:
            ciu_flux = [(f*np.exp(1.96*d))-f for f,d in zip(self.flux,self.dlogflux)]
            cil_flux = [f-(f*np.exp(-1.96*d)) for f,d in zip(self.flux,self.dlogflux)]
        except:
            raise AttributeError('Error: Flux and dLogFlux must be successfully calculated')
            
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')
        ax.errorbar(self.Eout, self.flux, yerr=[ciu_flux, cil_flux])
        ax.set_ylabel('Flux [#/cm$^2$.s.sr.keV]')
        ax.set_xlabel('Energy [keV]')
        ax.set_ylim(ymin=min(self.flux))
        if 'title' in kwargs:
            ax.set_title(kwargs['title'])
        plt.show()
            
        return fig


class AngInv(InvBase):
    """interface to omni2uni and wide2uni functions from invlib"""
    def __init__(self, **kwargs):
        super(AngInv, self).__init__(self, **kwargs) 
        #set null params for later population
        self._int_params = [None]*5
        self._real_params = [None]*3
        self.omniflux = None
        self.domniflux = None
        self.wideflux = None
        self.dwideflux = None
        self.PAgrid = None
    
    def retCodes(self, retval):
        """Function to parse INVLIB return codes"""
        #maybe replace this with defined exceptions
        err_codes = {1: 'Success', 0: 'Unknown Error (BUG)',
            -101: 'Null passed where pointer expected',
            -102: 'One or fewer valid data points',
            -103: 'One or more invalid (NaN or Inf) in input arrays, or unphysical flux/energy input',
            -104: 'Negative or zero input [omniflux, wideflux] or their error estimates',
            -401: 'Invalid minimizer or root-finder selected',
            -402: 'Invalid iterations requested (N<=0)',
            -501: 'Invalid verbose flag or NULL outfile',
            -502: 'User-requested output file could not be opened',
            -601: 'Output pitch angle requested out of range'}
        
        return err_codes[retval]
            
    def testOmni2Uni(self):
        """Runs the test case given in omni2uni"""
        self._int_params[0] = 50 #; /* NA - number of angular gridpoints */
        self._int_params[2] = 1 #; /* 1 = verbose to standard out */
        self._int_params[3] = 3 #; /* minimizer, 0=BFGS, 3=NM */
        self._int_params[4] = 1000 #; /* maximumn # of iterations */
        self._real_params[0] = 300.0 #; /* 300 keV */
        self._real_params[1] = 40000.0/100.0 #; /* B/B0 */
        self._real_params[2] = 6.6 #; /* Lm */
        self.omniflux = 10000 #; /* typical value 1E+4 */
        self.domniflux = math.log(2)/2 #natural log
        
        #run case 1
        self.setParams()
        ret1 = self.omni2uni('TEM1')
        #run case 2
        self.setParams(method='VAM')
        ret2 = self.omni2uni()
        
        return ret1+ret2         
        
    def setParams(self, method='TEM1', alpha=5):
        """Method called to set parameters to pass to INVLIB"""
        if method.upper() not in ['TEM1', 'VAM']:
            raise ValueError('Invalid angular inversion method requested')       
        if method.upper() == 'TEM1':
            self._int_params[1]=-1
        elif method.upper() == 'VAM':
            self._int_params[1]=-2
        #test for presence of inputs
        if self.omniflux is None or self.domniflux is None:
            raise ValueError('Flux or Error on flux undefined')
        if None in self._int_params:
            raise ValueError('Integer parameters not defined')
        
        intp = ctypes.c_long*5
        self._int_params = [int(j) for j in self._int_params]
        intp = intp(*self._int_params)
        realp = ctypes.c_double*3
        realp = realp(*self._real_params)
        
        stub = ctypes.c_double
        if isinstance(self.omniflux, numbers.Real):
            Lof = 1
            self.omniflux = [self.omniflux]
        else:
            Lof = len(self.omniflux)
        if isinstance(self.domniflux, numbers.Real):
            Ldf = 1
            self.domniflux = [self.domniflux]
        else:
            Ldf = len(self.domniflux)
        omniflux = (stub*int(Lof))(*self.omniflux) #need to catch if unpopulated
        domniflux = (stub*int(Ldf))(*self.domniflux)
        
        #initialize outputs
        uflux = (Lof * ctypes.c_double)(0)
        dloguflux = (Ldf * ctypes.c_double)(0)
        
        self._paramsO = (omniflux, domniflux, intp, realp, None, 
            uflux, dloguflux)
        #self._paramsW = (wideflux, dwideflux, PAgrid, H, 
            #intp, realp, None, alpha, uflux, dloguflux)
        return None
    
    def omni2uni(self, method=None):
        '''Angular Inversion: Omnidirectional to unidirectional'''
        try:
            assert self._paramsO
        except:
            if method != None:
                self.setParams(method=method)
            else:
                self.setParams() 
        retval = self._invlib.omni2uni(*self._paramsO)
        print('omni2uni Angular Inversion: %s\n' % self.retCodes(retval))
        
        return retval

    def wide2uni(self, method=None):
        '''Angular Inversion: Wide angle to unidirectional'''
        try:
            assert self._paramsW
        except:
            if method != None:
                self.setParams(method=method)
            else:
                self.setParams()
        retval = self._invlib.wide2uni(*self._paramsW)
        print('wide2uni Angular Inversion: %s\n' % self.retCodes(retval))
        
        return retval
        

if __name__ == '__main__':
    #run test suite if called from command line
    exec(compile(open('invlib_test.py').read(), 'invlib_test.py', 'exec'))