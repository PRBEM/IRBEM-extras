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
Finalize read support for LANL standard response files
Add support for calling of ana_spec_inv_multi (within extant functions)
Add support for calling pc_spec_inv (and _multi)
Improve error trapping
Update unit tests for additional code
Documentation
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

##following two functions can be removed when 
##numpy1.5 compatibility in Python3 is tested
##currently only used in invlib_test.py
def ravel(mylist):
    newlist = []
    for row in mylist:
        # may have nested tuples or lists
        if type(row) in (list, tuple):
            newlist.extend(ravel(row))
        else:
            newlist.append(row)
    return newlist

def transpose(arr):
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
            
    def _readLANLResp(self, fname):
        """read 'standard' LANL format response function files"""
        #modified port of r_sdtfrmt_resp in papco
        #currently data is nested lists (row major)
        #need to extract the columns and tag with energy channel names
        #this could be put into numpy arrays (ensuring Py3k compliance)
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
        self.H = np.array(data, dtype=float).transpose().ravel().tolist()  ##need to select which energy channel... how to implement?
        self.Egrid = Earray 
        self.Hcols = col_hdr
        
        if self._verb:
            print('LANL response file read:\n%s' % header[:2])
        try:
            assert self.Eout#==self.Egrid
        except AttributeError:
            raise AttributeError('Energy grid for output not set')
        #except AssertionError:
            #print('Response function not on same energy grid as desired output: interpolating')
            ##automatically interpolate to a specified grid
            #tmpH = np.zeros((self.H.shape[0],len(self.Eout)))
            #for n in range(self.H.shape[0]):
                #tmpH[n,:] = np.interp(self.Eout, self.Egrid, self.H[n,:])
            #self.H = tmpH.ravel().tolist()

        return None
        
    def readRespFunc(self, fname=None, std='LANL'):
        #method stub for reading PRBEM and LANL format response functions
        if std.upper() not in ['LANL', 'PRBEM']:
            raise ValueError('Unknown response file format specified')
        if fname == None or not os.path.exists(fname):
            raise IOError('Requested file %s does not exist' % (fname))
        if std.upper() == 'LANL':
            rfun = self._readLANLResp(fname)
            #need to work out how to get this into format required by *_spec_inv
            #and wide2uni
            #currently hacked into readLANLResp - good spot??
        elif std.upper() == 'PRBEM':
            try:
                from spacepy import pycdf
                rfun = pycdf.CDF(fname)
                #now get appropriate numbers from CDF to pass to *_spec_inv
                #then set self.H
            except ImportError:
                raise ImportError('''SpacePy is not installed; CDF support is unavailable without SpacePy\n
                Please get SpacePy from http://spacepy.lanl.gov''')
        
        return None
        
    
class SpecInv(InvBase):
    """Spectral Inversion class using INVLIB
    
    Example use:
    
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
        self._default_Eout = [  0.1       ,   0.1274275 ,   0.16237767,   0.20691381,
         0.26366509,   0.33598183,   0.42813324,   0.54555948,
         0.6951928 ,   0.88586679,   1.12883789,   1.43844989,
         1.83298071,   2.33572147,   2.97635144,   3.79269019,
         4.83293024,   6.15848211,   7.8475997 ,  10.        ]

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
        datlist = [float(d.rstrip()) for d in infile if len(d)>2] #remove all trailing \n and store each non-empty line as a string in a list
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
        self._int_params[3] = func
        self._int_params[4] = minim
        self._int_params[5] = niter
        self._int_params[6] = self._verb
        self._int_params[7:] = [0, 0, 0]
        self.setParams()
        
        return None
        
    def setParams(self, fittype='ana', fnc=None):
        """Generic method for setting parameter inputs to either ana_spec* or pc_spec*"""
        try:
            assert self.counts
            assert len(self.dcounts)==len(self.counts)
            if not self._int_params[0]:
                self._int_params[0] = len(self.Hcols) #how to calc number of channels? NC
            if not self._int_params[1]:
                self._int_params[1] = len(self.Egrid)  #NE
            if not self._int_params[2]:
                self._int_params[2] = len(self.Eout)
        except (AssertionError, AttributeError):
            raise TypeError('Counts, Relative Error on Counts and Energy grid must be defined') 
        try:
            assert self.H
        except:
            raise AttributeError('Response functions must be specified')
        try:
            assert self.b
        except AttributeError:
            print('Background not specified, defaulting to zero')
            self.b = [0]*len(self.counts)
        
        #define outputs
        NC, NE = int(self._int_params[0]), int(self._int_params[1])
        NEout = int(self._int_params[2])
        c, dc = NC*ctypes.c_double, NC*ctypes.c_double
        Egrid = NE*ctypes.c_double
        H = NE*NC*ctypes.c_double
        Eout = NEout*ctypes.c_double
        b = NC*ctypes.c_double
        flux = (NEout * ctypes.c_double)(0)
        dlogflux = (NEout * ctypes.c_double)(0)
        try:
            Eout = Eout(*self.Eout)
            NEout = int(len(Eout))
        except AttributeError:
            NEout = int(NEout) 
            Eout = NEout*ctypes.c_double
            Eout = Eout(*self._default_Eout)  ##set default log grid of 20 energies
        
        flux = (NEout * ctypes.c_double)(0)
        dlogflux = (NEout * ctypes.c_double)(0)
        
        #check for absence of keywords and set defaults
        if fnc == None:
            self._int_params[3] = 1+2 #defaults to power law and exponential
        else:
            self._int_params[3] = fnc
        if self._int_params[4] == None: self._int_params[4] = 0 #default to BFGS minimizer
        if self._int_params[5] == None: self._int_params[5] = 10000
        if self._int_params[6] == None: self._int_params[6] = self._verb
        if self._int_params[7] == None: self._int_params[7] = 0
        if self._int_params[8] == None: self._int_params[8] = 0
        if self._int_params[9] == None: self._int_params[9] = 0
        
        intp = ctypes.c_long*10
        realp = ctypes.c_double*10
        if fittype == 'ana':
            intp = intp(*self._int_params)
            realp = realp(*[self.rme, 100, 345, 0, 0, 0, 0, 0, 0, 0])
        elif fittype == 'pc':
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
        '''Analytic spectral inversion'''
        assert self._params
        retval = self._invlib.ana_spec_inv(*self._params)
        print('Analytic Spectral Inversion: %s\n' % self.retCodeAna(retval))
        
        return retval
    
    def _writeAnaSpec(self, retval, fnameout="specinv_test_py.out1"):
        #currently not called
        if retval == 1:
            Eout = list(self._params[-5])
            flux = list(self._params[-4])
            dlogflux = list(self._params[-3])
            
        outlist=[]
        for eo, fl, dlf in zip(Eout, flux, dlogflux):
            #add time to 1st position in nested list
            row = ['%6f' % eo, '%6g' % fl, '%6g' % dlf]
            outlist.append(row)
        csv.writer(open(fnameout,'w'), delimiter=',').writerows(outlist)
        print('Output written to file: %s' % fnameout)
        
        return None


class AngInv(InvBase):
    #interface to omni2uni and wide2uni functions from invlib
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
        """Runs the test case given in omni2uni
        
        """
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
        if method.upper() not in ['TEM1', 'VAM']:
            raise ValueError('Invalid angular inversion method requested')       
        if method.upper()=='TEM1':
            self._int_params[1]=-1
        elif method.upper()=='VAM':
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
        if isinstance(self.omniflux,numbers.Real):
            Lof = 1
            self.omniflux = [self.omniflux]
        else:
            Lof = len(self.omniflux)
        if isinstance(self.domniflux,numbers.Real):
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
    
    #dum = pinv.SpecInv(verbose=1)
    #dum.counts = [4.15524e+02, 3.70161e+02, 2.42137e+02, 2.12097e+02, 1.47379e+02, 1.40524e+02, 9.37500e+01, 5.08064e+01, 1.93548e+01, 3.22581e+00]
    #dum.dcounts=[0.3466]*10
    #dum.Eout = [1.00000000e-02, 1.31795546e-02, 1.73700659e-02, 2.28929731e-02, 3.01719189e-02, 3.97652452e-02, 5.24088219e-02, 6.90724929e-02, 9.10344690e-02, 1.19979375e-01,   1.58127472e-01,   2.08404965e-01, 2.74668462e-01,   3.62000798e-01,   4.77100928e-01, 6.28797771e-01,   8.28727455e-01,   1.09222587e+00, 1.43950505e+00,   1.89720354e+00,   2.50042976e+00, 3.29545504e+00,   4.34326296e+00,   5.72422712e+00, 7.54427638e+00, 9.94302023e+00, 1.31044578e+01, 1.72710917e+01, 2.27625295e+01, 3.00000000e+01]
    #dum.readRespFunc(fname='../../projects/responses/sopa/mcp_output/version_1/electrons/sopa_1989-046_t2_HSP_elec.001')
    #dum.setParams()
    #dum.anaSpecInv()
    #import matplotlib.pyplot as plt
    #flux = list(dum._params[-4])
    #dlogflux = list(dum._params[-3])
    #plt.semilogy(dum.Eout[:30], dlogflux[:30], 'r-d')
    #plt.semilogy(dum.Eout[:30], flux[:30], 'b-d')
    #plt.show()