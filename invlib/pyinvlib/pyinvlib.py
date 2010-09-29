#!/usr/lib/env python
# -*- coding: utf-8 -*-

"""Python interface to the IRBEM InvLib library

Author: Steve Morley, Los Alamos National Laboratory (smorley@lanl.gov)
Date Created: 23 Sept. 2010

Running this wrapper from the command line calls a test to reproduce the specinv_test
where the setup and function call are from within Python

To Do:
Enable read support for PRBEM response files
Finalize read support for LANL standard response files
Generalize calling of ana_spec_inv
Add support for calling of ana_spec_inv_multi (within extant functions)
Add support for calling pc_spec_inv (and _multi)
Update unit tests for additional code
Documentation
Angular Inversion code
"""

import ctypes, ctypes.util, os, sys, csv, math, numbers

if sys.platform == 'linux2':
    ext = 'so'
elif sys.platform == 'darwin':
    ext = 'dylib'
else:
    raise NotImplementedError('Your platform is not currently supported')

floc = locals()['__file__'].split('/')
libpath = '/'.join(floc[:-1])+'/' #should use os funcs for this to maintain platform independence


##following two functions can be removed when numpy1.5 compatibility in Python3 is tested
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


class SpecInv(object):
    """Spectral Inversion class using INVLIB
    """
    def __init__(self, verbose=True, rme=0.511):
        try:
            self._invlib = ctypes.CDLL(libpath+'invlib.'+ext)
        except: #case for testing in IRBEM dist
            self._invlib = ctypes.CDLL('../invlib.'+ext)
        if verbose==True:
            self._verb = 1
        elif verbose==False:
            self._verb = 0
        else:
            self._verb = verbose
        
        #set some defaults
        self._int_params = [None]*10
        self._real_params = [None]*10 
        self.rme = rme #defaults to electron with all units in MeV
        self._default_Eout = [  0.1       ,   0.1274275 ,   0.16237767,   0.20691381,
         0.26366509,   0.33598183,   0.42813324,   0.54555948,
         0.6951928 ,   0.88586679,   1.12883789,   1.43844989,
         1.83298071,   2.33572147,   2.97635144,   3.79269019,
         4.83293024,   6.15848211,   7.8475997 ,  10.        ]
        
        return None
    
    def readLANLResp(self, fname):
        #read 'standard' LANL format response function files
        #modified port of r_sdtfrmt_resp in papco
        
        #currently data is nested lists (row major)
        #this should be put into numpy arrays (ensuring Py3k compliance)

        fobj = open(fname, 'r')
        header, data, Earray = [], [], []
        
        c_symb = ';'
        for line in fobj:
            if c_symb in line: #if prefixed by ';' comment symbol it's header
                header.append(line.rstrip())
            else:
                if line.rstrip()!='': #skip blank lines
                    dum = line.rstrip().split()
                    data.append(dum[1:])
                    Earray.extend(dum[0]) #get energies from first column
        
        #move first line to column header list
        col_hdr = data.pop(0)
        
        #set attributes
        #probably need to interpolate response function onto same energy grid as Egrid
        #which corresponds to the data
        self.H = data
        self.Hgrid = Earray
        self.Hcols = col_hdr

        if self._verb: print('LANL response file read:\n%s' % header[:2])
        
        try:
            assert self.Hgrid==self.Egrid
        except AttributeError:
            raise UserWarning('Energy grid (Egrid) for data not specified')
        except AssertionError:
            raise UserWarning('Response function not on same energy grid as data')

        return None

    def readTestInput(self, verbose=True, func=1+2, minim=0, niter=1000):
        #port of specinv_test.c
        self._verb = verbose
        try:
            fnamein = libpath+"specinv_test.in1"
            assert os.path.exists(fnamein)
        except AssertionError:
            try: #case for testing in IRBEM dist
                fnamein = "../specinv_test.in1"
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
        self._int_params[7:] = [0,0,0]

        self.setParams()
        
        return None
    
    def readRespFunc(self, fname=None, std='LANL'):
        #method stub for reading PRBEM and LANL format response functions
        if std.upper() not in ['LANL','PRBEM']:
            raise ValueError('Unknown response file format specified')
        if fname == None or not os.path.exists(fname):
            raise IOError('Requested file %s does not exist' % (fname))

        if std.upper()=='LANL':
            rfun = self.readLANLResp(fname)
            #need to work out how to get this into format required by *_spec_inv
            #currently hacked into readLANLResp - good spot??
        elif std.upper()=='PRBEM':
            try:
                from spacepy import pycdf
                rfun = pycdf.CDF(fname)
                #now get appropriate numbers from CDF to pass to *_spec_inv
                #then set self.H
            except ImportError:
                raise ImportError('''SpacePy is not installed; CDF support is unavailable without SpacePy\n
                Please get SpacePy from http://spacepy.lanl.gov''')
        
        return None
        
    def setParams(self, niter=10000, fittype='ana', NEout=20):
        """Generic method for setting parameter inputs to either ana_spec* or pc_spec*"""
        try:
            assert self.counts
            assert len(self.dcounts)==len(self.counts)
            #NC = 1 #how to calc number of channels?
            NE = int(len(self.Egrid))
        except (AssertionError, AttributeError):
            raise TypeError('Counts, Relative Error on Counts and Energy grid must be defined')
        
        try:
            assert self.H
        except:
            raise AttributeError('Response functions must be specified')
        
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
        if self._int_params[3]==None: self._int_params[3] = 1+2 #defaults to power law and exponential
        if self._int_params[4]==None: self._int_params[4] = 0
        if self._int_params[5]==None: self._int_params[5] = 10000
        
        intp = ctypes.c_long*10
        realp = ctypes.c_double*10
        if fittype=='ana':
            intp = intp(*self._int_params)
            #intp = intp(*[self.NC, NE, NEout, func, minim, niter, 
                #self._verb, 0, 0, 0])
            realp = realp(*[self.rme, 100, 345, 0, 0, 0, 0, 0, 0, 0])
        elif fittype=='pc':
            raise NotImplementedError
        
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
    
    def anaSpecInv(self, restmass=0.511):
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


class AngInv(object):
    #interface to omni2uni and wide2uni functions from invlib
    def __init__(self, verbose=True):
        try:
            self._invlib = ctypes.CDLL(libpath+'invlib.'+ext)
        except: #case for testing in IRBEM dist
            self._invlib = ctypes.CDLL('../invlib.'+ext)
        if verbose==True:
            self._verb = 1
        elif verbose==False:
            self._verb = 0
        else:
            self._verb = verbose
            
        self._int_params = [None]*5
        self._real_params = [None]*3
        
        #set null params for later population
        self.omniflux = None
        self.domniflux = None
        self.wideflux = None
        self.dwideflux = None
        self.PAgrid = None
        self.H = None
        
        return None
    
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
        #run the test case given in omni2uni     
        
        self._int_params[0] = 50#; /* NA - number of angular gridpoints */
        self._int_params[2] = 1#; /* 1 = verbose to standard out */
        self._int_params[3] = 3#; /* minimizer, 0=BFGS, 3=NM */
        self._int_params[4] = 1000#; /* maximumn # of iterations */

        self._real_params[0] = 300.0#; /* 300 keV */
        self._real_params[1] = 40000.0/100.0#; /* B/B0 */
        self._real_params[2] = 6.6#; /* Lm */

        self.omniflux = 10000#; /* typical value 1E+4 */
        self.domniflux = math.log(2)/2 #natural log
        
        #run case 1
        self.setParams()
        self.omni2uni('TEM1')
        #run case 2
        self.setParams(method='VAM')
        self.omni2uni()
        
        return None        
        
    def setParams(self, method='TEM1', alpha=5):
        if method.upper() not in ['TEM1', 'VAM']:
            raise ValueError('Invalid angular inversion method requested')
        
        if method.upper()=='TEM1':
            self._int_params[1]=-1
        elif method.upper()=='VAM':
            self._int_params[1]=-2
            
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
        return None

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
        return None
    

if __name__=='__main__':
    #run test suite if called from command line
    exec(compile(open('invlib_test.py').read(), 'invlib_test.py', 'exec'))