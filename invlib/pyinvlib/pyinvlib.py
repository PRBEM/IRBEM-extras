#!/usr/lib/env python
# -*- coding: utf-8 -*-

"""Python interface to the IRBEM InvLib library

Author: Steve Morley, Los Alamos National Laboratory (smorley@lanl.gov)
Date Created: 23 Sept. 2010

Running this wrapper from the command line calls a test to reproduce the specinv_test
where the setup and function call are from within Python

To Do:
Enable read support for PRBEM response files
Enable read support for LANL standard response files
Generalize calling of ana_spec_inv
Add support for calling of ana_spec_inv_multi (within extant functions)
Add support for calling pc_spec_inv (and _multi)
Update unit tests for additional code
Documentation
"""

import ctypes, ctypes.util, os, sys, csv, math

if sys.platform == 'linux2':
    ext = 'so'
elif sys.platform == 'darwin':
    ext = 'dylib'
else:
    raise NotImplementedError('Your platform is not currently supported')

floc = locals()['__file__'].split('/')
libpath = '/'.join(floc[:-1])+'/'

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
        header, data = [], []
        
        c_symb = ';'
        for line in fobj:
            if c_symb in line: #if prefixed by ';' comment symbol it's header
                header.append(line.rstrip())
            else:
                if line.rstrip()!='': #skip blank lines
                    data.append(line.rstrip().split())
        
        #move first line to column header list
        col_hdr = data.pop(0)

        return {'data': data, 'columns': col_hdr, 'header': header}


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
        
        #define vars
        c, dc = NC*ctypes.c_double, NC*ctypes.c_double
        Egrid = NE*ctypes.c_double
        H = NE*NC*ctypes.c_double
        Eout = NEout*ctypes.c_double
        b = NC*ctypes.c_double
        flux = (NEout * ctypes.c_double)(0)
        dlogflux = (NEout * ctypes.c_double)(0)

        #populate vars
        datlist = [float(d.rstrip()) for d in infile if len(d)>2] #remove all trailing \n and store each non-empty line as a string in a list
        c = c(*datlist[:NC])
        dc = dc(*datlist[NC:NC*2])
        Egrid = Egrid(*datlist[NC*2:(NC*2)+NE])
        dum = ((NC*2)+NE) #set counter for readability
        H = H(*datlist[dum:dum+(NE*NC)])
        dum = dum+(NE*NC)
        b = b(*datlist[dum:dum+NC])
        Eout = Eout(*datlist[dum+NC:])
        
        print('Python Analytic Spectral Inversion test: NC=%d, NE=%d, NEout=%d' % (NC, NE, NEout))
        
        intp = ctypes.c_long*10
        intp = intp(*[NC, NE, NEout, func, minim, niter, self._verb, 0, 0, 0])
        
        realp = ctypes.c_double*10
        realp = realp(*[0.511, 100, 345, 0, 0, 0, 0, 0, 0, 0])

        self._params = (c, dc, Egrid, H, b, intp, realp, None, Eout, flux, dlogflux, None, None)
        
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
        elif std.upper()=='PRBEM':
            try:
                from spacepy import pycdf
                rfun = pycdf.CDF(fname)
                #now get appropriate numbers from CDF to pass to *_spec_inv
            except ImportError:
                raise ImportError('''SpacePy is not installed; CDF support is unavailable without SpacePy\n
                Please get SpacePy from http://spacepy.lanl.gov''')
        
        return rfun
        
    def setParams(self, func=1+2, minim=0, niter=10000, fittype='ana', NEout=20):
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
        try:
            Eout = self.Eout
            NEout = int(len(Eout))
        except AttributeError:
            NEout = int(NEout)
            
            Eout = NEout*ctypes.c_double
            Eout = Eout(*self._default_Eout)  ##set default log grid of 20 energies
        
        flux = (NEout * ctypes.c_double)(0)
        dlogflux = (NEout * ctypes.c_double)(0)
        
        intp = ctypes.c_long*10
        realp = ctypes.c_double*10
        if fittype=='ana':
            intp = intp(*[self.NC, NE, NEout, func, minim, niter, 
                self._verb, 0, 0, 0])
            realp = realp(*[self.rme, 100, 345, 0, 0, 0, 0, 0, 0, 0])
        elif fittype=='pc':
            pass
        
        self._params = (self.counts, self.dcounts, self.Egrid, 
            self.H, self.b, intp, realp, None, Eout, flux, dlogflux, None, None)
        
        return None

    def retCodeAna(self, retval):
        #replace this with defined exceptions
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
            
        msg = err_codes[retval]
        
        return msg
    
    def anaSpecInv(self, restmass=0.511):
        '''Analytic spectral inversion'''
        try:
            params = self._params
        except AttributeError:
            self.setParams()

        retval = self._invlib.ana_spec_inv(*self._params)
        print('Analytic Spectral Inversion: %s\n' % self.retCodeAna(retval))\
        
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


class AngInv(object):
    def __init__(self):
        pass
    
    def omni2uni(self):
        pass
    
    def wide2uni(self):
        pass
    

if __name__=='__main__':
    execfile('invlib_test.py')
    