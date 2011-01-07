#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Testing suite for the Python interface to the IRBEM InvLib library

Author: Steve Morley, Los Alamos National Laboratory (smorley@lanl.gov)
Date Created: 23 Sept. 2010
"""
try:
    import unittest_pretty as utp
except:
    pass
import unittest
import math
import pyinvlib as pinv
from pyinvlib import ravel, transpose

def feq(x,y, precision=0.0000005):
    """
    compare two floating point values if they are equal

    param x: a number
    type x: float
    param y: float  or array of floats
    type y: float
    keyword precision: precision for equal (default 0.0000005)
    type precision: float
    return: True (equal) or False (not equal)
    rtype: boolean

    author: Josef Koller
    organization: Los Alamos National Lab
    contact: jkoller@lanl.gov
    
    version: V1: 20-Jan-2010
    version: V2: 18-May-2010: User-specified precision added
       
    >>> index = where( feq(Lpos,Lgrid) ) # use float point comparison
 
    """

    boolean = abs(x-y) <= (abs(x+y)*precision)

    return boolean
    

class SpecInvTests(unittest.TestCase):

    def setUp(self):
        super(SpecInvTests, self).setUp()        

    def tearDown(self):
        super(SpecInvTests, self).tearDown()

    def testSpecInvTest(self):
        """running ana_spec_inv with test input file should match C output
        
        Basic regression test for ana_spec_inv
        """
        
        dum = pinv.SpecInv(rme=511)
        dum.readTestInput(verbose=False)
        result = dum.anaSpecInv()
        
        eout = list(dum._params[-5])
        flux = list(dum._params[-4])
        dlogflux = list(dum._params[-3])
        pyoutput = transpose([eout, flux, dlogflux])
        
        coutput = [ [0.0100,4.1954e+05,8.4936e-01],
                    [0.6159,2.0872e+05,3.7389e-01],
                    [1.2217,1.0599e+05,3.1122e-01],
                    [1.8276,5.3918e+04,2.6611e-01],
                    [2.4334,2.7446e+04,2.2762e-01],
                    [3.0393,1.3976e+04,1.9519e-01],
                    [3.6452,7.1190e+03,1.7126e-01],
                    [4.2510,3.6267e+03,1.5971e-01],
                    [4.8569,1.8478e+03,1.6350e-01],
                    [5.4627,9.4155e+02,1.8200e-01],
                    [6.0686,4.7980e+02,2.1166e-01],
                    [6.6744,2.4452e+02,2.4874e-01],
                    [7.2803,1.2462e+02,2.9056e-01],
                    [7.8862,6.3513e+01,3.3548e-01],
                    [8.4920,3.2371e+01,3.8251e-01],
                    [9.0979,1.6500e+01,4.3103e-01],
                    [9.7037,8.4100e+00,4.8066e-01],
                    [10.3096,4.2868e+00,5.3112e-01],
                    [10.9155,2.1851e+00,5.8224e-01],
                    [11.5213,1.1139e+00,6.3390e-01],
                    [12.1272,5.6779e-01,6.8599e-01],
                    [12.7330,2.8944e-01,7.3846e-01],
                    [13.3389,1.4755e-01,7.9124e-01],
                    [13.9447,7.5215e-02,8.4429e-01],
                    [14.5506,3.8343e-02,8.9757e-01],
                    [15.1565,1.9547e-02,9.5107e-01],
                    [15.7623,9.9646e-03,1.0048e+00],
                    [16.3682,5.0799e-03,1.0586e+00],
                    [16.9740,2.5897e-03,1.1126e+00],
                    [17.5799,1.3202e-03,1.1667e+00],
                    [18.1858,6.7307e-04,1.2210e+00],
                    [18.7916,3.4314e-04,1.2754e+00],
                    [19.3975,1.7493e-04,1.3299e+00],
                    [20.0033,8.9184e-05,1.3844e+00],
                    [20.6092,4.5468e-05,1.4391e+00],
                    [21.2151,2.3180e-05,1.4939e+00],
                    [21.8209,1.1818e-05,1.5487e+00],
                    [22.4268,6.0250e-06,1.6036e+00],
                    [23.0326,3.0717e-06,1.6586e+00],
                    [23.6385,1.5661e-06,1.7136e+00],
                    [24.2443,7.9842e-07,1.7687e+00],
                    [24.8502,4.0706e-07,1.8239e+00],
                    [25.4561,2.0754e-07,1.8791e+00],
                    [26.0619,1.0581e-07,1.9344e+00],
                    [26.6678,5.3946e-08,1.9897e+00],
                    [27.2736,2.7504e-08,2.0451e+00],
                    [27.8795,1.4022e-08,2.1005e+00],
                    [28.4854,7.1493e-09,2.1559e+00],
                    [29.0912,3.6450e-09,2.2114e+00],
                    [29.6971,1.8584e-09,2.2669e+00]]

        self.assertEqual(result, 1)
        
        for pel, cel in zip(ravel(pyoutput), ravel(coutput)):
            self.assertTrue(feq(pel, cel, precision=0.00005))
            
    def testSpecInvErrors(self):
        """running ana_spec_inv with bad test inputs should fail"""
        
        dum = pinv.SpecInv()
        dum.readTestInput(verbose=False, niter=-1)
        result = dum.anaSpecInv()
        self.assertEqual(result, -402)
        
        dum = pinv.SpecInv()
        dum.readTestInput(verbose=False, func=0)
        result = dum.anaSpecInv()
        self.assertEqual(result, -201)
        
        dum = pinv.SpecInv()
        dum.readTestInput(verbose=False, minim=9)
        result = dum.anaSpecInv()
        self.assertEqual(result, -401)
    
    def testSpecInvNoInput(self):
        """should raise errors if no inputs given"""
        
        dum = pinv.SpecInv()
        foo = lambda:dum.anaSpecInv()
        self.assertRaises(AttributeError, foo)


class AngInvTests(unittest.TestCase):

    def setUp(self):
        super(AngInvTests, self).setUp()        

    def tearDown(self):
        super(AngInvTests, self).tearDown()
        
    def testAngInvMethod(self):
        """should raise errors if bad method given"""
        
        dum = pinv.AngInv()
        func = lambda:dum.omni2uni(method='bad')
        self.assertRaises(ValueError, func)
        func = lambda:dum.wide2uni(method='bad')
        self.assertRaises(ValueError, func)
        
    def testAngInvOmni_regress(self):
        '''Basic regression test for omni2uni'''  
        dum = pinv.AngInv()
        dum._int_params[0] = 50 #; NA - number of angular gridpoints
        dum._int_params[2] = 1 #; 1 = verbose to standard out
        dum._int_params[3] = 3 #; minimizer, 0=BFGS, 3=NM
        dum._int_params[4] = 1000 #; maximumn # of iterations
        dum._real_params[0] = 300.0 #; 300 keV
        dum._real_params[1] = 40000.0/100.0 #; B/B0
        dum._real_params[2] = 6.6 #; Lm
        dum.omniflux = 10000 #; typical value 1E+4
        dum.domniflux = math.log(2)/2 # natural log
        
        #run case 1
        dum.setParams(alpha=5)
        ret1 = dum.omni2uni('TEM1')
        self.assertEqual(ret1, 1) #Test for convergence in omni2uni
        self.assertAlmostEqual(dum.uniflux[0], 28289.428111933394, places=6)
        #run case 2
        dum.setParams(alpha=5, method='VAM')
        ret2 = dum.omni2uni()
        self.assertEqual(ret2, 1) #Test for convergence in omni2uni
        self.assertAlmostEqual(dum.uniflux[0], 29295.302704246446, places=6)

        
if __name__ == "__main__":
    try:
        utp.main()
    except:
        unittest.main()
