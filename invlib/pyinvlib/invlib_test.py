#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Testing suite for the Python interface to the IRBEM InvLib library

Author: Steve Morley, Los Alamos National Laboratory (smorley@lanl.gov)
Date Created: 23 Sept. 2010
"""

import unittest
import pyinvlib as pinv
from pyinvlib import ravel, transpose

class SpecInvTests(unittest.TestCase):

    def setUp(self):
        super(SpecInvTests, self).setUp()        

    def tearDown(self):
        super(SpecInvTests, self).tearDown()

    def testSpecInvTest(self):
        """running ana_spec_inv with test input file should match C output"""
        
        dum = pinv.SpecInv()
        dum.readTestInput(verbose=False)
        result = dum.anaSpecInv()
        
        Eout = list(dum._params[-5])
        flux = list(dum._params[-4])
        dlogflux = list(dum._params[-3])
        pyoutput = transpose([Eout,flux,dlogflux])
        
        coutput = [[0.01,419526,0.849366],
            [0.615859,208717,0.373893],
            [1.22172,105992,0.311219],
            [1.82758,53916,0.266113],
            [2.43343,27445.2,0.227624],
            [3.03929,13975.9,0.19519],
            [3.64515,7118.68,0.171259],
            [4.25101,3626.53,0.159711],
            [4.85687,1847.72,0.163502],
            [5.46273,941.506,0.182004],
            [6.06859,479.779,0.211666],
            [6.67444,244.504,0.248739],
            [7.2803,124.61,0.290561],
            [7.88616,63.5089,0.335486],
            [8.49202,32.3693,0.382516],
            [9.09788,16.4986,0.43104],
            [9.70374,8.4095,0.480663],
            [10.3096,4.28651,0.531127],
            [10.9155,2.18498,0.582252],
            [11.5213,1.11378,0.633909],
            [12.1272,0.567748,0.686005],
            [12.733,0.289415,0.738471],
            [13.3389,0.147534,0.79125],
            [13.9447,0.0752085,0.844301],
            [14.5506,0.0383396,0.897589],
            [15.1565,0.0195449,0.951086],
            [15.7623,0.00996374,1.00477],
            [16.3682,0.00507944,1.05862],
            [16.974,0.00258948,1.11262],
            [17.5799,0.00132012,1.16676],
            [18.1858,0.000673001,1.22102],
            [18.7916,0.000343101,1.2754],
            [19.3975,0.000174916,1.32989],
            [20.0033,8.91747e-05,1.38447],
            [20.6092,4.54627e-05,1.43914],
            [21.2151,2.31778e-05,1.4939],
            [21.8209,1.18165e-05,1.54873],
            [22.4268,6.02434e-06,1.60364],
            [23.0326,3.07137e-06,1.65862],
            [23.6385,1.56587e-06,1.71366],
            [24.2443,7.98327e-07,1.76877],
            [24.8502,4.07013e-07,1.82393],
            [25.4561,2.07509e-07,1.87915],
            [26.0619,1.05795e-07,1.93442],
            [26.6678,5.39384e-08,1.98973],
            [27.2736,2.74999e-08,2.0451],
            [27.8795,1.40205e-08,2.10051],
            [28.4854,7.14826e-09,2.15596],
            [29.0912,3.64449e-09,2.21145],
            [29.6971,1.85812e-09,2.26698]]
            
        self.assertEqual(result,1)
        
        for pel, cel in zip(ravel(pyoutput),ravel(coutput)):
            pel = '%6g' % pel
            self.assertAlmostEqual(float(pel), cel, places=6)
            
    def testSpecInvErrors(self):
        """running ana_spec_inv with bad test inputs should fail"""
        
        dum = pinv.SpecInv()
        dum.readTestInput(verbose=False, niter=-1)
        result = dum.anaSpecInv()
        self.assertEqual(result, -402)
        
        dum.readTestInput(verbose=False, func=0)
        result = dum.anaSpecInv()
        self.assertEqual(result, -201)
        
        dum.readTestInput(verbose=False, minim=9)
        result = dum.anaSpecInv()
        self.assertEqual(result, -401)
    
    def testSpecInvNoInput(self):
        """should raise errors if no inputs given"""
        
        dum = pinv.SpecInv()
        foo = lambda:dum.anaSpecInv()
        self.assertRaises(TypeError, foo)


class AngInvTests(unittest.TestCase):

    def setUp(self):
        super(AngInvTests, self).setUp()        

    def tearDown(self):
        super(AngInvTests, self).tearDown()
        
    def testAngInvMethod(self):
        """should raise errors if bad method given"""
        
        dum = pinv.AngInv()
        foo = lambda:dum.omni2uni(method='bad')
        self.assertRaises(ValueError, foo)
        foo = lambda:dum.wide2uni(method='bad')
        self.assertRaises(ValueError, foo)

if __name__ == "__main__":
    unittest.main()
