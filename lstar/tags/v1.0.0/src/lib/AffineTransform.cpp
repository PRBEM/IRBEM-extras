//
//  AffineTransform.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/1/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "AffineTransform.h"
#include <cassert>
#include <iostream>

namespace UBK {

    //////////////////////////////////////////////////////////////////////
    // Test block
    //////////////////////////////////////////////////////////////////////
#ifdef DEBUG
    void AffineTransform::test()
    {
        Point scale(.1, .01, .001);
        Point offset(-1., -1., -1.);
        AffineTransform t(scale, offset);

        Point pt(1., 1., 1.);
        Point nt = t.pointConvertedToNormal(pt);
        std::cout << "nx=" << nt.x << ", ny=" << nt.y << ", nz=" << nt.z << std::endl;
        pt = t.pointConvertedFromNormal(nt);
        std::cout << "px=" << pt.x << ", py=" << pt.y << ", pz=" << pt.z << std::endl;

        assert(0);
    }
#endif

    //
    // Constructor
    //
    AffineTransform::AffineTransform (Point scale, Point offset) : _scale(scale), _offset(offset)
    {
        assert(_scale.x>0. && _scale.y>0. && _scale.z>0.);
    }

    /*! Forward transform. Normal = (Point-offset) / scale.
     */
    void AffineTransform::convertPoints_toNormal_count (Point const* pt, Point *normOut, long count) const
    {
        for (long idx=0; idx<count; idx++) {
            normOut[idx].x = (pt[idx].x - _offset.x) / _scale.x;
            normOut[idx].y = (pt[idx].y - _offset.y) / _scale.y;
            normOut[idx].z = (pt[idx].z - _offset.z) / _scale.z;
        }
    }
    Point AffineTransform::pointConvertedToNormal (Point pt) const
    {
        Point norm;
        this->convertPoints_toNormal_count(&pt, &norm, 1);
        return norm;
    }

    /*! Backward transform. Point = Normal*scale + offset.
     */
    void AffineTransform::convertPoints_fromNormal_count (Point *ptOut, Point const* norm, long count) const
    {
        for (long idx=0; idx<count; idx++) {
            ptOut[idx].x = norm[idx].x*_scale.x + _offset.x;
            ptOut[idx].y = norm[idx].y*_scale.y + _offset.y;
            ptOut[idx].z = norm[idx].z*_scale.z + _offset.z;
        }
    }
    Point AffineTransform::pointConvertedFromNormal (Point norm) const
    {
        Point pt;
        this->convertPoints_fromNormal_count(&pt, &norm, 1);
        return pt;
    }

}