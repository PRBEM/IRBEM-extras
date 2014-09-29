//
//  FieldLineTable.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/2/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "FieldLineTable.h"
#include "Point.h"
#include "FieldModel.h"
#include "AffineTransform.h"
#include "FieldLine.h"
#include "ParticleTemplate.h"
#include "Thread.h"
#include <cassert>
#include <cstdio>
#include <stdexcept>

namespace UBK {
    using namespace std;

    //
    // Test
    //
#ifdef DEBUG
    void FieldLineTable::test()
    {
        {
            Point scale(.1, M_PI/180., 1.);
            Point offset;
            AffineTransform t(scale, offset);

            double IMFb1 = 0.;
            double compb2 = 4.;
            CompressedDipole fm(IMFb1, compb2);

            Point x1(2., 0.), x2(10., M_PI*2.);
            Point n1 = t.pointConvertedToNormal(x1);
            Point n2 = t.pointConvertedToNormal(x2);
            Key ub(n2.x, n2.y+1);
            Key lb(n1.x, n1.y-1);
            FieldLineTableCyl<FieldLineTable, FieldLine> tbl(lb, ub, &fm, &t);
            for (int i=60; i<=60; i++) {
                for (int j=20; j<=100; j++) {
                    Key key = Key(i, j, 0);
                    printf("%f\n", tbl[key]->magneticEquatorFieldMagnitude());
                }
                break;
            }
        }

        assert(0);
    }
#endif

    //
    // Implementation
    //

    class FieldLineProxy {
    public:
        FieldLine const* _fl;
        Mutex _key;
        FieldLineTable const* _table;

        ~FieldLineProxy() {
            delete _fl;
        };
        FieldLineProxy() : _fl(NULL), _key() {
        };

        FieldLine const* fieldLine(Key const key) {
            Locker l(_key);
            if (NULL == _fl) {
                _fl = _table->newFieldLineAtNode(key);
            }

            return _fl;
        };
    };

    // Reclaim memory
    FieldLineTable::~FieldLineTable()
    {
        // Free all map elements
        delete [] _table;
    }
    FieldLineTable::FieldLineTable (Key const lBound, Key const uBound) : _lBound(lBound), _uBound(uBound)
    {
        _N1 = (_uBound.i1-_lBound.i1+1);
        _N2 = (_uBound.i2-_lBound.i2+1);
        _N = _N1 * _N2;

        if (_N1<=0 || _N2<=0) {
            throw invalid_argument("Bounds are not incremental.");
        }

        _table = new FieldLineProxy[_N];
        for (long idx=0; idx<_N; idx++) {
            _table[idx]._table = this;
        }
    }

    FieldLine const* FieldLineTable::operator[](const Key key)
    {
        if (key.i1<_lBound.i1 || key.i1>_uBound.i1 ||
            key.i2<_lBound.i2 || key.i2>_uBound.i2) { // Out of bound
            return NULL;
        }

        long idx = (key.i1-_lBound.i1)*_N2 + (key.i2-_lBound.i2);
        return _table[idx].fieldLine(key);
    }

}