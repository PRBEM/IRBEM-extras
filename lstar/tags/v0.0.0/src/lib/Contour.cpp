//
//  Contour.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 3/30/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "Contour.h"
#include <cassert>
#include <cmath>
#include <set>
#include <cstdio>

// Tiny offset
#define EPSILON (1e-15)

namespace UBK {
    using namespace std;

    //////////////////////////////////////////////////////////////////////
    // Test block
    //////////////////////////////////////////////////////////////////////
#ifdef DEBUG
    void Contour::test()
    {
        class ContourTest : public ContourLevel, public ContourCallback {
        public:
            ~ContourTest () {};
            void getLevel_atNode (double *lvlOut, Key const xy) {
                double x = xy.x, y = xy.y;
                *lvlOut = sqrt(x*x + y*y);
            };
            void handleEventAtIndex_node_cellSide_coordinate_level_stop (long const idx, Key *keyInOut, CellSide const side, Point const coordinate, double const level, BOOL *stop) {
                printf("%f, %f; ",coordinate.x, coordinate.y);
                //printf("%d, %d; ",keyInOut->x, keyInOut->y);
                *stop = idx>100;
            };
        } t;
        Contour c(&t);

        Point pt(-2., -1., 0.);
        printf("c = [ ");
        c.enumerateCoordinatesStartingAtPoint_usingCallback(pt, &t);
        printf("];\n");

        assert(0);
    }
#endif

    //////////////////////////////////////////////////////////////////////
    // Class implementation
    //////////////////////////////////////////////////////////////////////
    Contour::Contour (ContourLevel *level) : _level(level)
    {
        assert(level);
    };

    ContourTermination Contour::enumerateCoordinatesStartingAtPoint_usingCallback(Point const pt, ContourCallback *callback) const
    {
        assert(callback);

        Key key((int)floor(pt.x), (int)floor(pt.y), 0); // Bottom-left
        ContourTermination closed = kContourOpen;

        //
        // Search for initial level by interpolating from the adjacent four nodes
        //
        double level0;
        {
            Key bottomLeft = key;
            Point dR(pt.x-bottomLeft.x, pt.y-bottomLeft.y, 0.);
            double blLevel, brLevel, tlLevel, trLevel;
            _level->getLevel_atNode(&blLevel, bottomLeft);
            _level->getLevel_atNode(&brLevel, Key(bottomLeft.x+1, bottomLeft.y, 0));
            _level->getLevel_atNode(&tlLevel, Key(bottomLeft.x, bottomLeft.y+1, 0));
            _level->getLevel_atNode(&trLevel, Key(bottomLeft.x+1, bottomLeft.y+1, 0));
            level0 = blLevel*(1.-dR.x)*(1.-dR.y) + brLevel*dR.x*(1.-dR.y) + tlLevel*(1.-dR.x)*dR.y + trLevel*dR.x*dR.y;

            if (isnan(level0)) {
                return kContourOpen;
            }
        }

        //
        // Select first side
        //
        BOOL stop = NO;
        long idx = 0;
        Point coordinate(0., 0., 0.);
        // At first try, the initial point at exact nodes (bottom-right, top-left and top-right corners) may be treated as being outside of the cell. Therefore it has to search for adjacent three more cells (right, top and top-right) if not found at first try.
        CellSide currentSide = this->firstSideOfLevel_cellOrigin_excludingSide_coordinate(level0, &key, kCellNone, &coordinate);
        if (currentSide == kCellNone) {
            key.x += 1;
            currentSide = kCellRightSide;
            currentSide = this->firstSideOfLevel_cellOrigin_excludingSide_coordinate(level0, &key, -currentSide, &coordinate);
            if (currentSide == kCellNone) {
                key.y += 1;
                currentSide = kCellTopSide;
                currentSide = this->firstSideOfLevel_cellOrigin_excludingSide_coordinate(level0, &key, -currentSide, &coordinate);
                if (currentSide == kCellNone) {
                    key.x -= 1;
                    currentSide = kCellLeftSide;
                    currentSide = this->firstSideOfLevel_cellOrigin_excludingSide_coordinate(level0, &key, -currentSide, &coordinate);
                    if (currentSide == kCellNone) {
                        return kContourOpen;
                    }
                }
            }
        }

        // Call user-supplied callback with first selected side
        callback->handleEventAtIndex_node_cellSide_coordinate_level_stop(idx++, &key, currentSide, coordinate, level0, &stop);

        //
        // Connect contour
        //
        {
            using namespace std;
            Key key0 = key;
            set<long> visitPool; // TODO: Use set with better performance
            while ( !(stop || visitPool.find(key.key)!=visitPool.end()) ) {
                visitPool.insert(key.key);
                if (kCellNone==(currentSide=this->firstSideOfLevel_cellOrigin_excludingSide_coordinate(level0, &key, -currentSide, &coordinate))) {
                    break;
                }
                callback->handleEventAtIndex_node_cellSide_coordinate_level_stop(idx++, &key, currentSide, coordinate, level0, &stop);// block(idx++, &key, currentSide, coordinate, level0, &stop, context);
            }
            closed = (key.key==key0.key);
        }

        return closed;
    }

    //
    // Private
    //
    UBK_INLINE Point AddXYToPoint (Point pt, double const x, double const y)
    {
        pt.x += x;
        pt.y += y;
        return pt;
    }

    CellSide Contour::firstSideOfLevel_cellOrigin_excludingSide_coordinate(const double lvl, Key *originInOut, const CellSide exSides, Point *coordOut) const
    {
        double lvl1, lvl2, dl1, dl2;
        Point origin((double)originInOut->x, (double)originInOut->y, 0.);
        Key originKey = *originInOut;

        //
        // Counterclockwise search from origin
        //
        if (kCellBottomSide!=exSides) {
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl1, &dl1, originKey, lvl); // origin
            originKey.x++;
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl2, &dl2, originKey, lvl); // lower right corner
            if (dl1*dl2 < 0.) {
                *coordOut = AddXYToPoint(origin, (dl1/(dl1-dl2)), 0.);
                originInOut->y -= 1; // Move to lower cell for next iteration
                return kCellBottomSide;
            }
        } else {
            originKey.x++; // Skip to right side
        }
        if (kCellRightSide!=exSides) {
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl1, &dl1, originKey, lvl); // lower right corner
            originKey.y++;
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl2, &dl2, originKey, lvl); // upper right corner
            if (dl1*dl2 < 0.) {
                *coordOut = AddXYToPoint(origin, 1., (dl1/(dl1-dl2)));
                originInOut->x += 1; // Move to right cell for next iteration
                return kCellRightSide;
            }
        } else {
            originKey.y++; // Skip to top side
        }
        if (kCellTopSide!=exSides) {
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl2, &dl2, originKey, lvl); // upper right corner
            originKey.x--;
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl1, &dl1, originKey, lvl); // upper left corner
            if (dl1*dl2 < 0.) {
                *coordOut = AddXYToPoint(origin, (dl1/(dl1-dl2)), 1.);
                originInOut->y += 1; // Move to upper cell for next iteration
                return kCellTopSide;
            }
        } else {
            originKey.x--; // Skip to left side
        }
        if (kCellLeftSide!=exSides) {
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl2, &dl2, originKey, lvl); // upper left corner
            this->getNodeLevel_deltaLevel_atNodePoint_contourLevel(&lvl1, &dl1, *originInOut, lvl); // origin
            if (dl1*dl2 < 0.) {
                *coordOut = AddXYToPoint(origin, 0., (dl1/(dl1-dl2)));
                originInOut->x -= 1; // Move to left cell for next iteration
                return kCellLeftSide;
            }
        }

        return kCellNone;
    }

    void Contour::getNodeLevel_deltaLevel_atNodePoint_contourLevel(double *nodeLevelOut, double *deltaOut, const Key nodePt, const double lvl) const
    {
        //
        // Check if the initial level is the same as the node level. If so, move node point slightly consistently.
        //
        _level->getLevel_atNode(nodeLevelOut, nodePt);
        if (*nodeLevelOut != lvl) {
            *deltaOut=*nodeLevelOut-lvl;
            return;
        }

        // x increment
        double nodeLevel2;
        _level->getLevel_atNode(&nodeLevel2, Key(nodePt.x+1, nodePt.y, 0));
        if (*nodeLevelOut != nodeLevel2) {
            *deltaOut=*nodeLevelOut-nodeLevel2;
            *deltaOut *= EPSILON;
            return;
        }
        // y increment
        _level->getLevel_atNode(&nodeLevel2, Key(nodePt.x, nodePt.y+1, 0));
        if (*nodeLevelOut != nodeLevel2) {
            *deltaOut=*nodeLevelOut-nodeLevel2;
            *deltaOut *= EPSILON;
            return;
        }
        
        *deltaOut = NAN;
    }

}