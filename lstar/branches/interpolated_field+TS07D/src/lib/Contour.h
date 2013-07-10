//
//  Contour.h
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

#ifndef UBJDevelopment_Contour_h
#define UBJDevelopment_Contour_h

#include "Definitions.h"
#include "Point.h"
#include "Key.h"

namespace UBK {

    //!
    //! Contour termination code.
    //!
    enum _ContourTermination {
        kContourOpen = 0, //!< Contour is disconnected.
        kContourClosed = 1 //!< Contour is connected.
    };
    //!
    //! Contour termination code type.
    //!
    typedef long ContourTermination;

    //!
    //! Cell side id.
    //!
    enum _CellSide {
        kCellLeftSide = -1, //!< Left side of cell
        kCellRightSide = 1, //!< Right side of cell
        kCellTopSide = 10, //!< Top side of cell
        kCellBottomSide = -10, //!< Bottom side of cell
        kCellNone = 0 //!< Unknown
    };
    //!
    //! Cell side id type.
    //!
    typedef long CellSide;

    //!
    //! Contour level callback abstract.
    //!
    class ContourLevel {
    public:
        virtual ~ContourLevel () {};
        //!
        //! Contour object requests for level (scalar value) at a node identified by a index.
        //! Contour class assumes regularly spaced grid space with any coordinate system. The caller should interprete this index into his/her coordinate system.
        //!
        virtual void getLevel_atNode (double *lvlOut, Key const xy) = 0;
    };

    //!
    //! Contour tracing event handler abstract.
    //! Is called at every iteration instance. It is solely caller's responsibility for what to do.
    //! @sa Particle
    //!
    class ContourCallback {
    public:
        virtual ~ContourCallback () {};
        //!
        //! Event handler.
        //! @param [in] idx Iteration count starting from 0.
        //! @param [in,out] keyInOut The cell index at which the tracing is about to be done as soon as it returns. A cell index is the node index at the lower-left corner. Note that this parameter is bidirectional. See Particle class for use.
        //! @param [in] side The edge that the contour touches. It is the side that just found.
        //! @param [in] coordinate Interpolated coordinate of the touching point. Caller should transform into his/her coordinate system.
        //! @param [in] level The contour level. Won't change during the tracing.
        //! @param [out] stop Set this to YES to terminate.
        //!
        virtual void handleEventAtIndex_node_cellSide_coordinate_level_stop (long const idx, Key *keyInOut, CellSide const side, Point const coordinate, double const level, BOOL *stop) = 0;
    };

    //!
    //! @brief Unidirectional contour tracer.
    //! @details Traces a single level contour starting at a point. The class merely wraps tracing algorithm and is dynamically configurable through level abstract class and contour callback abstract. The object is initialized with level abstract object, an object returning scalar level at a given coordinate. Once initialized, it enumerates every contour coordinate unless discontinuity is found.
    //! The algorithm assumes equi-spaced integer grid space (nodes are placed at integral location).
    //! The algorithm assumes infinite 2D space (of course it is limited by int type). The boundary should be set up through either level callback or user-supplied enumeration callback or both. In level callback, callers can assign NAN to level argument which will terminate tracing. In enumeration callback, set *stop to NO.
    //! @note In this version, the contour is traced in only one direction. Unless the contour is closed, opposite side of contour segment from the starting point will not be traced.
    //!
    class Contour {
#ifdef DEBUG
    public:
        static void test ();
#endif

    private:
        ContourLevel *_level;

    public:
        virtual ~Contour () {};
        //!
        //! Initialize with level object.
        //!
        Contour (ContourLevel *level);

        //!
        //! Initiates Contour tracing.
        //! Starting from a given point, iterates contour coordinates. At every coordinate which, the receiver calls user-supplied block. The caller should transform back to the original space.
        //! @param [in] pt Tracing starting point.
        //! @param [in] callback Callback object.
        //! @return ContourTermination parameter indicating the connectivity of the contour.
        //! @note The reciever does not impose any limitation on iteration count. The enumeration callback handler may want to set the maximum iteration count for safety.
        //!
        virtual ContourTermination enumerateCoordinatesStartingAtPoint_usingCallback (Point const pt, ContourCallback *callback) const;

    private:
        //
        // Helper functions. Internal use only.
        //
        CellSide firstSideOfLevel_cellOrigin_excludingSide_coordinate (double const lvl, Key *originInOut, CellSide const exSides, Point *coordOut) const;
        void getNodeLevel_deltaLevel_atNodePoint_contourLevel (double *nodeLevelOut, double *deltaOut, Key const nodePt, double const lvl) const;
    };

}

#endif
