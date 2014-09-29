//
//  AffineTransform.h
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

#ifndef UBJDevelopment_AffineTransform_h
#define UBJDevelopment_AffineTransform_h

#include "Definitions.h"
#include "Point.h"

namespace UBK {

    //!
    //! Affine transform.
    //!
    class AffineTransform {
#ifdef DEBUG
    public:
        static void test ();
#endif

    protected:
        Point _scale;
        Point _offset;

    public:
        //! Scale factor (or spatial resolution, dx, dy or dz).
        Point const& scale () const {return _scale;};
        //! Offset in the FROM coordinate system.
        Point const& offset () const {return _offset;};

        virtual ~AffineTransform () {};
        AffineTransform (Point scale, Point offset);

        //!
        //! Forward transform.
        //! Normal = (Point-offset) / scale.
        //! @note In-place operation.
        //! @param [in] pt Points to be transformed.
        //! @param [out] normOut Points being tramsformed.
        //! @param count The number of elements.
        //!
        virtual void convertPoints_toNormal_count (Point const* pt, Point *normOut, long count) const;
        //!
        //! Forward transform.
        //! @sa convertPoints_toNormal_count
        //! @param [in] pt A point to be transformed.
        //! @return A point being transformed.
        //!
        virtual Point pointConvertedToNormal (Point pt) const;

        //!
        //! Backward transform.
        //! Point = Normal*scale + offset.
        //! @note In-place operation.
        //! @param [out] ptOut Points being transformed.
        //! @param [in] norm Points to be tramsformed.
        //! @param count The number of elements.
        //!
        virtual void convertPoints_fromNormal_count (Point *ptOut, Point const* norm, long count) const;
        //!
        //! Backward transform.
        //! @sa convertPoints_fromNormal_count
        //! @param [in] norm A point to be transformed.
        //! @return A point being transformed.
        //!
        virtual Point pointConvertedFromNormal (Point norm) const;
    };

}

#endif
