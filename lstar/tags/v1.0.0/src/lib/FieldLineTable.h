//
//  FieldLineTable.h
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

#ifndef UBJDevelopment_FieldLineTable_h
#define UBJDevelopment_FieldLineTable_h

#include "Definitions.h"
#include "Key.h"
#include "AffineTransform.h"

namespace UBK {

    //
    // Forward declarations.
    //
    class FieldLine;
    class FieldLineProxy;

    //!
    //! Fixed size 2D Field line hash table abstract.
    //! Implements base functionality of field line hash table.
    //!
    //! Use template form instead.
    //!
    class FieldLineTable {
#ifdef DEBUG
    public:
        static void test();
#endif

    private:
        Key _uBound;
        Key _lBound;
        long _N1, _N2, _N;

    private:
        //
        // Copy operation is not allowed.
        //
        FieldLineTable (FieldLineTable const& flt);
        FieldLineTable &operator = (FieldLineTable const& flt);

    private:
        //
        // Uses proxy to reserve space in the table.
        //
        FieldLineProxy *_table;

    public:
        virtual ~FieldLineTable ();
        //!
        //! Initialize with table size.
        //! The lower and upper bound of the index space determine the table size which is used to pre-allocate a memory block.
        //!
        //! Exception will be thrown if the bounds are not incremental.
        //! @note Only i1 and i2 components are used (2D table).
        //! @param [in] lBound Lower bound index.
        //! @param [out] uBound Upper bound index.
        //!
        FieldLineTable (Key const lBound, Key const uBound);

        //!
        //! Accessor to field line object associated with key (or index).
        //! @param [in] key A node index.
        //! @return NULL when the key is out of table range. Pointer to an object to a FieldLine base class.
        //!
        virtual FieldLine const* operator [] (Key const key);

        //!
        //! Create (with new operator) and returns a FieldLine (or its subclass) object.
        //! @note Table takes the owership of the created object. Caller should not delete the created object.
        //! @param [in] key A node index at which the field line information is asked.
        //! @return Newly created a FieldLine object.
        //!
        virtual FieldLine const* newFieldLineAtNode (Key const key) const = 0;
    };

}

#endif
