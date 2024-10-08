README for var2hdf5

See LICENSE.txt

var2hdf5 is a package that provides ways to save structured/array data to 
HDF5 files and load them again. The intent is to facilitate data transport
between computer languages.

Each language will implement a var2hdf5 writer function and an 
hdf52var reader function

A global HDF5 file attribute VAR2HDF5_SPEC_VERSION (an integer, starting at 0)
will be used to track changes to the file spec and allow reverse compatibility.

A global HDF5 file attribute 'writer' (a string) will be used to record what
program wrote the HDF5 file. E.g., 'python/var2hdf5.py'

Strings, including attributes, may be encoded differently by different languages.
The reader will need to account for this when processing type and subtype attributes
and when processing variables of type str and date.

The variables are stored recursively under /var
Each variable has a type attribute (see table below)

type attr | python      | idl         | matlab      | (D)ataset or (G)roup, notes
int       | int         | LONG        | int64       | D signed long integer
float     | float       | DOUBLE      | double      | D real number
bool      | bool        | BOOLEAN     | logical     | D boolean True or False (HDF5 enum)
str       | str         | STRING      | char        | D string as UTF-8/ASCII encoded bytes
date      | datetime    | STRING      | char        | D date/time stored as ISO 8601
timedelta | timedelta   | DOUBLE      | double      | D time difference in seconds
null      | None        | list()      | [] or zeros | D indicates absent/empty
list      | list        | LIST        | cell        | D list of items of heterogeneous types
array     | np.ndarray  | (array)     | N-D matrix  | D (N-D) list of variables of the same type
dict      | dict        | HASH        | struct      | G key-value pairs
table     | np.ndarray  | ORDEREDHASH | struct      | G collection of named variables of same length

Additional notes

list: 
    list are stored item-by-item under /__list_#__/
    The value stored at the nominal location (e.g., /var/mylist) of a list 
      variable is the list number (#) in /__list_#__/.
    Attribute list_length specifies number of entries in list
    List item names are zero-based, zero-padded integers as strings with enough 
      depth to hold all entries. So a list with 10 entries will have items "0".."9"
      and a list with 100 entries will have items "00".."99"
    List items can be any type, including other lists
    
array:
    Attribute subtype specifies the data type for the array entries
    Array of size 0 is written with subtype null
    
dict:
    key value pairs
    Note that IDL, being case-insensitive, will write all structure tag names as upper
      case. Conversely, a hash is used when reading dicts to preserve case-sensitivity.
      IDL's !NULL type cannot be used because it is ignored by IDL in lists and structures.
      E.g., [!NULL,3,2,!NULL,1] looks to IDL like [3,2,1], and {a:!NULL} is not allowed.

table:
    For tabular data that has named columns
    Stored like a dict, but each variable's entry has a 'column' attribute
      which is a zero-based column index used to reconstruct key order
    Attribute nested specifies if any of the columns contain non-scalar types,
      Some quantities in these structure arrays are themselves arrays


