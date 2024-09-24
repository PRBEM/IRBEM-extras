"""
var2hdf5.py - write variables to and read variables from hdf5 files
Prinicpal author: Paul O'Brien

public functions:
    var2hdf5 - save a variable to an HDF5 file
    hdf52var - read a variable from an HDF5 file (for files created by var2hdf5)
    test - run several test cases that write, read, and compare
"""

import os
import h5py
import numpy as np
import datetime as dt

VAR2HDF5_SPEC_VERSION = 0 # integer

def var2hdf5(var,filename,converter=None):
    """
    var2hdf5(var,filename,converter=None)
    write variable to hdf5 file
    var - variable to write
    filename - HDF5 filename to write
    converter - function pointer that converts variables to HDF5-recognized types
    calls var = conveter(var) before trying to handle the variable
    see README for conversion spec
    """
    if converter and not callable(converter):
        raise Exception('converter must be a callable')
    
    lists = 0 # number of lists written so far
    
    def to_friendly(var):
        """
        (var,attrs) = to_friendly(var)
        convert var to an HDF5-friendly type
        attrs holds attributes
        var is None if fail
        """
        new_var = None
        attrs = {}
        # type cast down to supported basic types
        for base in [bool,int,float,str]:
            if isinstance(var,base):
                new_var = base(var)
                attrs['type'] = type(new_var).__name__
                break
            
        # convert datetime to ISO 8601 string
        if isinstance(var,(dt.datetime,dt.date)):
            new_var = var.isoformat()
            attrs['type'] = 'date'

        # convert timedelta to total seconds
        if isinstance(var,dt.timedelta):
            new_var = var.total_seconds()
            attrs['type'] = 'timedelta'
            
        return (new_var,attrs)

    def write_array(fp,group,var):
        """
        write_array(fp,group,var)
        write numpy array either as 'array' or as 'table'
        """
        if len(var.dtype)==0: # simple array
            for base in [int,float,bool]:
                if np.issubdtype(var.dtype,base):
                    fp[group] = var
                    fp[group].attrs['type'] = 'array'
                    fp[group].attrs['subtype'] = base.__name__
                    return
            if var.size == 0: # for any other empty types, treat as empty array of None
                fp[group] = np.array([])
                fp[group].attrs['type'] = 'array'
                fp[group].attrs['subtype'] = 'null'
                return
            
            sample = var.ravel()[0]
            
            if sample is None: # none has zero int values and subtype null
                fp[group] = np.zeros(var.shape,dtype=int)
                fp[group].attrs['type'] = 'array'
                fp[group].attrs['subtype'] = 'null'
                return

            if isinstance(sample,str): # strings are a little picky, have to convert to dtype=S first
                fp[group] = np.array(var,dtype='S')
                fp[group].attrs['type'] = 'array'
                fp[group].attrs['subtype'] = 'str'
                return

            if isinstance(sample,(dt.datetime,dt.date))                    :
                new_var = np.array([v.isoformat() for v in var.ravel()],dtype='S')
                new_var.shape = var.shape
                fp[group] = np.array(new_var,dtype='S')
                fp[group].attrs['type'] = 'array'
                fp[group].attrs['subtype'] = 'date'
                return

            if isinstance(sample,dt.timedelta):
                new_var = np.array([v.total_seconds() for v in var.ravel()])
                new_var.shape = var.shape
                fp[group] = new_var
                fp[group].attrs['type'] = 'array'
                fp[group].attrs['subtype'] = 'timedelta'
                return
            raise Exception('Unable to process numpy array of type %s' % var.dtype)
        # now handle tables
        # write as dict, but with extra variable attr 'column' indicating column number
        # group attrs are type=table and nested=nested?
        # nested is when a column has more than 1 dimension
        fp.create_group(group)
        fp[group].attrs['type'] = 'table'
        nested = False
        for (column,key) in enumerate(var.dtype.names):
            vgroup = group + '/' + key
            write_recursive(fp,vgroup,var[key])
            fp[vgroup].attrs['column'] = column
            nested = nested or (var[key].ndim>1)
        fp[group].attrs['nested'] = nested

    def write_list(fp,group,var):
        """
        write_list(fp,group,var)
        entries of list to /__list_#__/
        where # indicates which list w/in the file (starting at 0)
        """
        nonlocal lists
        fp[group] = lists # number of lists so far
        fp[group].attrs['type'] = 'list'
        fp[group].attrs['list_length'] = len(var)
        ref = '/__list_%d__' % lists
        fp.create_group(ref)
        lists = lists+1 # increment list
        if len(var) == 0:
            return # done, no list entries to write
        fmt = '%s/%%0%.0fd' % (ref,np.ceil(np.log10(len(var))))
        for (i,v) in enumerate(var):
            write_recursive(fp,fmt % i,v)

    def write_recursive(fp,group,var):
        """
        write_recursive(fp,group,var)
        write variable recursively to group (str)
        """
        nonlocal converter
        original_type = type(var)
        if converter is not None:
            var = converter(var)
            
        if isinstance(var,np.ndarray):
            write_array(fp,group,var)
            return

        # call built in converters tolist, todict if present
        for conv in ['tolist','to_list','todict','as_dict']:
            if hasattr(var,conv):
                conv_func = getattr(var,conv)
                var = conv_func()
        
        if var is None:
            fp[group] = 0
            fp[group].attrs['type'] = 'null'
            return

        (new_var,attrs) = to_friendly(var)
        if new_var: # write regular scalars
            fp[group] = new_var
            for a in attrs:
                fp[group].attrs[a] = attrs[a]
            return # done

        # now look at compound types

        # tuples and lists become new lists
        if isinstance(var,(tuple,list)):
            old = var
            var = []
            for v in old:
                var.append(v)
        
        # write dicts recursively
        if isinstance(var,dict):
            fp.create_group(group)
            fp[group].attrs['type'] = 'dict'
            for key in var:
                write_recursive(fp,group+'/'+key,var[key])
            return
        
        # write lists into '/__list_#__/'
        if isinstance(var,list):
            write_list(fp,group,var)
            return

        raise Exception('Unable to convert data of type %s to HDF5-compatible form' % original_type)

    if os.path.exists(filename):
        os.remove(filename)
    if os.path.exists(filename):
        raise Exception('Could not remove %s (required to write a new hdf5 file)' % filename)
    with h5py.File(filename,'w') as fp:
        fp.attrs['VAR2HDF5_SPEC_VERSION'] = VAR2HDF5_SPEC_VERSION
        fp.attrs['writer'] = 'python/' + os.path.basename(__file__)
        write_recursive(fp,'/var',var)

def hdf52var(filename,converter=None):
    """
    var = hdf52var(filename,converter=None)
    read variable from hdf5 file
    filename - HDF5 filename to read
    converter - function pointer that converts variables to HDF5-recognized types
    calls var = conveter(var,attrs) before trying to handle the variable
    var - variable that was read
    see README for conversion spec
    """
    def read_list(fp,list_num,size):
        """
        var = read_list(fp,list_num,size)
        fp = HDF5 file pointer
        list_num # in /__list_#__/
        size - number of entries in list
        var - the variable read from the file
        """
        
        var = []
        if size==0:
            return var
        ref = '/__list_%d__' % list_num # group for this list
        # build 0-padded format for list entries:
        fmt = '%s/%%0%.0fd' % (ref,np.ceil(np.log10(size)))
        for i in range(size):
            var.append(read_recursive(fp,fmt % i))
        return var

    def read_recursive(fp,group):
        """
        var = read_recursive(fp,parent,,key)
        retrieve value from hdf5 file pointer
        """
        
        def decode_str(var):
            """
            var = decode_str(var)
                decode various string types
                returns a str
            """
            if isinstance(var,str):
                return var
            return var.decode('utf-8')
        def decode_date(var):
            """
            decodes date in YYYY-MM-DD:HH:MM:SS...Z format
            a bit smarter than strptime: allows more precision on
            seconds, allows with or without 'Z' suffix
            """
            var = decode_str(var)
            parts = var.replace('Z','').replace('-',':').replace('T',':').split(':')
            ints = [int(v) for v in parts[:5]]
            parts[5] = float(parts[5])
            parts[5] = round(parts[5]*1e6)/1e6 # round to nearest microsecond
            s = int(np.floor(parts[5])) # seconds
            us = int(round((parts[5]-s)*1e6)) # microseconds
            date = dt.datetime(*ints,s,us)
            return date
        
        nonlocal converter
            
        ty = decode_str(fp[group].attrs['type'])
        if isinstance(fp[group],h5py.Dataset):
            # not a group, read & process dataset
            if hasattr(fp[group],'value'):
                var = fp[group].value
            else:
                var = fp[group][()]
            if converter:
                var = converter(var,fp[group]['attrs'])
                
            if ty in ['int','float','bool']:
                return var
            if ty == 'str':
                return decode_str(var)
            if ty == 'date':
                return decode_date(var)
            if ty == 'timedelta':
                return dt.timedelta(seconds=var)
            if ty == 'null':
                return None
            if ty == 'list':
                return read_list(fp,var,fp[group].attrs['list_length'])
            if ty == 'array':
                subtype = decode_str(fp[group].attrs['subtype'])
                if subtype in ['int','float','bool']:
                    return var
                if subtype == 'str':
                    new_var = [decode_str(v) for v in var.ravel()]
                    new_var = np.array(new_var,dtype='U')
                    new_var.shape = var.shape
                    return new_var
                if subtype == 'null':
                    new_var = np.array([None]*var.size)
                    if var.shape:
                        new_var.shape = var.shape
                    return new_var
                if subtype == 'date':
                    new_var = np.array([decode_date(v) for v in var.ravel()])
                    new_var.shape = var.shape
                    return new_var
                if subtype == 'timedelta':
                    new_var = np.array([dt.timedelta(seconds=v) for v in var.ravel()])
                    new_var.shape = var.shape
                    return new_var
                raise Exception('Unable to process hdf5.Dataset of type %s, subtype %s' % (ty,subtype))
            raise Exception('Unable to process hdf5.Dataset of type %s' % ty)
        else: # it's a group
            var = {} # make a dict, return that
            for key in fp[group]:
                var[key] = read_recursive(fp,group+'/'+key)
            if ty == 'dict':
                return var
            if ty == 'table':
                names = list(var.keys())
                # sort by column order
                order = [fp[group+'/'+key].attrs['column'] for key in names]                
                names = [n for o, n in sorted(zip(order, names))] # sort names by order
                dtype = []
                data = []
                for key in names:
                    data.append(var[key])
                    if var[key].ndim > 1: # nested, include shape after 1st dim
                        dtype.append((key,var[key].dtype,var[key].shape[1:]))
                    else:
                        dtype.append((key,var[key].dtype))
                var = np.fromiter(list(zip(*data)),dtype=dtype)
                return var
        
    with h5py.File(filename,'r') as fp:
        var = read_recursive(fp,'/var')
    return var

def test():
    
    def recursive_equals(a,b,path='/var'):
        if isinstance(a,np.ndarray):
            if not isinstance(b,np.ndarray):
                print('%s are not both numpy ndarrays ' % path)
                return False
            if len(a.dtype) != len(b.dtype):
                print('%s have different dtype structure' % path)
                return False
            if len(a.dtype) == 0: # simple array
                if np.array_equal(a,b,equal_nan=np.issubdtype(a.dtype,float)):
                    return True
                else:
                    print('arrays %s are not equal ' % path)
                    return False
            else: # structured array
                if a.dtype.names != b.dtype.names:
                    print('structured arrays %s have different keys ' % path)
                    return False
                for key in a.dtype.names:
                    if not np.array_equal(a[key],b[key],equal_nan=np.issubdtype(a[key].dtype,float)):
                        print(key,a[key].dtype,np.issubdtype(a[key].dtype,float))
                        print('arrays %s/%s are not equal' % (path,key))
                        return False
                return True
                
        if isinstance(a,dict):
            if not isinstance(b,dict):
                print('%s are not both dicts ' % path)
                return False
            if not recursive_equals(sorted(a.keys()),sorted(b.keys()),path):
                print('dicts %s have different keys ' % path)
                return False
            for key in a:
                if not recursive_equals(a[key],b[key],path+'/' + key):
                    return False
            return True
        
        if isinstance(a,list):
            if not isinstance(b,list):
                print('%s are not both lists ' % path)
                return False
            if len(a) != len(b):
                print('lists %s have different length ' % path)
                return False
            for (i,ai) in enumerate(a):
                if not recursive_equals(ai,b[i],path+'/[' + str(i) + ']'):
                    return False
            return True

        if a is None:
            if b is None:
                return True
            else:
                print('%s are not both None' % path)
                return False

        for base in (int,float,bool,str,dt.datetime,dt.timedelta):
            if np.issubdtype(type(a),base):
                if not np.issubdtype(type(b),base):
                    print('%s are not both %s' % (path,base.__name__))
                    return False
                if (not np.isscalar(a)) or not np.isscalar(b):
                    return False
                if (base == float) and np.isnan(a) and np.isnan(b):
                    return True
                if a == b:
                    return True
                else:
                    print('%s are unequal' % path)
                    return False        
        
        print('Unable to compare %s types %s and %s' % (path,type(a),type(b)))
        return False
            
    fails = 0 # number of failures
    # set up tests, each test has a tag and variable
    # variable will be written to test_<tag>.h5, then read and the results compared
    tests = {
            'dict' : {'string':'foo','integer':3,'real':4.5,'boolean':True,
                      'date':dt.datetime.now(),
                      'list':['foo',3,4.5,True,None],
                      'dictionary':{'a':1,'b':2.5,'c':'bar'}},
            'list' : ['foo',3,4.5,True,{'a':1,'b':2.5,'c':'bar'},None],
            'scalar' : 4.5,
            'string' : 'Hello World',
            'none' : None,
            'nan' : np.nan,
            'null' : np.array([None]*4),
            'empty_list' : [],
            'empty_array' : np.array([]),
            'array' : np.array([[1,np.nan,2],[3.0,4.1,5]]),
            'table' : np.array([(1,2.0,'foo',True),(-1,3.5,'bar',False)],
                                dtype=[('integer',int),('real',float),
                                       ('string','U3'),('boolean',bool)]),
            'nested' : np.array([(1,2.0,'foo',True,np.zeros(3)),(-1,3.5,'bar',False,np.ones(3))],
                                 dtype=[('integer',int),('real',float),
                                        ('string','U3'),('boolean',bool),
                                        ('subarray',float,(3,3))]),
            }
    for key,var in tests.items():
        filename = 'test_%s.h5' % key        
        var2hdf5(var,filename)
        var2 = hdf52var(filename)
        if recursive_equals(var,var2):
            print('%s: PASS' % key)
        else:
            print('%s1:' % key,var,'\n%s2:' % key,var2)
            print('%s: FAIL' % key)
            fails += 1

    print('TOTAL FAILS: %d' % fails)        

    # now do IDL tests    
    idl_count = 0
    for file in os.listdir():
        if not file.endswith('_idl.h5'): continue
        print('Reading IDL file %s' % file)
        var = hdf52var(file)
        idl_count += 1
    print('Total IDL files read: %d' % idl_count)

if __name__ == '__main__':
    test()    
