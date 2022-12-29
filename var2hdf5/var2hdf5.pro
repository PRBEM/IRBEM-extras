; var2hdf5.pro - write variables to and read variables from hdf5 files
; Prinicpal author: Paul O'Brien
; public main functions:
;    var2hdf5 - save a variable to an HDF5 file
;    hdf52var - read a variable from an HDF5 file (for files created by var2hdf5)
;    test - run several test cases that write, read, and compare
; public utility functions:
;    ismember - determine whether an item is in a list
;    pad_int - left zero pad an integer to a specified strlen
    
; placeholder for var2hdf5.pro which will provide two functions

pro _init_var2hdf5
  ; initialize VAR2HDF5_COMMON common block
  common VAR2HDF5_COMMON,VAR2HDF5_SPEC_VERSION,INT_TYPES,FLOAT_TYPES,SIMPLE_TYPES
  VAR2HDF5_SPEC_VERSION = 0 ; integer
  INT_TYPES = ['BYTE','INT','UINT','LONG','ULONG','LONG64','ULONG64']
  FLOAT_TYPES = ['FLOAT','DOUBLE']
  SIMPLE_TYPES = [INT_TYPES,FLOAT_TYPES,'STRING']
end

function ismember,entry,list
  ; bool = ismember(entry,list)
  ; returns true if list is found anywhere (exactly) in list
  return, where(list eq entry,/NULL) ne !NULL
end

function pad_int,i,padding
  ; str = pad_int(i,padding)
  ; creates string from integer i
  ; left pads with zeros until strlen(str) >= padding
  key = strtrim(i,2)
  while strlen(key) lt padding do key = '0'+key
  return,key
end

function _build_simple_var,variable,name,conv_func,lists
  ; var = _build_simple_var(variable,name,conv_func,lists)
  ; builds a structure that holds needed info for h5_create
  ; to create this simple variable
  ; simple variables are just scalars or arrays
  ; conv_func is the converter provided to var2hdf5
  ; lists tracks list variables
  _init_var2hdf5
  common VAR2HDF5_COMMON
  
  if conv_func then variable = conv_func(variable)
  tname = typename(variable)
  if isa(variable,/boolean) then type = 'bool' $
  else if tname eq 'STRING' then begin
    type = 'str'
    ; now test for date-time string in ISO8601 YYYY-MM-DDTHH:MM:SS... format
    iso8601 = '^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}'
    if where(stregex(variable,iso8601,/BOOLEAN) eq !FALSE,/NULL) eq !NULL then type = 'date'
  endif else if ismember(tname,INT_TYPES) then type = 'int' $
  else if ismember(tname,FLOAT_TYPES) then type = 'float' $
  else message,"Variable " + name + " has unknown simple type " + tname
  
  subtype = type
  if n_elements(variable) gt 1 then type = 'array'
  
  var = {_NAME:name,_DATA:variable,_TYPE:'DATASET',$
    type:{_NAME:'type',_DATA:type,_TYPE:'ATTRIBUTE'}}
    
  if type eq 'array' then begin
    ; add subtype attribute
    var = create_struct('subtype',{_NAME:'subtype',_TYPE:'ATTRIBUTE',_DATA:subtype},var)
  endif
    
  return, var
end

function _build_var,variable,name,conv_func,lists
  ; var = _build_var(variable,name,conv_func,lists)
  ; builds a structure that holds needed info for h5_create
  ; to create this simple or complex variable
  ; conv_func is the converter provided to var2hdf5
  ; lists tracks list variables
  
  _init_var2hdf5
  common VAR2HDF5_COMMON
  
  if conv_func then variable = conv_func(variable)
  
  dims = size(variable,/dimensions)
  if n_elements(dims) gt 1 then begin
    ; reverse order of multi-dimensional arrays for writing to file
    variable = transpose(variable,reverse(indgen(n_elements(dims))))
  endif
  
  tname1 = size(variable,/tname) ; simple type names
  if ismember(tname1,SIMPLE_TYPES) then begin
    return,_build_simple_var(variable,name,conv_func,lists)
  endif

  tname2 = typename(variable) ; a bit more verbose on complex types (tname1 OBJREF)
  
  ; Cannot handle: tname1=COMPLEX, DCOMPLEX, POINTER, OBJREF
  ; Can handle: tname1=STRUCT, tname2=LIST, HASH, DICTIONARY, ORDEREDHASH
  if tname2 eq 'LIST' then begin
    list_len = n_elements(variable)
    ; add list_length attribute
    var = {_NAME:name,_TYPE:'DATASET',_DATA:n_elements(lists), $
      type:{_NAME:'type',_TYPE:'ATTRIBUTE',_DATA:'list'}, $
      list_length:{_NAME:'list_length',_TYPE:'ATTRIBUTE',_DATA:n_elements(lists)}}
    val = list()
    padding = ceil(alog10(1>list_len)) ; required digits to express all entries (1> means at least 1)
    for j=0,list_len-1 do begin
      key = pad_int(j,padding) ; left pad with zeros
      val.add,_build_var(variable[j],key,conv_func,lists)
    endfor
    lists.add,val
  endif else if (tname1 eq 'STRUCT') or ismember(tname2,['DICTIONARY','HASH','ORDEREDHASH']) then begin
    var = {_NAME:name,_TYPE:'GROUP'}
    type = 'table'
    nested = !FALSE
    if tname2 ne 'ORDEREDHASH' then type = 'dict' ; without ordering, not a table
    if tname1 eq 'STRUCT' then tags = tag_names(variable) $
      else tags = variable.keys()
    for i = 0,n_elements(tags)-1 do begin
      if tname1 eq 'STRUCT' then val = variable.(i) $
        else val = variable[tags[i]]
      dims = size(val,/dimensions) ; not yet processed, so dimensions in original order
      if i eq 0 then common_length = dims[0]
      if dims[0] ne common_length then type = 'dict' ; 1st dim is not common_length long
      if (n_elements(dims) gt 1) or (~ismember(typename(val),SIMPLE_TYPES)) then nested = !TRUE ; multidimensional or not simple type
      val = _build_var(val,tags[i],conv_func,lists) ; convert to h5_create required struct
      var = create_struct(tags[i],val,var)
    endfor
    if type eq 'table' then begin
      ; add nested attribute
      newvar = {nested:{_NAME:'nested',_TYPE:'ATTRIBUTE',_DATA:nested}}
      tags2 = tag_names(var)
      for i =0,n_elements(tags2)-1 do begin
        val = var.(i)
        if (size(val,/tname) eq 'STRUCT') && (val._TYPE ne 'ATTRIBUTE') then begin
          ; add column attribute
          column = where(tags eq tags2[i])
          val = create_struct('column',{_NAME:'column',_TYPE:'ATTRIBUTE',_DATA:column[0]},val)
        endif
        newvar = create_struct(tags2[i],val,newvar)         
      endfor
      var = newvar
    endif
    ; add type attribute
    var = create_struct('type',{_NAME:'type',_TYPE:'ATTRIBUTE',_DATA:type},var)
  endif else message, "Unable to save " + name + " of type " + tname1 + '/' + tname2
  
  return, var
end

pro var2hdf5,variable,filename,converter=conv_func
  ; var2hdf5,var,filename
  ; var2hdf5,var,filename,converter=conv_func
  ; write variable to hdf5 file
  ; var - variable to write
  ; filename - HDF5 filename to write
  ; converter - function pointer that (conv_func) converts variables to HDF5-recognized types
  ; calls var = conv_func(var) before trying to handle the variable
  ; see README for conversion spec
  _init_var2hdf5
  common VAR2HDF5_COMMON
  
  if n_elements(conv_func) eq 0 then conv_func = !FALSE
  
  ; remove file if it exists
  if file_test(filename) then begin
    file_delete,filename
  endif
  
  lists = list()
  ; each entry in lists is a list of built-out variables (w/ metadata)
  struct = _build_var(variable,'var',conv_func,lists)
  
  top = {_TYPE:'GROUP',$
    VAR2HDF5_SPEC_VERSION:{_TYPE:'ATTRIBUTE',_DATA:VAR2HDF5_SPEC_VERSION},$
    VAR:struct}

  if n_elements(lists) gt 0 then begin
    for i = 0,n_elements(lists)-1 do begin
      lname = '__list_' + strtrim(i,2) + '__'
      listi = lists[i]
      list_len = n_elements(listi)
      s = {_NAME:lname,_TYPE:'GROUP'}
      for j=0,list_len-1 do begin
        s = create_struct('_'+listi[j]._NAME,listi[j],s) ; prepend _ to make valid IDL variable name
      endfor ; j - index into list[i]
      top = create_struct(lname,s,top)
    endfor ; i - index into lists    
  endif
  

  h5_create,filename,top
end

function _read_var,rawvar,h5data,conv_func
; var = _read_var(rawvar,h5data,conv_func)
; recursively read and process rawvar
; rawvar is (sub)structure returned by h5_parse
; h5data is the top level structure returned by h5_parse
; conv_func is the converter provided to var2hdf5
; can't do this for python HDF5 files because of an IDL bug
; reading string attributes from HDF5 files. This is critical because
; the string attribute "type" is used to tell the reader what
; type the variable is supposed to have
; but can read HDF5 files written by IDL

; types are: int, float, bool, str, date, timedelta, null, list, array, dict, table
  toptags = tag_names(h5data)
  subtype = rawvar.TYPE._DATA
  if subtype eq 'array' then subtype = rawvar.SUBTYPE._DATA
  if ismember(subtype,['int','float','str','date','timedelta']) then begin
    var = rawvar._DATA ; these types need no further conversion
  endif else if subtype eq 'bool' then begin
    var = boolean(rawvar._DATA)    
  endif else if subtype eq 'null' then begin
    if n_elements(rawvar._DATA) le 1 then var = list() $
      else var = make_array(size(rawvar._DATA,/dimensions),value=list())
  endif else if subtype eq 'list' then begin
    list_num = rawvar._DATA
    list_name = '__LIST_' + strtrim(list_num,2) + '__' ; __list_#__
    i = where(toptags eq list_name)
    i = i[0]
    listi = h5data.(i)
    list_len = rawvar.list_length._DATA
    if list_len eq 0 then begin
      var = list()
    endif else begin
      var = list(indgen(list_len),/extract) ; dummy entries
      padding = ceil(alog10(1>list_len)) ; required digits to express all entries (1> means at least 1)
      tags = tag_names(listi)
      for j=0,list_len-1 do begin
        key = '_'+pad_int(j,padding) ; left pad with zeros, prepend '_' as IDL does with numeric tag names
        k = where(tags eq key)
        var[j] = _read_var(listi.(k[0]),h5data,conv_func) ; replace dummy with parsed variable
      endfor
    endelse
  endif else if subtype eq 'dict' then begin
    vtags = tag_names(rawvar)
    var = hash()
    for i=0,n_elements(vtags)-1 do begin
      val = rawvar.(i)
      if (size(val,/tname) eq 'STRUCT') && (val._TYPE ne 'ATTRIBUTE') then begin
        var[val._NAME] = _read_var(val,h5data,conv_func)
      endif
    endfor
  endif else if subtype eq 'table' then begin
    keys = []
    order = []
    values = []
    vtags = tag_names(rawvar)
    for i=0,n_elements(vtags)-1 do begin
      val = rawvar.(i)
      if (size(val,/tname) eq 'STRUCT') && (val._TYPE ne 'ATTRIBUTE') then begin
        keys = [keys,val._NAME]
        order = [order,val.column._DATA]
        values = [values,_read_var(val,h5data,conv_func)] 
      endif
    endfor
    var = orderedhash(keys[order],values[order])    
  endif else begin
    message,'Unable to process (sub)type ' + subtype
  endelse

  if conv_func then var = conv_func(var,rawvar)  
  
  dims = size(var,/dimensions)
  if n_elements(dims) gt 1 then begin
    ; reverse order of multi-dimensional arrays after reading
    var = transpose(var,reverse(indgen(n_elements(dims))))
  endif
  
end


function hdf52var,filename,converter=conv_func
  ;var = hdf52var(filename,converter=conv_func)
  ;read variable from hdf5 file
  ;filename - HDF5 filename to read
  ;converter - function pointer that converts variables to HDF5-recognized types
  ;calls var = conveter(var,attrs) before trying to handle the variable
  ;var - variable that was read
  ; converter works as var = conv_func(var,rawvar) before returning
  ;see README for conversion spec
  
  if n_elements(conv_func) eq 0 then conv_func = !FALSE
  
  _init_var2hdf5
  common VAR2HDF5_COMMON
  h5data = h5_parse(filename,/read_data) ; read the entire hdf5 file
  if h5data.VAR2HDF5_SPEC_VERSION._DATA gt VAR2HDF5_SPEC_VERSION then begin
    filev = strtrim(h5data.VAR2HDF5_SPEC_VERSION._DATA,2) ; conver to string, w/o padding
    codev = strtrim(VAR2HDF5_SPEC_VERSION,2) ; conver to string, w/o padding
    message,'VAR2HDF5_SPEC_VERSION '+ filev + '>' + codev  +' not suported'
  endif
  return, _read_var(h5data.var,h5data,conv_func)
end

function _recursive_equals,var1,var2,path=in_path
  ; bool = _recursive_equals(var1,var2,path='/var')
  if n_elements(in_path) eq 0 then path = '/var'
  message,'Not implemented yet'
  
  print,'Unable to compare ',in_path,' types ',typename(var1),' and ', ypename(var2)
  return, !FALSE
end

pro test
  ; run a series of write/read tests for consistency

  tests = { $; each entry is test_name:expected_value
    struct : {string:'foo',integer:3,real:4.5,boolean:!TRUE, $
      date:timestamp(),list:list('foo',3,4.5,!TRUE,!VALUES.F_NAN), $
      dictionary:{a:1,b:2.5,c:'bar'}}, $
    dict : dictionary('string','foo','integer',3,'real',4.5,'boolean',!TRUE, $
      'date',timestamp(),'list',list('foo',3,4.5,!TRUE,!VALUES.F_NAN), $
      'dictionary',{a:1,b:2.5,c:'bar'}), $
    hash : hash({string:'foo',integer:3,real:4.5,boolean:!TRUE, $
      date:timestamp(),list:list('foo',3,4.5,!TRUE,!VALUES.F_NAN), $
      dictionary:{a:1,b:2.5,c:'bar'}}), $
    list : list('foo',3,4.5,boolean(0),{a:1,b:2.5,c:'bar'},!VALUES.F_NAN), $
    scalar : 4.5, $
    string : 'Hello World', $    
    none : list(), $ ; list() has n_elements eq 0, serves as python None 
    nan : !VALUES.F_NAN, $
    null : make_array(4,value=!VALUES.F_NAN), $
    empty_list : list(), $
    empty_array : list(), $ ; UNDEFINED not allowed as dict entry value
    array : [[1,!VALUES.F_NAN,2],[3.0,4.1,5]], $
    table : orderedhash({integer:[1,-1],real:[2.0,3.5],string:['foo','bar'],boolean:[!TRUE,!FALSE]}), $
    nested : orderedhash({integer:[1,-1],real:[2.0,3.5],string:['foo','bar'],boolean:[!TRUE,!FALSE],$
    subarray:[make_array(1,3,3,value=0),make_array(1,3,3,value=1)]}) $ 
  }
  
  tnames = tag_names(tests)
  fails = 0
  for i = 0,n_elements(tnames)-1 do begin
    ;if ismember(strlowcase(tnames[i]),['table','nested']) then continue
    key = tnames[i]
    var = tests.(i)
    filename = 'test_' + strlowcase(tnames[i]) + '_idl.h5'
    print,filename
    var2hdf5,var,filename
    var2 = hdf52var(filename)
    continue ; cannot complete this test because _reverse_equals not written
    if _recursive_equals(var,var2) then begin
      print,key,': PASS'
    endif else begin
      print,key,'1: ',var
      print,key,'2: ',var2
      fails = fails + 1
    endelse
  endfor
  print,'TOTAL FAILS: ',strtrim(fails,2)

end


