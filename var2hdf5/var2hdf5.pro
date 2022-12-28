; var2hdf5.pro - write variables to and read variables from hdf5 files
; Prinicpal author: Paul O'Brien
; public functions:
;    var2hdf5 - save a variable to an HDF5 file
;    hdf52var - read a variable from an HDF5 file (for files created by var2hdf5)
;    test - run several test cases that write, read, and compare
    
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
  key = strmid(i,2)
  while strlen(key) lt padding do key = '0'+key
  return,key
end

function _build_simple_var,variable,name,lists
  ; var = _build_simple_var(variable,name,lists)
  ; builds a structure that holds needed info for h5_create
  ; to create this simple variable
  ; simple variables are just scalars or arrays
  ; lists tracks list variables
  _init_var2hdf5
  common VAR2HDF5_COMMON
  
  tname = typename(variable)
  if isa(variable,/boolean) then type = 'bool' $
  else if tname eq 'STRING' then begin
    type = 'str'
    ; now test for date-time string in ISO8601 YYYY-MM-DDTHH:MM:SS... format
    iso8601 = '^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}'
    if where(stregex(variable,iso8601,/BOOLEAN) eq !FALSE,/NULL) eq !NULL then type = 'datetime'
  endif else if ismember(tname,INT_TYPES) then type = 'int' $
  else if ismember(tname,FLOAT_TYPES) then type = 'float' $
  else message,"Variable " + variable + " has unknown simple type " + tname
  
  var = {_NAME:name,_DATA:variable,_TYPE:'DATASET',$
    type:{_NAME:'type',_DATA:type,_TYPE:'ATTRIBUTE'}}
  return, var
end

function _build_var,variable,name,lists
  ; var = _build_var(variable,name,lists)
  ; builds a structure that holds needed info for h5_create
  ; to create this simple or complex variable
  ; lists tracks list variables
  
  _init_var2hdf5
  common VAR2HDF5_COMMON
  
  tname1 = size(variable,/tname) ; simple type names
  if ismember(tname1,SIMPLE_TYPES) then begin
    return,_build_simple_var(variable,'var',lists)
  endif
  
  tname2 = typename(variable) ; a bit more verbose on complex types (tname1 OBJREF)
  
  ; Cannot handle: tname1=COMPLEX, DCOMPLEX, POINTER, OBJREF
  ; Can handle: tname1=STRUCT, tname2=LIST, HASH, DICTIONARY, ORDEREDHASH
  if tname1 eq 'STRUCT' then begin
    var = {_NAME:variable,_TYPE:'GROUP'}
    tags = tag_names(variable)
    for i = 0,n_elements(tags)-1 do begin
      val = _build_var(variable.(i),tags[i],lists)
      var = create_struct(tags[i],val,var)
    endfor
  endif else if tname2 eq 'LIST' then begin
    val = list()
    list_len = n_elements(variable)
    padding = ceil(alog10(list_len)) ; required digits to express all entries
    for j=0,list_len-1 do begin
      key = pad_int(j,padding) ; left pad with zeros
      val.add,_build_var(variable[i],key,lists)
    endfor
    lists.add,val
  endif else if tname2 eq 'HASH' then begin
  endif else if tname2 eq 'DICTIONARY' then begin
  endif else if tname2 eq 'ORDEREDHASH' then begin
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
  
  ; should be able to simply set conv_func = lambda(x:x)
  ; but IDL just throws errors when calling lambda(5)
  ; because it doesn't recognize lambda as a function pointer
  ; executing examples from the help reveals exactly this error
  ; IDL is buggy
  has_converter = n_elements(conv_func) eq 0 
  
  ; remove file if it exists
  if file_test(filename) then begin
    file_delete,filename
  endif
  
  
  lists = list()
  ; each entry in lists is a list of built-out variables (w/ metadata)
  struct = _build_var(variable,'var',lists)
  
  if n_elements(lists) gt 0 then begin
    for i = 0,n_elements(lists)-1 do begin
      lname = '__list_' + strmid(i,2) + '__'
      listi = lists[i]
      list_len = n_elements(listi)
      s = {_NAME:lname,_TYPE:'GROUP',$
           type:{_TYPE:'ATTRIBUTE',_DATA:'list'},$
           list_length:{_TYPE:'ATTRIBUTE',_DATA:list_len}}
      for j=0,list_len-1 do begin
        s = create_struct('x'+listsi(j)._NAME,listi(j),s) ; prepend x to make valid IDL variable
      endfor ; j - index into list[i]
      struct = create_struct(lname,s,struct)
    endfor ; i - index into lists    
  endif
  
  top = {_TYPE:'GROUP',$
    VAR2HDF5_SPEC_VERSION:{_TYPE:'ATTRIBUTE',_DATA:VAR2HDF5_SPEC_VERSION},$
    VAR:struct}

  h5_create,filename,top
end

function _read_var,rawvar
; var = _read_var(rawvar)
; recursively read and process rawvar
; can't do this now because IDL < 8.8.1 has a bug that prevents
; reading string attributes from HDF5 files. This is critical because
; the string attribute "type" is used to tell the reader what
; type the variable is supposed to have
; IDL is so frustrating
end


function hdf52var,filename,converter
  ;var = hdf52var(filename,converter=None)
  ;read variable from hdf5 file
  ;filename - HDF5 filename to read
  ;converter - function pointer that converts variables to HDF5-recognized types
  ;calls var = conveter(var,attrs) before trying to handle the variable
  ;var - variable that was read
  ;see README for conversion spec
  
  _init_var2hdf5
  common VAR2HDF5_COMMON
  h5data = h5_parse(filename,/read_data) ; read the entire hdf5 file
  if h5data.VAR2HDF5_SPEC_VERSION._DATA gt VAR2HDF5_SPEC_VERSION then begin
    filev = strtrim(h5data.VAR2HDF5_SPEC_VERSION._DATA,2) ; conver to string, w/o padding
    codev = strtrim(VAR2HDF5_SPEC_VERSION,2) ; conver to string, w/o padding
    message,'VAR2HDF5_SPEC_VERSION '+ filev + '>' + codev  +' not suported'
  endif
  return, _read_var(h5data.var)
end

pro test_hdf5


  tests = { $; each entry is test_name:expected_value
    dict : {string:'foo',integer:3,real:4.5,boolean:!TRUE, $
      date:timestamp(),list:list('foo',3,4.5,!TRUE,!VALUES.F_NAN), $
      dictionary:{a:1,b:2.5,c:'bar'}}, $
    list : list('foo',3,4.5,boolean(0),{a:1,b:2.5,c:'bar'},!VALUES.F_NAN), $
    scalar : 4.5, $
    string : 'Hello World', $    
    none : list(), $ ; list() has n_elements eq 0, serves as python None 
    nan : !VALUES.F_NAN, $
    null : make_array(4,value=!VALUES.F_NAN), $
    empty_list : list(), $
    empty_array : list(), $ ; UNDEFINED not allowed as dict entry value
    array : [[1,!VALUES.F_NAN,2],[3.0,4.1,5]], $
    table : {integer:[1,-1],real:[2.0,3.5],string:['foo','bar'],boolean:[!TRUE,!FALSE]}, $
    nested : {integer:[1,-1],real:[2.0,3.5],string:['foo','bar'],boolean:[!TRUE,!FALSE],subarray:make_array(3,3,value=1)} $ 
  }
  
  ; must run python test first to generate test files to read
  ; python files are lowercase <tag>.h5 
  
  tnames = [tag_names(tests)]
  for i = 0,n_elements(tnames)-1 do begin
    filename = 'test_' + strlowcase(tnames[i]) + '.h5'
    print,filename
  endfor
  
end


