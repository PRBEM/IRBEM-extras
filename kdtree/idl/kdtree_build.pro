FUNCTION kdtree_build,X,LIB_PATH=inLIB_PATH
  ; tree = kdtree_build(X) - build kdtree

  if size(X,/tname) ne 'DOUBLE' then begin
    Xdbl = double(X) ; convert X to double, try again
    tree = kdtree(Xdbl,inLIB_PATH)
    return,tree
  endif
  
  lib_file = kdtree_load(inLIB_PATH) ; load the kdtree library
  
  Xsize = size(X)
  Nx = ulong(Xsize(1))
  Nc = ulong(Xsize(2))
  
  ITYPE = 13; unsigned long int = ULONG = type 13
  
  root = ulong(7)
  c = intarr(Nx)
  parent = make_array(Nx,TYPE=ITYPE)
  left = make_array(Nx,TYPE=ITYPE)
  right = make_array(Nx,TYPE=ITYPE)
    
  flags = fix(1) ; IDL is row major (fix = int)
  
  result = call_external(lib_file,'kdtree_build_idl',X,Nx,Nc,flags, $
    root,c,parent,left,right,value=bytarr(9),/unload)    
    
  tree = create_struct('root',root,'c',c,'parent',parent, $
    'left',left,'right',right)

  return,tree
  
END

