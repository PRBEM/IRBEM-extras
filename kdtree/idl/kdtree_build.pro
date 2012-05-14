FUNCTION kdtree_build,X,LIB_PATH=inLIB_PATH,STOREX=inSTOREX
  ; tree = kdtree_build(X) - build kdtree
  ; X is an Nx x Nc double matrix
  ;  of Nx points in an Nc-dimensional space
  ; tree is an anonymous structure that describes
  ;  the kdtree for fast nearest-neighbors look-up
  ; KEYWORDS:
  ; LIB_PATH, string, contains extra paths in 
  ;  which to look for kdtree.so or kdtree.dll
  ;  path_list is delimited by the system's path
  ;  separator (; or :)

  if size(X,/tname) ne 'DOUBLE' then begin
    Xdbl = double(X) ; convert X to double, try again
    tree = kdtree_build(Xdbl,inLIB_PATH)
    return,tree
  endif
  
  lib_file = kdtree_load(inLIB_PATH) ; load the kdtree library
  
  if !VERSION.MEMORY_BITS eq 64 then ITYPE = 15 else ITYPE = 13
  ; ITYPE = 13; unsigned long int = ULONG = type 13 on win32
  ; ITYPE = 15; ULONG64 on linux64
  
  Xsize = size(X)
  Nx = fix(Xsize(1),type=ITYPE)
  Nc = fix(Xsize(2),type=ITYPE)
  
  root = fix(0,type=ITYPE)
  c = intarr(Nx)
  parent = make_array(Nx,TYPE=ITYPE)
  left = make_array(Nx,TYPE=ITYPE)
  right = make_array(Nx,TYPE=ITYPE)
    
  flags = fix(0) ; IDL is effectively column major because it indexes [column,row] and calls itself row major (fix = int)
  
  result = call_external(lib_file,'kdtree_build_idl',X,Nx,Nc,flags, $
    root,c,parent,left,right,value=bytarr(9))    
    
  tree = create_struct('root',root,'c',c,'parent',parent, $
    'left',left,'right',right,'X',X,'Nx',Nx,'Nc',Nc)

  return,tree
  
END

