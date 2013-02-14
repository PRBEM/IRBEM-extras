FORWARD_FUNCTION kdtree_buid
FUNCTION kdtree_kNN,tree,X0,k,R2=inR2,XNN=inXNN,DistScale=inDistScale,LIB_PATH=inLIB_PATH,MAX_ABS=inMAX_ABS
  ; index = kdtree_kNN(tree,X0,k) - % retrieve index and k points nearest to X0  
  ; tree is an anonymous structure created by kdtree_build
  ; X0 is an Nx0 x Nc matrix of Nx0 query points in an Nc-dimensional space 
  ; k is the number of neighbors to find for each query point
  ; KEYWORDS:
  ; DistScale: (array of length Nc), defines weighting for each dimension
  ; XNN: provides a variable (must be initialized to 1) in which to store the actual locations
  ;  of the nearest neighbors. On return, the size will be Nx0 x Nc x k
  ; /inMAX_ABS: use max absolute distance rather than pythagorean
  ; LIB_PATH:  string, contains extra paths in 
  ;  which to look for kdtree.so or kdtree.dll
  ;  path_list is delimited by the system's path
  ;  separator (; or :)

  if size(X0,/tname) ne 'DOUBLE' then begin
    X0dbl = double(X0) ; convert X0 to double, try again
    index = kdtree_kNN(tree,X0dbl,k,XNN=inXNN,R2=inR2,DistScale=inDistScale,LIB_PATH=inLIB_PATH)
    return,index
  endif
  
  ; now sure X0 is double
  
  lib_file = kdtree_load(inLIB_PATH) ; load the kdtree library

  X0size = size(X0) ; size is [ndims, Nx0, Nc,...]
  Nx0 = ulong(X0size(1))
  
  if X0size(2) ne tree.Nc then message,"X0 must have same number of columns as tree.X"
    
  if !VERSION.MEMORY_BITS eq 64 then ITYPE = 15 else ITYPE = 13
  ; ITYPE = 13; unsigned long int = ULONG = type 13 on win32
  ; ITYPE = 15; ULONG64 on linux64
  
  flags = fix(0) ; IDL is effectively column major because it indexes [column,row] and calls itself row major (fix = int)
  
  if keyword_set(inMAX_ABS) then flags = flags + 2 ; max absolute distance (otherwise pythagorean)
  
  ce_k = fix(k,type=ITYPE)
  if n_elements(inDistScale) eq 0 then begin
    ce_DistScale = 1.0+dblarr(tree.Nc)
  endif else begin
    ce_DistScale = double(inDistScale)
  endelse
  
  index = make_array(Nx0,k,type=ITYPE)
  R2 = dblarr(Nx0,k)
  
  result = call_external(lib_file,'kdtree_kNN_idl',tree.X,tree.Nx,tree.Nc,flags, $
    tree.root,tree.c,tree.parent,tree.left,tree.right,X0,NX0,ce_k,ce_DistScale,index,R2,value=bytarr(15))
  index = index-fix(1,type=ITYPE) ; tree indexing is 1-based, idl is zero based
    
  if keyword_set(inR2) then begin
    inR2 = R2
  endif
    
  if keyword_set(inXNN) then begin
    inXNN = dblarr(Nx0,tree.Nc,k)
    for ik = 0,k-1 do begin
      inXNN[*,*,ik] = tree.X[index[*,ik],*] 
    endfor
  endif

  return,index
  
END

