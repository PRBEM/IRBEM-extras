FUNCTION kdtree_kNN,tree,X0,k,R2=inR2,treeX=inX,XNN=inXNN,DistScale=inDistScale,LIB_PATH=inLIB_PATH,MAX_ABS=inMAX_ABS
  ; index = kdtree_kNN(tree,X0,k) - % retrieve index and k points nearest to X0  
  ; tree is an anonymous structure created by kdtree_build
  ; KEYWORDS:
  ; treeX: a Nx x Nc matrix
  ;  of Nx points in an Nc-dimensional space (only needed 
  ;  if tree was built w/o /STOREX switch. Use when X is 
  ;  not inherently double precision)
  ; DistScale: (Nc x 1), defines weighting for each dimension
  ; /inMAX_ABS: use max absolute distance rather than pythagorean
  ; LIB_PATH:  string, contains extra paths in 
  ;  which to look for kdtree.so or kdtree.dll
  ;  path_list is delimited by the system's path
  ;  separator (; or :)

  if size(X0,/tname) ne 'DOUBLE' then begin
    X0dbl = double(X0) ; convert X0 to double, try again
    index = kdtree_kNN(tree,X0dbl,k,XNN=inXNN,R2=inR2,X=inX,DistScale=inDistScale,LIB_PATH=inLIB_PATH)
    return,index
  endif
  
  ; now sure X0 is double
  
  if keyword_set(inX) then begin
    if size(inX,/tname) ne 'DOUBLE' then begin
      inXdbl = double(inX) ; convert inX to double, try again
      index = kdtree_kNN(tree,X0,k,XNN=inXNN,R2=inR2,X=inXdbl,DistScale=inDistScale,LIB_PATH=inLIB_PATH)
      return,index
    endif
  endif
  
  ; now, either inX is double or tree.X is  
  
  lib_file = kdtree_load(inLIB_PATH) ; load the kdtree library

  Nx = n_elements(tree.c)
  X0size = size(X0)
  Nx0 = ulong(X0size(1))
  Nc = ulong(X0size(2))
    
  ITYPE = 13; unsigned long int = ULONG = type 13
  
  flags = fix(0) ; IDL is effectively column major because it indexes [column,row] and calls itself row major (fix = int)
  
  if keyword_set(inMAX_ABS) then flags = flags + 2 ; max absolute distance (otherwise pythagorean)
  
  if keyword_set(inX) then X = inX else X = tree.X
  
  ce_k = ulong(k)
  if n_elements(inDistScale) eq 0 then begin
    ce_DistScale = 1.0+dblarr(1,Nc)
  endif else begin
    ce_DistScale = double(inDistScale)
  endelse
  
  index = make_array(Nx0,k,type=ITYPE)
  R2 = dblarr(Nx0,k)
  
  result = call_external(lib_file,'kdtree_kNN_idl',X,Nx,Nc,flags, $
    tree.root,tree.c,tree.parent,tree.left,tree.right,X0,NX0,ce_k,ce_DistScale,index,R2,value=bytarr(15))
  index = index-ulong(1) ; tree indexing is 1-based, idl is zero based
    
  if keyword_set(inR2) then begin
    inR2 = R2
  endif
    
  if keyword_set(inXNN) then begin
    inXNN = dblarr(Nx0,Nc,k)
    for ik = 0,k-1 do begin
      inXNN[*,*,ik] = X[index[*,ik],*] 
    endfor
  endif

  return,index
  
END

