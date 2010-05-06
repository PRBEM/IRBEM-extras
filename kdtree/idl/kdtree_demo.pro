pro kdtree_demo
  ; demonstrate kdtree

  Nx = 1000
  Nc = 2
  Nx0 = 10
  k = 5
  
  X = randomu(randomu_seed,Nx,Nc,/double)
  tree = kdtree_build(X)
  X0 = randomu(randomu_seed,Nx0,Nc,/double)
  XNN = 1 ; set to true, will be returned as 3-d matrix
  index = kdtree_kNN(tree,X0,k,treeX=X,XNN=XNN)
  
  plot,X[*,0],X[*,1],PSYM=3
  oplot,X0[*,0],X0[*,1],PSYM=2
  for i = 0,k-1 do begin
    oplot,XNN[*,0,i],XNN[*,1,i],PSYM=6
  endfor
  
  return
end
