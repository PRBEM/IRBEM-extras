pro kdtree_demo
  ; demonstrate kdtree

  Nx = 1000 ; number of points in data cloud
  Nc = 2 ; number of dimensions of data cloud
  Nx0 = 10 ; number of query points for Nearest Neighbors search
  k = 5 ; number of neighbors to find for each query point
  
  X = randomu(randomu_seed,Nx,Nc,/double) ; random Nx x Nc point cloud
  tree = kdtree_build(X)
  X0 = randomu(randomu_seed,Nx0,Nc,/double); random query locations Nx0 x Nc
  XNN = 1 ; set to true, will be returned as 3-d matrix
  index = kdtree_kNN(tree,X0,k,XNN=XNN)
  
  ; make a plot: dots mark point cloud, stars mark query points, squares mark neighbors
  plot,X[*,0],X[*,1],PSYM=3
  oplot,X0[*,0],X0[*,1],PSYM=2
  for i = 0,k-1 do begin
    oplot,XNN[*,0,i],XNN[*,1,i],PSYM=6
  endfor
  
  return
end
