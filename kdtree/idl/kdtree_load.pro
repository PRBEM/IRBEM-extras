function kdtree_load,LIB_PATH
  ; lib_file = kdtree_load()
  ; lib_file = kdtree_load(LIB_PATH)
  common kdtree_load_common,lib_file
  
  if n_elements(lib_file) eq 0 then begin
    case !VERSION.OS_FAMILY OF
      'unix' : ext = '.so'
      'Windows' : ext = '.dll'
    endcase
    dist_type = getenv('DISTTYPE') ; .el5 or .el6, look for kdtree.el5.so or kdtree.el6.so
    if strlen(dist_type) ne 0 then ext = dist_type + ext
    base_file = 'kdtree' + ext
    filesep = path_sep() ; / or \
    pathsep = path_sep(/search_path) ; : or ;
    paths = strsplit(!PATH,pathsep,/extract) ; search IDL path
    
    if n_elements(LIB_PATH) ne 0 then begin
      ; also search LIB_PATH if supplied
      paths = [strsplit(LIB_PATH,pathsep,/extract),paths]
    endif
    lib_file = file_which(strjoin(paths,pathsep),base_file,/include_current_dir)
    if n_elements(lib_file) eq 0 then begin
      message, 'Unable to locate ' + base_file
    endif
  endif
  
  return,lib_file
end
