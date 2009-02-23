pro figure_out_structure, varNames, varType, struc1, struc2, type1, type2
kk = 0
struc1 ='' 
struc2 ='' 
for kk = 0, n_elements(varnames)-1 do begin
 name = varNames(kk)
 type = varType(kk)
 result =  strsplit(name, '.', count=count, /extract)
 if count eq 1 then begin
    struc1 = name
    type1  = type
 endif
 if count eq 2 then begin
    teststring = '*' + result(1) + '.*'
    test = strmatch(varNames, teststring)
    wtest = where(test eq 1,  wcount)
    if wcount ge 1 then begin
       struc2 = name
       type2 = type
    endif
    if wcount eq 0 then begin
       struc1 = [[struc1], name]
       type1  = [[type1], type]
    endif
 endif
 if count eq 3 then begin
    struc2 = [[struc2], name]
    type2  = [[type2], type]
 endif
endfor
return
end 
