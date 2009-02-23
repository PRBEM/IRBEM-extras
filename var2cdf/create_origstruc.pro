PRO CREATE_ORIGSTRUC, id, varName, varType, varElem, data
;*
;* The Zvariable names need to be decoded and then used to recreate the
;* original data structure used in either Matlab or IDL.
;* A period is used when defining subelements of a structure.  The number
;* of periods in the Zvariable names will be used in helping to determine the
;* format of the structure.
;*
;* For example:  if there is no period "." in the name of the variable, then
;* that is the name  of the main structure.  One "."  means that
;* it is a variable of the main structure. Two or more "."s refers to nested
;* structures and their elements within the main structure.  
;*
;* This method is currently assuming that only one main structure and one substructure
;* were written to the CDF. Modifications will be needed if more than one main structure
;* or substructure exists in the CDF.
;* 
FIGURE_OUT_STRUCTURE, varName, varType, main_struc, sub_struc, main_type, sub_type
n_main = N_ELEMENTS(main_struc)
n_sub  = N_ELEMENTS(sub_struc)

FOR mm = 1, n_main - 1 DO BEGIN
    name     = main_struc(mm)
    type     = main_type(mm)

    ;*
    ;* Convert name to a format that can be used in a structure
    result = STRSPLIT(name, '.', /EXTRACT)
    sname  = RESULT(1)

    ;* 
    ;* Load in the data for eacg variable and the put it into the structure
    ;* by using "create_struct"
    IF (type EQ 'CDF_UCHAR') OR (type EQ 'CDF_CHAR') THEN BEGIN
       CDF_VARGET, id, name, var, /STRING
       numElems = N_ELEMENTS(var)
    ENDIF  ELSE BEGIN
       CDF_VARGET, id, name, var
    ENDELSE

    IF mm EQ 1 THEN $
       data = CREATE_STRUCT(name = main_struc(0), sname, var)
    IF mm GT 1 THEN $
       data = CREATE_STRUCT(data, sname, var) 
ENDFOR

;*
;* If the number of elements in struc2 is only equal to one, then there is no substructure
;* in the original data that was used to create the CDF.  Therefore, this next portion of the
;* program can be skipped.
;*
IF n_sub EQ 1 THEN GOTO, skip_struc2
 
FOR ss = 1, n_sub - 1 DO BEGIN
    name = sub_struc(ss)
    type = sub_type(ss)

    ;*
    ;* Convert name to a format that can be used in a structure
    result = STRSPLIT(name, '.', /EXTRACT)
    sname  = RESULT(2)

    ;* 
    ;* Load in the data for each variable and the put it into the structure
    ;* by using "create_struct"
    IF (type EQ 'CDF_UCHAR') OR (type EQ 'CDF_CHAR') THEN BEGIN
       CDF_VARGET, id, name, var, /STRING
       numElems = N_ELEMENTS(var)
    ENDIF  ELSE BEGIN
       CDF_VARGET, id, name, var
    ENDELSE

    IF ss EQ 1 THEN $
       subdata = CREATE_STRUCT(name = sub_struc(0), sname, var)
    IF ss GT 1 THEN $
       subdata = CREATE_STRUCT(subdata, sname, var) 
ENDFOR

;* Once both the main structure and sub structures are created, add the sub structures
;* to the main structure by using create_struc
result = STRSPLIT(sub_struc(0), ".", /EXTRACT)
subname = result(1)
data = CREATE_STRUCT(data, subname, subdata)

skip_struc2:
RETURN
END
