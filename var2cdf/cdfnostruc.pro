PRO CDFNOSTRUC, id, numAtts, numZvars, attNames, varNames, varType, VarElem, Data
;*
;* Procedure cdfnostruc.pro
;*
;* Written by	Edith L. Mazur
;*		Aerospace Corporation
;*		Space Instrumentation Department
;*		November 9, 2008
;*
;* This routine is called by cdf2var.pro
;* It handles situations in which the data in the CDF was not originally
;* formatted as a structure.
;*

zz = 0
FOR zz = 0, numZvars - 1 DO BEGIN
    vname = varNames(zz)
    vtype = varType(zz)
    IF (vtype EQ 'CDF_UCHAR') OR (vtype EQ 'CDF_CHAR') THEN BEGIN
       CDF_VARGET, id, vname, var_cdf, /STRING
       numElems = N_ELEMENTS(var_cdf)
    ENDIF ELSE BEGIN
       CDF_VARGET, id, vname, var_cdf
    ENDELSE
    IF numZvars EQ 1 THEN BEGIN
       data = var_cdf
       RETURN
    ENDIF ELSE BEGIN
       IF zz EQ 0 THEN BEGIN
          commandString = 'data = CREATE_STRUCT("' + varNames(zz) + '", var_cdf)'
          result        = EXECUTE(commandString)
       ENDIF ELSE BEGIN
          commandString = 'data = CREATE_STRUCT(data, "' + varNames(zz) + '", var_cdf)'
          result        = EXECUTE(commandString)
       ENDELSE
    ENDELSE
ENDFOR

END
