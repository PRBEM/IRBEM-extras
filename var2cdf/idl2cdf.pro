;*****************
;* Program IDL2CDF
;*
;* Written by	Edith L. Mazur and Dr. Timothy Guild
;*		The Aerospace Corporation
;*		April 15, 2008
;*
;* This program takes IDL savesets where the variables are saved are 
;* in a structural format and the output is written to a CDF file so that
;* the data can be used in other languages such as Matlab.
;*
;* The naming convention used for the output variables allows the user of the
;* CDF to recreate the original structure.  This is necessary because CDF does
;* not support the use of structures.
;*
;* Modification history:	August 28, 2008 - Edith Mazur
;*				Program was modified to be a procedure with no user
;*				input.  The saveset name and variable choices are
;*				now included in the procedures's call from the main
;*				program.
;*
;*****************
PRO CREATE_CDF, saveset_name, var, strucNames, outNames, varTypes, varDims, dimensions
;*
;* Determine output filename based upon name of saveset
result  = STRSPLIT(saveset_name, '.', /EXTRACT)
outFile = result(0)
numVars = N_ELEMENTS(outNames)
maxDim  = MAX(dimensions)
id = CDF_CREATE(  outFile, [numVars, maxDim], /SINGLE_FILE, /COL_MAJOR, /CLOBBER)
;CDF_COMPRESSION,id,SET_COMPRESSION=2
;CDF_COMPRESSION,id,SET_COMPRESSION=5,SET_GZIP_LEVEL=5  ; set the compression 
; to '5' = GZIP, and set the GZIP level to '9', the max.  TBG, 3/17/09

;*
;* Make sure output names are lower case.  CDF is case sensitive.
outNames   = STRLOWCASE(outNames)
strucNames = STRLOWCASE(strucNames)

;* Write the structure names to a global attribute
;* Because a global attribute cannot be a string array, the attribute
;* is written as mutiple entries.
variable_dummy = CDF_ATTCREATE(id, 'STRUCTURE_NAMES', /GLOBAL)
FOR gg = 0, N_ELEMENTS(strucNames) - 1 DO $
    CDF_ATTPUT, id, 'STRUCTURE_NAMES', gg, strucNames(gg)

variable_dummy = CDF_ATTCREATE(id, 'CreatedBy',   /GLOBAL)
variable_dummy = CDF_ATTCREATE(id, 'CreatedAt',   /GLOBAL)
variable_dummy = CDF_ATTCREATE(id, 'CreatedFrom', /GLOBAL)

CDF_ATTPUT, id, 'CreatedBy',   0, 'Edith Mazur'
CDF_ATTPUT, id, 'CreatedAt',   0, 'The Aerospace Corporation, Chantilly'
;CDF_ATTPUT, id, 'CreatedFrom', 0,  saveset_name
; TBG:  this doesn't work, CDF_ATTPUT is looking for a data type instead of 'saveset_name'

;* Create an attribute associated with each variable containing the 
;* number of elements in the arrays.
attid = CDF_ATTCREATE(id, 'size', /VARIABLE_SCOPE)

FOR hh = 0, N_ELEMENTS(outNames) -1 DO BEGIN
    CASE varTypes(hh) OF
	1:  type = '/CDF_BYTE'
	2:  type = '/CDF_INT2'
	3:  type = '/CDF_INT4'
	4:  type = '/CDF_FLOAT'
        5:  type = '/CDF_DOUBLE'
        7:  type = '/CDF_UCHAR'
        8:  type = '/CDF_FLOAT'	;* This code really means structure but CDF doesn't support
				;* structures. A variable type is needed, though.
        12: type = '/CDF_UINT2'
        13: type = '/CDF_UINT4'
        14: type = '/CDF_INT4'
	15: type = '/CDF_UINT4'
    ENDCASE
    name = outNames(hh)
    ;* When creating the variable, IDL/CDF needs to know the number of dimensions in 
    ;* in the variable and if it varies withing that dimension.
    ndims = varDims(hh)
 
    IF ndims LT 1 THEN BEGIN
       vary_str = "['vary']"
       the_dims = 1
    ENDIF
    IF ndims GE 1 THEN BEGIN
       vary_beg  = '['
       vary_end  = ']'
       vary_mid  = "'vary'"

       FOR ii = 0, ndims - 1 DO BEGIN
           IF ii EQ 0 THEN vary_str = vary_mid
           IF ii GT 0 THEN  vary_str = vary_str + ',' + vary_mid
       ENDFOR
       vary_str = vary_beg + vary_str + vary_end

      ;* Determine the dimensions of the CDF variable being created.
      ;*
      the_dims  = dimensions(0:ndims - 1,hh)
   ENDIF
   IF  varTypes(hh) EQ 8 THEN $
      command   = 'attid = CDF_ATTCREATE(id, "' + name + '", /VARIABLE_SCOPE)'
      command   = 'varid = CDF_VARCREATE(id, "' + name + '", ' + vary_str + ', ' + type + $
		   ', /REC_VARY, /ZVARIABLE)'
   IF varTypes(hh) NE 8  THEN $
      command   = 'varid = CDF_VARCREATE(id, "' + name + '", ' + vary_str + ', ' + type + $
		   ', /REC_VARY, /ZVARIABLE, dimensions = the_dims)'
   result    = EXECUTE(command)
;   print, command
;   print, the_dims
;   print, vartypes(hh)
   CDF_ATTPUT, id, attid, varid, the_dims, /ZVARIABLE
   IF varTypes(hh) NE 8 THEN BEGIN
	command = 'CDF_VARPUT, id, varid, '+ name + ' ,/ZVARIABLE'
;        print, command
        result = EXECUTE(command)
        IF result EQ 1 THEN print, 'PUT ' + name + ' successful!'
   ENDIF
ENDFOR
print, ' '
CDF_CLOSE, id
PRINT, 'OUTPUT WRITTEN TO: ', outfile + '.cdf'
RETURN
END

PRO IDL2CDF, saveset_name, variables = variables
;**************************
;* Main Procedure IDL2CDF *
;*			  *
;**************************
variables = STRLOWCASE(variables)
variables = STRCOMPRESS(variables, /REMOVE_ALL)
variables = STRSPLIT(variables, ',', count = count, /EXTRACT)

RESTORE, saveset_name, /VERBOSE
;*
;* Check to make sure that variables requested are structures.  If they
;* are not structures, put them into one.
;*
;* First, determine the type of each variable requested.
command = 'dataInfo = SIZE(' + variables(0) + ',/STRUCTURE)'
result  = EXECUTE(command)
type    = dataInfo.type
FOR bb = 1, count - 1 DO BEGIN
    command = 'dataInfo = SIZE(' + variables(bb) + ',/STRUCTURE)'
    result  = EXECUTE(command)
    type    = [type, dataInfo.type]
ENDFOR

;*
;* Load all requested variables into a new structure.  This is necessary
;* because the requested variables can have a varying number of
;* array elements.
command = 'var = CREATE_STRUCT("' + variables(0) + '", '+ variables(0)+')'
result  = EXECUTE(command)
FOR cc = 1, count - 1 DO BEGIN
    command = 'var = CREATE_STRUCT(var, "' + variables(cc) + '", '+ variables(cc)+')'
    result  = EXECUTE(command)
    IF (result ne 1) THEN BEGIN
      STOP, 'Could not create one of the substructures!  This command failed: ' + command
    ENDIF  
ENDFOR
;*
;* Create an array that contains the tag names of all the variables. Also create
;* an array that contains the variable type of each tag name.  This is necessary
;* when using cdf_varcreate
;*

varNames = TAG_NAMES(var)
numTags  = N_TAGS(var)
level1   = 'var.'
varinfo  = SIZE(var, /STRUCTURE)
varType  = varinfo.type
varDim   = varinfo.n_dimensions
dimensions = varinfo.dimensions
strucNames = 'var'
outNames = level1 + varNames(0)
FOR dd = 0, numTags - 1 DO BEGIN
    elemName = level1 + varNames(dd)
    ;* Determine if first level tag name is a variable or nested structure
    command = 'varinfo = SIZE(' + elemName + ', /structure)'
    result  = EXECUTE(command)
    varType = [varType, varinfo.type]
    varDim  = [varDim,  varinfo.n_dimensions]
    dimensions = [ [dimensions], [varinfo.dimensions] ]
    IF dd GT 0 THEN outNames = [outNames, elemName]
    IF varinfo.type EQ 8 THEN BEGIN ;* structure tag is a nested structure
       level2 = level1 + varNames(dd)
       strucNames = [strucNames, level2]
       command = 'subNames = TAG_NAMES(' + level2 + ')'
       result  = EXECUTE(command)
       subNames = level2 + '.' + subNames
       outNames = [outNames, subNames]
       ;varType = [varType, varinfo.type]   ; TBG: seeing if this works.  
       ;*
       ;* Determine if any of the substructure variables are themselves structures.
       Nsubs = N_ELEMENTS(subNames)
       FOR ff = 0, Nsubs - 1 DO BEGIN
           command = 'varinfo = SIZE(' + subNames(ff) + ', /structure)'
           result  = EXECUTE(command)
           varType = [varType, varinfo.type]
           varDim  = [varDim,  varinfo.n_dimensions]
           dimensions = [ [dimensions], [varinfo.dimensions] ]
           IF varinfo.type EQ 8 THEN BEGIN
              level3  = subNames(ff)
              strucNames = [strucNames, level3]
              command = 'sub2Names = TAG_NAMES(' + level3 + ')'
              result  = EXECUTE(command)
              sub2Names = level3 + '.' + sub2Names
              outNames = [outNames, sub2Names]
              ;varType = [varType, varinfo.type]   
              ; TBG: Now need to loop over at least one more level of sub2Names
              N2subs = N_ELEMENTS(sub2Names)
              FOR gg = 0, N2subs-1 DO BEGIN
                 command = 'varinfo = SIZE(' + sub2Names(gg) + ', /structure)'
                 result  = EXECUTE(command)
                 varType = [varType, varinfo.type]
                 varDim  = [varDim,  varinfo.n_dimensions]
                 dimensions = [ [dimensions], [varinfo.dimensions] ]
                 IF varinfo.type EQ 8 THEN BEGIN
                    level4  = sub2Names(gg)
                    strucNames = [strucNames, level4]
                    command = 'sub3Names = TAG_NAMES(' + level4 + ')'
                    result  = EXECUTE(command)
                    sub3Names = level4 + '.' + sub3Names
                    outNames = [outNames, sub3Names]
                 ENDIF
              ENDFOR
           ENDIF
              
       ENDFOR
    ENDIF
ENDFOR
;*
;* Add var structure name to outNames
outNames = ['var', outNames]
CREATE_CDF, saveset_name, var, strucNames, outNames, varType, varDim, dimensions
END
