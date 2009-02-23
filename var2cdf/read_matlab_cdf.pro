;************************************************************************************************
;* PROGRAM read_matlab_cdf.pro
;*
;* Written by	Edith L. Mazur
;*		Space Instrumentation Department
;*		January 8, 2007
;*
;* This program takes a CDF that was created by Matlab and loads it into IDL.  The idea is that
;* both IDL and MATLAB use structure formatted data but an indendent platform is needed to 
;* read structures that were created in IDL and are now needed in Matlab and vice versa.  CDF is
;* used for the independent platform but it is necessary to recreate the structures that were used
;* in either IDL or Matlab.  Since CDF does not handle structures, a variable naming structure is
;* used to give clues to the original structure format.
;*
;* This program loads in a CDF that was created in Matlab, but I suppose it should work just as
;* well for a CDF originally created in IDL as well. It then recreates the data structure based 
;* upon the filenames.  How this is done is explained in CREATE_ORIGSTRUC.PRO.
;*
;**************************************************************************************************
filename = ' '
READ, 'Enter filename (including path): ', filename
;filename = '/export/home/emazur/PROGRAMS/CDF_STUFF/s33esa_ave_bs20.cdf'
;filename = 'example'

;*
;* Open the CDF and retrieve global information about the file
;*
id = CDF_OPEN(filename)
varinfo = CDF_INQUIRE(id)

;*
;* Retrieve the variable information including the variable names,
;* variable types, and the number of elements for each variable.
;*
numAtts  = varinfo.natts
numZvars = varinfo.nzvars
attNames = STRARR(3,numAtts)
varNames = STRARR(numZvars)
varType  = STRARR(numZvars)
varElem  = LONARR(numZvars)

FOR aa = 0, numAtts - 1 DO BEGIN
    CDF_ATTINQ, id, aa, name, scope, maxentry
    attNames(0,aa) = name
    attNames(1,aa) = scope
    attNames(2,aa) = maxentry
ENDFOR

FOR zz = 0, numZvars - 1 DO BEGIN
    result = CDF_VARINQ(id, zz, /zvariable)
    varNames(zz) = result.name
    varType(zz)  = result.dataType
    varElem(zz)  = result.numElem

ENDFOR

CREATE_ORIGSTRUC, id, varNames, varType, varElem, data
CDF_CLOSE, id
END
