PRO CDFSTRUC, id, numAtts, numZvars, attNames, varNames, varType, VarElem, Data
;*
;* Procedure CDFSTRUC.PRO
;*
;* Written by	Edith L. Mazur
;*		Aerospace Corporation
;*		Space Instrumentation Department
;*		November 9, 2008
;*
;* This procedure is called by program CDF2VAR.  It handles loading in the cdf
;* files where the data needs to be in a structureal format.
;*
;* Figure out # of structures
nV         = N_ELEMENTS(varNames)
CDF_CONTROL, id, attribute = 'STRUCTURE_NAMES', get_attr_info = ginfo

;*  Determine the name of the main structure and its substructures.
numStrucs  = ginfo.numgentries
strucNames = STRARR(numStrucs)


FOR jj = 0, numStrucs - 1 DO BEGIN
    CDF_ATTGET, id, 'STRUCTURE_NAMES', jj, value, cdf_type = type
    strucNames(jj) = value
ENDFOR

;* Put the structure names in order from the deepest nested structure
;* to the stop structure.
strucNames = strucNames[REVERSE(SORT(strucNames))]

;*
;* This portion of the program loads the data and creates the structure at the same
;* time.  The "EXECUTE" command is used to recreate exactly the original name of the 
;* structure.  
;*
strucFlag = 1
FOR kk = 0, numStrucs - 1 DO BEGIN
    result     = STRSPLIT(strucNames(kk), '.', COUNT = count, /EXTRACT)
    sname      = strucNames(kk)

    ;*
    ;* Determine the variable type for each element of the structure
    ;* This is necessary when creating the structure.
    wtype      = WHERE(strucNames(kk) EQ varNames)
    stype      = varType(wtype)

    teststring = sname + '.*'
    test       = STRMATCH(varNames, teststring)
    wtest      = WHERE(test eq 1, wcount)

    ;* 
    ;* Here is where the variables included in each structure are determined.
    FOR mm = 0, wcount - 1 DO BEGIN
        vname    = varNames(wtest(mm))
        vtype    = varType(wtest(mm))
        vnresult = STRSPLIT(vname, '.', COUNT = vncount)
        IF vncount EQ (count + 1) THEN BEGIN
           ;* 
           ;* check to make sure variable is not the name of one of the structures.
           wstrucName = WHERE(strucNames EQ vname, scount)

           IF scount GE 0 THEN BEGIN
              ;*
              ;* Load variable names into an array
              sname = [[sname], vname]
              stype = [[stype], vtype]
           ENDIF
        ENDIF
    ENDFOR
    ;*
    ;* This is where we create the structure and load the variables
    ;* into it.
    numVars      = N_ELEMENTS(sname)
    dummy        = strucNames(kk)
    dummy        = sname(0)
    result_dummy = STRSPLIT(dummy, '.', COUNT = dcount, /EXTRACT)
    dummy        = result_dummy(dcount - 1)
    
    IF kk EQ 0 THEN begin
       storeStructs = dummy    
       level1 = dummy
    ENDIF
    IF kk GT 0 THEN BEGIN
       storeStructs = [[storeStructs], dummy]
    ENDIF
    IF numVars EQ 1 THEN BEGIN
       strucFlag = 0
       GOTO, next_structure
    ENDIF
    elemNames = STRARR(numVars)
    elemTypes = STRARR(numVars)
    elemNames(0) = dummy
    FOR nn = 1, numVars - 1 DO BEGIN
        Ename = sname(nn)
        type = stype(nn)

        ;* convert name to a format that can be used as an element
        ;* in the structure
        result_element = STRSPLIT(Ename, '.', /EXTRACT)
        element_name   = result_element(dcount)
        elemNames(nn)  = element_name
        elemTypes(nn)  = type

        ;*
	;* check to see if the variable to be loaded in is actually a
	;* stucture or substructure.  If so, do not load the data.
        ;* This is important because once a structural elements is
	;* defined, it cannot be changed to another element type.
        w1 = WHERE(Ename EQ strucNames, w1_count)
        IF w1_count NE 0 THEN GOTO, DONT_LOAD_DATA

        ;* Use varget to load in the data for each variable of the structure.
        ;*
        IF (type EQ 'CDF_UCHAR') or (type EQ 'CDF_CHAR') THEN BEGIN
           CDF_VARGET, id, Ename, var_cdf, /STRING
           numElems = N_ELEMENTS(var_cdf)
        ENDIF ELSE BEGIN
           CDF_VARGET, id, Ename, var_cdf
        ENDELSE
        commandString = element_name + '=  + var_cdf'
        Eresult       = EXECUTE(commandString) 
    ENDFOR
    DONT_LOAD_DATA:
    ;*
    ;* Once all the variables within the  have been loaded, they can be put into
    ;* the structure.
    ;* There is a possible problem, though.  If an element of the structure has the same name
    ;* as the structure or substructure, this causes a problem with the data being lost once
    ;* the structure is created.  This problem is solved by temporarily renaming the structure
    ;* until the data elements are loaded.  Then the structure's name is restored.
    ;*
    mainName = elemNames(0)
    subnames = elemNames(1:numVars-1)
    wmatch   = WHERE(subnames EQ mainName, match_cnt)
    IF match_cnt EQ 1 THEN BEGIN
       elemNames(0) = mainName + '0'
       level1 = level1 + '0'
       wstore = where(storestructs EQ mainName)
       storeStructs(wstore) = elemNames(0)
    ENDIF
    FOR pp = 1, numVars - 1 DO BEGIN
       IF pp EQ 1 THEN BEGIN
         commandString = elemNames(0) + '= CREATE_STRUCT(elemNames(pp), ' + elemNames(pp) + ')'
         result = EXECUTE(commandString)
       ENDIF 
       IF pp GT 1 THEN BEGIN
          tagname = elemNames(pp)
          IF tagname NE level1 THEN BEGIN
              commandString = elemNames(0) + '= CREATE_STRUCT('+ elemNames(0)+ ',"' $
                            + elemNames(pp) +'", ' + elemNames(pp) + ')'
              result = EXECUTE(commandString)
          ENDIF
          IF tagname EQ level1 THEN BEGIN
              commandString = elemNames(0) + '= CREATE_STRUCT(' + elemNames(0) + ',"'$
	      		    + storeStructs(kk-1) + '", ' + storeStructs(kk-1) + ')'
              result = EXECUTE(commandString)
          ENDIF
       ENDIF
    ENDFOR

    ;* Restore original structure name if a temporary one had been used.
    IF match_cnt EQ 1 THEN begin
       commandstring = mainName + '=' + elemNames(0)
       result = EXECUTE(commandstring) 
    endif 

    level1 = elemNames(0)
    next_structure:
ENDFOR
;*
;* Combine main structure with its substructures.
mainName = storestructs(kk-1)
command_string =  'data = ' + mainName
eresult = EXECUTE(command_string)
END
