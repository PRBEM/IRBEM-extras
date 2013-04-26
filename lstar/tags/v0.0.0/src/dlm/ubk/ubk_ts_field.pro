;+
;Procedure:
;   UBK_TS_FIELD
;
;Purpose:
;   To calculate magnetic field using Tsyganenko field models
;
;Calling Sequence:
;   UBK_TS_FIELD, X, Y, Z, $
;       YDHMS, IOPT_PARMOD, EXTERNAL, INTERNAL, $
;       /GSM, /SM, N_THREADS=N_THREADS, $
;       BX=BX, BY=BY, BZ=BZ
;
;Arguments:
;   X, Y, Z:
;       Input coordinate vectors in RE. Required.
;   YDHMS:
;       5 element array of [YEAR, DOY, HOUR, MINUTE, SECOND]. Required.
;   IOPT_PARMOD: External input parameters.
;       IF EXTERNAL=='T89' THEN IOPT_PARMOD = IOPT:
;           T89 input equivalent to Kp.
;           IOPT= 1       2        3        4        5        6        7
;           KP  = 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  >=6-
;       ELSEIF EXTERNAL>='T96' THEN IOPT_PARMOD = PARMOD:
;           10 element array for T96 model and above.
;           PARMOD=[PDYN (nPa), DST (nT), BYIMF, BZIMF (nT), W1, W2, W3, W4, W5, W6]
;       ELSE IOPT_PARMOD is not used.
;   EXTERNAL:
;       External field model. Currently, "T89", "T96", "T02" and "TS05" models are built in. Pass "none" not to use the external field model.
;   INTERNAL:
;       Internal magnetic field model. Either "DIP" or "IGRF" (case-insensitive). Required.
;
;Keywords:
;   /GSM and /SM: Either GSM or SM coordinate to use. Only one should be set.
;   N_THREADS = Integer value of the number of threads to assign the load. Default is 8. If set, the value should be an integer of value greater than or equal to 1.
;
;   B[XYZ] = Named variables in which to return the calculated magnetic field vectors in nT. Same as X.
;
;Example:
;   IOPT = 2
;   UBK_TS_FIELD, [-6,6], [0,0], [0,0], $
;       [2001,1,1,1,1], IOPT, "t89", "dip", $
;       N_THREADS=2, /GSM, $
;       BX=BX, BY=BY, BZ=BZ
;
;Notes:
;   1. Please refer for the field models to http://geo.phys.spbu.ru/~tsyganenko/modeling.html.
;   2. Solar wind velocity is fixed to [-400, 0, 0] km/s.
;
;Modifications:
;
;Author:
;   Kyungguk Min
;-

;
; $Author$
; $LastChangedDate$
; $Revision$
; $Id$
;

PRO UBK_TS_FIELD, X, Y, Z, $
    YDHMS, IOPT_PARMOD, EXTERNAL, INTERNAL, $
    GSM=GSM, SM=SM, N_THREADS=N_THREADS, $
    BX=BX, BY=BY, BZ=BZ

    ; Varify inputs
    if (n_elements(X) ne n_elements(Y)) or $
        (n_elements(X) ne n_elements(Z)) then $
        message, "The input vector elements are not in equal length."

    if n_elements(YDHMS) ne 5 then $
        message, "The number of elements of YDHMS variable is incorrect."

    case strlowcase(INTERNAL) of
        "dip": INTERNAL1 = 0
        "igrf": INTERNAL1 = 1
        else: begin
            message, "Invalid INTERNAL parameter."
            end
    endcase

    switch strlowcase(EXTERNAL) of
        "none": begin
            EXTERNAL1 = 0
            break
            end
        "t89": begin
            EXTERNAL1 = 1
            break
            end
        "t96": begin
            EXTERNAL1 = 2
            break
            end
        "t02": begin
            EXTERNAL1 = 3
            break
            end
        "ts05": begin
            EXTERNAL1 = 4
            break
            end
        else: begin
            message, "Invalid EXTERNAL parameter."
            break
            end
    endswitch

    if (EXTERNAL1 eq 1) then begin
        if (n_elements(IOPT_PARMOD) ne 1) or $
            (IOPT_PARMOD lt 1 or IOPT_PARMOD gt 7) then $
            message, "Invalid IOPT parameter."
    endif else if (EXTERNAL1 ge 2) then begin
        if (n_elements(IOPT_PARMOD) ne 10) then $
            message, "Invalid PARMOD parameter."
    endif

    if not keyword_set(N_THREADS) then N_THREADS = 8
    if (n_elements(N_THREADS) ne 1 or N_THREADS lt 1) then $
        message, "N_THREADS should be an integer greater than 0."

    if (keyword_set(GSM) and keyword_set(SM)) or $
        (not keyword_set(GSM) and not keyword_set(SM)) then $
        message, "Only one of GSM and SM can be set."
    co_system = keyword_set(SM)

    ; Call DLM
    ubk_ts_field_c, BX, BY, BZ, X, Y, Z, $
	    YDHMS, IOPT_PARMOD, EXTERNAL1, INTERNAL1, $
	    co_system, N_THREADS

    ; Reform
    dim = size(x, /dim)
    if dim then begin
        bx = reform(bx,dim, /overwrite)
        by = reform(by,dim, /overwrite)
        bz = reform(bz,dim, /overwrite)
    endif

end
