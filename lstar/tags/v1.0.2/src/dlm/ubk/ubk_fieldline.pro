;+
;Procedure:
;   UBK_FIELDLINE
;
;Purpose:
;   To calculate field line parameters (SM coordinate only)
;
;Calling Sequence:
;   UBK_FIELDLINE, X0, Y0, Z0, $
;       YDHMS, IOPT_PARMOD, EXTERNAL, INTERNAL, $
;       IONOR=IONOR, DS=DS, N_THREADS=N_THREADS, $
;       USE_INTERP_FIELD_WITH_POLY_ORDER=USE_INTERP_FIELD_WITH_POLY_ORDER, $
;       XYZMEQ=XYZMEQ, BMEQ=BMEQ, XYZEQ=XYZEQ, BEQ=BEQ, $
;       K=K, BMIRROR=BM, XFL=XFL, YFL=YFL, ZFL=ZFL, $
;       XYZFOOT=XYZFOOT
;
;Arguments:
;   X0, Y0, Z0:
;       Input coordinate vectors in RE. Required.
;   YDHMS:
;       5 element array of [YEAR, DOY, HOUR, MINUTE, SECOND]. Required.
;   IOPT_PARMOD: External input parameters.
;       IF EXTERNAL=='T89' THEN IOPT_PARMOD = IOPT:
;           T89 input equivalent to Kp.
;           IOPT= 1       2        3        4        5        6        7
;           KP  = 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  >=6-
;       ELSEIF EXTERNAL>='T96' AND EXTERNAL<='TS05' THEN IOPT_PARMOD = PARMOD:
;           10 element array for T96 model and above.
;           PARMOD=[PDYN (nPa), DST (nT), BYIMF, BZIMF (nT), W1, W2, W3, W4, W5, W6]
;       ELSEIF EXTERNAL=='TS07' THEN IOPT_PARMOD = PARMOD:
;           102 element array for TS07D model.
;           PARMOD=[PDYN (nPa), 101 element fitting coefficients]
;       ELSE IOPT_PARMOD is not used.
;   EXTERNAL:
;       External field model. Currently, "T89", "T96", "T02", "TS05", "TS07" models are built in. Pass "none" not to use the external field model.
;   INTERNAL:
;       Internal magnetic field model. Either "DIP" or "IGRF" (case-insensitive). Required.
;
;Keywords:
;   IONOR = The radial distance of the ionospheric boundary in RE. Default is 1.015 RE (~ 100 km above Earth).
;   DS = Field line integration step size in RE. Default is 0.05 RE.
;   N_THREADS = Integer value of the number of threads to assign the load. Default is 8. If set, the value should be an integer of value greater than or equal to 1.
;   USE_INTERP_FIELD_WITH_POLY_ORDER = Either 1 (linear) or 2 (2nd order).
;       If set, caches the field vectors at discrete points and interpolates the
;       field vector at the given location from the cached vectors.
;       The interpolation may increase the calculation speed if the underlying
;       field model is expensive, but the result includes the numerical error
;       from the interpolation and the interpolated vector violates divergence
;       free condition. Default is unset (equivalently pass 0).
;
;   XYZMEQ = Named variable in which to return the coordinates of the magnetic equator (the magnetic equator is not necessary where B=Bmin) in RE. The dimensions are [3, n_elements(X0)].
;   BMEQ = Named variable in which to return the magnetic field magnitude at XYZMEQ in nT. n_elements(X0) element array.
;   XYZEQ = Named variable in which to return the coordinates of the coordinate equator where the field lines thread in RE. The dimensions are [3, n_elements(X0)].
;   BEQ = Named variable in which to return the magnetic field magnitude at XYZEQ in nT. n_elements(X0) element array.
;   K = Named variable in which to return the modified second invariants at given field lines (refer to Roederer 1970). The dimensions are [M, n_elements(X0)] where M is the maximum number of nodes along all field lines. NAN's are filled if the field line is shorter. nT^.5 RE.
;   BMIRROR = Named variable in which to return the magnetic mirror magnitudes corresponding to K in nT. Refer to the K keyword.
;   [XYZ]FL = Named variables in which to return the field line coordinates corresponding to K in RE. Refer to the K keyword.
;   XYZFOOT = Named variables in which to return the foot points of the field lines in RE. The dimensions are [3, n_elements(X0)].
;
;Example:
;   UBK_FIELDLINE, [-6,6], [0,0], [0,0], $
;       [2001,1,1,1,1], 0, "none", "dip", $
;       N_THREADS=2, $
;       XYZMEQ=XYZMEQ, BMEQ=BMEQ, XYZEQ=XYZEQ, BEQ=BEQ, $
;       K=K, BMIRROR=BM, XFL=XFL, YFL=YFL, ZFL=ZFL, $
;       XYZFOOT=XYZFOOT
;
;Notes:
;   Please refer for the field models to http://geo.phys.spbu.ru/~tsyganenko/modeling.html and http://geomag_field.jhuapl.edu/model/.
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

PRO UBK_FIELDLINE, X0, Y0, Z0, $
    YDHMS, IOPT_PARMOD, EXTERNAL, INTERNAL, $
    IONOR=IONOR, DS=DS, N_THREADS=N_THREADS, $
    USE_INTERP_FIELD_WITH_POLY_ORDER=USE_INTERP_FIELD_WITH_POLY_ORDER, $
    XYZMEQ=XYZMEQ, BMEQ=BMEQ, XYZEQ=XYZEQ, BEQ=BEQ, $
    K=K, BMIRROR=BM, XFL=XFL, YFL=YFL, ZFL=ZFL, $
    XYZFOOT=XYZFOOT

    ; Varify inputs
    if (n_elements(x0) ne n_elements(y0)) or $
        (n_elements(x0) ne n_elements(z0)) then $
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
        "ts07": begin
            EXTERNAL1 = 5
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
    endif else if (EXTERNAL1 ge 2 and EXTERNAL1 le 4) then begin
        if (n_elements(IOPT_PARMOD) ne 10) then $
            message, "Invalid PARMOD parameter."
    endif else if (EXTERNAL1 eq 5) then begin
        if (n_elements(IOPT_PARMOD) ne 102) then $
            message, "Invalid PARMOD parameter."
    endif

    if not keyword_set(IONOR) then IONOR = 1.015
    if (n_elements(IONOR) ne 1 or IONOR lt 1d) then $
        message, "IONOR should be greater than or equal to 1 RE."

    if not keyword_set(DS) then DS = 0.05
    if (n_elements(DS) ne 1 or DS le 0d) then $
        message, "DS should be greater than 0 RE."

    if not keyword_set(N_THREADS) then N_THREADS = 8
    if (n_elements(N_THREADS) ne 1 or N_THREADS lt 1) then $
        message, "N_THREADS should be an integer greater than 0."

    if not keyword_set(USE_INTERP_FIELD_WITH_POLY_ORDER) then USE_INTERP_FIELD_WITH_POLY_ORDER = 0
    if (n_elements(USE_INTERP_FIELD_WITH_POLY_ORDER) ne 1 or (USE_INTERP_FIELD_WITH_POLY_ORDER lt 0 or USE_INTERP_FIELD_WITH_POLY_ORDER gt 2)) then $
        message, "USE_INTERP_FIELD_WITH_POLY_ORDER should be between 0 and 2."

    ; Call DLM
    ubk_field_line_c, K, Bm, XYZMEQ, BMEQ, XYZEQ, BEQ, XYZFOOT, XFL, YFL, ZFL, $
	    X0,Y0,Z0,YDHMS, IOPT_PARMOD, EXTERNAL1, INTERNAL1, $
        IONOR, DS, N_THREADS, USE_INTERP_FIELD_WITH_POLY_ORDER

end
