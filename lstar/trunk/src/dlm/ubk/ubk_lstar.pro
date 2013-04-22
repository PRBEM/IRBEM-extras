;+
;Procedure:
;   UBK_LSTAR
;
;Purpose:
;   To calculate the generalized L value, L* (SM coordinate only)
;
;Calling Sequence:
;   UBK_LSTAR, XorR0, YorPHI0, PA0, $
;       YDHMS, IOPT_PARMOD, EXTERNAL, INTERNAL, $
;       IONOR=IONOR, DS=DS, DXorDR=DXorDR, DYorDPHI=DYorDPHI, $
;       N_PHI=N_PHI, N_THETA=N_THETA, /CARTESIAN_GRID, /PRESERVE_DRIFT_CONTOUR, N_THREADS=N_THREADS, $
;       LSTAR=LSTAR, K=K, PHI0=PHI0, $
;       DRIFTSHELL_CONTOUR=DRIFTSHELL, FOOTPOINTS_CONTOUR=FOOTPOINTS
;
;Arguments:
;   XorR0, YorPHI0:
;       Equatorial input coordinate vectors in the cartesian coordinate system if /CARTESIAN_GRID is set. Otherwise the inputs should be R and Phi. Required. RE and radian.
;   PA0:
;       Initial pitch angles in radian. the number of elements should be the same as the coordinate inputs. Required.
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
;   IONOR = The radial distance of the ionospheric boundary in RE. Default is 1.015 RE (~ 100 km above Earth).
;   DS = Field line integration step size in RE. Default is 0.1 RE.
;   DXorDR = R (or X if /CARTESIAN_GRID is set) grid resolution in RE. Default is 0.2 RE.
;   DYorDPHI = PHI (or Y if /CARTESIAN_GRID is set) grid resolution in radian (or RE). Default is 0.2 RE if /CARTESIAN_GRID is set or !pi*2/N_PHI.
;   N_PHI = The number of azimuthal grids to discretize Earth's surface to integrate the magnetic flux. Default is 360/5.
;   N_THETA = The number of polar grids to discretize Earth's surface to integrate the magnetic flux. Default is 180*2.
;   /CARTESIAN_GRID: If set, the routine uses the cartesian grid. Otherwise cylindrical grid. The inputs/outputs should follow the same coordinate system except for FOOTPOINTS_CONTOUR which is always in spherical coordinate system.
;   /PRESERVE_DRIFT_CONTOUR: If set, the routine calculates the drift contours and returns the result through FOOTPOINTS_CONTOUR and DRIFTSHELL_CONTOUR.
;   N_THREADS = Integer value of the number of threads to assign the load. Default is 8. If set, the value should be an integer of value greater than or equal to 1.
;
;   LSTAR = Named variable in which to return the calculated L*. The dimensions are the same as XorR0.
;   K = Named variable in which to return the modified second invariants. The dimensions are the same as XorR0. nT^.5 RE.
;   PHI0 = Named variable in which to return Earth's magnetic flux, i.e., 2pi B0 where B0 is Earth's magnetic moment. nT RE^2.
;   DRIFTSHELL_CONTOUR = Named variable in which to return the drift contours. Structure variable of {X:[M, n_elements(XorR0)], Y:[M, n_elements(YorPHI0)]}, where M is the maximum number of nodes of all contours. NAN's are filled for the shorter contours if /CARTESIAN_GRID is set. If not, the returned values are in the cylindrical coordinate system and "X" and "Y" tags are replaced with "R" and "PHI" tags.
;   FOOTPOINTS_CONTOUR = Named variable in which to return the foot points of the drift contours. Structure variable of {PHI:[M, n_elements(XorR0)], THETA:[M, n_elements(YorPHI0)]}, where M is the maximum number of nodes of all contours.
;
;Example:
;   UBK_LSTAR, [-6,6], [0,0], [90,90] * !pi/180, $
;       [2001,1,1,1,1], 0, "none", "dip", $
;       N_THREADS=2, $
;       LSTAR=LSTAR, DRIFTSHELL=CNTR
;   plot,CNTR.x, CNTR.y
;
;Notes:
;   Cautious about memory pressure when /PRESERVE_DRIFT_CONTOUR is set.
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

PRO UBK_LSTAR, XorR0, YorPHI0, PA0, $
    YDHMS, IOPT_PARMOD, EXTERNAL, INTERNAL, $
    IONOR=IONOR, DS=DS, DXorDR=DXorDR, DYorDPHI=DYorDPHI, $
    N_PHI=N_PHI, N_THETA=N_THETA, CARTESIAN_GRID=CARTESIAN_GRID, N_THREADS=N_THREADS, $
    PRESERVE_DRIFT_CONTOUR=PRESERVE_DRIFT_CONTOUR, $
    LSTAR=LSTAR, K=K, PHI0=PHI0, $
    DRIFTSHELL_CONTOUR=DRIFTSHELL, FOOTPOINTS_CONTOUR=FOOTPOINTS

    ; Varify inputs
    if (n_elements(XorR0) ne n_elements(YorPHI0)) or $
        (n_elements(XorR0) ne n_elements(PA0)) then $
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

    if not keyword_set(IONOR) then IONOR = 1.015
    if (n_elements(IONOR) ne 1 or IONOR lt 1d) then $
        message, "IONOR should be greater than or equal to 1 RE."

    if not keyword_set(DS) then DS = 0.1
    if (n_elements(DS) ne 1 or DS le 0d) then $
        message, "DS should be greater than 0 RE."

    if not keyword_set(N_PHI) then N_PHI = 360/5
    if (n_elements(N_PHI) ne 1 or N_PHI le 10) then $
        message, "N_PHI should be greater than 10."

    if not keyword_set(N_THETA) then N_THETA = 180*2
    if (n_elements(N_THETA) ne 1 or N_THETA le 10) then $
        message, "N_THETA should be greater than 10."

    if not keyword_set(DXorDR) then DXorDR = 0.2
    if (n_elements(DXorDR) ne 1 or DXorDR le 0d) then $
        message, "DXorDR should be greater than 0 RE."

    if not keyword_set(DYorDPHI) then $
        if keyword_set(CARTESIAN_GRID) then $
            DYorDPHI = 0.2 $
            else DYorDPHI = !pi*2 / N_PHI
    if (n_elements(DYorDPHI) ne 1 or DYorDPHI le 0d) then $
        message, "DYorDPHI should be greater than 0 RE or 0 radian."

    if not keyword_set(N_THREADS) then N_THREADS = 8
    if (n_elements(N_THREADS) ne 1 or N_THREADS lt 1) then $
        message, "N_THREADS should be an integer greater than 0."

    ; Call DLM
    ubk_lstar_c, LSTAR,K,PHI0,XorRc,YorPhic, PF,TF, $
            XorR0,YorPHI0,PA0, $
            YDHMS,IOPT_PARMOD,EXTERNAL1,INTERNAL1, $
            IONOR,DS,DXorDR,DYorDPHI,N_PHI,N_THETA, $
            keyword_set(PRESERVE_DRIFT_CONTOUR), $
            keyword_set(CARTESIAN_GRID), N_THREADS
        
    ; Reform
    dim = size(XorR0, /dim)
    if dim then begin
        LSTAR = reform(LSTAR,dim, /overwrite)
        K = reform(K,dim, /overwrite)
    endif

    if keyword_set(PRESERVE_DRIFT_CONTOUR) then begin
        if keyword_set(CARTESIAN_GRID) then DRIFTSHELL = {X:XorRc, Y:YorPhic} $
        else DRIFTSHELL = {R:XorRc, PHI:YorPhic}
        FOOTPOINTS = {PHI:PF, THETA: TF}
    endif

end
