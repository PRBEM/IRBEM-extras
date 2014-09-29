;+
;Procedure:
;   UBK_COTRANS
;
;Purpose:
;   Coordinate transform between GSM and SM using GEOPACK
;
;Calling Sequence:
;   UBK_COTRANS, XIN, YIN, ZIN, YDHMS, $
;       /GSM2SM, /SM2GSM, N_THREADS=N_THREADS, $
;       XOUT=XOUT, YOUT=YOUT, ZOUT=ZOUT
;
;Arguments:
;   XIN, YIN, ZIN:
;       Input coordinate vectors. Required.
;   YDHMS:
;       5 element array of [YEAR, DOY, HOUR, MINUTE, SECOND]. Required.
;
;Keywords:
;   /GSM2SM or /SM2GSM: Conversion from GSM to SM or from SM to GSM. Only one should be set.
;   N_THREADS = Integer value of the number of threads to assign the load. Default is 8. If set, the value should be an integer of value greater than or equal to 1.
;
;   [XYZ]OUT = Named variables in which to return the transformed vectors. The dimensions are the same as XIN.
;
;Example:
;   UBK_COTRANS, XO=XO, YO=YO, ZO=ZO, [6,-6], [0, 0], [0, 0], [2001, 1, 1, 1, 1], $
;       /GSM2SM, N_THREADS=2
;
;Notes:
;   This routine internally used GEOPACK further information of which can be found at http://geo.phys.spbu.ru/~tsyganenko/modeling.html.
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

PRO UBK_COTRANS, XIN, YIN, ZIN, YDHMS, $
    GSM2SM=GSM2SM, SM2GSM=SM2GSM, N_THREADS=N_THREADS, $
    XOUT=XOUT, YOUT=YOUT, ZOUT=ZOUT

    ; Varify inputs
    if (n_elements(xin) ne n_elements(yin)) or $
        (n_elements(xin) ne n_elements(zin)) then $
        message, "The input vector elements are not in equal length."

    if n_elements(YDHMS) ne 5 then $
        message, "The number of elements of YDHMS variable is incorrect."

    if (keyword_set(GSM2SM) and keyword_set(SM2GSM)) or $
        (not keyword_set(GSM2SM) and not keyword_set(SM2GSM)) then $
        message, "Only one of GSM2SM and SM2GSM can be set."
    to_co_system = keyword_set(GSM2SM)

    if not keyword_set(N_THREADS) then N_THREADS = 8
    if (n_elements(N_THREADS) ne 1 or N_THREADS lt 1) then $
        message, "N_THREADS should be an integer greater than 0."

    ; Call DLM
    ubk_cotrans_c, xout, yout, zout, xin, yin, zin, $
        YDHMS, to_co_system, N_THREADS

    ; Reform
    dim = size(xin, /dim)
    if dim then begin
        xout = reform(xout,dim, /overwrite)
        yout = reform(yout,dim, /overwrite)
        zout = reform(zout,dim, /overwrite)
    endif

end
