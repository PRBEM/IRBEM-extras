
;
; $Author$
; $LastChangedDate$
; $Revision$
; $Id$
;

x0 = [6, -6]
y0 = [0, 0]
z0 = [0, 0]

date = [2001,1,1,1,1]
internal = "igrf"

iopt = 3
external = "t89"
ionoR = 1.1
n_threads = 2

ubk_fieldline, x0, y0, z0, $
    date, iopt, external, internal, $
    ionor=ionor, n_threads=n_threads, $
    xyzmeq=xyzmeq, bmeq=bmeq, xyzeq=xyzeq, beq=beq, $
    k=k, bmi=bm, xfl=xfl, yfl=yfl, zfl=zfl, $
    xyzfoot=xyzfoot

;plot, k, bm
plot, xfl[*,0], zfl[*,0]
oplot, -xfl[*,1], zfl[*,1], color=2l^16-1

end
