
;
; $Author$
; $LastChangedDate$size(xin, /dim)
; $Revision$
; $Id$
;

xgsm = [6, -6, 0, 0]
ygsm = [0, 0, 6, -6]
zgsm = [0, 0, 0, 0]

date = [2001,1,1,1,1]

ubk_cotrans, xgsm, ygsm, zgsm, date, $
    /gsm2sm, n_threads=n_threads, $
    xout=xsm, yout=ysm, zout=zsm

ubk_cotrans, xsm, ysm, zsm, date, $
    /sm2gsm, n_threads=n_threads, $
    xout=xgsm1, yout=ygsm1, zout=zgsm1

print, "Difference of abs(gsm - gsm1): "
print, abs(xgsm1-xgsm)
print, abs(ygsm1-ygsm)
print, abs(zgsm1-zgsm)

end
