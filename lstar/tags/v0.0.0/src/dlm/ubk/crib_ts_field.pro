
;
; $Author$
; $LastChangedDate$
; $Revision$
; $Id$
;

x = [6, -6]
y = [0, 0]
z = [0, 0]

date = [2001,1,1,1,1]
internal = "dip"

parmod = [3,-10,1,-5, 1,1,1,1,1,1]
external = "ts05"
n_threads = 2

ubk_ts_field, x, y, z, $
date, parmod, external, internal, $
/gsm, n_threads=n_threads, $
bx=bx, by=by, bz=bz

print, "bx = " + string(bx)
print, "by = " + string(by)
print, "bz = " + string(bz)
print, "|b| = " + string(sqrt(bx^2 + by^2 + bz^2))

end
