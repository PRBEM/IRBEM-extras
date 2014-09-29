
;
; $Author$
; $LastChangedDate$
; $Revision$
; $Id$
;

x = dblarr(17) - 6;
y = x * 0
pa = (dindgen(17)/16*80+10) * !pi/180
date = [2001,1,1,1,1]
internal = 'igrf'

parmod = [3,-10,1,-5, 1,1,1,1,1,1]
external = 'ts05'

n_phi = 180

r = sqrt(x^2 + y^2)
phi = atan(y,x)
phi(where(phi<0)) = 2*!pi + phi(where(phi<0))

message,/continue, "L* Cyl"
ubk_lstar, r, phi, pa, $
date, parmod, external, internal, $
/preserve_drift, $
lstar=lstar, k=k, phi0=phi0, $
driftshell=c, footpoints=f

plot,c.r*cos(c.phi),c.r*sin(c.phi)

print, "L*=", lstar, ", K=", k, ", Phi0=", phi0

end
