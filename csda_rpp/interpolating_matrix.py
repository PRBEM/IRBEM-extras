"""
interpolating_matrix.py
  
interpolating_matrix - produce an interpolating matrix from one grid to another
trapz_dx  - provide trapezoidal integral weights for array x
"""

import numpy as np

def interpolating_matrix(x,y,out_of_bounds=0):
    """
    A = interpolating_matrix(x,y,out_of_bounds=0)
    x - input grid
    y - output grid
    out_of_bounds - control behavior outside boundaries of x
        0,'ZERO' - pad with zeros
        'EXTRAP' - linearly extrapolate
        nan,'Nan','NAN' - NaNs
        'EXCEPTION' - raise an exception
    A interpolating matrix such that u(y) ~ A*u(x)
    where u(x) is a scalar function evaluated at each x
    A has shape (len(y),len(x))
    A*u(x) is zero outside the bounds of x, * is dot product
    """

    if isinstance(out_of_bounds,str):
        out_of_bounds = out_of_bounds.upper()
    elif np.isnan(out_of_bounds):
        out_of_bounds = 'NAN'
    elif out_of_bounds == 0:
        out_of_bounds = 'ZERO'
    assert out_of_bounds in ['ZERO','EXTRAP','EXCEPTION','NAN'],'unknown value for out_of_bounds: %s' % str(out_of_bounds)
    
    
    x = np.array(x)
    s = np.sign(x[-1]-x[0])
    xs = x*s;
    ys = np.array(y)*s
    Nx = len(x)
    Ny = len(y)
    A = np.zeros((Ny,Nx))
    
    # xs is now increasing
    # ys is multiplied by s so it can index into xs
    
    for (iy,yi) in enumerate(ys):
        ix = 0
        while (ix<Nx) and (xs[ix] <= yi):
            ix += 1
        # xs[ix] > yi
        if ix == 0: # out of range yi < xs[0]
            if out_of_bounds == 'EXCEPTION':
                raise Exception('y value %g out of bounds %g-%g' %(s*yi,s*xs[0],s*xs[-1]))
            elif out_of_bounds == 'NAN':
                A[iy,:] = np.nan
                continue
            elif out_of_bounds == 'EXTRAP':
                ix = 1 # causes extrapolation from ix=0,ix=1
            else: # ZERO
                continue # leave A[iy,:] == 0
        elif ix == Nx:
            if xs[-1] == yi:
                A[iy,-1] = 1 # all weight on last grid point
                continue
            else:  # out of range yi > xs[-1] 
                if out_of_bounds == 'EXCEPTION':
                    raise Exception('y value %g out of bounds %g-%g' % (s*yi,s*xs[0],s*xs[-1]))
                elif out_of_bounds == 'NAN':
                    A[iy,:] = np.nan
                    continue
                elif out_of_bounds == 'EXTRAP':
                    ix = Nx-1 # causes extrapolation from ix=Nx-2,ix=Nx-1
                else: # ZERO
                    continue # leave A[iy,:] == 0
        c = (yi - xs[ix-1]) / (xs[ix]-xs[ix-1])
        A[iy,ix] = c
        A[iy,ix-1] = 1-c
            
    
    return A

def trapz_dx(x):
    """
    return dx weights for x values in x
    dx follows trapezoidal rule
    """
    x = np.array(x)
    dx = np.full(x.shape,np.nan)
    # edges get one-sided half weight
    dx[0] = (x[1]-x[0])/2
    dx[-1] = (x[-1]-x[-2])/2
    # middle has centered weight
    dx[1:-1] = (x[2:]-x[:-2])/2
    return dx

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # increasing
    x = np.linspace(3,10,30)**2
    y = np.linspace(4**2,9**2,50)
    u = lambda x : np.sin(x/10)
    plt.close('all')
    
    # first case: all y within range of x
    plt.figure()
    Axy = interpolating_matrix(x,y)
    plt.plot(y,u(y),'k-',y,np.dot(Axy,u(x)),'r--')
    plt.title('All y in range of x, u = sin')
    
    # second case: y exceeds range of x. test out_of_bounds options
    plt.figure()
    plt.plot(x,u(x),'k-',label='Truth')
    for out_of_bounds in [0,'ZERO','EXTRAP','NAN','NaN',np.nan,7,'EXCEPTION']:
        try:
            Ayx = interpolating_matrix(y,x,out_of_bounds=out_of_bounds)
        except Exception as err:
            print('Exception raised when out_of_bounds=%s.' % str(out_of_bounds),err)
            continue
        plt.plot(x,np.dot(Ayx,u(y)),'--',label=str(out_of_bounds))
    plt.title('Range of y exceeds range of x, u = sin')
    plt.legend()
    
    print('dx',trapz_dx(list(range(5))))
    