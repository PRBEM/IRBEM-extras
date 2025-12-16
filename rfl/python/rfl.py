"""
@file response_function_library.py
@brief Response Function Library.
@details
Primary Author: Paul O'Brien, The Aerospace Corporation (paul.obrien@aero.org)

AI Disclosure Statement: This codebase was developed with the assistance of 
generative AI tools: Astro, Open AI - GPT-4o, Nov 2025, for the purpose of 
modifying documentation to fit the Doxygen automated documentation standard. 
All AI-assisted content has been reviewed, verified, and revised by developers
to ensure accuracy, completeness, and functionality.

A note on units:
    The standard allows for L_UNIT (length unit) and E_UNIT (energy unit),
    but internally, all parameters are represented using cm and MeV.
    Units will therefore be converted on initialization/load.

@par Glossary:
    - E         : energy, MeV
    - alpha     : local pitch angle, degrees
    - beta      : gyrophase angle, degrees
    - theta     : sensor polar angle, degrees
    - phi       : sensor azimuth angle, degrees
    - alpha0    : local pitch angle of particle incident at theta=0, degrees
    - beta0     : gyrophase angle of particle incident at theta=0, degrees
    - phib      : sensor azimuth angle of magnetic field vector, degrees
    - Egrid     : energy grid, MeV
    - alphagrid : alpha grid, degrees, 1-d numpy array or scalar
    - betagrid  : beta grid, degrees, 1-d numpy array or scalar
    - thetagrid : theta grid, degrees, 1-d numpy array or scalar
    - phigrid   : phi grid, degrees, 1-d numpy array or scalar
    - tgrid     : time grid, seconds, 1-d numpy array or scalar
    - hE        : energy response, units of MeV (CROSSCALIB applied)
    - hA*       : angular response, units of cm^2(-s) (CROSSCALIB *not* 
                  applied)
    - h*        : energy-angle response, units of MeV-cm^2-sr(-s) (CROSSCALIB 
                  applied)

The ChannelResponse and its internal classes, EnergyResponse and AngleResponse,
 have implicit factory constructors via the FactoryConstructorMixin. That means
 passing the constructor a data tree that actually defines a subclass will 
return the appropriate subclass. See FactoryConstructorMixin for information on
 how to define the is_mine method in any subclasses derived from 
ChannelResponse, EnergyResponse, and AngleResponse.

@par Major Not-Yet-Implemented Issues:
    - Internal unit conversion - presently assumes cm and MeV.
    - to_dict() - ability to convert a channel response back to a dict (needed
                  for writing). Currently just returns the instantiation dict.
    - Limited ability to write response function nested dicts.
    - No ability to read .mat files.
"""

# NOTE: In this code where the docstring is *INHERIT* that means
#       the string "*INHERIT*" will be replaced by the docstring for the
#       parent class's corresponding method.

import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator

# Utilities

## @brief Degree-based trigonometric functions.
## These lambda functions compute common trigonometric operations where the
## input is in degrees.
sind    = lambda x: np.sin(np.radians(x))
cosd    = lambda x: np.cos(np.radians(x))
tand    = lambda x: np.tan(np.radians(x))
asind   = lambda x: np.degrees(np.arcsin(x))
acosd   = lambda x: np.degrees(np.arccos(x))
atand   = lambda x: np.degrees(np.arctan(x))
atan2d  = lambda x, y: np.degrees(np.arctan2(y, x))


def inherit_docstrings(cls=None, *, parent=None, do_specials=False):
    """
    @brief Propagates docstrings throughout the class inheritance hierarchy.

    This function walks through the class tree and replaces the placeholder 
    "*INHERIT*" found in the docstrings of class members with the corresponding
    docstring from the parent class.
    If no parent is explicitly provided, the function automatically ascends to 
    the top of the inheritance tree and then propagates docstrings down through
    all subclasses.

    <b>Usage:</b>
      - As a class decorator:
          @inherit_docstrings
      - Or by calling on the base class once its subclasses have been defined.
      - Can be invoked multiple times on the same class hierarchy if needed.

    @param cls
        The class whose docstrings are to be processed. When None, the function
        returns a decorator that processes the class upon definition.
    @param parent
        The parent class from which to inherit docstrings. When set to None, 
        the function traverses up to the top of the inheritance tree and then
        processes the entire subclass tree.
    @param do_specials
        A boolean flag indicating whether to process special attributes (those
        whose names both start and end with double underscores, e.g., 
        __dict__). Default value is False.

    @return
        The original class with updated docstrings.
    """
    if cls is None: # called w/ named arguments
        return lambda cls: inherit_docstrings(cls, parent=parent,
                                              do_specials=do_specials)
    if (parent is None) and cls.__bases__ and (cls.__bases__ != (object,)):
        # go up to the top of the tree and start there
        for base in cls.__bases__:
            inherit_docstrings(base, do_specials=do_specials)
        return cls

    if parent:
        # do this level in tree (cls itself)
        for name in dir(cls):
            if name in ['__doc__','__module__']: continue 
            if (not do_specials) and name.startswith('__') and \
               name.endswith('__'): 
                continue # don't mess with special items
            if not hasattr(parent, name): continue
            p = getattr(parent, name)
            if p is None: continue
            m = getattr(cls, name)
            if (p.__doc__ is not None) and hasattr(m, '__doc__') and \
               (m.__doc__ is not None):
                if type(m).__name__ == 'method':
                    # set the method's function's docstring
                    if not hasattr(p, '__func__'): continue
                    p = p.__func__  # from parent function's docstring
                    if p is None: continue
                    doc_owner = m.__func__
                else:
                    # set the object's docstring
                    doc_owner = m
                if '*INHERIT*' in doc_owner.__doc__:
                    # check avoids some read-only issues
                    doc_owner.__doc__ = doc_owner.__doc__.replace('*INHERIT*',
                                                                  p.__doc__)

    # do subclasses of cls, down the tree
    for c in cls.__subclasses__():
        inherit_docstrings(c, parent=cls, do_specials=do_specials)
        
    return cls


def get_list_neighbors(lst, q):
    """
    @brief Find the bounding neighbors of query value(s) in a sorted list.
    
    This function determines the indices of the two nearest values in a sorted
    list 'lst' that satisfy the condition: lst[index0] <= q < lst[index1]. For
    a scalar query 'q', it returns an array of shape (2,). For multiple query
    values, it returns an array of shape (len(q), 2). An index value of -1 
    indicates that the corresponding query is out of bounds.
    
    @param lst A list or array-like sequence of numeric values sorted in 
               ascending order.
    @param q A numeric scalar or an array-like sequence of query values.
    
    @return 
      - For a scalar q: a NumPy array of shape (2,), where the first element is
        the index such that lst[index0] <= q < lst[index1].
      - For an array of query values: a NumPy array of shape (len(q), 2) with
        the same semantics.
      - An index of -1 indicates that the query value q is not within the 
        bounds of lst.
    """
    if lst != sorted(lst):
        raise ValueError("Input list is not sorted.")
    lst = np.array(lst).ravel()
    N = len(lst)
    is_scalar = np.isscalar(q)
    q = np.atleast_1d(q).ravel()
    i = np.full((len(q), 2), -1)
    iq = np.isfinite(q)
    i[iq, 0] = interp1d(lst, np.arange(N, dtype=int), 'nearest',
                        bounds_error=False)(q[iq])
    f = (i[:, 0] >= 0)
    if any(f):
        # force to greatest list <= q
        i[f, 0] = i[f, 0] - (q[f] < lst[i[f, 0]])

        # next neighbor to right
        i[f, 1] = i[f, 0] + 1  

    # limit checks
    i[i >= N] = -1
    i[i < 0] = -1
    if is_scalar:
        i = i[0, :]
        
    return i


def make_deltas(grid, int_method='trapz', **kwargs):
    """
    @brief Compute integration weights for a 1D grid.

    Calculates default weights (deltas) for a given 1-dimensional grid based on
    the specified integration method. For grids with more than one element, the
    trapezoidal method is applied: 
      - the endpoints are weighted as half the distance to the adjacent point
      - the interior points are weighted as half the difference between the 
        following and preceding grid values.
    
    Special cases:
      - If grid is None or empty, returns 1.
      - If grid has a single element, returns that element (assuming it 
        represents a delta value).
    Note that additional keyword arguments are ignored.

    @param grid
        A 1D array-like object representing the grid. If None or empty, the 
        function returns 1.
    @param int_method
        Integration method to use. Currently, only 'trapz' (trapezoidal 
        integration) is supported.
    @param kwargs
        Additional keyword arguments (ignored).

    @return
        A NumPy array of weights with the same shape as grid for multi-element
        grids, a scalar value when grid has one element, or 1 for a None/empty
        grid.

    @exception ValueError
        Thrown when an unknown integration method is specified.
    """
    if grid is None:
        return 1
    grid = np.atleast_1d(grid).ravel()
    if len(grid) == 0:
        return 1
    elif len(grid) == 1:
        return grid[0]  # grid is actually dt

    d = np.full(grid.shape, np.nan)
    if int_method == 'trapz' or int_method == 'trapezoid':
        d[0] = (grid[1] - grid[0]) / 2
        d[1:-1] = (grid[2:] - grid[:-2]) / 2
        d[-1] = (grid[-1] - grid[-2]) / 2
    else:
        raise ValueError('Unknown int_method "%s"' % int_method)

    return d


def broadcast_grids(*args, mesh=False):
    """
    @brief Prepare multiple grids for broadcasting.
    
    This function reshapes two or more input grids so that they are mutually
    broadcastable. Two modes are provided:
    
      - When mesh is False, the input grids are reshaped by inserting singleton
        dimensions before and after the data dimension. For example, in a 2-D 
        case:
          - The first grid becomes shape (N1, 1)
          - The second grid becomes shape (1, N2)
        In a 3-D case:
          - The first grid becomes shape (N1, 1, 1)
          - The second grid becomes shape (1, N2, 1)
          - The third grid becomes shape (1, 1, N3)
          
      - When mesh is True, the function returns meshgrid arrays corresponding
        to the provided inputs, with each output array having the shape 
        (len(x1), len(x2), ..., len(xN)).
    
    @param args
        Two or more array-like grids that are to be broadcast.
    @param mesh
        Boolean flag specifying the output format:
          - If True, output arrays are created using numpy.meshgrid with 
            indexing set to "ij".
          - If False, each grid is reshaped so that the grids are mutually 
            broadcastable.
    
    @return A tuple of arrays with the appropriate reshaped or meshed grids.
    """
    if mesh:
        #return np.meshgrid([x.ravel() for x in args], indexing='ij')
        return np.meshgrid(*[np.ravel(x) for x in args], indexing='ij')
    else:
        NX = len(args)
        X = [None] * NX
        for i, x in enumerate(args):
            s = [1] * NX
            s[i] = np.array(x).size
            X[i] = np.reshape(x, s)
        return tuple(X)

    
def validate_grid(grid):
    """
    @brief Validate a 1-D grid.
    
    Checks if the provided grid is a one-dimensional NumPy array and that the
    values are strictly increasing.
    
    @param grid A NumPy array representing the grid.
    
    @return True if grid is a flat (1-D) array and strictly increasing; False
    otherwise.
    """
    if grid.shape != (grid.size,):
        return False
    if any(np.diff(grid) < 0):
        return False
    return True


def default_thetagrid(**kwargs):
    """
    @brief Generate a default theta grid.
    
    Returns a default theta grid as a 1-D NumPy array linearly spaced between 
    0.0 and 180.0 degrees with one degree spacing (closed endpoints). Keyword
    arguments are accepted for future modifications.
    
    @param kwargs Additional keyword arguments (currently ignored).
    
    @return A NumPy array representing the theta grid.
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0, 180.0, 181)  # 1 degree spacing, closed endpoints


def default_phigrid(**kwargs):
    """
    @brief Generate a default phi grid.
    
    Returns a default phi grid as a 1-D NumPy array linearly spaced between 0.0
    and 360.0 degrees with one degree spacing (closed endpoints). Keyword
    arguments are accepted for future modifications.
    
    @param kwargs Additional keyword arguments (currently ignored).
    
    @return A NumPy array representing the phi grid.
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0, 360.0, 361)  # 1 degree spacing, closed endpoints

def default_alphagrid(**kwargs):
    """
    @brief Generate a default alpha grid.
    
    Returns a default alpha grid as a 1-D NumPy array linearly spaced between
    0.0 and 180.0 degrees with one degree spacing (closed endpoints). Keyword
    arguments are accepted for future modifications.
    
    @param kwargs Additional keyword arguments (currently ignored).
    
    @return A NumPy array representing the alpha grid.
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0, 180.0, 181)  # 1 degree spacing, closed endpoints

def default_betagrid(**kwargs):
    """
    @brief Generate a default beta grid.
    
    Returns a default beta grid as a 1-D NumPy array linearly spaced between
    0.0 and 360.0 degrees with one degree spacing (closed endpoints). Keyword
    arguments are accepted for future modifications.
    
    @param kwargs Additional keyword arguments (currently ignored).
    
    @return A NumPy array representing the beta grid.
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0, 360.0, 361)  # 1 degree spacing, closed endpoints

def alphabeta2thetaphi(alpha, beta, alpha0, beta0, phib):
    """
    @brief Convert pitch angle and gyrophase to instrument angles.

    Given the pitch angle (alpha) and gyrophase (beta) along with instrument 
    orientation parameters (alpha0, beta0, phib), this function computes the 
    corresponding instrument angles theta and phi (all in degrees). When theta
    equals 0 or 180 degrees, phi is set to 0.

    The input angles alpha and beta must be mutually broadcastable, and 
    similarly, the parameters alpha0, beta0, and phib must be mutually 
    broadcastable. The resulting output shape is determined by broadcasting 
    alpha and beta with the broadcast shape of alpha0, beta0, and phib.

    @param alpha
      Pitch angle(s) in degrees.
    @param beta
      Gyrophase angle(s) in degrees.
    @param alpha0
      Orientation parameter alpha0 in degrees.
    @param beta0
      Orientation parameter beta0 in degrees.
    @param phib
      Orientation parameter phib in degrees.

    @return A tuple (theta, phi) where:
      - theta: Instrument angle theta in degrees.
      - phi: Instrument angle phi in degrees.
    The shapes of theta and phi correspond to the combined broadcasted shapes 
    of the input parameters.

    @exception ValueError
      Thrown if alpha and beta are not mutually broadcastable or if alpha0, 
      beta0, and phib are not mutually broadcastable.
    """
    
    try:
        _ = np.broadcast(alpha, beta)
    except ValueError:
        raise ValueError('alpha, and beta must be mutually broadcastable')
    try:
        _ = np.broadcast(alpha0, beta0, phib)
    except ValueError:
        raise ValueError("alpha0, beta0, and phib must be mutually" + \
                         "broadcastable")

    # define rotation matrix
    # columns are coefficients of c,d,b terms of instrument basis vectors
    # rows are coefficients of s1,s2,s0 terms of magnetic basis vectors
    # note: python inner lists are columns; in MATLAB a row of code is a row
    #       of the matrix, so this code is transposed compared to MATLAB.
    R = np.array([
        [-sind(beta0)*sind(phib) -cosd(alpha0)*cosd(beta0)*cosd(phib),  cosd(beta0) *sind(phib)-cosd(alpha0)*cosd(phib) *sind(beta0), cosd(phib)  *sind(alpha0)],
        [ cosd(phib) *sind(beta0)-cosd(alpha0)*cosd(beta0)*sind(phib), -cosd(beta0) *cosd(phib)-cosd(alpha0)*sind(beta0)*sind(phib),  sind(alpha0)*sind(phib)  ],
        [ cosd(beta0)*sind(alpha0),                                     sind(alpha0)*sind(beta0),                                     cosd(alpha0)             ]
    ])

    x = np.array([sind(alpha)*cosd(beta), sind(alpha)*sind(beta), cosd(alpha)])
    xpad = x.ndim - 1 # number of dimensions in x after the first
                      # dimension, which is 3
    # y = R*x (sum over the 2nd dimension of R and 1st dimension of x)
    
    # Order chosen to maintain alpha,beta shape before alpha0,beta0,phib shape
    y = np.tensordot(x, R, axes=((0,), (1,))) # do it in this order to get
                                              # alpha, beta shape before
                                              # alpha0, beta0, phib shape
    # y is in sensor Cartesian coordinates, but in the dimension after xpad

    ## Old version:
    #theta = acosd(y[3])
    #phi = atan2d(y[2], y[1])
    #phi[(theta == 0.0) | (theta == 180.0)] = 0

    ## New version:
    z = np.take(y, 2, xpad)
    if np.any(z < -1) or np.any(z > 1):
        ## Correct for last-digit floating point errors creating
        ## numbers out of bounds for acosd.
        cond1 = z > 1
        cond2 = z < -1
        z[cond1] = np.round(z[cond1], decimals=8)
        z[cond2] = np.round(z[cond2], decimals=8)
    theta = acosd(z) # third Cartesian coordinate (z)
    phi = atan2d(np.take(y, 1, xpad), np.take(y, 0, xpad)) # y, x as above
    phi[(theta==0.0)|(theta==180.0)] = 0
    if phi < 0:
        phi = np.mod(phi, 360)

    return theta, phi


def thetaphi2alphabeta(theta, phi, alpha0, beta0, phib):
    """
    @brief Convert instrument angles to pitch angle and gyrophase.

    This function performs the inverse transformation of instrument angles 
    (theta, phi) into pitch angle (alpha) and gyrophase (beta), all given in 
    degrees. When alpha equals 0 or 180 degrees, beta is set to 0.

    The input angles theta and phi must be mutually broadcastable, and the 
    parameters alpha0, beta0, and phib must be mutually broadcastable. The 
    computation leverages the similarity of the forward and inverse 
    transformations by swapping beta0 and phib.

    The alpha, beta to theta, phi transform is identical when alpha corresponds
    to theta, beta corresponds to phi, and beta0 corresponds to phib (makes
    R correspond to R'). Thus, we can just call alphabeta2thetaphi with phib
    and beta0 swapped.

    @param theta
      Instrument angle theta in degrees.
    @param phi
      Instrument angle phi in degrees.
    @param alpha0
      Orientation parameter alpha0 in degrees.
    @param beta0
      Orientation parameter beta0 in degrees.
    @param phib
      Orientation parameter phib in degrees.

    @return A tuple (alpha, beta) where:
      - alpha: Pitch angle in degrees.
      - beta: Gyrophase in degrees.
    The output shapes correspond to the broadcasted shapes of theta, phi and 
    alpha0, beta0, phib.

    @note The transformation exploits the equivalence between the forward and 
    inverse mappings by swapping beta0 and phib.
    """
    # phib <-> beta0
    return alphabeta2thetaphi(theta, phi, alpha0, phib, beta0)  

def vectors_to_euler_angles(B, C, S0, S1):
    """
    @brief Compute Euler angles from input vectors.

    Given four 3-vectors in a common Cartesian coordinate system (e.g., 
    Cartesian GEI), this function computes the Euler angles required to convert
    between instrument coordinates and magnetic coordinates. Specifically, it
    determines:
      - @b alpha0: the pitch angle of the boresight,
      - @b beta0: the gyrophase angle of the boresight,
      - @b phib: the longitude of the magnetic field vector B.

    The input vectors follow these conventions:
      - @b B is parallel to the magnetic field.
      - @b C defines the B-C plane in which beta = 0.
      - @b S0 points into the instrument (parallel to normally incident
           particles).
      - @b S1 defines the S0-S1 plane in which phi = 0.

    The vectors need not be normalized, and the planes (B-C and S0-S1) are not
    required to be orthogonal. Inputs may be provided either as a single 
    3-element vector ((3,)) or as an array of vectors with shape (N, 3). The
    output angles will be scalars for (3,) inputs or 1-D arrays of shape (N,)
    for (N,3) inputs. All angles are returned in degrees.

    @param B A 3-vector or an (N,3) array representing the magnetic field 
             direction.
    @param C A 3-vector or an (N,3) array that, together with B, defines the
             plane in which beta = 0.
    @param S0 A 3-vector or an (N,3) array representing the instrument pointing
              direction parallel to incident particles and opposite the 
              boresight direction.
    @param S1 A 3-vector or an (N,3) array that defines the S0-S1 plane in 
              which phi = 0.

    @return A tuple (alpha0, beta0, phib) where:
            - @b alpha0: Pitch angle of the boresight (in degrees).
            - @b beta0: Gyrophase angle of the boresight (in degrees).
            - @b phib: Longitude of B (in degrees).
            The output is scalar if the inputs are 3-element vectors, or a 1-D
            array of shape (N,) if the inputs are (N,3) arrays.

    @exception ValueError Thrown if the input vectors do not have a shape of 
    (3,) or (N,3) with a consistent N.
    """
    # Verify shape of input vectors.
    N = None
    isscalar = True
    for x in [B, C, S0, S1]:
        s = np.shape(x)
        if s == (3,): # fixed bug - old version had != instead of ==
            # Single vector with shape (3,).
            isscalar = False
            N = 1
        if len(s) == 2:
            # Series of vectors with shape (N, 3).
            N = s[0]
    if N is None:
        raise ValueError('B, C, S0, S1 inputs must be (3,), or (N, 3)' + \
                         'with common N')
    B = np.atleast_2d(B)
    C = np.atleast_2d(C)
    S0 = np.atleast_2d(S0)
    S1 = np.atleast_2d(S1)
    for x in [B, C, S0, S1]:
        if x.shape[0] not in (1, N):
            raise ValueError('B,C,S0,S1 inputs must be (3,), or (N,3)' + \
                             'with common N')
            
    hat = lambda x: x / np.sqrt((x**2).sum(axis=1, keepdims=True))
    dot = lambda x, y: (x * y).sum(axis=1)
    bhat = hat(B)
    dhat = hat(np.cross(B, C))
    chat = np.cross(dhat, bhat)
    S0hat = hat(S0)
    S2hat = hat(np.cross(S0, S1))
    S1hat = np.cross(S2hat, S0hat)
    cosa0 = max(-1.0, min(1.0, dot(bhat, S0hat)))  # bound to -1,1
    alpha0 = acosd(cosa0)
    
    phib = np.full((N,), np.nan)  # default phib is NaN
    i0 = alpha0 == 0
    if any(i0):
        phib[i0] = atan2d(dot(S1hat[i0], dhat[i0]), -dot(S2hat[i0], dhat[i0]))
    
    beta0 = np.zeros(N)  # beta=0 when alpha0=0
    beta0[np.logical_not(np.isfinite(alpha0))] = 0  # beta0=NaN when alpha0=NaN

    i0 = alpha0 > 0
    if any(i0):
        beta0[i0] = atan2d(dot(S0hat[i0], dhat[i0]), -dot(S0hat[i0], chat[i0]))
        phib[i0] = atan2d(dot(S2hat[i0], bhat[i0]), -dot(S1hat[i0], bhat[i0]))

    if isscalar:
        alpha0 = np.squeeze(alpha0, axis=0)
        beta0 = np.squeeze(beta0, axis=0)
        phib = np.squeeze(phib, axis=0)
    
    return alpha0, beta0, phib
       

# classes



class RFLError(Exception):
    """
    @brief Generic RFL error.

    This exception is raised for general errors in RFL. It accepts an optional
    message argument that describes the error.

    @param message Optional error message. Defaults to "Unknown RFL Error".
    """
    def __init__(self, message='Unknown RFL Error'):
        """
        @brief Construct a new RFLError instance.

        @param message Optional error message. Defaults to "Unknown RFL Error".
        """
        super().__init__(message)


class ArgSizeError(ValueError):
    """
    @brief Incorrect/Inconsistent Argument Size.

    This exception is raised when the provided arguments have incorrect or 
    inconsistent sizes. It accepts an optional message argument describing the
    error.

    @param message Optional error message. Defaults to "Incorrect/Inconsistent
    Argument Size".
    """
    def __init__(self, message='Incorrect/Inconsistent Argument Size'):
        """
        @brief Construct a new ArgSizeError instance.

        @param message Optional error message. Defaults to 
        "Incorrect/Inconsistent Argument Size".
        """
        super().__init__(message)


class KeywordError(RFLError):
    """
    @brief Missing or Invalid Keyword.

    This exception is raised when a required keyword is missing or invalid.
    It accepts an optional message argument describing the error.

    @param message Optional error message. Defaults to 
    "Missing/Invalid Keyword".
    """
    def __init__(self, message='Missing/Invalid Keyword'):
        """
        @brief Construct a new KeywordError instance.

        @param message Optional error message. Defaults to 
        "Missing/Invalid Keyword".
        """
        super().__init__(message)


# RFL helper functions
def keyword_check(*keywords, **kwargs):
    """
    @brief Check for required string keywords in a set of keyword arguments.

    This function verifies that each specified keyword exists in the provided
    kwargs dictionary and that its associated value is a string. If any keyword
    is missing or its value is not a string, a KeywordError is raised.

    @param keywords One or more keyword names to verify.
    @param kwargs A dictionary of keyword arguments.
    
    @return True if all specified keywords are present and are strings.
    
    @exception KeywordError If any keyword is missing or its value is not a 
    string.
    """
    for key in keywords:
        if key not in kwargs:
            raise KeywordError('Keyword %s required, absent' % key)
        if not isinstance(kwargs[key], str):
            raise KeywordError('%s must be a string(str)' % key)
    return True


def keyword_check_bool(*keywords, **kwargs):
    """
    @brief Check for required string keywords, returning a boolean result.

    This function checks whether each specified keyword exists in the kwargs
    dictionary and that its associated value is a string. Instead of raising
    an exception when a condition fails, it returns False.

    @param keywords One or more keyword names to verify.
    @param kwargs A dictionary of keyword arguments.
    
    @return True if all specified keywords are present and are strings; False
    otherwise.
    """
    for key in keywords:
        if key not in kwargs:
            return False
        if not isinstance(kwargs[key], str):
            return False
    return True


def keyword_check_numeric(*keywords, **kwargs):
    """
    @brief Check for required numeric keywords in a set of keyword arguments.

    This function verifies that each specified keyword exists in the kwargs 
    dictionary and that its associated value, or the first element if the value
    is a list or array, is numeric. If any keyword is missing or its value is
    non-numeric, a KeywordError is raised.

    @param keywords One or more keyword names to verify.
    @param kwargs A dictionary of keyword arguments.
    
    @return True if all specified keywords are present and numeric.
    
    @exception KeywordError If any keyword is missing or its value is not 
    numeric.
    """
    for key in keywords:
        if key not in kwargs:
            raise KeywordError('Keyword %s required, absent' % key)
        val = kwargs[key]
        while isinstance(val, (list, np.ndarray)):
            val = val[0]
        if not np.issubdtype(type(val), np.number):
            raise KeywordError('%s must be numeric' % key)
    return True


def keyword_check_numeric_bool(*keywords, **kwargs):
    """
    @brief Check for required numeric keywords, returning a boolean result.

    This function checks whether each specified keyword exists in the kwargs
    dictionary and verifies that its associated value (or the first element if
    the value is a list or array) is numeric. If a keyword is missing or its
    value is not numeric, the function returns False.

    @param keywords One or more keyword names to verify.
    @param kwargs A dictionary of keyword arguments.
    
    @return True if all specified keywords are present and numeric; False
    otherwise.
    """
    for key in keywords:
        if key not in kwargs:
            return False
        val = kwargs[key]
        while isinstance(val, (list, np.ndarray)):
            val = val[0]
        if not np.issubdtype(type(val), np.number):
            return False
    return True


def squeeze(val):
    """
    @brief Remove trailing singleton dimensions from a value.

    This function applies NumPy's squeeze to the input value to remove trailing
    singleton dimensions. If the resulting array contains a single element, it
    is converted to a scalar.

    @param val A NumPy array or value to be squeezed.
    
    @return The squeezed value, or a scalar if the resulting array has only one
    element.
    """
    val = np.squeeze(val)
    if val.size == 1:
        val = val.item()  # Reduce singleton to scalar (replaces np.asscalar)
        
    return val


def interp_weights_1d(xgrid, xhat, extrap_left=False, extrap_right=False,
                      period=None):
    """
    @brief Compute one-dimensional interpolation weights.

    This function computes interpolation weights for mapping values from a 
    one-dimensional grid onto specified query points. Linear interpolation is
    used between grid points, and several options are available to control the
    behavior for query points outside the grid range as well as periodic 
    boundary conditions.

    @param xgrid A list or array of grid values with shape (Nx,). The grid 
                 values should be in ascending order.
    @param xhat A scalar or list/array of query values where interpolation is
                desired.
    @param extrap_left Controls extrapolation for xhat values left of xgrid[0]:
           - False or 0: Assign a zero weight for xhat < xgrid[0].
           - True: Linearly extrapolate for xhat < xgrid[0].
           - "fixed": For xhat < xgrid[0], use the value at xgrid[0] (i.e., 
                      assign full weight to xgrid[0]).
    @param extrap_right Controls extrapolation for xhat values right 
                        of xgrid[-1]:
           - False or 0: Assign a zero weight for xhat > xgrid[-1].
           - True: Linearly extrapolate for xhat > xgrid[-1].
           - "fixed": For xhat > xgrid[-1], use the value at xgrid[-1] (i.e.,
                      assign full weight to xgrid[-1]).
    @param period Optional scalar that specifies the period over which x values
                  wrap around.
                  If provided, query points are reduced modulo period.

    @return A NumPy array of interpolation weights from the grid to the query 
            points.
            - If xhat is an array with shape (N,), the output has 
              shape (N, Nx).
            - If xhat is a scalar, the output is a one-dimensional array with
              shape (Nx,).
    """
    xgrid = np.array(xgrid).ravel()
    isscalar = np.isscalar(xhat)
    xhat = np.atleast_1d(xhat).ravel()
    if period:
        xhat = np.remainder(xhat, period)
    Nx = xgrid.size
    if Nx < 2:
        raise ValueError('xgrid must have at least 2 entries')
    N = xhat.size
    v = np.zeros((N, Nx))
    I = get_list_neighbors(xgrid, xhat)  # (N,2)
    # Prepare to broadcast.
    
    # Calculate weights for grid intervals via linear interpolation.
    for (k, (ileft, iright)) in enumerate(I):
        if (ileft >= 0) & (ileft < Nx - 1):  # xgrid[i] <= x < xgrid[i+1]
            v[k, ileft] = (xgrid[ileft + 1] - xhat[k]) / \
                (xgrid[ileft + 1] - xgrid[ileft])
        if (iright > 0) & (iright < Nx):  # xgrid[i-1] <= x < xgrid[i]
            v[k, iright] = (xhat[k] - xgrid[iright - 1]) / \
                (xgrid[iright] - xgrid[iright - 1])

    if extrap_left == 'fixed':
        f = xhat < xgrid[0]
        v[f, 0] = 1.0
    elif extrap_left:
        dx = (xgrid[1] - xgrid[0])
        f = (xhat < xgrid[0]) & (xhat > xgrid[0] - dx)
        # Linear extrapolation to the left
        v[f, 0] = 1.0 - (xgrid[0] - xhat[f]) / dx

    if period:
        dx = xgrid[0] + period - xgrid[-1]
        f = xhat < xgrid[0]
        v[f, -1] = (xgrid[0] - xhat[f]) / dx
        v[f, 0] = 1 - v[f, -1]
        f = xhat > xgrid[-1]
        v[f, 0] = (xhat[f] - xgrid[-1]) / dx
        v[f, -1] = 1 - v[f, 0]
    elif extrap_right == 'fixed':
        f = xhat > xgrid[-1]
        v[f, -1] = 1.0
    elif extrap_right:
        dx = xgrid[-1] - xgrid[-2]
        f = (xhat >= xgrid[-1]) & (xhat < xgrid[-1] + dx)
        v[f, -1] = 1.0 - (xhat[f] - xgrid[-1]) / dx
    if isscalar:
        v = v[0, :]  # Convert from (1, Nx) to (Nx,)
        
    return v


def interp_weights_3d(xgrid, xhat, ygrid, yhat, zgrid=None, zhat=None,
                      xopts={}, yopts={}, zopts={}):
    """
    @brief Compute three-dimensional interpolation weights.

    This function calculates interpolation weights for mapping values from a
    three-dimensional grid (constructed from separate 1-D grids for x, y, and
    optionally z) onto query points. The function internally calls 
    interp_weights_1d for each dimension and then combines the results via 
    outer multiplication. If the z-dimension grid is omitted (i.e., both zgrid
    and zhat are None), the function produces two-dimensional interpolation
    weights.

    @param xgrid A list or array of grid values along the x-axis with 
                 shape (Nx,).
    @param xhat A scalar or list/array of query values along the x-axis.
    @param ygrid A list or array of grid values along the y-axis with shape 
                 (Ny,).
    @param yhat A scalar or list/array of query values along the y-axis.
    @param zgrid Optional. A list or array of grid values along the z-axis with
                           shape (Nz,).
    @param zhat Optional. A scalar or list/array of query values along the 
                          z-axis.
    @param xopts Optional. A dictionary of options to be passed to 
                           interp_weights_1d for the x-axis.
    @param yopts Optional. A dictionary of options to be passed to 
                           interp_weights_1d for the y-axis.
    @param zopts Optional. A dictionary of options to be passed to
                           interp_weights_1d for the z-axis.

    @return A NumPy array of interpolation weights:
            - If zgrid and zhat are provided, the output shape is 
              (N, Nx, Ny, Nz), where N is the number of query points after 
              broadcasting.
            - If zgrid and zhat are None, the output has shape (N, Nx, Ny).
            - If the query points are given as scalars, appropriate dimensions
              are squeezed.
    
    @exception ValueError Raised if one of zgrid or zhat is provided without
    the other, or if the query points are not mutually broadcastable.
    """
    if (zgrid is None) != (zhat is None):
        raise ValueError('if either zgrid or zhat is None, both must be')
    
    isscalar = np.isscalar(xhat) and np.isscalar(yhat)
    xhat = np.atleast_1d(xhat).ravel()
    yhat = np.atleast_1d(yhat).ravel()
    N = max(len(xhat), len(yhat), len(zhat))
    if zhat is not None:
        isscalar = isscalar and np.isscalar(zhat)
        zhat = np.atleast_1d(zhat).ravel()
        N = max(N, len(zhat))
    
    if len(xhat) not in [1, N]:
        raise ValueError('xhat must be broadcastable to same size as' + \
                         'yhat, zhat')
    wx = interp_weights_1d(xgrid, xhat, **xopts)  # (len(xhat) x Nx)
    
    if len(yhat) not in [1, N]:
        raise ValueError('yhat must be broadcastable to same size as' + \
                         'xhat, zhat')
    wy = interp_weights_1d(ygrid, yhat, **yopts)  # (len(yhat) x Ny)
    
    if zhat is not None:
        if len(zhat) not in [1, N]:
            N = max(N, len(zhat))
            raise ValueError('zhat must be broadcastable to same size as' + \
                             'xhat, yhat')
        wz = interp_weights_1d(zgrid, zhat, **zopts)  # (len(zhat) x Nz)
    
    if zhat is not None:  # Make weights broadcastable to (N, Nx, Ny, Nz)
        wy = np.expand_dims(np.expand_dims(wy, 1), 3)
        wz = np.expand_dims(np.expand_dims(wz, 1), 1)
    else:  # Make weights broadcastable to (N, Nx, Ny)
        wy = np.expand_dims(wy, 1)
        wz = 1  # Dummy multiplier for 2D interpolation
    
    w = wx * wy * wz
    if isscalar:  # Squeeze the first dimension if query points are scalars
        w = np.squeeze(w, axis=0)
        
    return w


def interp_weights_2d(xgrid, xhat, ygrid, yhat, xopts={}, yopts={}):
    """
    @brief Compute two-dimensional interpolation weights.

    This function calculates interpolation weights for mapping values from a 
    two-dimensional grid, defined by separate 1-D grids for x and y, onto query
    points. It leverages interp_weights_3d by setting the z-dimension inputs to
    None, thereby generating 2-D interpolation weights.

    @param xgrid A list or array of grid values along the x-axis with shape
                 (Nx,).
    @param xhat A scalar or list/array of query values along the x-axis.
    @param ygrid A list or array of grid values along the y-axis with
                 shape (Ny,).
    @param yhat A scalar or list/array of query values along the y-axis.
    @param xopts Optional. A dictionary of options to be passed to 
                           interp_weights_1d for the x-axis.
    @param yopts Optional. A dictionary of options to be passed to 
                           interp_weights_1d for the y-axis.

    @return A NumPy array of interpolation weights from the grid to the query 
            points.
            - If xhat and yhat are arrays with shape (N,), the output has 
              shape (N, Nx, Ny).
            - If xhat and yhat are scalars, the output is a two-dimensional 
              array with shape (Nx, Ny).
    """
    return interp_weights_3d(xgrid, xhat, ygrid, yhat, xopts, yopts)


class FactoryConstructorMixin(object):
    """
    @brief A mixin to support factory-style object construction.

    Enables an object constructor to traverse its subclass tree in order to 
    instantiate the correct subclass based on the initialization arguments.

    Classes that inherit from FactoryConstructorMixin should implement the 
    classmethod is_mine() in every subclass to “claim” the provided 
    initialization data and thereby instantiate that subclass rather than the
    parent.

    @classmethod 
    def is_mine(cls,*args,**kwargs):
        # search args, kwargs dict/tree to determine if this subclass should 
        # handle the data, and if so, return True     

    Any class that inherits directly from FactoryConstructorMixin should 
    implement is_mine() so that it raises an exception instead of returning
    False. 

    A subclass that does not override is_mine will simply replace its
    parent in the factory hierarchy. This approach will modify the behavior of
    a subclass. A subclass that defines its own is_mine will extend its parent
    into a distinct new subclass in the hierarchy. This approach is appropriate
    when the parent and child are alternatives that should coexist.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """
        bool = is_mine(cls, *args, **kwargs)
        @brief Return true if the provided arguments describe an object of 
               this class.
        """
        raise Exception("Class '%s' failed to" % cls.__name__ + \
                        " implement is_mine classmethod or didn't claim case")

    def __new__(cls, *args, **kwargs):
        """
        @brief Factory method for creating an instance of the correct subclass.

        Traverses the subclass tree and returns the first matching subclass 
        instance based on the provided arguments.
        """
        # Check children first
        if args or kwargs:
            for c in cls.__subclasses__():
                res = c.__new__(c, *args, **kwargs)
                if res:
                    return res
        if cls.is_mine(*args, **kwargs):
            # Python will pass kwargs to init.
            return super().__new__(cls)
        return None

    def __init__(self, *args, **kwargs):
        """
        @brief Object initialization.

        Calls the parent initializer without passing any additional arguments.
        """
        # Avoid trying to pass arguments to object initializer.
        super().__init__()  


class EnergyResponse(FactoryConstructorMixin):
    """
    @brief Abstract class representing energy responses.

    Initialization accepts the keywords CROSSCALIB and EPS (which default to 1
    if absent).
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        keyword_check('E_TYPE', **kwargs)
        raise KeywordError('Supplied keywords do not define a ' + \
                           'recognized EnergyResponse')

    def __init__(self, **kwargs):
        """
        @brief Initialize energy response properties.

        Handles the keywords 'CROSSCALIB' and 'EPS', defaulting to 1 if not 
        provided.
        Also sets additional default values such as the energy unit.
        """
        for arg in ['CROSSCALIB', 'EPS']:
            if keyword_check_numeric_bool(arg, **kwargs):
                setattr(self, arg, squeeze(kwargs[arg]))
            else:
                print("%s not found or " % arg + \
                      "not numeric in channel description, assuming unity")
                setattr(self, arg, 1.0)
        for arg, val in {'E_UNIT': 'MeV'}.items():
            if arg in kwargs:
                setattr(self, arg, kwargs[arg])
            else:
                setattr(self, arg, val)

    def RE(self, E):
        """
        @brief Evaluate the channel response at energy E (without CROSSCALIB).

        @param E Energy value or array of energy values.
        @return The channel response (dimensionless).
        @exception NotImplementedError if not overloaded by a subclass.
        """
        raise NotImplementedError("Class %s " % self.__class__.__name__ + \
                                  "did not overload of RE method, as required")

    def hE(self, Egrid, **kwargs):
        """
        @brief Compute energy response weights on the given energy grid.

        @param Egrid A 1-D numpy array representing energies.
        @param kwargs Optional additional parameters.
        @return A 1-D numpy array of energy response weights with CROSSCALIB 
                applied.
        @exception NotImplementedError if not overloaded by a subclass.
        """
        raise NotImplementedError('Class %s '  % self.__class__.__name__ + \
                                  'did not overload of hE method, as required')

    def hE0(self, Egrid=None, **kwargs):
        """
        @brief Return the integrated energy response.

        @param Egrid Optional 1-D numpy array over which to integrate. If not
                     supplied, the stored _E0 value is returned.
        @return A scalar energy response (in MeV) after integration.
        @exception ValueError if Egrid is None and _E0 is not defined.
        """
        if Egrid is None:
            if self._E0 is None:
                raise ValueError('Egrid must be supplied ' + \
                                 'for %s.hE0' % self.__class__.__name__)
            return self._E0
        else:
            hE = self.hE(Egrid, **kwargs)
            return hE.sum() # integrate over energy


class ER_Diff(EnergyResponse):
    """
    @brief Class representing differential energy channel response.

    Initialization requires EPS, E0, and DE, with EPS being optional (defaults
    to 1).
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'DIFF')

    def __init__(self, **kwargs):
        """
        @brief Initialize a differential energy response.

        Processes numeric keywords 'E0' and 'DE'.
        """
        super().__init__(**kwargs) # handles CROSSCALIB and EPS
        keyword_check_numeric('E0', 'DE', **kwargs)
        self.E0 = squeeze(kwargs['E0'])  # TODO: convert to MeV
        self.DE = squeeze(kwargs['DE'])  # TODO: convert to MeV
        self._hE0 = self.DE * self.EPS / self.CROSSCALIB # hE for flat spectrum

    def RE(self, E):
        """*INHERIT*"""
        # DE gets ignored, which is bad, but unavoidable
        return (E == self.E0) * self.EPS 

    def hE(self, Egrid, **kwargs):
        """*INHERIT*"""
        hE = interp_weights_1d(Egrid, self.E0) * self.dE * self.EPS
        return hE / self.CROSSCALIB
        

class ER_Int(EnergyResponse):
    """
    @brief Class representing integral energy channel response.

    Initialization requires EPS and E0; provides an integrated energy response
    for a flat spectrum.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'INT')

    def __init__(self, **kwargs):
        """
        @brief Initialize the integral energy response.

        Processes the numeric keyword 'E0'.
        """
        super().__init__(**kwargs) # handles CROSSCALIB and EPS
        keyword_check_numeric('E0', **kwargs)            
        self.E0 = squeeze(kwargs['E0'])  # TODO: convert to MeV
        self._hE0 = None # not defined for integral channel

    def RE(self, E):
        """*INHERIT*"""
        return (E >= self.E0) * self.EPS

    def hE(self, Egrid, **kwargs):
        """*INHERIT*"""
        Egrid = np.array(Egrid)
        NE = len(Egrid)
        hE = np.zeros(Egrid.shape)
        I = get_list_neighbors(Egrid, self.E0)
        # Egrid[I[0]] <= inst_info.E0 < Egrid[I[1]]
        # or I[1] = -1 or both I = -1
        
        dE = make_deltas(Egrid, **kwargs)
        hE = dE # default, eventually only kept for Egrid > Egrid[I[1]]
        hE[0:(I[0] + 1)] = 0 # zero out below E0

        # left side
        i = I[0]
        if (i >= 0) and (i < NE - 1): # E[i] <= E0 < E[i+1]
            hE[i] = (Egrid[i+1] - self.E0) ** 2 / 2 / (Egrid[i+1] - Egrid[i])

        # right side
        i = I[1]
        if (i >= 1) and (i < NE): # E[i-1] < E0 < E[i]
            hE[i] = (self.E0 - Egrid[i-1]) ** 2 / 2 / (Egrid[i] - Egrid[i-1])
            
        return hE * self.EPS / self.CROSSCALIB
    

class ER_Wide(EnergyResponse):
    """
    @brief Class representing wide differential energy channel response.

    Initialization requires EPS, E0, and E1 (with EPS optional and defaulting
    to 1). The wide-channel response is modeled as the difference between two
    integral energy channels.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'WIDE')

    def __init__(self, **kwargs):
        """
        @brief Initialize the wide differential energy response.

        Processes numeric keywords 'E0' and 'E1' and ensures E1 > E0. 
        Constructs two integral channels for computing the difference.
        """
        super().__init__(**kwargs) # handles CROSSCALIB and EPS
        keyword_check_numeric('E0', 'E1', **kwargs)            
        self.E0 = squeeze(kwargs['E0'])  # TODO: convert to MeV
        self.E1 = squeeze(kwargs['E1'])  # TODO: convert to MeV
        if self.E1 <= self.E0:
            raise ValueError('E1 must be greater than E0 for wide channel')

        # For flat spectrum:
        self._hE0 = (self.E1 - self.E0) * self.EPS / self.CROSSCALIB
        
        # build two integral channels to difference them
        tmp = {**kwargs, 'E_TYPE': 'INT'} # integral channel with same E0
        self.low = ER_Int(**tmp)
        tmp['E0'] = self.E1 # integral channel above this one
        self.high = ER_Int(**tmp)

    def RE(self, E):
        """*INHERIT*"""
        return ((E >= self.E0) & (E <= self.E1)) * self.EPS

    def hE(self, Egrid, **kwargs):
        """*INHERIT*"""
        # treat as difference between two integral channels
        return self.low.hE(Egrid) - self.high.hE(Egrid)


class ER_Table(EnergyResponse):
    """
    @brief Class representing tabular energy channel response.

    This class represents a tabular energy channel response. Initialization 
    requires the keywords EPS and E_GRID, where EPS is optional (defaults to 
    1).
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'TBL')

    def __init__(self, **kwargs):
        """
        @brief Initialize the ER_Table instance.

        Processes the keyword 'E_GRID' to define a valid 1-D energy grid (in 
        MeV) and the keyword 'EPS'. If EPS is provided as a scalar, it is 
        expanded to match the shape of E_GRID; if provided as an array, its
        shape must match that of E_GRID. Also calculates the integrated energy
        response for a flat spectrum (_hE0) and initializes the cached 
        interpolation function (_RE).

        @param kwargs Keyword arguments including:
               - E_GRID: A 1-D numeric grid of energy values.
               - EPS: A numeric value or array representing the response 
                      (optional, defaults to 1).
        @return None
        @exception ValueError If E_GRID is not a valid grid (must be 1-D and
                              contain unique values).
        @exception ArgSizeError If the shapes of E_GRID and EPS do not match.
        """
        super().__init__(**kwargs)  # handles CROSSCALIB and EPS
        keyword_check_numeric('E_GRID', **kwargs)
        self.E_GRID = squeeze(kwargs['E_GRID'])  # TODO: convert to MeV
        if not validate_grid(self.E_GRID):
            raise ValueError('E_GRID is not a valid grid: 1-d, unique')

        keyword_check_numeric('EPS', **kwargs)
        self.EPS = squeeze(kwargs['EPS'])  # TODO: convert to MeV
        if np.isscalar(self.EPS):
            self.EPS = np.full(self.E_GRID.shape, self.EPS)
        elif self.EPS.shape != self.E_GRID.shape:
            raise ArgSizeError('E_GRID and EPS are not the same shape')

        # For flat spectrum:
        self._hE0 = np.sum(self.E_GRID * self.EPS * \
                           make_deltas(self.E_GRID)) / self.CROSSCALIB
        
        self._RE = None

    def RE(self, E):
        """*INHERIT*"""
        if self._RE is None:
            # extrapolate beyond grid with nearest edge (first/last) value
            self._RE = interp1d(self.E_GRID, self.EPS, 'linear',
                                bounds_error=False,
                                fill_value=(self.EPS[0], self.EPS[-1]),
                                assume_sorted=True)
        return self._RE(E)

    def hE(self, Egrid, **kwargs):
        """*INHERIT*"""
        dE = make_deltas(Egrid, **kwargs)
        hE = self.RE(Egrid) * dE / self.CROSSCALIB
        return hE


class AngleResponse(FactoryConstructorMixin):
    """
    @brief Class representing angular responses.

    This class models angular responses (effective areas) that are forward-only
    by default. Initialization accepts the keyword BIDIRECTIONAL to indicate 
    whether the response is bidirectional.
    Properties include:
      - bidirectional: a boolean that specifies if the angular response is 
                       bidirectional.
      - hA0: the nominal geometric factor (in cm^2 sr), accounting for 
             bidirectional effects.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""     
        if not kwargs:
            return True # null response
        else:
            raise KeywordError('The data provided did not define a ' + \
                               'recognized AngleResponse')

    def __init__(self, **kwargs):
        """
        @brief Initialize the AngleResponse instance.

        Determines if the response is bidirectional based on the BIDIRECTIONAL
        keyword. Also sets default values for additional properties, such as
        the length unit (L_UNIT).

        @param kwargs Keyword arguments that may include:
               - BIDIRECTIONAL: A boolean or string ('TRUE') indicating 
                                bidirectionality.
               - L_UNIT: A string representing the length unit (defaults 
                         to 'cm').
        @return None
        """
        super().__init__(**kwargs)
        self.bidirectional = ('BIDIRECTIONAL' in kwargs) and \
            (kwargs['BIDIRECTIONAL'] in [True, 'TRUE'])
        
        # Set others to default values.
        for arg, val in {'L_UNIT': 'cm'}.items():
            if arg in kwargs:
                setattr(self, arg, kwargs[arg])
            else:
                setattr(self, arg, val)  # default

    def A(self, theta, phi):
        """
        @brief Evaluate the angular response (effective area) at specified 
               theta and phi.

        Computes the effective area in cm^2 given polar (theta) and azimuthal 
        (phi) angles. The input angles must be broadcast-compatible, and the 
        returned value or array will have the same broadcast shape as theta and
        phi.

        @param theta A numeric value or array representing the polar angle in
                     degrees.
        @param phi A numeric value or array representing the azimuthal angle
                   in degrees.
        @return The effective area (in cm^2) as a scalar or numpy array.
        @exception NotImplementedError If the method is not overridden in a 
                                       subclass.
        """
        raise NotImplementedError('Class %s ' % self.__class__.__name__ + \
                                  'did not overload of A method, as required')

    def hAthetaphi(self, thetagrid, phigrid, **kwargs):
        """
        @brief Compute angular response weights on a theta-phi grid.

        Returns the angular response weights (in cm^2 sr) on a grid defined by
        thetagrid and phigrid. Differential elements are computed using the
        deltas of -cos(theta) and phi (converted to radians), and the response
        A is evaluated on the broadcasted grid.

        @param thetagrid A 1-D numpy array of theta values in degrees.
        @param phigrid A 1-D numpy array of phi values in degrees.
        @param kwargs Additional keyword arguments to be passed to make_deltas.
        @return A scalar or a 2-D numpy array with shape (len(thetagrid), 
                len(phigrid)) representing the response weights.
        """
        dcostheta = np.abs(make_deltas(-cosd(thetagrid), **kwargs))
        dphi = make_deltas(np.radians(phigrid), **kwargs)

        ## Old version - no mesh
        #thetaI, phiI = broadcast_grids(thetagrid, phigrid)
        #dcosthetaI, dphiI = broadcast_grids(dcostheta, dphi)

        ## New version - mesh=True
        thetaI, phiI = broadcast_grids(thetagrid, phigrid, mesh=True)
        dcosthetaI, dphiI = broadcast_grids(dcostheta, dphi, mesh=True)
        
        # Includes backward response if present.
        return self.A(thetaI, phiI) * dcosthetaI * dphiI

    def hAtheta(self, thetagrid, phigrid=None, **kwargs):
        """
        @brief Compute angular response weights on a theta grid integrated 
               over phi.

        If phigrid is not provided, it defaults to a full 0–360 degree range.
        The response weights (in cm^2 sr) are integrated over the phi 
        dimension.

        @param thetagrid A 1-D numpy array of theta values in degrees.
        @param phigrid Optional 1-D numpy array of phi values in degrees; 
                       defaults to full range if omitted.
        @param kwargs Additional keyword arguments.
        @return A scalar or a 1-D numpy array with length equal to 
                len(thetagrid) representing the integrated response weights.
        """
        if phigrid is None:
            phigrid = default_phigrid(**kwargs)
        h = self.hAthetaphi(thetagrid, phigrid, **kwargs)

        return h.sum(axis=1)  # integrate over phi

    def hAalphabeta(self, alpha0, beta0, phib,
                    alphagrid, betagrid, tgrid=None, **kwargs):
        """
        @brief Compute angular response weights on an alpha-beta grid.

        This function returns angular response weights based on an alpha-beta
        grid. It computes differential elements using the deltas of -cos(alpha)
        and beta (converted to radians), converts alpha and beta to theta and
        phi using alphabeta2thetaphi, and evaluates the response A. If a time
        grid (tgrid) is provided, the result is integrated over the time
        dimension. The output units are cm^2·sr or cm^2·sr·s if time 
        integration is performed.

        @param alpha0 Numeric value representing the parameter alpha0.
        @param beta0 Numeric value representing the parameter beta0.
        @param phib Numeric value representing the parameter phib.
        @param alphagrid A 1-D numpy array of alpha values in degrees.
        @param betagrid A 1-D numpy array of beta values in degrees.
        @param tgrid Optional time grid; can be a scalar or 1-D numpy array.
        @param kwargs Additional keyword arguments to be passed to make_deltas.
        @return A scalar or a 2-D numpy array with shape (len(alphagrid), 
                len(betagrid)) representing response weights.
        """
        dcosa = make_deltas(-cosd(alphagrid), **kwargs)  # (Na,)
        db = np.radians(make_deltas(betagrid, **kwargs))  # (Nb,)
        dt = make_deltas(tgrid, **kwargs)  # (Nt,) or scalar

        ## Old version - no mesh
        #dcosa, db = broadcast_grids(dcosa, db)  # (Na, Nb)

        ## New version - mesh=True
        dcosa, db = broadcast_grids(dcosa, db, mesh=True)  # (Na, Nb)
        
        theta, phi = alphabeta2thetaphi(alphagrid, betagrid, alpha0, beta0,
                                        phib)  # (Na, Nb, Nt) or (Na, Nb)
        A = self.A(E, theta, phi)  # (Na, Nb, *Nt)
        if np.isscalar(dt):
            h = A * dt
        else:
            # dot product along the Nt dimension
            h = np.tensordot(A, dt, axes=((2,), (0,)))
            
        # now h is (Na,Nb)
        h = h * dcosa * db / self.CROSSCALIB
        
        return h

    def hAalpha(self, alpha0, beta0, phib, alphagrid, betagrid=None,
                tgrid=None, **kwargs):
        """
        @brief Compute angular response weights on an alpha grid after 
               integrating over beta.

        This function computes the angular response weights based on an alpha
        grid by first adjoining a beta grid. If betagrid is not provided, it is
        assumed to cover the full range of 0–360 degrees. If a time grid is 
        supplied, the result is integrated over time. The output units are 
        cm^2·sr or cm^2·sr·s if time-integrated.

        @param alpha0 Numeric value representing the parameter alpha0.
        @param beta0 Numeric value representing the parameter beta0.
        @param phib Numeric value representing the parameter phib.
        @param alphagrid A 1-D numpy array of alpha values in degrees.
        @param betagrid Optional 1-D numpy array of beta values in degrees; 
                        defaults to full range if omitted.
        @param tgrid Optional time grid for integration.
        @param kwargs Additional keyword arguments.
        @return A scalar or a 1-D numpy array of response weights with length
                equal to len(alphagrid).
        """
        if betagrid is None:
            betagrid = default_betagrid(**kwargs)
        h = self.hAalphabeta(alpha0, beta0, phib,
                             alphagrid, betagrid, tgrid, **kwargs)
        return h.sum(axis=1)  # integrate over beta

    @property
    def hA0(self):
        """
        @brief Return the nominal geometric factor.

        The nominal geometric factor (in cm^2·sr) is determined from the 
        property G. If the response is bidirectional, the geometric factor is
        doubled to include both hemispheres.

        @return The nominal geometric factor in cm^2·sr.
        """
        hA0 = self.G
        if self.bidirectional:
            hA0 *= 2
        return hA0

class AR_csym(AngleResponse):
    """
    @brief Virtual base class representing phi-symmetric angular responses.

    This abstract class represents angular responses that are symmetric in phi
    (i.e., independent of phi) for forward response only.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        return False  # abstract class, never the right answer
                
    def hAtheta(self, thetagrid, phigrid=None, **kwargs):
        """*INHERIT*"""
        dcostheta = np.abs(make_deltas(-cosd(thetagrid), **kwargs))
        if phigrid is None:
            dphi = 2*np.pi
        else:
            dphi = make_deltas(np.radians(phigrid), **kwargs).sum()
        phigrid = 0.0  # dummy value to save calculation
        # includes backward response if present
        return self.A(thetagrid, phigrid) * dcostheta * dphi


class AR_Omni(AR_csym):
    """
    @brief Class representing omnidirectional angular responses.

    This class represents omnidirectional angular responses. Initialization 
    requires G (None OK for derived classes). An omni response is modeled as a
    sphere (i.e., it has the same apparent area from all directions). For flat,
    single-element detectors, use AR_Disk or AR_Slab.

    Properties:
      - area: Detector area in cm^2.
      - G: Geometric factor for forward hemisphere particles (theta <= 90) 
           in cm^2 sr.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        return keyword_check_bool('RESP_TYPE', **kwargs) and \
               (kwargs['RESP_TYPE'] == '[E]')
               
    def __init__(self, **kwargs):
        """
        @brief Initialize the omnidirectional angular response.

        Sets up the response using the provided geometric factor G. The
        detector area is computed as G/(2pi) for a half omni and is adjusted if
        the response is bidirectional.
        
        @param kwargs Expected to include:
               - G: Geometric factor (optional).
               - BIDIRECTIONAL: Indicates bidirectional response (optional).
        """
        super().__init__(**kwargs)  # handles BIDIRECTIONAL
        G = None
        if ('G' in kwargs) and kwargs['G'] is not None:
            # derived class may provide G as None
            keyword_check_numeric('G', **kwargs)
            # Setup for half omni: full G over 2pi
            G = squeeze(kwargs['G'])  # TODO: convert to cm^2 sr
            self.G = G
            self.area = self.G / 2 / np.pi  # equivalent sphere's cross-section
            if self.bidirectional:
                self.area = self.area / 2
                
    def A(self, theta, phi):
        """*INHERIT*"""
        if self.bidirectional:
            return self.area
        else:
            return self.area * (theta <= 90)


class AR_SingleElement(AR_Omni):
    """
    @brief Abstract class representing single-element detectors.

    This class represents single-element detectors (i.e., with cos projection
    effect). Derived classes must define self.area (the detector area) and
    self.G (the geometric factor).
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        return False  # abstract class, never the right answer

    def __init__(self, **kwargs):
        """
        @brief Initialize the single-element detector response.

        Supplies a fake G to the superclass; the derived class should set area
        and G appropriately.
        
        @param kwargs Keyword arguments (handles BIDIRECTIONAL, etc.).
        """
        # Supply fake G to super; area and G will be set in derived class.
        super().__init__(G=None, **kwargs)  # handles BIDIRECTIONAL

    def A(self, theta, phi):
        """*INHERIT*"""
        return super().A(theta, phi) * np.abs(cosd(theta))

    @property
    def G(self):
        """
        @brief Geometric factor for forward-going particles.

        Computed as pi multiplied by the detector area.
        
        @return Geometric factor in cm^2 sr.
        """
        return self.area * np.pi


class AR_Disk(AR_SingleElement):
    """
    @brief Class representing disk angular responses.

    This class models disk angular responses. Initialization requires:
      - R1: Disk radius (cm).
      - R2: (Unused in this calculation but provided for compatibility.)
      - D: (Unused in this calculation but provided for compatibility.)
    The detector area is computed as pi * R1^2.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return keyword_check_bool('TH_TYPE', **kwargs) and \
               (kwargs['TH_TYPE'] == 'DISK')
               
    def __init__(self, **kwargs):
        """
        @brief Initialize the disk angular response.

        Processes the numeric keyword 'R1' to determine the detector area.
        
        @param kwargs Expected to include:
               - R1: Disk radius.
        """
        super().__init__(**kwargs)  # handles BIDIRECTIONAL
        keyword_check_numeric('R1', **kwargs)
        self.R1 = squeeze(kwargs['R1'])  # TODO: convert to cm
        self.area = np.pi * self.R1**2


class AR_Slab(AR_SingleElement):
    """
    @brief Class representing rectangular slab angular responses.

    This class models angular responses for rectangular slab detectors. 
    Initialization requires:
      - W1: Width (cm).
      - H1: Height (cm).
    The detector area is computed as W1 * H1.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return keyword_check_bool('TH_TYPE', **kwargs) and \
               (kwargs['TH_TYPE'] == 'SLAB')
               
    def __init__(self, **kwargs):
        """
        @brief Initialize the slab angular response.

        Processes numeric keywords 'W1' and 'H1' to compute the detector area.
        
        @param kwargs Expected to include:
               - W1: Width of the slab.
               - H1: Height of the slab.
        """
        super().__init__(**kwargs)  # handles BIDIRECTIONAL
        keyword_check_numeric('W1', 'H1', **kwargs)
        self.W1 = squeeze(kwargs['W1'])  # TODO: convert to cm
        self.H1 = squeeze(kwargs['H1'])  # TODO: convert to cm
        self.area = self.W1 * self.H1


class AR_Tele_Cyl(AR_csym):
    """
    @brief Class representing cylindrically-symmetric telescope angular 
           responses.

    Using Sullivan's 1971 paper (updated ~2010), this class models the angular
    response of a cylindrical telescope. Initialization requires:
      - R1: First radius (cm).
      - R2: Second radius (cm).
      - D: Distance (cm).
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return keyword_check_bool('TH_TYPE', **kwargs) and \
               (kwargs['TH_TYPE'] == 'CYL_TELE')
               
    def __init__(self, **kwargs):
        """
        @brief Initialize the cylindrical telescope angular response.

        Processes numeric keywords 'R1', 'R2', and 'D'. Calculates:
          - Rs: Minimum of R1 and R2.
          - thetac: Critical angle computed from |R1 - R2|/D.
          - thetam: Maximum angle computed from (R1 + R2)/D.
          - area: Effective area at normal incidence.
          - G: Geometric factor (computed via Equation 8 from Sullivan's paper)
        If the response is bidirectional, G is doubled.
        
        @param kwargs Expected to include:
               - R1, R2, D: Numeric values (in cm).
        """
        super().__init__(**kwargs)  # handles BIDIRECTIONAL
        keyword_check_numeric('R1', 'R2', 'D', **kwargs)
        self.R1 = squeeze(kwargs['R1'])  # TODO: convert to cm
        self.R2 = squeeze(kwargs['R2'])  # TODO: convert to cm
        self.D = squeeze(kwargs['D'])    # TODO: convert to cm
        self.Rs = min(self.R1, self.R2)
        self.thetac = atand(np.abs(self.R1 - self.R2) / self.D)
        self.thetam = atand((self.R1 + self.R2) / self.D)
        self.area = np.pi * self.Rs**2  # area at normal incidence
        tmp = self.R1**2 + self.R2**2 + self.D**2
        self.G = np.pi**2 / 2 * (tmp - np.sqrt(tmp**2 - 4 * self.R1**2 \
                                               * self.R2**2))  # eqn 8
        if self.bidirectional:
            self.G *= 2

    def A(self, theta, phi):
        """*INHERIT*"""
        # Uses Sullivan's updated paper, equation 10.
        A = np.zeros(theta.shape)
        if self.bidirectional:
            # treat backward as forward
            theta = np.minimum(theta, 180.0 - theta)
        thetarads = np.radians(theta)
        f = theta <= self.thetac
        if any(f):
            A[f] = np.pi * self.Rs**2 * np.cos(thetarads[f])
        f = (theta > self.thetac) & (theta < self.thetam)
        if any(f):
            tantheta = self.tan(thetarads[f])
            Psi1 = np.arccos((self.R1**2 + self.D**2 * tantheta**2 - \
                              self.R2**2) / (2 * self.D * self.R1 * tantheta))
            Psi2 = np.arccos((self.R2**2 + self.D**2 * tantheta**2 - \
                              self.R1**2) / (2 * self.D * self.R2 * tantheta))
            # Psi1 and Psi2 are both in radians.
            
            A[f] = np.cos(thetarads[f]) / 2 * \
                (self.R1**2 * (2 * Psi1 - np.sin(2 * Psi1)) + \
                 self.R2**2 * (2 * Psi2 - np.sin(2 * Psi2)))
        return np.broadcast_to(A, np.broadcast(theta, phi).shape)

class AR_Table_sym(AngleResponse):
    """
    @brief Class representing phi-symmetric tabular angular responses.

    This class represents a phi-symmetric tabular angular response. 
    Initialization requires TH_GRID and A. It ignores bidirectional effects; 
    for a bidirectional tabular response, provide A defined on TH_GRID spanning
    0 to 180 degrees.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return keyword_check_bool('TH_TYPE', **kwargs) and \
            (kwargs['TH_TYPE'] == 'TBL')
            
    def __init__(self, **kwargs):
        """
        @brief Initialize the AR_Table_sym instance.

        Processes the keywords 'TH_GRID' and 'A'. TH_GRID is validated to be a
        1-D, unique grid and must start at 0. The provided A is squeezed (and 
        should be converted to cm^2 if needed). The geometric factor G is 
        computed via numerical integration over -cos(theta).

        @param kwargs Keyword arguments including:
               - TH_GRID: 1-D numeric grid for theta (must start at 0).
               - A: Numeric response values corresponding to TH_GRID.
        @exception ValueError Raised if TH_GRID is not a valid grid or does not
                              start at 0.
        @exception ArgSizeError Raised if TH_GRID and A do not have the same 
                                shape.
        """
        super().__init__(**kwargs)
        keyword_check_numeric('TH_GRID', 'A', **kwargs)
        self.TH_GRID = squeeze(kwargs['TH_GRID'])
        if not validate_grid(self.TH_GRID):
            raise ValueError('TH_GRID is not a valid grid: 1-d, unique')
        if (self.TH_GRID[0] != 0):
            raise ValueError('TH_GRID should start at 0')
        self._A = squeeze(kwargs['A'])  # TODO: convert to cm^2
        self._Ainterpolator = None
        if self.TH_GRID.shape != self._A.shape:
            raise ArgSizeError('TH_GRID and A are not the same shape')
        self.G = np.trapz(self._A, x=-cosd(self.TH_GRID)) * 2 * np.pi

    def A(self, theta, phi):
        """*INHERIT*"""
        if self._Ainterpolator is None:
            # Extrapolate beyond grid with nearest edge with zero.
            self._Ainterpolator = interp1d(self.TH_GRID, self.A, 'linear',
                                           bounds_error=False, fill_value=0.0,
                                           assume_sorted=True)
        A = self._Ainterpolator(theta)
        return np.broadcast_to(A, np.broadcast(theta, phi).shape)
        

class AR_Pinhole(AngleResponse):
    """
    @brief Class representing pinhole angular responses (delta function at 
           THETA = 0).

    Initialization requires the geometric factor G.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):            
        return keyword_check_bool('TH_TYPE', **kwargs) and \
            (kwargs['TH_TYPE'] == 'PINHOLE')
            
    def __init__(self, **kwargs):
        """
        @brief Initialize the AR_Pinhole instance.

        Processes the keyword 'G' (converting it to a numeric value in cm^2·sr
        as needed) and sets the detector area to zero.
        
        @param kwargs Keyword arguments including:
               - G: Geometric factor.
        """
        super().__init__(**kwargs)  # handles BIDIRECTIONAL
        keyword_check_numeric('G', **kwargs)
        self.G = squeeze(kwargs['G'])  # TODO: convert to cm^2 sr  
        self.area = 0

    def A(self, theta, phi):
        """*INHERIT*"""
        if self.bidirectional:
            A = self.G * ((theta == 0) | (theta == 180)) / 2
        else:
            A = self.G * (theta == 0)
        return np.broadcast_to(A, np.broadcast(theta, phi).shape)

    def hAthetaphi(self, thetagrid, phigrid, **kwargs):
        """*INHERIT*"""
        # Compute interpolation weight in -cos(theta) so that integration over
        # dcos(theta) gives unity.
        w = interp_weights_1d(-cosd(thetagrid), -1.0)  # -cos(0) = -1
        if self.bidirectional:
            w += interp_weights_1d(-cosd(thetagrid), 1.0)  # -cos(180) = 1
        w = np.broadcast_to(w, np.broadcast(thetagrid, phigrid).shape)
        tmp = w / w.sum() * self.G  # normalize so the integral equals G.
        return tmp

    def hAtheta(self, thetagrid, phigrid=None, **kwargs):
        """*INHERIT*"""
        # Compute interpolation weight in -cos(theta) so that integration over
        # dcos(theta) gives unity.
        w = interp_weights_1d(-cosd(thetagrid), -1.0)  # -cos(0) = -1
        if self.bidirectional:
            w += interp_weights_1d(-cosd(thetagrid), 1.0)  # -cos(180) = 1
        tmp = w / w.sum() * self.G  # normalize to sum to G.
        return tmp

    def hAalphabeta(self, alpha0, beta0, phib,
                    alphagrid, betagrid, tgrid=None, **kwargs):
        """*INHERIT*"""
        dt = np.atleast_1d(make_deltas(tgrid, **kwargs))
        alpha0 = np.broadast_to(alpha0, dt.shape)
        beta0 = np.broadast_to(beta0, dt.shape)
        
        # Compute interpolation weight in -cos(alpha) so that the integration
        # over dcos(alpha) gives unity.
        dha = interp_weights_1d(-cosd(alphagrid), -cosd(alpha0))
        # (len(tgrid) x len(alphagrid))
        
        dhb = interp_weights_1d(betagrid, beta0)
        # (len(tgrid) x len(betagrid))
        
        if self.bidirectional:
            dha = (dha + interp_weights_1d(-cosd(alphagrid),
                                           -cosd(180 - alpha0))) / 2
            dhb = (dhb + interp_weights_1d(betagrid,
                                           (beta0 + 180.0) % 360.0)) / 2
            
        # Reshape to enable broadcasting.
        dha.shape = (len(dt), len(alphagrid), 1)
        dhb.shape = (len(dt), 1, len(betagrid))
        dt.reshape = (len(dt), 1, 1)
        h = (dt * dha * dhb).sum(axis=0)  # integrate over time.
        h = h * self.G * sum(dt) / sum(h)  # force to integrate to G·dt.
        return h

    def hAalpha(self, alpha0, beta0, phib,
                alphagrid, betagrid=None, tgrid=None, **kwargs):
        """*INHERIT*"""
        dt = np.atleast_1d(make_deltas(tgrid, **kwargs))
        alpha0 = np.broadast_to(alpha0, dt.shape)
        
        # Compute interpolation weight in -cos(alpha) so that integration over
        # dcos(alpha) gives unity.
        dh = interp_weights_1d(-cosd(alphagrid), -cosd(alpha0))
        # (len(tgrid) x len(alphagrid))
        
        if self.bidirectional:
            dh = (dh + interp_weights_1d(-cosd(alphagrid),
                                         -cosd(180 - alpha0))) / 2
        h = np.dot(dt, dh)  # integrate over time.
        h = h * self.G * sum(dt) / sum(h)  # force to integrate to G·dt.
        return h


class AR_Tele_Rect(AngleResponse):
    """
    @brief Class representing a telescope with rectangular elements.

    This class models the angular response of a telescope with rectangular 
    elements using Sullivan's 1971 paper (updated circa 2010), primarily based
    on equation (13). Initialization requires W1, H1, W2, H2, and D.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return keyword_check_bool('TP_TYPE', **kwargs) and \
            (kwargs['TP_TYPE'] == 'RECT_TELE')
            
    def __init__(self, **kwargs):
        """
        @brief Initialize the AR_Tele_Rect instance.

        Processes numeric keywords 'W1', 'H1', 'W2', 'H2', and 'D'. Computes
        several intermediate parameters (alpha, beta, gamma, delta) and then
        uses Sullivan's equation (11) to compute the geometric factor G.

        @param kwargs Keyword arguments including:
               - W1, H1: Dimensions of the first rectangular element.
               - W2, H2: Dimensions of the second rectangular element.
               - D: Distance parameter.
        """
        super().__init__(**kwargs)  # handles BIDIRECTIONAL
        keyword_check_numeric('W1', 'H1', 'W2', 'H2', **kwargs)
        self.W1 = squeeze(kwargs['W1'])  # TODO: convert to cm
        self.H1 = squeeze(kwargs['H1'])  # TODO: convert to cm
        self.W2 = squeeze(kwargs['W2'])  # TODO: convert to cm
        self.H2 = squeeze(kwargs['H2'])  # TODO: convert to cm
        self.D = squeeze(kwargs['D'])    # TODO: convert to cm
        
        # Initialize some Sullivan variables (in radians).
        alpha = (self.H1 + self.H2) / 2
        beta = (self.W1 + self.W2) / 2
        gamma = (self.H1 - self.H2) / 2
        delta = (self.W1 - self.W2) / 2
        
        # Sullivan's equation (11):
        self.G = self.D**2*np.log(
                ((self.D**2+alpha**2+delta**2)/(self.D**2+alpha**2+beta**2))  
                *(self.D**2+gamma**2+beta**2)/(self.D**2+gamma**2+delta**2)) \
                +2*alpha*np.sqrt(self.D**2+beta**2) \
                *np.arctan(alpha/np.sqrt(self.D**2+beta**2)) \
                +2*beta*np.sqrt(self.D**2+alpha**2) \
                *np.arctan(beta/np.sqrt(self.D**2+alpha**2)) \
                -2*alpha*np.sqrt(self.D**2+delta**2) \
                *np.arctan(alpha/np.sqrt(self.D**2+delta**2)) \
                -2*beta*np.sqrt(self.D**2+gamma**2) \
                *np.arctan(beta/np.sqrt(self.D**2+gamma**2)) \
                -2*gamma*np.sqrt(self.D**2+beta**2) \
                *np.arctan(gamma/np.sqrt(self.D**2+beta**2)) \
                -2*delta*np.sqrt(self.D**2+alpha**2) \
                *np.arctan(delta/np.sqrt(self.D**2+alpha**2)) \
                +2*gamma*np.sqrt(self.D**2+delta**2) \
                *np.arctan(gamma/np.sqrt(self.D**2+delta**2)) \
                +2*delta*np.sqrt(self.D**2+gamma**2) \
                *np.arctan(delta/np.sqrt(self.D**2+gamma**2))


        
    def _X(self, zeta, a1, a2):
        """
        @brief Compute Sullivan's X(zeta,a1,a2).

        @param zeta Numeric value.
        @param a1 Numeric parameter.
        @param a2 Numeric parameter.
        @return The computed X value.
        """
        return min(zeta + a1 / a2, a2 / 2) - max(zeta - a1 / a2, -a2 / 2)

    def _A(self, theta, phi):
        """
        @brief Compute the forward-only A.

        This function computes A (forward response only) based on the input 
        angles.

        @param theta Polar angle (in degrees).
        @param phi Azimuthal angle (in degrees).
        @return A as an array with the broadcast shape of theta and phi.
        """
        A = np.zeros(np.broadcast(theta, phi).shape)
        thetarads = np.radians(theta)
        phirads = np.radians(phi)
        tantheta = np.tan(thetarads)
        zeta = self.D * tantheta * np.cos(phirads)
        eta = self.D * tantheta * np.sin(phirads)
        X = self._X(zeta, self.W1, self.W2)
        Y = self._X(eta, self.H1, self.H2)
        costheta = np.cos(thetarads)

        # apply Heaviside implicitly and resolve tand(90)=inf
        f = (costheta > 0) & (X > 0) & (Y > 0)
        
        if any(f):
            # Equation (13) without the Heaviside functions.
            A[f] = costheta[f] * X[f] * Y[f]
        

    def A(self, theta, phi):
        """*INHERIT*"""
        A = self._A(theta, phi)
        if self.bidirectional:
            A += self._A(90 - theta, (phi + 180.0) % 360.0)
        return A


class AR_Table_asym(AngleResponse):
    """
    @brief Class representing phi-asymmetric tabular angular responses.

    This class represents a phi-asymmetric tabular angular response. 
    Initialization requires TH_GRID, PH_GRID, and A. Bidirectional effects are
    ignored; for a bidirectional tabular response, provide A on TH_GRID that
    spans 0 to 180 degrees.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return keyword_check_bool('TP_TYPE', **kwargs) and \
            (kwargs['TP_TYPE'] == 'TBL')
            
    def __init__(self, **kwargs):
        """
        @brief Initialize the AR_Table_asym instance.

        Processes the keywords 'TH_GRID', 'PH_GRID', and 'A'. Both TH_GRID and
        PH_GRID are validated to be 1-D, unique grids; TH_GRID must start at 0
        and PH_GRID must span [0, 360]. The shape of A must match 
        (len(TH_GRID), len(PH_GRID)). The geometric factor G is computed via 
        double trapezoidal integration.

        @param kwargs Keyword arguments including:
               - TH_GRID: 1-D numeric grid for theta.
               - PH_GRID: 1-D numeric grid for phi.
               - A: 2-D numeric response array corresponding to (TH_GRID, 
                    PH_GRID).
        @exception ValueError Raised if TH_GRID or PH_GRID is invalid or not 
                              spanning required ranges.
        @exception ArgSizeError Raised if the shape of A does not match 
                                (len(TH_GRID), len(PH_GRID)).
        """
        super().__init__(**kwargs)
        keyword_check_numeric('TH_GRID', 'PH_GRID', 'A', **kwargs)
        for arg in ['TH_GRID', 'PH_GRID']:
            setattr(self, arg, squeeze(kwargs[arg]))
            if not validate_grid(getattr(self, arg)):
                raise ValueError('%s is not a valid grid: 1-d, unique' % arg)
        self._A = squeeze(kwargs['A'])  # TODO: convert to cm^2
        self._Ainterpolator = None
        if self._A.shape != (len(self.TH_GRID), len(self.PH_GRID)):
            raise ArgSizeError('A does not have shape (TH_GRID, PH_GRID)')
        if (self.TH_GRID[0] != 0):
            raise ValueError('TH_GRID should start at 0')
        if (self.PH_GRID[0] != 0) or (self.PH_GRID[-1] != 360):
            raise ValueError('PH_GRID should span [0,360]')
        self.G = np.trapz(np.trapz(self._A, x=np.radians(self.PH_GRID),
                                   axis=1), x=-cosd(self.TH_GRID), axis=0)

    def A(self, theta, phi):
        """*INHERIT*"""
        if self._Ainterpolator is None:
            ## PH_GRID over [0, 360] restricted at init.
            ## Extrapolate beyond grid at nearest edge with fill zero.
            self._Ainterpolator = \
                RegularGridInterpolator((self.TH_GRID, self.PH_GRID),
                                        self._A, 'linear', bounds_error=False,
                                        fill_value=0.0)
            
        ## Restrict phi to [0, 360].
        if phi < 0 or phi > 360:
            phi = np.mod(phi, 360)

        return self._Ainterpolator((theta, phi))

class ChannelResponse(FactoryConstructorMixin):
    """
    @brief Virtual base class for channel responses.
    
    This class serves as an abstract base for channel responses. Its 
    constructor acts as a factory, returning the appropriate subclass based on
    the supplied response dictionary.
    
    The following methods are defined:
      - Working in (E, alpha, beta, time) coordinates (results integrated over
        time):
        * hEiso      - Weights for numerical integration over E (isotropic).
        * hEalpha    - Weights for double numerical integration over E and 
                       alpha (gyrotropic).
        * hEalphabeta- Weights for triple numerical integration over E, alpha, 
                       and beta.
        * halpha     - Weights for numerical integration over alpha (flat 
                       spectrum, gyrotropic).
        * halphabeta - Weights for double numerical integration over alpha and
                       beta (flat spectrum).
      
      - Working in (E, theta, phi) coordinates:
        * R          - Response function (effective area, cm^2).
        * hE         - Weights for numerical integration over E (isotropic).
        * hEtheta    - Weights for double numerical integration over E and 
                       theta (ignoring phi dependence).
        * hEthetaphi - Weights for triple numerical integration over E, theta,
                       and phi.
        * htheta     - Weights for numerical integration over theta (flat 
                       spectrum, ignoring phi dependence).
        * hthetaphi  - Weights for double numerical integration over theta and 
                       phi (flat spectrum).
    """
    
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """*INHERIT*"""
        keyword_check('RESP_TYPE', **kwargs)
        raise KeywordError('Supplied keywords to not define a recognized' + \
                           'ChannelResponse')
    
    def __init__(self, *args, **kwargs):
        """
        @brief Construct a ChannelResponse instance.
        
        The constructor captures the original keyword dictionary and assigns 
        any additional key-value pairs (not already set as attributes) to the 
        instance.
        
        @param args Positional arguments.
        @param kwargs Keyword arguments used to build the response.
        """
        super().__init__(*args, **kwargs)
        self._original_kwargs = {**kwargs}  # copy original dictionary
        
        # Capture any remaining key-value pairs as attributes.
        for key, val in kwargs.items():
            if not hasattr(self, key):
                setattr(self, key, val)
    
    def copy(self):
        """
        @brief Return a copy of this ChannelResponse.
        
        @return A new ChannelResponse object constructed from the original
                keyword dictionary.
        """
        return ChannelResponse(**self.to_dict())
    
    def to_dict(self):
        """
        @brief Return the original dictionary of keywords.
        
        @return The dictionary of keywords used to build this object.
        """
        return self._original_kwargs
    
    def R(self, E, theta, phi):
        """
        @brief Evaluate the 3-D response function.
        
        Returns the response (effective area multiplied by efficiency) in cm^2
        at the specified energy, polar angle, and azimuthal angle. The inputs
        E, theta, and phi must broadcast together.
        
        @param E Energy value(s) (MeV).
        @param theta Polar angle(s) in degrees.
        @param phi Azimuthal angle(s) in degrees.
        @return The response function with shape matching the broadcast of E,
                theta, and phi.
        @exception NotImplementedError If this method is not implemented in 
                                       a subclass.
        """
        raise NotImplementedError('Class %s ' % self.__class__.__name__ + \
                                  'did not overload of R method, as required')
    
    def hEthetaphi(self, Egrid, thetagrid, phigrid, **kwargs):
        """
        @brief Compute response weights on an E x theta x phi grid.
        
        Calculates response weights (in cm^2 sr MeV) using differential grid
        elements. If Egrid is None, the object will attempt to obtain its own 
        (e.g., from E_GRID attribute).
        
        @param Egrid 1-D numpy array of energy grid values (MeV).
        @param thetagrid 1-D numpy array of theta values (degrees).
        @param phigrid 1-D numpy array of phi values (degrees).
        @param kwargs Additional keyword arguments for integration (e.g., 
                      delta options).
        @return A scalar or a numpy array with shape (len(Egrid), 
                len(thetagrid), len(phigrid)).
        """
        if Egrid is None:
            if hasattr(self, 'E_GRID'):
                Egrid = self.E_GRID
            else:
                raise ValueError('Egrid input cannot be None for object ' + \
                                 'without its own E_GRID')
        dE = make_deltas(Egrid, **kwargs)
        dcostheta = np.abs(make_deltas(-cosd(thetagrid), **kwargs))
        dphi = make_deltas(np.radians(phigrid), **kwargs)

        ## New version: added mesh=True in the next two lines.
        EI, thetaI, phiI = broadcast_grids(Egrid, thetagrid, phigrid,
                                           mesh=True)
        dEI, dcosthetaI, dphiI = broadcast_grids(dE, dcostheta, dphi,
                                                 mesh=True)
        
        return self.R(EI, thetaI, phiI) * dEI * dcosthetaI * dphiI
    
    def hEtheta(self, Egrid, thetagrid, phigrid=None, **kwargs):
        """
        @brief Compute response weights on an E x theta grid.
        
        Integrates the full E x theta x phi response over phi to yield weights
        (in cm^2 sr MeV) on an E x theta grid. If phigrid is not provided, it 
        defaults to the object's PH_GRID or a default phi grid.
        
        @param Egrid 1-D numpy array of energy values (MeV).
        @param thetagrid 1-D numpy array of theta values (degrees).
        @param phigrid Optional 1-D numpy array of phi values (degrees). 
                       Defaults to 0–360.
        @param kwargs Additional integration options.
        @return A scalar or numpy array with shape (len(Egrid),
                len(thetagrid)).
        """
        if phigrid is None:
            if hasattr(self, 'PH_GRID'):
                phigrid = self.PH_GRID
            else:
                phigrid = default_phigrid(**kwargs)
        h = self.hEthetaphi(Egrid, thetagrid, phigrid, **kwargs)

        ## New version: changed axis=1 to axis=2 (fixed bug)
        return h.sum(axis=2)  # integrate over phi
    
    def hE(self, Egrid, thetagrid=None, phigrid=None, **kwargs):
        """
        @brief Compute response weights on an E grid.
        
        Obtains response weights (in cm^2 sr MeV) by integrating over theta 
        (0–180) and phi (0–360). If thetagrid is not provided, it defaults to
        the object's TH_GRID attribute or a default grid.
        
        @param Egrid 1-D numpy array of energy values (MeV).
        @param thetagrid Optional 1-D numpy array of theta values (degrees).
        @param phigrid Optional 1-D numpy array of phi values (degrees).
        @param kwargs Additional options.
        @return A scalar or 1-D numpy array (length equals len(Egrid)).
        """
        if thetagrid is None:
            if hasattr(self, 'TH_GRID'):
                thetagrid = self.TH_GRID
            else:
                thetagrid = default_thetagrid(**kwargs)
        h = self.hEtheta(Egrid, thetagrid, phigrid, **kwargs)
        return h.sum(axis=1)  # integrate over theta
    
    def hEalphabeta(self, alpha0, beta0, phib, Egrid,
                    alphagrid, betagrid, tgrid=None, **kwargs):
        """
        @brief Compute response weights on an E x alpha x beta grid.
        
        Returns weights (in cm^2·sr·MeV or cm^2·sr·MeV·s if tgrid is supplied)
        computed via a triple numerical integration. If Egrid is None, the 
        object will attempt to supply its own.
        
        @param alpha0 The alpha0 parameter.
        @param beta0 The beta0 parameter.
        @param phib The phib parameter.
        @param Egrid 1-D numpy array of energy values (MeV).
        @param alphagrid 1-D numpy array of alpha values (degrees).
        @param betagrid 1-D numpy array of beta values (degrees).
        @param tgrid Optional time grid; if provided, integration over time is
                     performed.
        @param kwargs Additional integration options.
        @return A scalar or numpy array with shape (len(Egrid), len(alphagrid),
                len(betagrid)).
        """
        dE = make_deltas(Egrid, **kwargs)  # (NE,)
        dcosa = make_deltas(-cosd(alphagrid), **kwargs)  # (Na,)
        db = np.radians(make_deltas(betagrid, **kwargs))  # (Nb,)
        dt = make_deltas(tgrid, **kwargs)  # (Nt,) or scalar

        ## Old version:
        #dE, dcosa, db = broadcast_grids(dE, dcosa, db)  # shape (NE, Na, Nb)
        #theta, phi = alphabeta2thetaphi(alphagrid, betagrid, alpha0, beta0,
        #                                phib)  # (Na, Nb, Nt) or (Na, Nb)
        #E = np.reshape(Egrid, (len(Egrid),
        #                       *np.ones(theta.ndim)))  # (NE, 1, 1, *1)

        ## New version:
        dE, dcosa, db, dt = broadcast_grids(dE, dcosa, db, dt,
                                            mesh=True) # (NE,Na,Nb, dt)
        alphagrid, betagrid = broadcast_grids(alphagrid, betagrid, mesh=True)
        theta, phi = alphabeta2thetaphi(alphagrid, betagrid, alpha0,
                                        beta0, phib) # (Na,Nb,Nt) or (NA,Nb)
        E = np.reshape(Egrid, (len(Egrid),
                               *np.ones(theta.ndim, dtype=int))) # (NE,1,1,*1)
        
        theta = np.expand_dims(theta, axis=0)  # (1, Na, Nb, *Nt)
        phi = np.expand_dims(phi, axis=0)  # (1, Na, Nb, *Nt)

        ## Old way (may contain bugs).
        #R_val = self.R(E, theta, phi)  # (NE, Na, Nb, *Nt)
        #if np.isscalar(dt):
        #    h = R_val * dt
        #else:
        #    # Take the dot product along Nt dimension.
        #    # AI replaced h in the function call below with R_val, since h is
        #    # not declared in this function previously. Is this a bug from the
        #    # original versin of rfl.py?
        #    h = np.tensordot(h, dt, axes=((3,), (0,)))
        #h = h * dE * dcosa * db / self.CROSSCALIB
        # now h is (NE,Na,Nb)
        
        ## New way (untested).
        R_val = self.R(E, theta, phi) # (NE, Na, Nb, *Nt) 
        h = R_val * dt * dE * dcosa * db / self.CROSSCALIB
        h = h.sum(axis=3) # integrate over time
        
        return h
    
    def hEalpha(self, alpha0, beta0, phib,
                Egrid, alphagrid, betagrid=None, tgrid=None, **kwargs):
        """
        @brief Compute response weights on an E x alpha grid (integrated over 
               beta).
        
        This integrates the triple integral of hEalphabeta over beta.
        If betagrid is None, a default numeric beta grid is provided.
        
        @param alpha0 The alpha0 parameter.
        @param beta0 The beta0 parameter.
        @param phib The phib parameter.
        @param Egrid 1-D numpy array of energy values (MeV).
        @param alphagrid 1-D numpy array of alpha values (degrees).
        @param betagrid Optional 1-D numpy array of beta values (degrees).
        @param tgrid Optional time grid.
        @param kwargs Additional integration options.
        @return A scalar or 1-D numpy array with length equal to 
                len(alphagrid).
        """
        if betagrid is None:
            betagrid = default_betagrid(**kwargs)
        #print(betagrid[0], betagrid[-1])
        
        h = self.hEalphabeta(alpha0, beta0, phib,
                             Egrid, alphagrid, betagrid, tgrid, **kwargs)
        return h.sum(axis=2)  # integrate over beta
    
    def hEiso(self, alpha0, beta0, phib,
              Egrid, alphagrid, betagrid=None, tgrid=None, **kwargs):
        """
        @brief Compute isotropic response weights on an alpha grid (integrated
               over beta).
        
        This is designed to yield isotropic weights, integrating hEalpha over 
        alpha. If alphagrid is None, a default alpha grid is used.
        
        @param alpha0 The alpha0 parameter.
        @param beta0 The beta0 parameter.
        @param phib The phib parameter.
        @param Egrid 1-D numpy array of energy values (MeV).
        @param alphagrid 1-D numpy array of alpha values (degrees). Defaults 
                         to a default if None.
        @param betagrid Optional 1-D numpy array of beta values (degrees).
        @param tgrid Optional time grid.
        @param kwargs Additional options.
        @return A scalar or 1-D numpy array (function of alphagrid) 
                representing the weights.
        """
        if alphagrid is None:
            alphagrid = default_alphagrid(**kwargs)
        #print(alphagrid[0], alphagrid[-1])
        
        h = self.hEalpha(alpha0, beta0, phib,
                         Egrid, alphagrid, betagrid, tgrid, **kwargs)
        return h.sum(axis=1)  # integrate over alpha
    
    def hthetaphi(self, Egrid, thetagrid, phigrid, **kwargs):
        """
        @brief Compute angular response weights on a theta x phi grid for a 
               flat spectrum.
        
        This method returns weights (cm^2 sr MeV) by integrating over energy 
        (for a flat spectrum). If Egrid is None, the object will supply its 
        own.
        
        @param Egrid 1-D numpy array of energy values (MeV).
        @param thetagrid 1-D numpy array of theta values (degrees).
        @param phigrid 1-D numpy array of phi values (degrees).
        @param kwargs Additional integration options.
        @return A scalar or 2-D numpy array with shape (len(thetagrid), 
                len(phigrid)).
        """
        h = self.hEthetaphi(Egrid, thetagrid, phigrid, **kwargs)
        return h.sum(axis=0)  # integrate over E
    
    def htheta(self, Egrid, thetagrid, phigrid=None, **kwargs):
        """
        @brief Compute response weights on a theta grid for a flat spectrum.
        
        Returns weights (cm^2 sr MeV) integrated over E and phi. If phigrid is
        None, it defaults to a full 0–360 range.
        
        @param Egrid 1-D numpy array of energy values (MeV).
        @param thetagrid 1-D numpy array of theta values (degrees).
        @param phigrid Optional 1-D numpy array of phi values (degrees).
        @param kwargs Additional options.
        @return A scalar or 1-D numpy array with length equal to 
                len(thetagrid).
        """
        h = self.hEtheta(Egrid, thetagrid, phigrid, **kwargs)
        return h.sum(axis=0)  # integrate over E
    
    def halphabeta(self, alpha0, beta0, phib,
                   Egrid, alphagrid, betagrid, tgrid=None, **kwargs):
        """
        @brief Compute response weights on an alpha x beta grid for a flat 
               spectrum.
        
        The integration is performed over energy to yield weights with units 
        cm^2·sr·MeV (or cm^2·sr·MeV·s if tgrid is provided).
        
        @param alpha0 The alpha0 parameter.
        @param beta0 The beta0 parameter.
        @param phib The phib parameter.
        @param Egrid 1-D numpy array of energy values (MeV).
        @param alphagrid 1-D numpy array of alpha values (degrees).
        @param betagrid 1-D numpy array of beta values (degrees).
        @param tgrid Optional time grid.
        @param kwargs Additional options.
        @return A scalar or numpy array with shape (len(alphagrid), 
                len(betagrid)).
        """
        h = self.hEalphabeta(alpha0, beta0, phib,
                             Egrid, alphagrid, betagrid, tgrid, **kwargs)
        return h.sum(axis=0)  # integrate over E
    
    def halpha(self, alpha0, beta0, phib,
               Egrid, alphagrid, betagrid=None, tgrid=None, **kwargs):
        """
        @brief Compute response weights on an alpha grid for a flat spectrum.
        
        This returns weights (cm^2·sr·MeV or with time integration) by 
        integrating over energy and beta.
        
        @param alpha0 The alpha0 parameter.
        @param beta0 The beta0 parameter.
        @param phib The phib parameter.
        @param Egrid 1-D numpy array of energy values (MeV).
        @param alphagrid 1-D numpy array of alpha values (degrees).
        @param betagrid Optional 1-D numpy array of beta values (degrees).
        @param tgrid Optional time grid.
        @param kwargs Additional options.
        @return A scalar or 1-D numpy array with length equal to
                len(alphagrid).
        """
        h = self.halphabeta(alpha0, beta0, phib,
                            Egrid, alphagrid, betagrid, tgrid, **kwargs)
        return h.sum(axis=1)  # integrate over beta

class CR_Esep(ChannelResponse):
    """
    @brief Abstract base class for all energy-separable responses.
    """
    
    @classmethod
    def is_mine(cls, *args, **kwargs):
        """
        @brief Return True if the provided initialization data claims this 
               subclass.
        
        For CR_Esep the abstract base class returns False because it is never
        the correct answer.
        """
        return False  # abstract class is never the right answer

    def __init__(self, **kwargs):
        """
        @brief Construct a CR_Esep instance.
        
        Initializes the constituent energy response (er) and angular response
        (ar). Copies any properties from the constituents to this object if not
        already defined.
        
        @param kwargs Keyword arguments used by the EnergyResponse and 
                      AngleResponse constructors.
        """
        super().__init__(**kwargs)
        self.er = EnergyResponse(**kwargs)
        self.ar = AngleResponse(**kwargs)
        # Copy properties from constituent objects, if present.
        for sub in [self.er, self.ar]:
            for key in sub.__dict__:
                if key not in self.__dict__:
                    setattr(self, key, getattr(sub, key))

    def R(self, E, theta, phi):
        """*INHERIT*"""
        return self.er.RE(E) * self.ar.A(theta, phi)

    def _merge_hEhA(hE, hA):
        """
        @brief Merge hE (size NE,) and hA (size N1, or (N1,N2)) into a 
               combined array.
        
        The resulting array hEA is of size (NE, N1) or (NE, N1, N2) and is 
        computed as the product hE*hA.
        
        @param hE 1-D array with NE elements.
        @param hA Array with shape (N1,) or (N1,N2).
        @return The merged array of shape (NE, N1) or (NE, N1, N2).
        """
        # Reshape hE to (hE.size, 1, 1, ... ) with number of dimensions
        # matching hA.ndim.
        hE = np.reshape(hE, (hE.size, np.ones(hA.ndim, dtype=int)))
        hA = np.reshape(hE, (1, *hA.shape))
        return hE * hA

    def hEthetaphi(self, Egrid, thetagrid, phigrid, **kwargs):
        """*INHERIT*"""
        if Egrid is None:
            if hasattr(self, 'E_GRID'):
                Egrid = self.E_GRID
            else:
                raise ValueError('Egrid input cannot be None for object' + \
                                 'without its own E_GRID')
        hE = self.er.hE(Egrid, **kwargs)
        hA = self.ar.hAthetaphi(thetagrid, phigrid, **kwargs)
        return self._merge_hEhA(hE, hA)

    def hEtheta(self, Egrid, thetagrid, phigrid=None, **kwargs):
        """*INHERIT*"""
        if phigrid is None:
            if hasattr(self, 'PH_GRID'):
                phigrid = self.PH_GRID
            else:
                phigrid = default_phigrid(**kwargs)
        hE = self.er.hE(Egrid, **kwargs)
        hA = self.ar.hAtheta(thetagrid, phigrid, **kwargs)
        return self._merge_hEhA(hE, hA)

    def hE(self, Egrid, thetagrid=None, phigrid=None, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE(Egrid, **kwargs)
        return hE * self.ar.hA0

    def hEalphabeta(self, alpha0, beta0, phib,
                    Egrid, alphagrid, betagrid, tgrid=None, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE(Egrid, **kwargs)
        hA = self.ar.hAalphabeta(alpha0, beta0, phib,
                                 alphagrid, betagrid, tgrid, **kwargs)
        return self._merge_hEhA(hE, hA)

    def hEalpha(self, alpha0, beta0, phib,
                Egrid, alphagrid, betagrid=None, tgrid=None, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE(Egrid, **kwargs)
        hA = self.ar.hAalpha(alpha0, beta0, phib,
                             alphagrid, betagrid, tgrid, **kwargs)
        return self._merge_hEhA(hE, hA)

    def hEiso(self, alpha0, beta0, phib,
              Egrid, alphagrid, betagrid=None, tgrid=None, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE(Egrid, **kwargs)
        hA = self.ar.hA0 * sum(make_deltas(tgrid, **kwargs))
        return hE * hA

    def hthetaphi(self, Egrid, thetagrid, phigrid, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid, **kwargs)
        hA = self.ar.hAthetaphi(thetagrid, phigrid, **kwargs)
        return hE * hA

    def htheta(self, Egrid, thetagrid, phigrid=None, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid, **kwargs)
        hA = self.ar.hAtheta(thetagrid, phigrid, **kwargs)
        return hE * hA

    def halphabeta(self, alpha0, beta0, phib,
                   Egrid, alphagrid, betagrid, tgrid=None, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid, **kwargs)
        hA = self.ar.hAalphabeta(alpha0, beta0, phib,
                                 alphagrid, betagrid, tgrid, **kwargs)
        return hE * hA

    def halpha(self, alpha0, beta0, phib,
               Egrid, alphagrid, betagrid=None, tgrid=None, **kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid, **kwargs)
        hA = self.ar.hAalpha(alpha0, beta0, phib,
                             alphagrid, betagrid, tgrid, **kwargs)
        return hE * hA


class CR_Omni(CR_Esep):
    """
    @brief Omnidirectional channel response.
    
    This subclass represents an omnidirectional channel response.
    """
    
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return ('RESP_TYPE' in kwargs) and (kwargs['RESP_TYPE'] == '[E]')
    # No further customization required.

    
class CR_Esep_sym(CR_Esep):
    """
    @brief Energy-separable, phi-symmetric channel response.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return ('RESP_TYPE' in kwargs) and (kwargs['RESP_TYPE'] == '[E],[TH]')
    # No further customization required.

    
class CR_Esep_asym(CR_Esep):
    """
    @brief Energy-separable, phi-asymmetric channel response.
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return ('RESP_TYPE' in kwargs) and \
            (kwargs['RESP_TYPE'] == '[E],[TH,PH]')
    # No further customization required.

    
class CR_Table_sym(ChannelResponse):
    """
    @brief Energy-inseparable, phi-symmetric channel response.
    
    The response is defined by a 2-D table. The initializer requires ET_TYPE,
    E_GRID, TH_GRID, and R.
    """
    
    @classmethod
    def is_mine(cls, *args, **kwargs):
        keyword_check('RESP_TYPE', **kwargs)
        return ('RESP_TYPE' in kwargs) and \
               (kwargs['RESP_TYPE'] == '[E,TH]') and \
               ('ET_TYPE' in kwargs) and \
               (kwargs['ET_TYPE'] == 'TBL')
    
    def __init__(self, **kwargs):
        """
        @brief Construct a CR_Table_sym instance.
        
        Checks that the grids E_GRID and TH_GRID are valid (1-D and unique) and
        squeezes the response table R. If R does not have shape (E_GRID.size,
        TH_GRID.size), an ArgSizeError is raised.
        
        @param kwargs Expected keywords include 'E_GRID', 'TH_GRID', and 'R'.
        @exception ValueError if any grid is invalid.
        @exception ArgSizeError if R does not match grid dimensions.
        """
        super().__init__(**kwargs)
        keyword_check_numeric('E_GRID', 'TH_GRID', 'R', **kwargs)
        for arg in ['E_GRID', 'TH_GRID']:
            setattr(self, arg, squeeze(kwargs[arg]))
            if not validate_grid(getattr(self, arg)):
                raise ValueError('%s is not a valid grid: 1-d, unique' % arg)
        # TODO: Convert E_GRID to MeV.
        self._R = squeeze(kwargs['R'])  # TODO: Convert to cm^2.
        if self._R.shape != (self.E_GRID.size, self.TH_GRID.size):
            raise ArgSizeError('R is not shape (E_GRID x TH_GRID)')
        self._Rinterpolator = None

    def R(self, E, theta, phi):
        """*INHERIT*"""
        if self._Rinterpolator is None:
            # Extrapolate beyond grid with fill value of zero.
            self._Rinterpolator = \
                RegularGridInterpolator((self.E_GRID, self.TH_GRID), self._R,
                                        'linear', bounds_error=False,
                                        fill_value=0.0)
        return self._Rinterpolator((E, theta))

class CR_Table_asym(ChannelResponse):
    """
    @brief Energy-inseparable, phi-asymmetric channel response defined by a
           3-D table.

    This class represents an energy-inseparable, phi-asymmetric channel
    response where the response is defined by a three-dimensional table. The
    initializer requires the following keywords:
      - ETP_TYPE must be 'TBL'
      - E_GRID: a 1-D numeric energy grid
      - TH_GRID: a 1-D numeric theta grid
      - PH_GRID: a 1-D numeric phi grid (must span [0,360])
      - R: a 3-D response table with shape (len(E_GRID), len(TH_GRID), 
           len(PH_GRID))
    """
    @classmethod
    def is_mine(cls, *args, **kwargs):
        return ('RESP_TYPE' in kwargs) and \
                (kwargs['RESP_TYPE'] == '[E,TH,PH]') and \
                ('ETP_TYPE' in kwargs) and \
                (kwargs['ETP_TYPE'] == 'TBL')

    def __init__(self, **kwargs):
        """
        @brief Initialize a CR_Table_asym instance.

        Validates and stores the energy grid (E_GRID), theta grid (TH_GRID), 
        phi grid (PH_GRID), and the 3-D response table R. Raises errors if any
        grid is invalid or if the response table does not match the dimensions 
        of the grids.

        @param kwargs Keyword arguments including:
          - E_GRID: a 1-D numeric grid for energy.
          - TH_GRID: a 1-D numeric grid for theta.
          - PH_GRID: a 1-D numeric grid for phi.
          - R: a 3-D response table.
        @exception ValueError If any grid is not valid (i.e. not 1-D, unique)
                              or if PH_GRID does not span [0, 360].
        @exception ArgSizeError If R does not have the shape (E_GRID.size, 
                                TH_GRID.size, PH_GRID.size).
        """
        super().__init__(**kwargs)
        keyword_check_numeric('E_GRID', 'TH_GRID', 'PH_GRID', 'R', **kwargs)
        for arg in ['E_GRID', 'TH_GRID', 'PH_GRID']:
            setattr(self, arg, squeeze(kwargs[arg]))
            if not validate_grid(getattr(self, arg)):
                raise ValueError('%s is not a valid grid: 1-d, unique' % arg)
        # TODO: convert E_GRID to MeV
        self._R = squeeze(kwargs['R'])  # TODO: convert to cm^2
        
        if (self._R.shape == (self.PH_GRID.size, self.TH_GRID.size,
                              self.E_GRID.size)) and \
                              (self.PH_GRID.size != self.E_GRID.size):
            ## Transpose automatically if E and PH dimensions are switched.
            print("Warning: R is shape (PH_GRID x TH_GRID x E_GRID), " + \
                  "which is reversed from the expected " + \
                  "(E_GRID x TH_GRID x PH_GRID). Automatically transposing.")
            self._R = np.transpose(self._R, (2, 1, 0))
        if self._R.shape != (self.E_GRID.size, self.TH_GRID.size,
                             self.PH_GRID.size):
            ## Raise exception for mismatched dimensions.
            raise ArgSizeError("R is not shape (E_GRID x TH_GRID x PH_GRID)" +
                               "\nR.shape = {}".format(self._R.shape) +
                               "\n(E, TH, PH) = " +
                               "({}, {}, {})".format(self.E_GRID.size,
                                                     self.TH_GRID.size,
                                                     self.PH_GRID.size))
        self._Rinterpolator = None
        if (self.PH_GRID[0] != 0) or (self.PH_GRID[-1] != 360):
            raise ValueError('PH_GRID should span [0,360]')

    def R(self, E, theta, phi):
        """*INHERIT*"""
        if self._Rinterpolator is None:
            ## Extrapolate beyond the grid with zero.
            self._Rinterpolator = \
                RegularGridInterpolator((self.E_GRID, self.TH_GRID,
                                         self.PH_GRID), self._R, 'linear',
                                        bounds_error=False, fill_value=0.0)
            
        ## Restrict theta to [0, 180] and phi to [0, 360].
        if theta < 0 or theta > 180:
            theta = np.mod(theta, 180)
        if phi < 0 or phi > 360:
            phi = np.mod(phi, 360)

        return self._Rinterpolator((E, theta, phi))


## Inherit docstrings for everything descended from FactoryConstructorMixin,
## including ChannelResponse and its internal classes EnergyResponse and
## AngleResponse.
inherit_docstrings(FactoryConstructorMixin)


def recursive_rename(var, renames):
    """
    @brief Recursively rename dictionary keys and list entries.

    Returns a new variable built from the input variable where all dictionary
    keys and list entries are renamed according to the provided renames
    dictionary. If an element is a string and appears in renames, it is 
    replaced with the new value.
    
    @param var The input variable (can be a list, dict, or any other type).
    @param renames A dictionary mapping old names to new names.
    @return A new variable with keys and list entries renamed.
    """
    if isinstance(var, list):
        out = []
        for x in var:
            if isinstance(x, str) and (x in renames):
                out.append(renames[x])
            else:
                out.append(recursive_rename(x, renames))
    elif isinstance(var, dict):
        out = {}
        for key, val in var.items():
            if key in renames:
                key = renames[key]
            out[key] = recursive_rename(val, renames)
    else:
        out = var        
    return out
    

def load_inst_info(inst_info):
    """
    @brief Load instrument info and construct channel response objects.

    The input inst_info is expected to be a dictionary compliant with the
    response function file format standard data model. This function converts
    individual channel data into corresponding ChannelResponse objects, 
    propagates higher-level data downward, and supports input as a filename.
    
    @param inst_info A dictionary or a filename (assumed to be .json or 
                     .h5/.hdf5; HDF5 not implemented yet).
    @return The instrument info with individual channels replaced by 
            ChannelResponse objects.
    @exception RFLError If the file type is unknown or the input is not a 
                        dictionary.
    """
    if isinstance(inst_info, str):
        if inst_info.endswith('.json'):
            inst_info = read_JSON(inst_info)
        elif inst_info.endswith('.h5') or inst_info.endswith('.hdf5') \
             or inst_info.endswith('.mat'):
            inst_info = read_h5(inst_info)
        else:
            raise RFLError('Unknown file type %s' % inst_info)
    
    if not isinstance(inst_info, dict):
        raise RFLError('Only dictionary supported for inst_info argument')

    renames = {'XCAL': 'CROSSCALIB', 'XCAL_RMSE': 'CROSSCALIB_RMSE'}
    # XCAL to CROSSCALIB rename following upgrade from RFL v1.0.0 to v1.1.0
    # in response to change from format specification v1.0.1 to v1.1.0

    inst_info = recursive_rename(inst_info, renames)

    # Define fields that can propagate down to channel/species responses.
    prop_flds = ['L_UNIT', 'E_UNIT', 'DEAD_TIME_PER_COUNT', 'DEAD_TYPE',
                 'COUNTS_MAX', 'CROSSCALIB', 'CROSSCALIB_RMSE', 'RESP_TYPE',
                 'ETP_TYPE', 'ET_TYPE', 'E_TYPE', 'TP_TYPE', 'TH_TYPE',
                 'E_GRID', 'TH_GRID', 'PH_GRID', 'EPS', 'R', 'A', 'G', 'E0',
                 'E1', 'DE', 'R1', 'R2', 'W1', 'W2', 'H1', 'H2', 'D',
                 'BIDIRECTIONAL']

    for chan in inst_info['CHANNEL_NAMES']:
        # Copy inst_info to channel.
        if 'SPECIES' not in inst_info[chan]:
            inst_info[chan]['SPECIES'] = inst_info['SPECIES']
        for fld in prop_flds:
            if (fld not in inst_info[chan]) and (fld in inst_info):
                inst_info[chan][fld] = inst_info[fld]

        # Copy channel info to species.
        for sp in inst_info[chan]['SPECIES']:
            if isinstance(inst_info[chan][sp], ChannelResponse):
                continue  # already initialized
            for fld in prop_flds:
                if (fld not in inst_info[chan][sp]) and \
                   (fld in inst_info[chan]):
                    inst_info[chan][sp][fld] = inst_info[chan][fld]
                # copy inst_info.(sp) to inst_info.(chan).(sp)
                if (fld not in inst_info[chan][sp]) and (sp in inst_info) and \
                   (fld in inst_info[sp]):
                    inst_info[chan][sp][fld] = inst_info[sp][fld]
    
    for chan in inst_info['CHANNEL_NAMES']:
        for sp in inst_info[chan]['SPECIES']:
            if isinstance(inst_info[chan][sp], ChannelResponse):
                continue  # already initialized
            inst_info[chan][sp] = ChannelResponse(**inst_info[chan][sp])

    return inst_info        


def CRs_to_dicts(inst_info):
    """
    @brief Convert ChannelResponse objects to dictionaries within inst_info.

    Returns a shallow copy of the inst_info dictionary structure where any 
    ChannelResponse objects are converted to dictionaries using their to_dict()
    method. This facilitates writing the data to files.
    
    @param inst_info The instrument info dictionary potentially containing 
                     ChannelResponse objects.
    @return A shallow copy of inst_info with ChannelResponse objects replaced
            by their dict representations.
    """
    inst_dict = {**inst_info}
    sp0 = None
    if 'SPECIES' in inst_info:
        sp0 = inst_info['SPECIES']
    for chan in inst_info['CHANNEL_NAMES']:
        inst_dict[chan] = {**inst_info[chan]}
        sp = sp0
        if 'SPECIES' in inst_info[chan]:
            sp = inst_info[chan]['SPECIES']
        if sp:
            for s in sp:
                if isinstance(inst_info[chan][s], ChannelResponse):
                    inst_dict[chan][s] = inst_info[chan][s].to_dict()
    return inst_dict


def write_JSON(inst_info, jsonfile):
    """
    @brief Write instrument info to a JSON file.

    Converts the instrument info dictionary (with ChannelResponse objects 
    converted to dictionaries) to a JSON string and writes it to the specified
    file.
    
    @param inst_info The instrument info dictionary.
    @param jsonfile The filename for the output JSON file.
    """
    import json
    def encoder(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        elif isinstance(o, (np.integer,)):  # catches numpy.int64, np.int32...
            return int(o)
        elif isinstance(o, (np.floating,)):  # catches numpy.float64, etc.
            return float(o)
        print("Unhandled type:", type(o), "value:", o)
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(o)
    inst_dict = CRs_to_dicts(inst_info)
    with open(jsonfile, 'wt') as f:
        json.dump(inst_dict, f, default=encoder, indent=1)

    return

def read_JSON(jsonfile):
    """
    @brief Read instrument info from a JSON file.

    Reads a JSON file and returns its contents as a dictionary. This function
    does not convert the data into ChannelResponse objects.
    
    @param jsonfile The JSON filename to read.
    @return The instrument info dictionary read from the file.
    """
    import json
    with open(jsonfile, 'r') as f:
        inst_dict = json.load(f)
    return inst_dict

def write_h5(inst_info, filename):
    """
    @brief Write instrument info dict/structure to an HDF5 file.
    
    This function writes the instrument info data structure (which may include
    nested dictionaries, lists, and NumPy arrays) to an HDF5 file. It adds
    attributes to each dataset or group to indicate if the stored value is a
    boolean, structure, array, or string.
    
    @param inst_info The instrument info dictionary (or structure) to write.
    @param filename The file path of the HDF5 file to be created.
    @exception Exception Raised if an existing file cannot be removed.
    """
    import os
    import h5py
    import numpy as np

    def defaultAtts(var):
        for key in ['isBoolean', 'isStruct', 'isArray', 'isString']:
            var.attrs[key] = False

    def writeVar(fp, prefix, var):
        """
        @brief Recursively write a variable to an HDF5 file.
        
        Converts the variable using built-in converters if available, then
        writes dictionaries as groups and lists as arrays with keys based on
        their index.
        
        @param fp The HDF5 file object.
        @param prefix Current file path (group) in the HDF5 file.
        @param var The variable to write.
    
        """
        ## Call built-in converters on var if present.
        for conv in ['tolist', 'to_list', 'todict', 'as_dict']:
            if hasattr(var, conv):
                conv_func = getattr(var, conv)
                var_converted = conv_func()
                if isinstance(var, np.ndarray):
                    ## Don't convert numpy arrays - should be saved as singular
                    ## object instead of recursively broken down.
                    break
                else:
                    ## Convert to list or dict.
                    var = var_converted
                    break
                
        if isinstance(var, dict):
            ## If var is a dictionary, recursively call the writing function.
            if prefix != '/':
                fp.create_group(prefix)
            defaultAtts(fp[prefix])
            fp[prefix].attrs['isStruct'] = True
            for key, val in var.items():
                writeVar(fp, prefix + '/' + key, val)
        elif isinstance(var, list):
            ## If var is homogeneous list representing an array,
            ## convert to np.array.
            try:
                arr = np.array(var)
                if arr.dtype != np.object:
                    ## homogeneous numeric or string array
                    ## Emulate file structure by creating a group and a
                    ## dataset inside it.
                    fp.create_group(prefix)
                    defaultAtts(fp[prefix])
                    fp[prefix].attrs['isArray'] = True
                    ds_name = prefix + '/data'
                    fp.create_dataset(ds_name, data=arr)
                    ## Set similar attributes on inner dataset if required.
                    defaultAtts(fp[ds_name])
                    return
            except Exception:
                ## Fallback to recursive treatment if conversion fails.
                ## Most common option.
                pass
            if prefix != '/':
                fp.create_group(prefix)
            defaultAtts(fp[prefix])
            fp[prefix].attrs['isArray'] = True
            n = np.ceil(np.log10(len(var)))
            fmt = '%s/%0' + ('%.0f' % n) + 'd'
            for i, val in enumerate(var):
                key = fmt % (prefix, i)
                writeVar(fp, key, val)
        else:
            fp[prefix] = var
            defaultAtts(fp[prefix])
            fp[prefix].attrs['isString'] = isinstance(var, str)
            fp[prefix].attrs['isBoolean'] = isinstance(var, bool)

    inst_dict = CRs_to_dicts(inst_info)

    if os.path.exists(filename):
        os.remove(filename)
    if os.path.exists(filename):
        raise Exception('Could not remove %s ' % filename + \
                        '(required to write a new hdf5 file)')
    with h5py.File(filename, 'w') as fp:
        writeVar(fp, '/', inst_dict)

        
def read_h5(filename):
    """
    @brief Read instrument info dict/structure from an HDF5 file.
    
    Reads the specified HDF5 file and reconstructs the instrument info data
    structure. This function does not convert data into ChannelResponse
    objects.
    
    @param filename The HDF5 file to read.
    @return The instrument info dictionary reconstructed from the HDF5 file.
    """
    import h5py

    def read_recursive(fp, prefix):
        """
        @brief Recursively read a variable from an HDF5 file.
        
        Reconstructs nested dictionaries or lists based on attributes stored
        in the file.
        
        @param fp The HDF5 file object.
        @param prefix The current path in the HDF5 file.
        @return The reconstructed variable.
        """
        # Only checks if attribute exists, not the value.
        if fp[prefix].attrs.get('isStruct', False):
            var = {}  # Create a dict result.
            for key in fp[prefix]:
                var[key] = read_recursive(fp, prefix + '/' + key)
        elif fp[prefix].attrs.get('isArray', False):
            keys = [k for k in fp[prefix]]
            var = []
            for key in sorted(keys):
                var.append(read_recursive(fp, prefix + '/' + key))
        else:
            ## Old version:
            #var = fp[prefix][()]
            #if fp[prefix].attrs['isBoolean']:
            #    var = bool(var)
            #elif fp[prefix].attrs['isString']:
            #    var = var.decode('utf-8')

            ## New version:
            if type(fp[prefix]) == h5py._hl.group.Group:
                var = {}
                for key in fp[prefix]:
                    var[key] = read_recursive(fp,prefix+'/'+key)
            else:   
                var = fp[prefix][()]
                if fp[prefix].attrs.get('isBoolean', False):
                    var = bool(var)
                elif fp[prefix].attrs.get('isString', False):
                    if not isinstance(var, np.ndarray):
                        var = var.decode('utf-8')
                    else:
                        var = "".join([chr(i) for i in var[:,0]])

        return var

    with h5py.File(filename, 'r') as fp:
        inst_dict = read_recursive(fp, '/')
    return inst_dict


if __name__ == '__main__':

    # Factory constructor demo.
    class A(FactoryConstructorMixin):  # Abstract base class.
        @classmethod
        def _is_mine(cls, *args, **kwargs):
            # Abstract base class never the right answer.
            raise Exception("Provided initialization data did not " + \
                            "correspond to any implemented subclass " + \
                            "of '%s'" % cls.__name__)
        def __init__(self, *args, **kwargs):
            if 'type' in kwargs:
                self.type = kwargs['type']
            else:
                self.type = None
        def __str__(self):
            return 'Object of class %s type=%s' % (self.__class__.__name__,
                                                   self.type)
            
    class B(A):
        @classmethod
        def is_mine(cls, *args, **kwargs):
            return ('type' in kwargs) and (kwargs['type'] == 'B')
    class C(A):
        @classmethod
        def is_mine(cls, *args, **kwargs):
            return ('type' in kwargs) and (kwargs['type'] == 'C')
    class D(C):  # Extends C. C and D coexist.
        @classmethod
        def is_mine(cls, *args, **kwargs):
            return ('type' in kwargs) and (kwargs['type'] == 'D')
    
    class E(C):
        pass  # Replaces C in hierarchy.
    
    try:    
        print('A()', A())
    except Exception as e:
        print(e)
    try:    
        print('A', A(type='A'))
    except Exception as e:
        print(e)
    print('B', A(type='B'))
    print('C', A(type='C'))
    print('D', A(type='D'))
    try:    
        print('E', A(type='E'))
    except Exception as e:
        print(e)
    
    # Very simple test of RFL class instantiation.
    print('Omni int', ChannelResponse(RESP_TYPE='[E]', E_TYPE='INT', E0=2))
    print('Tele_sym wide', ChannelResponse(RESP_TYPE='[E],[TH]',
                                           E_TYPE='WIDE', E0=2, E1=3,
                                           TH_TYPE='CYL_TELE', R1=1, R2=2,
                                           D=3.0))
    print('Tele_asym diff', ChannelResponse(RESP_TYPE='[E],[TH,PH]',
                                            E_TYPE='DIFF', E0=2, DE=1,
                                            TP_TYPE='RECT_TELE', W1=1.0,
                                            H1=2.0, W2=0.5, H2=1.0, D=3.0))
