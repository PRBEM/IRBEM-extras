"""
RFL Generic Telescope

A generic telescope is meant for combinations of circular and rectangular 
elements, allowed to be at offsets and tilted. However, it can handle any 
convex shape whose edges are a combination of lines and elliptical arcs

nomenclature notes:
A segment is the part of an ellipse between a chord and it's arc. 
A minor segment is the smaller part, major segment is the larger 
part which includes the center of the ellipse.

A sector is the part of an ellipse between two radii and their arc.

The area of a segment is the sector area minus the triangle formed by the 
radii and their chord. A point is inside a segment of it is inside the sector 
but not the triangle.

A point p is on an arc if it's on the ellipse and r1 x p is parallel 
to r1 x r2, and parallel to p x r2.

A patch is made up of straight edges and ellipse arcs. And it is convex. 
The overlap of two patches is defined by the intersection points of the 
patch boundaries and the vertices of either patch that are inside the other 
patch. Together these points define a new patch.

A Patch can be broken into ellipse segments and triangles. And these can be 
used for determining area and whether a point is inside the patch. 
A Patch is made up of Edges.

A generic CurvedEdge should exist which has .area and .inside methods for 
the area between the curve and its chord. Then getting an area of a patch 
with such features would be to get the area of the curved edges and then 
the area of the left over polygon.

So Polygon would be a special case of Patch with only straight edges. 
And it would consist of triangles.

EllipseArc would be a special case of CurvedEdge.

Triangle would be a special case of Polygon with only three edges. 
This is where we would implement the inside and area functions.

Edge would be a base class for StraightEdge and CurvedEdge. We would have 
intersection rules for curved and straight edges. We would also have 
rotation rules for edges.

Edge
. StraightEdge 
. CurvedEdge
..EllipseArc
...FullEllipse
....FullCircle

.intersect returns None if class cannot do intersection with other class. 
Returns list of intersecting points otherwise. Algorithm:
self knows how to do intersection? If no, check super. At Edge level, 
check other.intersect(self, False). False prevents endless recursion. 
If still no results raise error.

Patch
. Polygon
.. Triangle
. Rectangle
. Segment
. EllipsePatch
.. CirclePatch

.xy2p converts from local coords x,y to 3-d vector, e.g origin+x*xhat+y*yhat.
Triangles generate their own xy2p
Segments, Circles, and Ellipses get their xy2p from the edge

.overlap reruns a new patch that is the overlap of self with another patch.

Need to break curved edges into pieces <180 degrees around
This will occur when a patch is built from edges


module settings:
REFERENCE_LENGTH - 1.0. A "typical" length in the units being used
SMALL_FRACTION - 1e-10 - de minimis distance offset. Points closer than this
  fraction of REFERENCE_LENGTH are assumed to be the same point. 
  A value of 1e-10 means that if length units are cm, then any points 
  within a pm are "the same". Since we normally work in cm (or inches) this 
  should be sufficient. Also defines the smallest dot product considered to be 
  zero. Equivalent to two 1 cm long unit vectors whose tips differ by 1 pm
  Note: when comparing cosines to unity, comparison is |cos-1| to larger of
      SMALL_FRACTION**2/2 or 10*EPS, where EPS is machine precision. This
      is because cos(x) ~ 1-x^2/2, for "small" x

"""

# TODO: check spec to see what it says about orientation of phi relative
# to W,H in rectangular telescope case

import numpy as np
from scipy.integrate import dblquad

from rfl import AngleResponse,inherit_docstrings, \
    default_thetagrid,default_phigrid, broadcast_grids, \
    RFLError, sind, cosd

norm = np.linalg.norm # shorthand

REFERENCE_LENGTH = 1.0 # "typical" length in units being used
SMALL_FRACTION = 1e-10

DEBUG = False
DEBUG_STATE = {} # keeps debug variables to print out in exceptions

def between_rays(a,b,p,n=None):
    """
    bool = between_rays(a,b,p,n=None)
    is p between rays connecting origin to a and b
    a and b - points defining rays
    p - test point
    n - normal vector. if None, n=cross(a,b)
    """
    if n is None:
        n = np.cross(a,b)

    # check colinearity
    norma = norm(a)
    normp = norm(p)    
    SMALL_COS = max(SMALL_FRACTION**2/2,10*np.finfo(normp).eps) # when comparing cosine to unity, use this
    # when comparing sine to zero, use SMALL_FRACTION
    if DEBUG: print('1 norma',norma,'normp',normp)
    if norma < SMALL_FRACTION*REFERENCE_LENGTH: return False # not a ray
    if DEBUG: print('2 a dot p',np.dot(a,p))
    if abs(np.dot(a,p)-norma*normp) <= SMALL_COS*norma*normp: return True # p along a
    normb = norm(b)
    if DEBUG: print('3 normb',normb)
    if normb < SMALL_FRACTION*REFERENCE_LENGTH: return False # not a ray
    if DEBUG: print('4 (b dot p)/bp-1',np.dot(b,p)/normb/normp-1)
    if abs(np.dot(b,p)-normb*normp) <= SMALL_COS*normb*normp: return True # p along b
    
    # both a x p and p x b must go in same direction as n = a x b
    pn = normp*norm(n)
    if DEBUG: print('5 normp*normn',pn)
    if DEBUG: print('6 (axp).n)',np.dot(np.cross(a,p),n))
    if np.dot(np.cross(a,p),n) < -SMALL_FRACTION*norma*pn: return False
    if DEBUG: print('7 (pxb).n)',np.dot(np.cross(p,b),n))
    if np.dot(np.cross(p,b),n) < -SMALL_FRACTION*normb*pn:return False

    if DEBUG: print('8 between')
    return True

def points_equal(a,b):
    """
    bool = points_equal(a,b)
    a,b are 3-vectors
    returns true if a and b are effectively equal.
    Formally, returns ||a-b|| < REFERENCE_LENGTH*SMALL_FRACTION
    """
    return (norm(a-b) < REFERENCE_LENGTH*SMALL_FRACTION)
    
def make_axes(xyz=(None,)*3):
    """
    make unit axes
    (xhat,yhat,zhat) = make_axes((x,y,z))
    input is a tuple or other iterable of 3-vectors
    None is alos allowed in place of one or more vectors
    if all xyz are None(default): returns identity basis functions
    if one xyz is None: uses cross product to find third (right hand rule)
    if two xyz are None: generates other two axes perpendicular to non-None one
    xhat, yhat, and zhat are unit vectors
    """
    if len(xyz) != 3:
        raise ValueError('make_axes requires a 3-element input')
    ivec = []
    for i,x in enumerate(xyz):
        if x is None: continue
        if np.size(x) != 3:
            raise ValueError('All inputs to make_axes must be 3-vectors or None')
        ivec.append(i)

    if len(ivec) == 3:
        unitize = lambda x : np.array(x).ravel()/norm(x)
        return tuple(unitize(x) for x in xyz)
    if len(ivec) == 2:
        xnew = np.cross(xyz[ivec[0]],xyz[ivec[1]])
        if ivec[1]-ivec[0] == 2: # x,z provided
            xnew = -xnew # y = z cross x, not x cross z
        tmp = [xnew]*3
        for i in ivec: tmp[i] = xyz[i]
        return make_axes(tmp)
    if len(ivec) == 1:
        z = xyz[ivec[0]]
        x = np.roll(z,1)
        y = np.cross(z,np.roll(z,1)) # use roll to create non-parallel vector
        x = np.cross(y,z)
        return make_axes(np.roll([z,x,y],ivec[0],axis=0)) # use roll to get z in right spot
    # len(ivec) == 0
    return tuple(np.eye(3))
"""        
print('xyz',make_axes())
print('xyz',make_axes(([2,0,0],None,None)))
print('xtz',make_axes((None,None,[0,0,3])))
print('xyz',make_axes(([5,0,0],[0,4,0],None)))
print('xyz',make_axes(([6,0,0],None,[0,0,7])))
print('xyz',make_axes(([8,0,0],[0,9,0],[0,0,1])))
raise Exception('stop')
"""

def project(p,nhat,offset=0):
    """
    q = project(p,nhat,offset=0)
    project p into plane perpendicular to nhat
    by removing part of p along nhat
    q = p-p.nhat*nhat + offset*nhat
    accepts p as shape (3,) or (N,3)
    if offset is given as a 3-vector, then offset is taken to be a point in 
        the plane perpendicular to nhat:
        offset -> (p-offset)*offset.nhat
    """
    if p.ndim > 1:
        nhat = np.reshape(nhat,(1,3))
        if not np.isscalar(offset):
            offset = ((p-np.reshape(offset,1,3))*nhat).sum(1,keepdims=1) # (N,1) sum of product along 2nd axis
        return p-(p*nhat).sum(1,keepdims=1)*nhat + offset*nhat
    else:
        if not np.isscalar(offset):
            offset = np.dot(p-offset,nhat)
        return p-np.dot(p,nhat)*nhat + offset*nhat

def sort_points(points,nhat=None,return_index=False):
    """
    sort points such that angles increase in nhat plane
    points = sort_points(points,nhat=None,return_index=False)
    points,isort,ireverse = sort_points(points,return_index=True,...)
    nhat - ignored if sorted, otherwise used to determine how to sort points
        defines normal vector to polygon. Points will be anticlockwise about nhat
    points will be sorted in order of increasing angle around centroid (0,2pi)
    Note: identical points will be consolidated
    """
    # algorithm:
    # - replace near-duplicate points with reference to prior
    # - find centroid (center) of all unique points
    # - if nhat not provided, get it from vectors to centroid for first two unique points
    # - get angle made around nhat by each vector from centroid to all points
    # - sort/unique points by that angle
    
    points = list(points) # make mutable
    # to address near-duplicates, we replace them with the first instance. 
    # We do this to preserve book keeping for the call to np.unique
    uniques = [] # unique points previously seen
    # also compute centroid
    center = np.zeros(3) # average point / centroid
    for i,p in enumerate(points):
        if np.size(p) != 3:
            raise ValueError('points must be list of 3-d vectors')
        p = np.array(p).ravel()
        # use for-else syntax. else clause only executed if loop completes
        for q in uniques:
            if points_equal(p,q):
                # point is nearly identical to a prior point                    
                p = q # replace it with a reference to the prior point
                break
        else: # only get here if no break in loop
            uniques.append(p) # point is unique
            center += p
        points[i] = p # overwrite w/ cleaned up point
    if len(uniques) == 1:
        if return_index:
            return uniques,np.array(0,dtype=int),np.zeros(len(points),dtype=int)
        else:
            return uniques
    center = center/len(uniques) # average
    # find nhat if needed
    zhat = nhat
    xhat = None
    if zhat is None:
        zhat = 0
        for i,p in enumerate(uniques):
            zhat = np.cross(p-center,uniques[(i+1) % len(uniques)]-center) # p x next
            if norm(zhat)>SMALL_FRACTION*REFERENCE_LENGTH: 
                break # good zhat found
        else:  # loop completed, no zhat found
            # points are colinear
            # find an arbitrary direction to use as zhat
            for p in uniques:
                if points_equal(p,center): continue # p==center
                dp = p-center
                (xhat,yhat,zhat) = make_axes((dp,np.roll(dp,1),None))
                if norm(zhat)>SMALL_FRACTION*REFERENCE_LENGTH: 
                    break # suitable zhat found
            else:
                print(DEBUG_STATE)
                raise ValueError('Unable to create normal vector (zhat)')

    if xhat is None:
        # generate x,y,nhat system
        (zhat,xhat,yhat) = make_axes((zhat,uniques[0]-center,None))
    # compute angles about nhat relative to centroid
    angles = []
    for i,p in enumerate(points): # get angles for all points
        dp = p-center
        x = np.dot(dp,xhat)
        y = np.dot(dp,yhat)
        z = np.dot(dp,zhat)
        if abs(z) > SMALL_FRACTION*max(REFERENCE_LENGTH,norm(dp)):
            print(DEBUG_STATE)
            raise ValueError('points not coplanar')
        points[i] = p = center + x*xhat+y*yhat # ensure coplanar
        angle = np.arctan2(y,x) % (2*np.pi) # put on 0,2*pi
        angles.append(angle)
    # now do sort an unique on angles
    angles,isort,ireverse= np.unique(angles,return_index=True,return_inverse=True) # remove duplicates
    # get corresponding unique points in angle order
    points = [points[i] for i in isort]
    if return_index:
        return points,isort,ireverse
    else:
        return points

# classes
    
class LocalCoordsMixin(object):
    """
    LocalCoordsMixin - mixin to provide local coordinate system info
    properties:
        origin - 3-d origin of system
        xhat,yhat,zhat - 3-d unit vectors that define system
        (x1,x2) = xlim - lower and upper limit of x in object's coordinate system
        x1 - lower limit of x in object's coordinate system
        x2 - upper limit of x in object's coordinate system
    Properties are read only. Managed by private _ variables that are 
        initilized to None in the Edge base class constructor. 
        Only getters are provided.
    methods:
        (y1,y2) = ylim(x) - y limits in object's coordinate sytem at x
        y = y1(x) - lower limit of y in object's coordinate sytem at x
        y = y2(x) - upper limit of y in object's coordinate sytem at x
        p = xy2p(x,y) convert x,y from object's coordinates to point in 3-d frame
        (x,y) = p2xy(p) convert point from 3-d frame to object's local coordinate system
        xy2p and p2xy will accept scalar or 1-d input, adding the extra dimension to
        the output of p,x,y. So, e.g., if x,y are (N,) then p will be (N,3), but 
        if x,y are scalar, then p will be(3,)
        
    Derived class must populate _origin, _xhat, _yhat, _zhat, and _xlim
        and overload __init__ and ylim
    
    """
    def __init__(self,*args,**kwargs):
        """
        Initialzies private underscore variables to None
        """
        self._origin = self._xhat = self._yhat = self._zhat = self._xlim = None
    @property
    def origin(self):
        return self._origin
    @property
    def xhat(self):
        return self._xhat
    @property
    def yhat(self):
        return self._yhat
    @property
    def zhat(self):
        return self._zhat
    @property
    def xlim(self):
        return self._xlim
    @property
    def x1(self):
        return self._xlim[0]
    @property
    def x2(self):
        return self._xlim[1]
    def ylim(self,x):
        """(y1,y2) = ylim(x) - y limits in object's coordinate sytem at x"""
        raise NotImplementedError('Method ylim must be overloaded in derived classes')
    def y1(self,x):
        """y = y1(x) - lower limit of y in object's coordinate sytem at x"""
        return self.ylim(x)[0]
    def y2(self,x):
        """y = y2(x) - upper limit of y in object's coordinate sytem at x"""
        return self.ylim(x)[1]
    def xy2p(self,x,y):
        """p = xy2p(x,y) convert x,y from object's coordinates to point in 3-d frame
        if x or y is (N,) then output will be (N,3), otherwise output is (3,)
        """
        isscalar = np.isscalar(x) and np.isscalar(y)
        if isscalar:
            return self.origin+x*self.xhat+y*self.yhat
        else:
            # reshape x,y to (N,1) or (1,1)
            x = np.reshape(np.array(x),(np.size(x),1)) # (N,1) or (1,1)
            y = np.reshape(np.array(y),(np.size(y),1)) # (N,1) or (1,1)
            if (x.size>1) and (y.size>1) and (y.size != x.size):
                raise ValueError('x and y must be scalar or same length')
            # reshape origin, xhat, yhat to (1,3)
            origin = np.expand_dims(self.origin,0) # (3,) -> (1,3)
            xhat= np.expand_dims(self.xhat,0) # (3,) -> (1,3)
            yhat= np.expand_dims(self.yhat,0) # (3,) -> (1,3)
            # now do the formula and let it broadast
            return origin+x*xhat+y*yhat # broadcasts to = (N,3)
            
    def p2xy(self,p):
        """(x,y) = p2xy(p) convert point from 3-d frame to object's local coordinate system
        if p is (N,3) x and y will be (N,). Otherwise, x,y will be scalar
        """
        dp = p-self.origin
        return (np.dot(dp,self.xhat),np.dot(dp,self.yhat)) # self-broadcasting
            

class Edge(object):
    """
    Abstract Base Class for edges
    properties:
        points = identity - tuple of points that uniquely identify a class instance
            Note that order matters. So, e.g., a straight edge from a to b is not the
            same as a straight edge from b to a. Compare to edge.reverse to check for that.
        a,b = endpoints returns Edges' end points
            or None for closed edges (e.g., full ellipses)
        len = length - length of edge
        area - area of edge for curved edes, zero for straight
        nhat - unit vector normal to Edge's plane
        pieces - tuple of Edges that comprise this edge (used to break up curved edges into <180 degree chunks)
        reversed - edge with endpoints reversed
        Properties are read only. Managed by private _ variables that are 
            initilized to None in the Edge base class constructor. 
            Only getters are provided.
    methods:
        edge = from_identity(identity) (class method) returns new object based on identity vectors
        points = intersect(other) returns a tuple of points of intersection between self and other
        edge = project(nhat,offset) returns the edge projected along nhat, with offset
        bool = inside(point) returns True if a point is inside the edge
            Always False for straight edges
            True for curved eges if point lies between curve and its chord
            also available as "in" operator
        edge = cut(p1,p2) - return edge cut down to have endpoints p1,p2
        h = plot(ax,...) - plot on axes
    overloaded operators:
        in (__contains__) - "p in edge" is the same as edge.inside(p)
        == (__eq__) - chesk that both edges are same class and 
            have same identity property, using points_equal function to test
    """
    def __init__(self,*args,**kwargs):
        """edge = Edge()
        initializes underscore cache variables for all properties
        """
        self._identity = ()
        self._endpoints = ()
        self._length = None
        self._area = None
        self._nhat = None
        self._pieces = None
        self._reversed = None
    @property
    def identity(self):
        """points = identity returns tuple of points that uniquely identify Edge"""
        return self._identity
    @classmethod
    def from_identity(cls,identity):
        """edge = from_identity(identity) (class method) returns new object based on identity vectors"""
        raise NotImplementedError('Method from_identity must be overloaded in derived classes')
    def __eq__(self,other):
        if other is self: return True
        if self.__class__ != other.__class__: return False
        for s,o in zip(self.identity,other.identity):
            if not points_equal(s,o): return False
        return True
    @property
    def endpoints(self):
        """a,b = endpoints returns end points of Edge"""
        return self._endpoints
    @property
    def length(self):
        """length of Edge"""
        return self._length
    @property
    def area(self):
        """area of Edge"""
        return self._area
    @property
    def nhat(self):
        """nhat - unit vector normal to Edge's plane"""
        return self._nhat
    @property
    def reversed(self):
        """return edge flipped so that endpoints are reversed (but still represents the same set of points)"""
        raise NotImplementedError('Method reversed must be overloaded in derived classes')
    @property
    def pieces(self):
        """return list of edges that make up this one. Used to break curved edges into <180 degree chunks"""
        # for most things, just retuns [self] so that's OK
        if self._pieces is None:
            self._pieces = (self,)
        return self._pieces
    def intersect(self,other,try_other=True):
        """points = intersect(other,try_other=True) 
        returns a typle of points of intersection between self and Edge other
        try_other: try other's intersect method if this one doesn't support other's type
        """
        # we got here because self did not implement intersect with other        
        if try_other:
            return other.intersect(self,try_other=False) # let other try
        else: # oops, coding error, neither class has intersect for the other
            raise NotImplementedError('No intersection method found for %s,%s' %
                                      (self.__class__.__name__,other.__class__.__name__))
    def project(self,nhat,offset=0):
        """
        edge = project(nhat,offset=0) returns the edge projected into nhat plane with offset
        see function project for how points are projected.
        """
        return self.from_identity(project(np.array(self.identity),nhat,offset))
        
    def inside(self,point):
        """
        bool = inside(point) returns True if a point is inside the edge
        can also accept list of points and return list of bools
        """
        raise NotImplementedError('Method inside must be overloaded in derived classes')
    def cut(self,p1,p2):
        """edge = cut(p1,p2) - return edge cut down to have endpoints p1,p2"""
        raise NotImplementedError('Method cut must be overloaded in derived classes')
    def __contains__(self,point): # in operator calls inside
        return self.inside(point)
    def plot(self,ax,*args,**kwargs):
        """
        h = plot(self,ax,*args,**kwargs)
        call ax.plot on axis with args and kwargs
        returns result of ax.plot
        """
        raise NotImplementedError('Method plot must be overloaded in derived classes')
        
    
class StraightEdge(Edge):
    """StraightEdge implements Edge for edges that are line segments"""
    def __init__(self,p1,p2,*args,**kwargs):
        """straight = StraightEdge(p1,p2,...)
        p1 - 3-d location of start of line segment
        p2 - 3-d location of end of line segment
        """
        super().__init__(*args,**kwargs)
        self._area = 0.0
        for p in [p1,p2]:
            if np.size(p) != 3:
                raise ValueError('Inputs p1 and p2 to StraightEdge must be 3-vectors')
        p1 = np.array(p1).ravel()
        p2 = np.array(p2).ravel()
        self._endpoints = (p1,p2)
        self._identity = self._endpoints
    @property
    def length(self):
        """*INHERIT*"""
        if self._length is None:
            self._length = norm(self._endpoints[1]-self._endpoints[0])
        return self._length
    @property
    def reversed(self):
        """*INHERIT*"""
        if self._reversed is None:
            self._reversed = self.__class__(self._endpoints[1],self._endpoints[0])
        return self._reversed
    @classmethod
    def from_identity(cls,identity):
        """*INHERIT*"""
        return cls(*identity) # identity is endpoints
    def _intersect_with_straight(self,other):
        """
        same as intersect, but when other is known to be a StraightEdge
        if edges are co-linear, then up to two points will be returned
            indicating the endpoints of the overlap
        """
        empty = tuple() # no intersection
        # self is: p(t) = p1+(p2-p1)*t = p1 + dp*t
        # other is: q(s) = q1+(q2-q1)*s = q1 + dq*s
        # find t,s such that p(t) = q(s)
        
        (p1,p2) = self.endpoints
        (q1,q2) = other.endpoints
        dp = p2-p1
        dq = q2-q1
        normdp = norm(dp)
        normdq = norm(dq)
        SMALL_COS = max(SMALL_FRACTION**2/2,10*np.finfo(normdp).eps) # when comparing cosine to unity, use this
        # when comparing sine to zero, use SMALL_FRACTION
        
        # if dp is parallel to dq only colinear case possibly applies
        if (abs(abs(np.dot(dp,dq))-normdp*normdq) <= SMALL_COS*normdp*normdq):
            # p1->p2 is (anti)parallel to q1->q1
            # check colinear case: intersection is endpoints of overlap
            # is q1 on p1->p2?
            dq1 = q1-p1
            normdq1 = norm(dq1)
            dq1dp = np.dot(dq1,dp)
            if (abs(abs(dq1dp)-normdq1*normdp) <= SMALL_COS*normdq1*normdp):
                # q1 is on pq->p2, colinear case
                t = [0.0,1.0] # p1,p2
                t.append(dq1dp/normdp**2) # q1
                t.append(np.dot(q2-p1,dp)/normdp**2) # q2
                t = np.sort(t) #  middle two points [1,2] are intersection
                if (t[1] < -SMALL_FRACTION) or (t[2]>1.0+SMALL_FRACTION):
                    return empty # line segments don't overlap
                if abs(t[1]-t[2]) < SMALL_FRACTION:
                    return (p1+dp*t[1],) # line segments are contiguous
                return (p1+dp*t[1],p1+dp*t[2]) # return middle two points 
            else:
                return empty # parallel lines that are not colinear do not cross

        # possibly crossing lines case        
        # solve A*x = y
        # A = [p2-p1,q2-q1] (stack as two rows)
        # y = q1-p1
        # x = (t,s)
        # use least squares since A is overspecified
        A = np.column_stack([dp,-dq])
        y = q1-p1
        (t,s) = np.linalg.lstsq(A,y,rcond=None)[0] # 1st output is solution
        if (t<-SMALL_FRACTION) or (t>1.0+SMALL_FRACTION):
            return empty # intersection is not between p1,p2
        if (s<-SMALL_FRACTION) or (s>1.0+SMALL_FRACTION):
            return empty # intersection is not between q1,q2
        if not points_equal(p1+dp*t,q1+dq*s):
            return empty # points do not intersect
        return (p1+dp*t,) # return tuple of single intersecting point        
    def intersect(self,other,try_other=True):
        """*INHERIT*"""
        # try intersect with types self knows about
        # then hand off to parent class
        if self == other:
            return self.endpoints # these will generally be trimmed
        if isinstance(other,StraightEdge):
            return self._intersect_with_straight(other)
        else:
            return super().intersect(other,try_other) # pass up the chain
    def inside(self,point):
        """*INHERIT*
        always false for straight edges
        """
        return False
    def cut(self,p1,p2):
        """*INHERIT*"""
        return StraightEdge(p1,p2)
    def plot(self,ax,*args,**kwargs):
        """*INHERIT*"""
        # make the colors match
        h = ax.plot(*np.array(self.endpoints).transpose(),*args,**kwargs)
        for hi in h[1:]:
            hi.set_color(h[0].get_color())
        return h

class CurvedEdge(Edge,LocalCoordsMixin):
    """
    CurvedEdge base subclass of Edge for curved edges
    Also implements LocalCoordsMixin to represent local coordinate system
    """
    @property
    def area(self):
        """*INHERIT*
        area between curve and its chord
        """
        raise NotImplementedError('Method area must be overloaded in derived classes')
    @property
    def reversed(self):
        """*INHERIT*"""
        raise NotImplementedError('Method reversed must be overloaded in derived classes')

class EllipseArc(CurvedEdge):
    """EllipseArc Edge, subclass of CurvedEdge, for elliptical arcs
    in addition to Edge, adds the following:
    properties:
        center - 3-d location of circle's center
        a - semimajor axis 3-d vector
        b - semiminor axis 3-d vector
        angles = tuple of start and end angle, always angles[1]>angles[0]
        extent = angular extend in degrees = angles[1]-angles[0]
        ra - length of semimajor axis
        rb - length of semiminor axis
    methods:
        theta2p(theta) = return 3-d point given angle in degrees
    """
    def __init__(self,center,r1,r2,p1,p2,*args,**kwargs):
        """ellpiseArc = Ellipse(center,r1,r2,p1,p2)
        center- 3-d location of ellipse's center
        r1 - 3-d location of point on ellipse
        r2 - 3-d location of point on ellipse along conjugated diamater to r1
        (if r1 perpendicular to r2, then they are the semimajor and semiminnor axes)
        p1 - theta or 3-d location of first point on arc
        p2 - theta or 3-d location of second point on arc
        if p1 and p2 are given as angles then they are taken to be degrees, such that:
            p = r1*cos(theta) + r2*sin(theta)
        arc spans from p1 to p2 anticlockwise (around r1 x r2)
        """
        # stub for testing. Does not compute semiaxes
        for p in [center,r1,r2]:
            if np.size(p) != 3:
                raise ValueError('Inputs center,r1,r2 to Ellipse must be 3-vectors')

        center = np.array(center).ravel()
        r1 = np.array(r1).ravel()
        r2 = np.array(r2).ravel()
        if np.isscalar(p1): # theta in degrees
            p1 = center+r1*cosd(p1) + r2*sind(p1)
        else:
            if np.size(p1) != 3:
                raise ValueError('Input p1 to Ellipse must be 3-vector')
            p1 = np.array(p1)
        if np.isscalar(p2): # theta in degrees
            p2 = center+r1*cosd(p2) + r2*sind(p2)
        else:
            if np.size(p2) != 3:
                raise ValueError('Input p2 to Ellipse must be 3-vector')
            p2 = np.array(p2)
        super().__init__()
        self._center = center
        self._endpoints = (p1,p2)

        # build semiaxes from r1, r2
        # https://en.wikipedia.org/wiki/Rytz%27s_construction#Computer_aided_solution
        
        t0 = np.arctan2(2*np.dot(r1,r2),(norm(r1)**2-norm(r2**2)))
        ax1 = center + r1*np.cos(t0) + r2*np.cos(t0)
        ax2 = center + r1*np.cos(t0+np.pi/2) + r2*np.cos(t0+np.pi/2)
        if norm(ax2) > norm(ax1): # first axis is semiminor
            ax1 = ax2
            ax2 = center + r1*np.cos(t0+np.pi) + r2*np.cos(t0+np.pi)
        self._a = ax1-center
        self._b = ax2-center
        self._ra = None
        self._rb = None
        angles = [np.degrees(np.arctan2(np.dot(p-center,self.b),np.dot(p-center,self.a))) for p in self._endpoints]
        if angles[1] <= angles[0]:
            angles[1] += 360 # ensure angles[1] > angles[0]
        self._angles = tuple(angles)
        self._identity = (self.center,self.a,self.b,*self._endpoints)
        
        # set up local coords mixin
        self._origin = self._center
        self._xhat = self._a/self.ra
        self._yhat = self._b/self.rb
        self._zhat = np.cross(self._xhat,self._yhat)
    @classmethod
    def from_identity(cls,identity):
        """*INHERIT*"""
        return cls(*identity)
    @property
    def center(self):
        return self._center
    @property
    def a(self):
        return self._a
    @property
    def b(self):
        return self._b
    @property
    def ra(self):
        if self._ra is None:
            self._ra = norm(self.a)
        return self._ra
    @property
    def rb(self):
        if self._rb is None:
            self._rb = norm(self.b)
        return self._rb
    @property
    def angles(self):
        """*INHERIT*"""
        return self._angles
    @property
    def extent(self):
        """*INHERIT*"""
        return self.angles[1]-self.angles[0]
    @property
    def length(self):
        """
        *INHERIT*
        length is computed along arc
        """
        # https://math.stackexchange.com/questions/433094/how-to-determine-the-arc-length-of-ellipse
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipeinc.html#scipy.special.ellipeinc
        raise NotImplementedError('TODO') # TODO
    @property
    def area(self):
        """*INHERIT*"""
        raise NotImplementedError('TODO') # TODO
        # area of elipse arc is area between arc and chord
        # so, area of pie slice, minus area of inscribed trangle
    @property
    def pieces(self):
        """*INHERIT*"""
        if self._pieces is None:
            N = np.floor(self.extent/180)+1 # number of pieces
            if N == 1:
                self._pieces = (self,)
            else: # break up
                theta = np.linspace(self.angles[0],self.angles[1],N+1)
                points = self.theta2p(theta)
                pieces = []
                for i in range(N):
                    p1 = points[i]
                    p2 = points[i+1]
                    pieces.append(self.__class__(self.center,self.center+self.a,self.center+self.b,p1,p2))
                self._pieces = tuple(pieces)
                
        return self._pieces
    @property
    def reversed(self):
        """*INHERIT*"""
        # reverse n=axb by swapping -b
        # reverse endpoints by swapping endpoints
        # yields same segment, but on flipped elipse
        return self.__class__(self.center,self.center+self.a,self.center-self.b,self.endpoints[1],self.endpoints[0])
    @property
    def xlim(self):
        """*INHERIT"""
        # xhat points along a
        if self._xlim is None:
            xlim = sorted([np.dot(self.endpoints[0]-self.center,self.a),np.dot(self.endpoints[1]-self.center,self.a)])
            if (self.angles[0] <= 0 <= self.angles[1]) or (self.angles[0] <= 360 <= self.angles[1]):
                xlim.append(self.ra) # semimajor axis is in elipse
            if (self.angles[0] <= 180 <= self.angles[1]) or (self.angles[0] <= 180+360 <= self.angles[1]):
                xlim.append(-self.ra) # anti-semimajor axis is in elipse
            xlim = sorted(xlim)
            self._xlim = (xlim[0],xlim[-1])
        return self._xlim            
    def ylim(self,x):
        """*INHERIT"""
        # yhat points along b
        # (x/a)^2 + (y/b)^2 = 1
        # y = +/- b*sqrt(1-(x/a)^2)
        if (x <= -self.ra) or (x >= self.ra):
            return (0.0,0.0)
        y = self.rb*np.sqrt(1-(x/self.ra)**2)
        ylim = [0.0,0.0]
        # check -b side
        theta = np.degrees(np.arctan2(-y,x))
        if (self.angles[0] <= theta <= self.angles[1]) or (self.angles[0] <= theta+360 <= self.angles[1]):
            ylim[0] = -y
        # check +b side
        theta = np.degrees(np.arctan2(y,x))
        if (self.angles[0] <= theta <= self.angles[1]) or (self.angles[0] <= theta+360 <= self.angles[1]):
            ylim[1] = y
        return tuple(ylim)
    
    
    def theta2p(self,theta):
        """p = theta2p(theta)
           convert angle in degrees to 3-d point
           theta - angle or list of angles in degrees
           p - (3,) 3-d point or list of 3-d points if theta is list
        """
        scalar = np.isscalar(theta)
        if scalar:
            theta = [theta]
        points = [None]*len(theta)
        for i,th in enumerate(theta):
            points[i] = self.center + self.a*cosd(th) + self.b*sind(th)
        if scalar:
            return points[0]
        else:
            return points
    def intersect(self,other,try_other=True):
        """*INHERIT*"""
        # try intersect with types self knows about
        # then hand off to parent class
        if isinstance(other,EllipseArc):
            raise NotImplementedError('TODO') # TODO
        elif isinstance(other,StraightEdge):
            raise NotImplementedError('TODO') # TODO
        else:
            super().intersect(other,try_other) # pass up the chain
    def inside(self,point):
        """*INHERIT*"""
        raise NotImplementedError('TODO') # TODO
    def cut(self,p1,p2):
        """*INHERIT*"""
        if (p1 not in self) or (p2 not in self):
            raise ValueError('Cannot cut to points outside %s' % self.__class__.__name__)
        return EllipseArc(self.center,self.a,self.b,p1,p2) # explicit class b/c inherited by FullEllipse

class FullEllipse(EllipseArc):
    """FullEllipse Edge, subclass of EllipseArc, full ellipse"""
    def __init__(self,center,r1,r2,*args,**kwargs):
        """fullEllpise = FullEllipse(center,r1,r2)
        center- 3-d location of ellipse's center
        r1 - 3-d location of point on ellipse
        r2 - 3-d location of point on ellipse along conjugated diamater to r1
        (if r1 perpendicular to r2, then they are the semimajor and semiminnor axes)
        """
        super().__init__(center,r1,r2,0.0,360.0)
        self._identity = (self.center,self.a,self.b)
    @classmethod
    def from_identity(cls,identity):
        """*INHERIT*"""
        return cls(*identity) # identity is endpoints        
    @property
    def area(self):
        """*INHERIT*
        area inside closed edge
        """
        return np.pi*self.ra*self.rb
    @property
    def xlim(self):
        return (-self.ra,self.ra)
    def ylim(self,x):
        """*INHERIT"""
        if (x <= -self.ra) or (x >= self.ra):
            return (0.0,0.0)
        y = self.rb*np.sqrt(1-(x/self.ra)**2)
        return (-y,y)
    def inside(self,point):
        """*INHERIT*"""
        (x,y) = self.p2xy(point)
        return (x/self.ra)**2 + (y/self.rb)**2 <= 1+SMALL_FRACTION

class CircleArc(EllipseArc):
    """CircleArc Edge, subclass of EllipseArc, for circular arcs"""
    def __init__(self,center,p1,p2):
        """arc = CircleArc(center,p1,p2)
        center- 3-d location of circle's center
        p1 - 3-d location of first point on arc
        p2 - 3-d location of second point on arc
        arc spans from p1 to p2 anticlockwise (around p1 x p2)
        """
        for p in [center,p1,p2]:
            if np.size(p) != 3:
                raise ValueError('Inputs center,p1,p2 to CircleArc must be 3-d vectors')
        center = np.array(center).ravel()
        p1 = np.array(p1).ravel()
        p2 = np.array(p2).ravel()
        r1 = p1-center
        r = norm(r1)
        r2 = p2-center
        if np.abs(r-norm(r2)) > REFERENCE_LENGTH*SMALL_FRACTION:
            raise ValueError('points p1 and p2 are not equidistant from center')
        nhat = np.cross(r1,r2)
        r2 = np.cross(nhat,r1)
        r2 = r2/norm(r2)*r
        super().__init__(center,center+r1,center+r2,p1,p2)
    @classmethod
    def from_identity(cls,identity):
        """*INHERIT*"""
        center,a,b,p1,p2 = identity
        return cls(center,p1,p2)
    def cut(self,p1,p2):
        """*INHERIT*"""
        if (p1 not in self) or (p2 not in self):
            raise ValueError('Cannot cut to points outside %s' % self.__class__.__name__)
        return self.__class__(self.center,p1,p2)
    @property
    def reversed(self):
        """*INHERIT*"""
        return self.__class__(self.center,p2,p1)

class FullCircle(FullEllipse):
    """FullCircle Edge, subclass of Ellipse, full circle
    In addition to FullEllipse, has the following:
    properties:
        center - 3-d location of circle's center
        r - radius of circle (scalar)
        nhat - 3-vector normal to plane of circle
    """
    def __init__(self,center,r,nhat):
        """circle = FullCircle(center,r)
        center - 3-d location of circle's center
        r - radius of circle (scalar)
        nhat - 3-vector normal to plane of circle
        """
        for p in [center,nhat]:
            if np.size(p) != 3:
                raise ValueError('Inputs center, nhat to FullCircle must be 3-vectors')
        if not np.isscalar(r):
            raise ValueError('Input r to FullCircle must be scalar')
        center = np.array(center).ravel()
        nhat = np.array(nhat).ravel()
        nhat = nhat/norm(nhat) # just to be sure
        self.nhat = nhat
        self.r = r
        p1 = np.cross(nhat,np.roll(nhat))
        p1 = p1/norm(p1)
        p2 = np.cross(nhat,p1)
        super.__init__(center,center+r*p1,center+r*p2)
    @classmethod
    def from_identity(cls,identity):
        """*INHERIT*"""
        center,a,b = identity
        nhat = np.cross(a,b)
        nhat = nhat/norm(nhat)
        r = norm(a)        
        return cls(center,r,nhat)
    @property
    def reversed(self):
        """*INHERIT*"""
        return self.__class__(self.center,self.r,-self.nhat) # same circle opposite sense
    def cut(self,p1,p2):
        """*INHERIT*"""
        if (p1 not in self) or (p2 not in self):
            raise ValueError('Cannot cut to points outside %s' % self.__class__.__name__)
        return CircleArc(self.center,p1,p2)


inherit_docstrings(Edge)
# end of Edges

class Patch(object):
    """
    Convex Patch object consisting of straight and curved edges
    convex patches only
    properties:
        edges = tuple of Edges that make up the Patch. May not be same as 
            those provided at initialization due to "pieces" breakup
        corners - tuple of points that make up the Patch (edge endpoints)
          edges and corners are ordered anticlockwise
        area - area of patch
        origin - vector position of reference point in plane of patch
        nhat - unit vector normal to Pacth's plane
        components - constituent component Patches
        chord_poly - a polygon made up of the endpoints for all edges (may be None)
        centroid - average location of corners
        Properties are read only. Managed by private _ variables and getters
    methods:
        bool = inside(point) returns True if a point is inside the Patch
            also available as "in" operator
        patch = project(nhat,offset=0) returns the patch projected along nhat, with offset
        h = plot(ax,...) - plot on axes
        patch = overlap(other) returns patch representing overlap of self and other
        from_edges(edges) - class method to create new object given edges
        
    Note: only component Patches (Triangles and Curved Edges) implement xlim, ylim, and xy2p
    """
    def __init__(self,edges):
        """
        patch = Patch(edges)
        edges - array of Edge objects
        """
        
        for edge in edges:
            if not isinstance(edge,Edge):
                raise ValueError('Patch should be initialized with list of Edge objects')

        # break up >=180 degree edges if needed
        pieces = []
        for edge in edges:
            pieces += edge.pieces

        self._edges = edges = tuple(pieces)
        
        corners = [*edges[0].endpoints]

        for edge in self._edges[1:]:
            if not points_equal(edge.endpoints[0],corners[-1]):
                raise ValueError('Supplied edges do not connect')
            corners.append(edge.endpoints[1])
        if points_equal(corners[0],corners[-1]):
            corners.pop() # remove duplicate 1st/last point
        else:
            raise ValueError('Supplied edges do not connect')
        self._corners = tuple(corners) # corners (Edge endpoints) that make up this patch        
        # initialize cached properties
        self._chord_poly = None
        self._components = None
        self._area = None
        self._centroid = None
        self._nhat = None
    @classmethod
    def from_edges(cls,edges):
        """
        obj = from_edges(edges)
        create new object given edges
        """
        return cls(edges)
    def __str__(self):
        return 'Patch %s' % (','.join([e.__class__.__name__ for e in self.edges]))
    @property
    def edges(self):
        """a list of Edges that make up the Patch"""
        return self._edges
    @property
    def corners(self):
        """a list of points that make up the Patch (from each Edge.endpoints)"""
        return self._corners
    @property
    def chord_poly(self):
        """a Polygon made from the Patch's corners. May be None"""
        if (self._chord_poly is None) and (len(self._corners)>2):
            self._chord_poly = Polygon(self._corners)
        return self._chord_poly
    @property
    def centroid(self):
        """average of corners"""
        if self._centroid is None:
            center = np.zeros(3)
            for p in self.corners:
                center += p
            self._centroid = center/len(self.corners)
        return self._centroid
                
    @property
    def area(self):
        """area of Patch"""
        if self._area is None:
            A = 0.0
            for c in self.components:
                A += c.area
            self._area = A
        else:
            A = self._area
        return A
    @property
    def origin(self):
        """reference point in Patch's plane"""
        if isinstance(self,LocalCoordsMixin):
            return super().origin
        else:
            return self.components[0].origin
    @property
    def nhat(self):
        """unit vector normal to Patch's plane"""
        return self.components[0].nhat
    @property
    def components(self):
        """constinuent components (Triangles and CurvedEdges)"""
        if self._components is None:
            components = []
            origin = self.edges[0].endpoints[0]
            for i,edge in enumerate(self.edges):
                add_triangle = False
                if isinstance(edge,CurvedEdge):
                    components.append(edge)
                elif isinstance(edge,StraightEdge):
                    # don't add triangle if first edge is straight
                    add_triangle = (i not in [0,len(self.edges)-1])
                else:
                    raise ValueError("Patch has edges that's neither straight nor curved")
                if add_triangle:
                    components.append(Triangle((origin,*edge.endpoints)))
            self._components = tuple(components)
        return self._components
    def inside(self,point):
        """
        bool = inside(point) returns True if a point is inside the patch
        can also accept list of points and return list of bools
        """
        for c in [*self.edges,*self.components]:
            if point in c:
                return True
        return False
    def __contains__(self,point): # in operator calls inside
        return self.inside(point)
    def project(self,nhat,offset=0):
        """
        patch = project(nhat,offset=0) returns the patch projected along nhat, with offset
        see project function definition for how offset works
        """
        edges = []
        for edge in self.edges:
            edges.append(edge.project(nhat,offset))
        return self.from_edges(edges)
    def plot(self,ax,*args,**kwargs):
        """
        h = plot(self,ax,*args,**kwargs)
        plots all edges, returns concatenated results as list
        """
        h = []
        for edge in self.edges:
            h = h+edge.plot(ax,*args,**kwargs)
        # make the colors match
        for hi in h[1:]:
            hi.set_color(h[0].get_color())
        return h
    def overlap(self,other):
        """
        patch = overlap(other) returns patch representing overlap of 
            self and other, or None if no overlap
        """
        # will need to overload this in CurvedEdge and in Triangle
        intersections = [] # entries are [point,[edge,...]]
        # add all intersections

        for e1 in self.edges:
            for e2 in other.edges:
                intersect = e1.intersect(e2)
                for i in intersect:
                    intersections.append([i,[e1,e2]])
        # add any endpoints that are inside other
        for (slf,oth) in [(self,other),(other,self)]:
            for i,e in enumerate(slf.edges):
                p = e.endpoints[1] # only need to check 2nd endpoint for each edge
                if p in oth:
                    e2 = slf.edges[(i+1)%len(slf.edges)] # next edge starts with p
                    intersections.append([p,[e,e2]])

        if len(intersections)>0:
            # get unique,sorted points
            points,isort,ireverse = sort_points([x[0] for x in intersections],nhat=self.nhat,return_index=True)
        else:
            isort = ireverse = points = []

        if len(points)<=1:
            if self.centroid in other:
                return self # self entirely in other
            elif other.centroid in self:
                return other # other entirely in self
            else:
                return None # overlap is a line or a point

        # assemble unique intersections
        uix = [] # unique intersections
        for i,j in enumerate(isort):
            p = intersections[j][0]
            edges = [] # accumulate unique edges associated with p
            for k in np.where(ireverse==i)[0]:
                for edge in intersections[k][1]:
                    if edge not in edges:
                        edges.append(edge) # combine edges tied to duplicate points
            uix.append([p,edges])
        
        N = len(uix)
        # find edges that span unique intersections
        edges = []
        for i in range(N):
            j = (i+1)%N
            p1 = uix[i][0]
            p2 = uix[j][0]
            # find edge both points have in common
            edges1 = uix[i][1]
            edges2 = uix[j][1]
            common_edge = None
            for edge in edges1:
                if edge in edges2:
                    common_edge = edge;
                    break
            if common_edge is None:
                
                # make a 2-d plot we can zoom in on
                plt.figure()
                xhat = self.edges[0].endpoints[1]-self.edges[0].endpoints[0]
                (zhat,xhat,yhat) = make_axes((self.nhat,xhat,None))
                for edge in self.edges+other.edges:
                    x1 = np.dot(edge.endpoints[0],xhat)
                    y1 = np.dot(edge.endpoints[0],yhat)
                    x2 = np.dot(edge.endpoints[1],xhat)
                    y2 = np.dot(edge.endpoints[1],yhat)
                    plt.plot([x1,x2],[y1,y2],'.-')
                plt.xlabel('X')
                plt.ylabel('Y')
                

                # intersect every edge with every other edge
                for (i1,e1) in enumerate(self.edges):
                    for (i2,e2) in enumerate(other.edges):
                        ax = plt.figure().add_subplot(projection='3d')
                        self.plot(ax,'ks-')
                        other.plot(ax,'gs-')
                        e1.plot(ax,'r-',lw=2)
                        e2.plot(ax,'b:',lw=2)
                        intersect = e1.intersect(e2)
                        for p in intersect:
                            ax.plot([p[0]],[p[1]],[p[2]],'ro')
                        ax.set_title('edge crossing check %dx%d = %d points' % (i1,i2,len(intersect)))
                # add any endpoints that are inside other
                for (slf,oth) in [(self,other),(other,self)]:
                    ax = plt.figure().add_subplot(projection='3d')
                    slf.plot(ax,'ks-')
                    oth.plot(ax,'gs-')
                    for i,e in enumerate(slf.edges):
                        p = e.endpoints[1] # only need to check 2nd endpoint for each edge
                        if p in oth:
                            ax.plot([p[0]],[p[1]],[p[2]],'ro')
                        else:
                            ax.plot([p[0]],[p[1]],[p[2]],'rx',lw=2)                            
                        ax.text(p[0],p[1],p[2],str(i),horizontalalignment='left')
                    ax.set_title('inside corner check')


                ax = plt.figure().add_subplot(projection='3d')
                self.plot(ax,'ks-')
                other.plot(ax,'gs-')
                ax.plot([p1[0]],[p1[1]],[p1[2]],'ro')
                ax.plot([p2[0]],[p2[1]],[p2[2]],'bo')
                for edge in edges1: edge.plot(ax,'r-',lw=2)
                for edge in edges2: edge.plot(ax,'b:',lw=2)
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')
                ax.set_title('All edges for these two points')
                
                for p,edges in intersections:
                    ax = plt.figure().add_subplot(projection='3d')
                    self.plot(ax,'ks-')
                    other.plot(ax,'gs-')
                    ax.plot([p[0]],[p[1]],[p[2]],'ro')
                    for edge in edges: edge.plot(ax,lw=2)
                    ax.set_xlabel('X')
                    ax.set_ylabel('Y')
                    ax.set_zlabel('Z')
                    ax.set_title('Edges for point at %.2f,%.2f,%.2f' % tuple(p.tolist()))
                
                # one of the points is incorrectly labeled "inside"
                p0 = self.edges[0].endpoints[1]
                c0 = other.components[0]
                global DEBUG
                DEBUG = True
                tmp = (p0 in c0)
                print(DEBUG_STATE)
                raise RFLError('Adjacent points have no common edge in overlap')
            # cut the edge to only span the intersection points
            edges.append(common_edge.cut(p1,p2))
        if (len(edges)<=2) and np.sum([e.area for e in edges]):
            # for a valid patch with less than two eges, at least one must 
            # have area (e.g., not be a straight edge)
            return None # overlap has no area
        return Patch(edges)

class Polygon(Patch):
    """Patch composed entirely of StraightEdges"""
    def __init__(self,points,*args,sorted=False,nhat=None,**kwargs):
        """poly = Polygon(points,sorted=False,nhat=None)
        points - list of points that define the polygon
        sorted - if True, points will be assumed to already form a convex polygon
          otherwise, they will be sorted to form a convex polygon
        nhat - ignored if sorted, otherwise used to determine how to sort points
            defines normal vector to polygon. Points will be anticlockwise about nhat
        """
        if not sorted:
            points = sort_points(points,nhat)

        edges = []
        for i,p in enumerate(points):
            p2 = points[(i+1)%len(points)]
            edges.append(StraightEdge(p,p2))
            
        super().__init__(edges)
        self._chord_poly = self # a polygon is its own chord_poly, prevents infinite recursion        
    @classmethod
    def from_edges(cls,edges):
        """*INHERIT*"""
        points = [*edges[0].endpoints]
        for edge in edges[1:-1]:
            points.append(edge.endpoints[1])
        return cls(points)

class Triangle(Polygon,LocalCoordsMixin):
    """Patch composed of a single triangle
    Also implements LocalCoordsMixin to represent local coordinate system
    """
    def __init__(self,*args,**kwargs): # TODO -- add other arguments
        """*INHERIT*
        A Triangle has exactly 3 points
        """
        
        super().__init__(*args,**kwargs)
        if len(self.corners) != 3:
            raise ValueError('Attempt to initialize Triangle with %d points' % len(self.corners))
        self._nhat = None

        ilongest = np.argmax([edge.length for edge in self._edges])
        longest = self._edges[ilongest]
        
        # define local: x along first segment, y = nhat cross x, z = nhat
        # call it in z,x,y order to ensure self.nhat->self._zhat
        (self._zhat,self._xhat,self._yhat) = make_axes((self.nhat,longest.endpoints[1]-longest.endpoints[0],None))
        self._apex3 = self.corners[(ilongest+2)%3] # 3-d apex, edge goes from corners[i] to corners[i+1]
        self._origin = longest.endpoints[0] # a 3-d vector
        self._base = longest # an edge
        self._apex2 = self.p2xy(self._apex3)# (x,y) in local coordinates, y may be negative if triangle entered left-handed
        self._xlim = (0,self._base.length)
    @property
    def components(self):
        if self._components is None:
            self._components = (self,)
        return self._components
    @property
    def area(self):
        """area of Patch"""
        if self._area is None:
            self._area = self._base.length*abs(self._apex2[1])/2 # base*height/2
        return self._area
    @property
    def nhat(self):
        """*INHERIT*"""
        # this is a "component" so it must compute nhat
        if self._nhat is None:
            # unit vector defined by 1st cross 2nd edge
            a = self.corners[1]-self.corners[0] # 1st edge vector
            b = self.corners[2]-self.corners[1] # 2nd edge vector            
            n = np.cross(a,b)
            self._nhat = n/norm(n)
        return self._nhat
    def ylim(self,x):
        """*INHERIT*"""        
        if (x <= 0.0) or (x>=self._base.length):
            ymax = 0.0
        elif x <= self._apex2[0]:
            ymax = self._apex2[1]*x/self._apex2[0]
        elif x < self._base.length:
            ymax = self._apex2[1]*(self._base.length-x)/(self._base.length-self._apex2[0])
        if ymax < 0.0: # shouldn't happen but could
            return (-ymax,0.0)
        else:
            return (0.0,ymax)
    def inside(self,point):
        """*INHERIT*"""
        for i in range(3):
            origin = self.corners[i]
            a = self.corners[(i+1) % 3]-origin
            b = self.corners[(i+2) % 3]-origin
            c = point-origin
            if DEBUG:
                plt.figure()
                xhat,yhat,zhat = make_axes((None,None,self.nhat))
                ax = np.dot(a,xhat)
                ay = np.dot(a,yhat)
                bx = np.dot(b,xhat)
                by = np.dot(b,yhat)
                cx = np.dot(c,xhat)
                cy = np.dot(c,yhat)
                plt.plot([0,ax],[0,ay],'.-')
                plt.plot([0,bx],[0,by],'.-')
                if between_rays(a,b,c,n=self.nhat):
                    plt.plot(cx,cy,'o')
                else:
                    plt.plot(cx,cy,'x')
                plt.title('inside %d' % i)
            if not between_rays(a,b,c,n=self.nhat): return False
        return True
"""    
t = Triangle([[0,0,0],[1,0,0],[0,1,0]])
print('area',t.area) # 0.5 w/in SMALL_FRACTION**2
print('inside',t.inside([0.25,0.25,0])) # True
print('inside',t.inside([0.75,0.75,0])) # False
"""

class Rectangle(Polygon):
    """Patch composed of a rectangle"""
    def __init__(self,points,*args,**kwargs):
        """*INHERIT*
        Rectangle(points,...,sorted=False)
        A Rectangle is specified by three points: (corner,a,b)
        where (a-corner) and (b-corner) must make a right angle
        cross(a,b) should point in the rectangle's normal direciton
        forth corner is found by projecting corner across diagonal ab
        """
        if len(points) != 3:
            raise ValueError('A Rectangle must be initialized with 3 points')
        corner,a,b = np.array(points)
        da = a-corner
        db = b-corner
        if np.dot(da,db)>SMALL_FRACTION*norm(da)*norm(db):
            raise ValueError('A Rectangle must be initialized with the first point being a right angle: corner,a,b')
        points = [corner,a,corner+da+db,b]
        kwargs['sorted'] = True # for sorted arg to be true
        super().__init__(points,*args,**kwargs)
    @classmethod
    def from_edges(cls,edges):
        """*INHERIT*"""
        try:
            points = [*edges[0].endpoints,edges[-1].endpoints[0]]
            return cls(points)
        except ValueError:
            # a projected rectangle may not still be a Rectangle, but is still a polygon
            return Polygon.from_edges(edges)

class Segment(Patch):
    """Patch composed of a single segment of an ellipse segment"""
    def __init__(self,edge,*args,**kwargs):
        """
        seg = Segment(curved,...)
        curved = a CurvedEdge instance
        """
        if not isinstance(edge,CurvedEdge):
            raise ValueError('edge input expected to be a CurvedEdge instance')
        super().__init__((edge,))
    @property
    def components(self):
        if self._components is None:
            self._components = (self,)
        return self._components
    @classmethod
    def from_edges(cls,edges):
        """*INHERIT*"""        
        return cls(edges[0]) # lose subclass-ness

class EllipsePatch(Segment):
    """Patch composed of a full ellipse"""
    def __init__(self,*args,**kwargs):
        """
        seg = EllipsePatch(...)
        produce an EllipsePatch by passing parameters to FullElipse Edge
        Also allows
        seg = EllipsePatch(edge,...) for edge of type FullEllipse
        """
        if isinstance(args[0],FullEllipse):
            super().__init__(*args,**kwargs)
        else:
            super().__init__(FullEllipse(*args,**kwargs))

class CirclePatch(Segment):
    """Patch composed of a full circle"""
    def __init__(self,*args,**kwargs):
        """
        seg = CirclePatch(...)
        produce a CirclePatch by passing parameters to FullCircle Edge
        Also allows
        seg = CirclePatch(edge,...) for edge of type FullCircle
        """
        if isinstance(args[0],FullCircle):
            super().__init__(*args,**kwargs)
        else:
            super().__init__(FullCircle(*args,**kwargs))

inherit_docstrings(Patch)
# end of patches

"""
# Test code for EllipsePatch from_edges
p = EllipsePatch([0,0,0],[1,0,0],[0,1,0])
q = p.from_edges(p.edges)
print(p.edges[0] == q.edges[0]) # True
"""

class AR_Tele_Generic(AngleResponse):
    """
    Generic telescope is meant for combinations of
    circular and rectangular elements, allowed to be at
    offsets and tilted. However, it can handle any convex
    shape whose edges are a combination of lines and elliptical arcs (Patches)
    Initialize with an array of Patch objects that define coincidence geometry
    Call method compute_G to override default integration options for computing
        G==hA0
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        return ('patches' in kwargs)
            
    def __init__(self,patches=None,**kwargs):
        """
        ar = AR_Tele_Generic(patches=None,**kwargs)
        patches - a tuple of Patch objects that must be in coincidence to be considered
        part of the angular response
        Can initialize with a list of patches or using
        conventional keywords and args from the standard
        """
        super().__init__(**kwargs) # handles BIDIRECTIONAL
        if patches is None:
            raise NotImplementedError('keyword standard for initializatoin not defined/implemented yet')
        self.patches = tuple(patches)
        for patch in patches:
            if not isinstance(patch,Patch):
                raise ValueError("Every entry in patches must be a Patch object")
        self.G = None
    def A(self,theta,phi):
        """*INHERIT*"""
        # algorithm:
        # - broadcaast to working shape of at least (1,1)
        # - define nhat from theta,phi (this is particle's momentum axis)
        # - project all patches onto nhat plane
        # - overlap all projected patches
        # - compute area of overlap
        sz = sz0 = np.broadcast(theta,phi).shape
        if len(sz) == 0:
            sz = (1,1)
        elif len(sz)==1:
            sz = (sz[0],1)
        theta = np.atleast_1d(theta)
        if theta.ndim == 1:
            theta = np.expand_dims(theta,1)
        theta = np.broadcast_to(theta,sz)
        phi = np.atleast_1d(phi)
        if phi.ndim == 1:
            phi = np.expand_dims(phi,1)
        phi = np.broadcast_to(phi,sz)
        A = np.zeros(sz)
        it = np.nditer(theta,flags=['multi_index']) # fancy any-shape iterator
        for th in it: # it.multi_index gives corresponding index into same-shape arrays
            if (not self.bidirectional) and (th>90):
                continue
            ph = phi[it.multi_index]
            global DEBUG_STATE
            DEBUG_STATE['A(th,ph)'] =(th,ph)
            nhat = np.array([cosd(ph)*sind(th),sind(ph)*sind(th),cosd(th)])
            projected = []
            # for-else clause executes else clause only if no break
            for p in self.patches:
                # check for grazing incidence
                if abs(np.dot(p.nhat,nhat)/norm(p.nhat)) < SMALL_FRACTION:
                    break # don't do grazing incidence
                # project patch into nhat plane
                projected.append(p.project(nhat))
            else: # only get here if no break
                # progressively shrink overlap by applying overlap of
                # all previous patches with next one

                p = projected[0] # start with full first patch
                for q in projected[1:]:
                    p = p.overlap(q) # overlap what's left with next patch
                    if p is None: break # no overlap
                else: # only get here if no break
                    # compute area of projected overlap
                    A[it.multi_index] = p.area 

        A.shape = sz0 # re-assert original shape
        return A
    def compute_G(self,**kwargs):
        """
        G = compute_G(**kwargs)
        (re)compute G with options passed to scipy dblquad
        kwargs will be passed to dblquad        
        sets self.G, which is used by self.hA0
        if not supplied, the following will be supplied to dblquad
            epsabs = 2e-4*REFERENCE_LENGTH**2
            epsrel = 1e-4
            if only one is supplied, the other will be inferred
            such that epsabs = 2*REFERENCE_LENGTH**2*epsrel
            factor of 2 arises from (1+x*2) = 1+2x+x**2 ~ 1+2x, for x<<1
        """
        
        if ('epsabs' not in kwargs) and ('epsrel' not in kwargs):
            kwargs['epsrel'] = 1e-4 # if both are missing, relative error is known
        
        if 'epsabs' not in kwargs:
            kwargs['epsabs'] = 2*REFERENCE_LENGTH**2*kwargs['epsrel']
        if 'epsrel' not in kwargs:
            kwargs['epsrel'] = kwargs['epsabs']/(2*REFERENCE_LENGTH**2)
        
        # double integral of first and last detector
        # integrate patch-by-patch
        def firstfunc(firsty,firstx,first):
            firstp = first.xy2p(firstx,firsty)
            def lastfunc(lasty,lastx,last):
                lastp = last.xy2p(lastx,lasty)
                # nhat is line connecting (firstx,firsty) and (lastx,lasty)
                r = firstp-lastp
                normr = norm(r)
                rhat = r/normr
                # now check middle patches
                for p in self.patches[1:-2]:
                    # project point onto patch
                    q = project(p,rhat,p.origin)
                    if q not in p: # check inside
                        return 0.0 # not inside, no coincidence
                    # no need to adjust for distance to q
                    # because this filter effectively reduces
                    # the area of the last patch at its real distance
                # compute cos factors for planes relative to line 
                firstcos = abs(np.dot(rhat,first.nhat))
                lastcos = abs(np.dot(rhat,first.nhat))
                return firstcos*lastcos/normr**2
             # sum over components of last detector
            innertmp = 0.0
            for lastc in self.patches[-1].components:
                innertmp += dblquad(lastfunc,lastc.x1,lastc.x2,lastc.y1,lastc.y2,
                                    args=(lastc,),**kwargs)[0]
            return innertmp
        # sum over components of first detector
        tmp = 0.0
        for firstc in self.patches[0].components:
            tmp += dblquad(firstfunc,firstc.x1,firstc.x2,firstc.y1,firstc.y2,
                           args=(firstc,),**kwargs)[0]
        if self.bidirectional:
            tmp *= 2
        self.G = tmp
        return self.G
    @property
    def hA0(self):
        """*INHERIT*
        hA0 (self.G) is computed via a double integral over the areas of the first
        and last patch, applying cosine factors and coincidence conditions on
        any intermediate patches (i.e., a 3-element telescope).
        result is cached in self.G, since it's a slow calculation
        to control the double integral, call compute_G
        """
        if self.G is None:
            self.G = self.compute_G()
        return self.G


inherit_docstrings(AngleResponse) # re-inherit tree

def test_G_thetaphi(ar,G):
    """    
    test calculation of geometric factor using theta,phi double integral and weights
    compare it to ar.hA0 and (if supplied G)
    G = test_G_thetaphi(ar,G=None)
    """
    from tictoc import tic,toc
    
    # typical performance: 
    # ~2-3 minutes for area method regardless of bidirectional
    # 5-10 minutes for angle method, depending on bidirectional
    # area method is faster and more accurate
    # Note that area integral runtime grows with detector area
    # because it has an absolute error requirement

    print('Computing G')

    tic()
    areaG = ar.hA0
    toc()
    print('areaG',areaG)
    
    if G is None:
        G = areaG
    else:
        print('Refernce G',G)
        print('error',areaG/G-1)        

    tic()
    # set up theta,phi grids
    thetagrid = default_thetagrid()
    if not ar.bidirectional:
        thetagrid = thetagrid[thetagrid<=90.0] # acute only
    phigrid = default_phigrid()
    thetagrid,phigrid = broadcast_grids(thetagrid,phigrid)

    # cmopute h for theta,phi grids. Sum over h. apply bidirectional factor
    thetaG = ar.hAthetaphi(thetagrid,phigrid).sum()
    toc()
    print('thetaG',thetaG)
    
    print('error',thetaG/G-1)

    return thetaG

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d # enables 3d projection
    plt.close('all')


    r = Rectangle(([0,0,0],[1,0,0],[0,1,0]))
    print('area',r.area) # 1.0 w/in SMALL_FRACTION**2
    points = [([-1,0.5,0],False),([0,0.5,0],True),([0.5,0.5,0],True),([0.75,0.5,0],True)]
    for p,answer in points:
        guess = r.inside(p)
        if guess != answer:
            print(p,'FAIL','right answer',answer,'guess',guess)
            raise Exception('Inside failed')
    print('PASS: inside tests')
    # test edge intersection
    tests = [(StraightEdge([0,0,0],[1,1,0]),StraightEdge([0,1,0],[1,0,0]),[np.array([0.5,0.5,0])]),
             (StraightEdge([0,0,0],[1,0,0]),StraightEdge([0,0,0],[1,0,0]),[np.array([0,0,0]),np.array([1,0,0])]),
             (StraightEdge([0,0,0],[1,0,0]),StraightEdge([0,0,0],[2,0,0]),[np.array([0,0,0]),np.array([1,0,0])]),
             (StraightEdge([-1,0,0],[2,0,0]),StraightEdge([0,0,0],[1,0,0]),[np.array([0,0,0]),np.array([1,0,0])]),
             (StraightEdge([0,0,0],[1,0,0]),StraightEdge([0.5,0,0],[1.5,0,0]),[np.array([0.5,0,0]),np.array([1,0,0])]),
             (StraightEdge([0,0,0],[1,0,0]),StraightEdge([1,0,0],[2,0,0]),[np.array([1,0,0])]),
             (StraightEdge([0,0,0],[1,0,0]),StraightEdge([2,0,0],[3,0,0]),[]),
             ]
    for edge1,edge2,answer in tests:
        guess = edge1.intersect(edge2)
        if len(guess) != len(answer):
            print('FAIL',edge1.endpoints,'intersect',edge2.endpoints,'right answer',answer,'guess',guess)
            raise Exception('Intersect failed on length check')
        for p1 in answer:
            for p2 in answer:
                if points_equal(p1,p2):
                    break # match found
            else:
                print('FAIL',edge1.endpoints,'intersect',edge2.endpoints,'right answer',answer,'guess',guess)
                raise Exception('Intersect failed on point equality check')
    print('PASS: intersect tests')
            
    r2 = Rectangle(([0,0,2],[1,0,2],[0,1,2]))
    
    # reference rectangle
    tele = AngleResponse(TP_TYPE='RECT_TELE',W1=1,H1=2,W2=3,H2=4,D=5)
    print('Analytical G',tele.hA0,'Matlab',0.804961)

    # off axis for debugging
    first = Rectangle(([0,0,0],[tele.W1,0,0],[0,tele.H1,0]))
    last = Rectangle(([0,0,tele.D],[tele.W2,0,tele.D],[0,tele.H2,tele.D]))
    ar = AR_Tele_Generic(patches=[first,last])
    print('A',ar.A(1.0,359.0))
    print('A',ar.A(1.0,91.0))
    print('A',ar.A(30.0,300.0))
    print('A',ar.A(90.0,0.0))
    print('A',ar.A(np.array([0,1,89,90,91,179,180]),3.0))
    print('A',ar.A(30.0,30.0))
    print('A',ar.A(1.0,90.0))
    print('A',ar.A(1.0,270.0))

    #on axis for comparison to analytical
    first = Rectangle(([-tele.W1/2,-tele.H1/2,0],[tele.W1/2,-tele.H1/2,0],[-tele.W1/2,tele.H1/2,0]))
    last = Rectangle(([-tele.W2/2,-tele.H2/2,tele.D],[tele.W2/2,-tele.H2/2,tele.D],[-tele.W2/2,tele.H2/2,tele.D]))

    ax = plt.figure().add_subplot(projection='3d')
    for rect in [first,last]:
        rect.plot(ax)
        for c in rect.components: # now draw components
            c.plot(ax,linestyle='--')
    for bidir in ['FALSE','TRUE']:
        print('BIDIRECTIONAL',bidir)
        ar = AR_Tele_Generic(patches=[first,last],BIDIRECTIONAL=bidir)
        print('A',ar.A(18.0,52.0))
        test_G_thetaphi(ar,tele.hA0*(1+ar.bidirectional))


"""
arc has property extent which is angular extent in degrees

Break up arcs into pieces floor(extent/180)+1
No need to merge adjacent arcs from same parent curve: their chords will never
be colinear.

Throw error for inside, xlim etc when extent >= 180

Projecting an ellipse and recovering its axes.
https://en.wikipedia.org/wiki/Rytz%27s_construction#Computer_aided_solution
 
First, project the ellipse into the plane by removing all vector components parallel to the normal direction:
x = x-(x.n)n
the original semiaxes are a and b, and the original center is c. The projected axes are a, b, and the projected center is c.
The new ellipse is given by:
r(t) = c+a*cos(t)+b*sin(t)
find t0 such that tan(2t0) = (2*a.b)/(|a|^2-|b|^2)
The semiaxes occur at t0, t0+pi/2, t0+pi, and t0+3pi/2
Test the first two semiaxes at t0, t0+pi/2 . The longer one is the semimajor axis. The semiminor axis is pi/2 later.
 
Computing ellipse arc length:
https://en.wikipedia.org/wiki/Ellipse#Arc_length 
The length of an arc from 0 to alpha for an ellipse with semiaxes of length a,b is:
L(alpha,a,b) = -a*ellipseinc(pi/2-a,sqrt(1-b^2/a^2))
Note Wikipedia reverses a,b.
The arc length from theta1 to theta2 is L(theta2,a,b)-L(theta1,a,b)

"""