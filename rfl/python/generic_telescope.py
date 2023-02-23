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

Will need a way to merge patches with a common edge. We will use 
REFERENCE_LENGTH and SMALL_FRACTION to identify points that are the same.

module settings:
REFERENCE_LENGTH - 1.0. A "typical" length in the units being used
SMALL_FRACTION - 1e-6 - de minimis distance offset. Points closer than this
  fraction of REFERENCE_LENGTH are assumed to be the same point. 
  A value of 1e-6 means that if length units are mm, then any points 
  within a nm are "the same". Since we normally work in cm (or inches) this 
  should be sufficient. Also defines the smallest dot product considered to be 
  zero. Equivalent to two 1 mm long unit vectors whose tips differ by 1 nm

"""

#TODO Creat an edge "cut" method that produces a new edge that only has part of the old edge in it.

#TODO Pull out sort_points.

import numpy as np
from rfl import AngleResponse, inherit_docstrings, sind, cosd

norm = np.linalg.norm # shorthand

REFERENCE_LENGTH = 1.0 # "typical" length in units being used
SMALL_FRACTION = 1e-6 # nm if using mm

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
    # both a x p and p x b must go in same direction as n = a x b
    pn = norm(p)*norm(n)
    if np.dot(np.cross(a,p),n) < -SMALL_FRACTION*norm(a)*pn: return False
    if np.dot(np.cross(p,b),n) < -SMALL_FRACTION*norm(b)*pn:return False
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
        offset -> offset.nhat
    """
    if not np.scalar(offset):
        offset = np.dot(offset,nhat)
    return p-np.dot(p,nhat)*nhat + offset*nhat

# classes
    
class IdentityMixin(object):
    """
    IdentityMixin - mixin to provide identity, from_identity, and __eq__
    functions
    
    """

class LocalCoordsMixin(object):
    """
    LocalCoordsMixin - mixin to provide local coordinate system info
    properties:
        origin - 3-d origin of system
        xhat,yhat,zhat - 3-d unit vectors that define system
        (x1,x2) = xlim - lower and upper limit of x in object's coordinate system
    Properties are read only. Managed by private _ variables that are 
        initilized to None in the Edge base class constructor. 
        Only getters are provided.
    methods:
        (y1,y2) = ylim(x) - y limits in object's coordinate sytem
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
    def ylim(self,x):
        """(y1,y2) = ylim(x) - y limits in object's coordinate sytem"""
        raise NotImplementedError # Abstract method
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
        a,b = endpoints returns Edges' end points
            or None for closed edges (e.g., full ellipses)
        len = length - length of edge
        area - area of edge for curved edes, zero for straight
        nhat - unit vector normal to Edge's plane
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
    @property
    def identity(self):
        """points = identity returns tuple of points that uniquely identify Edge"""
        return self._identity
    @classmethod
    def from_identity(cls,identity):
        raise NotImplementedError # abstract method
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
    def intersect(self,other,try_other=True):
        """points = intersect(other,try_other=True) 
        returns a typle of points of intersection between self and Edge other
        try_other: try other's intersect method if this one doesn't support other's type
        """
        # we got here because self did not implement intersect with otheer        
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
        return self.from_id(project(np.array(self.identity),nhat,offset))
        
    def inside(self,point):
        """
        bool = inside(point) returns True if a point is inside the edge
        can also accept list of points and return list of bools
        """
        return NotImplementedError # Abstract base class
    def __contains__(self,point): # in operator calls inside
        return self.inside(point)
    def plot(self,ax,*args,**kwargs):
        """
        h = plot(self,ax,*args,**kwargs)
        call ax.plot on axis with args and kwargs
        returns result of ax.plot
        """
        return NotImplementedError # Abstract base class
        
    
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
    @classmethod
    def from_identity(cls,id):
        """*INHERIT*"""
        return cls(*id) # identity is endpoints
    @property
    def length(self):
        """length of Edge"""
        if self._length is None:
            self._length = norm(self._endpoints[1]-self._endpoints[0])
        return self._length
    def _intersect_with_straight(self,other):
        """
        same as intersect, but when other is known to be a StraightEdge
        """
        # self is: p(t) = p1+(p2-p1)*t
        # other is: q(s) = q1+(q2-q1)*s
        # find t,s such that p(t) = q(s)
        # solve A*x = y
        # A = [p2-p1,q2-q1] (stack as two rows)
        # y = q1-p1
        # x = (t,s)
        # use least squares since A is overspecified
        (p1,p2) = self.endpoints[0]
        (q1,q2) = other.endpoints[1]
        A = np.column_stack([p2-p1,q1-q2])
        y = q1-p1
        (t,s) = np.linalg.solve(A,y,rcond=None)[0] # 1st output is solution
        empty = tuple()
        if (t<-SMALL_FRACTION) or (t>1.0+SMALL_FRACTION):
            return empty # intersection is not between p1,p2
        if (s<-SMALL_FRACTION) or (s>1.0+SMALL_FRACTION):
            return empty # intersection is not between q1,q2
        if not points_equal(p1+(p2-p1)*t,q1+(q1-q2)*s):
            return empty # points do not intersect
        return (p1+(p2-p1)*t,) # return tuple of single intersecting point        
    def intersect(self,other,try_other=True):
        """*INHERIT*"""
        # try intersect with types self knows about
        # then hand off to parent class
        if isinstance(other,StraightEdge):
            return self._intersect_with_straight(other)
        else:
            return super().intersect(other,try_other) # pass up the chain
    def inside(self,point):
        """*INHERIT*
        always false for straight edges
        """
        return False
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
        raise NotImplementedError # Abstract method

class EllipseArc(CurvedEdge):
    """EllipseArc Edge, subclass of CurvedEdge, for elliptical arcs"""
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
        super().__init__()
        r1 = np.array(r1)
        r2 = np.array(r2)
        if np.isscalar(p1): # theta in degrees
            p1 = center+r1*cosd(p1) + r2*sind(p1)
        if np.isscalar(p2): # theta in degrees
            p2 = center+r1*cosd(p2) + r2*sind(p2)
        self._endpoints = (p1,p2)
        self._identity = self._endpoints
        # raise NotImplementedError # TODO
    @classmethod
    def from_identity(cls,id):
        """*INHERIT*"""
        raise NotImplementedError # TODO
    @property
    def length(self):
        """
        *INHERIT*
        length is computed along arc
        """
        # https://math.stackexchange.com/questions/433094/how-to-determine-the-arc-length-of-ellipse
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipeinc.html#scipy.special.ellipeinc
        raise NotImplementedError # TODO
    @property
    def area(self):
        """*INHERIT*"""
        raise NotImplementedError # TODO
    @property
    def xlim(self):
        """*INHERIT"""
        return NotImplementedError # TODO
    def ylim(self,x):
        """*INHERIT"""
        return NotImplementedError # TODO    
    def xy2p(self,x,y):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def intersect(self,other,try_other=True):
        """*INHERIT*"""
        # try intersect with types self knows about
        # then hand off to parent class
        if isinstance(other,EllipseArc):
            raise NotImplementedError # TODO
        elif isinstance(other,StraightEdge):
            raise NotImplementedError # TODO
        else:
            super().intersect(other,try_other) # pass up the chain
    def rotate(self,rot,project=True):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def inside(self,point):
        """*INHERIT*"""
        raise NotImplementedError # TODO

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
    @classmethod
    def from_identity(cls,id):
        """*INHERIT*"""
        raise NotImplementedError # TODO
    @property
    def area(self):
        """*INHERIT*
        area inside closed edge
        """
        raise NotImplementedError # TODO
    @property
    def xlim(self):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def ylim(self,x):
        """*INHERIT"""
        return NotImplementedError # TODO    
    def intersect(self,other,try_other=True):
        """*INHERIT*"""
        # try intersect with types self knows about
        # then hand off to parent class
        if isinstance(other,EllipseArc):
            raise NotImplementedError # TODO
        elif isinstance(other,StraightEdge):
            raise NotImplementedError # TODO
        else:
            super().intersect(other,try_other) # pass up the chain
    def rotate(self,rot,project=True):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def inside(self,point):
        """*INHERIT*"""
        raise NotImplementedError # TODO

class CircleArc(EllipseArc):
    """CircleArc Edge, subclass of EllipseArc, for circular arcs"""
    def __init__(self,center,p1,p2):
        """arc = CircleArc(center,p1,p2)
        center- 3-d location of circle's center
        p1 - 3-d location of first point on arc
        p2 - 3-d location of second point on arc
        arc spans from p1 to p2 anticlockwise (around p1 x p2)
        """
        raise NotImplementedError # TODO
    @classmethod
    def from_identity(cls,id):
        """*INHERIT*"""
        raise NotImplementedError # TODO

class FullCircle(FullEllipse):
    """FullCircle Edge, subclass of Ellipse, full circle"""
    def __init__(self,center,r,nhat):
        """circle = FullCircle(center,r)
        center - 3-d location of circle's center
        r - radius of circle (scalar)
        nhat - 3-vector normal to plane of circle
        """
        raise NotImplementedError # TODO
    @classmethod
    def from_identity(cls,id):
        """*INHERIT*"""
        raise NotImplementedError # TODO

inherit_docstrings(Edge)
# end of Edges

class Patch(object):
    """
    Convex Patch object consisting of straight and curved edges
    convex patches only
    properties:
        edges = tuple of Edges that make up the Patch
        corners - tuple of points that make up the Patch (edge endpoints)
          edges and corners are ordered anticlockwise
        area - area of patch
        nhat - unit vector normal to Pacth's plane
        components - constituent component Patches
        chord_poly - a polygon made up of the endpoints for all edges
        Properties are read only. Managed by private _ variables and getters
    methods:
        bool = inside(point) returns True if a point is inside the Patch
            also available as "in" operator
        patch = project(nhat,offset) returns the patch projected along nhat, with offset
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

        edges = tuple(edges)
        self._edges = tuple(edges)
        
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
        """a Polygon made from the Patch's corners"""
        if self._chord_poly is None:
            self._chord_poly = Polygon(self._corners)
        return self._chord_poly
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
    def overlap(self,other,try_other=True):
        """patch = overlap(other) returns patch representing overlap of self and other"""
        # we got here because self did not implement intersect with otheer        
        return NotImplementedError # TODO

class Polygon(Patch):
    """Patch composed entirely of StraightEdges"""
    def __init__(self,points,*args,sorted=False,nhat=None,**kwargs):
        """poly = Polygon(points,sorted=False)
        points - list of points that define the polygon
        sorted - if True, points will be assumed to already form a convex polygon
          otherwise, they will be sorted to form a convex polygon
        nhat - ignored if sorted, otherwise used to determine how to sort points
            defines normal vector to polygon. Points will be anticlockwise about nhat
        """
        if not sorted:
            # sort points such that angles increase in nhat plane
            center = np.zeros(3) # average point
            points = list(points) # make mutable
            for i,p in enumerate(points):
                if np.size(p) != 3:
                    raise ValueError('Polygon requires points input to be list of 3-d vectors')
                points[i] = p = np.array(p).ravel()
                center += p
            center = center/len(points) # average
            zhat = nhat
            if zhat is None:
                zhat = 0
                for i,p in enumerate(points):
                    zhat = np.cross(p-center,points[(i+1) % len(points)]-center)
                    if norm(zhat)>SMALL_FRACTION*REFERENCE_LENGTH: break
            if norm(zhat) <= REFERENCE_LENGTH*SMALL_FRACTION:
                raise ValueError('points provided to Polygon are degenerate (co-linear)')
            (zhat,xhat,yhat) = make_axes((zhat,points[0]-center,None))
            angles = []
            for (i,p) in enumerate(points):
                dp = p-center
                x = np.dot(dp,xhat)
                y = np.dot(dp,yhat)
                z = np.dot(dp,zhat)
                if abs(z) > SMALL_FRACTION*norm(dp):
                    raise ValueError('Polygon input points not coplanar')
                points[i] = p = center + x*xhat+y*yhat # ensure coplanar
                angles.append(np.arctan2(y,x))
            angles,isort = np.unique(angles,return_index=True) # remove duplicates
            if len(angles) < 3:
                raise ValueError('All points in polygon are coplanar')
            points = [points[i] for i in isort]

        edges = []
        for i,p in enumerate(points):
            p1 = points[i]
            p2 = points[(i+1)%len(points)]
            edges.append(StraightEdge(p1,p2))
            
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
            ymax = self._apex2[1]*x/self.apex2[0]
        elif x < self._base.length:
            ymax = self._apex2[1]*(self._base.length-x)/self._base.length
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
        points = [*edges[0].endpoints,edges[-1].endpoints[0]]
        return cls(points)

r = Rectangle(([0,0,0],[1,0,0],[0,1,0]))
print('area',r.area) # 0.5 w/in SMALL_FRACTION**2
print('inside',r.inside([0.2,0.1,0])) # True
print('inside',r.inside([1.75,0.75,0])) # False
print('inside',r.inside([0.1,0.1,0])) # True
r2 = Rectangle(([0,0,2],[1,0,2],[0,1,2]))
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d # enables 3d projection
plt.close('all')
ax = plt.axes(projection='3d')
for r_ in [r,r2]:
    r_.plot(ax)
    for c in r_.components: # now draw components
        c.plot(ax,linestyle='--')

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

#"""
# Test code for EllipsePatch from_edges
p = EllipsePatch([0,0,0],[1,0,0],[0,1,0])
q = p.from_edges(p.edges)
print(p.edges[0] == q.edges[0]) # True
#"""

class AR_Tele_Generic(AngleResponse):
    """
    Generic telescope is meant for combinations of
    circular and rectangular elements, allowed to be at
    offsets and tilted. However, it can handle any convex
    shape whose edges are a combination of lines and elliptical arcs
    Initialize with an array of Patch objects that define coincidence geometry
    """

inherit_docstrings(AngleResponse) # re-inherit tree

"""
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