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
the area between the curve and it's chord. Then getting an area of a patch 
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
...Ellipse
....Circle

.intersect returns None if class cannot do intersection with other class. 
Returns list of intersecting points otherwise. Algorithm:
self knows how to do intersection? If no, check super. At Edge level, 
check other.intersect(self, False). False prevents endless recursion. 
If still no results raise error.

Patch
. Polygon
.. Triangle
. Rectangle
.. Square
. Segment

.xy2v converts from local coords x,y to 3-d vector, e.g p+x*a+y*b for ellipse.
Triangles generate their own xy2v
Segments, Circles, and Ellipses get their xy2v from the edge

.overlap reruns a new patch that is the overlap of self with another patch.

Will need a way to merge patches with a common edge. Will work best if can 
do pointer compare (is) vs numeric compare. Can hopefully compare edge type 
value and pointer compare chord ends. Will use a Rotation object to rotate
points, which will cache any points it's already seen by pointer and return 
same pointer for any other rotations of that same point by other edges. 
This should enable pointer comparisons of points without worrying about 
floating point equality issues. Rotation can cache with dict of points' 
tobytes(). That's an efficient way to make a unique hashable value from a 
numpy array.

"""

import numpy as np
from rfl import AngleResponse, inherit_docstrings

class Rotation(object):
    """
    rotate a point from one 3-d system to another
    caches points to so that if the same point is rotated
    again, its prior result is returned (same object)
    propertis:
        R - (3,3) rotation matrix
        Rt - (3,3) transpose of R (reverse rotation matrix)
        cache - cached points that were prevously rotated
        reverse_cache - chace points that were previously inverse rotated
        (cache is used to ensure points are the same object if rotated again later)
    methods:
        rotate - rotate from old to new coordinates
        reverse_rotate - rotate from new to old coordinates
        inverse - return new rotation
    using cache is actually probably slower on large rotations, but it
        preserves object identity for "is" vs nearly-equal comparisons
    """
    def __init__(self,xhat,yhat,zhat):
        """
        Rotation(xhat,yhat,zhat)
        unit vectors of new system in old system
        """
        for hat in (xhat,yhat,zhat):
            if np.size(hat) != 3:
                raise ValueError('xhat, yhat, and zhat must each have 3 elements')
        self.cache = {}
        self.reverse_cache = {}
        self.R = np.row_stack((np.array(xhat).ravel(),np.array(yhat).ravel(),np.array(zhat).ravel()))
        self.Rt = self.R.transpose()
    def _rotate(self,x,R,use_cache,cache,reverse_cache):
        """
        y = _rotate(x)
        rotate x using rotation R
        with bool flag use_cache
        cache stores the x->y map
        reverse_cache stores the y->x map        
        """
        shape = np.shape(x)
        x = np.atlesat_2d(x) # (N,3)
        if use_cache:
            y = []
            for ix in x:
                k = x.tobytes()
                if k in self.cache:
                    iy = self.cache[k]
                else:
                    iy = self.dot(self.R,ix)
                    cache[k] = iy
                    reverse_cache[iy.tobytes()] = ix
                y.append(iy)
            y = np.array(y)
        else:
            y = np.tensordot(x,R,axes=(1,1)) # 2nd dim of R time 2nd dim of x
        y.shape = shape
        return y
    def rotate(self,x,use_cache=False):
        """
        y = rotate(x,use_cache=True)
        rotate x into the new coordinate system
        x can be (3,) or (N,3)
        y will be the same shape
        use_cache = use cache or not
        """
        return self._rotate(x,self.R,use_cache,self.cache,self.reverse_cache)
    def reverse_rotate(self,y,use_cache=True):
        """
        x = reverse_rotate(y,use_cache=True)
        rotate y into the old coordinate system
        y can be (3,) or (N,3)
        x will be the same shape
        use_cache = use cache or not
        """
        return self._rotate(y,self.Rt,use_cache,self.reverse_cache,self.cache)
    def inverse(self):
        """return the inverse rotation as a new Rotation object"""
        return Rotation(self.Rt[0,:],self.Rt[1,:],self.Rt[2,:])
        

class Edge(object):
    """
    Abstract Base Class for edges
    properties:
        a,b = chord returns StraightEdge that closes edge
            or None for closed edges (e.g., full ellipses)
        area - area of edge for curved edes, zero for straight
        nhat - unit vector normal to Edge's plane
        (x1,x2) = xlim - lower and upper limit of x in Patch's coordinate system
    methods:
        (y1,y2) = ylim(x) - y limits in Edge's coordinate sytem
        v = xy2v(x,y) convert x,y from Edge's coordinates to point in 3-d frame
        points = intersect(other) returns a list of points of intersection between self and other
        edge = rotate(rot) returns the edge rotated by Rotation rot
        bool = inside(point) returns True if a point is inside the edge
            Always False for straight edges
            True for curved eges if point lies between curve and chord
    """
    @property
    def chord(self):
        """a,b = chord() returns end points of Edge's chord"""
        return NotImplementedError # Abstract base class
    @property
    def area(self):
        """area of Edge"""
        return NotImplementedError # Abstract base class
    @property
    def nhat(self):
        """nhat - unit vector normal to Edge's plane"""
        return None
    @property
    def xlim(self):
        """(x1,x2) = xlim - lower and upper limit of x in Edge's coordinate system"""
        return NotImplementedError # Abstract base class
    def ylim(self,x):
        """(y1,y2) = ylim(x) - y limits in Edge's coordinate sytem at x"""
        return NotImplementedError # Abstract base class
    def xy2v(self,x,y):
        """v = xy2v(x,y) convert x,y from Edge's coordinates to point in 3-d frame"""
        return NotImplementedError # Abstract base class
    def intersect(self,other,try_other=True):
        """points = intersect(other,try_other=True) 
        returns a list of points of intersection between self and Edge other
        try_other: try other's intersect method if this one doesn't support other's type
        """
        # we got here because self did not implement intersect with otheer        
        if try_other:
            return other.intersect(self,try_other=False) # let other try
        else: # oops, coding error, neither class has intersect for the other
            raise NotImplementedError('No intersection method found for %s,%s' %
                                      (self.__class__.__name__,other.__class__.__name__))
    def rotate(self,rot):
        """
        edge = rotate(rot) returns the edge rotated by Rotation rot
        """
        return NotImplementedError # Abstract base class
    def inside(self,point):
        """
        bool = inside(point) returns True if a point is inside the edge
        can also accept list of points and return list of bools
        """
        return NotImplementedError # Abstract base class
    
class StraightEdge(Edge):
    """StraightEdge implements Edge for edges that are line segments"""
    def __init__(self):
        raise NotImplementedError # TODO
    @property
    def chord(self):
        """*INHERIT*
        returns self for straight edges
        """
        return self
    @property
    def area(self):
        """*INHERIT*
        returns zero for straight edges
        """
        return 0.0
    @property
    def xlim(self):
        """*INHERIT"""
        return NotImplementedError # TODO
    def ylim(self,x):
        """*INHERIT"""
        return NotImplementedError # TODO    
    def xy2v(self,x,y):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def intersect(self,other,try_other=True):
        """*INHERIT*"""
        # try intersect with types self knows about
        # then hand off to parent class
        if isinstance(other,StraightEdge):
            raise NotImplementedError # TODO
        else:
            super().intersect(other,try_other) # pass up the chain
    def rotate(self,rot):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def inside(self,point):
        """*INHERIT*
        always false for straight edges
        """
        return False

class CurvedEdge(Edge):
    """CurvedEdge base subclass of Edge for curved edges"""
    @property
    def chord(self):
        """*INHERIT*
        Returns a straight edge that closes the curve
        """
        raise NotImplementedError # Abstract method
    @property
    def area(self):
        """*INHERIT*
        area between curve and chord
        """
        raise NotImplementedError # Abstract method

class EllipseArc(CurvedEdge):
    """EllipseArc Edge, subclass of CurvedEdge, for elliptical arcs"""
    def __init__(self):
        raise NotImplementedError # TODO
    @property
    def chord(self):
        """*INHERIT*"""
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
    def xy2v(self,x,y):
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
    def rotate(self,rot):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def inside(self,point):
        """*INHERIT*"""
        raise NotImplementedError # TODO

class FullEllipse(EllipseArc):
    """FullEllipse Edge, subclass of EllipseArc, full ellipse"""
    def __init__(self):
        raise NotImplementedError # TODO
    @property
    def chord(self):
        """*INHERIT*
        closed curve, has no chord (always None)
        """
        raise None
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
    def rotate(self,rot):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def inside(self,point):
        """*INHERIT*"""
        raise NotImplementedError # TODO

class FullCricle(FullEllipse):
    """FullCircle Edge, subclass of Ellipse, full circle"""
    def __init__(self):
        raise NotImplementedError # TODO

inherit_docstrings(Edge)
# end of Edges

class Patch(object):
    """
    Convex Patch object consisting of straight and curved edges
    convex patches only
    properties:
        is_component - bool if this Patch is a component (Triangle or CurvedEdge)
        edges = list of Edges that make up the Patch
        corners - list of points that make up the Patch (chord ends)
          edges and corners are ordered anticlockwise
        area - area of patch
        nhat - unit vector normal to Pacth's plane
        components - constituent component Patches
        (x1,x2) = xlim - lower and upper limit of x in Patch's coordinate system
    methods:
        (y1,y2) = ylim(x) - y limits in Edge's coordinate sytem
        v = xy2v(x,y) convert x,y from Edge's coordinates to point in 3-d frame
        bool = inside(point) returns True if a point is inside the Patch
        patch = rotate(rot) returns the Patch rotated by Rotation rot
        patch = overlap(other) returns patch representing overlap of self and other
        
    Note: only component Patches (Triangles and Curved Edges) implement xlim, ylim, and xy2v
    """
    def __init__(self):
        self._components = None
        self._area = None
        self._edges = [] # list of edges that make up this Patch
        self._corners = [] # list of corners (chord ends) tha tmake up this patch
    def __str__(self):
        return 'Patch %s' % (','.join([e.__class__.__name__ for e in self.edges]))
    @property
    def is_component(self):
        return False
    @property
    def edges(self):
        """a list of Edges that make up the Patch"""
        return self._edges
    @property
    def corners(self):
        """a list of points that make up the Patch (from each Edge.chord)"""
        return self._corners
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
        return NotImplementedError # Abstract base class
    @property
    def xlim(self):
        """(x1,x2) = xlim - lower and upper limit of x in Patch's coordinate system"""
        return NotImplementedError # Abstract base class
    def ylim(self,x):
        """(y1,y2) = ylim(x) - y limits in Patch's coordinate sytem at x"""
        return NotImplementedError # Abstract base class
    def xy2v(self,x,y):
        """v = xy2v(x,y) convert x,y from Patch's coordinates to point in 3-d frame"""
        return NotImplementedError # Abstract base class
    def inside(self,point):
        """
        bool = inside(point) returns True if a point is inside the patch
        can also accept list of points and return list of bools
        """
        for c in self.components:
            if c.inside(point):
                return True
        return False
    def rotate(self,rot):
        """
        patch = rotate(rot) returns the patch rotated by Rotation rot
        """
        return NotImplementedError # Abstract base class
    def overlap(self,other,try_other=True):
        """patch = overlap(other) returns patch representing overlap of self and other"""
        # we got here because self did not implement intersect with otheer        
        return NotImplementedError # TODO

class Polygon(Patch):
    """Patch composed entirely of StraightEdges"""
    pass # TODO

class Triangle(Polygon):
    """Patch composed of a single triangle"""
    @property
    def is_component(self):
        return True
    @property
    def nhat(self):
        """*INHERIT*"""
        # unit vector defined by 1st cross 2nd edge
        a = self.corners[1]-self.corners[0] # 1st edge vector
        b = self.corners[2]-self.corners[1] # 2nd edge vector
        n = np.cross(a,b)
        return n/np.sum(n**2)
    @property
    def components(self):
        """*INHERIT*"""
        return [self]
    @property
    def xlim(self):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def ylim(self,x):
        """*INHERIT*"""
        return NotImplementedError # TODO
    def xy2v(self,x,y):
        """*INHERIT*"""
        return NotImplementedError # TODO
    pass # TODO
    
class Rectangle(Polygon):
    """Patch composed of a rectangle"""
    pass # TODO

class Square(Rectangle):
    """Patch composed of a square"""
    pass # TODO

class Segment(Patch):
    """Patch composed of a segment of an ellipse segment"""
    @property
    def is_component(self):
        return True
    pass # TODO

class EllipsePatch(Segment):
    """Patch composed of a segment of an ellipse"""
    pass # TODO

class CirclePach(EllipsePatch):
    """Patch composed of a segment of a circle"""
    pass # TODO

inherit_docstrings(Patch)
# end of patches

class AR_Tele_Generic(AngleResponse):
    """
    Generic telescope is meant for combinations of
    circular and rectangular elements, allowed to be at
    offsets and tilted. However, it can handle any convex
    shape whose edges are a combination of lines and elliptical arcs
    """

inherit_docstrings(AngleResponse) # re-inherit tree
