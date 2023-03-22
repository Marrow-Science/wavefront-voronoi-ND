import math
import numpy
import scipy
#from scipy.optimize import root_scalar as solve
from scipy.optimize import brentq as solve
from numpy.polynomial import Polynomial as poly
# For pairs of points
from itertools import combinations
# To exit cleanly
import sys

# We create a weighted Voronoi cell in ND using a wavefront algorithm

# Globals are bad, but hey...
EPS = 0.001
BOUND = 10.0

def kwargDef(arg, args, default):
	if arg in args: return args[arg]
	else: return default

# The form determines what the (hyper)plane of intersection is, and its
# orthogonal directions and weights
def formDir(dirs, N):
	ret = [0.0] * N
	for x in dirs: ret = vsum(ret, x)
	if magnitude(ret) < EPS: return ret
	return scale(ret, 1.0 / magnitude(ret))

# The main algorithm, add new waves on each collision
def runWevoNd(wavefront, eps = EPS):
	event, eventHandler = wavefront.nextEvent(eps = eps)
	# Handle events untill none are left
	while not event == None:
		eventHandler(event)
		event, eventHandler = wavefront.nextEvent(eps = eps)
	verticies = wavefront.getVerticies()
	base = {}
	for b in wavefront.base:
		wid = wavefront.ids[b]
		base[wid] = b
	return partition(base, verticies)

# Debug printing statements
def waveTreePrint(WF, wave,idx):
	wid = waveID(WF,wave)
	if wave in WF.alias: wid = WF.alias[wave] + "(" + str(wid) + ")"
	print("\t"*idx, wid, wave.span)
	if wave.leaf(): return
	p1, p2 = wave.parents()
	waveTreePrint(WF, p1, idx + 1)
	waveTreePrint(WF, p2, idx + 1)

def debugPrint(WF, nw):
	print("-------- DEBUG PRINT ---------")
	p1, p2 = nw.parents()
	print("PARENT N:",p1.N(), p2.N())
	waveTreePrint(WF, nw,1)
	print("------ END DEBUG PRINT -------")

# A subspace 
class subspace:
	def __init__(self, space, eps = EPS):
		self.vec, self.basis = self.cleanInit(space, eps)
		self.dim = len(self.vec)
	# Clean and initialize the space, calculate orthonormal basis
	def cleanInit(self, space, eps):
		# Clean zero vectors
		nozero = [s for s in space if magnitude(s) > eps]
		if len(nozero) == 0: return [], []
		vec, basis = [nozero[0]], [norm(nozero[0])]
		for x in range(1,len(nozero)):
			rej = reject(basis,[nozero[x]])[0]
			if magnitude(rej) > eps:
				vec.append(nozero[x])
				basis.append(norm(rej))
		return vec, basis
	def dim(self): return self.dim
	def X(self):
		ret = None
		for x in self.vec:
			ret = vsum(x,ret)
		return ret

class form:
	# TODO: Forms needs to take a centerline and build around an X vector
	# TODO: Forms need linear dependence code
	def __init__(self, center, I, p = None):
		if center == None:
			print("NO CENTER!", base, center)
			return None
		self.base = [center] if p == None else p[0].base + p[1].base
		wvec = [vec(b, center) for b in self.base]
		v = None
		for x in wvec: v = vsum(x, v)
		self.X = norm(v)
		space = [I] if p == None else [I] + p[0].dir + p[1].dir
		self.space = subspace(space)
		self.dir = self.space.vec
		self.dim = self.space.dim

	def N(self): return self.dim
	# To make the orthogonal dual do a series of rejections
	# Rejections maintain orthogonality over orthogonal bases
	def comp(self, compN, eps = EPS):
		N = self.N()
		comp = []
		i = 0
		while len(comp) + N < compN:
			ih = [hat(compN,i)]
			nxt = reject(self.space.basis + comp, ih)[0]
			i += 1
			if magnitude(nxt) < eps: continue
			comp.append(norm(nxt))
		return comp		

def interPlane(f1, f2):
	#print("INTERPLANE",f1.N(), f2.N())
	x1, x2 = (f1.X, f2.X)
	r1, r2 = (reject([x1],[x2]), reject([x2],[x1]))
	return subspace(r1 + r2)

# Strategy: find subspace of intersection, project X's onto it, find R's
def formInter(f1, f2):
	# Zero forms have no intersections
	if f1.N == 0 or f2.N == 0: return None
	plane = interPlane(f1, f2)
	p1, p2 = (project(plane, [f1.X]), project(plane, [f2.X]))
	return plane, p1, p2

# A wave starts from point (baseWave) or intersection (wave)
# A wave is always orthogonal to a form (0-form in the case of a point)
# A wave can be either convex or flat depending on weight ratios
class baseWave:
	def __init__(self, center, weight, N):
		# TODO: BETTER DEBUG!
		self.debug = []
		# TODO: BETTER DEBUG!
		self.weight = weight
		self.center = (center, None)
		self.span = (0.0, None)
		# TODO: split into function
		self.base = [center]
		self.form = form(center,[])
		self.dim = N - self.form.N() - 1
		# TODO: this?
		self.x1p = lambda T: center
		self.x2p = lambda T: center
	# If this is a leaf node
	def leaf(self): return True
	def valid(self): return True
	def N(self): return self.dim
	def mag(self): return self.weight
	def theta(self, T): return None
	def C(self, T): return [0.0]
	def L(self, T, clamp = True): return self.center[0]
	def R(self, T): return self.weight * T

class wave:
	# Basic initialization from a pair of intersecting parents
	def __init__(self, p1, p2, N):
		# TODO: BETTER DEBUG
		self.debug = []
		# TODO: BETTER DEBUG
		self.p1 = p1
		self.p2 = p2
		self.P, self.T, self.I, self.D, self.C2L = self.interInit(N)
		self.span, self.center, self.val = self.spanInit()
		# If the initialization is not valid just leave
		if not self.val: return
		# Otherwise the span is valid
		'''
		#TODO: better centerline handling for 4D+.
		Needs I to be between span[0] and span[1], make X direct from
		'''
		start = self.span[0]
		center = self.center[0]
		I = norm(self.I(self.span[0]))
		''''''
		self.form = form(center, I, (p1.form, p2.form))
		self.dim = N - (self.form.N() + 1)
		self.val = True
	# Init various intersection values needed to calculate the waveform
	def interInit(self, N, eps = EPS):
		# Local variable initialization
		p1, p2 = self.parents()
		c1, c2 = (None, None)
		c1, c2 = (p1.center[0], p2.center[0])
		f1, f2 = (p1.form, p2.form)
		# Lambdas for caching parent wave locations
		self.x1p = lambda T: p1.L(T)
		self.x2p = lambda T: p2.L(T)
		# 3 cases for wave intersection, all derive from same formula
		left, right = (f1.N() == 0, f2.N() == 0)
		if left and right:
			return self.interInit0x0(c1, c2)
		elif left or right:
			return self.interInit0xN(p1, p2, c1, c2, f1, f2)
		else:
			raise Exception("NxN intersections not needed")

	# Base degenerate case of two undirected wavefronts
	# The intersection plane is 1D
	# The POI intersects both parent centers
	def interInit0x0(self, c1, c2):
		Iv = lambda T: vec(c1, c2)
		Is = lambda T: magnitude(Iv(T))
		P = lambda T: [norm(vec(c1, c2))]
		theta = None
		C2L = None
		D = (None, None)
		return P, theta, Iv, D, C2L

	# Partial case of one directed and one undirected wavefront
	# In this case, the intersection plane is 1D
	# But the POI does not intersect both parent centers
	def interInit0xN(self, p1, p2, c1, c2, f1, f2):
		# Xp, the centers of each wave (lambdas)
		x1p, x2p = (self.x1p, self.x2p)
		# I, the vector from x1 to x2
		i1v = lambda T: vec(x1p(T),x2p(T))
		i2v = lambda T: vec(x2p(T),x1p(T))
		# D, vector from each X to the plane of intersection
		L, R = (f1.N() == 0, f2.N() == 0)
		d1 = None if R else lambda T: project([f2.X], [i1v(T)])[0]
		d2 = None if L else lambda T: project([f1.X], [i2v(T)])[0]
		D = (d1, d2)
		# P, the intersection plane
		F = f1.X if R else f2.X
		B = lambda T: reject([F], [i1v(T)])[0]
		P = lambda T: [norm(B(T))]
		Ip = lambda T: B(T)
		theta = None
		C2L = None
		return P, theta, Ip, D, C2L

	def spanInit(self):
		time,val = timeSpan(self)
		cen,_ = centerSpan(self, time)
		return time, cen, val

	# If this is a leaf node
	def leaf(self): return False
	# If the parents actually intersect at some point
	def valid(self): return self.val
	# Return the dimensionality of the wave
	def N(self): return self.dim
	# To make things simpler
	def parents(self): return self.p1, self.p2
	# The theta between centers of the two parent forms
	def theta(self, T): return None if self.T == None else self.T(T)
	# The distance along the parent form at time T!
	def C(self, T, clamp = True):
		# Clamp the output if C is called outside span boundaries
		if clamp and not inSpan(T, self.span):
			ret,valid = clampSpan(self, T, EPS)
			if not valid: print("Impossible situation occured")
			return ret
		# Calculate the form intersection
		theta = self.theta
		p1, p2 = self.parents()
		# Project waves onto the intersection subspace
		r1p = projectPOI(p1.R, self.D[0])
		r2p = projectPOI(p2.R, self.D[1])
		# Linear equation for two spheres
		I = lambda t: magnitude(self.I(t))
		f0 = lambda t: 2.0 * I(t)
		f1 = lambda t: r2p(t)**2.0 - I(t)**2.0 - r1p(t)**2.0
		return linear(f0,f1)(T)

	# The current location of the wave center
	def L(self, T, clamp = True):
		d = None if self.D[0] == None else self.D[0](T)
		C = scale(self.P(T)[0], self.C(T, clamp))
		ret = vsum(self.p1.L(T, clamp), vsum(d, C))
		return ret

	# The radius perpendicular to the form at time T
	def R(self, T):
		# Find projected C and R
		Cs = self.C(T)
		Rp = projectPOI(self.parents()[0].R,self.D[0])(T)
		# Find radius
		radius = abs((Rp ** 2.0) - (Cs ** 2.0)) ** 0.5
		return radius

# A vertex of a cell, computed by either wave cessation or intersection
# The set of all vertexes can be used to quickly reconstruct full cells
def waveVertexThunk(WF, wav, eps = EPS):
	# All vertexes have points ids and parents
	vid = WF.ids[wav]
	pt = wav.center[0]
	parents = waveID(WF, wav)
	# Some vertexes are just points in space
	if wav.N() == -1: return vertex(vid, pt, None, parents, None)
	# All others have parents and forms
	space = wav.form.space
	# Some vertexes are flat (either by stuff or things)
	if wav.center[1] == None: return vertex(vid, pt, None, parents, space)
	# If spans are too small flatten the vertex? TODO: is this the right location?
	span = magnitude(vec(wav.center[0],wav.center[1]))
	if span < eps: return vertex(vid, pt, None, parents, space)
	# Others are curved
	center = scale(vsum(wav.center[0],wav.center[1]), 0.5)
	return vertex(vid, pt, center, parents, space)

class vertex:
	def __init__(self, vid, vertex, center, parents, space):
		self.vid = vid
		self.vertex = vertex
		self.center = center
		self.parent = parents
		self.space = space

# The voronoi partition is composed of all vertexes and their relationships
class partition:
	# The core part of the partition is a set of vertexes
	def __init__(self, base, vertex, **kwargs):
		# Limit the number of intersection combinations we can cache
		self.intercachelimit = kwargDef("intercachelimit", kwargs, 100)
		# Base ID map
		self.base = base
		# Map for easy vertex relationships
		self.intermap = {}
		self.intermap[1] = {}
		self.intermap[2] = {}
		# Cells have one point in common, faces have two
		self.cellmap = self.intermap[1]
		self.facemap = self.intermap[2]
		# Add all the verticies
		self.vertex = []
		for vrt in vertex: self.addVertex(vrt)
	# Methods to extract base wave information
	def baseIDs(self): return [b for b in self.base]
	def baseCenter(self, wid): return self.base[wid].center[0]		
	def baseWeight(self, wid): return self.base[wid].weight

	# Cache cells and faces
	def addVertex(self, vertex):
		self.vertex.append(vertex)
		self.interAllCache(vertex)

	# Cache all combinations for this vertex
	def interAllCache(self, vertex):
		for N in range(1, len(vertex.parent)):
			if not self.interNCache(vertex, N): return False
		return True

	# Create a map for every combination of points. Within reason.
	def interNCache(self, vertex, N):
		parent = vertex.parent
		if not N in self.intermap: self.intermap[N] = {}
		numComb = 0
		# Use a generator so our cache limit applies during computation
		for faceID in combinations(parent,N):
			numComb += 1
			if numComb > self.intercachelimit: return False
			if faceID in self.intermap[N]:
				self.intermap[N].append(vertex)
			else:
				self.intermap[N] = [vertex]
		return True
		
	# Faces and cells are pre-computed
	def createFace(self, headID, tailID):
		pass
	def createCell(self, cellID):
		pass
	def createGraph(self, cellID):
		pass

# Whether a time is within a span
def inSpan(T, span, eps = EPS):
	return (((span[0] == None) or (T + EPS) > span[0]) and
		((span[1] == None) or (T - EPS) < span[1]))

# Clamp a wave's C between extremes
def clampSpan(wave, T, eps = EPS):
	if (T + EPS) < wave.span[0]:
		return [0.0], True
	if (T - EPS) > wave.span[1]:
		# TODO: potential bug here
		# TODO: fix by clamping spans in windowing, not selection
		#print("Potential bug! Need better span clamps...")
		return [magnitude(vec(wave.center[0], wave.center[1]))], True
	return None, False

def timeSpan(wave, eps = EPS, bound = BOUND):
	# Find the zeros of the wave intersection
	time = sorted(rzero(wave))
	time = [t for t in time if t < bound]
	# If the waves do not intersect
	if len(time) == 0: return (None, None), False
	# Return the span(s)
	if len(time) == 1: return (time[0], None), True
	return (time[0], time[1]), True

def centerSpan(wave, time, eps = EPS, bound = BOUND):
	if time[0] == None: return (None, None), False
	L = lambda T: wave.L(T, clamp = False)
	c0 = L(time[0])
	c1 = None if time[1] == None else L(time[1])
	return (c0, c1), True

# The timespans during which a wave is an exterior one that can intersect
def exteriorSpan(WF, wave):
	inf, sup = wave.span
	# Check breakouts for infimum calculation
	if wave in WF.breakout and len(WF.breakout[wave]) > 0:
		breakout = WF.breakout[wave]
		inf = max([i.span[0] for i in breakout])
	# Check breakins for supremum calculation
	if wave in WF.breakin and len(WF.breakin[wave]) > 0:
		breakin = WF.breakin[wave]
		sup = min([i.span[1] for i in breakin])
	return (inf, sup)

def linear(a, b, eps = EPS):
	p = lambda T: (-1.0 * b(T)) / a(T)
	cP = lambda T: 0.0 if abs(a(T)) < eps else p(T)
	return cP

def quadratic(a, b, c, eps = EPS):
	left = lambda T: -1.0 * b(T)
	right = lambda T: (b(T) ** 2.0 - (4.0 * a(T) * c(T))) ** 0.5
	denom = lambda T: 2.0 * a(T)
	p = lambda T: (left(T) + right(T)) / denom(T)
	n = lambda T: (left(T) - right(T)) / denom(T)
	rP = lambda T: p(T) if numpy.isreal(p(T)) else T
	rN = lambda T: n(T) if numpy.isreal(n(T)) else T
	cP = lambda T: linear(b,c)(T) if abs(a(T)) < eps else rP(T)
	cN = lambda T: linear(b,c)(T) if abs(a(T)) < eps else rN(T)
	return (cN,cP)

# Clip away times where the wave is not intersecting the plane of intersection
def clipWindow(p, d, window):
	if d == None: return window, False
	projectionDifference = lambda T: p.R(T) - magnitude(d(T))
	x0,x1 = window
	head, tail = projectionDifference(x0), projectionDifference(x1)
	if head < 0.0 and tail < 0.0: return None, True
	hasClip = head * tail < 0.0
	sign = head < 0.0
	window = (x0,x1)
	if hasClip:
		clipLocation = solve(projectionDifference, a = x0, b = x1)
		window = (clipLocation, x1) if sign else (x0, clipLocation)
	return window, hasClip

# Each window is the location of a potential minimum
def window(p1, p2, d1, d2, bound = BOUND):
	inf1 = p1.span[0]
	sup1 = bound if p1.span[1] == None else p1.span[1]
	inf2 = p2.span[0]
	sup2 = bound if p2.span[1] == None else p2.span[1]
	inf = max([inf1, inf2])
	sup = min([sup1, sup2])
	m1 = None if p1.span[1] == None else (p1.span[0] + p1.span[1]) / 2.0
	m2 = None if p2.span[1] == None else (p2.span[0] + p2.span[1]) / 2.0
	span = sorted([x for x in [inf,sup,m1,m2] if not x == None])
	clipped = []
	for x in range(len(span) - 1):
		clip1,has1 = clipWindow(p1,d1,(span[x],span[x+1]))
		clip2,has2 = clipWindow(p2,d2,(span[x],span[x+1]))
		if clip1 == None or clip2 == None: continue
		win = [max([clip1[0],clip2[0]]),min([clip1[1],clip2[1]])]
		# if has1 or has2:
			# print("CLIP:",x,span[x],span[x+1],clip1,clip2,win)
		clipped += [win]
	return clipped

def projectPOI(radius, poiNormal):
	M = lambda v: magnitude(v)
	sign = lambda t, R, D: 1.0 if R(t) > M(D(t)) else -1.0
	Rp = lambda t, R, D: sign(t,R,D) * (abs(R(t)**2.0 - M(D(t))**2.0)**0.5)
	projL = lambda t, R, D: R(t) if D == None else Rp(t, R, D)
	return lambda t: projL(t, radius, poiNormal)

#Zero when when r1 = c1, such as in 1d when abs(r1 - r2) = I or r1 + r2 = I
def rzero(wave):
	# Initialize local names for variables
	p1, p2 = wave.parents()
	d1, d2 = wave.D
	span = window(p1, p2, d1, d2)
	# Project waves onto the intersection subspace
	r1p = projectPOI(p1.R, d1)
	r2p = projectPOI(p2.R, d2)
	# Find point of interection between p1 and p2
	I = lambda t: magnitude(wave.I(t))
	# When the plane of intersection is 1D
	n1 = lambda t: r1p(t) - r2p(t) - I(t)
	n2 = lambda t: r2p(t) - r1p(t) - I(t)
	p = lambda t: r1p(t) + r2p(t) - I(t)
	func = [p, n1, n2]

	# For each possible root do root finding
	roots = []
	for inf,sup in span:
		for fun in func:
			if(fun(inf) * fun(sup)) < 0.0:
				root = solve(fun, a = inf, b = sup)
				roots.append(root)
	return roots

def waveEq(w1, w2):
	ret = True
	for h,t in zip(w1.center[0], w2.center[0]):
		ret = ret and (abs(h - t) < 0.0001)
	return ret

def hasParent(w1,w2):
	if w1 == w2: return True
	if w1.leaf(): return False
	ret = False
	for p in w1.parents(): ret = ret or hasParent(p,w2)
	return ret

def negWID(wid1,wid2):
	ret = []
	change = False
	for wid in wid2:
		if wid in wid1:
			change = True
			continue
		ret += [wid]
	return change, ret

def subWID(wid1, wid2):
	# All parts of id1 must be in id2
	for wid in wid1:
		if not wid in wid2:
			return False
	return True

def siblingWID(wid1, wid2):
	for wid in wid1:
		if wid in wid2: return True
	return False

# Siblings share a wave and have already intersected
def waveSibling(WF, w1, w2):
	# If their IDs share a number
	wid1, wid2 = waveID(WF, w1), waveID(WF,w2)
	return siblingWID(wid1, wid2)

# Wave cousins are when a 3rd party wave shares these waves as sub-waves
# They are candidates for interior waves that don't intersect
def waveCousin(WF, parent, w1, w2):
	wid = waveID(WF, parent)
	wid1, wid2 = waveID(WF, w1), waveID(WF, w2)
	sub1, sub2 = subWID(wid, wid1), subWID(wid, wid2)
	print(wid, wid1, wid2, sub1, sub2)
	return sub1 and sub2

# TODO: efficiency
def interiorCandidates(WF, wav):
	w1, w2 = wav.parents()
	wid1, wid2 = waveID(WF, w1), waveID(WF, w2)
	neg1 = []
	neg2 = []
	for wid in [WF.ids[x] for x in WF.wave]:
		ch1, rem1 = negWID(wid1, wid)
		ch2, rem2 = negWID(wid2, wid)
		# Interiors are always relative to a single base wave
		if ch1 and len(rem1) == 1: neg1 += [rem1]
		if ch2 and len(rem2) == 1: neg2 += [rem2]
	union = []
	for neg in neg1:
		if neg in neg2: union.append(tuple(neg))
	return union

# Debug vectors to print
class debugvec:
	def __init__(self, span, vec):
		self.span = span
		self.vec = vec

# Identify a wave uniquely in the wavefront
def waveID(WF, wav):
	if wav == None: return [-1]
	if wav in WF.ids.keys(): return WF.ids[wav]
	if wav.leaf(): return [-1]
	p1, p2 = wav.parents()
	return tuple(sorted(waveID(WF,p1) + waveID(WF,p2)))

def waveID2(WF, w1, w2):
	ret = [w for w in waveID(WF,w1)]
	for w in waveID(WF,w2):
		if not w in ret: ret.append(w)
	return sorted(ret)

# return if two waves are equivalent
def waveIDEQ(id1,id2):
	if not len(id1) == len(id2): return False
	for x,y in zip(id1,id2):
		if not x == y: return False
	return True

def waveIDIN(overcache, wid):
	for wid2 in overcache:
		if waveIDEQ(wid, wid2): return True
	return False

# A wavefront is composed of many waves
# TODO: encapsulation. Wavefronts need bounds they are provably correct in
# TODO: encapsulation. By proving that "addWave" always correctly adds a wave
# TODO: within the current bounds (by finding further verticies) we can do this 
class wavefront:
	def __init__(self, points, weights, ids={}):
		# Metadata for debug
		self.debugMD = False
		# Base data
		self.point = points
		self.weight = weights
		self.alias = {}
		self.dim = None
		self.base = []
		self.wave = []
		self.ended = []
		# Keep track of wave IDs for debug and readability purposes
		self.debug = []
		self.ids = {}
		self.rev = {}
		self.nextID = 0
		# There are no start verticies
		self.vertex = []
		# Keep track of interior / exterior crossovers
		self.breakout = {}
		self.breakin = {}
		# Keep track of children for purposes of wave endings
		self.children = {}
		# Add the base waves from input data
		self.start(ids)
	def setDebug(self,flag = True): self.debugMD = flag
	def start(self,ids={}):
		self.wave = []
		num = 0
		for point, weight in zip(self.point, self.weight):
			# Ensure points are equidimensional
			if self.dim == None: self.dim = len(point)
			if not len(point) == self.dim: continue
			# Calculate alias if needed
			alias = None
			if num in ids: alias = ids[num]
			num += 1
			# Add the wave
			self.addWave(baseWave(point,weight,self.dim),alias)
	def N(self): return self.dim
	# Bind a wave ID to a wave
	def bindID(self, wav, wid):
		self.ids[wav] = tuple(wid)
		self.rev[tuple(wid)] = wav
		return

	# Return the next wid for the given wave
	def ensureWID(self, wav):
		wid = None
		if wav.leaf():
			wid = tuple([self.nextID])
			self.nextID += 1
		else:
			wid = waveID(self, wav)
		self.bindID(wav, wid)

	# Ensure that the wave has correct child relations
	def ensureChildren(self, wav):
		if not wav in self.children:
			self.children[wav] = []
		if wav.leaf(): return
		for cand in self.wave:
			if wav in self.children[cand]: continue
			candWID = waveID(self, cand)
			waveWID = waveID(self, wav)
			if len(candWID) == len(waveWID): continue
			if subWID(candWID, waveWID):
				self.children[cand].append(wav)
		return

	# Ensure that the wave has the correct relationships
	def ensureRelations(self, wav):
		# Ensure the wave has a WID
		self.ensureWID(wav)
		# Ensure the wave has right child relationships
		self.ensureChildren(wav)
		# Ensure default breakin and breakout waves
		if not wav in self.breakin: self.breakin[wav] = []
		if not wav in self.breakout: self.breakout[wav] = []
		# Leaf waves have a special cache
		if wav.leaf(): self.base.append(wav)
		
	# Ensure that the given wave has correct breakout behavior
	def ensureBreakout(self, wav):
		if wav.leaf(): return
		# Interior calculations TODO: cache, clean
		# TODO: move outside of this function
		if not wav in self.breakout:
			self.breakout[wav] = []
		interiorCand = interiorCandidates(self, wav)
		for interior in interiorCand:
			# Store breakout wave(s) when confirmed
			if self.checkInterior(wav, interior):
				if(self.debugMD): print("CHECKED")
				out = self.rev[interior]
				BW = wave(wav, out, self.dim)
				self.breakout[wav].append(BW)
	# Add the specified vertex to the wavefront
	def addVertex(self, vert):
		self.ensureWID(vert)
		self.ensureChildren(vert)
		wid = waveID(self, vert)
		if(self.debugMD): print("ADDING VERTEX", wid)
		self.vertex.append(waveVertexThunk(self, vert))
	# Add the specified wave to the wavefront, return if it was a success
	# TODO: encapsulation: this needs to guarantee a proper waveform
	def addWave(self,wav,alias = None):
		if wav == None: return False
		# Debug vectors
		if not wav.debug == None:
			for vec in wav.debug:
				self.debug.append(debugvec(wav.span,vec))
		# Store alias if needed
		if not alias == None: self.alias[wav] = alias
		# Waves that barely span space have special relationships
		if wav.N() == -1: return self.addVertex(wav)
		# Otherwise the wave has standard relationships
		self.ensureRelations(wav)
		self.wave.append(wav)
		return True
	# TODO: encapsulation to assure good waveform
	def endWave(self, wav):
		# Base waves do not end
		if wav.leaf(): raise Exception("Tried to end a base wave")
		self.ended.append(wav)
		# Store the wave vertex associated with the end
		self.addVertex(wav)
		# If this is a 1-form mark the smaller wave as interior
		p1, p2 = wav.parents()
		if p1.leaf() and p2.leaf():
			if p1.weight > p2.weight:
				self.breakin[p2].append(wav)
			else:
				self.breakin[p1].append(wav)
			return
		# Remove the wave from parents' children and mark lonelies
		lonely = []
		for cand in self.wave:
			if wav in self.children[cand]:
				self.children[cand].remove(wav)
				if len(self.children[cand]) == 0:
					lonely.append(cand)
		# If parents are loney, mark them as breakin if needed
		for loner in lonely:
			if len(self.breakout[loner]) > 0: continue
			self.breakin[loner].append(wav)
		
	def widList(self): return [x for x in self.ids.values()]
	def popWave(self): self.wave.pop()
	# Return all waves active at t = T
	def cut(self, T):
		ret = []
		for wav in self.wave:
			s = wav.span
			if T > s[0] and (s[1] == None or T < s[1]):
				ret.append(wav)
		return ret
	# Return all verticies in the currenet wavefront, including remnants
	def getVerticies(self):
		remnants = []
		for wav in self.wave:
			if wav.leaf(): continue
			if wav in self.ended: continue
			remnants.append(waveVertexThunk(self, wav))
		return self.vertex + remnants
	# Calculate if a wave is interior to another
	def checkInterior(self, wav, interiorID):
		# Just get the wave center and compare it to the wid
		start = wav.span[0]
		interWave = self.rev[interiorID]
		disvec = vec(wav.center[0], interWave.L(start))
		dismag = magnitude(disvec)
		# The wave is interior if the center lies within the outer wave
		return dismag < interWave.R(start)
	# TODO: need massive speedups (Don't compare all waves baka!)
	# TODO: encapsulation. This should be a function over wavefronts.
	# TODO: encapsulation: there should only be an "addNext" function?
	def nextWaveMerge(self, eps = EPS):
		col = None
		minT = 0.0
		# Cache to prevent symmetric overlapp
		overcache = []
		# Check for wave mergings
		for w1 in self.base:
			for w2 in self.wave:
				# ID candidates to prevent overlapp.
				wid = waveID2(self, w1, w2)
				if waveIDIN(overcache, wid): continue
				if waveIDIN(self.widList(), wid): continue
				overcache.append(wid)
				# Wave siblings have already intersected
				if waveSibling(self, w1, w2): continue
				# Waves need to span space to intersect
				if w1.N() + w2.N() < self.dim - 1: continue
				# Perform wave merge calculations
				merge = wave(w1, w2, self.dim)
				# Perform breakout calculations
				self.ensureBreakout(merge)
				# An invalid merge only happens in rare cases
				if not merge.valid(): continue
				# Find the next event, either merge or join
				T = exteriorSpan(self, merge)[0]
				if col == None or T < minT:
					col = merge
					minT = T
		return col, minT
		
	def nextWaveEnd(self, eps = EPS):
		end = None
		minT = 0.0
		# Check for (unmarked) wave endings
		for w in self.wave:
			if w in self.ended: continue
			T = w.span[1]
			if T == None: continue
			if end == None or T < minT:
				end = w
				minT = T
		return end, minT
		
	def nextEvent(self, eps = EPS):
		waveEnd, eT = self.nextWaveEnd(eps)
		waveMerge, mT = self.nextWaveMerge(eps)
		if(self.debugMD):
			print("END CANDIDATE", eT, waveID(self, waveEnd))
			print("MERGE CANDIDATE", mT, waveID(self,waveMerge))
		# Events are just a wave and what to do (add or remove)
		endEvent = (waveEnd, lambda x: self.endWave(x))
		mergeEvent = (waveMerge, lambda x: self.addWave(x))
		# Wave ends are prioritized over new merges in edge cases
		if waveEnd == None and waveMerge == None: return None, None
		if waveEnd == None: return mergeEvent
		if waveMerge == None: return endEvent
		return endEvent if (eT - eps) < mT else mergeEvent

# TODO: move these into libraries
# '''''''''''
# VECTOR MATH
# '''''''''''
def vec(head, tail):
	return [t - h for h,t in zip(head, tail)]

def hat(N, index):
	ret = [0.0] * N
	ret[index] = 1.0
	return ret

def norm(vector, eps = EPS):
	if magnitude(vector) < eps: return vector
	return scale(vector, 1.0 / magnitude(vector))

def magnitude(vector):
	return dot(vector, vector) ** 0.5

def scale(vector, scale):
	if vector == None: return vector
	return [x * scale for x in vector]

def toSize(vector, size):
	return scale(norm(vector), size)

def dis(head, tail):
	return magnitude(vec(head, tail))

def vsum(v1, v2):
	if v1 == None and v2 == None: return None
	if v1 == None: return [x for x in v2]
	if v2 == None: return [x for x in v1]
	return [h + t for h,t in zip(v1, v2)]

def dot(head, tail):
	dot = 0.0
	for x,y in zip(head,tail):
		dot += x * y
	return dot

def angle(head, tail):
	cosang = cosangle(head,tail)
	if cosang > 1.0:
		print("Cosangle Innacuracy", cosang)
		cosang = 1.0
	if cosang < -1.0:
		print("Cosangle Inaccuracy", cosang)
		cosang = -1.0
	return math.acos(cosang)

def cosangle(head, tail, eps = EPS):
	if magnitude(head) * magnitude(tail) < eps: return 0.0
	return dot(head, tail) / (magnitude(head) * magnitude(tail))

# '''''''''''''
# SUBSPACE MATH
# '''''''''''''
# Project a subspace onto another subspace TODO: linear dependence?
# TODO: for higher dimensions we need orthonormal basis, otherwise the projection will also lead to a scaling operation
def project(basis, vector):
	if basis == None or len(basis) == 0 or len(vector) == 0: return vector
	ret = []
	for vec in vector:
		proj = None
		for b in basis:
			base = norm(b)
			#TODO: Does this work?
			mag = magnitude(vec)
			comp = scale(base, mag * cosangle(base, vec))
			proj = vsum(comp, proj)
		ret.append(proj)
	return ret

# Reject from a subspace onto another subspace
def reject(basis, vector):
	if basis == None or len(basis) == 0 or len(vector) == 0: return vector
	projection = project(basis, vector)
	return [vec(p,v) for p,v in zip(projection, vector)]

# Find the union of two (sub)spaces
def span(s1, s2, eps = EPS): return clean(s1 + s2, eps)


if __name__ == "__main__":
	# The points for the voronoi cells!
	points = [[1.0,1.0],[1.0,0.5],[0.5,1.0]]
	weights = [1.0,2.0,3.0]

	#Initalize the wavefront, set debug mode
	WF = wavefront(points, weights)
	WF.setDebug()

	#Initialize the waves
	wav = []
	for point, weight in zip(points, weights):
		wav.append(baseWave(point, weight,len(points[0])))

	# Run the Algorithm
	partition = runWevoNd(WF,0.001)
	# Return success
	sys.exit(0)
