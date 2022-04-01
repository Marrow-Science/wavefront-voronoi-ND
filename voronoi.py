import math
import numpy
import scipy
#from scipy.optimize import root_scalar as solve
from scipy.optimize import brentq as solve
from numpy.polynomial import Polynomial as poly

# We create a weighted Voronoi cell in ND using a wavefront algorithm

# Globals are bad, but hey...
EPS = 0.001
BOUND = 100.0

# The form determines what the (hyper)plane of intersection is, and its
# orthogonal directions and weights
def formDir(dirs, N):
	ret = [0.0] * N
	for x in dirs: ret = vsum(ret, x)
	if magnitude(ret) < EPS: return ret
	return scale(ret, 1.0 / magnitude(ret))

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
		while len(comp) < compN and i < N:
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
			return self.interInitNxN(p1, p2, c1, c2, f1, f2)

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
	# Only case with D vectors TODO: polish
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

	# Full case of two directed wavefronts
	# In this case the intersection plane is 2D
	# The POI intersects both parent centers
	def interInitNxN(self, p1, p2, c1, c2, f1, f2):
		# Calculate I
		x1p = lambda T: p1.L(T, clamp = False)
		x2p = lambda T: p2.L(T, clamp = False)
		I = lambda T: vec(x1p(T), x2p(T))
		Ir = lambda T: vec(x2p(T), x1p(T))
		# The POI is defined by I and parent forms
		orthP = interPlane(f1, f2)
		Ip = lambda T: project(orthP.basis,[I(T)])[0]
		H = lambda T: reject([Ip(T)],[orthP.X()])[0]
		P = lambda T: subspace([I(T),H(T)])
		# Calculate projections
		x1 = lambda T: norm(project(P(T).basis,[f1.X])[0])
		x2 = lambda T: norm(project(P(T).basis,[f2.X])[0])
		d1, d2 = None, None
		# Calculate raw angles (exterior or interior)
		t1a = lambda T: angle(x1(T), I(T))
		t2a = lambda T: angle(x2(T), Ir(T))
		# Ensure the angle is interior
		hC = lambda th, T: th(T) < (math.pi / 2.0)
		rC = lambda th, T: math.pi - th(T)
		interior = lambda th, T: rC(th,T) if hC(th,T) else th(T)
		# Calculate final angle lambdas
		t1 = lambda T: interior(t1a,T) - (math.pi / 2.0)
		t2 = lambda T: interior(t2a,T) - (math.pi / 2.0)
		Th = lambda T: math.pi - (t1(T) + t2(T))
		# Parameterize C2 in terms of C1
		sin = lambda th: math.sin(th)
		C2L = (lambda T: 0.0, lambda T: sin(t1(T)) / sin(t2(T)))
		# Store the plane of intersection TODO: polish
		Pl = lambda T: P(T).vec
		# Return the setup information
		self.x1 = lambda T: x1(T)
		self.x2 = lambda T: x2(T)
		self.v1 = lambda T: P(T).basis[0]
		self.v2 = lambda T: P(T).basis[1]
		self.t1a = t1a
		self.t2a = t2a
		self.t1 = t1
		self.t2 = t2
		self.Ip = Ip
		self.H = H
		return Pl, Th, I, (d1, d2), C2L

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
		# Solve the quadratic (or linear) equation for two spheres
		I = lambda t: magnitude(self.I(t))
		solution = []
		# Linear equation (0x0, 0xN)
		if self.C2L == None:
			f0 = lambda t: 2.0 * I(t)
			f1 = lambda t: r2p(t)**2.0 - I(t)**2.0 - r1p(t)**2.0
			solution = [linear(f0,f1)]
		# Quadratic equation (NxN)
		else:
			cos = lambda t: math.cos(t)
			c2l = lambda x, t: self.C2L[x](t)
			f0 = lambda t: 2.0 * (1.0 - cos(theta(t)) * c2l(1,t))
			f1 = lambda t: -2.0 * cos(theta(t)) * c2l(0,t)
			f2 = lambda t: r2p(t)**2.0 - I(t)**2.0 - r1p(t)**2.0
			n, p = quadratic(f0, f1, f2)
			# Use the solution that matches I
			sign = lambda t: dot(self.I(t),self.P(t)[0])
			signed = lambda t: n(t) if sign(t) < 0.0 else p(t)
			# TODO: need to use negative when appropriate?
			#solution = [n,p]
			solution = [signed]
		# We only use the smallest solution TODO: possible NxN bug here
		sol = sorted([s(T) for s in solution])
		return sol

	# The current location of the wave center
	def L(self, T, clamp = True):
		d = None if self.D[0] == None else self.D[0](T)
		C = scale(self.P(T)[0], self.C(T, clamp)[0])
		#return vsum(self.p1.center[0], vsum(d, C))
		ret = vsum(self.p1.L(T, clamp), vsum(d, C))
		return ret

	# The radius perpendicular to the form at time T
	def R(self, T):
		# Find projected C and R
		Cs = self.C(T)[0]
		Rp = projectPOI(self.parents()[0].R,self.D[0])(T)
		# Find radius
		radius = abs((Rp ** 2.0) - (Cs ** 2.0)) ** 0.5
		#print("RAD -- C -- Rp", radius, Cs, Rp)
		return radius

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
		print("Potential bug! Need better span clamps...")
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

#Zero when r1 + r2 = I or abs(r1 - r2) = I! (Prev: when r1 = c1)
def rzero(wave):
	# Initialize local names for variables
	p1, p2 = wave.parents()
	d1, d2 = wave.D
	span = window(p1, p2, d1, d2)
	# Project waves onto the intersection subspace
	r1p = projectPOI(p1.R, d1)
	r2p = projectPOI(p2.R, d2)
	# Find point of interection between p1 and p2
	p,n = (None, None)
	fA, fB = (None, None)
	I = lambda t: magnitude(wave.I(t))
	# When the plane of intersection is 1D
	if wave.C2L == None:
		n = lambda t: abs(r1p(t) - r2p(t)) - I(t)
		p = lambda t: abs(r1p(t) + r2p(t)) - I(t)
	# When the POI is 2D
	else:
		mx = lambda t: r1p(t)**2.0 + r2p(t)**2.0
		cosR = lambda t: 2.0 * r1p(t) * r2p(t) * math.cos(wave.T(t))
		n = lambda t: mx(t) - cosR(t) - (I(t)**2.0)
		p = None

	# For each possible root do root finding
	roots = []
	for inf,sup in span:
		r1, r2 = None, None
		if n(inf) * n(sup) < 0.0:
			r1 = solve(n, a = inf, b = sup)
			roots.append(r1)
		if not p == None and p(inf) * p(sup) < 0.0:
			r2 = solve(p, a = inf, b = sup)
			roots.append(r2)
	#TODO: remove debug statements
	if p1.N() + p2.N() == 2:
		print("PARSPAN", p1.span, p2.span)
		print("ROOTS!",roots, inf, sup, n(inf), n(sup))
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

def waveSibling(wave, w1, w2):
	# If their IDs share a number
	id1, id2 = waveID(wave, w1), waveID(wave,w2)
	isSibling = False
	for x in id1:
		if x in id2:
			isSibling = True
			break
	return isSibling
	# TODO: is this old way the correct one?
	'''
	for w in [x for x in wave if not x.leaf()]:
		l, r = w.parents()
		if l == w1 and r == w2: return True
		if r == w1 and l == w2: return True
	return hasParent(w1,w2) or hasParent(w2,w1)
	'''

# Debug vectors to print
class debugvec:
	def __init__(self, span, vec):
		self.span = span
		self.vec = vec

# Identify a wave uniquely in the wavefront
def waveID(WF, wave):
	if wave in WF.ids.keys(): return WF.ids[wave]
	if wave.leaf(): return [-1]
	return sorted(waveID(WF,wave.p1) + waveID(WF,wave.p2))

def waveID2(WF, w1, w2):
	return sorted(waveID(WF,w1) + waveID(WF,w2))

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
class wavefront:
	def __init__(self, points, weights):
		self.point = points
		self.weight = weights
		self.dim = len(points[0])
		self.wave = []
		self.debug = []
		# Keep track of wave IDs for debug purposes
		self.ids = {}
		self.nextID = 0
		# Add the base waves from input data
		self.start()
	def start(self):
		self.wave = []
		for point, weight in zip(self.point, self.weight):
			self.addWave(baseWave(point,weight,self.dim))
	def N(self): return self.dim
	# Add the specified wave to the wavefront, return if it was a success
	def addWave(self,wave):
		if wave == None: return False
		if not wave.debug == None:
			for vec in wave.debug:
				self.debug.append(debugvec(wave.span,vec))
		self.wave.append(wave)
		if wave.leaf(): 
			self.ids[wave] = [self.nextID]
			self.nextID += 1
		else:
			self.ids[wave] = waveID(self,wave)
		return True
	def widList(self): return [x for x in self.ids.values()]
	def popWave(self): self.wave.pop()
	# Return all waves active at t = T
	def cut(self, T):
		ret = []
		for wave in self.wave:
			s = wave.span
			if T > s[0] and (s[1] == None or T < s[1]):
				ret.append(wave)
		return ret
	# TODO: need massive speedups (Don't compare all waves baka!)
	def nextCollision(self):
		col = None
		minT = 0.0
		endT = 0.0
		overcache = []
		for w1 in self.wave:
			for w2 in self.wave:
				# ID candidates to prevent overlapp.
				wid = waveID2(self, w1, w2)
				if waveIDIN(overcache, wid): continue
				if waveIDIN(self.widList(), wid): continue
				else:
					special = (w1.N() == 2 and w2.N() == 0) or (w1.N() == 0 and w2.N() == 2)
					if not special:
						overcache.append(wid)
				# Wave siblings have already intersected
				if waveSibling(self, w1, w2): continue
				# PRINT WAVE PLUS
				# Waves need to span space to intersect
				if w1.N() + w2.N() < self.dim - 1: continue
				# Waves that barely span space make a vertex
				if w1.N() + w2.N() == self.dim - 1:
					print(w1.N(), w2.N())
					if not (w1.N() == 1 and w2.N() == 1):
						continue
					# TODO: Vertex intersections!
					print("POTVRT:", w1.N(), w2.N())
					print("WIDS",self.ids[w1],self.ids[w2])
					m = wave(w1, w2, self.dim)
					print("WAVE DIMS", w1.N(), w2.N(), m.N())
					inf, sup = m.span
					p1,p2 = m.parents()
					r1p = projectPOI(p1.R, m.D[0])
					r2p = projectPOI(p2.R, m.D[1])
					theta = m.theta
					cos = math.cos
					print("Span inf sup", inf, sup)
					print("D", m.D)
					if not inf == None:
						print("Th:",m.theta(inf))
						print("t1:",m.t1(inf))
						print("t2:",m.t2(inf))
						print("t1a:",m.t1a(inf))
						print("t2a:",m.t2a(inf))
					if not m.valid():
						print("INVALID!")
					else:
						print("VRT:",m,m.span)
					print(waveID(self,p1), waveID(self,p2))
					# Add debug vectors to parent(s)
					# TODO: SIGH this is still borked >->
					s = m.span
					s1 = p1.span
					s2 = p2.span
					c1 = lambda t: p1.L(t)
					c2 = lambda t: p2.L(t)
					l1 = lambda t: vec(c1(t),c2(t))
					l2 = lambda t: vec(c2(t),c1(t))
					cD = lambda t: vsum(c1(t),scale(l1(t),0.5))
					x1 = m.x1
					x2 = m.x2
					self.debug = []
					# X's
					self.debug.append(debugvec(s1,(c1,x1)))
					self.debug.append(debugvec(s2,(c2,x2)))
					# I debug
					self.debug.append(debugvec(s1,(c1,l1)))
					self.debug.append(debugvec(s1,(cD,m.v1)))
					self.debug.append(debugvec(s1,(cD,m.v2)))
					self.debug.append(debugvec(s1,(cD,m.H)))
					continue
				merge = wave(w1, w2, self.dim)
				if not merge.valid(): continue
				T = merge.span[0]
				if col == None or T < minT:
					col = merge
					minT = T
		return col

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

	#Initalize the wavefront
	WF = wavefront(points, weights)

	#Initialize the waves
	wav = []
	for point, weight in zip(points, weights):
		wav.append(baseWave(point, weight,len(points[0])))

	# Sort intersections by temporal order
	# Advance time to the next step, clean the wavefront by removing
	# unused waves, then repeat until done
	# All steps are parallelizable

	collision = WF.nextCollision()
	print(collision.span, collision.center)
