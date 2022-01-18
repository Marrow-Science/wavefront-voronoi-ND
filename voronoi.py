import math
import numpy
import scipy
#from scipy.optimize import root_scalar as solve
from scipy.optimize import brentq as solve
from numpy.polynomial import Polynomial as poly

# We create a weighted Voronoi cell in ND using a wavefront algorithm

# Globals are bad, but hey...
EPS = 0.001
BOUND = 10.0

# The form determines what the (hyper)plane of intersection is, and its
# orthogonal directions and weights
def formDir(dirs, N):
	ret = [0.0] * N
	for x in dirs: ret = vsum(ret, x)
	if magnitude(ret) < EPS: return ret
	return scale(ret, 1.0 / magnitude(ret))

class form:
	def __init__(self, center, I, p = None):
		if center == None:
			print("NO CENTER!", base, center)
			return None
		#print("FORM INIT")
		self.base = [center] if p == None else p[0].base + p[1].base
		wvec = [vec(b, center) for b in self.base]
		#print("wvec:", wvec)
		v = None
		for x in wvec: v = vsum(x, v)
		# TODO: THIS SIGN IS FUCKING UP!
		self.X = norm(v)
		#self.X = scale(norm(v), -1.0)
		#I = project(P, [self.X])
		#print("I", I)
		#print("X:", self.base, center, self.X)
		form = I if p == None else span([I] + p[0].dir, p[1].dir)
		self.dir = form
		self.dim = len(self.dir)

	def N(self): return self.dim
	# To make the orthogonal dual do a series of rejections
	# Rejections maintain orthogonality over orthogonal bases
	def comp(self, compN):
		N = self.N()
		comp = []
		i = 0
		while len(comp) < compN and i < N:
			nxt = reject(self.dir + comp, [hat(N,i)])[0]; i += 1
			if magnitude(nxt) < 0.01: continue
			comp.append(nxt)
		return comp		

def interPlane(f1, f2):
	print("INTERPLANE",f1.N(), f2.N())
	x1, x2 = (f1.X, f2.X)
	r1, r2 = (reject([x1],[x2]), reject([x2],[x1]))
	return span(r1, r2)

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
		self.debug = None
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
	def theta(self): return None
	def C(self, T): return [0.0]
	def L(self, T): return self.center[0]
	def R(self, T): return self.weight * T

class wave:
	# Basic initialization from a pair of intersecting parents
	def __init__(self, p1, p2, N):
		# TODO: BETTER DEBUG
		self.debug = None
		# TODO: BETTER DEBUG
		self.p1 = p1
		self.p2 = p2
		self.P, self.T, self.I, self.D, self.C2L = self.interInit(N)
		self.span, self.center, self.val = timeSpan(self)
		#TODO: BETTER DEBUG
		#TODO: BETTER DEBUG
		# TODO: split into function
		if self.val:
			print("CENTER!", self.span, self.center)
			start = self.span[0]
			center = self.center[0]
			I = [vec(p1.x1p(start), center)]
			I = project(self.P(start), [self.I(start)])
			#I = self.P(self.span[0])
			I = norm(self.I(self.span[0]))
			self.form = form(center, I, (p1.form, p2.form))
			print("BASE", self.form.base, self.center)
			self.dim = N - self.form.N() - 1
			#print("FORM DIMS:", N, self.form.N(), self.dim)
			if not(self.dim == p1.dim-1 or self.dim == p2.dim-1):
				self.val = False
			# TODO: DEBUG
			#if self.form.N() == 2:
			#	self.debug = [(lambda T: p2.x1p(self.span[0]), lambda T: I[0])]
			# TODO: DEBUG
		else: self.dim = N - 1
		# TODO: REMOVE DEBUGS
		if self.form.N() == 2:
			pC = lambda T: self.p1.center[0]
			D = None if self.D[0] == None else self.D[0]
			dp = lambda T: vsum(pC(T), D(T))
			Cvec = lambda T: scale(self.P(T)[0], self.C(T)[0])
			#self.debug += [(dp,Cvec)]
			self.debug += [(pC,D)]
		# TODO: REMOVE DEBUGS
	# Init various intersection values needed to calculate the waveform
	def interInit(self, N, eps = EPS):
		p1, p2 = self.parents()
		c1, c2 = (None, None)
		#c1 = p1.center[0] if p1.leaf() else p1.parents()[0].center[0]
		#c2 = p2.center[0] if p2.leaf() else p2.parents()[0].center[0]
		c1, c2 = (p1.center[0], p2.center[0])
		f1, f2 = (p1.form, p2.form)
		left, right = (f1.N() == 0, f2.N() == 0)
		#self.x1p = lambda T: vsum(c1, scale(f1.X, p1.C(T)[0]))
		#self.x2p = lambda T: vsum(c2, scale(f2.X, p2.C(T)[0]))
		self.x1p = lambda T: p1.L(T)
		self.x2p = lambda T: p2.L(T)
		if left and right:
			return self.interInit0x0(c1, c2)
		elif left or right:
			return self.interInit0xN(p1, p2, c1, c2, f1, f2)
		else:
			return self.interInitNxN(p1, p2, c1, c2, f1, f2)

	# Base degenerate case of two undirected wavefronts
	def interInit0x0(self, c1, c2):
		theta = 0.0
		Iv = lambda T: vec(c1, c2)
		Is = lambda T: magnitude(Iv(T))
		P = lambda T: [norm(vec(c1, c2))]
		C2L = (	(lambda T: Is(T), lambda T: -1.0),
			(lambda T: Is(T), lambda T: 1.0))
		D = (None, None)
		# TODO: BETTER DEBUG
		#self.debug = [(lambda T: c1, lambda T:vec(c1,c2))]
		# TODO:
		return P, theta, Iv, D, C2L

	# Partial case of one directed and one undirected wavefront
	def interInit0xN(self, p1, p2, c1, c2, f1, f2):
		# Cases, which one is 0?
		L, R = (f1.N() == 0, f2.N() == 0)
		# Xp
		x1p, x2p = (self.x1p, self.x2p)
		# I
		i1v = lambda T: vec(x1p(T),x2p(T))
		i2v = lambda T: vec(x2p(T),x1p(T))
		# F
		f1, f2 = (p1.form.X, p2.form.X)
		# D
		d1 = None if R else lambda T: project([f2], [i1v(T)])[0]
		d2 = None if L else lambda T: project([f1], [i2v(T)])[0]
		D = (d1, d2)
		# Special Case
		F = f1 if R else f2
		B = lambda T: reject([F], [i1v(T)])[0]
		P = lambda T: [norm(B(T))]
		Ip = lambda T: B(T)
		theta = 0.0
		C2L = (	(lambda T: Ip(T), lambda T: -1.0),
			(lambda T: Ip(T), lambda T: 1.0))
		self.debug = []
		#self.debug += [(x1p, i1v)]
		dp1 = lambda T: vsum(x1p(T),d1(T))
		#dp2 = lambda T: vsum(x2p(T),d2(T))
		#self.debug += [(x1p, d1)]
		#self.debug += [(x2p, d2)]
		self.debug += [(dp1, B)]
		#self.debug += [(dp1, lambda T: F)]
		self.debug += [(dp1, lambda T: P(T)[0])]
		return P, theta, Ip, D, C2L

	# Full case of two directed wavefronts
	def interInitNxN(sekf, p1, p2, c1, c2, f1, f2):
		P = lambda T: interPlane(f1, f2)
		x1 = lambda T: toSize(project(P, f1.X)[0], p1.C(T)[0])
		x2 = lambda T: toSize(project(P, f2.X)[0], p2.C(T)[0])
		x1p = lambda T: vsum(c1, x1(T))
		x2p = lambda T: vsum(c2, x2(T))
		I = lambda T: vec(x1p(T), x2p(T))
		d1 = lambda T: reject([P], x1(T))
		d2 = lambda T: reject([P], x2(T))
		c1a = lambda T: angle(I(T), x1(T))
		c2a = lambda T: angle(I(T), x2(T))
		c1t = lambda T: (math.pi / 2.0) - c1a(T)
		c2t = lambda T: (math.pi / 2.0) - c2a(T)
		t1 = angle(vec(c1,c2), project(P, f1.X)[0])
		t2 = angle(vec(c2,c1), project(P, f2.X)[0])
		nT = (t1 + t2)
		Th = math.pi - (t1 + t2)
		sin = math.sin
		C2L = ( (lambda T: 0.0, lambda T: sin(c2a(T)) / sin(c1a(T))),
			(lambda T: 0.0, lambda T: sin(c2t(T)) / sin(c1t(T))))
		return P, Th, I, D, C2L

	# If this is a leaf node
	def leaf(self): return False
	# If the parents actually intersect at some point
	def valid(self): return self.val
	# Return the dimensionality of the wave
	def N(self): return self.dim
	# To make things simpler
	def parents(self): return self.p1, self.p2
	# The theta between centers of the two parent forms
	def theta(self): return self.T
	# The distance along the parent form at time T!
	def C(self, T):
		#if not inspan(T, self.span): return [0.0]
		if (T + EPS) < self.span[0]:
			return [0.0]
		if not self.span[1] == None and (T - EPS) > self.span[1]:
			return [magnitude(vec(self.center[0], self.center[1]))]
		theta = self.theta()
		p1, p2 = self.parents()
		d1, d2 = self.D
		# Project waves onto the intersection subspace TODO: R > D!
		sign = lambda t, R, D: 1.0 if R(t) > magnitude(D(t)) else -1.0
		Rp = lambda t, R, D: abs(R(t)**2.0 - magnitude(D(t))**2.0)**0.5
		p1s = lambda t: 1.0 if d1 == None else sign(t,p1.R,d1)
		p2s = lambda t: 1.0 if d2 == None else sign(t,p2.R,d2)
		Ps = lambda t: p2s(t) if d1 == None else p1s(t)
		r1p = lambda t: p1.R(t) if d1 == None else p1s(t)*Rp(t,p1.R,d1)
		r2p = lambda t: p2.R(t) if d2 == None else p1s(t)*Rp(t,p2.R,d2)
		# Solve the quadratic (or linear) equation for two spheres
		I = lambda t: magnitude(self.I(t))
		solution = []
		# TODO: NEED A WAY TO HANDLE NEGATIVE RADII
		if abs(theta) < EPS:
			f0 = lambda t: 2.0 * I(t)
			f1 = lambda t: r2p(t)**2.0 - I(t)**2.0 - r1p(t)**2.0
			solution = [linear(f0,f1)]
		else:
			f0 = lambda t: 2.0 * (1.0 - math.cos(theta) * self.C2L[0][1](t))
			f1 = lambda t: -2.0 * math.cos(theta) * self.C2L[0][0](t)
			f2 = lambda t: r1p(t)**2.0 - I(t)**2.0 - r2p(t)**2.0
			n, p = quadratic(f0, f1, f2)
			solution = [n,p]
		# Reject them back to the big subspace
		Cr = lambda t, C, D: C(t) if D == None else (C(t)**2.0 + magnitude(D(t))**2.0)**0.5
		#Cr = lambda t, C, D: C(t) if D == None else Ps(t) * (C(t)**2.0 + magnitude(D(t))**2.0)**0.5
		sol = sorted([Cr(T, s, d1) for s in solution])
		return sol
	# The current location of the wave center
	def L(self, T):
		d = None if self.D[0] == None else self.D[0](T)
		C = scale(self.P(T)[0], self.C(T)[0])
		return vsum(self.p1.center[0], vsum(d, C))
	# The radius perpendicular to the form at time T (Just do a thingy)
	def R(self, T):
		p1, p2 = self.parents()
		r1 = p1.R(T)
		c1 = self.C(T)[0]
		#if self.form.N() == 2:
		#	print("R -- C1!",r1, c1)
		#c1 = magnitude(self.selfC(T))
		return abs((r1 ** 2.0) - (c1 ** 2.0)) ** 0.5

# Whether a time is within a span
def inspan(T, span, eps = EPS):
	return (((span[0] == None) or (T + EPS) > span[0]) and
		((span[1] == None) or (T - EPS) < span[1]))

# When r = 0 we got ourselves a start / end point!
def timeSpan(wave, eps = EPS, bound = BOUND):
	time = sorted(rzero(wave))
	time = [t for t in time if t < bound]
	#print(time)
	if len(time) == 0: return (None, None), (None, None), False
	p1, p2 = wave.parents()
	c1, c2 = (p1.center[0], p2.center[0])
	d1, d2 = wave.D
	# Projection onto the intersection plane
	Rp = lambda t, R, D: abs(R(t)**2.0 - magnitude(D(t))**2.0)**0.5
	#cen = lambda T: vsum(wave.x1p(T), toSize(wave.P(T)[0], -1.0 * p1.R(T)))
	Pb = lambda T: c1 if d1 == None else vsum(c1, d1(T))
	Cp = lambda T: p1.R(T) if d1 == None else Rp(T,p1.R,d1)
	# TODO: this is the problem! P(T) does not track correctly! Need B?
	#cen = lambda T: vsum(Pb(T), toSize(wave.P(T)[0], p1.R(T)))
	Cv = lambda T: scale(wave.P(T)[0], Cp(T))
	cen = lambda T: vsum(Pb(T), Cv(T))
	#wave.debug = [(Pb,Cv)]
	#cen = lambda T: vsum(wave.x1p(T), toSize(wave.form.X, wave.C(T)[0]))
	#cen = lambda T: vsum(wave.x1p(T), toSize(wave.P(T)[0], wave.C(T)[0]))
	if len(time) == 1: return (time[0], None), (cen(time[0]), None), True
	return (time[0], time[1]), (cen(time[0]), cen(time[1])), True

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

def clipWindow(p, d, window):
	if d == None: return window, False
	projectionDifference = lambda T: p.R(T) - magnitude(d(T))
	x0,x1 = window
	head, tail = projectionDifference(x0), projectionDifference(x1)
	print("HEAD TAIL", head, tail)
	if head < 0.0 and tail < 0.0: return None, True
	hasClip = head * tail < 0.0
	sign = head < 0.0
	window = (x0,x1)
	if hasClip:
		clipLocation = solve(projectionDifference, a = x0, b = x1)
		print("CLIPLOC:", clipLocation, projectionDifference(clipLocation))
		window = (clipLocation, x1) if sign else (x0, clipLocation)
	return window, hasClip

def window(p1, p2, d1, d2, MAX = 10.0):
	inf1 = p1.span[0]
	sup1 = MAX if p1.span[1] == None else p1.span[1]
	inf2 = p2.span[0]
	sup2 = MAX if p2.span[1] == None else p2.span[1]
	inf = max([inf1, inf2])
	sup = min([sup1, sup2])
	m1 = None if p1.span[1] == None else (p1.span[0] + p1.span[1]) / 2.0
	m2 = None if p2.span[1] == None else (p2.span[0] + p2.span[1]) / 2.0
	span = sorted([x for x in [inf,sup,m1,m2] if not x == None])
	clipped = []
	# TODO: d1 and d2 are a fucking mess, why are they mixed up this way!
	# TODO: FIX IT! (This has been verified to be the correct ordering btw.
	for x in range(len(span) - 1):
		clip1,has1 = clipWindow(p1,d1,(span[x],span[x+1]))
		clip2,has2 = clipWindow(p2,d2,(span[x],span[x+1]))
		if clip1 == None or clip2 == None: continue
		win = [max([clip1[0],clip2[0]]),min([clip1[1],clip2[1]])]
		if has1 or has2:
			print("CLIP:", x, span[x], span[x+1], clip1, clip2, win)
		clipped += [win]
	return clipped

#Zero when r1 + r2 = I or abs(r1 - r2) = I! (Prev: when r1 = c1)
def rzero(wave):
	theta = wave.theta()
	p1, p2 = wave.parents()
	d1, d2 = wave.D
	span = window(p1, p2, d1, d2)
	print("D1,D2:", d1, d2)
	# Project waves onto the intersection subspace
	sign = lambda t, R, D: 1.0 if R(t) > magnitude(D(t)) else -1.0
	Rp = lambda t, R, D: sign(t,R,D) * (abs(R(t)**2.0 - magnitude(D(t))**2.0)**0.5)
	r1p = lambda t: p1.R(t) if d1 == None else Rp(t, p1.R, d1)
	r2p = lambda t: p2.R(t) if d2 == None else Rp(t, p2.R, d2)
	#Cr = lambda t, C, D: C(t) if D == None else (C(t)**2.0 + magnitude(D(t))**2.0)**0.5
	#c1r = lambda C, t: C(t) if d1 == None else Cr(t, C, d1)
	#c2r = lambda C, t: C(t) if d2 == None else Cr(t, C, d2)
	p,n = (None, None)
	fA, fB = (None, None)
	I = lambda t: magnitude(wave.I(t))
	if abs(wave.theta()) < EPS:
		n = lambda t: abs(r1p(t) - r2p(t)) - I(t)
		p = lambda t: abs(r1p(t) + r2p(t)) - I(t)
		
	else:
		Th = wave.theta()
		nT = math.pi - wave.theta() 
		mx = lambda t: r1p(t)**2.0 + r2p(t)**2.0
		cosR = lambda t,th: 2.0 * r1p(t) * r2p(t) * math.cos(th)
		n = lambda t: mx(t) - cosR(t,nT) - (I(t)**2.0)
		p = lambda t: mx(t) - cosR(t,Th) - (I(t)**2.0)
	

	
	'''
	print("RANGE!")
	for x in range(100):
		t = spans[0] + (x * ((spans[-1] - spans[0]) / 100.0))
		if d1 == None and d2 == None: break
		if d1 == None:
			print()
			print(t, "-----")
			print("D2:", magnitude(d2(t)))
			print("R1:", p1.R(t))
			print("R2:", p2.R(t))
			print(sign(t,p1.R,d2), Rp(t,p1.R,d2))
			print(sign(t,p2.R,d2), Rp(t,p2.R,d2))
		if d2 == None:
			print(t, "-----")
			print("D1:", magnitude(d1(t)))
			print("R1:", p1.R(t))
			print("R2:", p2.R(t))
			print(sign(t,p1.R,d1), Rp(t,p1.R,d1))
			print(sign(t,p2.R,d1), Rp(t,p2.R,d1))
		#print(n(t), p(t))
	'''
	print("SPANS!", span)
	
	roots = []
	# For each possible root do root finding
	for inf,sup in span:
		#inf,sup = spans[x], spans[x+1]
		r1, r2 = None, None
		if n(inf) * n(sup) < 0.0:
			r1 = solve(n, a = inf, b = sup)
			roots.append(r1)
		if p(inf) * p(sup) < 0.0:
			r2 = solve(p, a = inf, b = sup)
			roots.append(r2)
		
	print("ROOTS", roots)

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
	for w in [x for x in wave if not x.leaf()]:
		l, r = w.parents()
		if l == w1 and r == w2: return True
		if r == w1 and l == w2: return True
	return hasParent(w1,w2) or hasParent(w2,w1)

# Debug vectors to print
class debugvec:
	def __init__(self, span, vec):
		self.span = span
		self.vec = vec

# A wavefront is composed of many waves
class wavefront:
	def __init__(self, points, weights):
		self.point = points
		self.weight = weights
		self.dim = len(points[0])
		self.wave = []
		self.debug = []
		self.start()
	def start(self):
		self.wave = []
		for point, weight in zip(self.point, self.weight):
			self.wave.append(baseWave(point, weight, self.dim))
	def N(self): return self.dim
	def addWave(self,wave):
		if not wave.debug == None:
			for vec in wave.debug:
				self.debug.append(debugvec(wave.span,vec))
		self.wave.append(wave)
	# Return all waves active at t = T
	def cut(self, T):
		ret = []
		for wave in self.wave:
			s = wave.span
			if T > s[0] and (s[1] == None or T < s[1]):
				ret.append(wave)
		return ret
	# TODO: make this recognize lower-form collisions
	# TODO: need massive speedups (Don't compare all waves baka!)
	def nextCollision(self):
		col = None
		minT = 0.0
		endT = 0.0
		for w1 in self.wave:
			for w2 in self.wave:
				#print(w1.N(),w2.N())
				if waveSibling(self.wave, w1, w2): continue
				if w1.N() + w2.N() < self.dim: continue
				# TODO: REMOVE THIS DEBUG!
				#if w1.N() == 1 or w2.N() == 1: continue
				merge = wave(w1, w2, self.dim)
				print("MERGE",
					w1.N(), w2.N(),
					merge.N(), merge.valid())
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
	return math.acos(cosangle(head, tail))

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
			#print(mag)
			comp = scale(base, mag * cosangle(base, vec))
			proj = vsum(comp, proj)
		ret.append(proj)
	return ret

# Reject from a subspace onto another subspace TODO: linear dependence?
def reject(basis, vector):
	if basis == None or len(basis) == 0 or len(vector) == 0: return vector
	projection = project(basis, vector)
	return [vec(p,v) for p,v in zip(projection, vector)]

# Find the union of two (sub)spaces
def span(s1, s2, eps = EPS): return clean(s1 + s2, eps)

# clean a space by removing linearly dependent vectors
def clean(space, eps = EPS):
	if len(space) < 2: return [norm(s) for s in space if magnitude(s)>EPS]
	ret = []
	for vec in space:
		if magnitude(reject(ret, [norm(vec)])[0]) > eps:
			ret.append(norm(vec))
	return ret

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