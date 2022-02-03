# Local Library Modules
import voronoi
# Math modules
import numpy
# Graphics modules
#from tkinter import *
from panda3d.core import GeomVertexFormat, GeomVertexData, GeomVertexWriter
from panda3d.core import GeomTriangles
from panda3d.core import Geom, GeomNode, ModelNode
from panda3d.core import TransparencyAttrib
from panda3d.core import Vec3, Quat
from direct.gui.DirectGui import DirectSlider

# The "default" runtime system
from direct.showbase.ShowBase import ShowBase
from direct.showbase import DirectObject
from direct.task import Task

size = 1200
height = size
width = size

# Smooth a geometry of triangles over a geodesic midpoint equation
def smooth(vert, face, midpoint, itr = 1):
	if itr == 0: return vert, face
	newvert = []
	newface = []
	for f in face:
		place = len(vert + newvert)
		i0, i1, i2 = f
		i3, i4, i5 = range(place, place + 3)
		v0, v1, v2 = (vert[i0], vert[i1], vert[i2])
		v3 = midpoint(v0,v1)
		v4 = midpoint(v1,v2)
		v5 = midpoint(v2,v0)
		newvert += [v3,v4,v5]
		newface += [
				[i0,i5,i3],
				[i1,i3,i4],
				[i2,i4,i5],
				[i3,i5,i4]
			]
	return smooth(vert + newvert, newface, midpoint, itr - 1)

def translate(vert, point):
	return [voronoi.vsum(v,point) for v in vert]

def sphere2midpoint(x, y, r):
	return voronoi.toSize(voronoi.vsum(x,y),r)

def torus2midpoint(x, y, R, r):
	RX = voronoi.toSize((x[0],x[1],0.0),R)
	RY = voronoi.toSize((y[0],y[1],0.0),R)
	rX, rY = (voronoi.reject([RX],[x])[0], voronoi.reject([RY],[y])[0])
	nR = voronoi.toSize(voronoi.vsum(RX,RY),R)
	nr = voronoi.toSize(voronoi.vsum(rX,rY),r)
	return voronoi.vsum(nR, nr)

# Create a hexagon in the given (hyper-) plane
def createHexagon(r, form):
	sin60 = (3 ** 0.5) / 2.0
	cos60 = 0.5
	# The two flat points
	flat = [voronoi.toSize(form[0], r), voronoi.toSize(form[0], -1.0 * r)]
	# The four other points
	rect = []
	for x,y in [(1.0,1.0),(-1.0,1.0),(-1.0,-1.0),(1.0,-1.0)]:
		f0 = voronoi.toSize(form[0], x * cos60 * r)
		f1 = voronoi.toSize(form[1], y * sin60 * r)
		point = voronoi.vsum(f0,f1)
		rect.append(point)
	# Return the points in counterclockwise order
	return [flat[0]] + rect[:2] + [flat[1]] + rect[2:]
	
# Create a hexagonal torus which perserves equilateral triangles
def createHexagonalTorus(R,r):
	formBase = [voronoi.hat(3,0), voronoi.hat(3,1)]
	hexbase = createHexagon(R, formBase)
	# Create hexagons from each base point
	zh = voronoi.hat(3,2)
	allhex = []
	for h in hexbase:
		allhex += [voronoi.vsum(h,v) for v in createHexagon(r,[h,zh])]
	# Join hexagons with weaving triangles
	tri = []
	for x in range(6):
		for y in range(6):
			X,Y = (x * 6, y)
			tailX = ((x + 1) % 6) * 6
			tailY = (y + 1) % 6
			face = (X + Y, tailX + y, X + tailY)
			comp = (X + tailY, tailX + y, tailX + tailY)
			tri.append(face)
			tri.append(comp)
	rev = [(z,y,x) for x,y,z in tri]		
	return allhex, tri + rev

# Creates an icosahedron with circumradius R
def createIcosahedron(r):
	base = []
	GR = (1 + (5 ** 0.5)) / 2.0  # The golden ratio
	RS = r / ((GR + 2.0) ** 0.5) # The unit icosa has weird circumradius
	for y in [-1.0,1.0]:
		for z in [-1.0,1.0]:
			base.append([0.0, y * RS, z * GR * RS])
	VRT = []
	# Create vertexes from the base GR rectangles
	for x in range(12):
		new = [[b[(x+0)%3],b[(x+1)%3],b[(x+2)%3]] for b in base]
		VRT += new

	# Create the 20 triangles, using 6 pennants and 8 individuals
	base = (0,1,3)
	TRI = []

	# Pennants
	TRI += [(0,2,8),(0,2,9),(1,3,10),(1,3,11)]
	TRI += [(4,6,0),(4,6,1),(5,7,2),(5,7,3)]
	TRI += [(8,10,4),(8,10,5),(9,11,6),(9,11,7)]

	# Individuals
	TRI += [(0,4,8),(0,6,9),(1,4,10),(1,6,11),
		(2,5,8),(2,7,9),(3,7,11),(3,5,10)]

	# Need to reverse so the triangles become double sided
	rev = [(z,y,x) for x,y,z in TRI]
	TRI += rev

	return VRT,TRI

def create0Sphere(r):
	# Start with two small spheres to simulate points
	posv, posf, _ = create2Sphere(0.1)
	negv, negf, _ = create2Sphere(0.1)
	# Move them to r radius
	# posv = translate(posv[0],(r, 0.0, 0.0))
	# negv = translate(negv[0],(-1.0 * r, 0.0, 0.0))
	# Need to translate the indices as well
	# shift = len(posv)
	# negf = [(f[0] + shift, f[1] + shift, f[2] + shift) for f in negf]
	return posv + negv, posf + negf, 2

def create1Sphere(R):
	r = 0.03
	# Start with a hexagonal torus
	ivrt, iface = createHexagonalTorus(R,r)
	# Make a geodescic midpoint function
	eq = lambda x,y: torus2midpoint(x,y,R,r)
	# Smooth the hextorus into a real torus
	vert, face = smooth(ivrt, iface, eq, 3)
	return [vert], [face], 1

def create2Sphere(r):
	# Start with an icosahedron
	ivrt, iface = createIcosahedron(r)
	# Make a geodesic midpoint function
	eq = lambda x,y: sphere2midpoint(x,y,r)
	# Smooth the icosahedron into a 2sphere
	vert, face = smooth(ivrt, iface, eq, 3)
	return [vert], [face], 1

# Create an N-sphere of radius 1.0 on the subspace comp
def createSphereGeometry(comp, N, tag = 'def'):
	mesh = [create0Sphere, create1Sphere, create2Sphere]
	name = ['0Sphere', '1Sphere', '2Sphere']
	# Colors: Blue, Orange, Purple
	color = [(0.5,0.1,0.8,0.2),(0.8,0.4,0.1,0.2),(0.0,0.0,0.7,0.2)]
	smvrt, smfce, parts = mesh[N](comp)
	# Put into a geometry thingy
	P = range(parts)
	form = [GeomVertexFormat.get_v3c4() for x in P]
	vdat = [GeomVertexData(name[N] + tag, f, Geom.UHStatic) for f in form]
	for v,s in zip(vdat, smvrt): v.setNumRows(len(s))
	# Add vertexes
	wvrt = [GeomVertexWriter(dat, 'vertex') for dat in vdat]
	wcol = [GeomVertexWriter(dat, 'color') for dat in vdat]
	ptri = [GeomTriangles(Geom.UHStatic) for x in P]
	# Add vertexes
	for s, v, c in zip(smvrt, wvrt, wcol):
		for x, y, z in s:
			v.addData3(x,y,z)
			c.addData4(color[N])
	# Add faces
	for p, f in zip(ptri, smfce):
		for a,b,c in f: p.addVertices(a,b,c)
	# Return the mesh
	geom = [Geom(d) for d in vdat]
	for g, p in zip(geom, ptri): g.addPrimitive(p)
	return geom

def createVector(num):
	rv, rf = [], []
	for x in range(num):
		vert, face, _ = create2Sphere(0.02)
		rv += vert
		rf += face
	return rv, rf

# Create a line of evenly spaced spheres
def createVectorGeometry(tag = 'def'):
	parts = 60
	smvrt, smfce = createVector(parts)
	# SIGH
	P = range(parts)
	form = [GeomVertexFormat.get_v3c4() for x in P]
	vdat = [GeomVertexData('debug' + tag, f, Geom.UHStatic) for f in form]
	for v,s in zip(vdat, smvrt): v.setNumRows(len(s))
	# Add vertexes
	wvrt = [GeomVertexWriter(dat, 'vertex') for dat in vdat]
	wcol = [GeomVertexWriter(dat, 'color') for dat in vdat]
	ptri = [GeomTriangles(Geom.UHStatic) for x in P]
	# Color, a dark red
	color = (0.5,0.1,0.1,0.8)
	# Add vertexes
	for s, v, c in zip(smvrt, wvrt, wcol):
		for x, y, z in s:
			v.addData3(x,y,z)
			c.addData4(color)
	# Add faces
	for p, f in zip(ptri, smfce):
		for a,b,c in f: p.addVertices(a,b,c)
	# Return the mesh
	geom = [Geom(d) for d in vdat]
	for g, p in zip(geom, ptri): g.addPrimitive(p)
	return geom

class VoronoiVis(ShowBase):
	def __init__(self, WF):
		ShowBase.__init__(self)
		# Store the wavefront here
		self.wavefront = WF
		# Initiate hotkeys
		self.init_GUI()
		self.init_camera()
		self.init_keybinds()
		# Initialize the geometry
		self.geom0 = createSphereGeometry(1.0,0)
		self.geom1 = createSphereGeometry(1.0,1)
		self.geom2 = createSphereGeometry(1.0,2)
		self.geom = [self.geom0, self.geom1, self.geom2]
		self.vector = createVectorGeometry()
		self.debug = []
		self.wave = []

	def init_GUI(self):
		self.T = DirectSlider(
			range=(0,10.0),
			value = 5.0,
			pageSize = 0.0003,
			pos = (0.0,0.0,0.9),
			command = self.show)
		self.T.resetFrameSize()
		
	def init_camera(self):
		keymap = base.win.get_keyboard_map()
		self.key = {
			# Camera Keys
			"w": keymap.get_mapped_button('w'),
			"s": keymap.get_mapped_button('s'),
			"a": keymap.get_mapped_button('a'),
			"d": keymap.get_mapped_button('d'),
			"q": keymap.get_mapped_button('q'),
			"e": keymap.get_mapped_button('e'),
			"u": keymap.get_mapped_button('u'),
			"j": keymap.get_mapped_button('j'),
			"h": keymap.get_mapped_button('h'),
			"k": keymap.get_mapped_button('k'),
			"alt": keymap.get_mapped_button("lalt"),
		}
		self.taskMgr.add(self.camera_task, 'CameraControl')
		self.disableMouse()
	
	def init_keybinds(self):
		self.accept('n', self.addWave)
		self.accept('b', self.removeWave)
		self.accept('i', self.isolateWave)
		self.accept('p', self.printWave)

	def camera_task(self,task):
		isDown = self.mouseWatcherNode.is_button_down
		speed = 20.0
		pos = speed * globalClock.dt
		neg = -1.0 * speed * globalClock.dt
		c = self.camera
		M = lambda i: c.getNetTransform().getMat().getRow3(i)
		D = lambda s,i: c.setPos(c.getPos() + M(i) * s)
		R = lambda h,p,r: c.setHpr(c.getH()+h, c.getP()+p, c.getR()+r)
		if isDown(self.key["alt"]): pos,neg = (pos * 0.1, neg * 0.1)
		if isDown(self.key["q"]): D(pos,1)
		if isDown(self.key["e"]): D(neg,1)
		if isDown(self.key["w"]): R(0.0,pos * 2.0,0.0)
		if isDown(self.key["s"]): R(0.0,neg * 2.0,0.0)
		if isDown(self.key["a"]): R(pos * 2.0,0.0,0.0)
		if isDown(self.key["d"]): R(neg * 2.0,0.0,0.0)
		return Task.cont

	def addWave(self):
		nw = self.wavefront.nextCollision()
		if not nw == None: self.wavefront.addWave(nw)
		self.updateWavefront()
		for i,vec in enumerate(debugVectors(nw)):
			handle = self.debugNode(-1 * i)
			self.debug.append((vec, handle))
		debugPrint(nw)
		self.show()

	def removeWave(self):
		self.wavefront.popWave()
		self.updateWavefront()
		self.show()

	def isolateWave(self):
		w = self.wave[-1][0]
		tree = []
		nxt = [w]
		while len(nxt) > 0:
			nw = nxt.pop()
			if not nw.leaf():
				p1,p2 = nw.parents()
				nxt += [p1,p2]
			tree.append(nw)
		# Update wave visibility
		for w,h in self.wave:
			if w in tree: continue
			else: h.hide()
		self.show()

	def printWave(self):
		w = self.wave[-1][0]
		T = self.T['value']
		print("-------- WAVE DIAGNOSTIC -------")
		print("Wave", w)
		print("C:", w.C(T, clamp = False))
		print("L:", w.L(T, clamp = False))
		print("------ END WAVE DIAGNOSTIC -----")
		
		pass

	def show(self):
		#print(self.T['value'])
		self.setTime(self.T['value'])
	
	def clearWaveCache(self):
		for w,h in self.wave:
			h.removeNode()
		for v,h in self.debug:
			h.removeNode()
		self.wave = []
		self.debug = []

	def updateWavefront(self):
		wavefront = self.wavefront
		self.clearWaveCache()
		N = self.wavefront.N()
		# Visualize the wavefronts
		for i,wave in enumerate(wavefront.wave):
			# Render the node
			geom = self.geom[wave.N()]
			node = ModelNode('wave' + str(i))
			for n,g in enumerate(geom):
				geomNode = GeomNode('geom' + str(n))
				geomNode.addGeom(g)
				node.addChild(geomNode)
			handle = self.render.attachNewNode(node)
			handle.setTransparency(TransparencyAttrib.MDual)
			# Store the actor
			self.wave.append((wave,handle))
		# Register the debug vectors
		for i,vec in enumerate(wavefront.debug):
			handle = self.debugNode(i)
			self.debug.append((vec, handle))

	def debugNode(self, i):
		node = ModelNode('debug'+str(i))
		for n,g in enumerate(self.vector):
			geomNode = GeomNode('geom' + str(n))
			geomNode.addGeom(g)
			node.addChild(geomNode)
		handle = self.render.attachNewNode(node)
		handle.setTransparency(TransparencyAttrib.MDual)
		return handle

	def vecplace(self, handle, vec, T):
		place = vec[0](T)
		vec = [0.0]*len(place) if vec[1] == None else vec[1](T)
		#print("placeVec", vec)
		siz = len(handle.children)
		for i,c in enumerate(handle.children):
			per = (1.0 * i) / siz
			pos = voronoi.vsum(place, voronoi.scale(vec,per))
			c.setPos(pos[0],pos[1],pos[2])
			c.setScale(1.0 + (2.0 * per))

	def scale(self, handle, N, r):
		if N == 2:
			handle.setSz(r)
			handle.setSy(r)
			handle.setSx(r)
		if N == 1:
			handle.setSy(r)
			handle.setSx(r)
		if N == 0:
			# TODO: why?
			#r /= 2.0
			handle.getChild(0).setX(r)
			handle.getChild(1).setX(-1.0 * r)


	# Morph the wavefronts to the time
	def setTime(self, T):
		#print("TIME", T)
		# Debug vectors
		for v,h in self.debug:
			if not voronoi.inSpan(T, v.span):
				h.detachNode()
				continue
			h.reparentTo(self.render)
			self.vecplace(h, v.vec, T)
		# Wave placement
		for w,h in self.wave:
			if not voronoi.inSpan(T, w.span):
				h.detachNode()
				continue
			h.reparentTo(self.render)
			N = w.N()
			center = w.L(T)
			p1, p2 = (None,None)
			if not w.leaf(): p1, p2 = w.parents()
			# Set the center position of the form
			cX, cY, cZ = center
			h.setPos(cX, cY, cZ)
			if w.form.N() == 0:
				self.scale(h,N,w.R(T))
			# Set the rotations of the geom
			if w.form.N() == 1:
				look = voronoi.vsum(center, w.P(T)[0])
				lX, lY, lZ = look
				# Use Look At to set a base orientation
				h.lookAt(lX, lY, lZ)
				h.setP(h.getP() + 90.0)
				self.scale(h,N,w.R(T))
			# HACK: I give up, will just set location not orient
			# TODO: make a general vis algorithm
			if w.form.N() == 2:
				compV = w.form.comp(3)[0]
				compVp = voronoi.toSize(compV, w.R(T))
				compVn = voronoi.toSize(compV, -1.0 * w.R(T))
				pos = voronoi.vsum(center, compVp)
				neg = voronoi.vsum(center, compVn)
				h.getChild(0).setPos(toV3(compVp))
				h.getChild(1).setPos(toV3(compVn))

def toV3(basic):
	x, y, z = basic
	return Vec3(x,y,z)

def planeAngle(plane, idx):
	basis = [voronoi.hat(3,x) for x in range(3)][:idx+1:-1]
	prj = voronoi.project(basis, plane)
	return numpy.rad2deg(voronoi.angle(prj[idx], voronoi.hat(3,idx)))

# The canvas goes from -4.0 to 4.0
def canvCoords(coords, r = 0.01):
	center = [(size / 2) + x * (size / 8) for x in coords]
	R = (size / 8) * r
	return [c - R for c in center] + [c + R for c in center] 

def drawwaves(root, canvas, scale, WF):

	canvas.delete("all")

	time = scale.get() / 500.0

	print("TIME", time)

	# Print the centers, X's, and R's for all active waves
	for wave in WF.cut(time):
		center = wave.center[0]
		canvas.create_oval(canvCoords(center), fill = "red")

		delR = wave.R(time)
		canvas.create_oval(canvCoords(center, delR))

	col = WF.nextCollision()
	T, eT = col.span
	c, eC = col.center

	#print(T, eT)

	if eT == None: eT = 10.0

	for x in range(400):
		dT = x / 140.0
		if dT > T and dT < eT:
			delX = voronoi.toSize(col.form.X, col.C(dT)[0])
			delR = col.R(dT)
			C = voronoi.vsum(col.p1.center[0], delX)
			O = col.form.comp(2)
			pF = voronoi.vsum(C, voronoi.scale(O[0], delR))
			pB = voronoi.vsum(C, voronoi.scale(O[0], -1.0 * delR))
			canvas.create_oval(canvCoords(pF), fill = "black")
			canvas.create_oval(canvCoords(pB), fill = "black")

	if time > T and time < eT:
		for c in col.center:
			if c == None: continue
			canvas.create_oval(canvCoords(c), fill = "red")
		delX = voronoi.toSize(col.form.X, col.C(time)[0])
		delR = col.R(time)
		C = voronoi.vsum(col.p1.center[0], delX)
		O = col.form.comp(2)
		pF = voronoi.vsum(C, voronoi.scale(O[0], delR))
		pB = voronoi.vsum(C, voronoi.scale(O[0], -1.0 * delR))
		canvas.create_oval(canvCoords(C), fill = "green")
		canvas.create_oval(canvCoords(pF), fill = "orange")
		canvas.create_oval(canvCoords(pB), fill = "orange")

	root.after(20, lambda: drawwaves(root, canvas, scale, WF))


def printwave(wave):
	if wave == None:
		print("NO WAVE!")
		return
	if not wave.leaf():
		print("PARENTS")
		for p in wave.parents():
			print("---")
			printwave(p)
		print("---")
	print("DIM", wave.N())
	print("CENTER", wave.center)
	print("SPAN", wave.span)

def debugVectors(nw):
	ret = []
	lamcen = lambda T: nw.L(T)
	lamvec = lambda T: voronoi.scale(nw.P(T)[0], nw.C(T)[0])
	ret.append((lamcen,lamvec))
	ret.append((lamcen,lambda T: nw.P(T)[0]))
	#ret.append((lamcen,lambda T: nw.form.comp(3)[0]))
	return [voronoi.debugvec(nw.span, x) for x in ret]

def waveTreePrint(wave,idx):
	print("\t"*idx, wave, wave.span)
	if wave.leaf(): return
	p1, p2 = wave.parents()
	waveTreePrint(p1, idx + 1)
	waveTreePrint(p2, idx + 1)

def debugPrint(nw):
	print("-------- DEBUG PRINT ---------")
	p1, p2 = nw.parents()
	print("PARENT N:",p1.N(), p2.N())
	waveTreePrint(nw,1)
	print("------ END DEBUG PRINT -------")

if __name__ == "__main__":

	'''
	data = [
	(0.0732764158707837, 0.48585568272536983, 0.5781378249494682),
	(0.3939588778405376, 0.7757835671885716, 0.08201036150890328),
	(0.9042783856470057, 0.6139615023676342, 0.03324267155880689),
	(0.48039019498094027, 0.7366494230107247, 0.16819997100066375),
	(0.7805904208855332, 0.8857318080514044, 0.5330266843849021),
	(0.3968644499216717, 0.06525707960187233, 0.9648408682413959),
	(0.4666398805949532, 0.7113641279878318, 0.2547308933969078)]
	data = [[x*10,y*10,z*10] for x,y,z in data]
	weight = [1.0, 2.0, 3.0, 4.0, 1.0, 0.5, 1.5]
	'''

	#data = [[2.0,3.0,0.0],[-1.0,0.0,2.0],[0.0,2.0,-1.0],[-1.0,1.0,3.0]]
	#weight = [1.0, 2.5, 2.0, 2.1]
	data = [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0],[5.0,0.0,5.0]]
	weight = [2.5, 2.0, 2.0, 2.2]
	#data = [[0.0,2.0,-1.0],[-1.0,1.0,3.0]]
	#weight = [2.0, 2.0]


	WF = voronoi.wavefront(data, weight)
	root = VoronoiVis(WF)
	root.setTime(0.5)
	root.run()
