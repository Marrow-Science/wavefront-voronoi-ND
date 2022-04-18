import vis
import voronoi
import sys
import argparse

def createArgumentParser():
	des = '''
		A wavefront algorithm for computing Nd voronoi cells.
	'''
	parser = argparse.ArgumentParser(description=des)
	parser.add_argument(
		"-b", "--binary",
		action='store_true',
		help="set binary input/output mode")
	parser.add_argument(
		"-v", "--vis",
		action='store_true',
		help="visualize with panda3d")
	parser.add_argument(
		"-i", "--input",
		nargs='?',
		type=argparse.FileType('r'),
		default=sys.stdin)
	parser.add_argument(
		"-o", "--output",
		nargs='?',
		type=argparse.FileType('w'),
		default=sys.stdout)
	return parser

class humanInput:
	def __init__(self, fileI = None):
		self.len = 0
		self.points = []
		self.weights = []
		self.ids = {}
		self.readFromFile(fileI)

	def toWavefront(self):
		print(self.points,self.weights,self.ids)
		return voronoi.wavefront(self.points, self.weights, self.ids)

	def readFromFile(self, fileI):
		if fileI == None: return True
		for line in fileI: self.parseLine(line)
		return True

	def parseLine(self, line):
		# Check for 'w=' and '!'
		newID = self.len
		weight = 1.0
		point = []
		idmatch = lambda t: t[0] == '!'
		weightmatch = lambda t: t[0] == 'w' and t[1] == '='
		for tok in line.strip().split():
			if idmatch(tok): newID = str(tok[1:])
			elif weightmatch(tok): weight = float(tok[2:])
			else: point.append(float(tok))
		self.points.append(point)
		self.weights.append(weight)
		if not newID == self.len: self.ids[self.len] = newID
		self.len += 1

class humanOutput:
	def __init__(self, partition):
		self.raw = partition
	def writeToFile(self, fileO):
		pass

def parse(parser, args):
	try: return parser.parse_args(args)
	except: sys.exit(2)

def readInput(fileI, binary):
	try: return readBinary(fileI) if binary else readHumanReadable(fileI)
	except Exception as e:
		print("Exception:", e)
		sys.exit(2)

def writeOutput(fileO, result, binary):
	wB = writeBinary
	wHR = writeHumanReadable
	try: return wB(fileO, result) if binary else wHR(fileO, result)
	except:
		print("Exception:", e)
		sys.exit(2)

def readHumanReadable(fileI):
	# Human readable has stuff and things
	read = humanInput(fileI)
	return read

def writeHumanReadable(fileO, result):
	pass

def readBinary(fileI):
	raise RuntimeError('Binary I/O function(s) not implemented')

def writeBinary(file0):
	raise RuntimeError('Binary I/O function(s) not implemented')

def wavefrontSetup(read): return read.toWavefront()

def visualize(wavefront, result):
	raise RuntimeError('Visualization not integrated')

if __name__ == "__main__":
	# Create a parser for command line arguments
	parser = createArgumentParser()
	# Allright, command line arguments parsed!
	parsed = parse(parser,sys.argv[1:])
	# Now to read the input from file
	readinput = readInput(parsed.input,parsed.binary)
	# Pass the input to the wavefront algorithm
	wavefront = wavefrontSetup(readinput)
	# Run the algorithm
	partition = voronoi.runWevoNd(wavefront)
	# Do visualization if wanted
	if parsed.vis: visualize(wavefront, partition)
	# Pass the result to the output file
	writeOutput(parsed.output, partition, parsed.binary)
	# We are done!
	sys.exit(0)
