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

class fileInput:
	def __init__(self, points, weights, ids):
		self.points = points

def readHumanReadable(fileI):
	# Human readable has stuff and things
	ret = {'points':[],'weights':[],'ids':{}}
	return ret

def writeHumanReadable(fileO, result):
	pass

def readBinary(fileI):
	raise RuntimeError('Binary I/O function(s) not implemented')

def writeBinary(file0):
	raise RuntimeError('Binary I/O function(s) not implemented')

def wavefrontSetup(read):
	points, weights, ids = read['points'], read['weights'], read['ids']
	return voronoi.wavefront(points, weights, ids=ids)

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
	print(wavefront.wave)
	# Run the algorithm
	partition = voronoi.runWevoNd(wavefront)
	# Do visualization if wanted
	if parsed.vis: visualize(wavefront, partition)
	# Pass the result to the output file
	writeOutput(parsed.output, partition, parsed.binary)
	# We are done!
	sys.exit(0)
