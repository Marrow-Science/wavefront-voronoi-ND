import wevoNd

testdir = "test/"
testfile = [
	"simple2.dat",
	"simple3.dat",
	"simple4.dat"
]

def makeArgs(test):
	# We want to visualize the run
	return ["wevoTest","-v","-i",test]

if __name__ == "__main__":
	for test in testfile:
		wevoNd.main(makeArgs(testdir+test))
