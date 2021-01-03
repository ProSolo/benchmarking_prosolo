import sys
import os
import subprocess
import math

EXEC = "./perfectphy" + "_mintree"

class Result:
    def __init__(self):
        self.numSmallest = 0
        self.numTrials = 0
        
    def add(self, sizes):
        sizes = map(int,sizes.split())
        
        if sizes[0] == sizes[1]:
            self.numSmallest = self.numSmallest + 1
        
        self.numTrials = self.numTrials + 1
    
    def val(self):
        return "%d/%d" % (self.numSmallest, self.numTrials)

def printResults(outf, resultmatrix, nmRange, kRange, TRIALS):
    outf.write("Testing KW Tree Size vs Smallest Minimal Tree Size\nk,");
    outf.write(','.join(str(k) for k in kRange))
    
    for n,m in nmRange:
        outf.write( "\n%d %d" % (n,m) )
        for k in kRange:
            entry = resultmatrix[ (n,m,k) ]
            outf.write( ", " + entry.val() )
    outf.write("\n")

def runTests():
    TRIALS = int(sys.argv[1])
    PARALLEL = int(sys.argv[2])
    
    if TRIALS < PARALLEL or TRIALS % PARALLEL != 0:
        print PARALLEL, "does not divide", TRIALS, "evenly"
        return
    
    for f in (EXEC, "ms", "multextract.pl"):
        if not os.path.isfile(f):
            print 'Cannot find file "%s" in current directory' % f
            return
    
    nmRange = [ (x,x) for x in [10, 50, 100, 500, 1000, 2000] ]
    kRange = [4, 10, 20]
    
    # create matrix of results for every (N,M,K) triplet
    resultmatrix = dict( [ ((n,m,k), Result()) for n,m in nmRange for k in kRange ] )
    
    try:
        for k in kRange:
            for n,m in nmRange:
                params = (n,m,k)
                
                # generate and interpret datasets
                gendata = subprocess.Popen("./ms %d %d -s %d -r 0.0 5000 > mybigdata ; perl multextract.pl mybigdata 0 %d > /dev/null ; rm -r o%dstate* > /dev/null" % (n, TRIALS, m * (k-1), k, k), shell = True )
                
                print "(%d,%d,%d) Generating..." % params,
                sys.stdout.flush()
                
                entry = resultmatrix[ params ]
                
                assert(0 == gendata.wait())
                
                print "Running...",
                sys.stdout.flush()
                
                for t in xrange(1,TRIALS+1,PARALLEL):
                    
                    # start PARALLEL instances of the program  
                    kws = [ subprocess.Popen("%s -f %dstate%d.%d.%d -enum" % (EXEC, k, n, m, ext),
                            executable = "/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True) for ext in xrange(t,t+PARALLEL) ]
                    
                    # as they finish, assert the nonexistence of errors
                    #for ext in xrange(t,t+PARALLEL):
                    #    assert 0 == os.wait()[1]
                    #    
                    #    print ext,
                    #    sys.stdout.flush()
                    
                    # all are finished, thus extract and record runtimes
                    for kw in kws:
                        resultoutput, erroroutput = kw.communicate()
                        
                        try:
                            sizes = resultoutput.splitlines()[2]
                            entry.add( sizes )
                        except Exception as e:
                            print "Error:", [c for c in resultoutput], [c for c in erroroutput]
                            print e
                
                subprocess.Popen("rm -r *state*", shell = True).wait()
                
                print "(%s/%s small)" % (entry.numSmallest, TRIALS)
    finally:
        subprocess.Popen("rm -r *state* mybigdata* *datalist seedms", shell = True).wait()
        
        outf = open("results_mintree.txt", 'a')
        
        outf.write("Iterations = %s\n" % sys.argv[1])
        
        printResults(outf, resultmatrix, nmRange, kRange, TRIALS)
        
        outf.write("\n\n")
        
        outf.close()
        
        print "Results saved"

if __name__ == '__main__':
    runTests()
