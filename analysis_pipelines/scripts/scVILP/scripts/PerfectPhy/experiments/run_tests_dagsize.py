import sys
import os
import subprocess
import math

EXEC = "./perfectphy"

class Result:
    def __init__(self):
        self.data = []
        
    def add(self, size):
        self.data.append(int(size.split()[0]))
    
    def avg(self):
        if len(self.data) == 0:
            return -1.0
        else:
            return sum(self.data)/float(len(self.data))

def printResults(outf, resultmatrix, nmRange, kRange, TRIALS):
    outf.write("Testing DAG output file size\nk,");
    outf.write(','.join(str(k) for k in kRange))
    
    for n,m in nmRange:
        outf.write( "\n%d %d" % (n,m) )
        for k in kRange:
            entry = resultmatrix[ (n,m,k) ]
            outf.write( ", " + str(entry.data) )
    outf.write("\n")
    
    for n,m in nmRange:
        outf.write( "\n%d %d" % (n,m) )
        for k in kRange:
            entry = resultmatrix[ (n,m,k) ]
            outf.write( ", " + str(entry.avg()) )
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
    
    nmRange = [ (x,x) for x in [10, 50, 100, 500]]#, 1000, 2000] ]
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
                    kws = [ subprocess.Popen("%s -f %dstate%d.%d.%d -enum dag%d" % (EXEC, k, n, m, ext, ext),
                            executable = "/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True) for ext in xrange(t,t+PARALLEL) ]
                    
                    for kw in kws:
                        kw.wait()
                    
                    for ext in xrange(t,t+PARALLEL):
                        try:
                            resultoutput, erroroutput = subprocess.Popen("du -b dag%d" % (ext), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True).communicate()
                            
                            entry.add(resultoutput)
                        except Exception as e:
                            print e
                    
                subprocess.call("rm -r *state* dag*", shell = True)
                
                print "(%s)" % (entry.avg())
    finally:
        subprocess.call("rm -r *state* dag* mybigdata* *datalist seedms", shell = True)
        
        outf = open("results_dagsize.txt", 'a')
        
        outf.write("Iterations = %s\n" % sys.argv[1])
        
        printResults(outf, resultmatrix, nmRange, kRange, TRIALS)
        
        outf.write("\n\n")
        
        outf.close()
        
        print "Results saved"

if __name__ == '__main__':
    runTests()
