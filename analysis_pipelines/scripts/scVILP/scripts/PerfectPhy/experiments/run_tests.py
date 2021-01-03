import sys
import os
import subprocess

EXEC = "./perfectphy"

class Result:
    def __init__(self):
        self.count = 0
        self.ms = 0
    
    def add(self, usert):
        self.count += 1
        
        # usert = '14.460'
        part = usert.partition('.')
        self.ms += 1000 * int(part[0])
        self.ms += int(part[2])
    
    def avg(self):
        if self.count > 0:
            return int(self.ms / self.count)
        else:
            return -1

def printResults(outf, resultmatrix, nmRange, kRange, TRIALS):
    outf.write("Testing KW\nk,");
    
    outf.write(','.join(str(k) for k in kRange))
    
    for n,m in nmRange:
        outf.write( "\n%d %d" % (n,m) )
        
        for k in kRange:
            entry = resultmatrix[ (n,m,k) ]
            if entry.count == TRIALS:
                outf.write(",%d" % entry.avg() )
            else:
                outf.write(",-1")
    
    outf.write("\n\n")


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
                    kws = [ subprocess.Popen("TIMEFORMAT='%%3U' ; time %s -f %dstate%d.%d.%d" % (EXEC, k, n, m, ext),
                            executable = "/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True) for ext in xrange(t,t+PARALLEL) ]
                    
                    # as they finish, assert the nonexistence of errors
                    for ext in xrange(t,t+PARALLEL):
                        assert 0 == os.wait()[1]
                        
                        print ext,
                        sys.stdout.flush()
                    
                    # all are finished, thus extract and record runtimes
                    for kw in kws:
                        resultoutput, timeoutput = kw.communicate()
                        
                        assert resultoutput[0] == '1'
                        
                        # timeoutput format:
                        # timeoutput[0] = 'xxx.xxx'
                        
                        entry.add(timeoutput)
                
                subprocess.Popen("rm -r *state*", shell = True).wait()
                
                print "(%d ms)" % entry.avg()
    finally:
        subprocess.Popen("rm -r *state* mybigdata* *datalist seedms", shell = True).wait()
        
        outf = open("results.txt", 'a')
        
        outf.write("Iterations = %s\n" % sys.argv[1])
        
        printResults(outf, resultmatrix, nmRange, kRange, TRIALS)
        
        outf.write("\n\n")
        
        outf.close()
        
        print "Results saved"

if __name__ == '__main__':
    runTests()
