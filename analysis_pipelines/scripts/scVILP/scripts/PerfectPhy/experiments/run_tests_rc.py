import sys
import os
import subprocess
import time

class Result:
    def __init__(self):
        self.count = 0
        self.ms = 0
        self.chars = 0
    
    def add(self, usert, c):
        self.count += 1
        self.chars += int(c)
        
        # usert = '14.460'
        part = usert.partition('.')
        self.ms += 1000 * int(part[0])
        self.ms += int(part[2])
    
    def avg(self):
        if self.count > 0:
            return int(self.ms / self.count), int(self.chars / self.count)
        else:
            return -1, -1

def printResults(outf, resultmatrix, nmRange, kRange, TRIALS):
    outf.write("Testing Character Removal\nk,");
    
    outf.write(','.join(str(k) + ' ' + md for k,md in kRange))
    
    for n,m in nmRange:
        outf.write( "\n%d %d" % (n,m) )
        
        for k,md in kRange:
            entry = resultmatrix[ (n,m,k,md) ]
            if entry.count == TRIALS:
                outf.write(",(%d ms, %d ch)" % entry.avg() )
            else:
                outf.write(",-")
    
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
    
    nmRange = [ (50,50) ]#m) for m in xrange(10,41,10) ]
    kRange = [ (4,r) for r in ["0.0"] ] #["0.01", "0.1", "0.25"] ]
    
    # create matrix of results for every (N,M,K) triplet
    resultmatrix = dict( [ ((n,m,k,r), Result()) for n,m in nmRange for k,r in kRange ] )
    
    try:
        for k, r in kRange:
            for n,m in nmRange:
                params = (n,m,k,r)
                entry = resultmatrix[ params ]
                
                #print entry.count
                #continue
                
                while entry.count < TRIALS:
                    # generate and interpret datasets
                    gendata = subprocess.Popen("./ms %d %d -s %d -r %s 5000 > mybigdata ; perl multextract.pl mybigdata 0 %d > /dev/null ; rm -r o%dstate* > /dev/null" % (n, TRIALS, m * (k-1), r, k, k), shell = True )
                    
                    print ""
                    print "(%d,%d,%d,%s) Generating..." % params,
                    sys.stdout.flush()
                    
                    assert(0 == gendata.wait())
                    
                    print "Running...",
                    sys.stdout.flush()
                    
                    for t in xrange(1,TRIALS+1,PARALLEL):
                        
                        # start PARALLEL instances of the program
                        kws = [ subprocess.Popen("TIMEFORMAT='%%3R' ; time python charremoval.py < %dstate%d.%d.%d" % (k, n, m, ext),
                                executable = "/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True) for ext in xrange(t,t+PARALLEL) ]
                        
                        # as they finish, assert the nonexistence of errors
                        for kw in kws:
                                resultoutput, timeoutput = kw.communicate()
                                
                                # resultoutput format:
                                #    1 [num chars] ...
                                # or 0 ...
                                
                                # timeoutput format:
                                # timeoutput[0] = 'xxx.xxx'
                                
                                resultoutput = resultoutput.split(' ')
                                
                                if resultoutput[0] == "1":
                                    if int(resultoutput[1]) < m:
                                        entry.add(timeoutput, resultoutput[1])
                                elif resultoutput[0] == "0":
                                    entry.add(timeoutput, "0")
                                elif resultoutput[0] != "x":
                                    print resultoutput
                                    raise Exception(str(params))
                                
                                print entry.count,
                                sys.stdout.flush()
                    
                    subprocess.Popen("rm -r *state*", shell = True).wait()
                    
                    #print "(%d/%d)" % (entry.count, TRIALS),
                
                print "(%d ms, %d ch)" % entry.avg()
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
