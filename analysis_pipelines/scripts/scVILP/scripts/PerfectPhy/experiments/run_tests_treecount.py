import sys
import os
import subprocess
import math

EXEC = "./perfectphy"

def avg(l):
    if len(l) == 0:
        return -1.0
    return sum(l) / float(len(l))

def sd(l):
    if len(l) == 0:
        return -1.0
    if len(l) == 1:
        return 0.0
    
    a = avg(l)
    return math.sqrt( sum( (t - a)**2 for t in l ) / len(l) )

def median(l):
    if len(l) == 0:
        return -1.0
    if len(l) == 1:
        return l[0]
    
    l = sorted(l)
    centerIndex = (len(l) - 1)/2
    
    if len(l) % 2 == 0:
        return (l[centerIndex] + l[centerIndex + 1]) / 2.0
    else:
        return float(l[centerIndex])

class Result:
    def __init__(self):
        self.trees = []
        self.times = []
        
    def add(self, tcount, usert):
        self.trees.append( int(tcount) )
        # usert = '14.460'
        part = usert.partition('.')
        ms = 1000 * int(part[0]) + int(part[2])
        self.times.append(ms)
    
    def avg(self):
        return "(%s, %s ms)" % (avg(self.trees), avg(self.times))

    def sd(self):
        return "(%s, %s ms)" % (sd(self.trees), sd(self.times))
    
    def median(self):
        return "(%s, %s ms)" % (median(self.trees), median(self.times))

def printResults(outf, resultmatrix, nmRange, kRange, TRIALS):
    outf.write("Testing KW Tree Count\nk,");
    outf.write(','.join(str(k) for k in kRange))
    
    for n,m in nmRange:
        outf.write( "\n%d %d" % (n,m) )
        for k in kRange:
            entry = resultmatrix[ (n,m,k) ]
            outf.write( ", " + str(zip(entry.trees, entry.times)) )
    outf.write("\n")
    
    for n,m in nmRange:
        outf.write( "\n%d %d" % (n,m) )
        for k in kRange:
            entry = resultmatrix[ (n,m,k) ]
            outf.write(", (%s avg, %s sd, %s med)" % (entry.avg(), entry.sd(), entry.median()) )
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
                    kws = [ subprocess.Popen("TIMEFORMAT='%%3U' ; time %s -f %dstate%d.%d.%d -enum" % (EXEC, k, n, m, ext),
                            executable = "/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True) for ext in xrange(t,t+PARALLEL) ]
                    
                    # as they finish, assert the nonexistence of errors
                    #for ext in xrange(t,t+PARALLEL):
                    #    assert 0 == os.wait()[1]
                    #    
                    #    print ext,
                    #    sys.stdout.flush()
                    
                    # all are finished, thus extract and record runtimes
                    for kw in kws:
                        resultoutput, timeoutput = kw.communicate()
                        
                        try:
                            tcount = resultoutput.splitlines()[1].split()[0]
                            entry.add( tcount, timeoutput )
                        except Exception as e:
                            print "Error:", [c for c in resultoutput], [c for c in timeoutput]
                            print e
                
                subprocess.Popen("rm -r *state*", shell = True).wait()
                
                print "(%s trees)" % str(entry.avg())
    finally:
        subprocess.Popen("rm -r *state* mybigdata* *datalist seedms", shell = True).wait()
        
        outf = open("results_trees.txt", 'a')
        
        outf.write("Iterations = %s\n" % sys.argv[1])
        
        printResults(outf, resultmatrix, nmRange, kRange, TRIALS)
        
        outf.write("\n\n")
        
        outf.close()
        
        print "Results saved"

if __name__ == '__main__':
    runTests()
