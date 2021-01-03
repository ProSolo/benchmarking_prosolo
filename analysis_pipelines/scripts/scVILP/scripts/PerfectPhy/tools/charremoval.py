import sys
import subprocess
import itertools

EXEC = "./perfectphy"

def charremoval(inf):
    s = inf.read().split()
    
    INDEX, N, M = s[0], s[1], s[2]
    n = int(N)
    m = int(M)
    matrix = [ s[3 + i*m : 3 + i*m + m] for i in xrange(n) ]
    
    # loop over number of characters to remove
    for rem in xrange(0,m):
        count = m - rem
        
        # check every subset
        for p in itertools.combinations(range(m), count): #subsets(m, count):
            kw = subprocess.Popen([EXEC], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            pipe = kw.stdin
            w = pipe.write
            
            w("%s\n%s %d" % (INDEX, N, count) )
            
            for row in matrix:
                w('\n')
                for i in p:
                    w( row[i] )
                    w(' ')
            
            output, err = kw.communicate()
            
            if output[0] == "1": # return column indices with a pp
                return p
    
    return [] # return no columns

if __name__ == '__main__':
    try:
        result = charremoval(sys.stdin)
        
        if len(result) > 0:
            print "1 %d characters: %s" % (len(result), ' '.join(str(i) for i in result))
        else:
            print "0 No Perfect Phylogeny for any subset"
        
    except Exception as e:
        print type(e), e
