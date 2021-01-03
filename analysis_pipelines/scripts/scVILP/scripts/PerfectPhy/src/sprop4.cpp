/*
    Michael Coulombe 2012
    
    sprop4 changes:
        - Extended to C++ (fixed memory bugs)
        - Implemented full Kannan-Warnow improvements to DP algorithm to determine/construct perfect phylogenies
            > Extended the above to determine uniqueness (minimal)
        - Implemented KW algorithm for enumeration of minimal perfect phylogenies in a DAG
            > Created algorithms, a tree iterator, and simplified Dynamic Program method for the DAG
*/


/*
 * Written by Dan Gusfield copyright 2008
 * August 11, 2008 Permission is granted to use but not distribute this code or extensions or modifications of it.
 *
 *
 * Aug 2 increased the defined value of S, which is an upper bound on the number of proper clusters found.
 * July 30 create tracebackadj that produces an adjacency list of the PP created, where there is one.
 * July 28 create statetraceback to include the state labels for the nodes of the PP.
 * July 27 create taxatraceback to print out the proper cluster at each node in the PP.
 * July 26 begin programming the traceback to actually describe a PP if there is one. This  only shows
 *  the index of each proper cluster.
 * July 17, changed sdt[K] to sdt[M] in findpropclusters, which seems to have fixed a bug where some of the
 * data entries were changed (on homer but not imap) at each iteration in findpropclusters.
 * July 15, 2008 fixed a bug in findpropclusters - k had been used where statecount[j] was needed
 * DG June 7, not done - need to find the proper nesting of the tests of the connected components of the
 * complement of G' in G, and test each one to see if it is good

 * June 4, 2008 this modifies sprop1.c in the unfinished section for case 2 - building and searching the
relation when G - G' is not a proper cluster. In the unfinished section of sprop1.c I tried to be thoughtful
about making that code efficient in practice. But it became too messy. Here I do the most straightforward thing
of building the relation and recording it in a relation matrix, and then I take transitive closer of the
relation matrix. Later I can implement algorithm, not programming, ideas to make this part more efficient *

 * April 23, 2008 sprop1.c differs from spropercl.c in that it removes duplicate rows from the data array before additional
  processing

 * April 18, implemented the idea that we define every row to be a good proper cluster
 * April 19 Added in the full set with Sv vector all -1, as a final
step to finaly determine if there is a pp for the data;
 */

#include "sprop4.h"

// initialize global configuration variables (modified by command line args)
bool gbl::PRINTDATA       = false;
bool gbl::CHECKUNIQUE     = false;
bool gbl::ENUMERATE       = false;
bool gbl::DAGSTATISTICS   = false;
bool gbl::ALLOWINTERNAL   = false;
bool gbl::NEWICK          = false;

int gbl::OUTGROUPTAXON    = -1;

const char *gbl::DAGPRETTYF = NULL;
const char *gbl::DAGDOTOUTF = NULL;
//const char *gbl::DPDOTOUTF  = NULL;

int main(int argc, char * argv[]) {
    int dataindex, k, maxk = 0, n, m, properclustercount;
    const char* infilename = NULL;
    
    // define common datasets to be initialized throughout program
    vector<vector<int> > datam(0), properclusterm(0), Svm(0), SSvG(0);
    vector<int> statecountm(0), clustersizem(0);
    
    Dict pcdict;
    
    // read command line arguments
    ArgError err = handleArguments(argc, argv, infilename);
    if (err != NONE)
        return err;
    
    // read in dataset of taxon from input file (or stdin)
    if (infilename) {
        ifstream inf(infilename, ios_base::in);
        readmatrix(inf, datam, dataindex, n, m);
    }
    else
        readmatrix(cin, datam, dataindex, n, m);
    
    if (gbl::OUTGROUPTAXON == -1)
        gbl::OUTGROUPTAXON = 0;
    else if (gbl::OUTGROUPTAXON >= n) {
        printf("Bad argument: outgroup index (%d) >= input size (%d)\n", gbl::OUTGROUPTAXON, n);
        return -1;
    }
    
    if (DEBUGTRACE)
        printf("Read matrix...\n");
    
    if (PRINTREAD)
        prettyprintarray(datam, "datam", PPAF_DUMMY);
    
    // count the number of states in each column of data
    int skipChar = countstates (datam, statecountm, n, m);
    
    if (skipChar >= 0) {
        printf("Bad matrix: character %d is missing a character-state.\n", skipChar);
        return -1;
    }
    
    if (DEBUGTRACE) {
        printf("Counted states...\n");
        //prettyprintarray(statecountm, "statecount");
    }
    
    for (int j = 0; j < m; j++) {
        k = statecountm[j]; /* k is the  number of states seen in column j, which is one more
                               than the highest state label, since the states are numbered 
                               starting at 0. Jan. 16, 2009 */
        if (maxk < k) maxk = k; /* maxk is the largest number of states in any column, and
                                   hence is one larger than the highest state label. Jan. 16, 2009 */
    }
    
    // removes duplicate taxa in data and sorts them
    removeduprows (datam, statecountm, n, m);
    
    if (DEBUGTRACE || DEBUG2) {
        printf("Removed duplicate rows in data...\n");
        prettyprintarray(datam, "data");
        prettyprintarray(statecountm, "statecount");
    }
    
    // finds all proper clusters and their splitting vectors
    findpropclusters(datam, properclusterm, Svm, statecountm, n, m, maxk); 
    
    if (properclusterm.size() > (unsigned) numeric_limits<int>::max()) {
        printf("Number of proper clusters has overflown signed int (%zu > %d)\n", properclusterm.size(), numeric_limits<int>::max());
        
        return -1;
    }
    
    properclustercount = (int) properclusterm.size();
    
    if (DEBUGTRACE)
        printf("Found %d proper clusters...\n", properclustercount);

    if (CLUSTER) { 
        printf("\nThere are %d proper clusters. They are defined as follows:\n", properclustercount);
        prettyprintarray(properclusterm, "pc", PPAF_PC);
    }
    
    // sorts and removes duplicate proper clusters using radix sort
    // also constructs pcdict
    sort(properclusterm, Svm, clustersizem, pcdict);
    
    // cache new size
    properclustercount = (int) properclusterm.size();
    
    
    if (DEBUGTRACE) {
        printf("Sorted them into %d...\n", properclustercount);
        prettyprintarray(properclusterm, "pc", PPAF_PC);
        printf("And here is Sv...\n");
        prettyprintarray(Svm, "Sv", PPAF_DUMMY);
    }

    if (CLUSTER) {
        printf("\nThere are %d distinct proper clusters. \n", properclustercount);
        prettyprintarray(properclusterm, "pc", PPAF_PC);
        prettyprintarray(clustersizem, "Their sizes");
    }
    
    if (DEBUGTRACE)
        printf("n= %d m= %d k= %d |pc|= %d\n", n, m, maxk, properclustercount);
    
    if (DEBUGTRACE)
        printf("Constructed pcdict...Building S/Sv(G)...\n");
    
    findSSvG(SSvG, datam, statecountm, Svm, pcdict, properclusterm);
    
    if (DEBUGTRACE) {
        printf("Found S/Sv(G)...\n");
        prettyprintarray(SSvG, "S/Sv");
    }
    
    if (gbl::DAGPRETTYF && gbl::DAGSTATISTICS && !gbl::ENUMERATE) {
        // when specifying dagstats and an input file for it
        // load dag from file and go straight to statistics
        PPDAG dag(gbl::DAGPRETTYF);
        takeDagStatistics(dag, properclusterm, pcdict, datam, statecountm);
    }
    else {
        // run determining and construction dynamic program
        // as well as uniqueness, enumeration, and DAG algorithms if specified
        perf(datam, properclusterm, pcdict, Svm, SSvG, clustersizem, statecountm, properclustercount, n, m, maxk);
    }
}

/*
sfile | efile | stats | enum | result
  T   | T->F  |   T   | T->F | no enumeration, take statistics on dag in DAGPRETTYF
  T   |   T   |   T   |   F  | x
  T   |   T   |   F   |   T  | x
  T   |   T   |   F   |   F  | x
  
  T   |   F   |   T   | T->F | no enumeration, take statistics on dag in DAGPRETTYF
  T   |   F   |   T   |   F  | no enumeration, take statistics on dag in DAGPRETTYF
  T   |   F   |   F   |   T  | x
  T   |   F   |   F   |   F  | x

  F   |   T   |   T   |   T  | enumerate new dag, take statistics, and save to DAGPRETTYF
  F   |   T   |   T   |   F  | x
  F   |   T   |   F   |   T  | enumerate new dag, save to DAGPRETTYF
  F   |   T   |   F   |   F  | x

  F   |   F   |   T   |   T  | enumerate new dag, take statistics, and print to cout
  F   |   F   |   T   |   F  | nothing
  F   |   F   |   F   |   T  | enumerate new dag, print to cout
  F   |   F   |   F   |   F  | nothing
*/

ArgError handleArguments(int argc, char* argv[], const char* &infilename) {
    const char* helpmessage = 
"Usage: %s [options]\n"
"Options:\n"
"-h, --help       Displays this message then quits.\n"
"-f FILE          Read input matrix from FILE instead of stdin.\n"
"-newick          Prints the perfect phylogeny in the Newick tree format.\n"
"-outgroup T      Forces row index T as the outgroup taxon of the Newick tree.\n"
"\n"
"-unique          Determines if input has a unique minimal perfect phylogeny.\n"
"                 With -newick, output tree will be minimal.\n"
"-enum [FILE]     Runs enumeration algorithms to construct DAG of all minimal\n"
"                 perfect phylogenies. Optional: save DAG to FILE.\n"
"-dagstats [FILE] Runs statistical algorithms on DAG generated by enumeration.\n"
"                 If given FILE, will instead load the DAG from there.\n"
"-graphdag FILE   Outputs dot FILE representing the DAG.\n"
"\n"
"-showdata        Prints significant computed datasets and information.\n"
"-internal        Does not force taxa to be leaves during tree minimization.\n";
    
//-graphdp FILE    Outputs dot FILE of the initial phylogeny found
    
    string arg;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            printf("Bad argument: %s is not a flag\n", argv[i]);
            return ERROR;
        }
        
        arg = argv[i];
        
        if (arg.compare("-f") == 0) {
            if (argc > i + 1 && argv[i+1][0] != '-') {
                infilename = argv[i+1];
                i++;
            }
            else {
                printf("Bad argument: %s must precede filename\n", argv[i]);
                return ERROR;
            }
        }
        else if (arg.compare("-showdata") == 0)     gbl::PRINTDATA = true;
        else if (arg.compare("-unique")   == 0)   gbl::CHECKUNIQUE = true;
        else if (arg.compare("-internal") == 0) gbl::ALLOWINTERNAL = true;
        else if (arg.compare("-newick")   == 0)        gbl::NEWICK = true;
        else if (arg.compare("-outgroup") == 0) {
            if (argc > i + 1 && argv[i+1][0] != '-') {
                stringstream s(argv[i+1]);
                
                s >> gbl::OUTGROUPTAXON;
                
                if (s.fail()) {
                    printf("Bad argument: %s must precede a int >= 0\n", argv[i]);
                    return ERROR;
                }
                
                i++;
            }
            else {
                printf("Bad argument: %s must precede a int >= 0\n", argv[i]);
                return ERROR;
            }
        }
        else if (arg.compare("-dagstats") == 0) {
            if (argc > i + 1 && argv[i+1][0] != '-') {
                if (gbl::DAGPRETTYF && gbl::ENUMERATE)
                    gbl::ENUMERATE = false;
                
                gbl::DAGPRETTYF = argv[i+1];
                i++;
            }
            
            gbl::DAGSTATISTICS = true;
        }
        else if (arg.compare("-enum") == 0){
            if (argc > i + 1 && argv[i+1][0] != '-') {
                if (gbl::DAGPRETTYF && gbl::DAGSTATISTICS) {
                    i++;
                    continue;
                }
                
                gbl::DAGPRETTYF = argv[i+1];
                i++;
            }
            else if (gbl::DAGPRETTYF && gbl::DAGSTATISTICS) {
                continue;
            }
            
            gbl::ENUMERATE = true;
        }
        
        else if (arg.compare("-graphdag") == 0) {
            if (argc > i + 1 && argv[i+1][0] != '-')
                gbl::DAGDOTOUTF = argv[i+1];
            else {
                printf("Bad argument: %s must precede filename\n", argv[i]);
                return ERROR;
            }
            
            i++;
        }
        /*
        else if (arg.compare("-graphdp") == 0) {
            if (argc > i + 1 && argv[i+1][0] != '-')
                gbl::DPDOTOUTF = argv[i+1];
            else {
                printf("Bad argument: %s must precede filename\n", argv[i]);
                return false;
            }
            
            i++;
        }
        */
        else if (arg.compare("-h") == 0 || arg.compare("--help") == 0) {
            printf(helpmessage, argv[0]);
            return QUIT;
        }
        else {
            printf("Bad argument: %s is invalid flag\n", argv[i]);
            return ERROR;
        }
    }
    
    return NONE;
}


/**************************/
// if any column skips a number, returns index of that column
// else fills statecount with the number of states per column and returns -1
int countstates (const vector<vector<int> > &data, vector<int> &statecount, int n, int m) {
    statecount.resize(m);
    
    vector<int> observedStates; // used to check that no states are skipped
    
    for (int j = 0; j < m; j++) {
        observedStates.clear();
        
        int max = -1;
        for (int i = 0; i < n; i++) {
            int k = data[i][j];
            
            if (k > max) {
                max = k;
                observedStates.resize(max + 1, 0);
            }
            observedStates[k] = 1;
        }
        
        statecount[j] = max + 1;
        
        for(int k = 0; k <= max; k++) {
            if (observedStates[k] == 0) {
                return j;
            }
        }
    }
    
    return -1;
}

/**************************/

void readmatrix (istream &inf, vector<vector<int> > &data, int &dataindex, int &n, int &m) {
    if (!inf) {
        cerr << "Input read error" << endl;
        exit(1);
    }
    
    double cols; // multextract.pl may output a non-integral column number
    
    inf >> dataindex >> n >> cols;
    m = (int) cols; // m is truncated
    
    // read in data normally
    data.assign(n, vector<int>(m, 0));
    
    for (int i = 0; i < n; i++) {
        vector<int> &datai = data[i];
        
        for (int j = 0; j < m; j++) {
            int& cell = datai[j];
            inf >> cell;
        }
    }
}

/****************/
// removes duplicate taxa in data via radix sort
// O( m n k )
int removeduprows (vector<vector<int> > &data, const vector<int> &statecount, int &n, int m)
{
    vector<int> nextpermvector(n);
    vector<int> permvector(n);
    vector<vector<int> > savedata(n);
    
    // nextpermvector and permvector are used to implement radixsort without moving rows until all columns have been processed
    for (int p = 0; p < n; p++)
        permvector[p] = nextpermvector[p] = p;

    for (int j = m-1; j >= 0; j--) {               /* start of radix sort, processing the cols right to left */
        swap(permvector, nextpermvector);

        if (DEBUG)
            printf("The permutation vectors for radix sort are:\n");
        
        int p = 0;
        for (int k = 0; k < statecount[j]; k++) {
            for (int q = 0; q < n; q++) {  /* find the locations in col. j of the ks */
                if (data[ permvector[q] ][j] == k) {
                    nextpermvector[p] = permvector[q];
                    p++;
                }
            }
        }
        
        if (DEBUG)
            prettyprintarray(nextpermvector, "");
    }
    
    if (VERBOSE)
        printf("\nThe distinct (copies removed) data rows in radix sorted order are:\n");
    
    // swap table for permutation => data's rows become empty
    swap(data, savedata);
    
    int nn = -1, outgroupos = -1;
    
    // move the remaining data rows into their radix sorted order looking for and eliminating duplicate rows
    for (int q = 0; q < n; q++) {
        int i = nextpermvector[q];
        
        if (q == 0 || data[nn] != savedata[i] ) {
            nn++; // nn is an index, not a count, so it is one less than the number of distinct rows found so far
            swap(data[nn], savedata[i]);
        }
        
        if (i == gbl::OUTGROUPTAXON) // nn is the new position of the outgroup row at this point
            outgroupos = nn;
    }
    
    gbl::OUTGROUPTAXON = outgroupos;
    
    // this should change the value of n in the main to be the number of distinct rows, 
    // and hence the number of rows in the modified data array 
    n = nn+1;
    
    if (VERBOSEDATA) {
        prettyprintarray(data, "data");
        printf("There are %d distinct rows\n", n);
    }
    
    data.resize(n);
    
    return n;
}

/*
// traverses data to fill in missing data (-1) by columns in O(nm)
void fillmissingdata(vector<vector<int> > &data, vector<int> &statecount, int n, int m) {
    for (int j = 0; j < m; j++)
        for(int i = 0; i < n; i++)
            if (data[i][j] == -1)
                data[i][j] = statecount[j]++;
}
*/

/************************/

// increments a bitset in O(k)
// returns whether result is all zero or not
bool next_bitset(vector<bool>& num, int k) {
    bool carry = true;
    
    for(int i = 0; carry && i < k; i++) {
        vector<bool>::reference bit = num[i];
        
        carry = bit;
        bit.flip();
    }
    
    return !carry;
}

// increments a bitset in O(k)
// returns whether result is all zero or not
bool next_bitset(vector<int>& num, int k) {
    bool carry = true;
    
    for(int i = 0; carry && i < k; i++) {
        carry = num[i];
        num[i] ^= 1;
    }
    
    return !carry;
}

// find all proper clusters and splitting vectors in O(2^k m^2 n)

void findpropclusters (const vector<vector<int> > &data, vector<vector<int> > &propercluster, vector<vector<int> > &Sv, const vector<int> &statecount, int n, int m, int kStates)
{
    int state, sharedstate;
    
    bool provenproper;
    vector<int> includedrows(n), singleton(n,false), rowstates(kStates) /*includedrowstates(kStates), excludedrowstates(kStates)*/, splitter(m);
    
    vector<int> select(kStates, 0);
    
    // singleton will record which single rows are true proper clusters, because we want to *declare* every singleton row to be a proper cluster. That 
    // handles the case that an input row is actually an interior node in the prefect phylogeny
    
    propercluster.reserve(2*n);
    Sv.reserve(2*n);
    
    if (DEBUG)
        printf("INSIDE PROPCL\n");
    
    // insert pseudo proper cluster of all taxa
    propercluster.push_back(vector<int>(n,1));
    Sv.push_back(vector<int>(m,-1));
    
    // O(m)
    for (int col = 0; col < m; col++) {
        const int k = statecount[col];
        
        // enumerate the appropriate selector vectors in array sets, the first and last rows specify the trivial sets, so ignore them
        // O(2^k) iterations
        while( next_bitset(select, k - 1) ) {
            // increments in O(k) from [1,0,..,0] to [1,..,1,0], terminates at [0,..,0,0]
            
            if (VERBOSE1) {
                printf("\n Starting next selector vector for column %d and k = %d \n", col, k);
                printf("The states included in this selector = { ");
                for(int i = 0; i < k; i++)
                    printf(select[i] ? "1 " : "0 ");
                printf("}\n");
            }
            
            // we start with the optimistic expectation that the cluster is proper, and look for a disproof of that expectation
            provenproper = true; 
            
            if (DEBUG2) {
                printf ("In findpropclusters, the data array is\n"); 
                prettyprintarray(data, "data");
            } 
            
            //for a given selector vector, look down column col to see which rows should be included in the potential proper cluster
            
            int includedcount = 0, lastincludedrow = -1, lastexcludedrow = -1;
            
            for (int i = 0; i < n; i++) { 
                state = data[i][col];
                
                if (select[state]) {
                    if (VERBOSE1) {
                        printf("row %d is included in the potential cluster\n", i);
                    }
                    includedcount++;
                    includedrows[i] = 1;
                    lastincludedrow = i;
                }
                else {
                    includedrows[i] = 0;
                    lastexcludedrow = i;
                }
            }
            
            for (int j = 0; j < m; j++) {   
                // for a fixed col and a fixed state partition, we may find cols. j which share a state across
                // the induced partition of data rows, and splitter will accumulate that information
                
                splitter[j] = -1;
                
                if (j == col)
                    continue;
                
                sharedstate = statecount[j];  // changed k to statecount[j] here and places below. July 15, 2008
                
                enum { NOTSEEN = -1, EXCLUDED = 0, INCLUDED = 1, SHARED = 2 };
                
                rowstates.assign(sharedstate, NOTSEEN); // changed K to sharedstate July 17, 2008
                
                for (int i = 0; i < n; i++) {
                    state = data[i][j];
                    
                    if (DEBUG)
                        printf("Working down column %d, at row %d, seeing state %d\n", j, i, state);
                    
                    int side = includedrows[i];
                    
                    if (rowstates[state] == NOTSEEN)
                        rowstates[state] = side;
                    
                    else if (rowstates[state] == 1 - side) {
                        if (DEBUG)
                            printf("Here we catch a shared state with state %d and rowstates[state] = %d\n", state, rowstates[state]);
                        if (VERBOSE)
                            printf("The intersection contains the character-state pair (%d, %d)\n", j, state);
                        
                        if (sharedstate == statecount[j]) {
                            // state is the first shared state seen
                            sharedstate = state;
                            splitter[j] = state;
                            rowstates[state] = SHARED;
                        }
                        else {
                            if (DEBUG)
                                printf("Case %d: Selector vector for column %d fails because of column %d\n", 2 - includedrows[i], col, j);
                            
                            // state is the second shared state seen
                            provenproper = false;
                            break;
                        }
                    }
                }
                
                if (!provenproper)
                    break;
                
                if (DEBUG)
                    printf("W\n");
            } // end of the for j loop. The break above breaks out of the for j loop

            if (provenproper) {
                if (VERBOSE)
                    printf("Partition of states in column %d, induces a proper cluster\n", col);
                
                // insert proper cluster and its complement
                
                if (includedcount == 1) {
                    singleton[lastincludedrow] = true;
                }
                else if (includedcount == n-1) {
                    singleton[lastexcludedrow] = true;
                }
                propercluster.push_back(includedrows);
                
                for (int i = 0; i < n; i++) {
                    includedrows[i] ^= 1;
                }
                propercluster.push_back(includedrows);
                
                // insert splitting vector (equal for pc and comp)
                Sv.push_back(splitter); /*Dec. 30, 2008 add these two lines, and delete a line below */
                Sv.push_back(splitter);
                
                if (VERBOSE)
                    prettyprintarray(Sv.back(), "The Sv vector is: ", PPAF_DUMMY);
                
                
            }
            else if (VERBOSE)
                printf("Partition of states in column %d, does NOT induce a proper cluster\n", col);
            
        } /* end of the for rowcount loop */
    } /* end of the for col loop */
    

    for (int i = 0; i < n; i++) {
        // here we look through the rows to spot any single row that is not a true proper cluster, in order to *declare* it to be a proper cluster,
        // in order to handle the case of rows that should be interior nodes - does this correctly handle interiors?
        
        if (! singleton[i]) {
            if (DEBUGTRACE)
                printf("Forcing %d to be proper cluster...", i);
            
            propercluster.push_back(vector<int>(n,0));
            propercluster.back()[i] = 1;
            
            Sv.push_back(data[i]);
        }
    }
}


/***********************/
/*
long expt ( int a, int b)
{
    int i;
    long result = 1;
    
    for (i = 0; i < b; i++)
        result = result * a;
    
    return result;
}
long exp2(int b)
{
    long result = 1;
    return result << b;
}
*/
/********************/

unsigned int digits(unsigned int n) {
    unsigned int i, result = 1;
    
    for (i = 10; i <= n; i *= 10)
        result++;
    
    return result;
}

void prettyprintarray(const vector<int> &array, const char name[], PrettyPrintFormat format) {
    unsigned int i, size = array.size();
    
    switch (format) {
        case PPAF_DUMMY: {
            printf("%s = { ", name);
            for (i = 0; i < size; i++) {
                if (array[i] < 0)
                    printf("* ");
                else
                    printf("%d ", array[i]);
            }
            printf("}\n");
            break;
        }
        case PPAF_PC: {
            unsigned int elemWidth = digits(size);
            
            printf("%s = { ", name);
            for (i = 0; i < size; i++) {
                if (array[i] == 0)
                    printf("- ");
                else
                    printf("%*u ", elemWidth, i);
            }
            printf("}\n");
            break;
        }
        case PPAF_NONE: {
            printf("%s = { ", name);
            for (i = 0; i < size; i++)
                printf("%d ", array[i]);
            printf("}\n");
            break;
        }
    }
}

void prettyprintarray(const vector<vector<int> > &array, const char name[], PrettyPrintFormat format)
{
    unsigned int i, size = array.size(), indexWidth = digits(size);
    
    for (i = 0; i < size; i++) {
        printf("%s[%*u", name, indexWidth, i);
        prettyprintarray(array[i], "]", format);
    }
}

/*****************/


// radix sort of parallel vectors propercluster and Sv
// O( 2^k m n )
void sort(vector<vector<int> > &propercluster, vector<vector<int> > &Sv, vector<int> &clustersize, Dict &qtable)
{
    unsigned int properclustercount = propercluster.size();
    const unsigned int n = propercluster[0].size();
    
    vector<vector<int> > savepropercluster(properclustercount), saveSv(properclustercount); // these are swap tables
    vector<int> nextpermvectors(properclustercount), permvectors(properclustercount),
                saveclustersize(properclustercount),
                nextofsize(n+1, 0),
                countclustersize(n+1, 0);
    
    // nextpermvector and permvectors are used to implement radixsort without moving rows until all columns have been processed
    for (unsigned int p = 0; p < properclustercount; p++) {
        permvectors[p] = nextpermvectors[p] = p;
    }
    
    // O( 2^k m n )
    for (int j = n-1; j >= 0; j--) {               // start of radix sort, processing the cols right to left
        swap(permvectors, nextpermvectors);
        
        if (DEBUG)
            printf("The permutation vectors for radix sort are:\n");
        
        unsigned int p = 0;
        for (unsigned int q = 0; q < properclustercount; q++) {  // find the locations in col. j of the 0s
            int i = permvectors[q];
            
            if (propercluster[i][j] == 0) {
                nextpermvectors[p] = permvectors[q];
                p++;
            }
        }
        
        for (unsigned int q = 0; q < properclustercount; q++) { // find the locations in col. j of the 1s
            int i = permvectors[q];
            
            if (propercluster[i][j] == 1) {
                nextpermvectors[p] = permvectors[q];
                p++;
            }
        }

        if (DEBUG)
            prettyprintarray(nextpermvectors, "nextpermvectors");
    }

    if (VERBOSE)
        printf("\n\nThe distinct (copies removed) properclusters in radix sorted order are:\n");
    
    // move tables for permutating
    swap(propercluster, savepropercluster);
    swap(Sv, saveSv);
    
    int end = -1;
    // O(2^k m n)
    for (unsigned int q = 0; q < properclustercount; q++) { // move the propercluster rows into their radix sorted order
                                               // looking for and eliminating duplicate rows
        int i = nextpermvectors[q];
        
        if (q == 0 || propercluster[end] != savepropercluster[i]) {
            end++; // end is an index, not a count, so it is one less than the number of distinct proper clusters found so far
            
            swap(propercluster[end], savepropercluster[i]);
            
            if (VERBOSE) {
                printf("propercluster[%d", end);
                prettyprintarray(propercluster[end], "]", PPAF_PC);
            }
            
            unsigned int count = accumulate(propercluster[end].begin(), propercluster[end].end(), 0);
            
            saveclustersize[end] = count;
            countclustersize[count]++;

            swap(Sv[end], saveSv[i]); // parallel the move of propercluster rows for Sv rows
        }
    }
    properclustercount = end+1;
    
    if (VERBOSE) {
        printf("There are %u distinct properclusters, which have sizes", properclustercount);
        prettyprintarray(saveclustersize, "");
    }
    
    unsigned int totalcount = 0;
    for (unsigned int count = 1; count <= n; count++) {
        nextofsize[count] = totalcount;
        totalcount += countclustersize[count];
        
        if (VERBOSE) {
            printf ("The number of proper clusters of size %u is %d\n", count, countclustersize[count]);
            printf ("The starting location for proper clusters of size %u is %d\n", count, nextofsize[count]);
        }

    }
    
    // resize tables to true count
    propercluster.resize(properclustercount);
    Sv.resize(properclustercount);
    
    // clear swap tables
    savepropercluster.assign(properclustercount, vector<int>(0));
    saveSv.assign(properclustercount, vector<int>(0));
    
    //printf("Radix sorted proper cluster table:\n");
    //prettyprintarray(propercluster, "pc", PPAF_PC);
    
    #if UseQTable
        qtable.construct(propercluster);
    #endif
    
    swap(propercluster, savepropercluster);
    swap(Sv, saveSv);
    clustersize.assign(properclustercount, -1);
    
    // O( 2^k m n )
    for (unsigned int i = 0; i < properclustercount; i++) { // Now sort the distinct properclusters by the number of 1s (their sizes)
        unsigned int count = saveclustersize[i];
        unsigned int q = nextofsize[count];
        nextofsize[count]++;
        
        clustersize[q] = saveclustersize[i];
        swap(propercluster[q], savepropercluster[i]);
        swap(Sv[q], saveSv[i]);   // parallel the move of propercluster rows for Sv rows
        
        #if UseQTable
            qtable.setIndexMap(i, q); // since propercluster is resorted, qtable must know the true indices
        #endif
    }
    
    #if 1 - UseQTable
        qtable.construct(propercluster);
    #endif
    
    if (VERBOSE)  {
        printf("\nThe distinct properclusters sorted by size are:\n");
        for (unsigned int i = 0; i < properclustercount; i++) {
            printf("propercluster[%u",i);
            prettyprintarray(propercluster[i], "]", PPAF_PC);
            
            printf("\nWhich has size %d\n", clustersize[i]);
            prettyprintarray(Sv[i], "The associated Sv vector", PPAF_DUMMY);
            printf("\n");
        }
    }
}

/*********************/

template <class DataMatrix>
void createDot_findEqClassesMulti(int maxC, const vector<int> &eqclasses, const DataMatrix &data, const vector<vector<int> > &prevTaxa) {
    const unsigned int n = data.size(), m = data[0].size();
    
    cout << "digraph findEqClassesMulti {" << endl;
    
    for (unsigned int curTaxa = 0; curTaxa < n; curTaxa++) {
        cout << curTaxa << ';' << endl;
    }
    
    for (unsigned int curChar = 0; curChar < m; curChar++) {
        for (unsigned int curTaxa = 0; curTaxa < n; curTaxa++) {
            const int prev = prevTaxa[curChar][curTaxa];
            
            if (prevTaxa[curChar][curTaxa] != -1) {
                cout << prev << " -> " << curTaxa << ';' << endl;
            }
        }
    }
    
    cout << '}' << endl;
}



// creates linear two-way adjacency lists for each character, grouping them by common states, O(nm)
template <class DataMatrix>
void createStateGraph(const DataMatrix &data, const vector<int> &statecount, vector<vector<int> > &prevTaxa, vector<vector<int> > &nextTaxa) {
    const unsigned int n = data.size(), m = data[0].size();
    
    prevTaxa.assign(m, vector<int>(n,-1));
    nextTaxa.assign(m, vector<int>(n,-1));
    vector<int> prevStateTaxa;
    
    // prevTaxa is the "outward" edges for each taxa in each state class
    // nextTaxa is the "inward" edges for each taxa in each state class
    // prevStateTaxa points to the end taxon of each eq class
    
    for (unsigned int curChar = 0; curChar < m; curChar++) {
        prevStateTaxa.assign(statecount[curChar], -1);
        
        for (unsigned int curTaxa = 0; curTaxa < n; curTaxa++) {
            int curState = data[curTaxa][curChar];
            int &prev = prevStateTaxa[curState];
            
            // -1 means first taxa in this class
            if (prev != -1) {
                // set it as the taxa at the end of this state class
                nextTaxa[curChar][prev] = curTaxa;
                prevTaxa[curChar][curTaxa] = prev;
            }
            prev = curTaxa;
        }
    }
}



// finds the equivalence classes of G/splitter in O(nm)
template <class SplitVector, class DataMatrix>
int findEqClassesMulti(const SplitVector &splitter, vector<int> &eqclasses, const DataMatrix &data, const vector<vector<int> > &prevTaxa, const vector<vector<int> > &nextTaxa) {
    const unsigned int n = data.size(), m = data[0].size();
    
    vector<int> seen(n,false), queue;
    int maxC = 0;
    
    eqclasses.assign(n,-1);
    queue.reserve(n);
    
    for (unsigned int i = 0; i < m; i++) {
        if (splitter[i] == -1) {
            continue; // when Sv is in a dummy state, that character doesn't influence the classes
        }
        
        // check every unseen taxon, O(n) total per proper cluster
        for (unsigned int t = 0; t < n; t++) {
            if (seen[t]) {
                continue; // taxon already exists in an eq class
            }
            
            if (splitter[i] == data[t][i]) {
                continue; // states in Sv don't influence the classes
            }
            
            // new eq class possibly found, find connected components in O(nm) total per proper cluster
            
            queue.assign(1, t);
            seen[t] = true;
            int curClass = maxC;
            
            // check every edge, O(m) once per taxa total
            for (unsigned int front = 0; front < queue.size(); front++) {
                int curTaxa = queue[front];
                
                for (unsigned int curChar = 0; curChar < m; curChar++) {
                    // skip if edge is restricted by Sv
                    if (splitter[curChar] == -1 || splitter[curChar] == data[curTaxa][curChar]) {
                        continue;
                    }
                    
                    int pt = prevTaxa[curChar][curTaxa];
                    
                    if (pt != -1) {
                        // found a real edge, add the taxa to eq class
                        if (!seen[pt]) {
                            // first time seeing taxa
                            seen[pt] = true;
                            queue.push_back(pt);
                        }
                        else if (eqclasses[pt] != -1) {
                            // found another class, change from making a new one to adding to it
                            curClass = eqclasses[pt];
                        }
                    }
                    
                    pt = nextTaxa[curChar][curTaxa];
                    
                    if (pt != -1) {
                        // found a real edge, add the taxa to eq class
                        if (!seen[pt]) {
                            // first time seeing taxa
                            seen[pt] = true;
                            queue.push_back(pt);
                        }
                        else if (eqclasses[pt] != -1) {
                            // found another class, change from making a new one to adding to it
                            curClass = eqclasses[pt];
                        }
                    }
                }
            } // end check edges
            
            // label the taxa in the class
            unsigned int queuesize = queue.size();
            for (unsigned int j = 0; j < queuesize; j++) {
                eqclasses[ queue[j] ] = curClass;
            }
            
            if (curClass == maxC) {
                maxC++;
            }
        }
    }
    
    // any straglers are put in their own eq classes
    for (unsigned int i = 0; i < n; i++) {
        if (!seen[i]) {
            eqclasses[i] = maxC;
            maxC++;
        }
    }
    
    if (DEBUGTRACE) {
        //createDot_findEqClassesMulti(maxC, eqclasses, data, prevTaxa);
    }
    
    return maxC;
}



// finds equivalence classes of every S/Sv(G) in O(Snm)
void findSSvG(vector<vector<int> > &SSvG, const vector<vector<int> > &data, const vector<int>& statecount, const vector<vector<int> > &Sv, const Dict &pcdict, const vector<vector<int> > &propercluster) {
    const unsigned int properclustercount = propercluster.size();
    
    vector<vector<int> > prevTaxa, nextTaxa;
    
    createStateGraph(data, statecount, prevTaxa, nextTaxa);
    
    // scan prevTaxa from them all to calculate connected components in O(nm) per proper cluster
    SSvG.resize(properclustercount);
    
    for (unsigned int x = 0; x < properclustercount; x++) {
        // if complement already computed this, then skip
        if (! SSvG[x].empty())
            continue;
        
        // otherwise compute S/Sv(G)
        findEqClassesMulti(Sv[x], SSvG[x], data, prevTaxa, nextTaxa);
        
        // then copy to S/Sv(S-G)
        int comp = pcdict.lookupComplement(propercluster[x]);
        
        if (comp != -1)
            SSvG[comp] = SSvG[x];
    }
}

class TransposedMatrix {
protected:
    typedef const vector<int>* R;
    vector<R> matrix;
    
    class TransposedRow {
    protected:
        const TransposedMatrix *rows;
        int col;
    public:
        TransposedRow(const TransposedMatrix *m, int j):rows(m), col(j) { }
        
        size_t size() const {
            return rows->matrix.size();
        }
        int operator[](size_t i) const {
            return (* (rows->matrix[i]) )[col];
        }
    };
    
public:
    void reserve(size_t cap) {
        matrix.reserve(cap);
    }
    void push_back(R e) {
        matrix.push_back(e);
    }
    size_t size() const {
        // all are the same length
        return matrix[0]->size();
    }
    TransposedRow operator[](size_t j) const {
        return TransposedRow(this,j);
    }
};

template<int value>
class InfiniteOneValueVector {
public:
    InfiniteOneValueVector() {}
    inline int operator[](size_t) const { return value; }
};

// merges two equivalence classes into eqclass in O(n) - calls findEqClassesMulti on two vectors
// returns the largest number of a class
int mergeEqClasses(vector<int>& eqclass, const vector<int>& v0, const vector<int>& v1, int n) {
    TransposedMatrix data; // the data must be a map from taxa to "states," thus is the transpose of { v0 , v1 }
    data.reserve(2);
    data.push_back(&v0);
    data.push_back(&v1);
    
    const InfiniteOneValueVector<-137> dummysplitter; // there is no true splitter to consider, so this should never conflict with any value in [-1,n)
    
    vector<vector<int> > prevTaxa, nextTaxa; // will be 2 by n
    vector<int> classcount(2);
    
    classcount[0] = 1 + *max_element(v0.begin(), v0.end());
    classcount[1] = 1 + *max_element(v1.begin(), v1.end());
    
    createStateGraph(data, classcount, prevTaxa, nextTaxa);
    
    return findEqClassesMulti(dummysplitter, eqclass, data, prevTaxa, nextTaxa);
}

// find the proper clusters (as equivalence classes) of G2/Sv(G,G1) in O(n)
// returns whether successful or not
bool findDecompSvGG1(const Dict &pcdict, const vector<vector<int> > &SSvG, const vector<int> &good, const vector<vector<int> > &propercluster, vector<int> &tempdecomp, int G, int G1, int n) {
    vector<int> eqclass;
    
    int maxC = mergeEqClasses(eqclass, SSvG[G], SSvG[G1], n); // O(n)
    
    return pcdict.lookupParallel(eqclass, good, tempdecomp, G, G1, n, maxC, propercluster); // O(n)
}

bool decompCoversTaxa(const vector<int>& decomp, const vector<int>& clustersize, int H) {
    int count = 0;
    for (unsigned int i = 0; i < decomp.size(); i++) {
        count += clustersize[ decomp[i] ];
    }
    
    return (count == clustersize[H]);
}



// find the proper clusters (as equivalence classes) of H/cl(decomp) in O(d n) where |decomp| = d <= n (O(kn) as d <= k if using k-MIS)
// returns whether successful or not
bool findExpandedDecomp(const Dict &pcdict, const vector<vector<int> > &SSvG, vector<int> &decomp, int H, int n, const vector<int> &clustersize, const vector<vector<int> > &verygood) {
    // treats each S/Sv(G) for G in decomp as a column in a data matrix
    // then merges the equivalence classes column-wise
    // and looks up the proper clusters in the result
    
    TransposedMatrix data;
    vector<int> classcount, finaleqclasses;
    vector<vector<int> > prevTaxa, nextTaxa;
    const InfiniteOneValueVector<-137> dummysplitter;
    
    const int datasize = decomp.size() + 2;
    data.reserve(datasize);
    classcount.reserve(datasize);
    
    for (unsigned int i = 0; i < decomp.size(); i++) {
        const vector<int>& eqclasses = SSvG[ decomp[i] ];
        
        classcount.push_back(1 + *max_element(eqclasses.begin(), eqclasses.end()));
        data.push_back(&eqclasses);
        
        if (DEBUGTRACE) {
            prettyprintarray(eqclasses, "partial eqclasses");
            /*
            vector<int> SvDecomp;
            
            if (!pcdict.lookupParallel(eqclasses, SvDecomp, n, classcount[i], clustersize)) {
                printf("ERROR: Not all pc's ???");
            }
            
            sort(SvDecomp.begin(), SvDecomp.end());
            printf("Sv(%d)", decomp[i]);
            prettyprintarray(SvDecomp, ": PC's represented");
            */
        }
    }
    
    // add an equivalence class of containment within members of the decomposition
    // in order to force the old decomp to be a subset of the new one
    vector<int> decompClasses(n, -1);
    {
        int numClasses = (int) decomp.size();
        
        for (int t = 0; t < n; t++) {
            for (unsigned int i = 0; i < decomp.size(); i++) {
                if (decomp[i] >= t && verygood[ decomp[i] ][t] != 0) {
                    decompClasses[n-t-1] = i; // assign eqclass index to decomp index
                    break;
                }
            }
            if (decompClasses[n-t-1] == -1) {
                decompClasses[n-t-1] = numClasses; // assign eqclass index to next largest unused index
                numClasses++;
            }
        }
        
        if (DEBUGTRACE) {
            prettyprintarray(decompClasses, "decomp eqclasses");
        }
        
        classcount.push_back(numClasses);
        data.push_back(&decompClasses);
    }
    
    if (true) {
        // do we need to add S/Sv(H) ?
        const vector<int>& eqclasses = SSvG[H];
        
        classcount.push_back(1 + *max_element(eqclasses.begin(), eqclasses.end()));
        data.push_back(&eqclasses);
        
        if (DEBUGTRACE) {
            prettyprintarray(eqclasses, "supset eqclasses");
            /*
            vector<int> SvDecomp;
            
            if (!pcdict.lookupParallel(eqclasses, SvDecomp, n, 1 + *max_element(eqclasses.begin(), eqclasses.end()), clustersize)) {
                printf("ERROR: Not all pc's ???");
            }
            
            //sort(SvDecomp.begin(), SvDecomp.end());
            printf("Sv(%d)", H);
            prettyprintarray(SvDecomp, ": PC's represented");
            */
        }
    }
    
    createStateGraph(data, classcount, prevTaxa, nextTaxa);
    
    int maxC = findEqClassesMulti(dummysplitter, finaleqclasses, data, prevTaxa, nextTaxa);
    
    if (DEBUGTRACE) {
        prettyprintarray(finaleqclasses, "unioned eqclasses");
    }
    
    // now eqclasses holds the union of all S/Sv(G)
    // thus verify the proper clusters they correspond to
    
    #if UseQTable
        bool success = pcdict.lookupParallel(finaleqclasses, decomp, n, maxC, clustersize); // O(n)
    #else
        bool success = pcdict.lookupParallel(finaleqclasses, decomp, n, maxC); // O(n)
    #endif
    
    // filter decomp for only verygood partners of H
    
    if (DEBUGTRACE) {
        prettyprintarray(decomp, "d before filter");
    }
    
    for(unsigned int i = 0; i < decomp.size(); i++) {
        if (H < decomp[i] || verygood[H][ decomp[i] ] != 2) {
            
            
            if (H > decomp[i] && verygood[H][ decomp[i] ] == 1) {
                // a not verygood subset cannot be in the decomp
                success = false;
                
                //if (DEBUGTRACE) {
                    printf("    %d : subset but not verygood! OH NO!\n", decomp[i]);
                //}
            }
            else if (DEBUGTRACE) {
                printf("    %d : not a subset\n", decomp[i]);
            }
            // remove from decomp by moving back() to position i
            decomp[i] = decomp.back();
            decomp.pop_back();
            i--;
            
        } else if (DEBUGTRACE) {
            printf("    %d : a-ok!\n", decomp[i]);
        }
    }
    
    if (DEBUGTRACE) {
        prettyprintarray(decomp, "d after filter");
    }
    
    // should always cover taxa if every subset is good
    if (! decompCoversTaxa(decomp, clustersize, H)) {
        prettyprintarray(decomp, "Expansion doesn't cover taxa");
        success = false;
    }
    
    return success;
}

/*
bool findExpandedDecomp(const Dict &pcdict, const vector<vector<int> > &SSvG, vector<int> &decomp, int H, int n, const vector<int> &clustersize, const vector<vector<int> > &verygood) {
    vector<int> eqclass = SSvG[ decomp[0] ], // initialize to eq classes of the first pc
                nexteqclass;                 // temporary contents
    
    // find initial maxC O(d)
    int maxC = 1 + *max_element(eqclass.begin(), eqclass.end());
    
    for (unsigned int i = 1; i < decomp.size(); i++) {
        maxC = mergeEqClasses(nexteqclass, eqclass, SSvG[ decomp[i] ], n); // O(n)
        
        eqclass.swap(nexteqclass);
    }
    
    // now eqclasses holds the union of all S/Sv(G)
    // thus verify the proper clusters they correspond to
    
    #if UseQTable
        bool success = pcdict.lookupParallel(eqclass, decomp, n, maxC, clustersize); // O(n)
    #else
        bool success = pcdict.lookupParallel(eqclass, decomp, n, maxC); // O(n)
    #endif
    
    for(unsigned int i = 0; i < decomp.size(); i++) {
        if (H < decomp[i] || verygood[H][ decomp[i] ] != 2) {
            decomp[i] = decomp.back();
            decomp.pop_back();
            
            if (verygood[H][ decomp[i] ] == 1) {
                // a not verygood subset cannot be in the decomp
                success = false;
            }
        }
    }
    
    return success;
}
*/

///////////////////////////////

// root = root (+) sub
void unionStates(vector<int> &root, const vector<int> &sub) {
    const unsigned int size = root.size();
    
    for (unsigned int i = 0; i < size; i++)
        if (root[i] == -1)
            root[i] = sub[i];
}


bool isDecompUnique(const vector<int> &decomp, const vector<int> &unique) {
    const unsigned int size = decomp.size();
    
    for (unsigned int i = 0; i < size; i++)
        if (unique[ decomp[i] ] == 0)
            return false;
    
    return true;
}

class SetEqChecker {
private:
    vector<unsigned long long> buckets;
    unsigned long long curHash;
    
public:
    // must be called if curHash would overflow (not likely if 64-bit...)
    void reset(int len) {
        buckets.assign(len, 0);
        curHash = 1;
    }
    
    /*
    Given that every element in d1,d2 are less than buckets.size(),
    and that d1 and d2 individually have no duplicates,
    returns whether or not they contain the same elements in O(n) time
    */
    bool compare(const vector<int>& d1, const vector<int> &d2) {
        const unsigned int n = d1.size();
        
        if (n != d2.size())
            return false;
        
        unsigned long long uniqueHash = curHash++;
        
        // add all pcs in d1 to the table
        for(unsigned int i = 0; i < n; i++)
            buckets[ d1[i] ] = uniqueHash;
        
        // then check if all pcs in d2 came from d1
        for(unsigned int i = 0; i < n; i++)
            if (buckets[ d2[i] ] != uniqueHash)
                return false;
        
        return true;
    }
};

// "bucket-sort" comparison of two decompositions in O(n) per call, with a one-time price of O(L) to make table
// works because every call to minimizeSubtrees is for a unique pair (H,G1), and table memory is persistent
// TODO: this only works if each H is only ever used once as an argument to this function
// must use pairs of H,G1 to get a reusable comparator for H
/*
bool decompEqual(const vector<int> &d1, const vector<int> &d2, int H, vector<int>& table) {
    const unsigned int n = d1.size();
    
    if (n != d2.size())
        return false;
    
    // add all pcs in d1 to the table
    for(unsigned int i = 0; i < n; i++)
        table[ d1[i] ] = H;
    
    // then check if all pcs in d2 came from d1
    for(unsigned int i = 0; i < n; i++)
        if (table[ d2[i] ] != H)
            return false;
    
    return true;
}
*/

// check whether a particular decomposition is minimal, if not minimize it, and check if it is unique for H
// time O(n + m)
int minimizeSubtrees(vector<vector<int> > &rootlabels, vector<vector<int> > &gooddecomp, vector<int> &newDecomp, const vector<vector<int> > &Sv, const vector<int> &clustersize, vector<int> &unique, /*vector<int> */SetEqChecker &table, bool first, int H, int G1) {
    int G2 = newDecomp[1];
    
    // enum variable to store how to minimize the decomp
    enum { NONE, G1_ONLY, G2_ONLY, EITHER, BOTH } numMinimized = NONE;
    
    vector<int> clH(Sv[H].size());
    
    // find the unminimized root label of H O(m)
    // if dummy states, thus find canonical labeling from the 3 edges
    // if fully set, thus label = Sv(H,G1)
    if (newDecomp.size() == 2)
        getCanonicalLabeling(newDecomp, Sv, clH, H);
    else {
        // G1 is in sorted order
        compatible(Sv[H], Sv[G1], clH);
    }
    
    // if any subtrees are not unique, then H is not unique O(n)
    bool isUnique;
    if (!isDecompUnique(newDecomp, unique))
        isUnique = false;
    else
        isUnique = unique[H];
    
    if (isUnique) {
        // if all subtrees are unique, then test newDecomp for minimal uniqueness
        
        if (newDecomp.size() == 2) {
            // decomposition is two proper clusters
            
            // At this point, both subtrees have unique decompositions.
            // The only way for minimization to make two different trees is for
            // them to be exclusively minimizable.
            // If taxa cannot be internal, discounts them
            
            if (compatible(clH, rootlabels[G1]) && (gbl::ALLOWINTERNAL || clustersize[G1] > 1) )
                numMinimized = G1_ONLY;
            
            if (compatible(clH, rootlabels[G2]) && (gbl::ALLOWINTERNAL || clustersize[G2] > 1)) {
                if (numMinimized == G1_ONLY) {
                    numMinimized = BOTH;
                    
                    if (!compatible( rootlabels[G1], rootlabels[G2])) {
                        // if both subtrees can be contracted but they aren't compatible with each other,
                        // then there are two possible decompositions for each case, thus not unique
                        isUnique = false;
                        numMinimized = EITHER;
                    }
                }
                else
                    numMinimized = G2_ONLY;
            }
            
            // minimize decomp in the ways specified and set the root of H
            
            switch (numMinimized) {
                case BOTH:
                case G2_ONLY:
                    // minimize G2
                    if (clustersize[G2] > 1) {
                        newDecomp.erase( --newDecomp.end() );
                        newDecomp.insert(newDecomp.end(), gooddecomp[G2].begin(), gooddecomp[G2].end());
                    }
                    unionStates(clH, rootlabels[G1]);
                    
                    if (numMinimized == G2_ONLY)
                        break;
                    // else continue to next case
                    
                case EITHER:
                case G1_ONLY:
                    // minimize G1
                    if (clustersize[G1] > 1) {
                        newDecomp.erase(newDecomp.begin());
                        newDecomp.insert(newDecomp.end(), gooddecomp[G1].begin(), gooddecomp[G1].end());
                    }
                    unionStates(clH, rootlabels[G1]);
                    
                    //isUnique = (numMinimized != EITHER);
                
                case NONE:
                    break;
            }
        }
        else {
            // decomposition is a proper cluster G1 and a set of other pc's which are necessarilly minimal
            // except the case where there are is an interior taxa identical to clH which is fully specified
            // so all that must be done is minimizing G1, uniqueness is preserved
            
            if (compatible(clH, rootlabels[G1]) && clustersize[G1] > 1) {
                newDecomp.erase(newDecomp.begin());
                newDecomp.insert(newDecomp.begin(), gooddecomp[G1].begin(), gooddecomp[G1].end());
            }
        }
        
        if (first || !isUnique) {
            // if this is the first decomp found or the first non-unique one found, then output it
            gooddecomp[H] = newDecomp;
            rootlabels[H] = clH;
        }
        else {
            // if this decomp is unique and not the first, compare it with the first one found
            //isUnique = decompEqual(gooddecomp[H], newDecomp, H, table);
            isUnique = table.compare(gooddecomp[H], newDecomp);
        }
    }
    else {
        // H is not unique already
        
        if (first) {
            // if this is the first decomp found and H is not unique, thus just minimize this decomp and output it
            
            if (clustersize[G1] > 1 && compatible(clH, rootlabels[G1])) {
                newDecomp.erase(newDecomp.begin());
                newDecomp.insert(newDecomp.begin(), gooddecomp[G1].begin(), gooddecomp[G1].end());
                
                unionStates(clH, rootlabels[G1]);
            }
            
            if (clustersize[G2] > 1 && compatible(clH, rootlabels[G2])) {
                newDecomp.erase(newDecomp.begin());
                newDecomp.insert(newDecomp.begin(), gooddecomp[G1].begin(), gooddecomp[G1].end());
            }
            
            gooddecomp[H] = newDecomp;
            rootlabels[H] = clH;
        }
    }
    
    return (int) isUnique;
}

int findcomplement(const vector<int> &sub, const vector<int> &sup, const Dict &pcdict, vector<int> &complement) {
    unsigned int n = sup.size();
    
    for (unsigned int i = 0; i < n; i++)
        complement[i] = (int) ( sup[i] > sub[i] );
    
    return pcdict.lookup(complement);
}

/*********************/

int perf(const vector<vector<int> > &data, const vector<vector<int> > &propercluster, const Dict &pcdict, const vector<vector<int> > &Sv, const vector<vector<int> > &SSvG, const vector<int> &clustersize, const vector<int> &statecount, int properclustercount, int n, int m, int k) {
    int G, G1, G2;
    
    // gooddecomp records how proper cluster G is decomposed into good proper subclusters
    vector<vector<int> > gooddecomp(properclustercount);
    vector<int> good(properclustercount, 0), complement(n), tempdecomp;
    
    // specialized datasets for gbl::ENUMERATE and gbl::CHECKUNIQUE
    vector<vector<int> > verygood, rootlabels;
    vector<int> unique;
    // reusable table for O(n) set comparison in minimizeSubtrees
    //vector<int /*pair<int,int>*/ > tableDecompEqual;
    SetEqChecker tableDecompEqual;
    
    if (gbl::ENUMERATE) {
        // verygood holds subsets and every very good pair, must be filtered for only good pc complement pairs at the end
        // 0 = not subset   1 = subset  2 = very good pair
        verygood.resize(properclustercount);
    }
    
    if (gbl::CHECKUNIQUE) {
        unique.assign(properclustercount, 1);
        rootlabels.resize(properclustercount);
        //tableDecompEqual.resize(properclustercount, -1 /*make_pair(-1,-1)*/ ); 
        tableDecompEqual.reset(properclustercount);
    }
    
    if (DEBUGTRACE || VERBOSE)
        printf("Starting perf:\n");
    
    for (G = 0; clustersize[G] == 1; G++) {  // this sets all the singletons to good
        if (VERBOSE)
            printf("%d, %d\n", G, clustersize[G]); 
        
        good[G] = 1;
        gooddecomp[G].push_back(G);
        
        if (gbl::ENUMERATE)
            verygood[G].assign(G + 1, 0);
        
        if (gbl::CHECKUNIQUE) {
            unique[G] = 1;
            rootlabels[G] = data[n-1 - G];
        }
    }
    
    for ( ; G < properclustercount; G++) { // now we consider the non-singleton proper clusters
        if (gbl::ENUMERATE)
            verygood[G].assign(G + 1, 0);
        
        if (VERBOSE) {
            printf("\nStarting the consideration of proper cluster %d\n", G);
            prettyprintarray(propercluster[G], "", PPAF_PC);
        }
        
        //for (G1 = G-1; G1 >= 0; G1--) {
        for (G1 = 0; G1 < G; G1++) {
            // correct???
            //if (propercluster[G1][gbl::OUTGROUPTAXON] == 1)
            //    continue;
            
            bool gcp = false; /* a flag to note if G1 is a good, compatible, proper subset of G */
            
            if (good[G1]) {
                if (VERBOSE)
                    printf("\nproper cluster %d is good, so we must determine if it is a compatible proper subset of proper cluster %d\n", G1, G);
                
                // require that G1 is a proper subset of G
                if (clustersize[G1] < clustersize[G] && subset(propercluster[G1], propercluster[G])) {
                    if (gbl::ENUMERATE)
                        verygood[G][G1] = 1;
                    
                    // G1 is compatible with G iff G1 is exactly a union of classes in S/Sv(G)
                    if (isClassUnion(SSvG[G], propercluster[G1])) {
                        gcp = true;
                        G2 = findcomplement(propercluster[G1], propercluster[G], pcdict, complement);
                        
                        if (VERBOSE) {
                            printf("%d is a good compatible proper subset of %d. Those clusters and complement are: \n", G1, G);
                            prettyprintarray(propercluster[G], " G", PPAF_PC);
                            prettyprintarray(propercluster[G1],"G1", PPAF_PC);
                            prettyprintarray(complement,       "G2", PPAF_PC);
                        }
                    }
                    else if (VERBOSE) 
                        printf("%d is a good but NOT compatible proper subset of %d.\n", G1, G);
                }
                else if (VERBOSE) 
                    printf("%d is good but NOT a proper subset of %d.\n", G1, G);

            }
            else if (VERBOSE) 
                printf("%d is NOT good\n", G1);  
            
            if (gcp) { /* Now determine if the complement set is a good proper cluster */
                if (G2 != -1) {   /* the complement set is a proper cluster */
                    if (VERBOSE)
                        printf("\nG2 is %d, which is a proper cluster\n", G2);
                    
                    if (good[G2] == 1) {
                        if (VERBOSE)
                            printf("proper cluster %d is GOOD and is decomposed into good proper clusters %d and %d\n", G, G1, G2);
                        
                        bool firstDecomp = !good[G];
                        
                        if (firstDecomp) {
                            good[G] = 1;
                            
                            gooddecomp[G].reserve(2);
                            gooddecomp[G].push_back(G1);  /* gooddecomp records how proper cluster G is decomposed into good proper subclusters */
                            gooddecomp[G].push_back(G2);
                        }
                        
                        if (gbl::CHECKUNIQUE) {
                            tempdecomp.clear();
                            tempdecomp.push_back(G1);
                            tempdecomp.push_back(G2);
                            
                            if (G != properclustercount - 1)
                                unique[G] = minimizeSubtrees(rootlabels, gooddecomp, tempdecomp, Sv, clustersize, unique, tableDecompEqual, firstDecomp, G, G1);
                            
                            else if (!isDecompUnique(tempdecomp, unique))
                                unique[G] = false;
                            
                            if (!gbl::ENUMERATE && !unique[G])
                                break; // no need to search for more decomps
                        }
                        
                        if (gbl::ENUMERATE) {
                            verygood[G][G1] = 2;
                            continue;
                        }
                        else
                            break; /* this breaks to the end of the for G1 block */
                    
                    }
                    else if (VERBOSE)
                        printf("The complement is a proper cluster but NOT good\n");
                }
                else {
                    /* the complement of G1 in G is not a proper cluster */
                    
                    if (VERBOSE)
                        printf("The complement is NOT a proper cluster\n");
                    
                    // find equivalence classes of G2 / Sv(G,G1)
                    tempdecomp.clear();
                    
                    if ( findDecompSvGG1(pcdict, SSvG, good, propercluster, tempdecomp, G, G1, n) ) {
                        // a good decomp was found
                        
                        bool firstDecomp = !good[G];
                        
                        if (firstDecomp) {
                            good[G] = 1;
                            gooddecomp[G] = tempdecomp;
                        }
                        
                        if (gbl::CHECKUNIQUE) {
                            if (G != properclustercount - 1)
                                unique[G] = minimizeSubtrees(rootlabels, gooddecomp, tempdecomp, Sv, clustersize, unique, tableDecompEqual, firstDecomp, G, G1);
                            
                            else if (unique[G] == -1) {
                                prettyprintarray(tempdecomp, "ERROR! Top-level split failure");
                                
                                if (!isDecompUnique(tempdecomp, unique))
                                    unique[G] = false;
                            }
                            
                            if (!gbl::ENUMERATE && !unique[G])
                                break; // no need to search for more decomps
                        }
                        
                        if (gbl::ENUMERATE) {
                            verygood[G][G1] = 2;
                            continue;
                        }
                        else
                            break;
                    }
                    
                } // end of the else block for when the complement of G1 inside G is not a proper cluster
            } // end of the if (gcp) block which means that G1 is a good proper cluster and subset of G

            if (gbl::ENUMERATE)
                continue;
            else if (good[G]) // no need to consider other subsets of G, so break out of the G1 loop
                break;
        
        } /* end of the for G1 loop that enumerates proper clusters inside proper cluster G. When the complement of
             G1 in G is a good proper cluster, the break is to here. */

        if (VERBOSE) {
            if (! good[G])
                printf("Proper cluster %d is NOT good\n", G);
            else {
                if (clustersize[G] == 1)
                    printf("cluster %d is an input row and does not decompose further\n", G);
                else {
                    printf("proper cluster %d decomposes into proper clusters", G);
                    prettyprintarray(gooddecomp[G], "");
                }
            }
        }
    } /* end of the for G loop that enumerates proper clusters to check which are good */
    
    if (DEBUGTRACE) {
        if (gbl::CHECKUNIQUE) {
            printf("Uniqueness and roots:\n");
            
            for (int i = 0; i < properclustercount; i++) {
                printf("unique[%3d] = %d, ", i, unique[i]);
                prettyprintarray(rootlabels[i], "root", PPAF_DUMMY);
            }
        }
        
        prettyprintarray(gooddecomp, "gooddecomp");
    }
    
    switch (good[properclustercount -1]) {
        case 0: printf("0 The data does NOT have a perfect phylogeny\n"); return 0;
        case 1: printf("1 The data DOES have a perfect phylogeny\n");
                if (VERBOSEDATA) {
                    traceback(gooddecomp, propercluster, properclustercount, n);
                    tracebackadj(gooddecomp, propercluster, properclustercount, n);
                    taxatraceback(gooddecomp, propercluster, clustersize, properclustercount, n);
                    statetraceback(data, gooddecomp, propercluster, Sv, properclustercount, n,m);
                }
                
                //tracebackDOT(gooddecomp, propercluster, properclustercount, n, rootlabels);
                
                if (gbl::CHECKUNIQUE) {
                    if (unique[properclustercount-1])
                        printf("1 Phylogeny is unique!\n");
                    else
                        printf("0 Phylogeny is not unique!\n");
                }
                
                if (gbl::NEWICK) {
                    int outgrouppc = n-1 - gbl::OUTGROUPTAXON;
                    int comp = pcdict.lookupComplement(propercluster[outgrouppc]);
                    
                    vector<int> &topdecomp = gooddecomp[properclustercount - 1];
                    tempdecomp.clear();
                    
                    if (comp != -1 && good[comp]) {
                        // outgrouppc is a proper cluster, thus the split is sufficient
                        topdecomp.clear();
                        topdecomp.push_back(comp);
                        topdecomp.push_back(outgrouppc);
                    }
                    else if ( findDecompSvGG1(pcdict, SSvG, good, propercluster, tempdecomp, properclustercount -1, outgrouppc, n) ) {
                        // outgrouppc is not a proper cluster, thus a multiway decomposition should work
                        topdecomp = tempdecomp;
                        
                        // move outgrouppc to the end of the decomp
                        if (outgrouppc != topdecomp.back())
                            swap( *find(topdecomp.begin(), topdecomp.end(), outgrouppc), topdecomp.back() );
                    }
                    else {
                        // theoretically unreachable
                        printf("Outgroup Goodness Error!");
                    }
                    if (DEBUGTRACE) {
                        prettyprintarray(topdecomp, "Outgroup gooddecomp.back()");
                    }
                    
                    newicktraceback(data, gooddecomp, Sv, properclustercount - 1);
                    printf("\n");
                }
                
                if (gbl::ENUMERATE) {
                    enumeration(verygood, propercluster, pcdict, good, Sv, SSvG, clustersize, data, statecount, k);
                    
                    if (RUNTESTMINTREE) {
                        printf(" %d\n", treesizetraceback(gooddecomp, properclustercount - 1));
                    }
                }
                break;
    }
    
    return 1;
}

/*********************/

// tests whether v1 is a subset of v2
bool subset(const vector<int> &v1, const vector<int> &v2) {
    const unsigned int n = v1.size();
    
    for (unsigned int i = 0; i < n; i++)
        if (v1[i] > v2[i])
            return false;
    return true;
}

// O(m) compatibility check
bool compatible(const vector<int> &v1, const vector<int> &v2) {
    const unsigned int n = v1.size();
    
    for (unsigned int i = 0; i < n; i++)
        if (v1[i] != -1 && v2[i] != -1 && v1[i] != v2[i])
            return false;
    return true;
}

// O(m) compatibility check + creating a state union of the two
bool compatible(const vector<int> &Sv1, const vector<int> &Sv2, vector<int> &SvGGprime) {
    if (!compatible(Sv1, Sv2))
        return false;
    
    const unsigned int m = Sv1.size();
    
    // if compatible, create combined splitting vector Sv(G,G')
    for (unsigned int j = 0; j < m; j++) {
        if (Sv1[j] < Sv2[j])
            SvGGprime[j] = Sv2[j];
        else
            SvGGprime[j] = Sv1[j];
    }
    return true;
}

// O(n) compatibility check: G1 must be a union of classes in S/Sv(G)
bool isClassUnion(const vector<int> &eqclasses, const vector<int> &pcluster) {
    const unsigned int n = eqclasses.size();
    vector<int> owns(n,-1); // partitioning of eqclasses
    
    for (unsigned int i = 0; i < n; i++) {
        if (owns[eqclasses[i]] == -1)
            owns[eqclasses[i]] = pcluster[i];
        
        else if (owns[eqclasses[i]] != pcluster[i])
            return false;
    }
    
    return true;
}

bool hasDummyStates(const vector<int> &x) {
    const unsigned int n = x.size();
    
    for (unsigned int i = 0; i < n; i++)
        if (x[i] == -1)
            return true;
    return false;
}

// checks if x = v1 (+) v2
bool equalsStateUnion(const vector<int> &x, const vector<int> &v1, const vector<int> &v2) {
    const unsigned int n = x.size();
    
    for (unsigned int i = 0; i < n; i++) {
        // if neither equals x, not true
        if (v1[i] != x[i] && v2[i] != x[i])
            return false;
        
        // if v1 and v1 aren't equal...
        else if (v1[i] != v2[i]) {
            
            // if x is dummy, one of them must be forced thus not true
            if (x[i] == -1) 
                return false;
            
            // if x is not dummy, one of them must be dummy, else not true
            else if (v1[i] != -1 && v2[i] != -1) 
                return false;
        }
        
    }
    return true;
}


// calculate canonical labeling of root with subtrees in decomp in O(m * |decomp|)
// return value is whether the labels are compatible thus labeling is valid
bool getCanonicalLabeling(const vector<int> &decomp, const vector<vector<int> > &Sv, vector<int> &clx, int H) {
    const unsigned int d = decomp.size(), m = Sv[H].size();
    
    // start with the splitting vector of H
    clx = Sv[H];
    
    for (unsigned int j = 0; j < d; j++) {
        int G = decomp[j];
        
        for (unsigned int i = 0; i < m; i++) {
            if (Sv[G][i] != clx[i]) {
                if (clx[i] == -1)
                    clx[i] = Sv[G][i];
                
                else if (Sv[G][i] != -1)
                    return false;
            }
        }
    }
    
    return true;
}

// calculate canonical labeling of root with subtrees in decomp in O(m * |decomp|)
// trusts that the labeling will be valid
void getCanonicalLabeling_nocheck(const vector<int> &decomp, const vector<vector<int> > &Sv, vector<int> &clx, int H) {
    const unsigned int d = decomp.size(), m = Sv[H].size();
    
    // start with the splitting vector of H
    clx = Sv[H];
    
    for (unsigned int j = 0; j < d; j++) {
        int G = decomp[j];
        
        for (unsigned int i = 0; i < m; i++) {
            if (Sv[G][i] != clx[i] && clx[i] == -1) {
                clx[i] = Sv[G][i];
                break;
            }
        }
    }
}

const bool printValidDecompErrors = false;

bool checkValidDecomp(vector<int>& decomp, vector<int>& clx, const vector<vector<int> >& propercluster, const vector<vector<int> > &Sv, const vector<vector<int> > &SSvG, const vector<int> &clustersize, const Dict &pcdict, const vector<vector<int> > &verygood, int H, int n, int k) {
    clx.clear();
    
    if (DEBUGTRACE && (int) decomp.size() > k) {
        printf("Warning! Too large Decomp for %d: k=%d size=%d\n", H, k, (int) decomp.size());
        /*
        prettyprintarray(decomp, "decomp");
        prettyprintarray(propercluster[H], "H");
        
        for(unsigned int j = 0; j < decomp.size(); j++) {
            prettyprintarray(propercluster[decomp[j]], (j == 0) ? " " : "+");
        }
        */
    }
    
    // find & verify canonical labeling
    if (! getCanonicalLabeling(decomp, Sv, clx, H)) {
        if (printValidDecompErrors || DEBUGTRACE) {
            prettyprintarray(decomp, "Bad labeling      ");
        }
        return false;
    }
    
    if (! decompCoversTaxa(decomp, clustersize, H)) {
        if (((int) decomp.size()) < k) {
            if (printValidDecompErrors || DEBUGTRACE) {
                prettyprintarray(decomp, "Small, not a cover");
            }
            return false;
        }
        else
        if (printValidDecompErrors || DEBUGTRACE) {
            prettyprintarray(decomp, "Doesn't cover taxa");
        }
        
        if (findExpandedDecomp(pcdict, SSvG, decomp, H, n, clustersize, verygood) ) {
            if (printValidDecompErrors || DEBUGTRACE) {
                prettyprintarray(decomp, " ! Expanded decomp");
            }
            if (! getCanonicalLabeling(decomp, Sv, clx, H)) {
                if (printValidDecompErrors || DEBUGTRACE) {
                    prettyprintarray(decomp, "Bad labeling      ");
                }
                
                return false;
            }
        }
        else {
            if (printValidDecompErrors || DEBUGTRACE) {
                printf(" ! ERROR: no true decomposition\n");
            }
            return false;
        }
    }
    
    return true;
}


/*********************/

// stack-based search for all proper cluster subsets of H-G in O(S^(k-1) * k m)
// Ext(H,G) is all the possible canonical labelings of the root of H given G as a subtree
// each decomposition is bounded by size k since labels are necessarilly fully determined at that point
void findExtBruteHG(const vector<int> &verygoodH, const vector<vector<int> > &propercluster, const Dict &pcdict, const vector<int> &clustersize, const vector<vector<int> > &Sv, const vector<vector<int> > &SSvG, const vector<vector<int> > &verygood, const vector<vector<int> > &data, set<vector<int> > &extHG, int H, int G, int k, int comp) {
    const int n = propercluster[0].size();
    vector<int> curCluster(n,0), curDecomp, clx;
    
    int pc = findcomplement(propercluster[G], propercluster[H], pcdict, curCluster);
    if (pc == -1)
        pc = H;
    else
        pc++; // compensates for decrement in loop
    
    int taxaCount = count(curCluster.begin(), curCluster.end(), 1);
    
    curDecomp.reserve(k);
    curDecomp.push_back(G);
    
    if (DEBUGTRACE) {
        printf("\nDecomps for (%d,%d):\n", H, G);
    }
    
    while (!curDecomp.empty()) {
        pc--;
        
        if (pc < 0) {
            // push the last pc off the stack since it has no more subsets
            pc = curDecomp[curDecomp.size() - 1];
            
            for (int i = 0; i < n; i++) {
                if (propercluster[pc][i] == 1) {
                    curCluster[i] = 1;
                }
            }
            
            taxaCount += clustersize[pc];
            
            curDecomp.pop_back();
            continue;
        }
        
        if (verygoodH[pc] == 2 && subset(propercluster[pc], curCluster)) {
            // pc is a subset of the current remaining taxa
            
            if (taxaCount == clustersize[pc]) {
                // subsets cover all taxa or all characters, so this is a partition of H-G
                
                curDecomp.push_back(pc);
                
                if (checkValidDecomp(curDecomp, clx, propercluster, Sv, SSvG, clustersize, pcdict, verygood, H, n, k)) {
                    
                    if (DEBUGTRACE) {
                        prettyprintarray(curDecomp, "   d");
                    }
                    
                    if (gbl::ALLOWINTERNAL) {
                        bool internal = false;
                        int i;
                        // since in descending order, all taxa will be at the end, will be G, or will be the parent of H
                        
                        for (i = (int) curDecomp.size() - 1; i >= 0 && curDecomp[i] < n; i--)
                            if (compatible(clx, data[n-1 - curDecomp[i]])) {
                                extHG.insert( data[n-1 - curDecomp[i]] );
                                internal = true;
                            }
                        
                        if (i != -1 && G < n && compatible(clx, data[n-1 - G]) ) {
                            extHG.insert( data[n-1 - G] );
                            internal = true;
                        }
                        
                        if (comp != -1 && compatible(clx, data[n-1 - comp])) {
                            extHG.insert( data[n-1 - comp] );
                            internal = true;
                        }
                        
                        if (!internal) {
                            extHG.insert(clx);
                        }
                    }
                    else {
                        extHG.insert(clx);
                    }
                }
                
                curDecomp.pop_back();
            }
            else /*if (!MISLIMITK || (int) curDecomp.size() < k)*/ {
                // if not the last one, subtract it from the current cluster and push onto stack
                
                for (int i = 0; i < n; i++) {
                    if (propercluster[pc][i] == 1) {
                        curCluster[i] = 0;
                    }
                }
                taxaCount -= clustersize[pc];
                
                curDecomp.push_back(pc);
            }
        }
    }
}

void findExtBrute(const vector<vector<int> > &propercluster, const vector<vector<int> > &Sv, const vector<vector<int> > &SSvG, const vector<int> &clustersize, const Dict &pcdict, const vector<vector<int> > &verygood, const vector<vector<int> > &data, vector<map<int, set<vector<int> > > > &Ext, int k) {
    printf("\n");
    int n = data.size();
    int properclusterend = propercluster.size() - 1;
    
    for (int H = 0; H < properclusterend; H++) {
        int comp = -1;
        
        if (gbl::ALLOWINTERNAL && clustersize[H] == n-1) {
            comp = pcdict.lookupComplement(propercluster[H]);
        }
        
        for (int G = 0; G < H; G++) {
            if (verygood[H][G] == 2) {
                findExtBruteHG(verygood[H], propercluster, pcdict, clustersize, Sv, SSvG, verygood, data, Ext[H][G], H, G, k, comp);
            }
        }
    }
}

void findExtMIS( const vector<vector<int> > &propercluster, const vector<vector<int> > &Sv, const vector<vector<int> > &SSvG, const vector<int> &clustersize, const Dict &pcdict, const vector<vector<int> > &verygood, const vector<vector<int> > &data, vector<map<int, set<vector<int> > > > &Ext, int k) {
    vector<vector<int> > partitions;
    vector<int> partners;
    const int properclusterend = (int) propercluster.size() - 1, n = propercluster[0].size();
    
    for (int H = 0; H < properclusterend; H++) {
        int comp;
        if (gbl::ALLOWINTERNAL && clustersize[H] == n-1) {
            comp = pcdict.lookupComplement(propercluster[H]);
        }
        else {
            comp = -1;
        }
        
        partitions.clear();
        partners.clear();
        
        // gather all very good partners of H
        for (int G = 0; G < H; G++) {
            if (verygood[H][G] == 2) {
                partners.push_back(G);
            }
        }
        
        if (DEBUGTRACE) {
            printf("\nDecomps for %d:\n", H);
            prettyprintarray(partners, "partners");
        }
        
        // find maximal independent sets of them
        enumMaxIndSets(partners, propercluster, Sv, clustersize, partitions, H, k);
        
        vector<int> clx;
        
        // for every partition found...
        unsigned int partitionsize = partitions.size();
        for (unsigned int p = 0; p < partitionsize; p++) {
            vector<int> decomp = partitions[p];
            
            if (! checkValidDecomp(decomp, clx, propercluster, Sv, SSvG, clustersize, pcdict, verygood, H, n, k) ) {
                continue;
            }
            
            // check for internal and insert labeling
            if (gbl::ALLOWINTERNAL) {
                bool internal = false;
                
                for (unsigned int i = 0; i < decomp.size(); i++) {
                    int G = decomp[i];
                    
                    if (G < n && compatible(clx, data[n-1 - G])) {
                        // G is internal, thus insert taxon as ext of all partitions
                        for (unsigned int t = 0; t < decomp.size(); t++) {
                            Ext[H][ decomp[t] ].insert( data[n-1 - G] );
                        }
                        
                        internal = true;
                    }
                }
                
                if (comp != -1 && compatible(clx, data[n-1 - comp])) {
                    // insert taxon as ext of all partitions
                    for (unsigned int t = 0; t < decomp.size(); t++) {
                        Ext[H][ decomp[t] ].insert( data[n-1 - comp] );
                    }
                    
                    internal = true;
                }
                
                if (!internal) {
                    // insert taxon as ext of all partitions
                    for (unsigned int t = 0; t < decomp.size(); t++) {
                        Ext[H][ decomp[t] ].insert( clx );
                    }
                }
            }
            else {
                for (unsigned int t = 0; t < decomp.size(); t++) {
                    Ext[H][ decomp[t] ].insert( clx );
                }
            }
            
            if (DEBUGTRACE) {
                prettyprintarray(decomp, "Good decomposition");
                prettyprintarray(clx,    " > Canonical Label");
            }
        }
    }
}

// filter verygood to remove pairs (H,G) where H or S-H are not good
// if H is very good, good[H] is set to 2
// O(2L + (L + n) * (number of not good pc's))
void identifyVeryGoodness(vector<vector<int> > &verygood, vector<int> &good, const Dict &pcdict, const vector<vector<int> > &propercluster) {
    int properclusterend = propercluster.size() - 1;
    
    // assume all good H's are verygood
    for (int H = 0; H < properclusterend; H++) {
        if (good[H]) {
            good[H] = 2;
        }
    }
    
    // look at not good pc's for counterexamples
    for (int H = 0; H < properclusterend; H++) {
        if (good[H] == 0) {
            int comp = pcdict.lookupComplement(propercluster[H]);
            
            if (comp != -1 && good[comp]) {
                // H is not good, so comp is not very good
                good[comp] = 1;
                
                for (int G = 0; G < comp; G++) {
                    if (verygood[comp][G] == 2) {
                        verygood[comp][G] = 1;
                    }
                }
            }
        }
    }
}

int enumeration(vector<vector<int> > &verygood, const vector<vector<int> > &propercluster, const Dict &pcdict, vector<int> &good, const vector<vector<int> > &Sv, const vector<vector<int> > &SSvG, const vector<int> &clustersize, const vector<vector<int> > &data, const vector<int> &statecount, int k) {
    
    const int n = data.size(), m = data[0].size();
    const int properclustercount = (int) propercluster.size(), properclusterend = properclustercount-1;
    
    vector< vector<vector<int> > > ConnVect, VpVect;
    
    // ConnVect and VpVect created in this scope:
    {
        int H;
        
        vector<map<int, set<vector<int> > > > Ext(properclustercount); // a vector over H's of maps from G to Ext(H,G) for only very good pairs
        vector< set<vector<int> > > Conn(properclustercount), Vp(properclustercount);
        vector<int> decomp, temp, SvHG(m);
        vector<vector<int> > prevTaxa, nextTaxa;
        
        //set<vector<int> > extHG; // unused
        set<vector<int> >::iterator x, y;
        
        identifyVeryGoodness(verygood, good, pcdict, propercluster);
        
        if (DEBUGTRACE) {
            prettyprintarray(verygood, "verygood");
            printf("Finding Ext(H,G) for every good pair:");
        }
        
        #if MISEXT
            findExtMIS(propercluster, Sv, SSvG, clustersize, pcdict, verygood, data, Ext, k);
        #else
            findExtBrute(propercluster, Sv, SSvG, clustersize, pcdict, verygood, data, Ext, k);
        #endif
        
        
        // create adj lists for the case where H/x must be computed
        createStateGraph(data, statecount, prevTaxa, nextTaxa);
        
        if (DEBUGTRACE)
            printf("\nFinding Conn(H) and Vp(G) ...\n");
        
        for (H = 0; clustersize[H] == 1; H++) {
            // Conn(H) for singletons are the base case of the dynamic program, compute separately
            // since Conn[H] is the set of all canonical labelings of the roots of minimal perfect subphylogenies for H
            // if H is a single taxa, the only canonical labeling is itself
            
            Conn[H].insert(data[n - 1 - H]);
            
            if (n == 2) {
                // corner case: a two-taxa tree doesn't have any proper cluster supersets
                int comp = 1 - H;
                Vp[H].insert(data[n - 1 - comp]);
                continue;
            }
            
            if (gbl::ALLOWINTERNAL) {
                int comp = pcdict.lookupComplement(propercluster[H]);
                
                if (comp == -1) {
                    // if no complement (aka is a pseudo proper cluster), then shortcut: Vp is only itself
                    Vp[H] = Conn[H];
                }
                else {
                    // Vp's must be incompatible or equal to itself
                    for (int G = H; G < properclusterend; G++)
                        if (verygood[G][H] == 2) 
                            for (x = Ext[G][H].begin(); x != Ext[G][H].end(); ++x)
                                if ( !compatible(*x, data[n-1 - H]) || *x == data[n-1 - H] )
                                    Vp[H].insert(*x);
                }
            }
            else {
                // taxon can't be internal, so all parents are valid
                for (int G = H; G < properclusterend; G++)
                    if (verygood[G][H] == 2)
                        Vp[H].insert(Ext[G][H].begin(), Ext[G][H].end());
            }
        }
        
        
        for ( ; H < properclusterend; H++) {
            if (good[H] != 2)
                continue;
            
            int G;
            
            // Conn(H)
            for (G = 0; G < H; G++) {
                if (verygood[H][G] != 2)
                    continue;
                
                // for every vector x in Ext(H,G)
                for (x = Ext[H][G].begin(); x != Ext[H][G].end(); ++x) {
                    
                    if (Conn[H].count(*x) && G >= n)
                        continue; // x is already in Conn(H), skip
                    
                    bool inConnH;
                    
                    if ( hasDummyStates(*x) ) {
                        int comp = pcdict.lookupComplement(propercluster[G], propercluster[H]);
                        
                        if (comp != -1) {
                            // when x has dummy states, can put x in Conn(H) iff:
                            // x in Vp(G)
                            // x in Vp(H-G) or x = Sv(H,G) + y for some y in Conn(H-G)
                            
                            // check if x in Vp(G) and Vp(H-G)
                            inConnH = Vp[G].count(*x);
                            
                            if (!inConnH)
                                continue;
                                
                            inConnH = Vp[comp].count(*x);
                            
                            // check if x = Sv(G,H) + y for some y in Conn(H-G)
                            if (!inConnH && compatible(Sv[G], Sv[H], SvHG)) {
                                inConnH = false;
                                for (y = Conn[comp].begin(); y != Conn[comp].end(); ++y) {
                                    inConnH = equalsStateUnion(*x, *y, SvHG);
                                    if (inConnH)
                                        break;
                                }
                            }
                            
                            if (inConnH)
                                Conn[H].insert(*x); // x satisfied all requirements
                        }
                    }
                    else {
                        // when x has no dummy states, x in Conn(H) iff:
                        // x in Vp(G) for every G in H/x
                        
                        int maxC = findEqClassesMulti(*x, temp, data, prevTaxa, nextTaxa);
                        
                        decomp.clear();
                        
                        #if UseQTable
                            inConnH = pcdict.lookupParallel(temp, decomp, n, maxC, clustersize);
                        #else
                            inConnH = pcdict.lookupParallel(temp, decomp, n, maxC);
                        #endif
                        
                        if (!inConnH)
                            continue;
                        
                        for (vector<int>::iterator t = decomp.begin(); t != decomp.end(); ++t) {
                            if (*t >= H || verygood[H][*t] == 0)
                                continue;
                            
                            if (Vp[*t].count(*x) == 0) {
                                inConnH = false;
                                break;
                            }
                        }
                        
                        if (inConnH)
                            Conn[H].insert(*x);
                    }
                }
                
            }
            
            // Vp(H)
            
            if (clustersize[H] == n-1) {
                // if complement is a singleton, Vp can only be that taxon
                int comp = pcdict.lookupComplement(propercluster[H]);
                
                if (gbl::ALLOWINTERNAL) {
                    x = Conn[comp].begin();
                    
                    for (y = Conn[H].begin(); y != Conn[H].end(); ++y)
                        if(!compatible(*y, *x)) {
                            Vp[H] = Conn[comp];
                            break;
                        }
                }
                else
                    Vp[H] = Conn[comp];
                
                continue;
            }
            
            //  for every very good pair (G,H)
            //      for every vector x in Ext(G,H)
            //          if x is incompatible with any vector y in Conn(H), x is in Vp(H)
            
            for ( ; G < properclusterend; G++)
                if (verygood[G][H] == 2)
                    for (x = Ext[G][H].begin(); x != Ext[G][H].end(); ++x)
                        for (y = Conn[H].begin(); y != Conn[H].end(); ++y)
                            if ( !compatible(*x, *y)) {
                                Vp[H].insert(*x);
                                break;
                            }
        }
        
        ConnVect.resize(properclustercount);
        VpVect.resize(properclustercount);
        
        // move contents of Conn and Vp into vectors for indexing
        for (H = 0; H < properclustercount; H++) {
            ConnVect[H].assign(Conn[H].begin(), Conn[H].end()); // unfortunately must copy all data, set does not allow for move-out
            VpVect[H].assign(Vp[H].begin(), Vp[H].end());
        }
        
        if (DEBUGTRACE) {
            for (H = 0; H < properclusterend; H++) {
                if (good[H] != 2)
                    continue;
                
                int Extsize = 0;
                for (int G = 0; G < H; G++) {
                    if (verygood[H][G] == 2 && Ext[H][G].size() > 0) {
                        Extsize += Ext[H][G].size();
                        
                        printf("Ext(%d,%d): \n", H, G);
                        int j = 0;
                        
                        for (x = Ext[H][G].begin(); x != Ext[H][G].end(); ++x) {
                            printf(" [%d] = { ", j++);
                            for (int i = 0; i < m; i++) {
                                if ((*x)[i] == -1)
                                    printf("* ");
                                else
                                    printf("%d ", (*x)[i]);
                            }
                            printf("}\n");
                        }
                        
                    }
                }
                
                if (VpVect[H].size() > 0) {
                    printf("Vp(%d): \n", H);
                    prettyprintarray(VpVect[H], " ", PPAF_DUMMY);
                }
                
                if (ConnVect[H].size() > 0) {
                    printf("Conn(%d): \n", H);
                    prettyprintarray(ConnVect[H], " ", PPAF_DUMMY);
                }
                
                printf("Ext(%3d) has %3d vectors.  Conn(%3d) has %3d vectors.  Vp(%3d) has %3d vectors.\n", H, Extsize, H, (int) ConnVect[H].size(), H, (int) VpVect[H].size());
                
            }
            fflush(stdout);
        }
    }
    
    // construct DAG, returns tree count
    return constructDAG(ConnVect, VpVect, verygood, good, pcdict, propercluster, Sv, data, statecount, clustersize, properclustercount, k);
}

// Find Maximal Independent Sets
void enumMaxIndSets(const vector<int> &partners, const vector<vector<int> > &propercluster, const vector<vector<int> > &Sv, const vector<int>& clustersize, vector<vector<int> > &output, int H, int k) {
    const unsigned int n = propercluster[0].size(), psize = partners.size();
    unsigned int edgeCount = 0;
    
    //prettyprintarray(partners, "partners");
    
    vector<vector<int> > graph(psize), rawmis;
    
    // create graph from partners and find maximal independent sets
    // edges == intersecting proper clusters
    
    for (unsigned int G1 = 0; G1 < psize; G1++) {
        
        if (DEBUGTRACE) {
            printf("%d;\n", partners[G1]);
        }
        
        for (unsigned int G2 = G1+1; G2 < psize; G2++) {
            bool addEdge = false;
            
            if (!compatible(Sv[ partners[G1] ],Sv[ partners[G2] ])) {
                addEdge = true;
            }
            else {
                for (unsigned int i = 0; i < n; i++)
                    if (propercluster[ partners[G1] ][i] == 1 && propercluster[ partners[G2] ][i] == 1) {
                        addEdge = true;
                        break;
                    }
            }
            
            if (addEdge) {
                graph[G1].push_back(G2);
                graph[G2].push_back(G1);
                edgeCount++;
                
                if (DEBUGTRACE) {
                    printf("%d -- %d;\n", partners[G1], partners[G2]);
                }
            }
        }
    }
    
    
    if (edgeCount > 0) {
        // find and store every MIS
        
        if (MISLIMITK) {
            KMAXIS(graph, k, rawmis);
        }
        else {
            MIS(graph, rawmis);
        }
        
        output.resize(rawmis.size());
        
        for (unsigned int i = 0; i < rawmis.size(); i++) {
            for (unsigned int j = 0; j < rawmis[i].size(); j++)
                if (rawmis[i][j] == 0)
                    output[i].push_back(partners[j]);
        }
    }
    else if (psize > 0) {
        // only one maximal independent set: every partner
        
        output.push_back(partners);
    }
}

// find Maximal Independent Sets of a graph of pairs
// NOT used to find limited-sized sets
void enumMaxIndSetsPairs(const vector<pair<int,int> > &partners, const vector<vector<int> > &propercluster, const vector<vector<int> > &Sv, const vector<int>& clustersize, vector<PPDAG::ProductNode> &output, int H, int k) {
    const unsigned int n = propercluster[0].size(), psize = partners.size();
    unsigned int edgeCount = 0;
    
    vector<vector<int> > graph(psize), rawmis;
    
    // create graph from partners and find maximal independent sets
    // edges == intersecting proper clusters
    
    for (unsigned int G1 = 0; G1 < psize; G1++) {
        for (unsigned int G2 = G1+1; G2 < psize; G2++) {
            bool addEdge = false;
            
            if (!compatible(Sv[ partners[G1].first ],Sv[ partners[G2].first ])) {
                addEdge = true;
            }
            else {
                for (unsigned int i = 0; i < n; i++) {
                    if (propercluster[ partners[G1].first ][i] == 1 && propercluster[ partners[G2].first ][i] == 1) {
                        addEdge = true;
                        break;
                    }
                }
            }
            
            if (addEdge) {
                graph[G1].push_back(G2);
                graph[G2].push_back(G1);
                edgeCount++;
            }
        }
    }

    
    if (edgeCount > 0) {
        // find and store every MIS
        MIS(graph, rawmis);
        
        output.resize(rawmis.size());
        
        for (unsigned int i = 0; i < rawmis.size(); i++) {
            for (unsigned int j = 0; j < rawmis[i].size(); j++)
                if (rawmis[i][j] == 0)
                    output[i].decomp.push_back(partners[j]);
        }
    }
    else if (psize > 0) {
        // only one maximal independent set: every partner
        
        PPDAG::ProductNode pn(partners);
        
        output.push_back(pn);
    }
}


int constructDAG(vector< vector<vector<int> > > &Conn, vector< vector<vector<int> > > &Vp, const vector<vector<int> > &verygood, const vector<int> &good, const Dict &pcdict, const vector<vector<int> > &propercluster,
                 const vector<vector<int> > &Sv, const vector<vector<int> > &data, const vector<int> &statecount, const vector<int> &clustersize, int properclustercount, int k) {
    const unsigned int n = data.size(), m = data[0].size();
    
    vector<vector< vector<PPDAG::ProductNode> > > pcGraph(properclustercount);
    // for each pc H, each x in Conn[H] has a list of Product Nodes that were maximal independent sets of the graph L(H,x)
    
    vector<int> SvHG(m), temp;
    vector<pair<int,int> > partners, partnerZ;
    map<int,int> partnerVp;
    
    PPDAG dag(propercluster, pcdict); // initialize DAG
    
    vector<bool> inDag(dag.maxPC + 1, false);
    
    inDag[dag.maxPC] = true;
    Vp[dag.maxPC].push_back( vector<int>(m,-2) );
    Conn[dag.maxPC].push_back( data[n-1-dag.rootTaxaPC] );
    
    // first, search for all proper clusters where rootTaxaPC is a valid parent
    // the MIS of the graph is the Product Node children of the root
    
    int H = dag.maxPC;
    dag.resizeRow(H, 1);
    
    for (int G1 = 0; G1 < H; G1++)
        if (verygood[H][G1] == 2)
            for (int y = 0; y < (int) Vp[G1].size(); y++)
                if (Vp[G1][y] == Conn[H][0]) {
                    partners.push_back( pair<int,int>(G1,y) );
                    inDag[G1] = true;
                    break;
                }
    
    enumMaxIndSetsPairs(partners, propercluster, Sv, clustersize, dag.sumNodes[H][0].children, H, k);
    
    H--;
    
    for ( ; H >= 0; H--) {
        if (!inDag[H])
            continue;
        
        pcGraph[H].resize(Conn[H].size());

        dag.resizeRow(H, Vp[H].size());
        
        // look for minimal subphylogenies of H with roots y in Vp(H)
        // by finding all connections x of H incompatible with y
        
        for (unsigned int y = 0; y < Vp[H].size(); y++) {
            for (unsigned int x = 0; x < Conn[H].size(); x++) {
                if ( compatible(Conn[H][x], Vp[H][y]) && (gbl::ALLOWINTERNAL || clustersize[H] < (int) n-1) ) // allows for rootTaxaPC to be compatible
                    continue;
                
                partners.clear();
                
                // find all very good partners of H and their Vp index of x
                for (int G1 = 0; G1 < H; G1++)
                    if (verygood[H][G1] == 2) {
                        
                        int vp = binary_find(Vp[G1], Conn[H][x]);
                        
                        if (vp != -1) {
                            partners.push_back( pair<int,int>(G1,vp) );
                            inDag[G1] = true;
                        }
                    }
                
                if (hasDummyStates(Conn[H][x])) {
                    // for each partner, check if it and its complement work as a decomposition
                    const int psize = partners.size();
                    
                    for (int G1 = 0; G1 < psize; G1++) {
                        int G2 = pcdict.lookupComplement(propercluster[partners[G1].first], propercluster[H]);
                        
                        if (G2 == -1)
                            continue;
                        
                        int vp = binary_find(Vp[G2], Conn[H][x]);
                        
                        if (vp != -1) {
                            // create product node off of (H,y) for (G1,G2;x)
                            
                            PPDAG::ProductNode pn;
                            if (partners[G1].first < G2) {
                                pn.decomp.push_back( partners[G1] );
                                pn.decomp.push_back( pair<int,int>(G2, vp) );
                            }
                            else {
                                pn.decomp.push_back( pair<int,int>(G2, vp) );
                                pn.decomp.push_back( partners[G1] );
                            }
                            
                            dag.insertProductNode(pn, H, y);
                            
                            inDag[G2] = true;
                        }
                        
                        compatible(Sv[ partners[G1].first ], Sv[H], SvHG);
                        
                        // look for z in Conn[G2] such that x = z (+) Sv(H,G1)
                        for (unsigned int z = 0; z < Conn[G2].size(); z++) {
                            if (equalsStateUnion(Conn[H][x], Conn[G2][z], SvHG)) {
                                
                                if (pcGraph[G2].size() <= z)
                                    pcGraph[G2].resize(z+1);
                                
                                if (pcGraph[G2][z].empty()) {
                                    // create new unenumerated MIS for (G2,z) and make product nodes labeled with x
                                    
                                    partnerZ.clear();
                                    for (int i = 0; i < G2; i++)
                                        if (verygood[G2][i] == 2) {
                                            int vp = binary_find(Vp[i], Conn[G2][z]);
                                            
                                            if (vp != -1) {
                                                partnerZ.push_back( pair<int,int>(i,vp) );
                                                inDag[i] = true;
                                            }
                                        }
                                    
                                    enumMaxIndSetsPairs(partnerZ, propercluster, Sv, clustersize, pcGraph[G2][z], G2, k);
                                }
                                
                                // create product node off of (H,y) for each precomputed maximal independent set in pcGraph[G2][z]
                                // must add G1 to each set and relabel their Vp's to x
                                
                                partnerVp.clear(); // partnerVp is a "cache" of preconverted valid parent values
                                PPDAG::ProductNode pn;
                                
                                for (unsigned int i = 0; i < pcGraph[G2][z].size(); i++) {
                                    pn = pcGraph[G2][z][i];
                                    
                                    // check everything to relabel it to x, compute if not seen yet
                                    for (unsigned int j = 0; j < pn.decomp.size(); j++) {
                                        if (partnerVp.count(pn.decomp[j].first) == 0) {
                                            pn.decomp[j].second = binary_find(Vp[ pn.decomp[j].first ], Conn[H][x]);
                                            partnerVp.insert(pn.decomp[j]);
                                        }
                                        else
                                            pn.decomp[j].second = (partnerVp.find(pn.decomp[j].first))->second;
                                    }
                                    
                                    // insert G1 into correct sorted position
                                    pn.decomp.insert( lower_bound(pn.decomp.begin(), pn.decomp.end(), partners[G1]), partners[G1]);
                                    
                                    dag.insertProductNode(pn, H, y);
                                }
                                
                            }
                        } // end for every z
                    }
                } // end if dummy states
                
                else {
                    if (pcGraph[H][x].empty()) // enumerate MIS if not precomputed
                        enumMaxIndSetsPairs(partners, propercluster, Sv, clustersize, pcGraph[H][x], H, k);
                    
                    dag.insertManyProductNodes(pcGraph[H][x], H, y);
                }
            }
        }
    }
    
    dag.finalize(Conn, Vp, pcdict, propercluster, data);
    
    // no longer require these datasets, clear to save memory
    pcGraph.clear();
    Conn.clear();
    Vp.clear();
    
    printf("%d Minimal Perfect Phylogeny Count\n", dag.getCount());
    
    
    // output the DAG in .dot format or custom format as specified by cmd args
    if (gbl::DAGDOTOUTF)
        dag.makeDOT();
    
    if (gbl::DAGPRETTYF) {
        ofstream outf(gbl::DAGPRETTYF);
        dag.prettyprint(outf);
    }
    else if (DEBUGTRACE)
        dag.prettyprint(cout);
    
    if (gbl::DAGSTATISTICS)
        takeDagStatistics(dag, propercluster, pcdict, data, statecount);
    
    if (RUNTESTMINTREE) {
        printf("%d ", dag.minNodeCount());
    }
    
    return dag.getCount();
}

void takeDagStatistics(const PPDAG &dag, const vector<vector<int> > &propercluster, const Dict &pcdict, const vector<vector<int> > &data, const vector<int> &statecount) {
    // at this point, the DAG should be fully constructed and ready for use
    
    PPDAG::tree_iterator mpp(dag);
    
    //if (gbl::DAGDOTOUTF)
    //    dag.makeDOTForTrees();
    
    if (DEBUGTRACE) {
        // verify that every tree is a minimal perfect phylogeny
        for (mpp.gotoIndex(0); !mpp.atend(); ++mpp) {
            if (!mpp.isMPP(data, statecount)) {
                printf("0 NOT all mpp are minimal perfect phylogenies:\n");
                mpp.printStack();
                prettyprintarray(dag.nodeLabels, "nodelabels", PPAF_DUMMY);
                
                ofstream badoutf("MPP_bad.dot", ios_base::out);
                mpp.makeTreeDOTLabeled(badoutf);
                
                mpp.printTreeDetailed();
                
                return;
            }
        }
        
        printf("1 All mpp are indeed minimal perfect phylogenies.\n");
        
        vector<int> curSplits;
        
        for (mpp.gotoIndex(0); !mpp.atend(); ++mpp) {
            mpp.getTreeSplits(curSplits);
            printf("mpp[%3d", *mpp);
            prettyprintarray(curSplits, "]");
        }
    }
    
    printf("%d Min node count\n", dag.minNodeCount() );
    printf("%d Max node count\n", dag.maxNodeCount() );
    
    // find supports for splits and nodes with dynamic programming
    
    vector<vector<int> > snSupport;
    vector<int> splitSupport, nodeSupport;
    
    dag.getSupportStats(snSupport, splitSupport, nodeSupport);
    
    // find the tree with maximum support
    pair<int,int> maxTreeDP = dag.getMaxSupportTree(splitSupport);
    
    mpp.gotoIndex(maxTreeDP.second);
    
    printf("Split Supports:\n");
    for (unsigned int i = 0; i < splitSupport.size() - 1; i++) {
        if (splitSupport[i] != 0)
            printf("  (%3u, %3d) = %d\n", i, pcdict.lookupComplement(propercluster[i]), splitSupport[i]);
    }
    
    printf("Node Label Supports:\n");
    for (unsigned int i = 0; i < nodeSupport.size(); i++) {
        printf("  nodeSupport[%2u] = %d", i, nodeSupport[i]);
        prettyprintarray(dag.nodeLabels[i], ", label", PPAF_DUMMY);
    }
    
    printf("Tree with maximum split support of %d is #%d.\n", maxTreeDP.first, maxTreeDP.second);
    
    mpp.printTreeDetailed();
    //mpp.makeTreeDOTLabeled("MPP_maxSupport.dot", &splitSupport);
    mpp.printNewick(&splitSupport);
    printf("\n");
    
    // calculate robinson-foulds distances
    dag.printRobinsonFoulds();
}


/**************************/

/*
void tracebackDOT(const vector<vector<int> > &gooddecomp, const vector<vector<int> > &propercluster, int properclustercount, int n, const vector<vector<int> > &rootlabels)
{
    int i, curPC, head = 0, tail = 0;
    vector<int> childNum(n,-1), stack;
    vector<int>::const_iterator j;
    
    ofstream outf(gbl::DPDOTOUTF);
    
    stack.reserve(n);
    stack.push_back(properclustercount - 1);
    
    outf << "digraph PerfectPhylogeny {" << endl;
    
    while ( !stack.empty() ) {
        curPC = stack[tail];
        if (gooddecomp[curPC].size() == 1) {
            // leaf, thus end of a path, so print it
            outf << stack[head];
            
            for (i = head+1; i < tail; i++)
                outf << " -> " << stack[i];
            
            outf << " -> " << curPC << ';' << endl;
            
            if (gbl::CHECKUNIQUE) {
                for (i = head+1; i < tail; i++) {
                    outf << stack[i] << " [ label = \"";
                    for (j = rootlabels[ stack[i] ].begin(); j != rootlabels[ stack[i] ].end(); j++)
                        outf << *j << ' ';
                    outf << "\" ];" << endl;
                }
                
                outf << curPC << " [ label = \"";
                for (j = rootlabels[ curPC ].begin(); j != rootlabels[ curPC ].end(); j++)
                    outf << *j << ' ';
                outf << "\" color=\"green\" ];" << endl;
            }    
            else
                outf << curPC << " [ color=\"green\" ];" << endl;
            
            outf << endl;
            
            stack.pop_back();
            tail--;
            head = tail;
        }
        else {
            // interior node
            childNum[tail]++;
            if ( childNum[tail] < (int) gooddecomp[curPC].size() ) {
                // has unexplored edges, push child onto stack
                stack.push_back( gooddecomp[curPC][ childNum[tail] ] );
                tail++;
            }
            else {
                // all children explored, fall back to parent
                childNum[tail] = -1;
                stack.pop_back();
                head = --tail;
            }
        }
        
    }
    
    outf << properclustercount - 1 << " [ label = \"root\" ];" << endl;
    
    outf << '}' << endl;
    
    outf.close();
}
*/

void traceback(const vector<vector<int> > &gooddecomp, const vector<vector<int> > &propercluster, int properclustercount, int n)
{
    int ii, head = 0, tail = 0, endoflevel;
    vector<int> queue;
    vector<int>::const_iterator j;

    printf("In traceback \n");
    queue.push_back(properclustercount - 1);
    endoflevel = queue[head]; /* endoflevel is the last propercluster in a level of the tree */
    printf("%d\n", properclustercount-1);

    while (head <= tail) {
        ii = queue[head]; /* propercluster ii needs to get expanded */
        printf("%d:", ii);
        head++;
        for (j = gooddecomp[ii].begin(); j != gooddecomp[ii].end(); ++j) { 
            printf(" %d,", *j);
            if (ii != *j) {
                tail++;
                queue.push_back(*j);
            }
        }
        printf("|  ");

        if (ii == endoflevel) {
            printf("\n\n");
            endoflevel = queue[tail];
        }
    }
}

/**************/

void tracebackadj(const vector<vector<int> > &gooddecomp, const vector<vector<int> > &propercluster, int properclustercount, int n)
{
    int ii, head = 0, tail = 0, endoflevel;
    vector<int> queue;
    vector<int>::const_iterator j;


    printf("(Non-symmetric) Adjacency list \n");
    queue.push_back(properclustercount - 1);
    endoflevel = queue[head]; /* endoflevel is the last propercluster in a level of the tree */

    while (head <= tail) {
        ii = queue[head]; /* propercluster ii needs to get expanded */
        printf("%d:", ii);
        head++;
        for (j = gooddecomp[ii].begin(); j != gooddecomp[ii].end(); ++j) { 
            printf(" %d,", *j);
            if (ii != *j) {
                tail++;
                queue.push_back(*j);
            }
        }
        printf("\n");

        if (ii == endoflevel) {
            printf("\n\n");
            endoflevel = queue[tail];
        }
    }
    printf("\n\n");
}



/*****************/

void taxatraceback(const vector<vector<int> > &gooddecomp, const vector<vector<int> > &propercluster, const vector<int> &clustersize, int properclustercount, int n)
{
    int ii, head = 0, tail = 0, endoflevel, k, properclusterindex;
    vector<int> queue;
    vector<int>::const_iterator j;

    printf("In taxatraceback \n");


    queue.push_back(properclustercount - 1);
    endoflevel = queue[head]; /* endoflevel is the last propercluster in a level of the tree */

    printf ("%d: ", endoflevel);
    for (k = 0; k < n; k++) {
        if (propercluster[endoflevel][k] == 1) {
            printf("%d ", k);
        }
    }
    printf("\n");

    while (head <= tail) {
        ii = queue[head]; /* propercluster ii needs to get expanded */
        /* if (clustersize[ii] > 1) { */
        printf("%d:", ii); 
        head++;
        for (j = gooddecomp[ii].begin(); j != gooddecomp[ii].end(); ++j) { 
            properclusterindex = *j;

            for (k = 0; k < n; k++) { /* print out the taxa that make up the good proper cluster */
                if (propercluster[properclusterindex][k] == 1) {
                    printf("%d ", k);
                }
            }
            printf ("; ");

            if (ii != *j) {
                tail++;
                queue.push_back(*j);
            }
        }
        printf("|  ");
        /* } */

        if (ii == endoflevel) {
            printf("\n\n");
            endoflevel = queue[tail];
        }
    }

}




/***************/

void statetraceback(const vector<vector<int> > &data, const vector<vector<int> > &gooddecomp, const vector<vector<int> > &propercluster, const vector<vector<int> > &Sv, int properclustercount, int n, int m)
{
    int head = 0, tail = 0, endoflevel, ii, jj, firstcluster;
    vector<int> queue;//, SvGGprime(m);
    vector<int>::const_iterator j;

    printf("In statetraceback \n");
    queue.push_back(properclustercount - 1);
    endoflevel = queue[head]; /* endoflevel is the last propercluster in a level of the tree */
    printf("%d\n", properclustercount-1);

    while (head <= tail) {
        ii = queue[head]; /* propercluster ii needs to get expanded */
        printf("%d", ii);
        firstcluster = gooddecomp[ii][0];
        if (firstcluster != ii) {
            if (!compatible(Sv[ii], Sv[firstcluster] /* , SvGGprime */)) {
                printf("Something wrong in statetraceback\n");
            }
            
            printf("(");
            for (jj = 0; jj < m; jj++) {
                /*if (SvGGprime[jj] != -1)  {
                printf("%d", SvGGprime[jj]);  */
                if (Sv[ii][jj] != -1) { /*   Change Aug. 3 to just print Sv not SvGGprime   */
                    printf("%d", Sv[ii][jj]);  
                }
                else printf("*");
            }
            printf("):");
        }
        else {   /* these are the singleton `proper clusters', i.e., the original data rows. Unfortunately, the ii'th
                proper cluster does not correspond to the ii'th data row, because of the way the radix sort of
                the proper clusters sorted them. In fact, the ordering is exactly reversed, so we use kk to
                compute the correct index into the data matrix. */
            printf("(");
            //kk = n - ii -1;

                for (jj = 0; jj < m; jj++) {
                    if (Sv[ii][jj] != -1) {    
                        printf("%d", Sv[ii][jj]);  
                    }
                    else printf("*");
                    /* printf("%d", data[kk][jj]); */
                }
                printf("):");
        }

        head++;
        for (j = gooddecomp[ii].begin(); j != gooddecomp[ii].end(); ++j) { 
            printf(" %d,", *j);
            if (ii != *j) {
                tail++;
                queue.push_back(*j);
            }
        }
        printf("|  ");

        if (ii == endoflevel) {
            printf("\n\n");
            endoflevel = queue[tail];
        }
    }
}

void newicktraceback(const vector<vector<int> > &data, const vector<vector<int> > &gooddecomp, const vector<vector<int> > &Sv, int H) {
    // http://evolution.genetics.washington.edu/phylip/newick_doc.html
    
    vector<int> cl;
    vector<int>::const_iterator j;
    
    if (H < (int) data.size()) {
        // singleton label is the data row
        cl = data[ data.size() - 1 - H ];
    }
    else {
        // recurse then get canonical labeling of internal node
        printf("(");
        
        j = gooddecomp[H].begin();
        newicktraceback(data, gooddecomp, Sv, *j);
        
        for (++j; j != gooddecomp[H].end(); ++j) {
            printf(",");
            newicktraceback(data, gooddecomp, Sv, *j);
        }
        
        printf(")");
        
        // get cannonical labeling
        getCanonicalLabeling(gooddecomp[H], Sv, cl, H);
    }
    
    if (H < (int) gooddecomp.size() - 1) {
        // print label
        j = cl.begin();
        if (*j == -1) printf("'*");
        else          printf("'%d", *j);
        
        for (++j; j != cl.end(); ++j) {
            if (*j == -1) printf(" *");
            else          printf(" %d", *j);
        }
        
        printf("'");
    }
}

int treesizetraceback(const vector<vector<int> > &gooddecomp, int H) {
    int nodeCount = 1;
    for(unsigned int i = 0; i < gooddecomp[H].size(); i++) {
        if (gooddecomp[H][i] != H) {
            nodeCount += treesizetraceback(gooddecomp, gooddecomp[H][i]);
        }
    }
    return nodeCount;
}
