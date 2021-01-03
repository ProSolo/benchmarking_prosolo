#include<vector>
#include<iostream>
#include<sstream>
#include<cstdlib>
#include<ctime>
using namespace std;

#define ADDMISSING 1

enum ArgError { ERROR = -1, QUIT = 0, NONE = 1 };

ArgError handleArguments(int argc, char** argv, bool &fillmissing, double& prob) {
    const char* helpmessage =
"Usage: ./missingdata.out [option] \n"
"\n"
"Options:\n"
"-h, --help       Displays this message then quits.\n"
#if ADDMISSING
"-add PERCENT     Randomly add missing data with PERCENT probability per cell\n"
"                 then fills them in with new states\n"
#endif
"-fill            Replace instances of -1 with new states\n";
    
    if (argc == 1) {
        cerr << "This program cannot be executed with no arguments." << endl;
        return ERROR;
    }
    
    string arg;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            cerr << "Bad argument: " << argv[i] << " is not a flag" << endl;
            return ERROR;
        }
        
        arg = argv[i];
        
        if (arg.compare("-fill") == 0) {
            fillmissing = true;
        }
        #if ADDMISSING
        else if (arg.compare("-add") == 0) {
            if (argc > i + 1) {
                stringstream ss(argv[i+1]);
                ss >> prob;
                
                if ( !ss || prob < 0.0 || prob > 1.0) {
                    cerr << "Bad argument: " << argv[i+1] << " must be a valid double in range [0.0, 1.0]" << endl;
                    return ERROR;
                }
                
                i++;
            }
            else {
                cerr << "Bad argument: " << arg << " must precede a double" << endl;
                return ERROR;
            }
        }
        #endif
        else if (arg.compare("-h") == 0 || arg.compare("--help") == 0) {
            cout << helpmessage;
            
            return QUIT;
        }
        else {
            cerr << "Bad argument: " << arg << " is not a valid flag" << endl;
            
            return ERROR;
        }
    }
    
    return NONE;
}

#define readinto(var) do {\
        cin >> var;\
        if (cin.bad()) {\
            cerr << "Input error: " << #var << " could not be read" << endl;\
            cerr << "Expected input: index n m (n*m matrix cells...)" << endl;\
            exit(-1);\
        }\
    } while(false)

#define readintoas(type, var) do {\
        type temp;\
        readinto(temp);\
        var = temp;\
    } while(false)

int main(int argc, char** argv) {
    bool fillmissing = false;
    double prob = 500.0;
    
    // improve io speed since not using C streams
    ios_base::sync_with_stdio(false);
    
    ArgError argerr = handleArguments(argc, argv, fillmissing, prob);
    if (argerr != NONE)
        return argerr;
    
    // read size of matrix
    int index, n, m;
    readinto(index);
    readintoas(double, n);
    readintoas(double, m);
    
    vector<vector<int> > data(n, vector<int>(m));
    vector<int> statecount(m,0);
    
    for (int i = 0; i < n; i++) {
        vector<int>& row = data[i];
        
        for (int j = 0; j < m; j++) {
            int& cell = row[j];
            
            readinto(cell);
            
            if (cell < 0 && cell != -1) {
                cerr << "Input error: matrix values must >= -1 : " << cell << endl;
                cerr << "Expected input: index n m (n*m matrix cells...)" << endl;
                exit(-1);
            }
            
            if (statecount[j] <= cell)
                statecount[j] = cell+1;
        }
    }
    
    // add and fill in missing data
    #if ADDMISSING
    srand(time(NULL));
    #endif
    
    for (int i = 0; i < n; i++) {
        vector<int>& row = data[i];
        
        for (int j = 0; j < m; j++) {
            int& cell = row[j];
            
            if ( cell == -1
                #if ADDMISSING
                || ( prob > 0.0 && prob >= (rand() / (double) RAND_MAX) )
                #endif
                )
                cell = statecount[j]++;
        }
    }
    
    #if ADDMISSING
    // adjust states to eliminate gaps
    vector<int> stateMap;
    
    for (int j = 0; j < m; j++) {
        stateMap.assign(statecount[j], -1);
        int maxState = 0;
        
        for (int i = 0; i < n; i++) {
            int& cell = data[i][j];
            
            if (stateMap[cell] == -1)
                stateMap[cell] = maxState++;
            
            cell = stateMap[cell];
        }
        statecount[j] = maxState;
    }
    #endif
    
    // write data to stdout
    cout << index << endl << n << ' ' << m << endl;
    
    for (int i = 0; i < n; i++) {
        const vector<int>& row = data[i];
        
        for (int j = 0; j < m; j++)
            cout << row[j] << ' ';
        
        cout << endl;
    }
    
    return 0;
}

#undef readinto
#undef readintoas
