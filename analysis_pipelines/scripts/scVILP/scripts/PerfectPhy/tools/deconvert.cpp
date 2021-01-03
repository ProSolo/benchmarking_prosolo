#include<vector>
#include<fstream>
#include<iostream>
#include<cctype>
using namespace std;

enum ArgError { ERROR = -1, QUIT = 0, NONE = 1 };

ArgError handleArguments(int argc, char** argv, string &mapfile ) {
    const char* helpmessage = 
"Usage: ./deconvert.out [option] \n"
"\n"
"Options:\n"
"-h, --help       Displays this message then quits.\n"
"-mapfile FILE    Read the state mapping from FILE (REQUIRED).\n";
    
    string arg;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            cerr << "Bad argument: " << argv[i] << " is not a flag." << endl;
            return ERROR;
        }
        
        arg = argv[i];
        
        if (arg.compare("-mapfile") == 0) {
            if (argc > i + 1) {
                mapfile = argv[i+1];
                
                i++;
            }
            else {
                cerr << "Bad argument: " << arg << " must precede a filename." << endl;
                
                return ERROR;
            }
        }
        else if (arg.compare("-h") == 0 || arg.compare("--help") == 0) {
            cout << helpmessage;
            
            return QUIT;
        }
        else {
            cerr << "Bad argument: " << arg << " is not a valid flag." << endl;
            
            return ERROR;
        }
    }
    
    if (mapfile.empty()) {
        cerr << "Bad argument: " << "no mapfile given." << endl;
        
        return ERROR;
    }
    
    return NONE;
}

int main(int argc, char** argv) {
    char dummychar = '\0';
    string mapfile;
    
    #define NAMEWIDTH 10
    
    // improve io speed since not using C streams
    ios_base::sync_with_stdio(false);
    
    ArgError argerr = handleArguments(argc, argv, mapfile);
    if (argerr != NONE)
        return argerr;
    
    ifstream inf(mapfile.c_str());
    
    // read number of characters
    int index, m;
    inf >> index >> m;
    
    // check for a dummy character
    inf >> noskipws;
    dummychar = inf.peek();
    inf >> skipws;
    
    if (dummychar != '\n') {
        inf >> dummychar; // consume actual dummy char
        
        int temp; inf >> temp; // temp must be -1
    }
    else
        dummychar = '\0'; // reset
    
    
    // read in rows of mapfile
    vector<vector<char> > allmaps(m);
    
    for (int j = 0; j < m; j++) {
        vector<char>& row = allmaps[j];
        
        int statecount; inf >> statecount;
        
        row.resize(statecount);
        
        for (int k = 0; k < statecount; k++) {
            char& cell = row[k];
            inf >> cell;
        }
    }
    
    
    // traverse stdin as newick format tree and replace characters as needed
    
    cin >> noskipws;
    while(! cin.eof()) {
        // get and echo next character
        char nextc;
        cin >> nextc;
        cout << nextc;
        
        if (nextc == '\'') {
            
            // found label
            for (int j = 0; j < m; j++) {
                while (cin.peek() == ' ') {
                    // skip and echo spaces
                    cin >> nextc;
                    cout << nextc;
                }
                
                if (cin.peek() == '*') {
                    cin >> nextc;
                    cout << '*';
                }
                else {
                    int num;
                    cin >> num;
                    
                    if (num < (int) allmaps[j].size()) {
                        cout << allmaps[j][num];
                    }
                    else {
                        cerr << "Mapping Error: " << num << " is outside range [0," << allmaps[j].size() << ')' << endl;
                        return -1;
                    }
                }
                
            }
            
            // echo last '\''
            cin >> nextc;
            cout << nextc;
            continue;
        }
    }
    
    cin >> skipws;
    
    return 0;
}
