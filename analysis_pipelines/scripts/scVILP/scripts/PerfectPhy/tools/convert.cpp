#include<vector>
#include<fstream>
#include<iostream>
#include<cctype>
using namespace std;

enum ArgError { ERROR = -1, QUIT = 0, NONE = 1 };

ArgError handleArguments(int argc, char** argv, char &dummychar, string &index, string &mapfile ) {
    const char* helpmessage = 
"Usage: ./convert.out [option] \n"
"\n"
"Options:\n"
"-h, --help       Displays this message then quits.\n"
"-dummy CHAR      Set CHAR as a dummy state (outputs -1).\n"
"-index I         Set index of dataset as I (default is 1).\n"
"-mapfile FILE    Writes the state mapping to FILE\n";
    
    string arg;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            cerr << "Bad argument: " << argv[i] << " is not a flag" << endl;
            return ERROR;
        }
        
        arg = argv[i];
        
        if (arg.compare("-dummy") == 0) {
            if (argc > i + 1) {
                dummychar = argv[i+1][0];
                
                if (!isgraph(dummychar)) {
                    cerr << "Bad argument: dummy character must be graphical" << endl;
                    return ERROR;
                }
                
                i++;
            }
            else {
                cerr << "Bad argument: " << arg << " must precede a character" << endl;
                
                return ERROR;
            }
        }
        else if (arg.compare("-index") == 0) {
            if (argc > i + 1) {
                index = argv[i+1];
                
                for (unsigned int j = 0; j < index.size(); j++)
                    if ( !( isdigit(index[j]) || (j == 0 && index[j] == '-' && index.size() > 1) ) ) {
                        cerr << "Bad argument: index must be an integer, not " << index << endl;
                        return ERROR;
                    }
                
                i++;
            }
            else {
                cerr << "Bad argument: " << arg << " must precede an integer" << endl;
                
                return ERROR;
            }
        }
        else if (arg.compare("-mapfile") == 0) {
            if (argc > i + 1) {
                mapfile = argv[i+1];
                
                i++;
            }
            else {
                cerr << "Bad argument: " << arg << " must precede a filename" << endl;
                
                return ERROR;
            }
        }
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

class Mapping {
private:
    char c;
    int n;
    
public:
    Mapping(char in, int out):c(in), n(out) {}
    
    void printTo(ostream& outf) const { outf << c << ' ' << n << ' '; }
};

#define readinto(var) do {\
        cin >> var;\
        if (!cin) {\
            cerr << "Input error: " << #var << " could not be read" << endl;\
            cerr << "Expected input:" << endl << "n m" << endl << "(n rows of: (10-char name) (m matrix cells...) )" << endl;\
            return -1;\
        }\
    } while(false)

int main(int argc, char** argv) {
    char dummychar = '\0';
    string index;
    string mapfile;
    
    #define NAMEWIDTH 10
    
    // improve io speed since not using C streams
    ios_base::sync_with_stdio(false);
    
    ArgError argerr = handleArguments(argc, argv, dummychar, index, mapfile);
    if (argerr != NONE)
        return argerr;
    
    if (index.empty())
        index = "1";
    
    bool outputmap = (mapfile.size() > 0);
    
    // read size of matrix
    int n, m;
    readinto(n);
    readinto(m);
    
    vector<vector<int> > data(n, vector<int>(m));
    
    for (int i = 0; i < n; i++) {
        vector<int>& row = data[i];
        char cell;
        
        // skip over the 10-character name
        cin >> noskipws;
        while (cin.get() != '\n');
        cin.ignore(NAMEWIDTH);
        cin >> skipws;
        
        for (int j = 0; j < m; j++) {
            int& mcell = row[j];
            
            // skip over spaces between cells
            //while (cin.peek() == ' ')
            //    cin.get(c);
            
            readinto(cell);
            
            if (cell == dummychar)
                mcell = -1;
            else
                mcell = cell;
        }
    }
    
    // convert to multistate numbers without gaps
    vector<int> stateMap;
    vector<vector<Mapping> > allmaps;
    
    if (outputmap) {
        allmaps.resize(m);
        
        for(int j = 0; j < m; j++)
            allmaps[j].reserve(n);
    }
    
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            int& cell = data[i][j];
            
            if (cell == -1)
                continue;
            
            int k, statecount = stateMap.size();
            
            for (k = 0; k < statecount; k++)
                if (stateMap[k] == cell) {
                    cell = k;
                    break;
                }
            
            if (k == statecount) {
                if (outputmap)
                    allmaps[j].push_back(Mapping(cell, k));
                
                stateMap.push_back(cell);
                cell = k;
            }
        }
        
        stateMap.clear();
    }
    
    if (outputmap) {
        // write mapping to mapfile
        ofstream outf(mapfile.c_str());
        
        outf << index << endl << m;
        
        if (dummychar != '\0') {
            outf << ' ';
            Mapping(dummychar, -1).printTo(outf);
        }
        outf << endl;
        
        for (int j = 0; j < m; j++) {
            const vector<Mapping>& row = allmaps[j];
            
            outf << (int) row.size() << "   ";
            
            for (unsigned int k = 0; k < row.size(); k++)
                row[k].printTo(outf);
            
            outf << endl;
        }
    }
    
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
