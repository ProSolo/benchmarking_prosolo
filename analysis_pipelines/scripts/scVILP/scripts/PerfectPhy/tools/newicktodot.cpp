#include<vector>
#include<iostream>
#include<cstdlib>
using namespace std;


enum ArgError { ERROR = -1, QUIT = 0, NONE = 1 };

ArgError handleArguments(int argc, char** argv, string &name) {
    const char* helpmessage = "\
Usage: ./newicktodot.out [option] \n\
\n\
Options:\n\
-h, --help       Displays this message then quits.\n\
-name NAME       Gives the graph the name NAME (default is \"PerfectPhylogeny\")\n";
    
    string arg;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            cerr << "Bad argument: " << argv[i] << " is not a flag" << endl;
            return ERROR;
        }
        
        arg = argv[i];
        
        if (arg.compare("-name") == 0) {
            if (argc > i + 1) {
                name = argv[i+1];
                
                i++;
            }
            else {
                cerr << "Bad argument: " << arg << " must precede a string" << endl;
                
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

void writenode(int depth);
void writesubtree(int depth);
void writechildren(int depth, vector<int> &children);

int nodeid;
vector<int> support;

void setSupport(int i, int s) {
    if (i >= (int) support.size())
        support.resize(i+1);
    
    support[i] = s;
}

void writeEdges(int parent, vector<int>& children) {
    for(unsigned int i = 0; i < children.size(); i++) {
        cout << "    " << nodeid << " -- " << children[i];
        
        if (! support.empty())
            cout << " [ label = \"" << support[ children[i] ] << "\" ];" << endl;
        else
            cout << ';' << endl;
    }
}

void iogood() {
    if (! cin.good() ) {
        cerr << "IO error! cin.rdstate() == " << cin.rdstate() << endl;
        
        exit(1);
    }
}

void mustbe(char c) {
    iogood();
    
    if (cin.peek() != c) {
        cerr << "Format error! Next char must be '" << c << "' not '" << (char) cin.peek() << '\'' << endl;
        exit(1);
    }
    
    cin.get();
    
    iogood();
}

// write a node definition and label
void writenode(int depth) {
    mustbe('\'');
    
    cout << "    " << nodeid << " [ label = \"";
    
    cin >> noskipws;
    while(cin.peek() != '\'') {
        char c;
        cin.get(c);
        cout << c;
    }
    cin >> skipws;
    
    cout << "\" ];" << endl;
    
    mustbe('\'');
    
    if (cin.peek() == ':') {
        mustbe(':');
        int s;
        
        cin >> s;
        
        setSupport(nodeid, s);
    }
}

// writes a subtree
void writesubtree(int depth) {
    vector<int> children;
    
    if (cin.peek() == '(') { // internal node
        writechildren(depth+1, children);
        cout << "    {rank=" << depth << "; " << nodeid << ";}" << endl;
    }
    else // leaf
        cout << "    {rank=max; " << nodeid << ";}" << endl;
    
    writenode(depth); // current node
    
    writeEdges(nodeid, children);
}

// writes a list of subtrees
void writechildren(int depth, vector<int> &children) {
    mustbe('(');
    
    while (cin.peek() != ')') {
        writesubtree(depth);
        children.push_back(nodeid);
        
        nodeid++;
        
        if (cin.peek() == ',')
            cin.get();
    }
    
    mustbe(')');
}

int main(int argc, char** argv) {
    string name;
    
    // improve io speed since not using C streams
    ios_base::sync_with_stdio(false);
    
    ArgError argerr = handleArguments(argc, argv, name);
    if (argerr != NONE)
        return argerr;
    
    if (name.empty())
        name = "PerfectPhylogeny";
    
    cout << "graph " << name << " {" << endl;
    cout << "    " << "splines = false;" << endl;
    
    nodeid = 0;
    vector<int> children;
    
    writechildren(1, children);
    
    #if 1
        cout << "    {rank=0; " << nodeid << ";}" << endl;    
        cout << "    " << nodeid << " [ label = \"\" shape = \"point\" ];" << endl;
        
        writeEdges(nodeid, children);
    #else
        int parent = children.back();
        children.pop_back();
        
        writeEdges(parent, children);
    #endif
    
    cout << '}' << endl;
    
    return 0;
}
