#ifndef __PPDAG_H__
#define __PPDAG_H__

class PPDAG {
public:
    // forward declarations
    class ProductNode;
    class SumNode;
    class tree_iterator;
    
    // maxPC is the pseudo proper cluster which represents all taxa
    // rootTaxaPC is the taxon which is at the "top" of the DAG
    // - all Product Nodes under maxPC are rooted by it
    
    int maxPC, rootTaxaPC, n;
    vector<vector<SumNode> > sumNodes;  // for each proper cluster H, for each valid parent y of H
    vector<vector<int> > nodeLabels;    // stores all possible node labels
    
    tree_iterator tree_begin();
    tree_iterator tree_end();
    tree_iterator tree_rbegin();
    tree_iterator tree_rend();
    
    // constructors
    PPDAG(const vector<vector<int> > &propercluster, const Dict &pcdict);
    PPDAG(const char* filename);
protected:
    bool makeFromFile(ifstream &inf);
public:
    
    // construction functions (for sprop4.cpp)
    void resizeRow(int H, int s);
    void insertProductNode(ProductNode &pn, int H, int y);
    void insertManyProductNodes(vector<ProductNode> &mis, int H, int y);
    bool handleInterior(SumNode *sn, vector<int> &interiorTaxaVp);
    void finalize(const vector< vector<vector<int> > > &Conn, const vector< vector<vector<int> > > &Vp, const Dict &pcdict, const vector<vector<int> > &propercluster, const vector<vector<int> > &data);
    
    // printing functions
    void makeDOT() const;
    void makeDOTForTrees() const;
    void writeProductNode(ofstream &outf, const ProductNode &pn, int H, int y, bool notseen) const;
    void prettyprint(ostream &outf) const;
    
    // accessor functions
    int getCount();
    int getCount() const;
    int getVpSize(int H) const;
    const SumNode& getSumNode(int H, int y) const;
    const SumNode& getSumNode(const pair<int,int> Hy) const;
    const vector<int>& getNodeLabel(int label) const;
    vector<int> getNodeLabel(int label);
    
    bool isSingleton(int H) const;
    
    // dynamic programming iteration functions
    void bottomDP(int &H, int &y) const;
    void topDP(int &H, int &y) const;
    bool nextSumNode(int &H, int &y) const;
    bool prevSumNode(int &H, int &y) const;
    
    // algorithms
    int minNodeCount() const;
    int maxNodeCount() const;
    void getSupportStats(vector<vector<int> > &snSupport, vector<int> &splitSupport, vector<int> &nodeSupport) const;
    pair<int,int> getMaxSupportTree(vector<int> &support) const;
    bool getRobinsonFoulds(vector<vector<int> > &RF) const;
    bool printRobinsonFoulds() const;
    
    /*
    // memoization-supporting, pre-order recursive apply
    // return false in sfunc/pfunc to stop downwards calls
    template<class T>
    void recursePreorder(bool (*sfunc)(T &output, const PPDAG &dag, pair<int,int> snIndex, const SumNode &sn), 
                         bool (*pfunc)(T &output, const PPDAG &dag, pair<int,int> snIndex, const ProductNode &pn),
                         T &output, pair<int,int> root) const
    {
        int i, j;
        const SumNode &sn = getSumNode(arg_split(root));
        
        // call sfunc and stop if returns false
        if (sfunc)
            if (!sfunc(output, *this, root, sn));
                return;
        
        // call pfunc on each product node and recurse 
        for (i = 0; i < sn.children.size(); i++) {
            if (!pfunc || pfunc(output, *this, root, sn.children[i]))
                for (j = 0; j < sn.children[i].decomp.size(); j++)
                    recursePreorder(sfunc, pfunc, output, sn.children[i].decomp[j]);
        }
        
        if (i == 0 && pfunc) {
            // singleton
            ProductNode pn;
            pn.count = 1;
            pn.rootLabel = root.first;
            
            pfunc(output, *this, root, pn);
        }
    }
    
    // plain post-order recursion
    template<class T>
    void recursePostorder(void (*sfunc)(T &output, const PPDAG &dag, pair<int,int> snIndex, const SumNode &sn), 
                          void (*pfunc)(T &output, const PPDAG &dag, pair<int,int> snIndex, const ProductNode &pn),
                          T &output, pair<int,int> root) const
    {
        int i, j;
        const SumNode &sn = getSumNode(arg_split(root));
        
        for (i = 0; i < sn.children.size(); i++) {
            for (j = 0; j < sn.children[i].decomp.size(); j++)
                recursePostorder(sfunc, pfunc, output, sn.children[i].decomp[j]);
            
            if (pfunc)
                pfunc(output, *this, root, sn.children[i]);
        }
        
        if (i == 0 && pfunc) {
            // singleton
            ProductNode pn;
            pn.count = 1;
            pn.rootLabel = root.first;
            
            pfunc(output, *this, root, pn);
        }
        
        if (sfunc)
            sfunc(output, *this, root, sn);
    }
    */
};

class PPDAG::ProductNode {
public:
    int rootLabel, count;
    vector<pair<int,int> > decomp; // 2D indeces in dag.sumNodes
    
    ProductNode();
    ProductNode(const vector<pair<int,int> > &d);
    ProductNode(const ProductNode &pn);
    
    int  getCount(PPDAG &dag);
    void getIValues(const PPDAG &dag, vector<int> &iValues, int i) const;
    int  collapseIValues(const PPDAG &dag, const vector<int> &iValues) const;
    bool hasTaxaRoot(int n) const;
    
    bool operator<(const ProductNode& rhs) const;
    bool operator==(const ProductNode &rhs) const;
    bool operator!=(const ProductNode &rhs) const;
};


class PPDAG::SumNode {
public:
    int count;
    vector<ProductNode> children;
    
    SumNode();
    SumNode(vector<ProductNode> &cld);
    
    int getCount(PPDAG &dag);
};

class PPDAG::tree_iterator {    
public:
    typedef struct {
        int pc, vp, // indeces in dag->sumNodes
            product,// which product node was chosen
            label,  // rootLabel of the node: if != pc then an internal taxa
            isum,   // i value at the sum node
            iprod,  // i value at the product node
            prev;   // parent's index in stack
        vector<int> next; // childrens' indeces in stack (maps to pn->decomp, -1 means internal)
    } Branch;
    
    // Branch b = { pc, vp, product, label, isum, iprod, prev, next[] }
    
protected:
    const PPDAG *dag;       // the dag to which this iterator belongs
    vector<Branch> stack;   // complete representation of the tree pointed to
    enum {
        UNDEF,
        MIDDLE,
        REND,
        END
    } boundaryState;
    
    // internal utility functions
    void getMPP(int pos);
    
    bool getNextMPP(int pos);
    void get0thMPP(int pos);
    
    bool getPrevMPP(int pos);
    void getLastMPP(int pos);
    
    void printTreeDetailed(int pos, int depth) const;
    void printNewick(int pos) const;
    
public:
    
    tree_iterator();
    tree_iterator(const PPDAG &d);
    tree_iterator(const PPDAG &d, int i);
    
    int getIndex() const;
    bool gotoIndex(int i);
    
    bool isValid() const; // checks if iteration is possible
    bool atend() const;   // checks if cannot iterate forward
    bool atrend() const;  // checks if cannot iterate backward
    
    tree_iterator& operator=(const tree_iterator &rhs);
    tree_iterator& operator++();
    tree_iterator& operator--();
    
    int operator*() const; // same as getIndex()
    
    bool operator!=(const tree_iterator &rhs) const;
    bool operator==(const tree_iterator &rhs) const;
    
    void getTreeSplits(vector<int> &output) const;
    void printTree() const;
    void printStack() const;
    void printTreeDetailed() const;
    void printNewick(vector<int> * const support = NULL, int pos = 0) const;
    
    // data extraction functions
    const vector<Branch>& getRawTree() const;
    int getNodeCount() const;
    void makeTreeDOT(ofstream &outf) const;
    void makeTreeDOTLabeled(ofstream &outf, vector<int> *support = NULL) const;
    void getLabels(vector<vector<int> > &output) const;
    bool isMPP(const vector<vector<int> > &data, const vector<int> &statecount) const;
};

#endif
