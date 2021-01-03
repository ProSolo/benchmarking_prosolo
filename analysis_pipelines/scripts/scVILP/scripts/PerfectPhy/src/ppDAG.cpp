#include "ppDAG.h"

// allows readable 2D indexing with pairs
#define pair_index(array,p) ( (array)[ (p).first ][ (p).second ] )

//////// ProductNode ////////

PPDAG::ProductNode::ProductNode():rootLabel(-1), count(-1) {}
PPDAG::ProductNode::ProductNode(const vector<pair<int,int> > &d):rootLabel(-1), count(-1), decomp(d) {}
PPDAG::ProductNode::ProductNode(const ProductNode &pn):rootLabel(pn.rootLabel), count(pn.count), decomp(pn.decomp) {}

// Pn.count = product of Sn.count below it
int PPDAG::ProductNode::getCount(PPDAG &dag) {
    if (count == -1) {
        count = 1;
        
        for (unsigned int i = 0; i < decomp.size(); i++)
            count *= pair_index(dag.sumNodes, decomp[i]).getCount(dag);
    }
    
    return count;
}

// returns whether the label of the node is a taxon from the original dataset
bool PPDAG::ProductNode::hasTaxaRoot(int n) const {
    return (rootLabel < n);
}

// get values of i for each sum node in decomp (linear time)
// example: i =           i3
//                 + i2 * c3
//            + i1 * c2 * c3
//       + i0 * c1 * c2 * c3
void PPDAG::ProductNode::getIValues(const PPDAG &dag, vector<int> &iValues, int i) const {
    int product = count;
    
    // find i values to give to each sum node in decomp
    for (unsigned int j = 0; j < decomp.size(); j++) {
        product /= dag.getSumNode(decomp[j]).count;
        iValues[j] = i / product;
        i -= iValues[j] * product;
    }
}

// inverse operation of getIValues()
int PPDAG::ProductNode::collapseIValues(const PPDAG &dag, const vector<int> &iValues) const {
    int i = 0, product = count;
    
    for (unsigned int j = 0; j < decomp.size(); j++) {
        product /= dag.getSumNode(decomp[j]).count;
        i += iValues[j] * product;
    }
    
    return i;
}


bool PPDAG::ProductNode::operator<(const ProductNode& rhs) const {
    return (decomp < rhs.decomp);
}

bool PPDAG::ProductNode::operator==(const ProductNode &rhs) const {
    return (decomp == rhs.decomp && rootLabel == rhs.rootLabel);
}

bool PPDAG::ProductNode::operator!=(const ProductNode &rhs) const {
    return !(*this == rhs);
}

//////// SumNode ////////

PPDAG::SumNode::SumNode():count(-1) {}
PPDAG::SumNode::SumNode(vector<ProductNode> &cld): count(-1), children(cld) {}

// Sn.count = sum of Pn.count below it
int PPDAG::SumNode::getCount(PPDAG &dag) {
    if (count == -1) {
        if (children.empty())
            count = 1;
        else  {
            count = 0;
            for (unsigned int i = 0; i < children.size(); i++)
                count += children[i].getCount(dag);
        }
    }
    
    return count;
}



//////// PPDAG ////////

PPDAG::PPDAG(const vector<vector<int> > &propercluster, const Dict &pcdict) { 
    int root;
    
    n = propercluster[0].size();
    
    // the root taxon of the dag must be a proper cluster, searches for the highest numbered one
    for (root = n - 1; root >= 0; root--)
        if (-1 != pcdict.lookupComplement(propercluster[root]))
            break;
    
    rootTaxaPC = root;
    maxPC = (int) propercluster.size() - 1;
    sumNodes.resize(propercluster.size());
}

// inputs the DAG from a file created by prettyprint()
PPDAG::PPDAG(const char* filename):maxPC(-1), rootTaxaPC(-1), n(-1) {
    ifstream inf(filename);
    
    if (!inf)
        printf("DAG input file not found.\n");
    else if (!makeFromFile(inf))
        printf("DAG input file formatted incorrectly.\n");
}

// if internal taxa are allowed, changes product nodes to reflect this
bool PPDAG::handleInterior(SumNode *sn, vector<int> &interiorTaxaVp) {
    vector<int> interiors;
    
    const int snChildSize = sn->children.size();
    
    // loop over all product nodes present at the start
    for (int i = 0; i < snChildSize; i++) {
        ProductNode* pn = &sn->children[i];
        
        // find all possible interior taxa
        int decompsize = (int) pn->decomp.size();
        for (int j = 0; j < decompsize; j++) {
            pair<int,int> &G = pn->decomp[j];
            
            if (G.first < n) {
                if (interiorTaxaVp[G.first] == G.second) {
                    interiors.push_back(j);
                }
                else if (interiorTaxaVp[G.first] == -2) {
                    interiors.assign(1, j);
                    break;
                }
            }
        }
        
        if (interiors.empty())
            continue;
        
        // copy product node
        ProductNode origpn = *pn;
        
        // modify the original pn for the first root
        vector<int>::const_iterator j = interiors.begin();
        
        pn->rootLabel = pn->decomp[*j].first;
        pn->decomp.erase( pn->decomp.begin() + *j);
        
        // then add new Product Nodes for the remaining roots (if they exist)
        for (++j; j != interiors.end(); ++j) {
            sn->children.push_back(origpn);
            pn = &sn->children.back();
            
            pn->rootLabel = pn->decomp[*j].first;
            pn->decomp.erase( pn->decomp.begin() + *j);
        }
        
        cout << "---" << endl;
        
        interiors.clear();
    }
    
    // return whether any new nodes were added
    return (snChildSize == (int) sn->children.size());
}

// in constructDAG(), this is a required call to make finishing touches on the DAG to allow for use
void PPDAG::finalize(const vector< vector<vector<int> > > &Conn, const vector< vector<vector<int> > > &Vp, const Dict &pcdict, 
                     const vector<vector<int> > &propercluster, const vector<vector<int> > &data) {
    int H, y, newSize, properclustercount = propercluster.size();
    SumNode *sn;
    vector<ProductNode>::iterator pn;
    vector<pair<int,int> >::iterator G;
    vector<vector<bool> > inDAG(properclustercount);
    
    // the root taxa has no true sum node in the DAG, so remove it
    resizeRow(rootTaxaPC,0);
    
    vector<int> interiorTaxaVp(n, -1);
    // interiorTaxaVp holds info for singletons that can be contracted:
    // -1 == exterior, -2 == must always be interior, otherwise can be interior
    
    if (gbl::ALLOWINTERNAL) {
        for (H = 0; H < n; H++) {
            if (-1 == pcdict.lookupComplement(propercluster[H]) )
                interiorTaxaVp[H] = -2;
            else
                interiorTaxaVp[H] = binary_find(Vp[H], Conn[H][0]);
        }
    }
    
    H = rootTaxaPC;
    
    if (interiorTaxaVp[H] == -1) {
        for (y = 0; y < getVpSize(H); y++) {
            if (compatible(Vp[H][y], Conn[H][0])) {
                interiorTaxaVp[H] = -2;
                break;
            }
        }
    }
    
    // initialize array to store whether nodes are in the DAG or not (not valid after handleInterior calls)
    for (H = 0; H < maxPC; H++)
        inDAG[H].assign(sumNodes[H].size(), false);
    
    inDAG[H].assign(1, true);
    
    // populate said array
    for ( ; H >= 0; H--) {
        int maxy = -1;
        for (y = 0; y < getVpSize(H); y++) {
            sn = &(sumNodes[H][y]);
            
            // disconnect nodes not connected to the root
            if (!inDAG[H][y]) {
                sn->children.clear();
                continue;
            }
            maxy = y;
            
            sn = &(sumNodes[H][y]);
            
            if (H < n) {
                // if singleton, add empty product node
                sn->children.resize(1);
            }
            else {
                // Remove duplicate product nodes in each sum node
                sort(sn->children.begin(), sn->children.end());
                newSize = unique(sn->children.begin(), sn->children.end()) - sn->children.begin();
                
                sn->children.resize(newSize);
                
                // mark children as in the DAG
                for (pn = sn->children.begin(); pn != sn->children.end(); ++pn) {
                    for (G = pn->decomp.begin(); G != pn->decomp.end(); ++G) {
                        inDAG[G->first][G->second] = true;
                    }
                }
            }
        }
        
        // shorten row to remove sum nodes not actually in the DAG
        resizeRow(H, maxy + 1);
    }
    
    if (gbl::ALLOWINTERNAL) {
        // Create/modify product nodes if interior nodes exist
        // resort if any new product nodes were created
        
        for (H = 0; H < maxPC; H++) {
            for (y = 0; y < getVpSize(H); y++) {
                sn = &(sumNodes[H][y]);
                
                if (!inDAG[H][y])
                    continue;
                
                if (handleInterior(sn, interiorTaxaVp))
                    sort(sn->children.begin(), sn->children.end());
            }
        }
    }
    // contract rootTaxaPC as product node roots at the top
    sn = &(sumNodes[maxPC][0]);
    handleInterior(sn, interiorTaxaVp);
    
    
    // the dag is structurally finished, so find tree counts (used by DP functions)
    getCount();
    
    
    // next create a set of all node labels
    
    set<vector<int> > sortedLabels;
    
    bottomDP(H,y);
    
    while (nextSumNode(H,y)) {
        sn = &(sumNodes[H][y]);
        
        for (pn = sn->children.begin(); pn != sn->children.end(); ++pn) {
            if (pn->decomp.empty())
                continue;
            
            G = pn->decomp.begin();
            
            if (pn->rootLabel == -1)
                sortedLabels.insert(Vp[G->first][G->second]);
        }
    }
    
    // treat labels which are taxa as special
    // they are put in reverse-sorted order at the front of nodeLabels
    // so their indeces are equal to the taxon's proper cluster index
    for (H = n - 1; H >= 0; H--) {
        sortedLabels.erase(data[H]);
        nodeLabels.push_back(data[H]);
    }
    
    // push the other node labels in sorted order
    nodeLabels.insert(nodeLabels.end(), sortedLabels.begin(), sortedLabels.end());
    
    
    // then look at each product node to give it the correct index
    vector<vector<int> >::iterator start = nodeLabels.begin() + n;
    
    bottomDP(H,y);
    
    while (nextSumNode(H,y)) {
        sn = &(sumNodes[H][y]);
        
        for (pn = sn->children.begin(); pn != sn->children.end(); ++pn) {
            if (pn->decomp.empty()) {
                pn->rootLabel = H;
                continue;
            }
            
            G = pn->decomp.begin();
            
            if (pn->rootLabel == -1) {
                pn->rootLabel = lower_bound(start, nodeLabels.end(), Vp[G->first][G->second]) - nodeLabels.begin();
                
                if (pn->rootLabel == (int) nodeLabels.size() || nodeLabels[ pn->rootLabel ] != Vp[G->first][G->second]) {
                    pn->rootLabel = n - 1 - (lower_bound(data.begin(), data.end(), Vp[G->first][G->second]) - data.begin());
                }
            }
        }
    }
    
}


////////////////


void PPDAG::resizeRow(int H, int s) {
    if ((int) sumNodes[H].size() < s)
        sumNodes[H].resize(s);
}

void PPDAG::insertProductNode(ProductNode &pn, int H, int y) {
    sumNodes[H][y].children.push_back(pn);
}

// add a product node for every pn in mis to the sum node for (H,y)
void PPDAG::insertManyProductNodes(vector<ProductNode> &mis, int H, int y) {
    ProductNode pn;
    
    for (unsigned i = 0; i < mis.size(); i++) {
        pn = mis[i];
        
        insertProductNode(pn, H, y);
    }
}

const vector<int>& PPDAG::getNodeLabel(int label) const {
    return nodeLabels[label];
}

vector<int> PPDAG::getNodeLabel(int label) {
    return nodeLabels[label];
}

// calculate the count for every sum and product node, returning the total count
int PPDAG::getCount() {
    return sumNodes[maxPC][0].getCount(*this);
}

// simply returns the count, must be computed already
int PPDAG::getCount() const {
    return sumNodes[maxPC][0].count;
}

const PPDAG::SumNode& PPDAG::getSumNode(int H, int y) const {
    return sumNodes[H][y];
}

const PPDAG::SumNode& PPDAG::getSumNode(const pair<int,int> Hy) const {
    return sumNodes[Hy.first][Hy.second];
}

// returns whether H is a proper cluster of a single taxon
bool PPDAG::isSingleton(int H) const {
    return (H < n);
}

int PPDAG::getVpSize(int H) const {
    return sumNodes[H].size();
}


// dynamic programming incrementing functions
void PPDAG::bottomDP(int &H, int &y) const {
    H = 0;
    y = -1;
}

void PPDAG::topDP(int &H, int &y) const {
    H = maxPC;
    y = 1;
}

bool PPDAG::nextSumNode(int &H, int &y) const {
    for (y++; H <= maxPC; H++) {
        for ( ; y < getVpSize(H); y++)
            if (getSumNode(H,y).count > 0)
                return true;
        y = 0;
    }
    
    return false;
}

bool PPDAG::prevSumNode(int &H, int &y) const {
    for (y--; H >= 0; H--) {
        for ( ; y >= 0; y--)
            if (getSumNode(H,y).count > 0)
                return true;
        
        if (H != 0)
            y = getVpSize(H-1) - 1;
    }
    
    return false;
}


////////////////


// Writes a dot representation of a product node to outf
void PPDAG::writeProductNode(ofstream &outf, const ProductNode &pn, int H, int y, bool notseen) const {
    unsigned int i = 0;
    stringstream pnName (stringstream::in | stringstream::out); // the dot file "name" of the product node
    
    pnName << "\"p";
    
    if (! pn.decomp.empty()) {
        for (i = 0; i < pn.decomp.size(); i++)
             pnName << '.' << pn.decomp[i].first << '.' << pn.decomp[i].second;
    
        if (pn.hasTaxaRoot(n))
            pnName << '.' << pn.rootLabel;
    }
    else
        pnName << pn.rootLabel;
    
    pnName << '\"';
    
    if (notseen) {
        // print visible label for the product node
        outf << "  " << pnName.str() << " [ label = \"";
        
        if (! pn.decomp.empty()) {
            for (i = 0; i < pn.decomp.size(); i++)
                outf << "<f" << i << "> " << pn.decomp[i].first << " | ";
        }
        
        outf << "{ <f" << i << "> r= " << pn.rootLabel;
        
        i++;
        outf << " | <f" << i << "> c= " << pn.count << " } \"";
        
        if (pn.hasTaxaRoot(n))
            outf << " color=\"green\"";
        
        outf << " shape = \"record\" ];" << endl;
    }
    
    // create edge from (H,y) to the product node
    outf << "  " << H << '.' << y << " -> " << pnName.str() << ';' << endl;
    
    if (notseen) {
        // create edges from each part of the decomp to their sum nodes
        
        for (i = 0; i < pn.decomp.size(); i++)
            outf << "    " << pnName.str() << ":f" << i << " -> " << pn.decomp[i].first << '.' << pn.decomp[i].second << ';' << endl;
    }
}

// creates a dot representation of the DAG, writing to file in DAGDOTOUTF
void PPDAG::makeDOT() const {
    int H, y;
    
    stringstream filename;
    filename << gbl::DAGDOTOUTF << ".dot";
    
    ofstream outf(filename.str().c_str(), ios_base::out);
    
    // these are used to write duplicate Product Nodes only once
    set<ProductNode> uniquePN;
    vector<ProductNode>::const_iterator pn;
    pair<set<ProductNode>::iterator, bool> pos;
    
    outf << "digraph DAG { " << endl;
    
    for (H = 0; H < (int) sumNodes.size(); H++) {
        uniquePN.clear();
        
        for (y = 0; y < (int) sumNodes[H].size(); y++) {
            const SumNode& sn = getSumNode(H,y);
            
            if (sn.count == -1)
                continue;
            
            // print node label of sum node
            outf << H << "." << y << " [ label = \"" << H << ", " << y << " c=" << sn.count << "\" ];" << endl;
            
            for (pn = sn.children.begin(); pn != sn.children.end(); ++pn) {
                pos = uniquePN.insert(*pn);
                
                writeProductNode(outf, *pn, H, y, pos.second);
                outf << endl;
            }
        }
    }
    
    // now print out a table of node labels
    
    int id = 0;
    
    outf << "properclusters [ label = \"{ <f0> r";
    
    for (unsigned int i = 0; i < nodeLabels.size(); i++)
        outf << " | <f" << ++id << "> " << i;
    
    outf << "} | { <f" << ++id << "> Node Label";
    
    for (unsigned int i = 0; i < nodeLabels.size(); i++) {
        
        outf << " | <f" << ++id << '>';
        
        for (unsigned int j = 0; j < nodeLabels[i].size(); j++) {
            if (nodeLabels[i][j] == -1)
                outf << " *";
            else
                outf << ' ' << nodeLabels[i][j];
        }
    }
    outf << "} \" shape = \"record\" ];" << endl;
    
    outf << "} " << endl;
}

// creates a directory containing dot files for all trees in the dag
void PPDAG::makeDOTForTrees() const {
    tree_iterator mpp(*this);
    
    for (mpp.gotoIndex(0); !mpp.atend(); ++mpp) {
        stringstream filename;
        
        filename << gbl::DAGDOTOUTF << "_tree" << mpp.getIndex() << ".dot";
        
        ofstream outf(filename.str().c_str(), ios_base::out);
        
        mpp.makeTreeDOTLabeled(outf);
    }
}

// outputs text representation of the DAG which can be used to recreate it, outputs to outf
void PPDAG::prettyprint(ostream &outf) const {
    unsigned int i, j;
    int H, y;
    // vector<int> curDecomp; // unused
    const ProductNode *pn;
    
    outf << "R " << rootTaxaPC << ' ' << maxPC << ' ' << n << endl;
    
    topDP(H,y);
    
    while( prevSumNode(H,y) ) {
        const SumNode& sn = getSumNode(H,y);
        
        outf << "S " << H << ' ' << y << " c " << sn.count << " s " << sn.children.size() << endl;
        
        for (i = 0; i < sn.children.size(); i++) {
            pn = &( sn.children[i] );
            
            outf << "P c " << pn->count << " x " << pn->rootLabel << " s " << pn->decomp.size() << " D ";
            
            for (j = 0; j < pn->decomp.size(); j++)
                outf << pn->decomp[j].first << ' ' << pn->decomp[j].second << ' ';
            
            outf << endl;
        }
    }
    
    outf << "N " << nodeLabels.size() << ' ' << nodeLabels[0].size() << endl;
    
    for (i = 0; i < nodeLabels.size(); i++) {
        
        for (j = 0; j < nodeLabels[i].size(); j++)
            outf << nodeLabels[i][j] << ' ';
        
        outf << endl;
    }
}

// constructs DAG from file created by prettyprint() 
bool PPDAG::makeFromFile(ifstream &inf) {
    int i, j, H, y, length, width;
    char type;
    
    #define assertNextCharIs(c) do { \
        inf >> type;                 \
        if (type != c)               \
            return false;            \
    } while(false)
    
    sumNodes.clear();
    nodeLabels.clear();
    
    /*
        specification:
        
        R { rootTaxaPC } { maxPC } { n }
        S { H } { y } c { sn.count } p { sn.children.size }
        P c { pn.count } x { pn.rootLabel } s { pn.decomp.size } D { pn.decomp[0].first } { pn.decomp[0].second } { pn.decomp[1].first } ...
        ...
        N { n* } { m }
        nodeLabels[0][0] ...                   ...
        ...              ...                   ...
        ...              ... nodeLabels[n*-1][m-1]
    */
    
    assertNextCharIs('R');
    
    inf >> rootTaxaPC >> maxPC >> n;
    sumNodes.resize(maxPC + 1);
    H = -1;
    
    // read in Sum Nodes
    while (!inf.eof()) {
        
        inf >> type;
        if (type == 'N') break;
        else if (type != 'S') return false;
        
        i = H;
        inf >> H >> y;
        
        // resize row if new H found
        if (i != H)
            sumNodes[H].resize(y+1);
        
        SumNode& sn = sumNodes[H][y];
        
        assertNextCharIs('c');
        inf >> sn.count;
        
        assertNextCharIs('s');
        inf >> length;
        sn.children.resize(length);
        
        // read in Product Node children
        for (i = 0; i < length; i++) {
            ProductNode& pn = sn.children[i];
            
            assertNextCharIs('P');
            assertNextCharIs('c');
            inf >> pn.count;
            
            assertNextCharIs('x');
            inf >> pn.rootLabel;
            
            assertNextCharIs('s');
            inf >> width;
            pn.decomp.resize(width);
            
            assertNextCharIs('D');
            
            // read in decomposition
            for (j = 0; j < width; j++)
                inf >> pn.decomp[j].first >> pn.decomp[j].second;
        }
    }
    
    inf >> length >> width;
    nodeLabels.resize(length);
    
    for (i = 0; i < length; i++) {
        nodeLabels[i].resize(width);
        
        for (j = 0; j < width; j++)
            inf >> nodeLabels[i][j];
    }
    return true;
    
    #undef assertNextCharIs
}


// finds the minimum number of nodes required to build a tree in the DAG
// uses dynamic programming on min node counts of sub-DAGs, O(s+p)
int PPDAG::minNodeCount() const {
    int H, y, curMin, minNodes;
    vector<ProductNode>::const_iterator pn;
    
    vector<vector<int> > nodeCount(maxPC+1);
    
    for (H = 0; H <= maxPC; H++)
        nodeCount[H].assign(getVpSize(H), 0);
    
    bottomDP(H,y);
    
    while( nextSumNode(H,y) ) {
        minNodes = n * 2;
        
        const SumNode& sn = getSumNode(H,y);
        
        for (pn = sn.children.begin(); pn != sn.children.end(); ++pn) {
            curMin = 1; // starts at 1 for the root
            
            for (unsigned int i = 0; i < pn->decomp.size(); i++)
                curMin += pair_index(nodeCount, pn->decomp[i]);
            
            if (curMin < minNodes)
                minNodes = curMin;
        }
        
        nodeCount[H][y] = minNodes;
    }
    
    return nodeCount[maxPC][0];
}

// same as above, but maximum instead
int PPDAG::maxNodeCount() const {
    int H, y, curMax, maxNodes;
    vector<ProductNode>::const_iterator pn;
    
    vector<vector<int> > nodeCount(maxPC+1);
    
    for (H = 0; H <= maxPC; H++)
        nodeCount[H].assign(getVpSize(H), 0);
    
    bottomDP(H,y);
    
    while( nextSumNode(H,y) ) {
        maxNodes = 0;
        
        const SumNode& sn = getSumNode(H,y);
        
        for (pn = sn.children.begin(); pn != sn.children.end(); ++pn) {
            curMax = 1; // starts at 1 for the root
            
            for (unsigned int i = 0; i < pn->decomp.size(); i++)
                curMax += pair_index(nodeCount, pn->decomp[i]);
            
            if (curMax > maxNodes)
                maxNodes = curMax;
        }
        
        nodeCount[H][y] = maxNodes;
    }
    
    return nodeCount[maxPC][0];
}

// top-down dynamic program algorithm for finding support statistics, O(s+p)
void PPDAG::getSupportStats(vector<vector<int> > &snSupport, vector<int> &splitSupport, vector<int> &nodeSupport) const {
    int H, y, curSupport, supportHy;
    vector<ProductNode>::const_iterator pn;
    
    // set default values
    snSupport.resize(maxPC+1);
    for (H = 0; H < maxPC; H++)
        snSupport[H].assign(getVpSize(H), 0);
    
    splitSupport.assign(maxPC+1, 0);
    nodeSupport.assign(nodeLabels.size(), 0);
    
    // set starting dp values: the maxPC Sum Node is in every tree
    snSupport[maxPC].assign(1, getCount());
    nodeSupport[ rootTaxaPC ] = getCount();
    
    topDP(H,y);
    
    while( prevSumNode(H,y) ) {
        supportHy = snSupport[H][y];
        
        // add to the support of H the support of (H,y)
        splitSupport[H] += supportHy;
        
        const SumNode& sn = getSumNode(H,y);
        
        for (pn = sn.children.begin(); pn != sn.children.end(); ++pn) {
            curSupport = supportHy / sn.count * pn->count;
            
            // add to the support for the children sum nodes
            for (unsigned int i = 0; i < pn->decomp.size(); i++)
                pair_index(snSupport, pn->decomp[i]) += curSupport;
            
            // add to the support of the root of pn (if taxa are not allowed to be internal, must ignore coincidental internal taxa)
            if (gbl::ALLOWINTERNAL || H < n || pn->rootLabel >= n)
                nodeSupport[ pn->rootLabel ] += curSupport;
        }
    }
    
    // root taxon will be counted once for each product node 
    nodeSupport[rootTaxaPC] = getCount();
}

// dynamic programming algorithm to compute the mpp with the maximum split support, O(s+p)
// returns < support, tree index >
pair<int,int> PPDAG::getMaxSupportTree(vector<int> &support) const {
    int i, H, y, sum, product;
    pair<int,int> maxSubtree, curSubtree;
    
    vector<vector<pair<int,int> > > maxNum(maxPC+1); // map from sn to (support, subtree index)
    vector<ProductNode>::const_iterator pn;
    
    for (H = 0; H <= maxPC; H++)
        maxNum[H].resize(getVpSize(H));
    
    bottomDP(H,y);
    
    while( nextSumNode(H,y) ) {
        /*
        if (isSingleton(H)) {
            curSubtree.first = support[H];
            curSubtree.second = 0;
            
            maxNum[H][y] = curSubtree;
            continue;
        }
        */
        
        const SumNode& sn = getSumNode(H,y);
        
        maxSubtree.first = 0;
        sum = 0;
        
        // for each product node, find its maximum support subtree
        for (pn = sn.children.begin(); pn != sn.children.end(); ++pn) {
            
            // calculate the index and support of this product node's subtree
            curSubtree.first = 0;
            curSubtree.second = 0;
            product = 1;
            
            for (i = (int) pn->decomp.size() - 1; i >= 0; i--) {
                curSubtree.first  += pair_index(maxNum, pn->decomp[i]).first;
                curSubtree.second += pair_index(maxNum, pn->decomp[i]).second * product;
                
                product *= getSumNode(pn->decomp[i]).count;
            }
            
            // then adjust index based on the ordering of product nodes as children of sn
            curSubtree.second += sum;
            sum += pn->count;
            
            // lastly check if it has greater support or not
            if (curSubtree.first > maxSubtree.first)
                maxSubtree = curSubtree;
        }
        
        maxSubtree.first += support[H];
        
        maxNum[H][y] = maxSubtree;
    }
    
    return maxNum[maxPC][0];
}

// find all pairwise Robinson-Foulds distances in O(n t^2) = O(s t)
// the Robinson-Foulds Distance is the average number of exclusive splits between two trees
bool PPDAG::getRobinsonFoulds(vector<vector<int> > &RF) const {
    int i, j, t = getCount(), rfd, size;
    
    // RF will be a lower-triangular matrix of values
    RF.resize(t);
    for (j = 0; j < t; j++)
        RF[j].assign(j,-1);
    
    // has[H] = latest baseMpp which has H as a split, O(nt) total size
    vector<int> has(sumNodes.size(), -1);
    PPDAG::tree_iterator baseMpp(*this), mpp(*this);
    
    // for every tree j, O(t)
    for (baseMpp.gotoIndex(0); !baseMpp.atend(); ++baseMpp) {
        j = *baseMpp;
        const vector<tree_iterator::Branch>& basetree = baseMpp.getRawTree();
        size = basetree.size();
        
        // for every split in j, mark that j has it, O(n)
        for (i = 0; i < size; i++)
            has[ basetree[i].pc ] = j;
        
        // for every tree indexed less than j, O(t)
        for (mpp.gotoIndex(0); mpp != baseMpp; ++mpp) {
            const vector<tree_iterator::Branch>& rawtree = mpp.getRawTree();
            size = rawtree.size();
            rfd = baseMpp.getRawTree().size();
            
            // for each of its splits, O(n)
            for (i = 0; i < size; i++) {
                if (has[ rawtree[i].pc ] == j)
                    rfd--; // this split is shared, thus not exclusive
                else
                    rfd++; // this split is a new exclusive split
            }
            
            RF[j][*mpp] = rfd / 2;
            
            if (rfd < 2) {
                printf("0 Two trees are equivalent! #%d == #%d\n", j, *mpp);
                
                baseMpp.printTreeDetailed();
                mpp.printTreeDetailed();
                
                return false;
            }
        }
    }
    
    printf("1 All trees are distict.\n");
    
    return true;
}

// same as above, except the RF matrix is printed out as it is computed
bool PPDAG::printRobinsonFoulds() const {
    int i, j, rfd, size;
    
    // RF will be a lower-triangular matrix of values
    
    // has[H] = latest baseMpp which has H as a split, O(nt) total size
    vector<int> has(sumNodes.size(), -1);
    PPDAG::tree_iterator baseMpp(*this), mpp(*this);
    
    // for every tree j, O(t)
    for (baseMpp.gotoIndex(0); !baseMpp.atend(); ++baseMpp) {
        j = *baseMpp;
        const vector<tree_iterator::Branch>& basetree = baseMpp.getRawTree();
        size = basetree.size();
        
        // for every split in j, mark that j has it, O(n)
        for (i = 0; i < size; i++)
            has[ basetree[i].pc ] = j;
        
        printf("RF[%d] = { ", j);
        
        // all splits are defaulted to exclusive to j, O(t)
        //RF[j].assign(j, baseMpp.getRawTree().size());
        
        // for every tree indexed less than j, O(t)
        for (mpp.gotoIndex(0); mpp != baseMpp; ++mpp) {
            const vector<tree_iterator::Branch>& rawtree = mpp.getRawTree();
            size = rawtree.size();
            rfd = baseMpp.getRawTree().size();
            
            // for each of its splits, O(n)
            for (i = 0; i < size; i++) {
                if (has[ rawtree[i].pc ] == j)
                    rfd--; // this split is shared, thus not exclusive
                else
                    rfd++; // this split is a new exclusive split
            }
            
            printf("%d ", rfd / 2);
            
            if (rfd < 2) {
                printf("0 Two trees are equivalent! #%d == #%d\n", j, *mpp);
                
                baseMpp.printTreeDetailed();
                mpp.printTreeDetailed();
                
                return false;
            }
        }
        
        printf("}\n");
    }
    
    printf("1 All trees are distict.\n");
    
    return true;
}


//////// PPDAG::tree_iterator ////////

PPDAG::tree_iterator PPDAG::tree_begin() {
    PPDAG::tree_iterator tmpp(*this, 0);
    return tmpp;
}

PPDAG::tree_iterator PPDAG::tree_end() {
    PPDAG::tree_iterator tmpp(*this, -1);
    return tmpp;
}

PPDAG::tree_iterator PPDAG::tree_rbegin() {
    PPDAG::tree_iterator tmpp(*this, getCount() - 1);
    return tmpp;
}

PPDAG::tree_iterator PPDAG::tree_rend() {
    PPDAG::tree_iterator tmpp(*this, getCount());
    return tmpp;
}


PPDAG::tree_iterator::tree_iterator():dag(NULL),boundaryState(PPDAG::tree_iterator::UNDEF) {}

PPDAG::tree_iterator::tree_iterator(const PPDAG &d):dag(&d),boundaryState(PPDAG::tree_iterator::UNDEF) {
    stack.reserve(d.n * 2);
}

PPDAG::tree_iterator::tree_iterator(const PPDAG &d, int i):dag(&d) {
    stack.reserve(d.n * 2);
    
    gotoIndex(i);
    
    if (i < 0)
        boundaryState = PPDAG::tree_iterator::REND;
    else if (i >= d.getCount())
        boundaryState = PPDAG::tree_iterator::END;
    else
        boundaryState = PPDAG::tree_iterator::MIDDLE;
}

// returns enumerated index of the current tree
int PPDAG::tree_iterator::getIndex() const {
    if (stack.empty())
        return -1;
    else
        return stack[0].isum;
}

// returns whether the iterator is pointing to a tree or an endpoint, or is undefined ("singular")
bool PPDAG::tree_iterator::isValid() const {
    return (boundaryState != PPDAG::tree_iterator::UNDEF);
}

// returns whether the iterator is past the end of the trees
bool PPDAG::tree_iterator::atend() const {
    return (boundaryState == PPDAG::tree_iterator::END);
}

// returns whether the iterator is before the start of the trees
bool PPDAG::tree_iterator::atrend() const {
    return (boundaryState == PPDAG::tree_iterator::REND);
}

// returns immutable reference to the raw data
const vector<PPDAG::tree_iterator::Branch>& PPDAG::tree_iterator::getRawTree() const {
    return stack;
}

PPDAG::tree_iterator& PPDAG::tree_iterator::operator=(const PPDAG::tree_iterator &rhs) {
    if (*this != rhs) {
        dag = rhs.dag;
        stack = rhs.stack;
        boundaryState = rhs.boundaryState;
    }
    
    return *this;
}


// set iterator to a specific MPP in the DAG
// if outside the range, sets it to the end or rend
// returns whether the result points to a tree or not
bool PPDAG::tree_iterator::gotoIndex(int i) {
    if (i == getIndex())
        return true;
    
    if (i < 0) {
        boundaryState = PPDAG::tree_iterator::REND;
        return false;
    }
    
    if (i >= dag->getCount()) {
        boundaryState = PPDAG::tree_iterator::END;
        return false;
    }
    
    stack.clear();
    
    Branch b;
    
    b.pc = dag->maxPC;
    b.vp = 0;
    b.product = -1;
    b.label = -1;
    b.isum = i;
    b.iprod = -1;
    b.prev = -1;
    
    stack.push_back(b);
    
    if (i == 0)
        get0thMPP(0);     // access the 0th tree, O(n)
    
    else if (i == dag->getCount() - 1)             
        getLastMPP(0);    // access the last tree, O(n)   

    else
        getMPP(0);        // Random access tree i, O(n+p)
    
    boundaryState = PPDAG::tree_iterator::MIDDLE;
    
    return true;
}


// move to the next MPP in the DAG, O(n)
PPDAG::tree_iterator& PPDAG::tree_iterator::operator++() {
    if (boundaryState == PPDAG::tree_iterator::MIDDLE) {
        getNextMPP(0);
        
        if (getIndex() == -1)
            boundaryState = PPDAG::tree_iterator::END;
    }
    else if (boundaryState == PPDAG::tree_iterator::REND) {
        gotoIndex(0);
        
        boundaryState = PPDAG::tree_iterator::MIDDLE;
    }
    
    return *this;
}

// move to the previous MPP in the DAG, O(n)
PPDAG::tree_iterator& PPDAG::tree_iterator::operator--() {
    if (boundaryState == PPDAG::tree_iterator::MIDDLE) {
        getPrevMPP(0);
        
        if (getIndex() == -1)
            boundaryState = PPDAG::tree_iterator::REND;
    }
    else if (boundaryState == PPDAG::tree_iterator::REND) {
        gotoIndex(dag->getCount() - 1);
        
        boundaryState = PPDAG::tree_iterator::MIDDLE;
    }
    
    return *this;
}


bool PPDAG::tree_iterator::operator!=(const PPDAG::tree_iterator &rhs) const {
    return (   dag != rhs.dag
            || getIndex() != rhs.getIndex()
            || boundaryState != rhs.boundaryState );
}


bool PPDAG::tree_iterator::operator==(const PPDAG::tree_iterator &rhs) const {
    return !(operator!=(rhs));
}


int PPDAG::tree_iterator::operator*() const {
    return getIndex();
}

// fills output with all pc splits in the tree, specifically the side NOT containing the rootTaxa
void PPDAG::tree_iterator::getTreeSplits(vector<int> &output) const {
    output.clear();
    output.reserve(stack.size());
    
    for (unsigned int j = 0; j < stack.size(); j++)
        output.push_back(stack[j].pc);
}

// prints the rootLabels of all product nodes in the tree
void PPDAG::tree_iterator::printTree() const {
    printf("mpp[%d] = ", **this);
    for (vector<Branch>::const_iterator b = stack.begin(); b != stack.end(); ++b)
        printf("%d ", b->label);
    
    printf("\n");
}

// prints each Branch in detail
void PPDAG::tree_iterator::printStack() const {
    printf("mpp[%d] = \n", **this);
    
    for (vector<Branch>::const_iterator b = stack.begin(); b != stack.end(); ++b)
        printf("[%3d] = { pc = %3d, vp = %3d, product = %3d, label = %3d, isum = %3d, iprod = %3d, prev = %3d, |next| = %3d }\n", (int) (b - stack.begin()), 
                       b->pc,    b->vp,    b->product,    b->label,    b->isum,    b->iprod,    b->prev, (int) b->next.size());
}

// prints the node labels of the current tree hierarchically & recursively
void PPDAG::tree_iterator::printTreeDetailed() const {
    printTreeDetailed(0,0);
}
void PPDAG::tree_iterator::printTreeDetailed(int pos, int depth) const {
    if (pos == 0)
        printf("Tree #%d:\n", getIndex());
    
    printf("%*c%3d", depth * 3, ' ', stack[pos].pc);
    prettyprintarray(dag->getNodeLabel(stack[pos].label), "", PPAF_DUMMY);
    
    for (unsigned int child = 0; child < stack[pos].next.size(); child++)
        printTreeDetailed(stack[pos].next[child], depth + 1);
}

// prints the current tree in newick format
void PPDAG::tree_iterator::printNewick(vector<int> * const support, int pos) const {
    bool internal = stack[pos].next.size() > 0;
    
    // print children
    if (internal) {
        printf("(");
        printNewick(support, stack[pos].next[0]);
    }
    
    for (unsigned int child = 1; child < stack[pos].next.size(); child++) {
        printf(",");
        printNewick(support, stack[pos].next[child]);
    }
    
    if (internal) {
        printf( (pos != 0) ? ")" : ",");
    }
    
    // print this node's label
    const vector<int>& label = dag->getNodeLabel(stack[pos].label);
    
    printf("\'%d", label[0]);
    for (unsigned int j = 1; j < label.size(); j++)
        printf(" %d", label[j]);
    printf("\'");
    
    if (support) {
        // adds support in the "branch length" section of newick standard
        printf(":%d", (*support)[stack[pos].pc] );
    }
    
    if (internal && pos == 0)
        printf(")");
}



// returns the number of [product] nodes in the current tree
int PPDAG::tree_iterator::getNodeCount() const {    
    return 1 + stack.size();
}

////// Private tree_iterator Methods //////

// append the ith MPP to the stack and recurse (pre-order traversal)
// from root node, O(n+p) for random access, O(n) to get the zeroth

void PPDAG::tree_iterator::getMPP(int pos) {
    unsigned int j;
    Branch& b = stack[pos];
    const SumNode& sn = dag->getSumNode(b.pc,b.vp);
    
    // incoming correct definitions for b
    // def = { pc, vp, isum, prev }
    // undef = { product, label, iprod, next[] }
    
    /*
    if (b.pc < n) {
        // external singleton
        b.label = b.pc;
        return;
    }
    */
    
    if (b.product == -1) {
        // choose which product node to follow - O(p) worst, O(1) zeroth
        int sum = 0;
        for (j = 0; j < sn.children.size(); j++) {
            sum += sn.children[j].count;
            if (sum > b.isum)
                break;
        }
        b.product = j;
        
        b.iprod = b.isum - (sum - sn.children[j].count);
    }
    
    const ProductNode& pn = sn.children[b.product];
    b.label = pn.rootLabel;
    
    vector<int> iValues(pn.decomp.size(), 0);
    
    // find i values to give each sum node in the decomposition - O(n)
    pn.getIValues(*(dag), iValues, b.iprod);
    
    b.next.assign(pn.decomp.size(), -1);
    
    Branch b2;
    b2.product = -1;
    b2.label = -1;
    b2.iprod = -1;
    b2.prev = pos;
    
    // recurse on subtree sum nodes
    for (j = 0; j < pn.decomp.size(); j++) {
        // create a branch for this sum node and recurse - O(n) iterations total
        
        b2.pc = pn.decomp[j].first;
        b2.vp = pn.decomp[j].second;
        b2.isum = iValues[j];
        
        int stackEnd = stack.size();
        b.next[j] = stackEnd;
        
        stack.push_back(b2);
        
        if (b2.isum == 0)
            get0thMPP(stackEnd);
        
        //else if (b2.isum == sn.count - 1)
        //    getLastMPP(stackEnd);
        
        else
            getMPP(stackEnd);
    }
}


// increments to the next MPP in the subtree at stack[pos] (post-order traversal)
// returns false if reaches the end (subtree will be removed from stack)
// from root node, O(n)

bool PPDAG::tree_iterator::getNextMPP(int pos) {
    Branch& b = stack[pos];
    const SumNode& sn = dag->getSumNode(b.pc,b.vp);
    
    b.isum++;
    b.iprod++;
    
    if (b.isum == sn.count) {
        // sum node has overflowed, no more MPP's in this subtree
        stack.resize(pos);
        return false;
    }
    
    const ProductNode& pn = sn.children[b.product];
    
    if (b.iprod == pn.count) {
        // product node has overflowed
        // start getMPP on the next product node to get a new subtree
        
        b.product++;
        b.label = -1;
        b.iprod = 0;
        b.next.clear();
        
        stack.resize(pos+1);
        
        getMPP(pos);
        
        return true;
    }
    else { // we're staying in the current product node
        
        unsigned int j, bNextSize = b.next.size();
        
        // attempt to get the next MPP of each subtree until one is found - O(n) iterations
        for (j = bNextSize - 1; j >= 0; j--) {
            if (getNextMPP(b.next[j])) {
                j++; // the next subtree was found
                break;
            }
            else
                b.next[j] = -1; // since that subtree is finished, remove from stack
        }
        
        if (j == bNextSize)
            return true; // least significant subtree didn't end, so nothing more to check
        
        // create 0th subtrees for each overflowed subtree - from root node, O(n) iterations
        Branch b2;
        b2.product = 0;
        b2.label = -1;
        b2.isum = 0;
        b2.iprod = 0;
        b2.prev = pos;
        
        for ( ; j < pn.decomp.size(); j++) {
            
            b2.pc = pn.decomp[j].first;
            b2.vp = pn.decomp[j].second;
            
            int stackEnd = stack.size();
            b.next[j] = stackEnd;
            
            stack.push_back(b2);
            
            get0thMPP(stackEnd);
        }
        
        return true;
    }
    
}

// append the 0th MPP to the stack and recurse (pre-order traversal)
// from root node, O(n) specialized for 0th trees
void PPDAG::tree_iterator::get0thMPP(int pos) {
    unsigned int j;
    Branch& b = stack[pos];
    const SumNode& sn = dag->getSumNode(b.pc,b.vp);
    
    // incoming correct definitions for b
    // def = { pc, vp, isum, prev }
    // undef = { product, label, iprod, next[] }
    
    /*
    if (sn.children.empty()) {
        // external singleton
        b.label = b.pc;
        return;
    }
    */
    
    const ProductNode& pn = sn.children[0];
    
    b.product = 0;
    b.iprod = 0;
    b.label = pn.rootLabel;
    b.next.assign(pn.decomp.size(), -1);
    
    Branch b2;
    b2.isum = 0;
    b2.prev = pos;
    
    // recurse on subtree sum nodes
    for (j = 0; j < pn.decomp.size(); j++) {
        // create a branch for this sum node and recurse - O(n) iterations total
        
        b2.pc = pn.decomp[j].first;
        b2.vp = pn.decomp[j].second;
        
        int stackEnd = stack.size();
        b.next[j] = stackEnd;
        
        stack.push_back(b2);
        
        get0thMPP(stackEnd);
    }
}



// increments to the previous MPP in the subtree at stack[pos] (post-order traversal)
// returns false if reaches the end (subtree will be removed from stack)
// from root node, O(n)

bool PPDAG::tree_iterator::getPrevMPP(int pos) {
    Branch& b = stack[pos];
    
    b.isum--;
    b.iprod--;
    
    if (b.isum == -1) {
        // sum node has underflowed, no more MPP's in this subtree
        stack.resize(pos);
        return false;
    }
    
    if (b.iprod == -1) {
        // product node has underflowed
        // start getMPP on the previous product node to get a new subtree
        
        const SumNode& sn = dag->getSumNode(b.pc,b.vp);
        
        b.product--;
        b.label = -1;
        b.iprod = sn.children[b.product].count - 1;
        b.next.clear();
        
        stack.resize(pos+1);
        
        getMPP(pos);
        
        return true;
    }
    
    else { // we're staying in the current product node
        
        unsigned int j, bNextSize = b.next.size();
        
        // attempt to get the next MPP of each subtree until one is found - O(n) iterations
        for (j = bNextSize - 1; j >= 0; j--) {
            if (getPrevMPP(b.next[j])) {
                // true means the next subtree was found
                j++;
                break;
            }
            else {
                // since that subtree is finished, remove from stack
                b.next.pop_back();
            }
        }
        
        if (j == bNextSize)
            return true; // least significant subtree didn't end, so nothing more to check
        
        // create (count-1)-th subtrees for each overflowed subtree - O(n) iterations
        
        Branch b2;
        b2.label = -1;
        b2.prev = pos;
        
        const ProductNode& pn = dag->getSumNode(b.pc,b.vp).children[b.product];
        
        for ( ; j < pn.decomp.size(); j++) {
            
            b2.pc = pn.decomp[j].first;
            b2.vp = pn.decomp[j].second;
            
            int stackEnd = stack.size();
            b.next.push_back(stackEnd);
            
            stack.push_back(b2);
            
            getLastMPP(stackEnd);
        }
        
        return true;
    }
}

// append the last MPP to the stack and recurse (pre-order traversal)
// from root node, O(n) specialized for the highest-indexed tree
void PPDAG::tree_iterator::getLastMPP(int pos) {
    unsigned int j;
    Branch& b = stack[pos];
    const SumNode& sn = dag->getSumNode(b.pc,b.vp);
    
    // incoming definitions for b
    // def = { pc, vp, prev }
    // undef = { product, label, isum, iprod, next[] }
    
    /*
    if (sn.children.empty()) {
        // external singleton
        b.product = -1;
        b.label = b.pc;
        b.isum = 0;
        b.iprod = 0;
        return;
    }
    */
    
    const ProductNode& pn = sn.children.back();
    
    b.product = (int) sn.children.size() - 1;
    b.label = pn.rootLabel;
    b.isum = sn.count - 1;
    b.iprod = pn.count - 1;
    
    b.next.assign(pn.decomp.size(), -1);
    
    Branch b2;
    b2.prev = pos;
    
    // create a branch for this sum node and recurse - O(n) iterations total
    for (j = 0; j < pn.decomp.size(); j++) {
        
        b2.pc = pn.decomp[j].first;
        b2.vp = pn.decomp[j].second;
        
        int stackEnd = stack.size();
        b.next[j] = stackEnd;
        
        stack.push_back(b2);
        
        getLastMPP(stackEnd);
    }
}

////////////////////////////


// creates a dot file for the current MPP tree, labeled by proper cluster indeces and taxa (green)
void PPDAG::tree_iterator::makeTreeDOT(ofstream& outf) const {
    unsigned int j, stackEnd = stack.size();
    int n = dag->n;
    
    outf << "graph MPP" << getIndex() << " {" << endl;
    
    for (j = 0; j < stackEnd; j++) {
        const Branch& b = stack[j];
        
        // print node for b
        outf << b.pc << " [ label = \"" << b.label;
        
        if (b.label < n)
            outf << "\" color=\"green";
        
        outf << "\" ];" << endl;
        
        // print upward edge
        if (b.prev != -1)
            outf << stack[ b.prev ].pc << " -- " << b.pc << ';' << endl;
    }
    
    outf << "}" << endl;
}

// fills output with the node labels in the current perfect phylogeny
void PPDAG::tree_iterator::getLabels(vector<vector<int> > &output) const {
    vector<Branch>::const_iterator b;
    
    output.reserve(stack.size());
    
    for (b = stack.begin(); b != stack.end(); ++b)
        output.push_back( dag->getNodeLabel(b->label) );
}

// creates a dot file for the current MPP tree with character states labeled on nodes
// also prints the support of edges if given
void PPDAG::tree_iterator::makeTreeDOTLabeled(ofstream& outf, vector<int> *support) const {
    unsigned int i, j, mppCount = dag->getCount(), stackEnd = stack.size();
    int n = dag->n;
    vector<vector<int> > labels;
    
    getLabels(labels);
    
    outf << "graph MPP" << getIndex() << " {" << endl;
    
    for (j = 0; j < stackEnd; j++) {
        const Branch& b = stack[j];
        
        // print node for b
        outf << b.pc << " [ label = \"";
        
        for (i = 0; i < labels[j].size(); i++) {
            if (labels[j][i] == -1)
                outf << "* ";
            else
                outf << labels[j][i] << ' ';
        }
        
        if (b.label < n)
            outf << "\" color=\"green";
        
        outf << "\" ];" << endl;
        
        // print upward edge
        if (b.prev != -1) {
            outf << stack[ b.prev ].pc << " -- " << b.pc;
            if (support)
                outf << " [ label = \"" << (*support)[b.pc] * 100 / mppCount << "%\" ]";
            
            outf << ';' << endl;
         }
    }
    
    outf << "}" << endl;
}

// verify that the current tree is a minimal perfect phylogeny in O(mk n log n)
bool PPDAG::tree_iterator::isMPP(const vector<vector<int> > &data, const vector<int> &statecount) const {
    int pos, curChar, curState, stackEnd = stack.size();
    bool provenFalse = false;
    vector<vector<int> > labels;
    vector<int> subtree(stackEnd, 0);
    
    getLabels(labels);
    
    for (curChar = 0; curChar < (int) statecount.size(); curChar++) {
        for (curState = 0; curState < statecount[curChar]; curState++) {
            
            pos = 0;
            
            // first find the root of the subtree with curChar on curState
            for (pos++; pos < stackEnd; pos++) {
                if (labels[pos][curChar] == curState) {
                    subtree[pos] = 1;
                    break;
                }
                else
                    subtree[pos] = 0;
            }
            
            // then make sure it is connected
            for (pos++; pos < stackEnd; pos++) {
                if (labels[pos][curChar] == curState) {
                    subtree[pos] = 1;
                    
                    if (subtree[ stack[pos].prev ] == 0) {
                        provenFalse = true;
                        break;
                    }
                }
                else
                    subtree[pos] = 0;
            }
            
            if (provenFalse)
                break;
        }
        
        if (provenFalse)
            break;
    }
    
    if (provenFalse) {
        printf("mpp is NOT a perfect phylogeny: state %d on char %d\n", curState, curChar);
        return false;
    }
    
    // all state-induced subtrees are connected, now check if minimal
    
    vector<int> hasTaxon(dag->n, 0);
    
    if (stack[0].label < dag->n)
        hasTaxon[stack[0].label] = 1;
    
    for (pos = 1; pos < stackEnd; pos++) {
        // not minimal iff compatible nodes, unless it's a leaf taxon not allowed to be internal (including rootTaxaPC)
        
        if (compatible(labels[pos], labels[ stack[pos].prev ]) && (gbl::ALLOWINTERNAL || (stack[pos].pc >= dag->n && pos != 1) ) ) {
            printf("mpp is NOT minimal, the following are compatible:\n");
            prettyprintarray(labels[pos], "child");
            prettyprintarray(labels[ stack[pos].prev ], "parent");
            
            return false;
        }
        
        if (stack[pos].label < dag->n)
            hasTaxon[stack[pos].label] = 1;
    }
    
    // lastly, make sure every taxon is present in the tree
    for (pos = 1; pos < dag->n; pos++) {
        if (hasTaxon[pos] == 0) {
            printf("mpp is NOT a perfect phylogeny, taxon of singleton %d is missing\n", pos);
            
            return false;
        }
    }
    
    return true;
}

#undef pair_index
