#ifndef __QTABLE_H__
#define __QTABLE_H__

// QTable wraps a L*n table with the following definition:
// Q[p][i] = the smallest x >= p such that T[x][i] == 1

/*
??? bucket sort output lists of pc's ???
*/

class QTable {
private:
    vector<vector<int> > table; // Q
    vector<int> indexmap;       // index mapping from lexicographic order to cluster size then lexicographic order
    
public:
    QTable();
    
    // constructs table from lexicographically-ordered proper cluster list
    void construct(const vector<vector<int> > &propercluster);
    
    // externally initialize index mapping
    // used to handle the fact that propercluster is resorted after construction
    void setIndexMap(int i, int pc);
    
    // find proper cluster index, O(n)
    int lookup(const vector<int> &pcluster) const;
    
    // finds complement of pcluster in superpcluster if it exists, O(n)
    int lookupComplement(const vector<int> &pcluster, const vector<int> &superpcluster) const;
    
    // finds full complement of pcluster if it exists in O(n)
    int lookupComplement(const vector<int> &pcluster) const;
    
    // find all proper cluster indeces in decomp, O(n)
    // returns whether a non-proper cluster is found
    bool lookupParallel(const vector<int> &classes, vector<int> &decomp, int n, int maxC, const vector<int> &clustersize) const;
    
    // find all proper cluster indeces in decomp, O(n)
    // returns whether a non-good cluster is found
    bool lookupParallel(const vector<int> &classes, const vector<int> &good, vector<int> &decomp, int G, int G1, int n, int maxC, const vector<vector<int> > &propercluster) const;
};

#endif
