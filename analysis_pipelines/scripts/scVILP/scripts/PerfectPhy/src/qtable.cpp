#include"qtable.h"

QTable::QTable() {}

void QTable::construct(const vector<vector<int> > &propercluster) {
    // propercluster must be sorted lexicographically, NOT also by cluster size
    int properclustercount = propercluster.size(), n = propercluster[0].size();
    vector<int> lastone(n);
    
    table.resize(properclustercount);
    indexmap.assign(properclustercount, -1);
    
    for(int pc = properclustercount - 1; pc >= 0; pc--) {
        for(int t = 0; t < n; t++)
            if (propercluster[pc][t] == 1)
                lastone[t] = pc;
        
        table[pc] = lastone;
    }
}

void QTable::setIndexMap(int i, int pc) {
    indexmap[i] = pc;
}

int QTable::lookup(const vector<int> &pcluster) const {
    int pc = 0;
    unsigned int t, n = pcluster.size();
    
    // jump through Q table to find sorted position of pcluster
    for(t = 0; t < n; t++)
        if (pcluster[t] == 1)
            pc = table[pc][t];
    
    // ensure that pcluster is actually present by checking that
    // pcluster[t] == 1 iff Q[pc][t] == pc
    for(t = 0; t < n; t++)
        if ((pcluster[t] == 1) != (table[pc][t] == pc))
            return -1;
    
    // return true index
    return indexmap[pc];
}

int QTable::lookupComplement(const vector<int> &pcluster, const vector<int> &superpcluster) const {
    int pc = 0;
    unsigned int t, n = pcluster.size();
    
    #define compHas(i) (pcluster[i] == 0 && superpcluster[i] == 1)
    
    // jump through Q table to find sorted position of comp
    for(t = 0; t < n; t++)
        if (compHas(t))
            pc = table[pc][t];
    
    // ensure that comp is actually present by checking that
    // comp[t] == 1 iff Q[pc][t] == pc
    for(t = 0; t < n; t++)
        if (compHas(t) != (table[pc][t] == pc))
            return -1;
    
    // return true index
    return indexmap[pc];
    
    #undef compHas
}

int QTable::lookupComplement(const vector<int> &pcluster) const {
    int pc = 0;
    unsigned int t, n = pcluster.size();
    
    #define compHas(i) (pcluster[i] == 0)
    
    // jump through Q table to find sorted position of comp
    for(t = 0; t < n; t++)
        if (compHas(t))
            pc = table[pc][t];
    
    // ensure that comp is actually present by checking that
    // comp[t] == 1 iff Q[pc][t] == pc
    for(t = 0; t < n; t++)
        if (compHas(t) != (table[pc][t] == pc))
            return -1;
    
    // return true index
    return indexmap[pc];
    
    #undef compHas
}

bool QTable::lookupParallel(const vector<int> &classes, vector<int> &decomp, int n, int maxC, const vector<int> &clustersize) const {
    bool allpc = true;
    vector<int> pc(maxC, 0), size(maxC, 0);
    
    // jump through Q table to find sorted position of each pc
    for(int t = 0; t < n; t++) {
        unsigned int curClass = classes[t];
        int& p = pc[curClass];
        
        p = table[p][t];
    }
    
    // we do not know if each pc is present, so ensure that:
    // 1. each class is a subset of its pc
    // 2. each class is the same size as its pc
    
    for(int t = 0; t < n; t++) {
        unsigned int curClass = classes[t];
        int& p = pc[curClass];
        
        if (p == -1)
            continue; // p is invalid
        else if (p != table[p][t]) {
            allpc = false;
            p = -1; // table entry indicates that p should not contain t, thus p is invalid
        }
        
        size[curClass]++;
    }
    
    decomp.clear();
    decomp.reserve(maxC);
    
    for(int i = 0; i < maxC; i++) {
        if (pc[i] == -1)
            continue;
        
        int p = indexmap[ pc[i] ];
        
        if (size[i] == clustersize[p])
            decomp.push_back(p); // p must be equal to its pc, thus add it to decomp
        else {
            allpc = false;
            pc[i] = -1; // p is a proper subset of its pc, thus not a proper cluster
        }
    }
    
    // return whether or not all members of decomp are pc's
    return allpc;
}

bool QTable::lookupParallel(const vector<int> &classes, const vector<int> &good, vector<int> &decomp, int G, int G1, int n, int maxC, const vector<vector<int> > &propercluster) const {
    vector<int> pc(maxC, 0); // pc indices in table
    
    // jump through Q table to find sorted position of each pc
    for(int t = 0; t < n; t++) {
        unsigned int curClass = classes[t];
        
        if (propercluster[G][t] == 1 && propercluster[G1][t] == 0)
            pc[curClass] = table[ pc[curClass] ][t]; // t is in G2
        else
            pc[curClass] = -1; // t is either in G1 or not in G at all
    }
    
    decomp.reserve(maxC);
    decomp.push_back(G1);
    
    // we know each pc is present, so must only ensure goodness
    for(int i = 0; i < maxC; i++) {
        if (pc[i] == -1)
            continue;
        
        int p = indexmap[ pc[i] ];
        
        if (good[p])
            decomp.push_back(p);
        else
            return false;
    }
    
    // return that all members of decomp are good pc's
    return true;
}
