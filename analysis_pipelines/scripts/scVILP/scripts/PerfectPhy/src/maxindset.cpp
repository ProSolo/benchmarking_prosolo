//////// Maximal Independent Set Algorithm by Tsukiyama et al. ////////

void BACKTRACK(vector<vector<int> > &graph, vector<int> &mis, int i, int n, int &count, vector<vector<int> > &output) {
    int c, x = i+1;
    bool f;
    vector<int> bucket;
    vector<int>::iterator y, z;
    
    if (x < n) { // if current scan is not at the end
        
        // count the number of neighbors of x which are in mis
        c = 0;
        for (y = graph[x].begin(); y != graph[x].end(); ++y) {
            if (*y < x && mis[*y] == 0) {
                c++;
            }
        }
        
        if (c == 0) {
            // if no neighbors in mis, add it
            for (y = graph[x].begin(); y != graph[x].end(); ++y) {
                if (*y < x) {
                    mis[*y]++;
                }
            }
            
            // recurse on vertices larger than x
            BACKTRACK(graph, mis, x, n, count, output);
            
            // then remove x
            for (y = graph[x].begin(); y != graph[x].end(); ++y) {
                if (*y < x) {
                    mis[*y]--;
                }
            }
        }
        else {
            // x has neighbors, so don't include it
            mis[x] = c;
            
            // recurse on vertices larger than x
            BACKTRACK(graph, mis, x, n, count, output);
            
            mis[x] = 0;
            
            f = true;
            bucket.reserve(graph[x].size());
            
            // now attempt to insert x
            for (y = graph[x].begin(); y != graph[x].end(); ++y) {
                if (*y < x) {
                    if (mis[*y] == 0) {
                        // for each neighbor of x in mis less than x
                        // remove their neighbors from mis
                        
                        bucket.push_back(*y);
                        for (z = graph[*y].begin(); z != graph[*y].end(); ++z) {
                            if (*z < x) {
                                mis[*z]--;
                                if (mis[*z] == 0) {
                                    f = false;
                                }
                            }
                        }
                    }
                    mis[*y]++;
                }
            }
            
            // if all 2-away neighbors <= i are not in mis, recurse on x
            if (f) {
                BACKTRACK(graph, mis, x, n, count, output);
            }
            
            // remove x
            for (y = graph[x].begin(); y != graph[x].end(); ++y) {
                if (*y < x) {
                    mis[*y]--;
                }
            }
            
            // add back 2-away neighbors of x
            for (y = bucket.begin(); y != bucket.end(); ++y) {
                for (z = graph[*y].begin(); z != graph[*y].end(); ++z) {
                    if (*z < x) {
                        mis[*z]++;
                    }
                }
            }
        }
    }
    else {
        // bottom of recursion, found a complete mis
        output.push_back(mis);
    }
}

void MIS(vector<vector<int> > &graph, vector<vector<int> > &output) {
    int count = 0, n = graph.size();
    
    vector<int> mis(n, 0);
    
    BACKTRACK(graph, mis, 0, n, count, output);
}


//////// Maximal Independent Set Algorithms by Kristian ////////

#ifndef KMAXISBRANCHTYPE
    #define KMAXISBRANCHTYPE 3
#endif
#define KMAXISCOND   (KMAXISBRANCHTYPE == 3)
#define KMAXISFIRST  (KMAXISBRANCHTYPE == 1)
#define KMAXISSECOND (KMAXISBRANCHTYPE == 2)

const bool kmaxisDebugPrint = false;

void KMAXISrecurse(vector<vector<int> > &graph, vector<int> &removed, int unseen, vector<int> &mis, int k, vector<vector<int> > &output) {
    if (kmaxisDebugPrint) {
        printf("{ k = %d\n", k);
    }
    
    if (0 == k || 0 == unseen) {//count(removed.begin(), removed.end(), false)) {
        //prettyprintarray(mis, "MIS");
        //prettyprintarray(removed, "rem");
        
        //check if (mis is actually an k-MIS)
        
        output.push_back(mis);
        
        if (kmaxisDebugPrint) {
            printf("} k = %d\n", k);
        }
        
        return;
    }
    
    unsigned int vmax = 0, maxdeg = 0,
                 vmin = 0, mindeg = -1;
    
    // find the minimum and maximum degree vertices
    for (unsigned int v = 0; v < graph.size(); v++) {
        if (removed[v]) {
            continue;
        }
        
        unsigned int vdeg = 0;
        
        for (unsigned int nbr = 0; nbr < graph[v].size(); nbr++) {
            if (!removed[ graph[v][nbr] ]) {
                vdeg++;
            }
        }
        
        if (vdeg >= maxdeg) {
            vmax = v;
            maxdeg = vdeg;
        }
        if (vdeg <= mindeg) {
            vmin = v;
            mindeg = vdeg;
        }
    }
    
    if ( KMAXISFIRST || (KMAXISCOND && maxdeg >= (unsigned) k - 1 )) {
        vector<int> neighbors;
        
        // first add vmax to the mis,
        mis[vmax] = 0; 
        
        if (kmaxisDebugPrint) {
            printf("+ %u (maxdeg %u)\n", vmax, maxdeg);
        }
        
        // remove it and its neighbors from the graph,
        removed[vmax] = true;
        unseen--;
        
        for (unsigned int nbr = 0; nbr < graph[vmax].size(); nbr++) {
            int w = graph[vmax][nbr];
            
            if (!removed[w]) {
                neighbors.push_back(w);
                removed[w] = true;
                unseen--;
                
                if (kmaxisDebugPrint) {
                    printf(" - %d\n", w);
                }
            }
        }
        
        // and recurse on the rest of the graph
        KMAXISrecurse(graph, removed, unseen, mis, k-1, output);
        
        // then try to find all the mis excluding vmax
        mis[vmax] = 1;
        
        if (kmaxisDebugPrint) {
            printf("- %u\n", vmax);
        }
        
        // by putting its neighbors back in
        for (unsigned int w = 0; w < neighbors.size(); w++) {
            removed[ neighbors[w] ] = false;
            unseen++;
            if (kmaxisDebugPrint) {
                printf(" + %u\n", w);
            }
        }
        
        neighbors.clear();
        
        // and recursing
        KMAXISrecurse(graph, removed, unseen, mis, k, output);
        
        removed[vmax] = false;
        unseen++;
    }
    else if (KMAXISSECOND || KMAXISCOND) {
        vector<int> Vneighbors, Wneighbors;
        
        // first add vmin to the mis
        mis[vmin] = 0;
        
        if (kmaxisDebugPrint) {
            printf("+ %u (mindeg %u)\n", vmin, mindeg);
        }
        
        // thus remove its neighbors from the graph
        removed[vmin] = true;
        unseen--;
        
        for (unsigned int nbr = 0; nbr < graph[vmin].size(); nbr++) {
            int w = graph[vmin][nbr];
            
            if (!removed[w]) {
                Vneighbors.push_back(w);
                removed[w] = true;
                unseen--;
                
                if (kmaxisDebugPrint) {
                    printf(" - %d\n", w);
                }
            }
        }
        
        // and recurse on the rest
        KMAXISrecurse(graph, removed, unseen, mis, k-1, output);
        
        mis[vmin] = 1;
        
        if (kmaxisDebugPrint) {
            printf("- %u\n", vmin);
        }
        
        removed[vmin] = false;
        unseen++;
        
        for (unsigned int w = 0; w < Vneighbors.size(); w++) {
            removed[ Vneighbors[w] ] = false;
            unseen++;
            
            if (kmaxisDebugPrint) {
                printf(" + %u\n", w);
            }
        }
        
        // then consider adding each of its neighbors into the mis
        for (unsigned int wi = 0; wi < Vneighbors.size(); wi++) {
            
            int w = Vneighbors[wi];
            
            // first add it to the mis,
            mis[w] = 0;
            
            if (kmaxisDebugPrint) {
                printf("+ %d\n (nbr)", w);
            }
            
            // remove its neighbors from the graph,
            removed[w] = true;
            unseen--;
            
            for (unsigned int nbr = 0; nbr < graph[w].size(); nbr++) {
                int u = graph[w][nbr];
                
                if (!removed[u]) {
                    Wneighbors.push_back(u);
                    removed[u] = true;
                    unseen--;
                    
                    if (kmaxisDebugPrint) {
                        printf(" - %d\n", u);
                    }
                }
            }
            
            // then recurse on the rest
            KMAXISrecurse(graph, removed, unseen, mis, k-1, output);
            
            mis[w] = 1;
            removed[w] = false;
            unseen++;
            
            
            if (kmaxisDebugPrint) {
                printf("- %d\n", w);
            }
            
            for (unsigned int u = 0; u < Wneighbors.size(); u++) {
                removed[ Wneighbors[u] ] = false;
                unseen++;
                
                if (kmaxisDebugPrint) {
                    printf(" + %u\n", u);
                }
            }
            
            Wneighbors.clear();
        }
    }
    
    if (kmaxisDebugPrint) {
        printf("} k = %d\n", k);
    }
}

void KMAXIS(vector<vector<int> > &graph, int k, vector<vector<int> > &output) {
    int n = graph.size();
    
    vector<int> mis(n,1), removed(n,false);
    
    KMAXISrecurse(graph, removed, n, mis, k /* - 1*/, output);
}
