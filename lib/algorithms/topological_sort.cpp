/******************************************************************************
 * topological_sort.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 * Modified and extended for O'Reach by
 *      * Kathrin Hanauer
 *      * Jonathan Trummer
 *****************************************************************************/

#include <algorithm>

#include "random_functions.h"
#include "topological_sort.h"

inline void levels_helper(graph_access & G, std::vector<unsigned> & degrees,
        std::vector<unsigned int> & levels, unsigned int level,
        std::vector<NodeID> & q1, std::vector<NodeID> & q2);
inline void ext_sort_dfs(NodeID node, graph_access & G,
                      std::vector<int>    & dfsnum,
                      int                 & dfscount,
                      std::vector<NodeID> & sorted_sequence,
                      std::vector<unsigned int> & low_index,
                      std::vector<unsigned int> & min_index);


topological_sort::topological_sort() {
                
}

topological_sort::~topological_sort() {
                
}

void topological_sort::sort( graph_access & G, std::vector<NodeID> & sorted_sequence) {
        std::vector<int> dfsnum(G.number_of_nodes(), -1);
        int dfscount = 0;

        std::vector<NodeID> nodes(G.number_of_nodes());
        random_functions::permutate_vector_good(nodes, true);

        forall_nodes(G, node) {
                NodeID curNode = nodes[node];
                if(dfsnum[curNode] == -1) {
                        sort_dfs(curNode, G, dfsnum, dfscount, sorted_sequence); 
                }
        } endfor

        std::reverse(sorted_sequence.begin(), sorted_sequence.end());
}

void topological_sort::sortStartingWith( graph_access & G, std::vector<NodeID> & startNodes, std::vector<NodeID> & sorted_sequence, bool reverse) {
        std::vector<int> dfsnum(G.number_of_nodes(), -1);
        int dfscount = 0;

        // std::vector<NodeID> nodes(G.number_of_nodes());
        // random_functions::permutate_vector_good(nodes, true);

        // start the topological sort from specific nodes
        for (auto curNode : startNodes) {
                if(dfsnum[curNode] == -1) {
                        sort_dfs(curNode, G, dfsnum, dfscount, sorted_sequence); 
                }
        }
        
        if (reverse) {
                std::reverse(sorted_sequence.begin(), sorted_sequence.end());
        }
}

void topological_sort::extendedSortStartingWith( graph_access & G, std::vector<NodeID> & startNodes, std::vector<NodeID> & sorted_sequence,
            std::vector<unsigned int> & high_index, std::vector<unsigned int> & max_index, bool reverse) {
        std::vector<int> dfsnum(G.number_of_nodes(), -1);
        high_index.assign(G.number_of_nodes(), 0);
        max_index.assign(G.number_of_nodes(), 0);
        int dfscount = 0;

        random_functions::permutate_vector_good(startNodes, false);

        // start the topological sort from specific nodes
        for (auto curNode : startNodes) {
                if(dfsnum[curNode] == -1) {
                        ext_sort_dfs(curNode, G, dfsnum, dfscount, sorted_sequence, high_index, max_index);
                }
        }
        // we assume that all vertices can be reached via one of the start nodes

        if (reverse) {
                std::reverse(sorted_sequence.begin(), sorted_sequence.end());
                for (auto i = 0U; i < high_index.size(); i++) {
                        high_index[i] = sorted_sequence.size() - 1 - high_index[i];
                        max_index[i] = sorted_sequence.size() - 1 - max_index[i];
                }
        }
}


void topological_sort::sort_dfs(NodeID node, graph_access & G, 
                                std::vector<int> & dfsnum, 
                                int & dfscount,
                                std::vector<NodeID> & sorted_sequence){ 

        dfsnum[node] = dfscount++;
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                //explore edge (node, target)
                if(dfsnum[target] == -1) {
                        sort_dfs(target, G, dfsnum, dfscount, sorted_sequence); 
                } 
        } endfor
       
        //return from call of node node
        sorted_sequence.push_back(node);
}

void ext_sort_dfs(NodeID node, graph_access & G,
                                std::vector<int> & dfsnum,
                                int & dfscount,
                                std::vector<NodeID> & sorted_sequence,
                                std::vector<unsigned int> & low_index,
                                std::vector<unsigned int> & min_index) {

        dfsnum[node] = dfscount++;
        low_index[node] = sorted_sequence.size();
        min_index[node] = low_index[node];
        std::vector<NodeID> targets;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                //explore edge (node, target)
                if(dfsnum[target] == -1) {
                        targets.push_back(target);
                } else {
                    if (min_index[target] < min_index[node]) {
                        min_index[node] = min_index[target];
                    }
                }
        } endfor

        random_functions::permutate_vector_fast(targets, false);

        for (auto target : targets) {
                if (dfsnum[target] == -1) {
                        ext_sort_dfs(target, G, dfsnum, dfscount, sorted_sequence, low_index, min_index);
                }
                if (min_index[target] < min_index[node]) {
                    min_index[node] = min_index[target];
                }
        }

        //return from call of node node
        sorted_sequence.push_back(node);
}


void topological_sort::levels(graph_access & G, graph_access & G_r, std::vector<unsigned int> & levels, const std::vector<NodeID> &sources) {
        std::vector<unsigned int> degrees(G.number_of_nodes(), 0);
        for (NodeID n = 0; n < G.number_of_nodes(); n++) {
                degrees[n] = G_r.getNodeDegree(n);
        }

        unsigned int level = 0;
        std::vector<NodeID> q1;
        std::vector<NodeID> q2;
        for (auto el : sources) {
                q1.push_back(el);
        }
        unsigned int sum = 0;

        while (q1.size() > 0 || q2.size() > 0) {
                sum += q1.size() + q2.size();
                if (level % 2 == 0) {
                        levels_helper(G, degrees, levels, level, q1, q2);
                } else {
                        levels_helper(G, degrees, levels, level, q2, q1);
                }
                level++;
        }
}

inline void levels_helper(graph_access & G, std::vector<unsigned> & degrees, std::vector<unsigned int> & levels, unsigned int level,
                std::vector<NodeID> & q1, std::vector<NodeID> & q2) {

        while(q1.size() > 0) {
                NodeID next = q1.back();
                q1.pop_back();
                levels[next] = level;
                forall_out_edges(G, e, next) {
                        NodeID target = G.getEdgeTarget(e);
                        degrees[target]--;
                        if (degrees[target] == 0) {
                                q2.push_back(target);
                        }
                } endfor
        }              
}
