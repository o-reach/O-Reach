/******************************************************************************
 * topological_sort.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 * Modified and extended for O'Reach by
 *      * Kathrin Hanauer
 *      * Jonathan Trummer
 *****************************************************************************/

#ifndef TOPOLOGICAL_SORT_GB9FC2CZ
#define TOPOLOGICAL_SORT_GB9FC2CZ

#include <stack>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"

class topological_sort {
public:
        topological_sort();
        virtual ~topological_sort();

        void sort( graph_access & SG, std::vector<NodeID> & sorted_sequence);
        void sortStartingWith( graph_access & SG, std::vector<NodeID> & startNodes, std::vector<NodeID> & sorted_sequence, bool reverse = true);
        void extendedSortStartingWith( graph_access & SG, std::vector<NodeID> & startNodes, std::vector<NodeID> & sorted_sequence,
            std::vector<unsigned int> & high_index, std::vector<unsigned int> & max_index, bool reverse = true);

        void sort_dfs(NodeID node, graph_access & G,
                      std::vector<int>    & dfsnum,
                      int                 & dfscount,
                      std::vector<NodeID> & sorted_sequence);

        void levels(graph_access & G, graph_access & G_r, std::vector<unsigned int> & levels, const std::vector<NodeID> & sources);
};


#endif /* end of include guard: TOPOLOGICAL_SORT_GB9FC2CZ */
