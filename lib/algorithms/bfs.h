/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Jonathan Trummer
 *****************************************************************************/

#pragma once

#include <unordered_set>
#include <vector>
#include <boost/circular_buffer.hpp>
#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "lib/reachability/Query.h"

class bfs {
    public:
        bfs();
        bfs(graph_access & G);
        ~bfs();

        unsigned int explore_neighbors(std::vector<char> & neighbors, NodeID start, std::vector<NodeID> & touched);
        unsigned int explore(NodeID start);
        std::vector<NodeID> & getTouched() {
            return touched;
        }
        bool standard_w_top_order(NodeID s, NodeID t, std::vector<NodeID> & node_order);
        bool standard(NodeID s, NodeID t);
        void doCleanup() {
            cleanup();
        }

    private:
        graph_access & G;
        std::vector<char> seen;
        std::vector<NodeID> touched;
        boost::circular_buffer<NodeID> queue;

        inline void cleanup() {
            for (auto el : touched) {
                seen[el] = false;
            }
            touched.clear();
            queue.clear();
        }
};

class bibfs {
    public:
        bibfs(graph_access * G, graph_access * G_r);
        ~bibfs();
        bool standard(NodeID s, NodeID t);
        bool standard_w_top_order(NodeID s, NodeID t, std::vector<NodeID> & node_order);

    private: 
        graph_access * G;
        graph_access * G_r;

        std::vector<char> seen_forward;
        std::vector<char> seen_backward;

        std::vector<NodeID> touched_forward;
        std::vector<NodeID> touched_backward;

        boost::circular_buffer<NodeID> queue_forward;
        boost::circular_buffer<NodeID> queue_backward;

        // resets seen_forward/backward
        inline void cleanup() {
            for (auto el : touched_backward) {
                seen_backward[el] = false;
            }

            for (auto el : touched_forward) {
                seen_forward[el] = false;
            }
            touched_forward.clear();
            touched_backward.clear();
            queue_forward.clear();
            queue_backward.clear();
        }

        inline bool bidirectional_step(
                graph_access * G, boost::circular_buffer<NodeID> & q, 
                std::vector<char> & seen_s, std::vector<char> & seen_t, std::vector<NodeID> & touched_s
            ) {

            if (q.empty()) {
                return false;
            }

            NodeID next = q.front();
            // std::cout << "starting from " << next << "\n";
            q.pop_front();
            forall_out_edges((*G), edge, next) {
                NodeID t = G->unsafe_getEdgeTarget(edge);
                if (seen_t[t]) {
                    // std::cout << "already saw vertex " << t << " in other direction\n";
                    return true;
                }
                if (!seen_s[t]) {
                    // std::cout << "found new vertex " << t << "\n";
                    seen_s[t] = true;
                    touched_s.push_back(t);
                    q.push_back(t);
                }
            } endfor 

            return false;
        }

        // bidirectional step taking into account the topological node orders when considering nodes
        inline bool bidirectional_step_w_top_order(
                graph_access * G, boost::circular_buffer<NodeID> & q, 
                std::vector<char> & seen_s, std::vector<char> & seen_t, 
                std::vector<NodeID> & touched_s, 
                std::pair<unsigned int, unsigned int> & range, std::vector<NodeID> & node_order
            ) {

            if (q.empty()) {
                return false;
            }

            NodeID next = q.front();
            // std::cout << "starting from " << next << "\n";
            q.pop_front();
            forall_out_edges((*G), edge, next) {
                NodeID t = G->unsafe_getEdgeTarget(edge);

                if (node_order[t] < range.first || node_order[t] > range.second) {
                    continue;
                }

                if (seen_t[t]) {
                    // std::cout << "already saw vertex " << t << " in other direction\n";
                    return true;
                }
                if (!seen_s[t]) {
                    // std::cout << "found new vertex " << t << "\n";
                    seen_s[t] = true;
                    touched_s.push_back(t);
                    q.push_back(t);
                }
            } endfor 

            return false;
        }
};
