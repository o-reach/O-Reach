/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Jonathan Trummer
 *****************************************************************************/

#include "bfs.h"

// REGULAR BFS
bfs::bfs(graph_access & G): G(G), queue(G.number_of_nodes()) {
    seen.assign(G.number_of_nodes(), false);
    touched.reserve(G.number_of_nodes()/2);
}
bfs::~bfs() {}


unsigned int bfs::explore_neighbors(std::vector<char> & neighbors, NodeID start, std::vector<NodeID> & touched) {
    queue.push_back(start);
    neighbors[start] = true;
    touched.push_back(start);

    unsigned int count = 0;
    
    while (!queue.empty()) {
        NodeID next = queue.front();
        queue.pop_front();
        forall_out_edges(G, edge, next) {
            NodeID t = G.unsafe_getEdgeTarget(edge);
            if (!neighbors[t]) {
                queue.push_back(t);
                neighbors[t] = true;
                touched.push_back(t);
                count++;
            }
        } endfor

    }
    queue.clear();
    return count;
}

unsigned int bfs::explore(NodeID start) {
    return explore_neighbors(seen, start, touched);
}

bool bfs::standard_w_top_order(NodeID start, NodeID target, std::vector<NodeID> & node_order) {
    if (start==target) {
        return true;
    }

    seen[start] = true;
    touched.push_back(start);
    queue.push_back(start);

    std::pair<unsigned int, unsigned int> range = {node_order[start], node_order[target]};
    bool reachable = false;
    while (!queue.empty() && !reachable) {
        NodeID next = queue.front();
        queue.pop_front();
        forall_out_edges(G, edge, next) {
            NodeID t = G.unsafe_getEdgeTarget(edge);
            if (t == target) {
                reachable = true;
                break;
            }
            if (node_order[t] < range.first || node_order[t] > range.second) {
                continue;
            }

            if (!seen[t]) {
                seen[t] = true;
                touched.push_back(t);
                queue.push_back(t);
            }
        } endfor

    }

    cleanup();

    return reachable;
}

bool bfs::standard(NodeID start, NodeID target) {
    if (start==target) {
        return true;
    }

    seen[start] = true;
    touched.push_back(start);
    queue.push_back(start);

    bool reachable = false;
    while (!queue.empty() && !reachable) {
        NodeID next = queue.front();
        queue.pop_front();
        forall_out_edges(G, edge, next) {
            NodeID t = G.unsafe_getEdgeTarget(edge);
            if (t == target) {
                reachable = true;
                break;
            }

            if (!seen[t]) {
                seen[t] = true;
                touched.push_back(t);
                queue.push_back(t);
            }
        } endfor

    }

    cleanup();

    return reachable;
}


// BIDIRECTIONAL BFS
bibfs::bibfs(graph_access * G, graph_access * G_r): G(G), G_r(G_r), queue_forward(G->number_of_nodes()), queue_backward(G_r->number_of_nodes()) {
    seen_forward.assign(G->number_of_nodes(), false);
    seen_backward.assign(G->number_of_nodes(), false);
    touched_forward.reserve(G->number_of_nodes()/2);
    touched_backward.reserve(G->number_of_nodes()/2);
}

bibfs::~bibfs(){}

bool bibfs::standard(NodeID s, NodeID t){
    if (s==t) {
        return true;
    }
    
    queue_forward.push_back(s);
    queue_backward.push_back(t);

    seen_forward[s] = true;
    seen_backward[t] = true;
    touched_forward.push_back(s);
    touched_backward.push_back(t);


    bool res = false;
    // std::cout << "query: " << s << " -> " << t << "\n";
    while (!queue_forward.empty() && !queue_backward.empty() && !res) {
        // std::cout << "forward step \n";
        res = bidirectional_step(G, queue_forward, seen_forward, seen_backward, touched_forward);
        if (!res) {
            // std::cout << "backward step: \n";
            res = bidirectional_step(G_r, queue_backward, seen_backward, seen_forward, touched_backward);
        }
    }

    cleanup();

    return res;
}


bool bibfs::standard_w_top_order(NodeID s, NodeID t, std::vector<NodeID> &  node_order) {
    if (s==t) {
        return true;
    }

    queue_forward.push_back(s);
    queue_backward.push_back(t);

    seen_forward[s] = true;
    seen_backward[t] = true;
    touched_forward.push_back(s);
    touched_backward.push_back(t);

    bool res = false;
    std::pair<unsigned int, unsigned int> range = {node_order[s], node_order[t]};
    // std::cout << "query: " << s << " -> " << t << "\n";
    while (!queue_forward.empty() && !queue_backward.empty() && !res) {
        // std::cout << "forward step \n";
        res = bidirectional_step_w_top_order(G, queue_forward, seen_forward, seen_backward, touched_forward, range, node_order);
        if (!res) {
            res = bidirectional_step_w_top_order(G_r, queue_backward, seen_backward, seen_forward, touched_backward, range, node_order);
        }
    }

    cleanup();

    return res;
}
