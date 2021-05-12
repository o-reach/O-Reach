/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Kathrin Hanauer
 *      Jonathan Trummer
 *****************************************************************************/

#pragma once
#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/algorithms/strongly_connected_components.h"
#include "lib/algorithms/topological_sort.h"
#include "lib/algorithms/bfs.h"
#include "extern/KaHIP/lib/tools/random_functions.h"
#include "extern/KaHIP/lib/data_structure/union_find.h"
#include "app/reachabilityconfig.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <boost/circular_buffer.hpp>
#include "lib/reachability/Query.h"

#ifdef RDEBUG
    #include <iomanip>
#endif

#ifndef NUM_SUPPORTS
#define NUM_SUPPORTS 8U
#endif

#define NUM_TOPSORTS 4U

#define PRUNE_WITH_LEVELS
//#define RESTRICT_TO_LARGE

enum FallbackMode {FB_PIBFS, FB_PDFS, FB_NONE};

// store information regarding node in *one* struct
// => query takes at most two memory lookups, one for
// NodeInfo of s, one for t
struct NodeInfo {
#if NUM_SUPPORTS > 32
    typedef uint64_t SVType;
#elif NUM_SUPPORTS > 16
    typedef uint32_t SVType;
#elif NUM_SUPPORTS > 8
    typedef uint16_t SVType;
#else
    typedef uint8_t SVType;
#endif
    const static unsigned int unset = std::numeric_limits<unsigned int>::max();
    const static unsigned int isolated = std::numeric_limits<unsigned int>::max() - 1;
    unsigned int component_number{unset};
    unsigned int pos[NUM_TOPSORTS]      =  { 0 };
    unsigned int high[(NUM_TOPSORTS)/2] =  { 0 };
    unsigned int low[(NUM_TOPSORTS)/2]  =  { 0 };
    unsigned int max[(NUM_TOPSORTS)/2]  =  { 0 };
    unsigned int min[(NUM_TOPSORTS)/2]  =  { 0 };
    unsigned int fw_level{ unset };
    unsigned int bw_level{ unset };
    SVType fw{0};
    SVType bw{0};

    NodeInfo() {}
    void setIsolated() {
        component_number = isolated;
        fw_level = 0;
        bw_level = 0;
    }

    inline
    bool isIsolated() const {
        return component_number == isolated;
    }

    void setSource() {
        fw_level = 0;
    }

    void setSink() {
        bw_level = 0;
    }

    bool isSource() const {
        return fw_level == 0;
    }

    bool isSink() const {
        return bw_level == 0;
    }

    unsigned int minOutReachability() const {
        return std::max({high[0] - pos[0], high[1] - pos[1], bw_level}) + 1;
    }

    unsigned int minInReachability() const {
        return std::max({ pos[2] - low[0], pos[3] - low[1], fw_level}) + 1;
    }

    unsigned int maxOutReachability() const {
        return std::min(max[0] - pos[0], max[1] - pos[1]) + 1;
    }

    unsigned int maxInReachability() const {
        return std::min(pos[2] - min[0], pos[3] - min[1]) + 1;
    }
};

inline
std::ostream &operator<<(std::ostream &os, const NodeInfo &n) {
    os << "(cc=" << n.component_number
        << ", source=" << (n.isSource() ? "y" : "n")
        << ", sink=" << (n.isSink() ? "y" : "n")
        << ", isolated=" << (n.isIsolated() ? "y" : "n")
        << ", fw level=" << n.fw_level
        << ", bw level=" << n.bw_level
        << ", ts 0: " << n.pos[0] << "(hi " << n.high[0] << ", max " << n.max[0] << ")"
        << ", ts 1: " << n.pos[1] << "(hi " << n.high[1] << ", max " << n.max[1] << ")"
        << ", ts 2: " << n.pos[2] << "(lo " << n.low[0] << ", min " << n.min[0] << ")"
        << ", ts 3: " << n.pos[3] << "(lo " << n.low[1] << ", min " << n.min[1] << ")"
//        << ", " << n.minOutReachability() << " <= r+ <= " << n.maxOutReachability()
//        << ", " << n.minInReachability()  << " <= r- <= " << n.maxInReachability()
        << ", #fw reach=" << __builtin_popcount(n.fw) //<< std::hex << static_cast<unsigned>(n.fw)
        << ", #bw reach=" << __builtin_popcount(n.bw) //<< std::hex << static_cast<unsigned>(n.bw)
        //<< std::dec
        << ")";
    return os;
}

class oreach {
    public:
        oreach(graph_access & G, ReachabilityConfig & config);
        ~oreach();

        void initialize();

        template<FallbackMode fm>
        void inline query(std::vector<Query> &queries) {
            for (auto &q : queries) {
#ifdef VERBOSE
              std::cout << "Query: " << q.s << " ->? " << q.t << std::endl;
#endif
               singleQuery<fm>(q);
               q.setBinaryFromResult();

#ifdef RDEBUG
               query_result_stats[q.resultToIndex()]+=1;
#endif
            }
        }

        template<FallbackMode fm, bool pbfs_subquery = false>
        inline
        void singleQuery(Query &q) /*const*/ {
#ifdef VERBOSE
              std::cout << "  Single query:\nsource: " << nodes[q.s] << "\ntarget: " << nodes[q.t] << std::endl;
#endif
            if (!pbfs_subquery && q.s == q.t) {
                #ifdef VERBOSE
                    std::cout << "  Observation: Same vertex" << std::endl;
                #endif
                q.result = Query::Result::SAME_VERTEX;
                return;
            }

            auto &s = nodes[q.s];
            auto &t = nodes[q.t];

            // s is sink or t is source
            /*if (nodes[q.s].isSink() || t.isSource()) {
                q.result = Query::Result::SRC_SNK;
            }*/

            // topological order of t is smaller than s
            // or topological level of t is smaller/equal to s
            // => there can't be a path from s to t
            if (!pbfs_subquery && s.fw_level >= t.fw_level) {
                #ifdef VERBOSE
                                    std::cout << "  Observation: Forward level of s is at least forward level of t\n";
                #endif
                q.result = Query::Result::TS_FW_LEVELS;
                return;
            }

            if(!pbfs_subquery && s.bw_level <= t.bw_level) {
                #ifdef VERBOSE
                                    std::cout << "  Observation: Backward level of s is at most backward level of t\n";
                #endif
                q.result = Query::Result::TS_BW_LEVELS;
                return;
            }

            // Observation 1: s in backward and t in forward => True
            if (!pbfs_subquery && (s.bw & t.fw))  {
#ifdef VERBOSE
                std::cout << "  Observation: Backward and forward reachable from support\n";
                const auto result =  s.bw & t.fw;
                std::cout << __builtin_popcount(result) << " matches for SV1:";
                for (auto i = 0U; i < NUM_SUPPORTS; i++) {
                    if (result & (1 << i)) {
                        std::cout << " " << i;
                    }
                }
                std::cout << std::endl;
#endif
                q.result = Query::Result::SUPPORT_POS;
                return;
            }


            // TS 0
            if(s.pos[0] >= t.pos[0]) {
#ifdef RESTRICT_TO_LARGE
            if (s.component_number >= large_components_end
                    || t.component_number >= large_components_end) {
#ifdef VERBOSE
                std::cout << "  One of s or t is in a small component.\n";
#endif
                goto small_component;
            }
#endif
#ifdef VERBOSE
                std::cout << "  Observation: s after t in TS 0\n";
#endif
                q.result = Query::Result::TS_NUM;
                return;
            }

            if (t.pos[0] <= s.high[0]) {
#ifdef RESTRICT_TO_LARGE
            if (s.component_number >= large_components_end
                    || t.component_number >= large_components_end) {
#ifdef VERBOSE
                std::cout << "  One of s or t is in a small component.\n";
#endif
                goto small_component;
            }
#endif
#ifdef VERBOSE
                std::cout << "  Observation: t within s.high in TS 0\n";
#endif
                q.result = Query::Result::TS_HI;
                return;
            }

            if (t.pos[0] > s.max[0]) {
#ifdef RESTRICT_TO_LARGE
            if (s.component_number >= large_components_end
                    || t.component_number >= large_components_end) {
#ifdef VERBOSE
                std::cout << "  One of s or t is in a small component.\n";
#endif
                goto small_component;
            }
#endif
#ifdef VERBOSE
                std::cout << "  Observation: t beyond s.max in TS 0\n";
#endif
                q.result = Query::Result::TS_MAX;
                return;
            }

            if (t.pos[0] == s.max[0]) {
#ifdef RESTRICT_TO_LARGE
            if (s.component_number >= large_components_end
                    || t.component_number >= large_components_end) {
#ifdef VERBOSE
                std::cout << "  One of s or t is in a small component.\n";
#endif
                goto small_component;
            }
#endif
#ifdef VERBOSE
                std::cout << "  Observation: t at s.max in TS 0\n";
#endif
                q.result = Query::Result::TS_EQMAX;
                return;
            }
            // END TS 0

            // Observation 2: s in forward but t not in forward => false
            if (s.fw & ~t.fw) {
#ifdef VERBOSE
                std::cout << "  Observation: s forward reachable by support but not t\n";
                const auto result = s.fw & ~t.fw;
                std::cout << __builtin_popcount(result) << " matches for SV2:";
                for (auto i = 0U; i < NUM_SUPPORTS; i++) {
                    if (result & (1 << i)) {
                        std::cout << " " << i;
                    }
                }
                std::cout << std::endl;
#endif
                q.result = Query::Result::SUPPORT_NEG_FW;
                return;
            }

            // Observation 3: t in backward but s not in backward => false
            // (s.bw NAND t.bw) AND t.bw
            //if ((~(s.bw & t.bw) & t.bw) > 0) {
            //if ((s.bw ^ t.bw) & t.bw) {
            if (~s.bw & t.bw) {
#ifdef VERBOSE
                std::cout << "  Observation: t backward reachable by support but not s\n";
                auto result = ~s.bw & t.bw;
                std::cout << __builtin_popcount(result) << " matches for SV3:";
                for (auto i = 0U; i < NUM_SUPPORTS; i++) {
                    if (result & (1 << i)) {
                        std::cout << " " << i;
                    }
                }
                std::cout << std::endl;
#endif
                q.result = Query::Result::SUPPORT_NEG_BW;
                return;
            }

#ifdef RESTRICT_TO_LARGE
            if (s.component_number >= large_components_end
                    || t.component_number >= large_components_end) {
#ifdef VERBOSE
                std::cout << "  One of s or t is in a small component.\n";
#endif
                goto small_component;
            }
#endif

            // TS 1
            if(s.pos[1] >= t.pos[1]) {
#ifdef VERBOSE
                std::cout << "  Observation: s after t in TS 1\n";
#endif
                q.result = Query::Result::TS_NUM;
                return;
            }

            if (t.pos[1] <= s.high[1]) {
#ifdef VERBOSE
                std::cout << "  Observation: t within s.high in TS 1\n";
#endif
                q.result = Query::Result::TS_HI;
                return;
            }

            if (t.pos[1] > s.max[1]) {
#ifdef VERBOSE
                std::cout << "  Observation: t beyond s.max in TS 1\n";
#endif
                q.result = Query::Result::TS_MAX;
                return;
            }

            if (t.pos[1] == s.max[1]) {
#ifdef VERBOSE
                std::cout << "  Observation: t at s.max in TS 1\n";
#endif
                q.result = Query::Result::TS_EQMAX;
                return;
            }
            // END TS 1

            // TS 2
            if(s.pos[2] >= t.pos[2]) {
#ifdef VERBOSE
                std::cout << "  Observation: s after t in TS 2\n";
#endif
                q.result = Query::Result::TS_NUM;
                return;
            }

            if (s.pos[2] >= t.low[0]) {
#ifdef VERBOSE
                std::cout << "  Observation: s within t.low in TS 2\n";
#endif
                q.result = Query::Result::TS_LO;
                return;
            }

            if (s.pos[2] < t.min[0]) {
#ifdef VERBOSE
                std::cout << "  Observation: s beyond t.min in TS 2\n";
#endif
                q.result = Query::Result::TS_MIN;
                return;
            }

            if (s.pos[2] == t.min[0]) {
#ifdef VERBOSE
                std::cout << "  Observation: s at s.min in TS 2\n";
#endif
                q.result = Query::Result::TS_EQMIN;
                return;
            }
            // END TS 2

            // TS 3
            if(s.pos[3] >= t.pos[3]) {
#ifdef VERBOSE
                std::cout << "  Observation: s after t in TS 3\n";
#endif
                q.result = Query::Result::TS_NUM;
                return;
            }

            if (s.pos[3] >= t.low[1]) {
#ifdef VERBOSE
                std::cout << "  Observation: s within t.low in TS 3\n";
#endif
                q.result = Query::Result::TS_LO;
                return;
            }

            if (s.pos[3] < t.min[1]) {
#ifdef VERBOSE
                std::cout << "  Observation: s beyond t.min in TS 3\n";
#endif
                q.result = Query::Result::TS_MIN;
                return;
            }

            if (s.pos[3] == t.min[1]) {
#ifdef VERBOSE
                std::cout << "  Observation: s at s.min in TS 3\n";
#endif
                q.result = Query::Result::TS_EQMIN;
                return;
            }
            // END TS 3

#ifdef RESTRICT_TO_LARGE
small_component:
#endif
            if (pbfs_subquery) {
                q.result = Query::Result::UNKNOWN;
                return;
            }

            // both vertices in different connected components
            if (s.component_number != t.component_number) {
                #ifdef VERBOSE
                    std::cout << "  Observation: Different connected components" << std::endl;
                #endif
                q.result = Query::Result::DIFF_WCC;
                return;
            }

            // if we come this far the query isn't answered by an observation...
            //return Query::Result::UNKNOWN;

            switch (fm) {
                case FB_PIBFS:
#ifdef VERBOSE
                    std::cout << "  Couldn't answer query via observation, falling back to pBiBFS.\n";
#endif
                    //if (checkEdge(q)) {
                    //    return;
                    //}
                    pruning_bfs(q);
                    break;
                case FB_PDFS:
#ifdef VERBOSE
                    std::cout << "  Couldn't answer query via observation, falling back to pDFS.\n";
#endif
                    pruning_dfs(q);
                    break;
                default:
#ifdef VERBOSE
                    std::cout << "  Couldn't answer query via observation, leaving query unanswered.\n";
#endif
                    break;
            }
        }

        inline
        bool checkEdge(Query &q) {
            if (condensation.getNodeDegree(q.s) <=
                    condensation_r.getNodeDegree(q.t)) {
                forall_out_edges(condensation, edge, q.s) {
                    if (condensation.unsafe_getEdgeTarget(edge) == q.t) {
                        q.result = Query::Result::EDGE;
                        return true;
                    }
                } endfor
                return false;
            } else {
                forall_out_edges(condensation_r, edge, q.t) {
                    if (condensation_r.unsafe_getEdgeTarget(edge) == q.s) {
                        q.result = Query::Result::EDGE;
                        return true;
                    }
                } endfor
                return false;
            }
        }

        void debug();

        std::vector<unsigned int> getObservationDistribution() {
            //std::vector<unsigned int> dist(Query::castResult(Query::Result::LAST),0);
            //for (unsigned int i = 0; i < Query::castResult(Query::Result::LAST); i++) {
            //    dist[i] = query_result_stats[i];
            //}
            //return dist;
            return std::vector<unsigned int>(std::begin(query_result_stats), std::end(query_result_stats));
        }

        void reset() {
            for (unsigned int i = 0; i < Query::castResult(Query::Result::LAST); i++) {
                query_result_stats[i] = 0;
            }
            num_pruned_search_steps = 0;
        }

        void setFallbackMode(FallbackMode mode) {
            fallback = mode;
        }


    private:
        void chooseSupportingVertex();
        void supportingVertexStep(
            unsigned int round,
            bfs & backward_bfs, std::pair<NodeID, NodeID> range,
            std::vector<NodeID> & sortedNodes
        );
        void sampleSupports( bfs & backward_bfs);
        void sampleSupportsFromLevels(bfs & backward_bfs);
        void sampleSupportsFromLevels2(bfs & backward_bfs);
        void pickSupports(bfs &backward_bfs);

        void doTopologicalSorts();
        void determineComponents();

         // resets seen_forward/backward
        void inline cleanupBfs() {
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

        // pruning bfs
        inline
        void pruning_bfs(Query &q) {
            // already tested
            assert(q.s != q.t);

            queue_forward.push_back(q.s);
            queue_backward.push_back(q.t);

            seen_forward[q.s] = true;
            seen_backward[q.t] = true;
            touched_forward.push_back(q.s);
            touched_backward.push_back(q.t);

            bool res = false;
#ifdef PRUNE_WITH_LEVELS
            auto fw_level_s = nodes[q.s].fw_level;
            auto bw_level_s = nodes[q.s].bw_level;

            auto fw_level_t = nodes[q.t].fw_level;
            auto bw_level_t = nodes[q.t].bw_level;

            auto cur_level_queue_fw_size = 1;
            auto cur_level_queue_bw_size = 1;
            #ifdef VERBOSE
                std::cout  << "Forward level of s is " << fw_level_s << std::endl;
                std::cout  << "Backward level of s is " << bw_level_s << std::endl;
                std::cout  << "Forward level of t is " << fw_level_t << std::endl;
                std::cout  << "Backward level of t is " << bw_level_t << std::endl;

                auto old_pruned_bfs_steps = num_pruned_search_steps;
            #endif
#endif

            while (!queue_forward.empty() && !queue_backward.empty() && !res) {
#ifdef PRUNE_WITH_LEVELS
                res = pruning_bibfs_step(condensation, queue_forward, seen_forward, seen_backward, touched_forward,
                        fw_level_t, bw_level_t, q.t);
                cur_level_queue_fw_size--;
                if (cur_level_queue_fw_size == 0) {
                    cur_level_queue_fw_size = queue_forward.size();
                    fw_level_s++;
                    bw_level_s--;
                    #ifdef VERBOSE
                        std::cout  << "Increased forward level of S to " << fw_level_s << std::endl;
                        std::cout  << "Decreased backward level of S to " << bw_level_s << std::endl;
                        std::cout  << "Size of forward queue is " << cur_level_queue_fw_size << std::endl;
                    #endif
                }
#else
                res = pruning_bibfs_step(condensation, queue_forward, seen_forward, seen_backward, touched_forward,
                        q.t);
#endif
                #ifdef RDEBUG
                  num_pruned_search_steps +=1;
                #endif
                if (!res && !queue_forward.empty()) {
#ifdef PRUNE_WITH_LEVELS
                    res = pruning_bibfs_step(condensation_r, queue_backward, seen_backward, seen_forward, touched_backward,
                            fw_level_s, bw_level_s, q.s, true);
                    cur_level_queue_bw_size--;
                    if (cur_level_queue_bw_size == 0) {
                        cur_level_queue_bw_size = queue_backward.size();
                        fw_level_t--;
                        bw_level_t++;
                        #ifdef VERBOSE
                            std::cout  << "Decreased forward level of T to " << fw_level_t << std::endl;
                            std::cout  << "Increased backward level of T to " << bw_level_t << std::endl;
                            std::cout  << "Size of backward queue is " << cur_level_queue_bw_size << std::endl;
                        #endif
                    }
#else
                    res = pruning_bibfs_step(condensation_r, queue_backward, seen_backward, seen_forward, touched_backward,
                            q.s, true);
#endif
                #ifdef RDEBUG
                    num_pruned_search_steps +=1;
                #endif
                }
            }
           #ifdef VERBOSE
              std::cout  << "#Pruned BiBFS steps: " << num_pruned_search_steps - old_pruned_bfs_steps << std::endl;
           #endif

            cleanupBfs();
            q.result = res ?
                Query::Result::FALLBACK_POS : Query::Result::FALLBACK_NEG;
        }

        bool inline pruning_bibfs_step(
                graph_access & G, boost::circular_buffer<NodeID> & q,
                std::vector<char> & seen_s, std::vector<char> & seen_t, std::vector<NodeID> & touched_s,
#ifdef PRUNE_WITH_LEVELS
                unsigned fw_level, unsigned bw_level,
#endif
                NodeID target, bool reversed = false
                ) {

            NodeID next = q.front();
            q.pop_front();
            Query query(next, target);
            if (reversed) {
                query.s = target;
                query.t = next;
            }

            forall_out_edges(G, edge, next) {
                NodeID t = G.unsafe_getEdgeTarget(edge);
#ifdef VERBOSE
                std::cout << "Looking at neighbor " << t << " of " << next
                    << " on levels " << nodes[t].fw_level << "/" << nodes[t].bw_level
                    << "..." << std::endl;
#endif
                if (seen_t[t]) {
#ifdef RDEBUG
                    query_result_stats[Query::castResult(Query::Result::SAME_VERTEX)]++;
#endif
#ifdef VERBOSE
                    std::cout << "  Already seen " << t << " from other side." << std::endl;
#endif
                    return true;
                }

                if (!seen_s[t]) {
                    seen_s[t] = true;
                    touched_s.push_back(t);
                    if (reversed) {
#ifdef PRUNE_WITH_LEVELS
                        if (nodes[t].fw_level <= fw_level) {
#ifdef RDEBUG
                            query_result_stats[Query::castResult(Query::Result::TS_FW_LEVELS)]++;
#endif
#ifdef VERBOSE
                            std::cout << "  Forward level of " << t << ": " << nodes[t].fw_level << " <= " << fw_level << ", can't reach S." << std::endl;
#endif
                            continue;
                        }
                        if (nodes[t].bw_level >= bw_level) {
#ifdef RDEBUG
                            query_result_stats[Query::castResult(Query::Result::TS_BW_LEVELS)]++;
#endif
#ifdef VERBOSE
                            std::cout << "  Backward level of " << t << ": " << nodes[t].bw_level  << " >= " << bw_level << ", can't reach S." << std::endl;
#endif
                            continue;
                        }
#endif
                        query.t = t;
                    } else {
#ifdef PRUNE_WITH_LEVELS
                        if (nodes[t].fw_level >= fw_level) {
#ifdef RDEBUG
                            query_result_stats[Query::castResult(Query::Result::TS_FW_LEVELS)]++;
#endif
#ifdef VERBOSE
                            std::cout << "  Forward level of " << t << ": " << nodes[t].fw_level  << " >= " << fw_level << ", can't reach T." << std::endl;
#endif
                            continue;
                        }
                        if (nodes[t].bw_level <= bw_level) {
#ifdef RDEBUG
                            query_result_stats[Query::castResult(Query::Result::TS_BW_LEVELS)]++;
#endif
#ifdef VERBOSE
                            std::cout << "  Backward level of " << t << ": " << nodes[t].bw_level  << " <= " << bw_level << ", can't reach T." << std::endl;
#endif
                            continue;
                        }
#endif
                        query.s = t;
                    }

                    // query oreach for t -> target
                    // on false, no further action
                    // on true => return true, no further steps required
                    // on maybe => add to queue for later processing
                    // query object is reused -> reset
                    query.result = Query::Result::UNKNOWN;
                    singleQuery<FB_NONE, true>(query);

                    // maybe branch
                    if (!query.answeredWithObservation()) {
                        q.push_back(t);
                    } else {
#ifdef RDEBUG
                        query_result_stats[query.resultToIndex()]++;
#endif
                        if (query.isPositive()) {
                          return true;
                        }
                    }
                }
#ifdef VERBOSE
                    else { std::cout << "  Already seen " << t << " from this side." << std::endl; }
#endif
            } endfor

            return false;
        }

        // pruning DFS
        inline
        void pruning_dfs(Query &q) {
#ifdef VERBOSE
            auto old_pruned_search_steps = num_pruned_search_steps;
#endif
            // already tested
            assert(q.s != q.t);

            auto res = pruning_dfs_step(q.s, q.t);
#ifdef VERBOSE
              std::cout  << "#Pruned DFS steps: " << num_pruned_search_steps - old_pruned_search_steps << std::endl;
#endif
            cleanupBfs();
            q.result = res ?
                Query::Result::FALLBACK_POS : Query::Result::FALLBACK_NEG;
        }

        inline
        bool pruning_dfs_step(NodeID n, NodeID target) {
            if (seen_forward[n]) {
                return false;
            }
#ifdef RDEBUG
            num_pruned_search_steps++;
#endif
            seen_forward[n] = true;
            touched_forward.push_back(n);

            Query query(n, target);
            auto &nt = nodes[target];
            forall_out_edges(condensation, edge, n) {
                NodeID x = condensation.unsafe_getEdgeTarget(edge);
                if (seen_forward[x]) {
                    continue;
                }
                if (x == target) {
#ifdef RDEBUG
                    query_result_stats[Query::castResult(Query::Result::SAME_VERTEX)]++;
#endif
                    return true;
                }

                if (nodes[x].fw_level >= nt.fw_level) {
#ifdef RDEBUG
                    query_result_stats[Query::castResult(Query::Result::TS_FW_LEVELS)]++;
#endif
                    continue;
                }

                if (nodes[x].bw_level <= nt.bw_level) {
#ifdef RDEBUG
                    query_result_stats[Query::castResult(Query::Result::TS_BW_LEVELS)]++;
#endif
                    continue;
                }

                query.s = x;
                query.result = Query::Result::UNKNOWN;
                singleQuery<FB_NONE, true>(query);
                if (!query.answeredWithObservation()) {
                    if (pruning_dfs_step(x, target)) {
                        return true;
                    }
                } else {
#ifdef RDEBUG
                    query_result_stats[query.resultToIndex()]++;
#endif
                    if (query.isPositive()) {
                        return true;
                    }
                }
            } endfor
            return false;
        }

        inline void setSupport(int num,
                const std::vector<NodeID> &out_reachable,
                const std::vector<NodeID> &in_reachable) {
            for (auto n : out_reachable) {
                nodes[n].fw += 1 << num;
            }
            for (auto n : in_reachable) {
                nodes[n].bw += 1 << num;
            }
        }

        inline void clearSupport(int num,
                const std::vector<NodeID> &out_reachable,
                const std::vector<NodeID> &in_reachable) {
            for (auto n : out_reachable) {
                nodes[n].fw &= ~(1 << num);
            }
            for (auto n : in_reachable) {
                nodes[n].bw &= ~(1 << num);
            }
        }

        inline void clearSupport(int num) {
            for (auto &node : nodes) {
                node.fw &= ~(1 << num);
                node.bw &= ~(1 << num);
            }
        }

        template<bool exclude_source_sink, bool forward = true>
        inline void countVerticesPerLevel(std::vector<unsigned> &vpl) {
            vpl.assign(max_ts_level+1, 0);
            for (auto &node : nodes) {
                if (node.component_number >= large_components_end
                        || (exclude_source_sink
                            && (node.isSource() || node.isSink()))) {
                    continue;
                }
                vpl[forward ? node.fw_level : node.bw_level]++;
            }
        }

        inline
        void filterVerticesInLargeComponents(std::vector<NodeID> &vs) {
            assert(!vs.empty());
            auto i = 0U;
            auto j = vs.size() - 1;
            while (i <= j) {
                if (nodes[vs[i]].component_number >= large_components_end) {
                    std::swap(vs[i], vs[j]);
#ifdef VERBOSE
                    std::cout << "Filtering out node " << nodes[vs[j]] << std::endl;
#endif
                    j--;
                } else {
                    i++;
                }
            }
#ifdef VERBOSE
                std::cout << "Filtering terminated with i = " << i << ", j = " << j << std::endl;
#endif
            if (i < vs.size()) {
#ifdef RDEBUG
            assert(i > 0);
            std::cout << "Last vertex to keep is " << nodes[vs[i-1]] << std::endl;
            std::cout << "First vertex to remove is " << nodes[vs[i]] << std::endl;
#endif
                vs.resize(i);
            }
        }

        static const int num_of_supports{NUM_SUPPORTS};

        std::vector<NodeInfo> nodes;

        std::vector<NodeID> sorted_nodes;

        graph_access & condensation;
        ReachabilityConfig & config;
        FallbackMode fallback;

        graph_access condensation_r;
        bfs * bfs_query;

        bool initialized {false};

        // components of the condensation, NOT SSCs!
        unsigned n_nonisolated;
        unsigned max_component_number;
        std::vector<unsigned> vertices_per_component;
        unsigned large_components_end;
        std::vector<NodeID> source_nodes;
        std::vector<NodeID> sink_nodes;

        unsigned int max_ts_level{0};

        std::vector<unsigned> vertices_per_level_fw;
        std::vector<unsigned> vertices_per_level_bw;
        //std::vector<unsigned> vertices_per_level_fw_non_sourcesink;
        //std::vector<unsigned> vertices_per_level_bw_non_sourcesink;

        // for bfs usages (ie when choosing supports, pruning_bibfs fallback)
        std::vector<char> seen_forward;
        std::vector<char> seen_backward;

        std::vector<NodeID> touched_forward;
        std::vector<NodeID> touched_backward;

        boost::circular_buffer<NodeID> queue_forward;
        boost::circular_buffer<NodeID> queue_backward;

        unsigned int query_result_stats[Query::castResult(Query::Result::LAST)] = {0};
        unsigned int num_pruned_search_steps = 0;
};
