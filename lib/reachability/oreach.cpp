/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Kathrin Hanauer
 *      Jonathan Trummer
 *****************************************************************************/

#include "oreach.h"
#include <algorithm>

#define SMALL_COMPONENT_RECIPROCAL 10U
#define SUPPORT_MIN_LEVEL_RECIPROCAL 4U
//#define ORDERED_SRC_SNK
#define MAX_SUPPORT_CANDIDATE_FACTOR 45

#ifndef MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS
//#define MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS 150
//#define MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS 100
#define MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS 75
#endif

//#define FREE_SPACE

struct SupportQuality {
    unsigned long sv1;
    unsigned long sv2;
    unsigned long sv3;
    SupportQuality(unsigned long q1 = 0, unsigned long q2 = 0, unsigned long q3 = 0)
        : sv1(q1), sv2(q2), sv3(q3) {}
    bool operator<(const SupportQuality &other) const {
        return sv1 < other.sv1
            && sv2 < other.sv2
            && sv3 < other.sv3;
    }

    void update(const SupportQuality &other) {
        if (other.sv1 > sv1) {
            sv1 = other.sv1;
        }
        if (other.sv2 > sv2) {
            sv2 = other.sv2;
        }
        if (other.sv3 > sv3) {
            sv3 = other.sv3;
        }
    }
};

std::ostream &operator<<(std::ostream &os, const SupportQuality &q) {
    os << q.sv1 << "/" << q.sv2 << "/" << q.sv3;
    return os;
}

oreach::oreach(
    graph_access & G,
    ReachabilityConfig & config): condensation(G), config(config), queue_forward(G.number_of_nodes()), queue_backward(G.number_of_nodes()) {
}

oreach::~oreach() {
    delete bfs_query;
}

// sample l*n supports, choose 8 best using three criteria
void oreach::sampleSupports(bfs & backward_bfs) {

    std::size_t n = nodes.size();
    std::vector<char> blacklist(n, false);

    //int middle = vertices_per_component[max_component_number] * 0.5;


    unsigned int fw_degree = 0;
    unsigned int bw_degree = 0;

    struct SupportData {
        NodeID node;
        unsigned int fw_degree;
        unsigned int bw_degree;
        std::size_t n;

        SupportData(NodeID node, unsigned fw, unsigned bw, std::size_t n ): node(node), fw_degree(fw), bw_degree(bw), n(n) {
        }
    };

    unsigned int samples = n*config.l;
    unsigned int non_source_sink_samples = samples * 0.4;

    std::vector<SupportData> candidates;
    candidates.reserve(n);

    std::vector<NodeID> found_fw;
    std::vector<NodeID> found_bw;
    int nextNode;
    unsigned int min_nodes_per_component = vertices_per_component[max_component_number] * 0.5;

    for (unsigned int i = 0; i < samples; i++) {
        nextNode = random_functions::nextInt(0, n-1);
        // skip vertices from tiny connected components
        while (vertices_per_component[nodes[nextNode].component_number] < min_nodes_per_component || blacklist[nextNode]) {
            nextNode = random_functions::nextInt(0, n-1);
        }
        /*while (nodes[nextNode].component_number != max_component_number) {
            nextNode = random_functions::nextInt(0, n-1);
        }*/

        // skip sources/sinks for a portion of the candidates => best for first maximation condition
        while (i < non_source_sink_samples && (nodes[nextNode].isSource() || nodes[nextNode].isSink())) {
            nextNode = random_functions::nextInt(0, n-1);
        }

        blacklist[nextNode] = true;

        fw_degree = bfs_query->explore_neighbors(seen_forward, nextNode, touched_forward);
        bw_degree = backward_bfs.explore_neighbors(seen_backward, nextNode, touched_backward);
        candidates.emplace_back(SupportData(nextNode, fw_degree, bw_degree, n));

        cleanupBfs();
        nextNode = random_functions::nextInt(0, n-1);
    }

    auto comp1 = [] (const SupportData & l, const SupportData & r) {
        return ((l.bw_degree * l.fw_degree) > (r.bw_degree * r.fw_degree));
    };
    auto comp2 = [] (const SupportData & l, const SupportData & r) {
        return ((l.bw_degree * (l.n - l.bw_degree)) > (r.bw_degree * (r.n - r.bw_degree)));
    };
    auto comp3 = [] (const SupportData & l, const SupportData & r) {
        return ((l.fw_degree * (l.n - l.fw_degree)) > (r.fw_degree * (r.n - r.fw_degree)));
    };

    NodeID supports[8];
    auto it = candidates.begin();
    std::advance(it, 4);
    std::nth_element(candidates.begin(), it, candidates.end(), comp1);
    supports[0] = candidates[0].node;
    supports[1] = candidates[1].node;
    supports[2] = candidates[2].node;
    supports[3] = candidates[3].node;

    std::nth_element(candidates.begin(), it, candidates.end(), comp2);
    supports[4] = candidates[0].node;
    supports[5] = candidates[1].node;

    std::nth_element(candidates.begin(), it, candidates.end(), comp3);
    supports[6] = candidates[0].node;
    supports[7] = candidates[1].node;

    for (unsigned int i = 0; i < NUM_SUPPORTS; i++) {
        fw_degree = bfs_query->explore_neighbors(seen_forward, supports[i], touched_forward);
        bw_degree = backward_bfs.explore_neighbors(seen_backward, supports[i], touched_backward);
#ifdef RDEBUG
        std::cout << "node " << supports[i] << " with fw " << fw_degree << " and bw: " << bw_degree << " on level " << nodes[supports[i]].fw_level << "\n";;
#endif

        for (auto n : touched_forward) {
            nodes[n].fw += 1<<i;
        }
        for (auto n : touched_backward) {
            nodes[n].bw += 1<<i;
        }

        cleanupBfs();
    }
}

// sample l*n supports, choose 8 best using three criteria
void oreach::sampleSupportsFromLevels(bfs & backward_bfs) {
    unsigned int middleDeviation = max_ts_level / 5;
    unsigned int minLevel = max_ts_level / 2 - middleDeviation;
    unsigned int maxLevel = minLevel + 2*middleDeviation + max_ts_level%2;

    struct SupportData {
        NodeID node;
        unsigned long long fw_degree;
        unsigned long long bw_degree;

        SupportData(NodeID node, unsigned fw, unsigned bw)
            : node(node), fw_degree(fw), bw_degree(bw) {
        }
    };

    unsigned int max_samples = NUM_SUPPORTS*MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS;
    std::vector<SupportData> candidates;
    candidates.reserve(max_samples);

    countVerticesPerLevel<true,true>(vertices_per_level_fw);
    countVerticesPerLevel<true,false>(vertices_per_level_bw);

    #ifdef RDEBUG
    std::cout << "Searching in levels " << minLevel << " - " << maxLevel << "\n";
    std::cout << "Certainly considering vertices on forward or backward level of size < " << config.t << std::endl;
    #endif
    std::vector<NodeID> pool;
    unsigned int fw_degree = 0;
    unsigned int bw_degree = 0;
    // sample nodes of largest component from the middle levels range
    for (std::size_t i = 0; i < nodes.size() && candidates.size() < max_samples; i++) {
        if ( nodes[i].component_number >= large_components_end
                || nodes[i].isSink()
                || nodes[i].isSource()) {
            continue;
        }
        if (vertices_per_level_fw[nodes[i].fw_level] < config.t
                || vertices_per_level_bw[nodes[i].bw_level] < config.t) {
            fw_degree = bfs_query->explore_neighbors(seen_forward, i, touched_forward);
            bw_degree = backward_bfs.explore_neighbors(seen_backward, i, touched_backward);
            candidates.emplace_back(SupportData(i, fw_degree, bw_degree));
            cleanupBfs();
        } else if (nodes[i].fw_level >= minLevel && nodes[i].fw_level <= maxLevel) {
            pool.push_back(i);
        }
    }

    #ifdef RDEBUG
    std::cout << "Sampling " << max_samples << " samples, "
        << " pool size: "
        << (pool.size()+candidates.size()) << " (" << pool.size() << " random, " << candidates.size() << " fixed samples)\n";
    #endif

    if (candidates.size() < max_samples) {
        random_functions::permutate_vector_fast(pool, false);
        auto diff = std::min(max_samples - candidates.size(), pool.size());
        for (auto i = 0U; i < diff; i++) {
            fw_degree = bfs_query->explore_neighbors(seen_forward, pool[i], touched_forward);
            bw_degree = backward_bfs.explore_neighbors(seen_backward, pool[i], touched_backward);
            candidates.emplace_back(SupportData(pool[i], fw_degree, bw_degree));
            cleanupBfs();
        }
    }

    auto comp1 = [] (const SupportData & l, const SupportData & r) {
        return ((l.bw_degree * l.fw_degree) > (r.bw_degree * r.fw_degree));
    };


    //NodeID supports[8];
    auto it = candidates.begin();
    std::advance(it, NUM_SUPPORTS);
    std::nth_element(candidates.begin(), it, candidates.end(), comp1);


#ifdef RDEBUG
    std::cout << "Picked " << NUM_SUPPORTS << " support(s) out of " << candidates.size() << " candidates:\n";
#endif
    for (unsigned int i = 0; i < NUM_SUPPORTS; i++) {
        NodeID nextNode = candidates[i].node;
        fw_degree = bfs_query->explore_neighbors(seen_forward, nextNode, touched_forward);
        bw_degree = backward_bfs.explore_neighbors(seen_backward, nextNode, touched_backward);

        for (auto n : touched_forward) {
            nodes[n].fw += 1<<i;
        }
        for (auto n : touched_backward) {
            nodes[n].bw += 1<<i;
        }

        #ifdef RDEBUG
        std::cout << "Picked node " << nextNode << nodes[nextNode]
            << " with r+ = " << fw_degree << " and r- = " << bw_degree << "\n";
        #endif

        cleanupBfs();
    }
}

void oreach::sampleSupportsFromLevels2(bfs & backward_bfs) {

    unsigned int median_level_fw;
    unsigned int median_level_bw;
    unsigned int n_half = vertices_per_component[max_component_number] * 0.5;
    unsigned int sum = 0;
    unsigned int i = 0;
    bool useFwMedian = false;
    bool useBwMedian = false;

    std::vector<unsigned> vpl_fw;
    std::vector<unsigned> vpl_nss_fw;
    std::vector<unsigned> vpl_bw;
    std::vector<unsigned> vpl_nss_bw;
    countVerticesPerLevel<false,true>(vpl_fw);
    countVerticesPerLevel<true,true>(vpl_nss_fw);
    countVerticesPerLevel<false,false>(vpl_bw);
    countVerticesPerLevel<true,false>(vpl_nss_bw);

    while (sum < n_half) {
        sum += vpl_fw[i];
        i++;
    }
    median_level_fw = i;
    if (vpl_fw[median_level_fw] <= NUM_SUPPORTS) {
        useFwMedian = true;
    }
    sum = 0;
    i = 0;
    while (sum < n_half) {
        sum += vpl_bw[i];
        i++;
    }
    median_level_bw = i;
    if (vpl_bw[median_level_bw] <= NUM_SUPPORTS) {
        useBwMedian = true;
    }
#ifdef RDEBUG
    std::cout << "median level fw: " << median_level_fw << " bw: " << median_level_bw << "\n";
#endif

    std::vector<char> blacklist(condensation.number_of_nodes(), false);
    // unsigned int max_num_candidates = 10*NUM_SUPPORTS;
    std::vector<NodeID> fixed_pool;
    // sample nodes of largest component from the middle levels range
    for (std::size_t i = 0; i < nodes.size(); i++) {
        if ( nodes[i].component_number == max_component_number) {
            if ( !(nodes[i].isSink() || nodes[i].isSource())) {
                if ((vpl_nss_fw[nodes[i].fw_level] < config.t) || (vpl_nss_bw[nodes[i].bw_level] < config.t)) {
                    fixed_pool.push_back(i);
                    blacklist[i] = true;
                }
            } else if (useFwMedian && nodes[i].fw_level == median_level_fw) {
                fixed_pool.push_back(i);
                blacklist[i] = true;
            } else if (useBwMedian && nodes[i].bw_level == median_level_bw) {
                fixed_pool.push_back(i);
                blacklist[i] = true;
            }
        }
    }

    struct LevelOrder {
        NodeID n;
        unsigned int min_level;
        unsigned int fw_level;
        unsigned int bw_level;
        LevelOrder(NodeID n, unsigned int fw_level, unsigned int bw_level): n(n), fw_level(fw_level), bw_level(bw_level) {
            min_level = std::min(fw_level, bw_level);
        }
    };
    std::vector<LevelOrder> level_order_data;
    //if (fixed_pool.size() < 4*NUM_SUPPORTS) {       
        level_order_data.reserve(condensation.number_of_nodes());
        for (std::size_t i = 0; i < nodes.size(); i++) {
            if (nodes[i].component_number == max_component_number && !blacklist[i]) {
                level_order_data.emplace_back(LevelOrder(i, nodes[i].fw_level, nodes[i].bw_level));
            }
        }

        auto comp1 = [median_level_bw, median_level_fw] (const LevelOrder & l, const LevelOrder & r) {
            if (l.min_level == r.min_level) {
                return 
                (std::abs((int)l.fw_level - (int) median_level_fw) * std::abs((int)l.bw_level - (int) median_level_bw)) 
                <
                (std::abs((int)r.fw_level - (int) median_level_fw) * std::abs((int)r.bw_level - (int) median_level_bw)) ;
            }
            return l.min_level > r.min_level ;
        };
        std::sort(level_order_data.begin(), level_order_data.end(), comp1);

        // first top 4*K "middle" nodes into consideration
        for (unsigned int i = 0; i < NUM_SUPPORTS*4; i++) {
            fixed_pool.push_back(level_order_data[i].n);
        }
    // }

    std::size_t n = nodes.size();

    unsigned int fw_degree = 0;
    unsigned int bw_degree = 0;

    struct SupportData {
        NodeID node;
        unsigned int fw_degree;
        unsigned int bw_degree;
        std::size_t n;

        SupportData(NodeID node, unsigned fw, unsigned bw, std::size_t n ): node(node), fw_degree(fw), bw_degree(bw), n(n) {
        }
    };

    std::vector<SupportData> candidates;
    candidates.reserve(condensation.number_of_nodes());

    std::vector<NodeID> found_fw;
    std::vector<NodeID> found_bw;
    NodeID nextNode;
#ifdef RDEBUG
    std::cout << " got " << fixed_pool.size() << " candidates\n";
#endif
    for (unsigned int i = 0; i < fixed_pool.size(); i++) {
        nextNode = fixed_pool[i];
        fw_degree = bfs_query->explore_neighbors(seen_forward, nextNode, touched_forward);
        bw_degree = backward_bfs.explore_neighbors(seen_backward, nextNode, touched_backward);
        candidates.emplace_back(SupportData(nextNode, fw_degree, bw_degree, n));

        cleanupBfs();
    }

    auto comp_support = [] (const SupportData & l, const SupportData & r) {
        return ((l.bw_degree * l.fw_degree) > (r.bw_degree * r.fw_degree));
    };


    auto it = candidates.begin();
    std::advance(it, NUM_SUPPORTS);
    std::nth_element(candidates.begin(), it, candidates.end(), comp_support);


    for (unsigned int i = 0; i < NUM_SUPPORTS; i++) {
        nextNode = candidates[i].node;
        fw_degree = bfs_query->explore_neighbors(seen_forward, nextNode, touched_forward);
        bw_degree = backward_bfs.explore_neighbors(seen_backward, nextNode, touched_backward);

        #ifdef RDEBUG
        std::cout << "node " << nextNode << " with fw " << fw_degree << " and bw: " << bw_degree << " on level " << nodes[nextNode].fw_level << "\n";
        #endif


        for (auto n : touched_forward) {
            nodes[n].fw += 1<<i;
        }
        for (auto n : touched_backward) {
            nodes[n].bw += 1<<i;
        }

        cleanupBfs();
    }

}

void oreach::supportingVertexStep(
        unsigned int round,
        bfs & backward_bfs, std::pair<NodeID, NodeID> range,
        std::vector<NodeID> & sortedNodes
    ) {

    std::size_t middle = (range.second - range.first) * 0.5 + range.first;

#ifdef RDEBUG
    unsigned int max_node = sortedNodes[middle];
#endif

    unsigned int fw_degree = 0;
    unsigned int bw_degree = 0;
    unsigned int closeness = 0;

    std::vector<NodeID> found_fw;
    std::vector<NodeID> found_bw;

    unsigned int max_fw = bfs_query->explore_neighbors(seen_forward, sortedNodes[middle], touched_forward);
    unsigned int max_bw = backward_bfs.explore_neighbors(seen_backward, sortedNodes[middle], touched_backward);

    // assign the touched nodes for possible further use later
    // when updating the member datastructure for supports
    found_fw = touched_forward;
    found_bw = touched_backward;

    cleanupBfs();
    unsigned int max_closeness = max_fw*max_bw;

    for (unsigned i = range.first; i < range.second; i++) {
        NodeID node = sortedNodes[i];

        fw_degree = bfs_query->explore_neighbors(seen_forward, node, touched_forward);
        bw_degree = backward_bfs.explore_neighbors(seen_backward, node, touched_backward);
        closeness = fw_degree * bw_degree;

        if (closeness > max_closeness) {
#ifdef RDEBUG
            max_node = node;
#endif
            found_fw = touched_forward;
            found_bw = touched_backward;
            max_fw = fw_degree;
            max_bw = bw_degree;
            max_closeness = closeness;
        }

        cleanupBfs();
    }

    for (auto n : found_fw) {
        nodes[n].fw += 1<<round;
    }
    for (auto n : found_bw) {
        nodes[n].bw += 1<<round;
    }


#ifdef RDEBUG
    std::cout << "max node: " << max_node << " with forward size " << max_fw << " backward size: " << max_bw << " closeness: " << max_closeness << "\n";
#endif
}

void oreach::pickSupports(bfs &backward_bfs) {
    const auto middle_level = max_ts_level / 2;
    const auto min_level = (max_ts_level + 1) / SUPPORT_MIN_LEVEL_RECIPROCAL;
    const auto num_stretched = max_ts_level - 2*min_level + 1 - (max_ts_level % 2);
#ifdef RDEBUG
    std::cout << "\nPicking supportive vertices around middle level " << middle_level << std::endl;
#endif
    std::vector<std::vector<NodeID>> stretched_levels(num_stretched);
    forall_nodes(condensation, n) {
        const auto node = nodes[n];
        if (node.isIsolated()
                || node.isSource()
                || node.isSink()
                || node.component_number >= large_components_end) {
            continue;
        }
        if (node.fw_level <= node.bw_level && node.fw_level >= min_level) {
            const auto index = middle_level == node.fw_level ? 0 : 2 * (middle_level - node.fw_level - 1) + 1;
            stretched_levels[index].push_back(n);
        } else if (node.bw_level < node.fw_level && node.bw_level >= min_level) {
            const auto index = middle_level == node.bw_level ? 0 : 2 * (middle_level - node.bw_level - 1) + 2;
            stretched_levels[index].push_back(n);
        } else {
            continue;
        }
    } endfor

#ifdef RDEBUG
    const auto num_print = 3UL;
    std::cout << "\nFirst " << num_print << " stretched levels:\n";
    for (auto i = 0U; i < std::min(stretched_levels.size(), num_print); i++) {
        std::cout << i << " (" << stretched_levels[i].size()  << "):";
        for (auto &n : stretched_levels[i]) {
            std::cout << " " << n << " (" << nodes[n].fw_level << "/" << nodes[n].bw_level << "),";
        }
        std::cout << std::endl;
    }
#endif

    std::stable_sort(stretched_levels.begin(), stretched_levels.end(),
            [](const std::vector<NodeID> &lhs, const std::vector<NodeID> &rhs) {
                return lhs.size() < rhs.size(); });
#ifdef RDEBUG
    std::cout << "\nFirst " << num_print << " sorted stretched levels:\n";
    for (auto i = 0U; i < std::min(stretched_levels.size(), num_print); i++) {
        std::cout << i << " (" << stretched_levels[i].size()  << "):";
        for (auto &n : stretched_levels[i]) {
            std::cout << " " << n << " (" << nodes[n].fw_level << "/" << nodes[n].bw_level << "),";
        }
        std::cout << std::endl;
    }
#endif

    std::vector<unsigned> vpl;
    std::vector<unsigned> fw_level_size_acc;
    std::vector<unsigned> bw_level_size_acc;
    fw_level_size_acc.assign(max_ts_level+1, 0);
    bw_level_size_acc.assign(max_ts_level+1, 0);

    countVerticesPerLevel<false,true>(vpl);
    std::partial_sum(vpl.cbegin(), vpl.cend(), fw_level_size_acc.begin());

    countVerticesPerLevel<false,false>(vpl);
    std::partial_sum(vpl.cbegin(), vpl.cend(), bw_level_size_acc.begin());

    std::vector<NodeID> supports;
    std::vector<SupportQuality> support_quality;
    SupportQuality best_quality;

    const auto num_supports = NUM_SUPPORTS;
    auto num_candidates = 0U;
    const auto max_candidates = MAX_SUPPORT_CANDIDATE_FACTOR * num_supports;

    for (const auto &level : stretched_levels) {
        if (num_candidates >= max_candidates) {
            break;
        }
        for (auto &candidate : level) {
#ifdef VERBOSE
            std::cout << "Checking support candidate " << candidate << nodes[candidate]
                << std::endl;
#endif

            if (supports.size() >= num_supports
                    && nodes[candidate].fw + nodes[candidate].bw == 0) {
#ifdef VERBOSE
                std::cout << "Can't replace earlier support and we already have enough."
                    << std::endl;
#endif
                continue;
            }

            const unsigned long out_reach = bfs_query->explore_neighbors(
                    seen_forward, candidate, touched_forward);
            const unsigned long in_reach = backward_bfs.explore_neighbors(
                    seen_backward, candidate, touched_backward);

            const auto max_out = bw_level_size_acc[nodes[candidate].bw_level - 1] + 1;
            const auto max_in = fw_level_size_acc[nodes[candidate].fw_level - 1] + 1;
            const auto sv1_quality = out_reach * in_reach;
            const auto sv2_quality = (max_out - out_reach) * out_reach;
            const auto sv3_quality = (max_in - in_reach) * in_reach;
            const auto quality = SupportQuality(sv1_quality, sv2_quality, sv3_quality);
#ifdef VERBOSE
            std::cout << "Candidate has r+=" << out_reach << ", r-=" << in_reach
                << ", max r+=" << max_out
                << ", max r-=" << max_in
                << ", r!+=" << max_out - out_reach
                << ", r!-=" << max_in - in_reach << " => " << quality << std::endl;
#endif
            num_candidates++;

            auto worst_quality = quality;
            auto worst_support = supports.size(); // that's me

            if (nodes[candidate].fw + nodes[candidate].bw > 0) {
                // the following two are only relevant if no
                // additional supports can be picked
                bool bestInOne = !(quality < best_quality);
                auto partially_worse = worst_support;

                for (auto s = 0U; s < supports.size(); s++) {
                    if ((nodes[candidate].fw & (1 << s)) > 0) {
                        if (support_quality[s] < worst_quality) {
                            // s-th support can reach candidate and is worse
                            worst_quality = support_quality[s];
                            worst_support = s;
                        } else if (bestInOne
                                && supports.size() == num_supports
                                && (support_quality[s] < best_quality)) {
                            partially_worse = s;
                        }
                    } else if ((nodes[candidate].bw & (1 << s)) > 0) {
                        if (support_quality[s] < worst_quality) {
                            // s-th support is reachable from candidate and this one is worse
                            worst_quality = support_quality[s];
                            worst_support = s;
                        } else if (bestInOne
                                && supports.size() == num_supports
                                && bestInOne
                                && (support_quality[s] < best_quality)) {
                            partially_worse = s;
                        }
                    }
                }
                if (partially_worse != supports.size()) {
                    worst_support = partially_worse;
#ifdef VERBOSE
                    std::cout << "Intersecting earlier support "
                        << supports[worst_support] << nodes[supports[worst_support]]
                        << " with quality " << support_quality[worst_support]
                        << " is partially worse and gets replaced." << std::endl;
#endif
                    clearSupport(worst_support);
                } else if (worst_support < supports.size()) {
#ifdef VERBOSE
                    std::cout << "Intersecting earlier support "
                        << supports[worst_support] << nodes[supports[worst_support]]
                        << " with quality " << support_quality[worst_support]
                        << " is worse and gets replaced." << std::endl;
#endif
                    clearSupport(worst_support);
                } else if (bestInOne && supports.size() < num_supports) {
#ifdef VERBOSE
                    std::cout<< "None worse, but candidate is best for at least one criterion. "
                        "Picking it additionally." << std::endl;
#endif
                } else {
#ifdef VERBOSE
                    std::cout << "No intersecting earlier support is worse. Discarding candidate."
                        << std::endl;
#endif
                    cleanupBfs();
                    continue;
                }
            }

            best_quality.update(quality);
#ifdef VERBOSE
            std::cout << "Candidate is picked." << std::endl;
#endif

            setSupport(worst_support, touched_forward, touched_backward);
            if (worst_support < supports.size()) {
                supports[worst_support] = candidate;
                support_quality[worst_support] = quality;
            } else {
                supports.push_back(candidate);
                support_quality.push_back(quality);
            }
            cleanupBfs();
        }
    }

#ifdef RDEBUG
    std::cout << "Picked " << supports.size() << " support(s) out of "
        << num_candidates << " candidates:";
    for (auto &p : supports) {
        std::cout << "\n" << p << nodes[p];
    }
    std::cout << std::endl;
#endif
}

void oreach::chooseSupportingVertex() {
    unsigned int n = condensation.number_of_nodes();
    bfs backward_bfs(condensation_r);
    seen_forward.assign(n, false);
    seen_backward.assign(n, false);
    unsigned int comp_n = vertices_per_component[max_component_number];
    touched_forward.reserve(comp_n);
    touched_backward.reserve(comp_n);

    if (config.support == 1) {
        // sort vertices of maximum component
        // ie iterate through topological order, add a vertex to max_component_top_sort
        // when vertex is part of maximum component
        std::vector<unsigned int> max_component_top_sort(comp_n);
        unsigned int count = 0;
        for (unsigned i = 0; i < n; i++) {
            NodeID node = sorted_nodes[i];
            if (nodes[node].component_number == max_component_number){
                max_component_top_sort[count++] = node;
            }
        }
        unsigned int vertices_per_range = (max_component_top_sort.size() / num_of_supports);
        #ifdef RDEBUG
        std::cout << "for each support testing " << (2*vertices_per_range*0.15) << " nodes\n";
        #endif
        for (std::size_t i = 0; i < num_of_supports; i++) {
            NodeID middle =  vertices_per_range * (0.5 + i);
            NodeID deviation = vertices_per_range * 0.15;

            std::pair<NodeID, NodeID> range = {middle-deviation, middle+deviation};
            // std::cout << "range: " << range.first << " to " << range.second << "\n";
            supportingVertexStep(i,backward_bfs, range, max_component_top_sort);
        }

    } else if (config.support == 2) {
        #ifdef RDEBUG
        std::cout << "SUPPORT CHOICE VARIANT 2 - random sampling\n"
            << "sampling " << n*config.l << " nodes for the supports \n";
        #endif
        sampleSupports(backward_bfs);
    } else if (config.support == 3) {
        #ifdef RDEBUG
        std::cout << "SUPPORT CHOICE VARIANT 3 - sampling from levels\n"
            << "Sampling at most " << NUM_SUPPORTS*MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS << " nodes for the supports\n";
        #endif
        sampleSupportsFromLevels(backward_bfs);
    } else if (config.support == 4) {
        #ifdef RDEBUG
        std::cout << "SUPPORT CHOICE VARIANT 4 - sampling from median levels\n";
        #endif
        sampleSupportsFromLevels2(backward_bfs);
    } else if (config.support == 5) {
        #ifdef RDEBUG
        std::cout << "SUPPORT CHOICE VARIANT 5 - sampling from middle small levels\n";
        #endif
        pickSupports(backward_bfs);
    }

}

void oreach::initialize() {
    unsigned int scc_n = condensation.number_of_nodes();

    nodes.assign(scc_n, NodeInfo());

    // reverse graph of comp graph
    condensation.buildReverseGraph(condensation_r);

    // initialize BFS classes
    bfs_query = new bfs(condensation);

    // determine components of the (condensation) graph
    determineComponents();

    doTopologicalSorts();

    // choose supporting vertice
    chooseSupportingVertex();

    #ifdef RDEBUG
        std::cout << "n,m: " << condensation.number_of_nodes() << ", " << condensation.number_of_edges() << "\n";
    #endif

    #ifdef FREE_SPACE
    source_nodes.clear();
    sink_nodes.clear();
    std::vector<NodeID>(0, 0).swap(source_nodes);
    std::vector<NodeID>(0, 0).swap(sink_nodes);
    // todo rest
    #endif

    initialized = true;
}

void oreach::doTopologicalSorts() {
        auto n = condensation.number_of_nodes();

        topological_sort top_sort;

        // topological levels
        std::vector<unsigned int> fw_levels(n, 0);
        std::vector<unsigned int> bw_levels(n, 0);
        top_sort.levels(condensation, condensation_r, fw_levels, source_nodes);
        top_sort.levels(condensation_r, condensation, bw_levels, sink_nodes);

        for (std::size_t i = 0; i < fw_levels.size(); i++) {
            nodes[i].fw_level = fw_levels[i];
            nodes[i].bw_level = bw_levels[i];
            if (nodes[i].component_number < large_components_end
                    && fw_levels[i] > max_ts_level) {
                max_ts_level = fw_levels[i];
            }
        }

        #ifdef RDEBUG
        std::cout << "Maximum topsort level (= path length) in large components: "
            << max_ts_level << "\n";
        #endif
        #ifdef VERBOSE
        std::vector<unsigned> vpl;
        std::vector<unsigned> vpl_nosourcesink;
        countVerticesPerLevel<false,true>(vpl);
        countVerticesPerLevel<true,true>(vpl_nosourcesink);
        std::cout << "Vertices per forward level (without sources and sinks)\n";
        for (unsigned i = 0; i <= max_ts_level; i++) {
            std::cout << "Level " << i << " : " << vpl[i]
                << " ("<< vpl_nosourcesink[i] << ")\n";
        }
        countVerticesPerLevel<false,false>(vpl);
        countVerticesPerLevel<true,false>(vpl_nosourcesink);
        std::cout << "\nVertices per backward level (without sources and sinks)\n";
        for (unsigned i = 0; i <= max_ts_level; i++) {
            std::cout << "Level " << i << " : " << vpl[i]
                << " ("<< vpl_nosourcesink[i] << ")\n";
        }
        #endif

#ifdef RESTRICT_TO_LARGE
#ifdef RDEBUG
        auto num_sources = source_nodes.size();
        auto num_sinks = sink_nodes.size();
#endif
        filterVerticesInLargeComponents(source_nodes);
        filterVerticesInLargeComponents(sink_nodes);
#ifdef RDEBUG
        std::cout << "Filtered out " << num_sources - source_nodes.size()
            << " sources and " << num_sinks - sink_nodes.size()
            << " sinks in small components." << std::endl;
#endif
#endif

        // extended topological orderings
        sorted_nodes.reserve(n);
        std::vector<unsigned int> high_indices (n, 0);
        std::vector<unsigned int> max_indices (n, 0);

        // topological sorting starting with sources and sinks
        for (unsigned i = NUM_TOPSORTS/2; i < NUM_TOPSORTS; i++) {
            sorted_nodes.clear();
            top_sort.extendedSortStartingWith(
                    condensation_r, sink_nodes, sorted_nodes, high_indices, max_indices, false);
            assert(sorted_nodes.size() == n_nonisolated);
            for (std::size_t j = 0; j < sorted_nodes.size(); j++) {
                nodes[sorted_nodes[j]].pos[i] = j;
            }
            for (auto &n : sorted_nodes) {
                nodes[n].low[i-(NUM_TOPSORTS)/2] = high_indices[n];
                nodes[n].min[i-(NUM_TOPSORTS)/2] = max_indices[n];
            }
        }

        for (unsigned i = 0; i < NUM_TOPSORTS/2; i++) {
            sorted_nodes.clear();
            top_sort.extendedSortStartingWith(condensation, source_nodes, sorted_nodes,
                    high_indices, max_indices);
            assert(sorted_nodes.size() == n_nonisolated);
            for (std::size_t j = 0; j < sorted_nodes.size(); j++) {
                nodes[sorted_nodes[j]].pos[i] = j;
            }
            for (auto &n : sorted_nodes) {
                nodes[n].high[i] = high_indices[n];
                nodes[n].max[i] = max_indices[n];
            }
        }
}

void oreach::determineComponents() {
    auto num_sources = 0U;
    auto num_sinks = 0U;
    auto num_isolated = 0U;

    union_find uf(condensation.number_of_nodes());
    forall_nodes(condensation, node) {
        if (condensation.getNodeDegree(node) == 0) {
            if (condensation_r.getNodeDegree(node) == 0) {
              nodes[node].setIsolated();
              num_isolated++;
              continue;
            }
            nodes[node].setSink();
            num_sinks++;
        } else if (condensation_r.getNodeDegree(node) == 0) {
            nodes[node].setSource();
            num_sources++;
        }
        forall_out_edges(condensation, e, node) {
            NodeID target = condensation.getEdgeTarget(e);
            uf.Union(node, target);
        } endfor
    } endfor

    unsigned int total_components = uf.n() - num_isolated;
    source_nodes.reserve(num_sources);
    sink_nodes.reserve(num_sinks);
    n_nonisolated = condensation.number_of_nodes() - num_isolated;

    if (total_components == 1) {
        #ifdef RDEBUG
            std::cout << "Just one component. Fast track." << std::endl;
        #endif
        forall_nodes(condensation, node) {
            if (nodes[node].isIsolated()) {
                continue;
            }
            nodes[node].component_number = 0;
            if (nodes[node].isSource()) {
              source_nodes.push_back(node);
            }
            if (nodes[node].isSink()) {
              sink_nodes.push_back(node);
            }
        } endfor
        vertices_per_component.push_back(n_nonisolated);
        large_components_end = 1;
        max_component_number = 0;
        return;
    }

    unsigned int counter = 0;
    vertices_per_component.assign(total_components, 0);

    // map large numbers to consecutive component numbers
    std::vector<unsigned> transl(condensation.number_of_nodes(),
            std::numeric_limits<unsigned>::max());

#ifdef ORDERED_SRC_SNK
#ifdef RDEBUG
    std::cout << "Ordering of sources and sinks is enabled." << std::endl;
#endif
    std::vector<NodeID> sources;
    sources.reserve(num_sources);
    std::vector<NodeID> sinks;
    sinks.reserve(num_sinks);
    std::vector<unsigned int> sources_per_component;
    sources_per_component.assign(total_components, 0);
    std::vector<unsigned int> sinks_per_component;
    sinks_per_component.assign(total_components, 0);
#endif

    unsigned num_large_components = 0;
    auto pivot = condensation.number_of_nodes() / SMALL_COMPONENT_RECIPROCAL ;

    forall_nodes(condensation, node) {
        if (nodes[node].isIsolated()) {
          continue;
        }
        unsigned int num = uf.Find(node);
        if (transl[num] == std::numeric_limits<unsigned>::max()) {
            transl[num] = counter++;
        }
        nodes[node].component_number = transl[num];
        if (nodes[node].isSource()) {
#ifdef ORDERED_SRC_SNK
            sources.push_back(node);
            sources_per_component[transl[num]]++;
#else
            source_nodes.push_back(node);
#endif
        }
        if (nodes[node].isSink()) {
#ifdef ORDERED_SRC_SNK
            sinks.push_back(node);
            sinks_per_component[transl[num]]++;
#else
            sink_nodes.push_back(node);
#endif
        }
        vertices_per_component[transl[num]]++;
        if (vertices_per_component[transl[num]] == pivot) {
          num_large_components++;
        }
    } endfor


    #ifdef RDEBUG
        auto minmax_size = std::minmax_element(
                vertices_per_component.cbegin(), vertices_per_component.cend());
        auto sum = std::accumulate(
                vertices_per_component.cbegin(), vertices_per_component.cend(), 0U);
        std::vector<unsigned int> sorted_components(vertices_per_component);
        std::sort(sorted_components.begin(), sorted_components.end(),
                std::greater<unsigned>());

        std::cout << "\n============== Graph Statistics =============="
                  << "\n #components (non-isolated): " << total_components
                  << "\n   threshold for large     : "
                  << "n/"<< SMALL_COMPONENT_RECIPROCAL << " = " << pivot
                  << "\n   #large                  : " << num_large_components
                  << "\n   min size                : " << *minmax_size.first
                  << "\n   max size                : " << *minmax_size.second
                  << "\n   avg size                : " << sum/vertices_per_component.size()
                  << "\n   first 5 sizes           : ";
        for (std::size_t i = 0; i < std::min((std::size_t)5, sorted_components.size()); i++) {
            std::cout << sorted_components[i] << " ";
        }
        std::cout << "\n #sources                  : " << num_sources
                  << "\n #sinks                    : " << num_sinks
                  << "\n #isolated                 : " << num_isolated
                  << "\n #non-isolated             : " << n_nonisolated
                  << "\n==============================================\n" << std::endl;
    #endif

#ifdef ORDERED_SRC_SNK
    // group sources and sinks by component
    // not stable currently: original order is reversed
#ifdef RDEBUG
    std::cout << "Grouping sources and sinks by component..." << std::endl;
#endif
    source_nodes.assign(num_sources, std::numeric_limits<unsigned int>::max());
    sink_nodes.assign(num_sinks, std::numeric_limits<unsigned int>::max());
    std::partial_sum(sources_per_component.begin(), sources_per_component.end(),
            sources_per_component.begin());
    std::partial_sum(sinks_per_component.begin(), sinks_per_component.end(),
            sinks_per_component.begin());
    assert(sources_per_component.back() == num_sources);
    assert(sinks_per_component.back() == num_sinks);
    for (auto &node : sources) {
        auto cc = nodes[node].component_number;
        assert(source_nodes[sources_per_component[cc] - 1]
                == std::numeric_limits<unsigned int>::max());
        source_nodes[sources_per_component[cc] - 1] = node;
        sources_per_component[cc]--;
    }
    for (auto &node : sinks) {
        auto cc = nodes[node].component_number;
        assert(sink_nodes[sinks_per_component[cc] - 1]
                == std::numeric_limits<unsigned int>::max());
        sink_nodes[sinks_per_component[cc] - 1] = node;
        sinks_per_component[cc]--;
    }
#endif

    max_component_number = std::distance(vertices_per_component.cbegin(),
        std::max_element(vertices_per_component.cbegin(), vertices_per_component.cend()));

    if (num_large_components == 0
            || (num_large_components == 1 && max_component_number == 0)) {
      #ifdef RDEBUG
            std::cout << "At most one large component, which is first. Nothing more to do."
                << std::endl;
      #endif
      large_components_end = num_large_components;
      return;
    }

    #ifdef RDEBUG
        std::cout
            << "Reordering components. Largest component currently has number "
            << max_component_number << "." << std::endl;
    #endif

    std::vector<unsigned int> cc_remap(total_components, 0);
    std::iota(cc_remap.begin(), cc_remap.end(), 0);

    if (max_component_number != 0) {
      cc_remap[max_component_number] = 0;
      cc_remap[0] = max_component_number;
      std::swap(vertices_per_component[0], vertices_per_component[max_component_number]);
      max_component_number = 0;
    }

    large_components_end = 1;

    for (auto i = total_components - 1;
            i > large_components_end
            && large_components_end < num_large_components;) {
        while (vertices_per_component[large_components_end] >= pivot
                && i > large_components_end) {
            large_components_end++;
        }
        if (large_components_end >= num_large_components) {
          break;
        }
        while (vertices_per_component[i] < pivot
                && i > large_components_end) {
            i--;
        }
        if (i > large_components_end) {
            cc_remap[large_components_end] = i;
            cc_remap[i] = large_components_end;
            std::swap(vertices_per_component[i],
                    vertices_per_component[large_components_end]);
            large_components_end++;
            i--;
        }
    }
    if (vertices_per_component[large_components_end] >= pivot) {
        large_components_end++;
    }

#ifdef RDEBUG
    std::cout << "Large components : [0, " << large_components_end << ")" << std::endl;
#endif

    // update component numbers
    forall_nodes(condensation, node) {
        if (!nodes[node].isIsolated()) {
            nodes[node].component_number = cc_remap[nodes[node].component_number];
        }
    } endfor
}

void oreach::debug() {
    #ifdef RDEBUG
    auto max = std::max(*std::max_element(std::begin(query_result_stats),
                std::end(query_result_stats)), num_pruned_search_steps);
    auto width = static_cast<unsigned int>(log10(max)) + 1U;
    auto sum = std::accumulate(std::begin(query_result_stats),
            std::end(query_result_stats), 0);
    auto sum_negative = std::accumulate(std::begin(query_result_stats),
            &query_result_stats[Query::castResult(Query::Result::POSITIVE)], 0);
    auto sum_positive = std::accumulate(&query_result_stats[Query::castResult(
                Query::Result::POSITIVE)], std::end(query_result_stats), 0);
    auto print = [&](auto description, Query::Result r, auto among) {
        auto value = query_result_stats[Query::castResult(r)];
        std::cout << description << std::setw(width) <<  value << " ("
          << std::fixed << std::setw(6) << std::setprecision(2)
          << 100.0 * value / sum << "% / "
          << std::fixed << std::setw(6) << std::setprecision(2)
          << 100.0 * value / among << "%)\n";
    };
    std::cout << "\n\n#### Query Statistics ####\n\n"
      << "  positive: " << sum_positive << " ("
      << std::fixed << std::setw(6) << std::setprecision(2)
      << 100.0 * sum_positive / sum << "%), "
      << "negative: " << sum_negative << " ("
      << std::fixed << std::setw(6) << std::setprecision(2)
      << 100.0 * sum_negative / sum << "%),  "
      << "total: " << sum << "\n\n"

      << "  negative (%all, %negative)\n";
    print("    different component           : ", Query::Result::DIFF_WCC, sum_negative);
    print("    t source or s sink            : ", Query::Result::SRC_SNK, sum_negative);
    print("    contradicting forward levels  : ", Query::Result::TS_FW_LEVELS, sum_negative);
    print("    contradicting backward levels : ", Query::Result::TS_BW_LEVELS, sum_negative);
    print("    contradicting topsort numbers : ", Query::Result::TS_NUM, sum_negative);
    print("    topsort max range             : ", Query::Result::TS_MAX, sum_negative);
    print("    topsort min range             : ", Query::Result::TS_MIN, sum_negative);
    print("    negative support (forward)    : ", Query::Result::SUPPORT_NEG_FW, sum_negative);
    print("    negative support (backward)   : ", Query::Result::SUPPORT_NEG_BW, sum_negative);
    print("    negative fallback             : ", Query::Result::FALLBACK_NEG, sum_negative);
    if (fallback == FB_PREACH) {
        print("    PReaCH levels                 : ", Query::Result::PREACH_LEVELS, sum_negative);
        print("    PReaCH min forward            : ", Query::Result::PREACH_MIN_FWD, sum_negative);
        print("    PReaCH min backward           : ", Query::Result::PREACH_MIN_BWD, sum_negative);
        print("    PReaCH max backward           : ", Query::Result::PREACH_MAX_BWD, sum_negative);
        print("    PReaCH max forward            : ", Query::Result::PREACH_MAX_FWD, sum_negative);
        print("    PReaCH negative fallback      : ", Query::Result::PREACH_FALLBACK_NEG, sum_negative);
    }

    std::cout  << "\n  positive (%all, %positive)\n";
    print("    same vertex                   : ", Query::Result::SAME_VERTEX, sum_positive);
    print("    same SCC                      : ", Query::Result::SAME_SCC, sum_positive);
    print("    topsort high range            : ", Query::Result::TS_HI, sum_positive);
    print("    topsort low range             : ", Query::Result::TS_LO, sum_positive);
    print("    positive support              : ", Query::Result::SUPPORT_POS, sum_positive);
    print("    topsort max is t              : ", Query::Result::TS_EQMAX, sum_positive);
    print("    topsort min is s              : ", Query::Result::TS_EQMIN, sum_positive);
    print("    edge                          : ", Query::Result::EDGE, sum_positive);
    print("    positive fallback             : ", Query::Result::FALLBACK_POS, sum_positive);

    if (fallback == FB_PREACH) {
        print("    PReaCH same vertex            : ", Query::Result::PREACH_SAME_VERTEX, sum_positive);
        print("    PReaCH forward range          : ", Query::Result::PREACH_RANGE_FWD, sum_positive);
        print("    PReaCH backward range         : ", Query::Result::PREACH_RANGE_BWD, sum_positive);
        print("    PReaCH positive fallback      : ", Query::Result::PREACH_FALLBACK_POS, sum_positive);
    }

    print("\n  Unanswered                      : ", Query::Result::UNKNOWN, sum_negative + sum_positive);
    std::cout <<
        "\n    #Pruned search steps          : " << std::setw(width) << num_pruned_search_steps
        << std::endl;
    #endif
}

