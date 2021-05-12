/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Jonathan Trummer
 *      Kathrin Hanauer
 *****************************************************************************/


#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include "extern/KaHIP/lib/tools/chronotimer.h"
#include "extern/KaHIP/lib/io/graph_io.h"
#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/definitions.h"
#include "parse_reachability_parameters.h"
#include "app/reachabilityconfig.h"
#include "lib/reachability/oreach.h"
#include "lib/algorithms/bfs.h"
#include "extern/KaHIP/lib/tools/random_functions.h"
#include "extern/KaHIP/lib/algorithms/strongly_connected_components.h"
#include <algorithm>
#include "lib/reachability/matrix/full_reach.h"

#include "query_generator.h"

#include <time.h>

using namespace std;

bool file_exists(const std::string &filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}


int main(int argn, char **argv) {
    bool runSReach{true};

    ReachabilityConfig config;
    std::string graph_filename;

    auto ret_code = parse_reachability_parameters(argn, argv, graph_filename, config);
    if(ret_code) {
        return ret_code;
    }

    graph_access input_graph;
    ChronoTimer t;
    t.restart();

    if (graph_io::readDiGraph(input_graph, graph_filename)) {
        return 1;
    }


    std::cout << "input IO took " << t.elapsed<std::milli>() << " ms\n"
              << "n(input): " << input_graph.number_of_nodes() << "\n"
              << "m(input): " << input_graph.number_of_edges() << std::endl;


    random_functions::setSeed(config.seed);


    if (config.writeCondensationGraph) {
        graph_access condensation;
        strongly_connected_components scc;
        std::vector<int> component_assignment(input_graph.number_of_nodes());
        scc.build_scc_graph(input_graph, component_assignment, condensation);
        graph_io::writeGraph(condensation, config.outputCondensation, false);
        std::cout << "done writing condensation to " << config.outputCondensation << "\n";
        graph_filename = config.outputCondensation;
    }

     double time_scc_init = 0;

    std::ifstream inputgraphstream;
    inputgraphstream.open(graph_filename);
    if (!inputgraphstream.is_open()) {
        std::cout << "ERROR OPENING GRAPH FILE" << graph_filename << "\n";
        exit(1);
    }

    oreach r(input_graph, config);
    t.restart();
    r.initialize();
    time_scc_init = t.elapsed<std::milli>();
    std::cout << "init took " << time_scc_init << " ms\n";
    FullReach fullreach(input_graph);


    std::vector<std::vector<Query>> all_queries;
    std::vector<std::string> query_types;
    int N = 100000;

    std::vector<Query> queries;
    queries.reserve(N);

    if (config.writeOutputQueryFile) {
        QueryGenerator generator(config, input_graph, 100000, r);
        generator.generate_negative();
        generator.generate_random();
        generator.generate_positive2();
    }



    if (config.useQueryFile) {
        for (auto inputQueryFile : config.queryFiles) {
            ifstream myfile (inputQueryFile);
            if (!myfile.is_open()) {
                std::cout << "failed to open " << inputQueryFile << '\n';
                exit(1);
            }
            // convention:
            // - positive queries: *.query.t.txt
            // - negative queries: *.query.f.txt
            // - random queries : *.query.txt
            // - mixed the positive and negative queries: *.query.m.txt
            if (inputQueryFile.find(".t.txt") != std::string::npos) {
                query_types.push_back("true");
            } else if (inputQueryFile.find(".f.txt") != std::string::npos) {
                query_types.push_back("false");
            } else if (inputQueryFile.find(".m.txt") != std::string::npos) {
                query_types.push_back("mixed");
            } else {
                query_types.push_back("random");
            }

            std::cout << "reading queries from " << inputQueryFile << "\n";

            queries.clear();

            NodeID source, target;
            bool success = false;
            string line;
            while (getline(myfile, line)) {
                std::stringstream ss(line);
                ss >> source;
                ss >> target;
                ss >> success;
                queries.push_back(Query(source, target));
            }
            all_queries.push_back(queries);
        }
    } else {
        QueryGenerator generator(config, input_graph, 100000, r);
        generator.setWriteOutput(false);
        generator.generate_random();
        queries = generator.getQueries();
        all_queries.push_back(queries);
        query_types.push_back("random");
    }
    std::vector<FallbackMode> fallbacks;

    if (config.pruning_bibfs) {
        fallbacks.push_back(FB_PIBFS);
    }
    if (config.pruning_dfs) {
        fallbacks.push_back(FB_PDFS);
    }
    if (fallbacks.empty()) {
        fallbacks.push_back(FB_NONE);
    }


    double time_full_reach_init = 0;
    if (config.fullReach) {
        t.restart();
        fullreach.init();
        time_full_reach_init += t.elapsed();
        r.reset();
    }


    for (std::size_t q = 0; q < all_queries.size(); q++) {
        double time_scc_query_pbibfs = 0;
        double time_scc_query_pdfs = 0;
        double time_cur_query = 0;
        double time_fullreach_query = 0;
        auto num_queries = all_queries[q].size();

        std::cout << "\n======= " << query_types[q] << " ========\n";

        if (config.fullReach) {
            std::cout << "running full reach matrix\n";
            t.restart();
            fullreach.query(all_queries[q]);
            time_fullreach_query = t.elapsed<std::micro>();

            if (config.writeResultFile) {
                    ofstream outfile;
                    bool nfile = !file_exists(config.output);
                    outfile.open(config.output, ios::out | ios::app);
                    if (nfile) {
                        outfile << "%instance,n,m,q,seed,querytype,queryfile,algorithm,time_init_ms,time_queries_total_mus,time_queries_avg_mus\n";
                    }
                    outfile << graph_filename << "," << input_graph.number_of_nodes() << "," << input_graph.number_of_edges()
                        << "," << num_queries << "," << config.seed << "," << query_types[q] << "," << config.queryFiles[q] << ","
                        << "full_reach" << "," << time_full_reach_init  << "," << time_fullreach_query << "," << time_fullreach_query/num_queries << "\n";
                    outfile.close();
                }
        }

        int successful_queries = 0;

        if (runSReach) {
            //double elapsed = 0;
            for (auto fb : fallbacks) {
                r.setFallbackMode(fb);
                std::string algorithm;
                std::vector<Query> queries(all_queries[q]);
                switch(fb) {
                    case FB_PIBFS:
                        algorithm = "oreach-pbibfs";
                        std::cout << "\n\nRunning O'Reach with fallback pBiBFS ..." << std::endl;
                        t.restart();
                        //r.query<FB_PIBFS>(all_queries[q]);
                        r.query<FB_PIBFS>(queries);
                        time_scc_query_pbibfs = t.elapsed<std::micro>();
                        time_cur_query = time_scc_query_pbibfs;
                        break;
                    case FB_PDFS:
                        algorithm = "oreach-pdfs";
                        std::cout << "\n\nRunning O'Reach with fallback pDFS ..." << std::endl;
                        t.restart();
                        //r.query<FB_PDFS>(all_queries[q]);
                        r.query<FB_PDFS>(queries);
                        time_scc_query_pdfs = t.elapsed<std::micro>();
                        time_cur_query = time_scc_query_pdfs;
                        break;
                    case FB_NONE:
                    default:
                        algorithm="none";
                        std::cout << "\n\nRunning O'Reach without any fallback! ..." << std::endl;
                        t.restart();
                        r.query<FB_NONE>(queries);
                        time_cur_query = t.elapsed<std::micro>();

                        if (config.writeUnansweredQueries) {
                            std::string base_filename = graph_filename.substr(graph_filename.find_last_of("/\\") + 1);
                            std::string outfile = config.outputUnanswered + base_filename;

                            if (!config.queryFiles.empty()) {
                                base_filename = config.queryFiles[q].substr(config.queryFiles[q].find_last_of("/\\") + 1);
                                outfile = config.outputUnanswered + base_filename;
                            }

                            std::cout << "writing unanswered queries out to " << outfile << "\n";
                            ofstream queryfile;
                            queryfile.open(outfile, ios::out | ios::trunc);
                            for (auto i = 0u; i < queries.size(); i++) {
                                if (queries[i].answeredWithObservation()) {
                                    continue;
                                }
                                queryfile << queries[i].s << " " << queries[i].t << " "
                                << (queries[i].reachable ? "1" : "0") << " 0"
                                << "\n";
                            }
                            queryfile.close();
                        }
                        break;

                }

                if (config.writeResultFile) {
                    ofstream outfile;
                    bool nfile = !file_exists(config.output);
                    outfile.open(config.output, ios::out | ios::app);
                    if (nfile) {
                        outfile << "%instance,n,m,q,seed,querytype,queryfile,algorithm,time_init_ms,time_queries_total_mus,time_queries_avg_mus\n";
                    }
                    outfile << graph_filename << "," << input_graph.number_of_nodes() << "," << input_graph.number_of_edges()
                        << "," << num_queries << "," << config.seed << "," << query_types[q] << "," << config.queryFiles[q] << ","
                        << algorithm << "," << time_scc_init  << "," << time_cur_query << "," << time_cur_query/num_queries << "\n";
                    outfile.close();
                }

                if (successful_queries == 0) {
                    for (auto &res : queries) {
                        if (res.isPositive()) {
                            successful_queries++;
                        }
                    }
                }
#ifdef RDEBUG
                r.debug();
                if (config.writeObservationsFile) {
                    std::vector<unsigned int> observations = r.getObservationDistribution();
                    ofstream outfile;
                    bool nfile = !file_exists(config.observations);
                    outfile.open(config.observations, ios::out | ios::app);

                    if (nfile) {
                        outfile << "%instance,q,seed,querytype,queryfile,algorithm";
                        for (std::uint_fast8_t i = 0; i < observations.size(); i++) {
                            if (!Query::isDummy(i)) {
                                outfile << "," << Query::resultName(Query::castIndex(i));
                            }
                        }
                        outfile << "\n";
                    }

                    outfile << graph_filename << "," << num_queries << "," << config.seed << "," << query_types[q]
                        << "," << config.queryFiles[q] << "," << algorithm;
                    for (std::uint_fast8_t i = 0; i < observations.size(); i++) {
                        if (!Query::isDummy(i)) {
                            outfile << "," << observations[i];
                        }
                    }
                    outfile << "\n";
                    outfile.close();
                }
#endif

                r.reset();
            }
        }


        /*if (config.writeOutputQueryFile) {
            ofstream queryfile;
            queryfile.open(config.outputQueryFile, ios::out | ios::trunc);
            for (int i = 0; i < num_queries; i++) {
                queryfile << queries[i].s << " " << queries[i].t << " "
                << (queries[i].reachable ? "1" : "0") << " "
                << (queries[i].answeredWithObservation ? "1" : "0")
                << "\n";
            }
            queryfile.close();
        } */

        for (auto &res : all_queries[q]) {
            if (res.isPositive()) {
                successful_queries++;
            }
        }

        std::cout << "\n\nSUMMARY:\n";
        std::cout << "\n## number of queries: " << num_queries << "\n";
        std::cout << "## percentage of successful queries: " << ((100.0*successful_queries)/num_queries) << "%\n";
        std::cout << "## number of successful queries: " << (successful_queries) << "\n\n";

        if (config.pruning_bibfs) {
            std::cout << std::setprecision(5) << "O'Reach with Pruning BiBFS took total (ms) " << time_scc_init  + time_scc_query_pbibfs << "\n"
                << "   total for queries (µs): " << time_scc_query_pbibfs << "\n"
                << "   average per query (µs): " << (time_scc_query_pbibfs / num_queries) << "\n\n";
        }
        if (config.pruning_dfs) {
            std::cout << std::setprecision(5) << "O'Reach with Pruning DFS took total (ms) " << time_scc_init  + time_scc_query_pdfs << "\n"
                << "   total for queries (µs): " << time_scc_query_pdfs << "\n"
                << "   average per query (µs): " << (time_scc_query_pdfs / num_queries) << "\n\n";
        }

        std::cout << "===================\n";

    }

    return ret_code;
}

