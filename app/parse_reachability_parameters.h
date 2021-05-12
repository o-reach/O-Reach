/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Jonathan Trummer
 *      Kathrin Hanauer
 *****************************************************************************/

#pragma once
#ifndef KAHIP_PARSE_REACHABILITY_PARAMETERS_H
#define KAHIP_PARSE_REACHABILITY_PARAMETERS_H

#include <string>
#include <argtable3.h>
#include "app/reachabilityconfig.h"

int parse_reachability_parameters(int argn, char **argv, std::string &graph_filename, ReachabilityConfig &reachability_config) {

    const char *progname = argv[0];

    struct arg_lit *help = arg_lit0(NULL, "help", "Print help.");
    struct arg_str *filename = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
    struct arg_int *seed = arg_int0(NULL, "seed", NULL, "random seed to use");
    struct arg_lit *use_pbibfs = arg_lit0(NULL, "pbibfs", "activate pruning bibfs as fallback");
    struct arg_lit *use_pdfs = arg_lit0(NULL, "pdfs", "activate pruning DFS as fallback");
    struct arg_lit *run_edgecheck = arg_lit0(NULL, "check-edge", "run check-edge before running fallback");
    struct arg_lit *run_fullreach = arg_lit0(NULL, "full-reach", "run full-reach implementation");

    struct arg_lit *output_result_file = arg_lit0(NULL, "write-results", "write results to results.csv");
    struct arg_lit *output_observations_file = arg_lit0(NULL, "write-observations", "write observations to results.csv");
    struct arg_str *input_query_file = arg_strn(NULL, "queries", NULL, 0, 10, "path to queryfile, if non given, queries will be generated");
    struct arg_str *output_query_file = arg_strn(NULL, "queries-output", NULL, 0, 1, "path to output queryfile");
    struct arg_str *output_condensation = arg_strn(NULL, "condensation-output", NULL, 0, 1, "path to output condensation graph");
    struct arg_str *output_unanswered = arg_strn(NULL, "unanswered-output", NULL, 0, 1, "path to where unanswered queries should be written");
    struct arg_int *support_variant = arg_int0("s", "support", "1|2|3|4|5", "which support calculation type to use, default 3");
    struct arg_dbl *support_l = arg_dbl0("l", NULL, "0-1", "sampling size per support, activates support=2");
    struct arg_dbl *support_l2 = arg_dbl0("k", NULL, "0-1", "sampling size per support, activates support=3");
    struct arg_int *threshold_level_size = arg_int0("t", NULL, "natural number", "threshold for level size when support=3, default=10");
    struct arg_str *result_filename = arg_strn(NULL, "output", NULL, 0, 1, "override path output csv file. enables output of results");
    struct arg_str *observations_filename = arg_strn(NULL, "observations", NULL, 0, 1, "write observations to this file");

    // struct arg_int *seed = arg_int0(NULL, "seed", NULL, "Seed to use for PRNG.");
    struct arg_end *end = arg_end(100);

    void *argtable[] = {
            help,
            filename,
            seed,
            use_pbibfs,
            use_pdfs,
            run_edgecheck,
            run_fullreach,
            result_filename,
            observations_filename,
            input_query_file,
            output_query_file,
            output_result_file,
            output_observations_file,
            output_condensation,
            output_unanswered,
            support_variant,
            support_l,
            support_l2,
            threshold_level_size,
            end
    };

    // Parse arguments.
    int nerrors = arg_parse(argn, argv, argtable);

    // Catch case that help was requested.
    if (help->count > 0) {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        arg_print_glossary(stdout, argtable,"  %-40s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }


    if (nerrors > 0) {
        arg_print_errors(stderr, end, progname);
        printf("Try '%s --help' for more information.\n",progname);
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }


    if (filename->count > 0) {
        graph_filename = filename->sval[0];
        reachability_config.graph_input = filename->sval[0];
    }

    if (seed->count > 0) {
        reachability_config.seed = seed->ival[0];
    }

    if (use_pbibfs->count > 0) {
        reachability_config.pruning_bibfs = true;
    }

    if (use_pdfs->count > 0) {
        reachability_config.pruning_dfs = true;
    }

    if (run_edgecheck->count > 0) {
        reachability_config.edgeCheck = true;
    }

    if (run_fullreach->count > 0) {
        reachability_config.fullReach = true;
    }

    if (result_filename->count > 0) {
        reachability_config.output = result_filename->sval[0];
        reachability_config.writeResultFile = true;
    }

    if (observations_filename->count > 0) {
        reachability_config.observations = observations_filename->sval[0];
        reachability_config.writeObservationsFile = true;
    }

    if (input_query_file->count > 0) {
        reachability_config.useQueryFile = true;
        // reachability_config.inputQueryFile = input_query_file->sval[0];
        for (int i = 0; i < input_query_file->count; i++) {
            reachability_config.queryFiles.push_back(input_query_file->sval[i]);
        }
    }

    if (output_query_file->count > 0) {
        reachability_config.writeOutputQueryFile = true;
        reachability_config.outputQueryFile = output_query_file->sval[0];
    }

    if (output_condensation->count > 0) {
        reachability_config.writeCondensationGraph = true;
        reachability_config.outputCondensation = output_condensation->sval[0];
    }

    if (output_unanswered->count > 0) {
        reachability_config.writeUnansweredQueries = true;
        reachability_config.outputUnanswered = output_unanswered->sval[0];
    }

    if (output_result_file->count > 0) {
        reachability_config.writeResultFile = true;
    }

    if (output_observations_file->count > 0) {
        reachability_config.writeObservationsFile = true;
    }

    if (support_variant->count > 0) {
        reachability_config.support = support_variant->ival[0];
    }

    if (support_l->count > 0) {
        reachability_config.support = 2;
        reachability_config.l = support_l->dval[0];
    }

    if (support_l2->count > 0) {
        reachability_config.support = 3;
        reachability_config.l2 = support_l2->dval[0];
    }

    if (threshold_level_size->count > 0) {
        if (threshold_level_size->ival[0] >= 0) {
            reachability_config.t = threshold_level_size->ival[0];
        }
    }

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}

#endif //KAHIP_PARSE_REACHABILITY_PARAMETERS_H
