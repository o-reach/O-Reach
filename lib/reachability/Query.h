/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Kathrin Hanauer
 *      Jonathan Trummer
 *****************************************************************************/

#ifndef REACH_QUERY_H
#define REACH_QUERY_H

#include "extern/KaHIP/lib/definitions.h"

struct Query {
    enum struct Result : std::uint_fast8_t {
        UNKNOWN = 0,
        KNOWN,  // dummy, use as barrier

        // negative
        SRC_SNK,           // t source or s sink
        DIFF_WCC,          // different weak components
        TS_FW_LEVELS,      // contradicting forward levels
        TS_BW_LEVELS,      // contradicting backward levels
        TS_NUM,            // contradicting topsort numbers
        TS_MAX,            // t beyond max range of s
        TS_MIN,            // s beyond min range of t
        SUPPORT_NEG_FW,    // SV2 / aka 02
        SUPPORT_NEG_BW,    // SV3 / aka 03
        FALLBACK_NEG,      // used fallback, answer negative
        PREACH_LEVELS,
        PREACH_MIN_FWD,
        PREACH_MIN_BWD,
        PREACH_MAX_BWD,
        PREACH_MAX_FWD,
        PREACH_FALLBACK_NEG,

        POSITIVE,          // dummy, use as barrier

        // positive
        SAME_VERTEX,       // s == t
        SAME_SCC,          // same strong component
        TS_HI,             // t in high range of s
        TS_LO,             // s in low range of t
        SUPPORT_POS,       // SV1 / aka 01
        TS_EQMAX,          // max range of s ends at t
        TS_EQMIN,          // min range of t ends at s
        EDGE,              // edge (s, t)
        FALLBACK_POS,      // used fallback, answer positive
        PREACH_SAME_VERTEX,
        PREACH_RANGE_FWD,
        PREACH_RANGE_BWD,
        PREACH_FALLBACK_POS,

        LAST               // dummy, use for #possibilities
    };

    NodeID s;
    NodeID t;
    bool reachable;
    Result result = Result::UNKNOWN;

    Query(NodeID s, NodeID t): s(s), t(t) { }

    bool answered() const { return result > Result::KNOWN; }

    bool isPositive() const { return result > Result::POSITIVE; }

    bool isNegative() const { return answered() && !isPositive(); }

    bool answeredWithObservation() const {
        return (isPositive() && result != Result::FALLBACK_POS)
            || (isNegative() && result != Result::FALLBACK_NEG);
    }

    long unsigned int resultToIndex() const {
        return castResult(result);
    }

    void setResultFromBinary() {
        if (reachable == 1) {
            result = Result::FALLBACK_POS;
        } else {
            result = Result::FALLBACK_NEG;
        }
    }
    void setBinaryFromResult() {
        reachable = isPositive();
    }

    constexpr
    static long unsigned int castResult(Result result) {
        return static_cast<long unsigned int>(result);
    }

    constexpr
        static Result castIndex(std::uint_fast8_t &index) {
            return static_cast<Result>(index);
        }

    constexpr
    static const char* resultName(Result result) {
            switch (result) {
                case Result::UNKNOWN:
                    return "unknown";
                case Result::KNOWN:
                    return "known";
                case Result::SRC_SNK:
                    return "neg_source-sink";
                case Result::DIFF_WCC:
                    return "neg_diff-wcc";
                case Result::TS_FW_LEVELS:
                    return "neg_fw-level";
                case Result::TS_BW_LEVELS:
                    return "neg_bw-level";
                case Result::TS_NUM:
                    return "neg_ts-num";
                case Result::TS_MAX:
                    return "neg_ts-max";
                case Result::TS_MIN:
                    return "neg_ts-min";
                case Result::SUPPORT_NEG_FW:
                    return "neg_support-fw";
                case Result::SUPPORT_NEG_BW:
                    return "neg_support-bw";
                case Result::FALLBACK_NEG:
                    return "neg_fallback";
                case Result::PREACH_LEVELS:
                    return "neg_pr-levels";
                case Result::PREACH_MIN_FWD:
                    return "neg_pr-min-fw";
                case Result::PREACH_MIN_BWD:
                    return "neg_pr-min-bw";
                case Result::PREACH_MAX_BWD:
                    return "neg_pr_max-bw";
                case Result::PREACH_MAX_FWD:
                    return "neg_pr_max-fw";
                case Result::PREACH_FALLBACK_NEG:
                    return "neg_pr_fallback";
                case Result::POSITIVE:
                    return "positive";
                case Result::SAME_VERTEX:
                    return "pos_same-vertex";
                case Result::SAME_SCC:
                    return "pos_same-scc";
                case Result::TS_HI:
                    return "pos_ts-high";
                case Result::TS_LO:
                    return "pos_ts-low";
                case Result::SUPPORT_POS:
                    return "pos_support";
                case Result::TS_EQMAX:
                    return "pos_ts-max";
                case Result::TS_EQMIN:
                    return "pos_ts-min";
                case Result::EDGE:
                    return "pos_edge";
                case Result::FALLBACK_POS:
                    return "pos_fallback";
                case Result::PREACH_SAME_VERTEX:
                    return "pos_preach-same-vertex";
                case Result::PREACH_RANGE_FWD:
                    return "pos_preach-range-fw";
                case Result::PREACH_RANGE_BWD:
                    return "pos_preach-range-bw";
                case Result::PREACH_FALLBACK_POS:
                    return "pos_pr-fallback";
                default:
                    return "dummy";
            }
        }

    constexpr
    static bool isDummy(Result result) {
        return result == Result::KNOWN
            || result == Result::POSITIVE
            || result == Result::LAST;
    }

    constexpr
    static bool isDummy( std::uint_fast8_t index) {
        auto result = castIndex(index);
        return isDummy(result);
    }
};

#endif

