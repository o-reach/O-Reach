/******************************************************************************
 * definitions.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DEFINITIONS_H_CHR
#define DEFINITIONS_H_CHR

#include <limits>
#include <queue>
#include <vector>

#include <climits>
#include "macros_assertions.h"
#include <cstdio>


// allows us to disable most of the output during partitioning
#ifdef KAFFPAOUTPUT
        #define PRINT(x) x
#else
        #define PRINT(x) do {} while (false);
#endif

/**********************************************
 * Constants
 * ********************************************/
//Types needed for the graph ds
typedef unsigned int 	NodeID;
typedef double 		EdgeRatingType;
typedef unsigned int 	PathID;
typedef unsigned int 	PartitionID;
typedef unsigned int 	NodeWeight;
typedef int 		EdgeWeight;
typedef EdgeWeight 	Gain;
#ifdef MODE64BITEDGES
typedef uint64_t 	EdgeID;
#else
typedef unsigned int 	EdgeID;
#endif
typedef int 		Color;
typedef unsigned int 	Count;
typedef std::vector<NodeID> boundary_starting_nodes;
typedef long FlowType;

const EdgeID UNDEFINED_EDGE            = std::numeric_limits<EdgeID>::max();
const NodeID UNDEFINED_NODE            = std::numeric_limits<NodeID>::max();
const NodeID UNASSIGNED                = std::numeric_limits<NodeID>::max();
const NodeID ASSIGNED                  = std::numeric_limits<NodeID>::max()-1;
const PartitionID INVALID_PARTITION    = std::numeric_limits<PartitionID>::max();
const PartitionID BOUNDARY_STRIPE_NODE = std::numeric_limits<PartitionID>::max();
const int NOTINQUEUE 		       = std::numeric_limits<int>::max();
const int ROOT 			       = 0;



#endif

