#ifndef CLIQUE_MERGE_H
#define CLIQUE_MERGE_H

#include "cgraph.h"

/** @brief merge cliques in a MIP
 *
 *  Merge cliques to produce a stronger
 *  mixed integer programming formulation
 *
 *  @param osi an OsiSolverInterface object containing a Mixed-Integer Program (MIP)
 *  @param maxExtensions maximum number of larger cliques generated from a single clique
 *  @param maxItBk maximum number of iteration in the Bron-kerbosch algorithm to extend a clique
 **/

void merge_cliques( void *osi, CGraph *cgraph, int maxExtensions, int maxItBk );

#endif

