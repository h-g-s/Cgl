#ifndef SPATHS_H
#define SPATHS_H

typedef struct _ShortestPathsFinder ShortestPathsFinder;
typedef ShortestPathsFinder *ShortestPathsFinderPtr;

/*
 * Creates Shortest Path Finder.
 * arcStart[i] indicates the position in vector
 * toNode and dist where arcs of node i start
*/
ShortestPathsFinderPtr spf_create(const int nodes, const int arcs, const int *arcStart, const int *toNode, const int *dist);

/* updates just one arc
 **/
void spf_update_arc(ShortestPathsFinder *spf, const int tail, const int head, const int cost);

/* returns the distance of one arc
 **/
int spf_get_arc(ShortestPathsFinder *spf, const int tail, const int head);

/*
 * queries number of nodes
 */
int spf_nodes(ShortestPathsFinder *spf);

/*
 * queries number of arcs
 */
int spf_arcs(ShortestPathsFinder *spf);

/*
 * executes the shortest path finder
 * using the Dijkstra algorithm
 */
void spf_find(ShortestPathsFinder *spf, const int origin);
void spf_find(ShortestPathsFinder *spf, const int origin, const int destination);

/*
 * solution query: returns distance to a node after executing spf_find
 */
int spf_get_dist(const ShortestPathsFinderPtr spf, const int node);

/*
 * solution query: returns previous node and allows one to build a path after executing spf_find
 */
int spf_get_previous(const ShortestPathsFinderPtr spf, const int node);

int *spf_previous(const ShortestPathsFinder *spf);

/* returns all previous nodes
 * which should be steped
 * to arrive at a given node (this node is not included)
 * returns how many nodes were filles in indexes
 */
int spf_get_path(const ShortestPathsFinder *spf, const int toNode, int indexes[]);

/*
 * releases Shortest Path Finder object
 */
void spf_free(ShortestPathsFinderPtr *spf);

#endif /* ifndef SPATHS_H */
