#include <cstdio>
#include <cstring>
#include <cassert>
#include <limits>
#include <algorithm>
#include "spaths.h"
#include "node_heap.h"

#define SP_INFTY_DIST ((std::numeric_limits<int>::max()/2)-1)
#define NULL_NODE (-1)

typedef struct {
    int node;
    int distance;
} Neighbor;

struct _ShortestPathsFinder {
    int nodes;
    int arcs;

    Neighbor *neighs; // all neighbors
    Neighbor **startn; // Start of neighbors for node i. The neighbor ends at startn[i+1]

    // solution:
    int *dist;
    int *previous;
    int *path; // temporary storage for path

    NodeHeap *nh;
};

int compNeighs(const void *n1, const void *n2) {
    const Neighbor *pn1 = (const Neighbor *) n1;
    const Neighbor *pn2 = (const Neighbor *) n2;

    return pn1->node - pn2->node;
}

bool compare_neighbors(const Neighbor &n1, const Neighbor &n2) {
    return n1.node < n2.node;
}

ShortestPathsFinderPtr spf_create(const int nodes, const int arcs, const int *arcStart, const int *toNode, const int *dist) {
    ShortestPathsFinder *spf = new ShortestPathsFinder;

    spf->nodes = nodes;
    spf->arcs = arcs;
    spf->startn = new Neighbor*[nodes + 1];
    spf->nh = nh_create(nodes, SP_INFTY_DIST);
    spf->previous = new int[nodes];
    spf->dist = new int[nodes];
    spf->path = new int[nodes];
    spf->neighs = new Neighbor[arcs];

    for (int idx = 0; idx < arcs; idx++) {
        spf->neighs[idx].node = toNode[idx];
        spf->neighs[idx].distance = dist[idx];
#ifdef DEBUG
        assert(spf->neighs[idx].node < spf->nodes);
#endif
    }

    for (int n = 0; n <= spf->nodes; n++) {
        spf->startn[n] = spf->neighs + arcStart[n];
    }

    return spf;
}

void spf_find(ShortestPathsFinder *spf, const int origin) {
    NodeHeap *nh = spf->nh;

    nh_reset(nh);

    for (int i = 0; i < spf->nodes; i++) {
        spf->dist[i] = SP_INFTY_DIST;
        spf->previous[i] = NULL_NODE;
    }

    spf->dist[origin] = 0;
    nh_update(nh, origin, 0);

    int topCost, topNode;
    while ((topCost = nh_remove_first(nh, &topNode)) < SP_INFTY_DIST) {
        // updating neighbors distances by iterating in all neighbors
        for (Neighbor *n = spf->startn[topNode]; n < spf->startn[topNode+1]; n++) {
            const int toNode = n->node;
            const int dist = n->distance;
            const int newDist = topCost + dist;

            if (spf->dist[toNode] > newDist) {
                spf->previous[toNode] = topNode;
                spf->dist[toNode] = newDist;
                nh_update(spf->nh, toNode, newDist);
            } // updating heap if necessary
        } // going through node neighbors
    } // going through all nodes in priority queue
}

void spf_find(ShortestPathsFinder *spf, const int origin, const int destination) {
    NodeHeap *nh = spf->nh;

    nh_reset(nh);

    for (int i = 0; i < spf->nodes; i++) {
        spf->dist[i] = SP_INFTY_DIST;
        spf->previous[i] = NULL_NODE;
    }

    spf->dist[origin] = 0;
    nh_update(nh, origin, 0);

    int topCost, topNode;
    while ((topCost = nh_remove_first(nh, &topNode)) < SP_INFTY_DIST) {
        if(topNode == destination) {
            break;
        }
        // updating neighbors distances by iterating in all neighbors
        for (Neighbor *n = spf->startn[topNode]; n < spf->startn[topNode+1]; n++) {
            const int toNode = n->node;
            const int dist = n->distance;
            const int newDist = topCost + dist;

            if (spf->dist[toNode] > newDist) {
                spf->previous[toNode] = topNode;
                spf->dist[toNode] = newDist;
                nh_update(spf->nh, toNode, newDist);
            } // updating heap if necessary
        } // going through node neighbors
    } // going through all nodes in priority queue
}

int spf_nodes(ShortestPathsFinder *spf) {
    return spf->nodes;
}

int spf_arcs(ShortestPathsFinder *spf) {
    return spf->arcs;
}

int spf_get_dist(const ShortestPathsFinderPtr spf, const int node) {
#ifdef DEBUG
    assert(node < spf->nodes);
#endif
    return spf->dist[node];
}

int spf_get_previous(const ShortestPathsFinderPtr spf, const int node) {
#ifdef DEBUG
    assert(node < spf->nodes);
#endif
    return spf->previous[node];
}

int *spf_previous(const ShortestPathsFinder *spf) {
    return spf->previous;
}

int spf_get_path(const ShortestPathsFinder *spf, const int toNode, int indexes[]) {
    int n = 0;
    int currNode = toNode;

    spf->path[n++] = currNode;

    while((currNode = spf->previous[currNode]) != NULL_NODE) {
        spf->path[n++] = currNode;
    }

    for (int i = 0; i < n; i++) {
        indexes[i] = spf->path[n-i-1];
    }

    return n;
}

void spf_free(ShortestPathsFinderPtr *spf) {
    delete[] (*spf)->neighs;
    delete[] (*spf)->startn;
    nh_free(&((*spf)->nh));
    delete[] (*spf)->previous;
    delete[] (*spf)->dist;
    delete[] (*spf)->path;

    delete (*spf);
    (*spf) = nullptr;
}

void spf_update_arc(ShortestPathsFinder *spf, const int head, const int tail, const int cost) {
    const Neighbor *start = spf->startn[head];
    const Neighbor *end = spf->startn[head + 1];
    const Neighbor key = {tail, 0};
    Neighbor *result = (Neighbor *) bsearch(&key, start, end - start, sizeof(Neighbor), compNeighs);

#ifdef DEBUG
    if (!result) {
        fprintf(stderr, "Could not find arc (%d->%d) in graph.\n", head, tail);
        abort();
    }
    assert(result->node == tail);
#endif
    result->distance = cost;
}

int spf_get_arc(ShortestPathsFinder *spf, const int head, const int tail) {
    const Neighbor *start = spf->startn[head];
    const Neighbor *end = spf->startn[head + 1];
    const Neighbor key = {tail, 0};
    Neighbor *result = (Neighbor *) bsearch(&key, start, end - start, sizeof(Neighbor), compNeighs);

#ifdef DEBUG
    if (!result) {
        fprintf(stderr, "Could not find arc (%d->%d) in graph.\n", head, tail);
        abort();
    }
    assert(result->node == tail);
#endif
    return result->distance;
}
