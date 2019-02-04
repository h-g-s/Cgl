#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cassert>
#include <cstring>
#include <ctime>
#include <cstdio>
#include <limits>
#include <algorithm>

#include "cgraph.h"
#include "oddhs.h"
#include "spaths.h"

#define ODDHS_DEF_MIN_FRAC 1e-5 // minimum to a variable to be considered active
#define DIST_MAX 100000 // maximum  distance for nodes with conflicts
#define MAX_WHEEL_CENTERS 256

typedef struct {
  int node;
  double cost;
} NodeCost;

bool cmp_node_cost(const NodeCost &np1, const NodeCost &np2)
{
  if (fabs(np1.cost - np2.cost) > ODDHS_DEF_MIN_FRAC) {
    return np1.cost < np2.cost;
  }

  return np1.node < np2.node;
}

struct _OddHoleSep {
  const CGraph *cgraph;

  const double *x;
  const double *rc;

  // from integer interesting columns
  // those which are active in solution
  int icaCount;
  int *icaIdx; // original index
  int *icaActivity; // mapping of the fractional solution value to an
  // integer value to made further computations easier

  // shortest path data
  std::vector< int > spArcStart; // start index for arcs of each node
  std::vector< int > spArcTo; // destination of each arc
  std::vector< int > spArcDist; // distance for each arc

  int *spOrigPathIdx; // to store discovered shortest path

  // incidence vector to check for repeated entries
  char *ivreIdx;

  ShortestPathsFinder *spf;

  // discovered odd holes
  // indexes are stored related to icaIdx
  // (only considering active variables)
  std::vector< std::vector< int > > doh; // indexes of all odd holes
  std::vector< double > dohViol; //violation of all odd holes
  std::vector< double > dohRHS; //RHS of all odd holes
  std::vector< std::vector< int > > dohWC; // indexes of all wheel centers

  char *dohIV; // incidence vector for variables considered here
};

OddHoleSep *oddhs_create(const CGraph *cgraph)
{
  OddHoleSep *oddhs = new OddHoleSep;
  const int numCols = cgraph_size(cgraph);

  oddhs->cgraph = cgraph;

  oddhs->icaCount = 0;
  oddhs->icaIdx = new int[numCols];
  oddhs->icaActivity = new int[numCols];
  oddhs->spOrigPathIdx = new int[numCols];
  oddhs->dohIV = new char[numCols];
  oddhs->ivreIdx = new char[numCols];

  return oddhs;
}

// methods to work with current fractional solution
int oddhs_fill_active_intcols(OddHoleSep *oddhs);

void oddhs_prepare_dist_graph(OddHoleSep *oddhs);

/**
 * Odd Hole store management functions.
 * Odd-holes are processed using ica indexes.
 * Only after all processing that they are
 * translated to original variable indexes.
 **/

/** Tries to add a newly discovered add hole that generates a violated cut.
 *  If it is a repeated entry, ignores it and returns zero, if suceed in
 *  inserting, returns one.
 **/
bool oddhs_add_doh(OddHoleSep *oddhs, const int nz, const int _idx[]);

/** tests if a newly discovered odd hole is not a repeated entry **/
int oddhs_doh_already_exists(OddHoleSep *oddhs, const std::vector< int > &doh);

/**
 * Tries to find an Odd Hole using the conflict graph
 * this odd hole may correspond to a violated cut.
 * the search may be aggressive or not
 * returns how many new odd holes were found
 **/
int oddhs_find_odd_holes_with_node(OddHoleSep *oddhs, const int node);

int oddhs_vector_has_repeated_entries(OddHoleSep *oddhs, const int nz, const int idx[]);

/**
 * tries to find wheel centers so that Odd Hole inequelities may be lifted
 */
size_t oddhs_search_wheel_centers_all_dohs(OddHoleSep *oddhs);

/**
 * finds at most "maxCenters" wheel centers, for one odd hole
 * priority is:
 *   1. most fractional one
 *   2. those with smallest reduced cost
 **/
size_t oddhs_find_wheel_centers(OddHoleSep *oddhs, const int dohIdx);

// transform a double [0 , 1] into an integer
// ready to use in the distance graph
int oddhs_get_activity(const double val);

int oddhs_search_odd_holes(OddHoleSep *oddhs, const double x[], const double rc[])
{
  oddhs->x = x;
  oddhs->rc = rc;

  int icCount = oddhs_fill_active_intcols(oddhs);

  if (icCount <= 4)
    return 0;

  oddhs->doh.clear();
  oddhs->doh.reserve(1024);
  oddhs->dohViol.clear();
  oddhs->dohViol.reserve(1024);
  oddhs->dohRHS.clear();
  oddhs->dohRHS.reserve(1024);
  oddhs->dohWC.clear();
  oddhs->dohWC.reserve(1024);
  oddhs->spArcStart.clear();
  oddhs->spArcStart.reserve(1024);
  oddhs->spArcTo.clear();
  oddhs->spArcTo.reserve(2048);
  oddhs->spArcDist.clear();
  oddhs->spArcDist.reserve(2048);

  oddhs_prepare_dist_graph(oddhs);

  int result = 0;
  for (int i = 0; i < oddhs->icaCount; i++) {
    result += oddhs_find_odd_holes_with_node(oddhs, i);
  }

  oddhs_search_wheel_centers_all_dohs(oddhs);

#ifdef DEBUG
  assert(oddhs->doh.size() == oddhs->dohRHS.size());
  assert(oddhs->doh.size() == oddhs->dohViol.size());
#endif

  return result;
}

size_t oddhs_search_wheel_centers_all_dohs(OddHoleSep *oddhs)
{
  /* trying to extend odd holes by including wheel centers */
  size_t insertedWC = 0;

  for (size_t i = 0; i < oddhs->doh.size(); i++) {
    const size_t nwc = oddhs_find_wheel_centers(oddhs, i);
    insertedWC += nwc;
  }

  return insertedWC;
}

int oddhs_fill_active_intcols(OddHoleSepPtr oddhs)
{
  int cnCount = cgraph_size(oddhs->cgraph);

  if (cnCount <= 4) {
    return 0;
  }

  oddhs->icaCount = 0;
  for (int j = 0; j < cnCount; j++) {
    if (cgraph_degree(oddhs->cgraph, j) < 2) {
      continue;
    }
    if (oddhs->x[j] < ODDHS_DEF_MIN_FRAC) {
      continue;
    }

#ifdef DEBUG
    assert(oddhs->x[j] >= -0.000001 && oddhs->x[j] <= 1.000001);
#endif

    oddhs->icaIdx[oddhs->icaCount] = j;
    oddhs->icaActivity[oddhs->icaCount] = oddhs_get_activity(oddhs->x[j]);
    oddhs->icaCount++;
  }

  return oddhs->icaCount;
}

void oddhs_prepare_dist_graph(OddHoleSepPtr oddhs)
{
  int idxArc = 0;
  const int nodes = oddhs->icaCount * 2;

  //Conflicts: (x', y'')
  for (int i1 = 0; i1 < oddhs->icaCount; i1++) {
    oddhs->spArcStart.push_back(idxArc);
    const int idx1 = oddhs->icaIdx[i1];

    for (int i2 = 0; i2 < oddhs->icaCount; i2++) {
      const int idx2 = oddhs->icaIdx[i2];

      if (cgraph_conflicting_nodes(oddhs->cgraph, idx1, idx2)) {
        oddhs->spArcTo.push_back(oddhs->icaCount + i2);
        oddhs->spArcDist.push_back(oddhs->icaActivity[i2]);
        idxArc++;
      } // conflict found
    } // i2
  } // i1

  //Conflicts: (x'', y')
  for (int i = 0; i < oddhs->icaCount; i++) {
    oddhs->spArcStart.push_back(idxArc);

    for (int j = oddhs->spArcStart[i]; j < oddhs->spArcStart[i + 1]; j++) {
      const int arcTo = oddhs->spArcTo[j] - oddhs->icaCount;
#ifdef DEBUG
      assert(arcTo >= 0 && arcTo < oddhs->icaCount);
#endif
      const int arcDist = oddhs->spArcDist[j];
      oddhs->spArcTo.push_back(arcTo);
      oddhs->spArcDist.push_back(arcDist);
      idxArc++;
    }
  }

  oddhs->spArcStart.push_back(idxArc);

  oddhs->spf = spf_create(nodes, idxArc, &oddhs->spArcStart[0], &oddhs->spArcTo[0], &oddhs->spArcDist[0]);
}

int oddhs_get_activity(const double val)
{
#ifdef DEBUG
  assert(val >= -0.000001 && val <= 1.000001);
#endif
  int intval = (int)nearbyint(DIST_MAX * val);

  if (intval > DIST_MAX) {
    intval = DIST_MAX;
  } else if (intval < 0) {
    intval = 0;
  }

  return DIST_MAX - intval;
}

int oddhs_find_odd_holes_with_node(OddHoleSep *oddhs, const int node)
{
  const int icaCount = oddhs->icaCount;
  const int dest = icaCount + node;
  ShortestPathsFinder *spf = oddhs->spf;
  int *origPathIdx = oddhs->spOrigPathIdx;
  int result = 0;

  spf_find(spf, node, dest);
  int oddSize = spf_get_path(spf, dest, origPathIdx) - 1;

  if (oddSize < 5) {
    return 0;
  }

  // translating indexes
  for (int j = 0; j < oddSize; ++j) {
    origPathIdx[j] %= oddhs->icaCount;
  }

  // checking for repeated entries
  if (oddhs_vector_has_repeated_entries(oddhs, oddSize, origPathIdx)) {
    return 0;
  }

  /* checking if it is violated */
  double lhs = 0;
  for (int i = 0; i < oddSize; i++) {
#ifdef DEBUG
    assert(origPathIdx[i] < oddhs->icaCount);
#endif
    lhs += oddhs->x[oddhs->icaIdx[origPathIdx[i]]];
  }
  const double rhs = floor(oddSize / 2.0);
  const double viol = lhs - rhs;
  if (viol < ODDH_SEP_DEF_MIN_VIOL) {
    return 0;
  }

  if (oddhs_add_doh(oddhs, oddSize, origPathIdx)) {
    oddhs->dohRHS.push_back(rhs);
    oddhs->dohViol.push_back(viol);
    result++;
  }

  return result;
}

int oddhs_vector_has_repeated_entries(OddHoleSep *oddhs, const int nz, const int idx[])
{
  if (nz <= 1) {
    return 0;
  }

#ifdef DEBUG
  assert(oddhs->icaCount <= cgraph_size(oddhs->cgraph));
#endif

  std::fill(oddhs->ivreIdx, oddhs->ivreIdx + oddhs->icaCount, 0);

  for (int i = 0; i < nz; i++) {
#ifdef DEBUG
    assert(idx[i] >= 0 && idx[i] < oddhs->icaCount);
#endif
    if (oddhs->ivreIdx[idx[i]]) {
      return 1;
    }

    oddhs->ivreIdx[idx[i]] = 1;
  }

  return 0;
}

bool oddhs_add_doh(OddHoleSep *oddhs, const int nz, int const _idx[])
{
  std::vector< int > doh(nz);

  for (int i = 0; i < nz; i++) {
    doh[i] = oddhs->icaIdx[_idx[i]];
  }

  // checking for repeated entries
  if (oddhs_doh_already_exists(oddhs, doh)) {
    return false;
  }

  // inserting
  oddhs->doh.push_back(doh);

  return true;
}

int oddhs_doh_already_exists(OddHoleSep *oddhs, const std::vector< int > &doh)
{
  const int cols = cgraph_size(oddhs->cgraph);
  const size_t dohSize = doh.size();

  std::fill(oddhs->dohIV, oddhs->dohIV + cols, 0);

  for (int idx : doh) {
    oddhs->dohIV[idx] = 1;
  }

  for (size_t idxDOH = 0; idxDOH < oddhs->doh.size(); idxDOH++) {
    // checking size
    const size_t otherSize = oddhs->doh[idxDOH].size();
    if (dohSize != otherSize) {
      continue;
    }

    bool isEqual = true;
    for (size_t j = 0; j < dohSize; j++) {
      if (!(oddhs->dohIV[oddhs->doh[idxDOH][j]])) {
        isEqual = false;
        break;
      }
    }
    if (isEqual) {
      return 1;
    }
  }

  return 0;
}

const std::vector< int > &oddhs_get_odd_hole(const OddHoleSep *oddhs, const int idx)
{
#ifdef DEBUG
  assert(idx >= 0 && idx <= oddhs->doh.size());
#endif
  return oddhs->doh[idx];
}

const size_t oddhs_get_odd_hole_size(const OddHoleSep *oddhs, const int idx)
{
#ifdef DEBUG
  assert(idx >= 0 && idx <= oddhs->doh.size());
#endif
  return oddhs->doh[idx].size();
}

const double oddhs_get_odd_hole_violation(const OddHoleSep *oddhs, const int idx)
{
#ifdef DEBUG
  assert(idx >= 0 && idx <= oddhs->doh.size());
#endif
  return oddhs->dohViol[idx];
}

const double oddhs_get_odd_hole_rhs(const OddHoleSep *oddhs, const int idx)
{
#ifdef DEBUG
  assert(idx >= 0 && idx <= oddhs->doh.size());
#endif
  return oddhs->dohRHS[idx];
}

const size_t oddhs_get_odd_hole_count(const OddHoleSep *oddhs)
{
  return oddhs->doh.size();
}

size_t oddhs_find_wheel_centers(OddHoleSep *oddhs, const int dohIdx)
{
  const int numCols = cgraph_size(oddhs->cgraph);
  char *iv = new char[numCols]();
  const std::vector< int > &oh = oddhs_get_odd_hole(oddhs, dohIdx);

  oddhs->dohWC.emplace_back();
  std::vector< int > &dohWC = oddhs->dohWC[oddhs->dohWC.size() - 1];
  dohWC.reserve(std::min(MAX_WHEEL_CENTERS, numCols));

  /* node with the smallest degree */
  int nodeSD = -1, minDegree = std::numeric_limits< int >::max();
  /* picking node with the smallest degree */
  for (const int &node : oh) {
    const int dg = cgraph_degree(oddhs->cgraph, node);
    if (dg < minDegree) {
      minDegree = dg;
      nodeSD = node;
    }

    iv[node] = 1; //incidence vector
  }

  int *neighs = new int[numCols];
  int nConfs = cgraph_get_all_conflicting(oddhs->cgraph, nodeSD, neighs, numCols);

  NodeCost *np = new NodeCost[nConfs];
  int npSize = 0;

  for (int i = 0; i < nConfs; i++) {
    const int idx = neighs[i];

    if (cgraph_degree(oddhs->cgraph, idx) < oh.size()) {
      continue;
    }

    /* checking if it has conflict with all nodes and if it is not included in oh */
    if (iv[idx]) {
      continue;
    }

    bool insert = true;
    for (const int &j : oh) {
      if (!cgraph_conflicting_nodes(oddhs->cgraph, idx, j)) {
        insert = false;
        break;
      }
    }

    if (!insert) {
      continue;
    }

    np[npSize].node = idx;

    if (oddhs->x[idx] >= ODDHS_DEF_MIN_FRAC) {
      np[npSize].cost = (oddhs->x[idx] * 1000.0);
    } else { /* avoiding negative numbers */
      np[npSize].cost = (1000000.0 + oddhs->rc[idx]);
    }
#ifdef DEBUG
    assert(np[npSize].cost >= 0);
#endif
    npSize++;
  }

  if (npSize) {
    if (npSize <= MAX_WHEEL_CENTERS) {
      std::sort(np, np + npSize, cmp_node_cost);
    } else {
      std::partial_sort(np, np + MAX_WHEEL_CENTERS, np + npSize, cmp_node_cost);
      npSize = MAX_WHEEL_CENTERS;
    }

    for (int i = 0; i < npSize; i++) {
      const int node = np[i].node;
      /* must have conflict with all other centers */
      bool insert = true;
      for (const int &c : dohWC) {
        if (!cgraph_conflicting_nodes(oddhs->cgraph, node, c)) {
          insert = false;
          break;
        }
      }

      if (!insert) {
        continue;
      }

      dohWC.push_back(node);
    }
  }

  delete[] iv;
  delete[] np;
  delete[] neighs;

  return dohWC.size();
}

size_t oddhs_get_nwc_doh(OddHoleSep *oddhs, const int doh)
{
#ifdef DEBUG
  assert(doh < oddhs->doh.size());
#endif
  return oddhs->dohWC[doh].size();
}

const std::vector< int > &oddhs_get_wc_doh(OddHoleSep *oddhs, const int doh)
{
#ifdef DEBUG
  assert(doh < oddhs->doh.size());
#endif
  return oddhs->dohWC[doh];
}

void oddhs_free(OddHoleSepPtr *oddhs)
{
  delete[](*oddhs)->icaIdx;
  delete[](*oddhs)->icaActivity;
  delete[](*oddhs)->dohIV;
  delete[](*oddhs)->spOrigPathIdx;
  delete[](*oddhs)->ivreIdx;

  spf_free(&((*oddhs)->spf));

  delete (*oddhs);
  *oddhs = nullptr;
}
