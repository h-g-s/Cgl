#include <cstdio>
#include <cassert>

#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "CoinTime.hpp"
#include "CglUsrOddHole.hpp"

extern "C" {
	#include "cut.h"
	#include "macros.h"
	#include "oddhs.h"
}

int CglUsrOddHole::sepCuts = 0;
double CglUsrOddHole::sepTime = 0.0;

CglUsrOddHole::CglUsrOddHole()
{
}

CglUsrOddHole::CglUsrOddHole(const CglUsrOddHole& rhs)
{
}

CglCutGenerator * CglUsrOddHole::clone() const
{
    CglUsrOddHole *aClq = new CglUsrOddHole();

    return static_cast<CglCutGenerator*>(aClq);
}

void CglUsrOddHole::fillColSolution(const OsiSolverInterface &si, double colSol[])
{
    const int numCols = si.getNumCols();
    const double* siColSol = si.getColSolution();

    for(int i = 0; i < numCols; i++)
    {
        colSol[i] = siColSol[i];
        colSol[i+numCols] = 1.0 - siColSol[i];
    }
}

void CglUsrOddHole::fillReducedCost(const OsiSolverInterface &si, double rCost[])
{
    const int numCols = si.getNumCols();
    const double* siRCost = si.getReducedCost();

    for(int i = 0; i < numCols; i++)
    {
        rCost[i] = siRCost[i];
        rCost[i+numCols] = -siRCost[i];
    }
}

void CglUsrOddHole::generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info ) {
    double startSep = CoinCpuTime();
    const int numCols = si.getNumCols();
    const CGraph *cg = si.getCGraph();
    double *x = new double[numCols*2];
    double *rc = new double[numCols*2];
    int *idx = new int[numCols*2];
    int *idxMap = new int[numCols*2];
    double *coef = new double[numCols*2];
    CutPool *cutPool = cut_pool_create(numCols);
    OddHoleSep *oddhs = oddhs_create(cg);
    const CoinAbsFltEq equal(1.0e-12);
    OsiRowCut osrc;
    
	if(numCols != cgraph_size(cg)/2) {
        fprintf(stderr, "Invalid conflict graph! Number of columns %d ... in graph %d\n",
                numCols, cgraph_size(cg)/2);
        exit(EXIT_FAILURE);
    }

    fillColSolution(si, x);
    fillReducedCost(si, rc);
    oddhs_search_odd_holes(oddhs, x, rc);

    /* adding odd holes */
    for(int j = 0; j < oddhs_get_odd_hole_count(oddhs); j++) {
        const std::vector<int> &oddEl = oddhs_get_odd_hole(oddhs, j);
        const size_t oddSize = oddEl.size();
        double viol = oddhs_get_odd_hole_violation(oddhs, j);
        double rhs = oddhs_get_odd_hole_rhs(oddhs, j);

        if(viol < ODDH_SEP_DEF_MIN_VIOL) {
            continue;
        }

        if(oddSize < 5) {
            fprintf(stderr, "Invalid size of cut: %lu\n", oddSize);
            exit(EXIT_FAILURE);
        }

        const size_t centerSize = oddhs_get_nwc_doh(oddhs, j);
        const std::vector<int> &centerIdx = oddhs_get_wc_doh(oddhs, j);

        const size_t cutSize = oddSize + centerSize;
        int realSize = 0;
        int dup = 0;

        std::fill(coef, coef + cutSize, 0.0);
        std::fill(idxMap, idxMap + numCols, -1);

        for(int k = 0; k < oddSize; k++) {
            if(oddEl[k] < numCols) {
                if(idxMap[oddEl[k]] == -1) {
                    idxMap[oddEl[k]] = realSize;
                    idx[realSize] = oddEl[k];
                    coef[realSize] = 1.0;
                    realSize++;
                }
                else {
                    coef[idxMap[oddEl[k]]] += 1.0;
                    dup++;
                }
            }
            else {
                if(idxMap[oddEl[k]-numCols] == -1) {
                    idxMap[oddEl[k]-numCols] = realSize;
                    idx[realSize] = oddEl[k] - numCols;
                    coef[realSize] = -1.0;
                    realSize++;
                }
                else {
                    coef[idxMap[oddEl[k]-numCols]] -= 1.0;
                    dup++;
                }
                rhs = rhs - 1.0;
            }
        }

        if(centerSize && fabs(rhs) > 1e-6) {
            const double oldRhs = rhs;
            for(int k = 0; k < centerSize; k++) {
                if(centerIdx[k] < numCols) {
                    if(idxMap[centerIdx[k]] == -1) {
                        idxMap[centerIdx[k]] = realSize;
                        idx[realSize] = centerIdx[k];
                        coef[realSize] = oldRhs;
                        realSize++;
                    }
                    else {
                        coef[idxMap[centerIdx[k]]] += oldRhs;
                        dup++;
                    }
                }
                else {
                    if(idxMap[centerIdx[k]-numCols] == -1) {
                        idxMap[centerIdx[k]-numCols] = realSize;
                        idx[realSize] = centerIdx[k] - numCols;
                        coef[realSize] = -1.0 * oldRhs;
                        realSize++;
                    }
                    else {
                        coef[idxMap[centerIdx[k]-numCols]] -= oldRhs;
                        dup++;
                    }
                    rhs = rhs - oldRhs;
                }
            }
        }

        if(dup) {
            int last = 0;
            for(int k = 0; k < realSize; k++) {
                if(fabs(coef[k]) > 1e-6) {
                    idx[last] = idx[k];
                    coef[last] = coef[k];
                    last++;
                }
            }
            realSize = last;
        }

        Cut *cut = cut_create(idx, coef, realSize, rhs, x);
        cut_pool_insert(cutPool, cut);
        cut_free(&cut);
    }

    cut_pool_update(cutPool);
    const int numberRowCutsBefore = cs.sizeRowCuts();
    for(int i = 0; i < cut_pool_size(cutPool); i++) {
        const Cut *cut = cut_pool_get_cut(cutPool, i);
        osrc.setRow(cut_size(cut) , cut_get_idxs(cut), cut_get_coefs(cut));
        osrc.setUb(cut_get_rhs(cut));
        cs.insertIfNotDuplicate(osrc, equal);
    }

    const int numberRowCutsAfter = cs.sizeRowCuts();
    CglUsrOddHole::sepCuts += numberRowCutsAfter - numberRowCutsBefore;

    if(!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
        int numberRowCutsAfter = cs.sizeRowCuts();
        for(int i = numberRowCutsBefore; i < numberRowCutsAfter; i++)
            cs.rowCutPtr(i)->setGloballyValid();
    }

    cut_pool_free(&cutPool);
    oddhs_free( &oddhs );
    delete[] x;
    delete[] rc;
    delete[] idx;
    delete[] idxMap;
    delete[] coef;
	
	CglUsrOddHole::sepTime += (CoinCpuTime() - startSep);
}

CglUsrOddHole::~CglUsrOddHole()
{

}
