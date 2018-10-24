#include <cstdio>
#include <cassert>

#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "CoinTime.hpp"
#include "CglAggressiveClique.hpp"
#include "bron_kerbosch.h"
#include "clique_separation.h"
#include "cut.h"

int CglAggressiveClique::sepCuts = 0;
double CglAggressiveClique::sepTime = 0.0;

CglAggressiveClique::CglAggressiveClique() : maxItBK(1000), maxItBKExt(100), extMethod(3)
{
}

CglAggressiveClique::CglAggressiveClique(const CglAggressiveClique& rhs)
{
    this->maxItBK = rhs.maxItBK;
    this->maxItBKExt = rhs.maxItBKExt;
    this->extMethod = rhs.extMethod;
}

CglCutGenerator * CglAggressiveClique::clone() const
{
    CglAggressiveClique *aClq = new CglAggressiveClique();

    aClq->maxItBK = this->maxItBK;
    aClq->maxItBKExt = this->maxItBKExt;
    aClq->extMethod = this->extMethod;

    return static_cast<CglCutGenerator*>(aClq);
}

void fillColSolution(const OsiSolverInterface &si, double colSol[])
{
    const int numCols = si.getNumCols();
    const double* siColSol = si.getColSolution();

    for(int i = 0; i < numCols; i++)
    {
        colSol[i] = siColSol[i];
        colSol[i+numCols] = 1.0 - siColSol[i];
    }
}

void fillReducedCost(const OsiSolverInterface &si, double rCost[])
{
    const int numCols = si.getNumCols();
    const double* siRCost = si.getReducedCost();

    for(int i = 0; i < numCols; i++)
    {
        rCost[i] = siRCost[i];
        rCost[i+numCols] = -siRCost[i];
    }
}

void CglAggressiveClique::generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info )
{
    double startSep = CoinCpuTime();

    const CoinAbsFltEq equal(1.0e-12);
    OsiRowCut osrc;

    const int numCols = si.getNumCols();
    const CGraph *cg = si.getCGraph();

    double *x = new double[numCols*2];
    double *rc = new double[numCols*2];
    int *idxs = new int[numCols*2];
    double *coefs = new double[numCols*2];
    CliqueSeparation *sep = clq_sep_create(cg);

    fillColSolution(si, x);
    fillReducedCost(si, rc);    
    clq_sep_set_rc(sep, rc);

    clq_sep_set_extend_method(sep, this->extMethod);
    clq_sep_set_max_it_bk(sep, this->maxItBK);
    clq_sep_set_max_it_bk_ext(sep, this->maxItBKExt);

    clq_sep_separate(sep, x);

    /* inserting cliques */
    const CliqueSet *clqSet = clq_sep_get_cliques(sep);
    CutPool *cutPool = cut_pool_create(numCols);
    for(int i = 0; i < clq_set_number_of_cliques(clqSet); i++)
    {
        int dup = 0;
        const int clqSize = clq_set_clique_size(clqSet, i);
        const int *el = clq_set_clique_elements(clqSet, i);
        double lhs = 0.0;
        double rhs = 1.0;
        char *iv = new char[numCols*2]();

        for(int j = 0; j < clqSize; j++)
        {
            if(iv[el[j]])
            {
                fprintf(stderr, "Variable already included!\n");
                exit(1);
            }

            iv[el[j]] = 1;
            if(el[j] < numCols)
            {
                coefs[j] = 1.0;
                idxs[j] = el[j];
                lhs += x[idxs[j]];

                if(iv[el[j]+numCols] == 1)
                    dup++;
            }
            else
            {
                coefs[j] = -1.0;
                idxs[j] = el[j]-numCols;
                rhs -= 1.0;
                lhs -= x[idxs[j]];

                if(iv[el[j]-numCols] == 1)
                    dup++;
            }
        }
        if(dup > 1)
        {
            fprintf(stderr, "Error at clique cut generation module\n");
            exit(1);
        }

        if(dup)
        {
            for(int j = 0; j < clqSize; j++)
            {
                int var = idxs[j];
                int complement = var + numCols;
                if(iv[var] && iv[complement])
                    continue;

                int iTmp[1]; iTmp[0] = var;
                double dTmp[1]; dTmp[0] = 1.0;

                if(coefs[j] == -1.0)
                {
                    dTmp[0] = -1.0;
                    Cut *cut = cut_create(iTmp, dTmp, 1, -1.0, si.getColSolution());
		            cut_pool_insert(cutPool, cut);
		            cut_free(&cut);
                }
                else
                {
                    Cut *cut = cut_create(iTmp, dTmp, 1, 0.0, si.getColSolution());
		            cut_pool_insert(cutPool, cut);
		            cut_free(&cut);
                }
            }

        }
        else
        {
            Cut *cut = cut_create(idxs, coefs, clqSize, rhs, si.getColSolution());
            cut_pool_insert(cutPool, cut);
		    cut_free(&cut);
        }

        delete[] iv;
    }

    cut_pool_update(cutPool);

    const int numberRowCutsBefore = cs.sizeRowCuts();
    for(int i = 0; i < cut_pool_size(cutPool); i++)
    {
        const Cut *cut = cut_pool_get_cut(cutPool, i);
        osrc.setRow(cut_size(cut) , cut_get_idxs(cut), cut_get_coefs(cut));
        osrc.setUb(cut_get_rhs(cut));
        cs.insertIfNotDuplicate(osrc, equal);
    }

    const int numberRowCutsAfter = cs.sizeRowCuts();
    CglAggressiveClique::sepCuts += numberRowCutsAfter - numberRowCutsBefore;

    if(!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
        int numberRowCutsAfter = cs.sizeRowCuts();
        for(int i = numberRowCutsBefore; i < numberRowCutsAfter; i++)
            cs.rowCutPtr(i)->setGloballyValid();
    }

    clq_sep_free(&sep);
    cut_pool_free(&cutPool);
    delete[] x;
    delete[] rc;
    delete[] idxs;
    delete[] coefs;

    CglAggressiveClique::sepTime += (CoinCpuTime() - startSep);
}

CglAggressiveClique::~CglAggressiveClique()
{

}

void CglAggressiveClique::setMaxItBK(int _maxItBK) {
    if(_maxItBK <= 0) {
        fprintf(stderr, "Invalid value for parameter maxItBK (%d).\n", _maxItBK);
        exit(1);   
    }

    this->maxItBK = _maxItBK;
}

void CglAggressiveClique::setMaxItBKExt(int _maxItBKExt) {
    if(_maxItBKExt <= 0) {
        fprintf(stderr, "Invalid value for parameter maxItBKExt (%d).\n", _maxItBKExt);
        exit(1);   
    }
    
    this->maxItBKExt = _maxItBKExt;
}

void CglAggressiveClique::setExtendingMethod(int _extMethod) {
    if(_extMethod < 0 || _extMethod > 4) {
        fprintf(stderr, "Invalid value for parameter extMethod (%d).\n", _extMethod);
        exit(1);   
    }
    
    this->extMethod = _extMethod;
}