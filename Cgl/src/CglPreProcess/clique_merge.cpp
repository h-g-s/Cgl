#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <clique.h>
#include "clique_merge.h"
#include "clique_extender.h"
#include "vint_set.h"
#include "memory.h"
#include <OsiSolverInterface.hpp>
#include <algorithm>

enum CliqueType
{
    NotAClique = 0,
    NotDominated = 1,
    Dominated = 2,
    MoreThanAClique = 3 // do not try to dominate this constraint
};

using namespace std;

class VecNames
{
public:
    VecNames() {
        starts.reserve(8192);
        charSpace.reserve(8192*61);
        starts.push_back( 0 );
    }

    const char *getName( int i ) {
        return &(charSpace[starts[i]]);
    }

    void add( const char *name ) {
        int size = strlen(name);
        charSpace.insert( charSpace.end(), name, name+size );
        charSpace.push_back('\0');
        starts.push_back( *starts.rbegin()+size+1 );
    }

    // should be called because using **
    void updateLines() {
        int nLines = starts.size()-1;
        if (lines) {
            delete[] lines;
            lines = new char*[nLines];
        }
        lines = new char*[nLines];
        for ( int i=0 ; (i<nLines) ; ++i )
            lines[i] = &(charSpace[starts[i]]);
    }

    // returns structure as char ** 
    char ** ptr() { 
        assert(lines);
        return lines; 
    }

    virtual ~VecNames() {
        if (lines)
            delete[] lines;
    }
private:
    vector< char > charSpace;
    vector< int > starts;
    char **lines;
};

// simple sparse matrix
class SpsMtx {
public:
    SpsMtx( int lines, int nzsToReserve = INT_MAX ) :
        starts( lines+1, INT_MAX ),
        nzs( lines, 0 ) {
        if (nzsToReserve != INT_MAX)
            elements.reserve( nzsToReserve );
    }

    void addRow( int line, int nz, const int els[], bool _sort = false ) {
        assert( line >= 0 && line<(int)nzs.size() );
        starts[line] = elements.size();
        nzs[line] = nz;
        elements.insert( elements.end(), els, els+nz );
        if (_sort) 
            std::sort( elements.begin()+starts[line], elements.begin()+starts[line]+nz );
    }

    const int *row( int line ) const {
        assert( line >= 0 && line<(int)nzs.size() );
        if (nzs[line] == 0)
            return NULL;
        return &elements[starts[line]];
    }

    int nz( int line ) const {
        assert( line >= 0 && line<(int)nzs.size() );
        return nzs[line];
    }

private:
    vector< int > starts;
    vector< int > nzs;
    vector< int > elements;
};

#define MAX_SIZE_CLIQUE_TO_BE_EXTENDED 256

// global configurations and stats
int clqMergeVerbose = 1;
double clqMergeSecsCheckClique = 0.0;
double clqMergeSecsExtendAndDominate = 0.0;
double clqMergeSecsAddAndRemove = 0.0;
double clqMergeSecsExtend = 0.0;
int clqMergeNExtended = 0;
int clqMergeNDominatedFull = 0;
int clqMergeNDominatedEqPart = 0;

static void *xmalloc( const size_t size );
static void *xcalloc( const size_t elements, const size_t size );
static void *xrealloc( void *ptr, const size_t size );

/* fills from start until the last element before end */
#define FILL( vector, start, end, value ) { \
    int i; \
    for ( i=start ; (i<end) ; ++i ) vector[i] = value; \
} \


/* allocate filling with zeros */
#define ALLOCATE_VECTOR_INI( ptr, type, nElements ) {\
    ptr = (type*) calloc( (nElements), sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define ALLOCATE_VECTOR( ptr, type, nElements ) {\
    ptr = (type*) malloc( sizeof(type)*(nElements) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }


static char dominates( int nClq1, const int clq1[], int nClq2, const int clq2[], char *iv )
{
    if ( nClq1 < nClq2 )
        return false;

    // filling incidence vector
    for ( int i=0 ; (i<nClq1) ; ++i )
    {
#ifdef DEBUG
        assert( iv[clq1[i]]==false );
#endif
        iv[clq1[i]] = true;
    }

    char res = true;
    // checking if clq2 is contained in clq1
    for ( int i=0 ; (i<nClq2) ; ++i )
    {
        if (iv[clq2[i]]==false)
        {
            res = false;
            break;
        }
    }


    // clearing iv
    for ( int i=0 ; (i<nClq1) ; ++i )
        iv[clq1[i]] = false;
    return res;
}

static void add_clique(
        OsiSolverInterface *mip,
        CliqueSet *newCliques,      // new cliques
        int size, const int el[],  // clique being added
        vector< enum CliqueType > &cliqueState,  // which are the clique rows and their state
        const SpsMtx &origClqs,
        char *iv, // temporary incidence vector
        char *ivrt, // temporary incidence vector
        int *rtc, // list of rows checked
        const int **colClqs, // column cliques
        const int nColClqs[] // number of cliques of a column
        )
{
#ifdef DEBUG
    for ( int i=0 ; (i<mip->getNumCols()) ; ++i )
    {
        assert( iv[i] == false );
    }
    for ( int i=0 ; (i<mip->getNumRows()) ; ++i )
    {
        assert( ivrt[i] == false );
    }
#endif

    int add = 0;
    int nc1 = clq_set_number_of_cliques( newCliques );
    clq_set_add( newCliques, el, size, size );
    int nc2 = clq_set_number_of_cliques( newCliques );
    add = nc2-nc1;

    assert(nc2>=nc1);
    if (add==0)
        return;
        
    int minColClq = INT_MAX;
    int maxColClq = INT_MIN;

    int nrtc = 0;
    for ( int i=0 ; (i<size) ; ++i )
    {
        int col = el[i];
        
        minColClq = min( minColClq, col );
        maxColClq = max( maxColClq, col );
        
        if (col >= mip->getNumCols()) // complimetary variable
            col -= mip->getNumCols();
#ifdef DEBUG            
        assert( col >= 0 && col<mip->getNumCols() );
#endif        
        // checking non-dominated cliques that column appear
        for ( int j=0 ; (j<nColClqs[col]) ; ++j )
        {
            int ir = colClqs[col][j];
#ifdef DEBUG
            assert( cliqueState[ir]!=NotAClique );
#endif
            
            // skipping already dominated, already included or larger rows (cannot be dominated by this one )
            if ( cliqueState[ir]==Dominated || ivrt[ir] || origClqs.nz(ir)>size )
                continue;

            ivrt[ir] = true;
            rtc[nrtc++] = ir;;
        }
    }

    for ( int i=0 ; (i<nrtc) ; ++i )
    {
        int rowClique = rtc[i];
        ivrt[rowClique] = false;
        
        // quick check before going to the slow check
        int nzRow = origClqs.nz(rowClique);
        const int *clqRow = origClqs.row(rowClique);
        if ( size == nzRow )
        {
            if ( ( clqRow[0] != minColClq ) ||
                 ( clqRow[nzRow-1] != maxColClq ) )
                     continue;
        }
        else
        {
            assert( nzRow<size );
            if ( (clqRow[0] < minColClq) ||
                 ( clqRow[nzRow-1] > maxColClq ) )
                     continue;
        }
        
        if (dominates(size, el, nzRow, clqRow, iv ))
        {
            cliqueState[rowClique] = Dominated;
            if (clqMergeVerbose>=2)
                printf("\t\tdominates %s\n", mip->getRowName(rowClique).c_str() );
        } // dominates 
    } // all rows to check
 
}

static int cmp_int( const void *v1, const void *v2 )
{
    const int *i1 = (const int *) v1;
    const int *i2 = (const int *) v2;

    return (*i1)-(*i2);
}

/* detect cliques */
static int check_cliques( 
        OsiSolverInterface *mip, 
        vector< enum CliqueType > &cliqueState,
        vector< int > &cliques, 
        SpsMtx &origClqs )
{
    const CoinPackedMatrix *cpmRow =  mip->getMatrixByRow();
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const double *Arhs = mip->getRightHandSide();
    const char *Asense = mip->getRowSense();
    const double *colLB = mip->getColLower();
    const double *colUB = mip->getColUpper();
    
    // used to fix column indexes
    vector< int > cidx( mip->getNumCols() );

    clock_t startcq = clock();
    if (clqMergeVerbose>=1)
        printf( "checking candidates for clique extension.\n" );
    int rows = mip->getNumRows();
    int nCliques = 0;

    int minClqSize = INT_MAX;
    int maxClqSize = 0;
    double avClgSize = 0.0;

    for ( int i=0 ; (i<rows) ; ++i )
    {
        int nz = cpmRow->getVectorLengths()[i];
        if ( nz<=1 || nz>MAX_SIZE_CLIQUE_TO_BE_EXTENDED )
            continue;

        const int *idx = cpmRow->getIndices() + starts[i];
        const double *coef = cpmRow->getElements() + starts[i];

        const char sense = Asense[i];
        double rhs = Arhs[i];
        // all constraints as <=
        double mult = (sense == 'G') ? -1.0 : 1.0;

        /* all variables should be positive, integers */
        char varsOk = true;
        double minCoef = DBL_MAX;

        // to check if original and complimentary variables are involved
        int nOnes = 0, nMinusOne = 0;

        for ( int j=0 ; (j<nz) ; ++j )
        {
            if ( (!mip->isInteger(idx[j])) || colLB[idx[j]]<=-1e-5 || ( coef[j]<=-1e8 && colUB[idx[j]]>=1.0+1e-8 ) )
            {
                varsOk = false;
                break;
            }

            // binaries, checking for original and complimentary variables
            if ( colLB[idx[j]] >= -1e-8 && colUB[idx[j]] <= 1.0+1e-8 )
            {
                if ( fabs(mult*coef[j]-1.0)<=1e-8 ) 
                    ++nOnes;
                else
                    if ( fabs(mult*coef[j]+1.0)<=1e-8 )
                        ++nMinusOne;
            }

            minCoef = min( mult * coef[j], minCoef );
        }

        if (!varsOk)
            continue;

        // first detect clique only with normal variables
        cliqueState[i] = ( 2*minCoef >= mult*rhs+1e-5 && rhs >= 1e-5 ) ? NotDominated : NotAClique;

        if (cliqueState[i]==NotDominated)
        {
            cliques[nCliques++] = i;

            origClqs.addRow( i, nz, idx, true );

            minClqSize = min( minClqSize, nz );
            maxClqSize = max( maxClqSize, nz );
            avClgSize += nz;
        } // found a clique candidate
        else
        {
            // checking for clique involving normal variables and complementary variables
            if ( (nOnes+nMinusOne==nz) && (fabs(rhs-(1.0-nMinusOne))<=1e-8) )
            {
                memcpy( &(cidx[0]), idx, sizeof(int)*nz );
                cliqueState[i]=NotDominated;
                
                cliques[nCliques++] = i;

                for ( int j=0 ; (j<nz) ; ++j )
                    if (coef[j]<-1e-8)
                        cidx[j] += mip->getNumCols();

                origClqs.addRow( i, nz, &cidx[0], true );

                minClqSize = min( minClqSize, nz );
                maxClqSize = max( maxClqSize, nz );
                avClgSize += nz;
            }
        }

        if ( clqMergeVerbose>=2 && cliqueState[i]==NotDominated )
        {
            printf( "\trow %s: ", mip->getRowName(i).c_str() );
            for ( int j=0 ; (j<nz) ; ++j )
            {
                printf(" %+g %s ", coef[j], mip->getColName(idx[j]).c_str() );
            }
            char strSense[3] = "";
            switch (sense)
            {
                case 'E' :
                {
                    strcpy( strSense, "=") ;
                    break;
                }
                case 'L' :
                {
                    strcpy( strSense, "<=") ;
                    break;
                }
                case 'G' :
                {
                    strcpy( strSense, ">=") ;
                    break;
                }
            } // sense

            printf("%s %g\n", strSense, rhs );
        } // if verbose

    } // all rows

    clqMergeSecsCheckClique = ((double)clock()-startcq) / ((double)CLOCKS_PER_SEC);

    if (clqMergeVerbose>=1)
        printf("model checked in %.4f seconds. %d candidate cliques for extension/merging. clique sizes range:[%d...%d], av %.2f.\n", clqMergeSecsCheckClique, nCliques, minClqSize, maxClqSize, avClgSize/((double)nCliques) );

    return nCliques;
}

/* to sort cliques per size */
struct CliqueSize
{
    int clqIdx;
    int size;
};

static int cmp_clq_size( const void *v1, const void *v2 )
{
    const struct CliqueSize *cs1 = (const struct CliqueSize *) v1;
    const struct CliqueSize *cs2 = (const struct CliqueSize *) v2;

    return cs2->size-cs1->size;
}

void addRow( 
        int *nrRows, int *nrCapRows, int *nrNz, int *nrCapNz, int **nrStart, 
        int **nrIdx, double **nrCoef, char **nrSense, double **nrRhs,
        int nz, const int idx[], const double coef[], char sense, double rhs,
        const char *rowName, 
        VecNames &rowNames )
{
    if ( (*nrRows)+2 >= (*nrCapRows) )
    {
        (*nrCapRows) = max( 2*(*nrCapRows),  (*nrCapRows)+32 );

        (*nrStart) = (int *) xrealloc( (*nrStart), sizeof(int)*((*nrCapRows)+1) );
        (*nrSense) = (char*) xrealloc( (*nrSense), sizeof(char)*((*nrCapRows)) );
        (*nrRhs) = (double *) xrealloc( (*nrRhs), sizeof(double)*(*nrCapRows) );
    }

    if ( ((*nrNz)+nz) >= (*nrCapNz))
    {
        (*nrCapNz) = max( (*nrCapNz)*2, ((*nrNz)+nz) );
        (*nrIdx) = (int *) xrealloc( (*nrIdx), sizeof(int)*(*nrCapNz) );
        (*nrCoef) = (double *) xrealloc( (*nrCoef), sizeof(double)*(*nrCapNz) );
    }

    /* adding */
    int *pidx = (*nrIdx) + (*nrNz);
    double *pcoef = (*nrCoef) + (*nrNz);
    memcpy( pidx, idx, sizeof(int)*nz );
    memcpy( pcoef, coef, sizeof(double)*nz );
    (*nrSense)[(*nrRows)] = sense;
    (*nrRhs)[(*nrRows)] = rhs;

    (*nrStart)[(*nrRows)+1] = (*nrStart)[(*nrRows)] + nz;
    (*nrNz) += nz;

    rowNames.add( rowName );
    
    ++(*nrRows);
}

/* tries to extend every clique in mip using
 * conflict graph cgraph, dominated cliques are removed */
void merge_cliques( void *osi, CGraph *cgraph, int maxExtensions, int maxItBk )
{
    OsiSolverInterface *mip = (OsiSolverInterface *) osi;

    const CoinPackedMatrix *cpmRow =  mip->getMatrixByRow();
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const double *Arhs = mip->getRightHandSide();
    const char *Asense = mip->getRowSense();
    //const double *colLB = mip->getColLower();
    //const double *colUB = mip->getColUpper();
    clock_t startExtend = clock(); // 
  
    double *rc = NULL;      // reduced costs (dummy values by now)

    int *idx = new int[mip->getNumCols()];        // temporary area
    double *coef = new double[mip->getNumCols()]; // to get indexes from rows
    
    // list of rows which are cliques
    vector< int > cliques( mip->getNumRows() ); 

    // if it is a clique or not and if it is dominated
    vector< enum CliqueType > cliqueState( mip->getNumRows(), NotAClique );

    CliqueSet *newCliques = clq_set_create(); // where new cliques will be stored

    int *nColClqs = NULL;    // number of cliques that a column is involved
    int *allCollClqs = NULL; // contiguous vector to store column cliques
    int **colClqs = NULL;    // cliques that a column is involved
    
    /* here cliques are stored only with column indexes
     * complimentary variables are stores with indexes
     * numCols + idxVar, so that coefficients are not necessary
     * to store */
    // clique elements per row
    SpsMtx origClqs( mip->getNumRows(), mip->getNumElements() );

    char *iv = new char[ mip->getNumCols()*2 ];  // incidence vector for columns (including compliment)
    memset( iv, 0, sizeof(char)*(mip->getNumCols()*2) );

    char *ivr = new char[ mip->getNumRows() ];   // incidence vector for rows
    memset( ivr, 0, sizeof(char)*(mip->getNumRows()) );

    int *rtc = new int[mip->getNumRows()];       // rows to check per thread

    int clqSizeCapIni = 4096;
    struct CliqueSize *clqsSize = NULL; // to sort per clique size
    ALLOCATE_VECTOR( clqsSize, struct CliqueSize, clqSizeCapIni );
    int clqSizeCap = clqSizeCapIni; // current capacity in number of cliques for each thread


    VecNames clqNames;

    /* new rows */
    int *nrStart = NULL;
    int *nrIdx = NULL;
    double *nrCoef = NULL;
    double *nrRhs = NULL;
    char *nrSense = NULL;
    int nrRows = 0;
    int nrNz = 0;
    int nrCapRows = 0;
    int nrCapNz = 0;
    
    /* row names */
    VecNames rowNames;

    const int nCliques = check_cliques( mip, cliqueState, cliques, origClqs );
    if ( nCliques == 0 )
        goto TERMINATE;

    /* filling cliques per col */
    {
        int totnz = 0;
        ALLOCATE_VECTOR_INI( nColClqs, int, mip->getNumCols() );
        for ( int i=0 ; (i<nCliques) ; ++i )
        {
            int ir = cliques[i];
            int nzr = origClqs.nz(ir);
            const int *row = origClqs.row(ir);
            assert( nzr >= 2 );
            for ( int j=0 ; (j<nzr) ; ++j )
            {
                int col = row[j];
                if (col >= mip->getNumCols())
                    col -= mip->getNumCols(); 
                ++nColClqs[col];
                ++totnz;
            }
        }
        ALLOCATE_VECTOR( allCollClqs, int, totnz );
        ALLOCATE_VECTOR( colClqs, int *, mip->getNumCols() );
        // start of each column
        colClqs[0] = allCollClqs;
        for ( int i=1 ; i<mip->getNumCols() ; ++i )
            colClqs[i] = colClqs[i-1] + nColClqs[i-1];

        FILL( nColClqs, 0, mip->getNumCols(), 0 );
        // filling cliques of each col
        for ( int i=0 ; (i<nCliques) ; ++i )
        {
            int ir = cliques[i];
            int nzr = origClqs.nz(ir);
            const int *row = origClqs.row(ir);
            assert( nzr >= 2 );
            for ( int j=0 ; (j<nzr) ; ++j )
            {
                int col = row[j];
                if (col >= mip->getNumCols())
                    col -= mip->getNumCols(); 
                colClqs[col][nColClqs[col]++] = ir;
            }
        }
    }


    nrCapRows = nCliques*min(maxExtensions,5) + 1024;
    nrCapNz = mip->getNumElements() * 2;
    ALLOCATE_VECTOR( nrStart, int, nrCapRows+1 );
    ALLOCATE_VECTOR( nrIdx, int, nrCapNz );
    ALLOCATE_VECTOR( nrCoef, double, nrCapNz );
    ALLOCATE_VECTOR( nrRhs, double, nrCapRows );
    ALLOCATE_VECTOR( nrSense, char, nrCapRows );
    nrStart[0] = 0;

    // just to fill parameter
    ALLOCATE_VECTOR( rc, double, mip->getNumCols()*2 );
    FILL( rc, 0, mip->getNumCols()*2, 1.0 );


    startExtend = clock();
    for ( int iclq=0 ; iclq<nCliques ; ++iclq )
    {
        int row = cliques[iclq];

        int nz = cpmRow->getVectorLengths()[row];

        const double *ccoef = cpmRow->getElements() + starts[row];
        const int *cidx = cpmRow->getIndices() + starts[row];

        memcpy( idx, cidx, sizeof(int)*nz );


        if ( cliqueState[row] != Dominated )
        {
            for ( int ii=0 ; (ii<nz) ; ++ii )
                if (ccoef[ii]<=-1e-5)
                    idx[ii] += mip->getNumCols();

            /* extending */
            {
                CliqueExtender *clqe = clqe_create(cgraph);
//                printf("aaa cgraph: %p\n", cgraph ); fflush(stdout); fflush(stderr);
//                printf("rc %p cgraph size %d\n", (void*) rc, cgraph_size(cgraph) ); fflush(stdout); fflush(stderr);
                clqe_set_max_it_bk( clqe, maxItBk );
                clqe_set_costs( clqe,  &rc[0], cgraph_size(cgraph) );

                clock_t startext = clock();
                int status = clqe_extend( clqe, idx, nz, mip->getNumCols(), CLQEM_EXACT );
                clqMergeSecsExtend += ((double)clock()-startext)/(double)CLOCKS_PER_SEC;
                if (status>0)
                {
                    if (clqMergeVerbose>=2)
                        printf("-> constraint %s extended: \n", mip->getRowName(row).c_str() );

                    ++clqMergeNExtended;

                    // to sort cliques found per size
                    struct CliqueSize *clqbs = clqsSize;

                    const CliqueSet *clqs = clqe_get_cliques( clqe );
                    assert(clq_set_number_of_cliques(clqs)>=1);

                    int nCliquesToInsert = min( maxExtensions, clq_set_number_of_cliques(clqs) );
                    if (nCliquesToInsert)
                    {
                        cliqueState[row] = Dominated;
                        if (nCliquesToInsert==1)
                        {
                            clqbs[0].clqIdx = 0;
                            clqbs[0].size = clq_set_clique_size( clqs, 0 );
                        }
                        else
                        {
                            /* resize */
                            if (clq_set_number_of_cliques(clqs)>clqSizeCap)
                            {
                                clqSizeCap = max( clqSizeCap*2, clq_set_number_of_cliques(clqs) );
                                clqbs = clqsSize = (CliqueSize*) xrealloc( clqsSize, sizeof(struct CliqueSize)*clqSizeCap );
                            }
                            /* sorting per size */
                            for ( int ic=0 ; (ic<clq_set_number_of_cliques(clqs)) ; ++ic )
                            {
                                clqbs[ic].clqIdx = ic;
                                clqbs[ic].size = clq_set_clique_size( clqs, ic );
                            }
                            /* maybe needed to sort cliques */

                            qsort( clqbs, clq_set_number_of_cliques(clqs), sizeof(struct CliqueSize), cmp_clq_size );
                            //printf( "from %d to %d\n", clqbs[0].size, clqbs[clq_set_number_of_cliques(clqs)-1].size );
                        }
                    }

                    for ( int ic=0 ; (ic<nCliquesToInsert) ; ++ic )
                    {
                        int idxclique = clqbs[ic].clqIdx;
                        int size = clq_set_clique_size( clqs, idxclique );
                        //printf("   adding %d\n", size );
                        const int *el = clq_set_clique_elements( clqs, idxclique );

                        add_clique( mip, newCliques, size, el, cliqueState, origClqs, iv, ivr, rtc, (const int **) colClqs, nColClqs );
                        char origrname[256] = "";
                        strncpy(origrname, mip->getRowName(row).c_str(), 256 );
                        char extn[64] = "";
                        if (nCliquesToInsert>1)
                        {
                            sprintf(extn, "xt(%d)", ic );
                            strcat( origrname, extn );
                        }
                        if (clqMergeVerbose>=2)
                        {
                            printf("\t%s : ", origrname );
                        }
                        clqNames.add( origrname );
                    
                        if (clqMergeVerbose>=2)
                        {
                            for ( int k=0 ; (k<size) ; ++k )
                            {
                                char strneg[3] = "";
                                int icol = el[k];
                                if (icol>=mip->getNumCols())
                                {
                                    strcpy( strneg, "!");
                                    icol -= mip->getNumCols();
                                }
                                printf("%s%s ", strneg, mip->getColName(icol).c_str() );
                            }
                            printf("\n");
                        }
                    }
                }
                clqe_free( &clqe );
            }
        } // non dominated clique constraints
    } // all clique constraints
    
    clqMergeSecsExtendAndDominate = ((double)clock()-(double)startExtend)/((double)CLOCKS_PER_SEC);

    {
        clock_t startRemAdd = clock();
        /* removing dominated cliques */
        {
            int nToRemove = 0;
            int *toRemove;
            ALLOCATE_VECTOR( toRemove, int, mip->getNumRows() );

            for ( int i=0 ; (i<mip->getNumRows()) ; ++i )
            {
                if ( cliqueState[i] == Dominated )
                {
                    if (mip->getRowSense()[i]=='E')
                    {
                        ++clqMergeNDominatedEqPart;
                        // adding >= part, <= part will be added separately
                        int nz = cpmRow->getVectorLengths()[i];
                        double rhs = Arhs[i];
                        char rname[256];
                        strncpy( rname, mip->getRowName(i).c_str(), 256 );
                        char nrname[512];
                        sprintf( nrname, "%sEp", rname );
                        
                        const double *crcoef = cpmRow->getElements() + starts[i];
                        const int *cridx = cpmRow->getIndices() + starts[i];

                        addRow(  &nrRows, &nrCapRows, &nrNz, &nrCapNz, &nrStart, &nrIdx, &nrCoef, &nrSense, &nrRhs,
                            nz, cridx, crcoef, 'G', rhs, nrname, rowNames );
                    }
                    else
                    {
                        ++clqMergeNDominatedFull;
                    }
                    toRemove[nToRemove++] = i;
                }
            }


            if (nToRemove)
                mip->deleteRows( nToRemove, toRemove );

            free( toRemove );
        }

        /* adding stronger cliques */
        if (clq_set_number_of_cliques(newCliques)>0)
        {
            int nNewCliques = clq_set_number_of_cliques(newCliques);

            for ( int ic=0 ; (ic<nNewCliques) ; ++ic )
            {
                int size = clq_set_clique_size( newCliques, ic );
                const int *el = clq_set_clique_elements( newCliques, ic );
                double rhs = 1.0;
                for ( int i=0 ; (i<size) ; ++i )
                {
                    if (el[i]>=mip->getNumCols())
                    {
                        coef[i] = -1.0;
                        rhs -= 1.0;
                        idx[i] = el[i]-mip->getNumCols();
                    }
                    else
                    {
                        idx[i] = el[i];
                        coef[i] = 1.0;
                    }
                }
                // adding
                addRow(  &nrRows, &nrCapRows, &nrNz, &nrCapNz, &nrStart, &nrIdx, &nrCoef, &nrSense, &nrRhs,
                    size, idx, coef, 'L', rhs,
                    clqNames.getName(ic), rowNames );
            }
        }
        
        clqMergeSecsAddAndRemove = ((double)clock()-(double)startRemAdd)/((double)CLOCKS_PER_SEC);
    }
    
    /* flushing rows */
    {
        vector< double > rowLB( nrRows ), rowUB( nrRows );
        for ( int i=0 ; (i<nrRows) ; ++i ) {
            switch (nrSense[i]) {
                case 'E':
                    rowLB[i] = rowUB[i] = nrRhs[i];
                    break;
                case 'L':
                    rowUB[i] = nrRhs[i];
                    rowLB[i] = -DBL_MAX;
                    break;
                case 'G':
                    rowLB[i] = nrRhs[i];
                    rowUB[i] = DBL_MAX;
                    break;
            }
        }
        mip->addRows( nrRows, nrStart, nrIdx, nrCoef, &rowLB[0], &rowUB[0] );
    }

    
    
    if (clqMergeVerbose)
    {
        printf("%d extended, %d dom full, %d dom eq.\n", clqMergeNExtended, clqMergeNDominatedFull, clqMergeNDominatedEqPart );
        fflush( stdout );
    }

TERMINATE:
    if (nrStart)
        free( nrStart );
    if (nrIdx)
        free( nrIdx );
    if (nrCoef)
        free( nrCoef );
    if (nrRhs)
        free( nrRhs );
    if (nrSense)
        free( nrSense );

    if (rc)
        free( rc );

    delete[] idx;
    delete[] coef;

    free( clqsSize );

    clq_set_free( &newCliques );

    delete[] iv;
    delete[] ivr;
    delete[] rtc;

    if (nColClqs)
        free(nColClqs);
    if (allCollClqs )
        free(allCollClqs);
    if (colClqs)
        free( colClqs);
}

static void *xmalloc( const size_t size )
{
   void *result = malloc( size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
      exit(1);
   }

   return result;
}

static void *xcalloc( const size_t elements, const size_t size )
{
   void *result = calloc( elements, size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to callocte %zu bytes.", size);
      exit(1);
   }

   return result;
}

static void *xrealloc( void *ptr, const size_t size )
{
   void *result = realloc( ptr, size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
      exit(1);
   }
   return result;
}
