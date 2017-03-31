//CHHS defines.h
//CHHS Store #define statements and other "global" preprocessor statements
//CHHS does not incorporate #include statements, conditional codes encapsulated by #ifdef, and safeguards
#ifndef _DEFINES_H
#define _DEFINES_H 1		//Safeguard to avoid multiple includes of the same things
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define FOR(i,n) for(i=0; i < (n); ++i)
#define FPN(file) fputc('\n', file)
#define F0 stdout
#define Pi  3.1415926535897932384626433832795
#define QuantileGamma(prob,alpha,beta) QuantileChi2(prob,2.0*(alpha))/(2.0*(beta))
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))
#define DGammaMean 0
#define FAST_RANDOM_NUMBER
#define m2NormalSym  0.95
#define MAXNFIELDS 10000
#define PAML_RELEASE      0
#define VerStr "slimcodeml version 2014_Feb_11"
#define NS            5000
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define MAXNSONS      100
#define NGENE         2000
#define LSPNAME       50
#define NCODE         64
#define NCATG         41
#define NBTYPE        17
#define NP            (NBRANCH*2+NGENE-1+2+NCODE+2)
#define CODEML 1
#ifdef  BASEML
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#define MINIMIZATION
#endif

#ifdef  CODEML
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#define MINIMIZATION
#endif

#ifdef  BASEMLG
#define REALSEQUENCE
#define NODESTRUCTURE
#define LSDISTANCE
#endif

#ifdef  RECONSTRUCTION
#define PARSIMONY
#endif

#ifdef  MCMCTREE
#define REALSEQUENCE
#define NODESTRUCTURE
#define LFUNCTIONS
#endif

#define gammap(x,alpha) (alpha*(1-pow(x,-1.0/alpha)))
#define NBESTANC  4		/* use 1 2 3 or 4 */
#define MAXNF2D  5
#endif //CHHS Safeguard
