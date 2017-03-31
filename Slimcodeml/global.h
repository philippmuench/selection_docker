// CHHS global.h, global variables declarations (cf. global.cpp)
#ifndef _GLOBAL_H
#define _GLOBAL_H 1		//Safeguard to avoid multiple includes of the same things
typedef enum
{ BASEseq = 0, CODONseq, AAseq, CODON2AAseq } DataTypes_typedef;	//CHHS see below
//CHHS extern enum {BASEseq=0, CODONseq, AAseq, CODON2AAseq} DataTypes; //CHHS causes problem with g++, see statements above and below
extern DataTypes_typedef DataTypes;	//CHHS see above
typedef enum
{ PrBranch = 1, PrNodeNum = 2, PrLabel = 4, PrAge = 8, PrOmega = 16 } OutTreeOptions_typedef;	//CHHS see blow
//CHHS extern enum {PrBranch=1, PrNodeNum=2, PrLabel=4, PrAge=8, PrOmega=16} OutTreeOptions; //CHHS causes problem with g++, see statements above and below
extern OutTreeOptions_typedef OutTreeOptions;	//CHHS see above
extern char BASEs[];
extern char AAs[];
extern char CODONs[256][4];
extern const char *EquateBASE[];	//CMV added const
extern char nChara[256], CharaMap[256][64];
extern char AA3Str[];
extern char BINs[];
extern int GeneticCode[][64];
extern int noisy, Iround, NFunCall, NEigenQ, NPMatUVRoot;
extern double SIZEp;		//CHHS This variable is also used as local in tools.c!
extern int AlwaysCenter;
extern double Small_Diff;
extern double prob_Quantile, *par_Quantile;	//CHHS dropped "static" property (should have no effect)
extern double (*cdf_Quantile) (double x, double par[]);	//CHHS This is a function pointer, dropped "static" property
extern unsigned int z_rndu;	//CHHS dropped "static" property (should have no effect)
extern int w_rndu;		//CHHS dropped "static" property (should have no effect)
extern time_t time_start;	//CHHS dropped "static" property (should have no effect)
extern struct TREEB
{
  int nbranch, nnode, root, branches[NBRANCH][2];
  double lnL;
} tree;
extern struct common_info
{
  char *z[NS], *spname[NS], seqf[96], outf[96], treef[96], daafile[96], cleandata;	//CHHS *z[NS] lost "unsigned" property to unify pointers
  char oldconP[NNODE];		/* update conP for nodes? to save computation */
  int seqtype, ns, ls, ngene, posG[NGENE + 1], lgene[NGENE], npatt, *pose,
    readpattern;
  int runmode, clock, verbose, print, codonf, aaDist, model, NSsites;
  int nOmega, nbtype, nOmegaType;	/* branch partition, AA pair (w) partition */
  int method, icode, ncode, Mgene, ndata, bootstrap;
  int fix_rgene, fix_kappa, fix_omega, fix_alpha, fix_rho, nparK, fix_blength,
    getSE;
  int np, ntime, nrgene, nkappa, npi, nrate, nalpha, ncatG, hkyREV;
  size_t sconP, sspace;
  double *fpatt, *space, kappa, omega, alpha, rho, rgene[NGENE];
  double pi[NCODE], piG[NGENE][64], fb61[64];
  double f3x4[NGENE][12], *pf3x4, piAA[20];
  double freqK[NCATG], rK[NCATG], MK[NCATG * NCATG], daa[20 * 20], *conP,
    *fhK;
  double (*plfun) (double x[], int np);
  double omega_fix;		/* fix the last w in the NSbranchB, NSbranch2 models
				   for lineages.  Useful for testing whether w>1 for some lineages. */
  int conPSiteClass;		/* conPSiteClass=0 if (method==0) and =1 if (method==1)?? */
  int NnodeScale;
  char *nodeScale;		/* nScale[ns-1] for interior nodes */
  double *nodeScaleF;		/* nScaleF[npatt] for scale factors */
  /* pomega & pkappa are used to communicate between SetParameters & ConditionalPNode
     & EigenQcodon.  Try to remove them? */
  double *pomega, pkappa[5], *ppi;
} com;
extern struct TREEN
{
  int father, nson, sons[MAXNSONS], ibranch, ipop;
  double branch, age, omega, *conP, label;
  char *nodeStr, fossil /*, usefossil */ ;	//CMV unused
} *nodes, **gnodes, nodes_t[2 * NS - 1];
/* for sptree.nodes[].fossil: lower, upper, bounds, gamma, inverse-gamma */
typedef enum
{ LOWER_F = 1, UPPER_F, BOUND_F } FOSSIL_FLAGS_typedef;	//CHHS see below
//CHHS extern enum {LOWER_F=1, UPPER_F, BOUND_F} FOSSIL_FLAGS; //CHHS g++ doesnt like it, typedef used; see above and below
extern FOSSIL_FLAGS_typedef FOSSIL_FLAGS;	//CHHS see above
extern char *fossils[];
extern struct SPECIESTREE
{
  int nbranch, nnode, root, nspecies, nfossil;
  struct TREESPN
  {
    char name[LSPNAME + 1], fossil, usefossil;	/* fossil: 0, 1, 2, 3 */
    int father, nson, sons[2];
    double age, pfossil[7];	/* lower and upper bounds or alpha & beta */
    //double *lnrates;          /* log rates for loci */ //CMV unused
  } nodes[2 * NS - 1];
} sptree;
/* all trees are binary & rooted, with ancestors unknown. */
extern struct DATA		/* locus-specific data and tree information */
{
  int ns[NGENE], ls[NGENE], npatt[NGENE], ngene /*, lgene[NGENE] */ ;	//CMV unused
  int root[NGENE + 1], /*BlengthMethod, */ fix_nu, nbrate[NGENE], icode[NGENE];	//CMV unused
  char *z[NGENE][NS], cleandata[NGENE];
  char /*idaafile[NGENE], */ daafile[NGENE][40];	//CMV unused
  double *fpatt[NGENE] /*, lnpT, lnpR, lnpDi[NGENE] */ ;	//CMV unused
  double /*Qfactor[NGENE], */ pi[NGENE][NCODE];	//CMV unused
  double /*rgene[NGENE], */ kappa[NGENE], alpha[NGENE], omega[NGENE];	//CMV unused
  int NnodeScale[NGENE];
  char *nodeScale[NGENE];	/* nScale[data.ns[locus]-1] for interior nodes */
} data;
extern int Nsensecodon, FROM61[64], FROM64[64], FourFold[4][4];
extern int ChangedInIteration;	/* 1: t changed, update P(t); 2: paras changed, update UVRoot */
extern double *PMat, *U, *V, *Root, *_UU[NBTYPE + 2], *_VV[NBTYPE + 2],
  *_Root[NBTYPE + 2];
/* 5 sets for branchsite models (YN2002); 6 sets for clade models */
extern double pcodon0[64], paa0[20], *pcodonClass;	/* for aaDist=FIT1 or FIT2 */
extern int BayesEB;		/* =1 for site models M2a & M8; =2 for branch-site models A & C */
extern int LASTROUND;
extern int IClass;
extern int OmegaAA[190], AA1STEP[190];
extern double _rateSite;
extern double Qfactor_NS, Qfactor_NS_branch[NBTYPE];
extern double AAchem[][20 + 1];
extern FILE *fout, *frub, *flnf, *frst, *frst1, *frst2, *finitials;
extern const char *ratef;	//CMV added const
typedef enum
{ Fequal, F1x4, F3x4, Fcodon, F1x4MG, F3x4MG, FMutSel0, FMutSel } CodonFreqs_typedef;	//CHHS see below
//CHHS extern enum {Fequal, F1x4, F3x4, Fcodon, F1x4MG, F3x4MG, FMutSel0, FMutSel} CodonFreqs; //CHHS g++ does not like it, see above and below
extern CodonFreqs_typedef CodonFreqs;	//CHHS see above
extern const char *codonfreqs[];	//CMV added const
typedef enum
{ NSbranchB = 1, NSbranch2, NSbranch3 } NSBranchModels_typedef;	//CHHS see below
//CHHS extern enum {NSbranchB=1, NSbranch2, NSbranch3} NSBranchModels; //CHHS g++ does not like it, see above and below
extern NSBranchModels_typedef NSBranchModels;	//CHHS see above
extern const char *NSbranchmodels[];	//CMV added const
typedef enum
{ Poisson, EqualInput, Empirical, Empirical_F, FromCodon = 6, REVaa_0 = 8, REVaa = 9 } AAModel_typedef;	//CHHS see below
//CHHS extern enum {Poisson, EqualInput, Empirical, Empirical_F,FromCodon=6, REVaa_0=8, REVaa=9} AAModel; //CHHS g++ does not like it, see above and below
extern AAModel_typedef AAModel;	//CHHS see above
extern const char *aamodels[];	//CMV added const
typedef enum
{ NSnneutral = 1, NSpselection, NSdiscrete, NSfreqs, NSgamma, NS2gamma, NSbeta, NSbetaw, NSbetagamma, NSbeta1gamma, NSbeta1normal, NS02normal, NS3normal } NSsitesModels_typedef;	//CHHS see blow
//CHHS extern enum {NSnneutral=1, NSpselection, NSdiscrete, NSfreqs, NSgamma, NS2gamma,NSbeta, NSbetaw, NSbetagamma, NSbeta1gamma, NSbeta1normal, NS02normal,NS3normal} NSsitesModels; //CHHS g++ does not like it, see above and below
extern NSsitesModels_typedef NSsitesModels;	//CHHS see above
extern const char *NSsitesmodels[];	//CMV added const
typedef enum
{ FIT1 = 11, FIT2 = 12 } SiteClassModels_typedef;	//CHHS see below
//CHHS extern enum {FIT1=11, FIT2=12} SiteClassModels; //CHHS g++ does not like it, see above and below
extern SiteClassModels_typedef SiteClassModels;	//CHHS see above
typedef enum
{ AAClasses = 7 } aaDistModels_typedef;	//CHHS see below
//CHHS extern enum {AAClasses=7 } aaDistModels; //CHHS g++ does not like it, see above and below
extern aaDistModels_typedef aaDistModels;	//CHHS see above
extern const char *clockstr[];
typedef enum
{ GlobalClock = 1, LocalClock, ClockCombined } ClockModels_typedef;	//CHHS see below
//CHHS extern enum {GlobalClock=1, LocalClock, ClockCombined} ClockModels; //CHHS g++ does not like it, see above and below
extern ClockModels_typedef ClockModels;	//CHHS see above
/* variables for batch run of site models */
extern int ncatG0, insmodel, nnsmodels, nsmodels[14];
/* used for sliding windows analysis */
extern int windowsize0, offset0, npositive;
extern double lnLmodel;
extern int _nestS;		/* 189= estimate the S elements, 0= use those from com.daa[] */
extern double *_Fij;		//CHHS "static" property removed
extern int ijAAref;		//CHHS "static" property removed
//CHHS The following variables moved from treesub.c
extern double *dfsites;		//CHHS lost "static" property
#ifdef LSDISTANCE		//CHHS We keep the #ifdef, although not very effective for defining single variables
extern double *SeqDistance;
extern int *ancestor;
#endif
#if(defined(BASEML) || defined(CODEML))	//CHHS We keep the #ifdef, although not very effective for defining single variables
extern double *AgeLow;
extern int NFossils, AbsoluteRate;
extern double ScaleTimes_TipDate, TipDate;
/* number of internal node times, usd to deal with known ancestors.  Broken? */
extern int innode_time;		//CHHS "static" property dropped (should have no effect)
#endif
extern long counter_PMatUVRoot;	//CHHS count number of invocations
extern int ONE;			//CHHS Mainly for BLAS and LAPACK function calls
extern double DONE;		//CHHS For BLAS
extern int ZERO;		//CHHS For BLAS
extern double DZERO;		//CHHS For BLAS
extern FILE *text_matrix;	//CHHS For debugging
#endif
