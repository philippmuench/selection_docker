//CHHS global.cpp contains definitions of global variables (with initializations), cf. global.h
#include "defines.h"
//typedef enum {PrBranch=1, PrNodeNum=2, PrLabel=4, PrAge=8, PrOmega=16} OutTreeOptions_typedef; //CHHS see above
typedef enum
{ LOWER_F = 1, UPPER_F, BOUND_F } FOSSIL_FLAGS_typedef;	//CHHS see above
typedef enum
{ Fequal, F1x4, F3x4, Fcodon, F1x4MG, F3x4MG, FMutSel0, FMutSel } CodonFreqs_typedef;	//CHHS see above
typedef enum
{ NSbranchB = 1, NSbranch2, NSbranch3 } NSBranchModels_typedef;	//CHHS see above
typedef enum
{ Poisson, EqualInput, Empirical, Empirical_F, FromCodon = 6, REVaa_0 = 8, REVaa = 9 } AAModel_typedef;	//CHHS see above
typedef enum
{ NSnneutral = 1, NSpselection, NSdiscrete, NSfreqs, NSgamma, NS2gamma, NSbeta, NSbetaw, NSbetagamma, NSbeta1gamma, NSbeta1normal, NS02normal, NS3normal } NSsitesModels_typedef;	//CHHS see above
typedef enum
{ FIT1 = 11, FIT2 = 12 } SiteClassModels_typedef;	//CHHS see above
typedef enum
{ AAClasses = 7 } aaDistModels_typedef;	//CHHS see above
typedef enum
{ GlobalClock = 1, LocalClock, ClockCombined } ClockModels_typedef;	//CHHS see above

//DataTypes_typedef DataTypes; //CHHS enum substituted in favor of typedef, see global.h
//OutTreeOptions_typedef OutTreeOptions; //CHHS enum substituted in favor of typedef, see global.h
char BASEs[] = "TCAGUYRMKSWHBVD-N?";
char AAs[] = "ARNDCQEGHILKMFPSTWYV-*?X";
char CODONs[256][4];
const char *EquateBASE[] = {
  "T", "C", "A", "G", "T", "TC", "AG", "CA", "TG", "CG", "TA",
  "TCA", "TCG", "CAG", "TAG", "TCAG", "TCAG", "TCAG"
};				//CMV added const

char nChara[256], CharaMap[256][64];
char AA3Str[] =
  { "AlaArgAsnAspCysGlnGluGlyHisIleLeuLysMetPheProSerThrTrpTyrVal***" };
char BINs[] = "TC";
int GeneticCode[][64] = {
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 0:universal */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, -1, -1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 1:vertebrate mt. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
   16, 16, 16, 16, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 2:yeast mt. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 3:mold mt. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 15, 15,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 4:invertebrate mt. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, 5, 5, 4, 4, -1, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 5:ciliate nuclear */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 2, 11, 15, 15, 15, 15,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 6:echinoderm mt. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 4, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 7:euplotid mt. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17,
   10, 10, 10, 15, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},
  /* 8:alternative yeast nu. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 7, 7,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 9:ascidian mt. */
  {
   13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, 5, 4, 4, -1, 17,
   10, 10, 10, 10, 14, 14, 14, 14, 8, 8, 5, 5, 1, 1, 1, 1,
   9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1,
   19, 19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},	/* 10:blepharisma nu. */
  {
   1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
   5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8,
   9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12,
   13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16}	/* 11:Ziheng's regular code */
};				/* GeneticCode[icode][#codon] */

int noisy = 0, Iround = 0, NFunCall = 0, NEigenQ, NPMatUVRoot;
double SIZEp = 0;		//CHHS This variable is also used as local in tools.c!
int AlwaysCenter = 0;
double Small_Diff = 1e-6;	/* reasonable values 1e-5, 1e-7 */
double prob_Quantile, *par_Quantile;	//CHHS dropped "static" property (should have no effect)
double (*cdf_Quantile) (double x, double par[]);	//CHHS This is a function pointer, property "static" dropped
unsigned int z_rndu = 1237;	//CHHS dropped "static" property (should have no effect)
int w_rndu = 1237;		//CHHS dropped "static" property (should have no effect)
time_t time_start;		//CHHS dropped "static" property (should have no effect)
struct TREEB
{
  int nbranch, nnode, root, branches[NBRANCH][2];
  double lnL;
} tree;
struct common_info
{
  char *z[NS], *spname[NS], seqf[96], outf[96], treef[96], daafile[96], cleandata;	//CHHS *z[NS] lost "unsigned" property to unify pointer data types
  char oldconP[NNODE];		/* update conP for nodes? to save computation */
  int seqtype, ns, ls, ngene, posG[NGENE + 1], lgene[NGENE], npatt,
    *pose, readpattern;
  int runmode, clock, verbose, print, codonf, aaDist, model, NSsites;
  int nOmega, nbtype, nOmegaType;	/* branch partition, AA pair (w) partition */
  int method, icode, ncode, Mgene, ndata, bootstrap;
  int fix_rgene, fix_kappa, fix_omega, fix_alpha, fix_rho, nparK,
    fix_blength, getSE;
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
struct TREEN
{
  int father, nson, sons[MAXNSONS], ibranch, ipop;
  double branch, age, omega, *conP, label;
  char *nodeStr, fossil /*, usefossil */ ;	//CMV unused
} *nodes, **gnodes, nodes_t[2 * NS - 1];
/* for sptree.nodes[].fossil: lower, upper, bounds, gamma, inverse-gamma */
FOSSIL_FLAGS_typedef FOSSIL_FLAGS;	//CHHS enum substituted in favor of typedef, see global.h
const char *fossils[] = { " ", "L", "U", "B" };	//CMV added const

struct SPECIESTREE
{
  int nbranch, nnode, root, nspecies, nfossil;
  struct TREESPN
  {
    char name[LSPNAME + 1], fossil, usefossil;	/* fossil: 0, 1, 2, 3 */
    int father, nson, sons[2];
    double age, pfossil[7];	/* lower and upper bounds or alpha & beta */
    //double *lnrates;      /* log rates for loci */ //CMV unused
  } nodes[2 * NS - 1];
} sptree;
/* all trees are binary & rooted, with ancestors unknown. */
struct DATA			/* locus-specific data and tree information */
{
  int ns[NGENE], ls[NGENE], npatt[NGENE], ngene /*, lgene[NGENE] */ ;	//CMV unused
  int root[NGENE + 1], /*BlengthMethod, */ fix_nu, nbrate[NGENE],	//CMV unused
    icode[NGENE];
  char *z[NGENE][NS], cleandata[NGENE];
  char /*idaafile[NGENE], */ daafile[NGENE][40];	//CMV unused
  double *fpatt[NGENE] /*, lnpT, lnpR, lnpDi[NGENE] */ ;	//CMV unused
  double /*Qfactor[NGENE], */ pi[NGENE][NCODE];	//CMV unused
  double /*rgene[NGENE], */ kappa[NGENE], alpha[NGENE], omega[NGENE];	//CMV unused
  int NnodeScale[NGENE];
  char *nodeScale[NGENE];	/* nScale[data.ns[locus]-1] for interior nodes */
} data;
int Nsensecodon, FROM61[64], FROM64[64], FourFold[4][4];
int ChangedInIteration;		/* 1: t changed, update P(t); 2: paras changed, update UVRoot */
double *PMat, *U, *V, *Root, *_UU[NBTYPE + 2], *_VV[NBTYPE + 2],
  *_Root[NBTYPE + 2];
/* 5 sets for branchsite models (YN2002); 6 sets for clade models */
double pcodon0[64], paa0[20], *pcodonClass;	/* for aaDist=FIT1 or FIT2 */
int BayesEB;			/* =1 for site models M2a & M8; =2 for branch-site models A & C */
int LASTROUND;
int IClass = -1;
int OmegaAA[190], AA1STEP[190];
double _rateSite = 1;
double Qfactor_NS, Qfactor_NS_branch[NBTYPE];
double AAchem[][20 + 1] =	/* last element is the max */
{
  {
   8.1, 10.5, 11.6, 13, 5.5, 10.5, 12.3, 9, 10.4, 5.2,
   4.9, 11.3, 5.7, 5.2, 8, 9.2, 8.6, 5.4, 6.2, 5.9, 13},	/* p */
  {
   31, 124, 56, 54, 55, 85, 83, 3, 96, 111,
   111, 119, 105, 132, 32.5, 32, 61, 170, 136, 84, 170},	/* v */
  {
   0, 0.65, 1.33, 1.38, 2.75, 0.89, 0.92, 0.74, 0.58,
   0, 0, 0.33, 0, 0, 0.39, 1.42, 0.71, 0.13, 0.2, 0, -999},	/* c */
  {
   -0.11, 0.079, -0.136, -0.285, -0.184, -0.067, -0.246, -0.073, 0.32,
   0.001,
   -0.008, 0.049, -0.041, 0.438, -0.016, -0.153, -0.208, 0.493, 0.381, -0.155}	/* a */
};				/* in the order p, v, c, a */

FILE *fout, *frub, *flnf, *frst, *frst1, *frst2 = NULL, *finitials;
const char *ratef = "rates";	//CMV added const
CodonFreqs_typedef CodonFreqs;	//CHHS enum substituted in favor of typedef, see global.h
const char *codonfreqs[] = {
  "Fequal", "F1x4", "F3x4", "Fcodon", "F1x4MG", "F3x4MG", "FMutSel0",
  "FMutSel"
};				//CMV added const

NSBranchModels_typedef NSBranchModels;	//CHHS enum substituted in favor of typedef, see global.h
const char *NSbranchmodels[] = { "One dN/dS ratio",
  "free dN/dS Ratios for branches", "several dN/dS ratios for branches",
  "NSbranch3"
};				//CMV added const

AAModel_typedef AAModel;	//CHHS enum substituted in favor of typedef, see global.h
const char *aamodels[] = {
  "Poisson", "EqualInput", "Empirical", "Empirical_F", "",
  "", "FromCodon", "", "REVaa_0", "REVaa"
};				//CMV added const

NSsitesModels_typedef NSsitesModels;	//CHHS enum substituted in favor of typedef, see global.h
const char *NSsitesmodels[] = {
  "one-ratio", "NearlyNeutral", "PositiveSelection", "discrete",
  "freqs",
  "gamma", "2gamma", "beta", "beta&w>1", "beta&gamma", "beta&gamma+1",
  "beta&normal>1", "0&2normal>0", "3normal>0"
};				//CMV added const

SiteClassModels_typedef SiteClassModels;	//CHHS enum substituted in favor of typedef, see global.h
aaDistModels_typedef aaDistModels;	//CHHS enum substituted in favor of typedef, see global.h
const char *clockstr[] = { "", "Global clock", "Local clock", "ClockCombined" };	//CMV added const

ClockModels_typedef ClockModels;	//CHHS enum substituted in favor of typedef, see global.h
/* variables for batch run of site models */
int ncatG0 = 10, insmodel = 0, nnsmodels = 1, nsmodels[14] = { 0 };

/* used for sliding windows analysis */
int windowsize0 = 20, offset0 = 1, npositive = 0;
double lnLmodel;
int _nestS = 0;			/* 189= estimate the S elements, 0= use those from com.daa[] */
double *_Fij;			//CHHS static property removed (should have no effect)
int ijAAref = 19 * 20 + 9;	//CHHS static propert removed (should have no effect)
//CHHS The following variables moved from treesub.c
double *dfsites;		//CHHS lost "static" property
#ifdef LSDISTANCE		//CHHS We keep the #ifdef, although not very effective for defining single variables
double *SeqDistance = NULL;
int *ancestor = NULL;
#endif
#if(defined(BASEML) || defined(CODEML))	//CHHS We keep the #ifdef, although not very effective for defining single variables
double *AgeLow = NULL;
int NFossils = 0, AbsoluteRate = 0;
double ScaleTimes_TipDate = 1, TipDate = 0;
/* number of internal node times, usd to deal with known ancestors.  Broken? */
int innode_time = 0;		//CHHS "static" property dropped (should have no effect)
#endif
int ONE = 1;			//CHHS Mainly for BLAS and LAPACK function calls
double DONE = 1.0;		//CHHS For BLAS
int ZERO = 0;			//CHHS For BLAS
double DZERO = 0.0;		//CHHS For BLAS
FILE *text_matrix;		//CHHS For debugging
