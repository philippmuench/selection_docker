//CHHS functions.h includes funtion headers
#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H 1		//CHHS Safeguard to avoid multiple includes of the same things
int ReadSeq (FILE * fout, FILE * fseq, int cleandata);
int ScanFastaFile (FILE * f, int *ns, int *ls, int *aligned);
void SetMapAmbiguity (void);
void ReadPatternFreq (FILE * fout, char *fpattf);
int Initialize (FILE * fout);
int MoveCodonSeq (int ns, int ls, char *z[]);
int PatternWeight (void);
int PatternWeightJC69like (FILE * fout);
int PatternWeightSimple (int CollapsJC);
void f_and_x_fromx (double x[], double f[], int n, int LastItem);	//CMV f_and_x is always called with fromf == 0
void SetSeed (unsigned int seed);
double rndu (void);
void rndu_vector (double r[], int n);
int SampleCat (double P[], int n, double space[]);
int AutodGamma (double Mmat[], double freqK[], double rK[], double *rho1,
		double alfa, double rho, int K);
int DiscreteGamma (double freqK[], double rK[], double alpha, double beta,
		   int K, int dgammamean);
double CDFBeta (double x, double p, double q, double lnbeta);
double QuantileBeta (double prob, double p, double q, double lnbeta);
double Quantile (double (*cdf) (double x, double par[]), double p, double x,
		 double par[], double xb[2]);
double CDFNormal (double x);
double LnGamma (double alpha);
double DFGamma (double x, double alpha, double beta);
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
double logLBinormal (double h, double k, double r);
void rainbowRGB (double temperature, int *R, int *G, int *B);
void GetIndexTernary (int *ix, int *iy, double *x, double *y, int itriangle,
		      int K);
int CodeChara (char b, int seqtype);
void printSeqs (FILE * fout, int *pose, char keep[], int format);
int printsma (FILE * fout, char *spname[], char *z[],	//CHHS unsigned dropped for second and third parameter to unify data types of pointers
	      int ns, int l, int lline, int gap, int seqtype,
	      int transformed, int simple, int pose[]);
int testXMat (double x[]);
double SeqDivergence (double x[], int model, double alpha, double *kapa);
int symtest (double freq[], int n, int nobs, double space[], double *chisym,
	     double *chihom);
int dSdNNG1986 (char *z1, char *z2, int lc, int icode, int transfed,
		double *dS, double *dN, double *Ssites, double *Nsites);
int difcodonNG (char codon1[], char codon2[], double *SynSite,
		double *AsynSite, double *SynDif, double *AsynDif,
		int transfed, int icode);
int testTransP (double P[], int n);
#if (defined(BASEML))
int PMatK80 (double P[], double t, double kapa);
#endif
int PMatTN93 (double P[], double a1t, double a2t, double bt, double pi[]);
void PMatUVRoot2 (double P[], double t, int n, double U[], double V[], double Root[]);	//CHHS Successor of above
int PMatCijk (double PMat[], double t);
int PMatQRev (double P[], double pi[], double t, int n, double space[]);
double DistanceIJ (int is, int js, int model, double alpha, double *kappa);
int DistanceMatNuc (FILE * fout, FILE * f2base, int model, double alpha);
int EigenQREVbase (FILE * fout, double kappa[], double pi[],
		   int *nR, double Root[], double Cijk[]);
int DistanceMatNG86 (FILE * fout, FILE * fds, FILE * fdn, FILE * dt,
		     double alpha);
int setmark_61_64 (void);
int rell (FILE * flnf, FILE * fout, int ntree);
int MultipleGenes (FILE * fout, FILE * fpair[], double space[]);
int lfunRates (FILE * fout, double x[], int np);
int AncestralSeqs (FILE * fout, double x[]);
int NucListall (char b, int *nb, int ib[4]);
char *getcodon (char codon[], int icodon);
char *getAAstr (char *AAstr, int iaa);
int Codon2AA (char codon[3], char aa[3], int icode, int *iaa);
int DNA2protein (char dna[], char protein[], int lc, int icode);
int printcu (FILE * f1, double fcodon[], int icode);
int printcums (FILE * fout, int ns, double fcodons[], int code);
#ifdef BASEML
int QtoPi (double Q[], double pi[], int n, double *space);
#endif
int PtoPi (double P[], double pi[], int n, double *space);
void starttimer (void);
char *printtime (char timestr[]);
char *strc (int n, int c);
int printdouble (FILE * fout, double a);
void strcase (char *str, int direction);
void error2 (const char *message);	//CMV made argument const
int indexing (double x[], int n, int index[], int descending, int space[]);
FILE *gfopen (const char *filename, const char *mode);	//CMV made arguments const
int appendfile (FILE * fout, const char *filename);	//CMV made filename const
int matby (double a[], double b[], double c[], int n, int m, int k);
int matIout (FILE * fout, int x[], int n, int m);
int matout (FILE * file, double x[], int n, int m);
int matout2 (FILE * fout, double x[], int n, int m, int wid, int deci);
int matinv (double x[], int n, int m, double space[]);
#ifdef BASEML
int matexp (double Q[], double t, int n, int TimeSquare, double space[]);
#endif
void eigenQREV (double Q[], double pi[], int n, double Root[], double U[], double V[], double spacesqrtpi[]);	//CMV changed to void
int comparedouble (const void *a, const void *b);
int splitline (char line[], int fields[]);
int scanfile (FILE * fin, int *nrecords, int *nx, int *ReadHeader,
	      char line[], int ifields[]);
double bound (int nx, double x0[], double p[], double x[],
	      int (*testx) (double x[], int nx));
int gradient (int n, double x[], double f0, double g[],
	      double (*fun) (double x[], int n), double space[], int Central);
int Hessian (int nx, double x[], double f, double g[], double H[],
	     double (*fun) (double x[], int n), double space[]);
int HessianSKT2004 (double xmle[], double lnLm, double g[], double H[]);
double LineSearch (double (*fun) (double x), double *f, double *x0,
		   double xb[2], double step, double e);
void xtoFreq (double x[], double freq[], int n);
int SetxBound (int np, double xb[][2]);
int ming2 (FILE * fout, double *f, double (*fun) (double x[], int n),
	   int (*dfun) (double x[], double *f, double dx[], int n),
	   double x[], double xb[][2], double space[], double e, int n);
int minB (FILE * fout, double *lnL, double x[], double xb[][2], double e,
	  double space[]);
int minB2 (FILE * fout, double *lnL, double x[], double xb[][2], double e,
	   double space[]);
int nls2 (FILE * fout, double *sx, double *x0, int nx,
	  int (*fun) (double x[], double y[], int nx, int ny),
	  int (*jacobi) (double x[], double J[], int nx, int ny),
	  int (*testx) (double x[], int nx), int ny, double e);
void NodeToBranch (void);
void BranchToNode (void);
int ReadTreeN (FILE * ftree, int *haslength, int *haslabel, int copyname,
	       int popline);
int OutTreeN (FILE * fout, int spnames, int printopt);
int OutTreeB (FILE * fout);
void PointconPnodes (void);
int SetBranch (double x[]);
int DistanceMat (FILE * fout, int ischeme, double alfa, double *kapa);
int StepwiseAddition (FILE * fout, double space[]);
int readx (double x[], int *fromfile);
int PopEmptyLines (FILE * fseq, int lline, char line[]);
int blankline (char *str);
void BranchLengthBD (int rooted, double birth, double death, double sample,
		     double mut);
int RandomLHistory (int rooted, double space[]);
int RootTN93 (int ischeme, double kapa1, double kapa2, double pi[],
	      double *scalefactor, double Root[]);
int EigenTN93 (int ischeme, double kapa1, double kapa2, double pi[],
	       int *nR, double Root[], double Cijk[]);
double MPScore (double space[]);
int MakeTreeIb (int ns, int Ib[], int rooted);
int NumberTrees (int ns, int rooted);
int CountLHistories (void);
void ReRootTree (int newroot);
int NeighborNNI (int i_tree);
int CountLHistory (char LHistories[], double space[]);
int GetSubSeqs (int nsnew);
int GenerateSeq (void);
void Evolve (int inode);
void EvolveJC (int inode);
int GetGtree (int locus);
int Forestry (FILE * fout);
int GetMemPUVR (int nc, int nUVR);
int testx (double x[], int np);
int SetxBound (int np, double xb[][2]);
int SetxInitials (int np, double x[], double xb[][2]);
int GetInitials (double x[], int *fromfile);
double *PointOmega (double xcom[], int igene, int inode, int isiteclass);
int SetParameters (double x[]);
int Set_UVR_BranchSite (int iclass, int branchlabel);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int SetPSiteClass (int iclass, double x[]);
int PMatJC69like (double P[], double t, int n);
int InitializeCodon (FILE * fout, double space[]);
int GetDaa (FILE * fout, double daa[]);
void EigenQcodon (int getstats, double blength, double *S, double *dS, double *dN, double Root[], double U[], double V[], double *meanrate, double kappa[], double omega, double Q[]);	//CMV made void
int EigenQaa (FILE * fout, double Root[], double U[], double V[],
	      double rate[]);
int TestModelQc (FILE * fout, double x[]);
int PairwiseCodon (FILE * fout, FILE * fds, FILE * fdn, FILE * dt,
		   double space[]);
int PairwiseAA (FILE * fout, FILE * f2AA);
double GetBranchRate (int igene, int ibrate, double x[], int *ix);
int GetPMatBranch (double Pt[], double x[], double t, int inode);
void ConditionalPNode (int inode, int igene, double x[]);
char GetAASiteSpecies (int species, int sitepatt);
void finishup (void);
int mergeSeqs (FILE * fout);
void Get4foldSites (void);
int AdHocRateSmoothing (FILE * fout, double x[NS * 3], double xb[NS * 3][2],
			double space[]);
void DatingHeteroData (FILE * fout);
void SimulateData2s61 (void);
void Ina (void);
void get_grid_para_like_M2M8 (double para[4][100], int n1d, int dim, int M2a,
			      int ternary, double p0b[], double p1b[],
			      double w0b[], double wsb[], double p_beta_b[],
			      double q_beta_b[], double x[], double *S);
void GetIndexTernary (int *ix, int *iy, double *x, double *y, int itriangle,
		      int K);
void get_grid_para_like_AC (double para[][100], int n1d, int dim,
			    double w0b[], double w2b[], double x[],
			    double *S);
double lfunAdG (double x[], int np);	//CHHS There was no signature found anywhere, added one here
double lfundG (double x[], int np);	//CHHS There was no signature found anywhere, added one here
double lfun (double x[], int np);	//CHHS There was no signature found anywhere, added one here
int InitializeBaseAA (FILE * fout);	//CHHS There was no signature found anywhere, added one here
int Perturbation (FILE * fout, int initialMP, double space[]);
int StarDecomposition (FILE * fout, double space[]);	//CHHS There was no signature found anywhere, added one here
void FreeMemPUVR (void);	//CHHS There was no signature found anywhere, added one here
int GetTreeFileType (FILE * ftree, int *ntree, int *pauptree, int shortform);
int OutputTimesRates (FILE * fout, double x[], double var[]);	//CHHS There was no signature found anywhere, added one here
int SetxBoundTimes (double xb[][2]);	//CHHS There was no signature found anywhere, added one here
int GetInitialsTimes (double x[]);	//CHHS There was no signature found anywhere, added one here
void InitializeNodeScale (void);	//CHHS There was no signature found anywhere, added one here
void NodeScale (int inode, int pos0, int pos1);	//CMV changed to void
void fx_r (double x[], int np);
void SaveToOctave (double *CVariable, char *OctaveVariable, FILE * FilePointer, int Rows, int Columns);	//CHHS new function to help comparing results with GNU Octave
//CMV Added file for all inline functions (like zero() etc.)
#include "inline.h"
#endif //CHHS Safeguard
