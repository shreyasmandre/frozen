#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <cblas.h>
#include <math.h>
long double scaleeig=1.0;

void dpttrf_(int *N, double *d, double *e, int *info);
void dpttrs_(int *N, int *nrhs, double *d, double *e, double *b, int *ldb, int *info);
void dgtsv_(int *N, int *nrhs, double *dl, double *d, double *du, double *vortnew, int *ldb, int *info);

inline double Czfunc(double x, double t) {return exp(-x*x/4/t)/sqrt(M_PI*t);};
double max(double a, double b) { if(a>b) return a; else return b; } ;
double min(double b, double a) { if(a<b) return a; else return b; } ;
int C2P(int N, double *p, double *c, double deta, double kx, double perm)
{
    int i;
    static double *d;
    static double *e;
    static int init=0, initkx=0;
    double tmp = 2/pow(deta,2) + kx*kx;
    double tmp2 = 1/pow(deta,2);
    //double zero = 0, mone=-1;
    static double kxold=-1;
	int info;
    kx = fabs(kx);
    if(!init)
    {
	d = (double*)malloc(N*sizeof(double));
	e = (double*)malloc((N-1)*sizeof(double));
	init = 1;
    }
    {
	for(i=0;i<N;i++) d[i] = tmp;
	for(i=0;i<N-1;i++) e[i] = -tmp2;
    d[0] = tmp2 + kx*perm/deta;
	//d[0]   = tmp2 + kx/deta; // -psi1/deta**2 + psi0  /deta**2 + kx psi0/deta = 0
	d[N-1] = tmp2 + (1e-4+fabs(kx))/deta; //  psin/deta**2 - psinm1/deta**2 + kx psin/deta = 0
    if(perm>1e5)
    {
        e[0] = 0; d[0] = 1;
    }
	dpttrf_(&N, d, e, &info);
	kxold = kx;
	initkx = 1;
	//fprintf(stderr, "Error code = %d\n", info);
    }
    //for(i=0;i<N-1;i++) printf("%lf\t%lf+I%lf\n", d[i], creal(e[i]), cimag(e[i]));
    for(i=1;i<N-1;i++) p[i] = -(c[i]-c[i-1])/deta;
    p[0] = 0; p[N-1] = 0;
    {
	int ldrhs=N, info, nrhs=1;
        dpttrs_(&N, &nrhs, d, e, p, &ldrhs, &info);
    }
    return 0;
}
int ConcStep(double ds, double s, double *conc, double *p, double *pnew, 
   double *concnew, const double *z, double dz, double kx, int N)
{
    static int init=0;
    static double *dl, *d, *du;
    static double *Cz, *Cznew;
    double nu=1;
    double lambda = nu*ds/dz/dz;
    int i;
    if(!init)
    {
	d  = (double*)malloc(N*sizeof(double));
	dl = (double*)malloc((N-1)*sizeof(double));
	du = (double*)malloc((N-1)*sizeof(double));
	Cz     = (double*)malloc(N*sizeof(double));
	Cznew  = (double*)malloc(N*sizeof(double));
	init = 1;
    }
    for(i=0;i<N;i++)
    {
        Cz[i]       = Czfunc(z[i], s);
        Cznew[i]    = Czfunc(z[i], s+ds);
    }

    for(i=0;i<N;i++) d[i] = 1 + lambda + nu*ds*kx*kx/2;
    for(i=0;i<N-1;i++) dl[i] = -lambda/2 ;
    for(i=0;i<N-1;i++) du[i] = -lambda/2 ;
    for(i=1;i<N-1;i++) concnew[i] = conc[i]*(1 - lambda - nu*ds*kx*kx/2 + ds*Cz[i])
	    - ds*( (p[i]-p[i-1])*Cz[i] )/dz
        +lambda/2*(conc[i+1]+conc[i-1]);
    d[0] = 1; du[0]= -1e-6; concnew[0] = 0;
    d[N-1] = 1+fabs(kx); dl[N-2] = -1; concnew[N-1] = 0;
    {
	int nrhs=1, ldb=N, info=0;
        dgtsv_(&N, &nrhs, dl, d, du, concnew, &ldb, &info);
    }
    return 0;
}
void step(double dt, double t, double *c, double *p, double *pnew, double *cnew, 
          const double *z, double dz, double kx, double perm, int N)
{
    C2P(N-1, p, c, dz, kx, perm);
    ConcStep(dt, t, c, p, p, cnew, z, dz, kx, N);
    //C2P(N-1, pnew, cnew, dz, kx, perm);
    //ConcStep(dt, t, c, p, pnew, cnew, z, dz, kx, N);
    return;
}

double errorcontrol(double ds, double s, double *c, double *p, double *pnew, 
	     double *cnew, const double *z, double dz, double kx, double perm, int N)
{
    static int init=0;
    static double *ptmp1, *ctmp1, *ptmp2, *ctmp2;
    double res;
    // double den;
    int i;
    if(!init)
    {
        ptmp1 = (double*)malloc(N*sizeof(double));
        ctmp1 = (double*)malloc(N*sizeof(double));
        ptmp2 = (double*)malloc(N*sizeof(double));
        ctmp2 = (double*)malloc(N*sizeof(double));
	init = 1;
    }
    //den = cblas_dnrm2(N, c, 1);
    //for(i=0;i<N;i++) c[i] /= den;
    //for(i=0;i<N;i++) p[i] /= den;
    step(ds,   s, c, p, ptmp1, ctmp1, z, dz, kx, perm, N);
    step(ds/2, s, c, p, ptmp2, ctmp2, z, dz, kx, perm,N);
    step(ds/2, s, ctmp2, ptmp2, pnew, cnew, z, dz, kx, perm, N);
    for(i=0,res=0;i<N;i++) res += pow(fabs(cnew[i]-ctmp1[i]),2);
    //for(i=0;i<N;i++) c[i] *= den;
    //for(i=0;i<N;i++) p[i] *= den;
    //for(i=0;i<N;i++) cnew[i] *= den;
    //for(i=0;i<N;i++) pnew[i] *= den;
    return sqrt(res);
}

double MatVec(const double *in, double *out, const double t1, const double t2, const double kx, double perm, const double *z, double dz, const int N)
{
    static double *p, *c, *pnew, *cnew;
    static int init=0;
    double t, dt, dtsuggest, derror, error=0, tolerance = 5e-6;
    long double factor=1;
    double den;
    int i, stepno;
    if(!init)
    {
        p = (double*)malloc(N*sizeof(double));
        c = (double*)malloc(N*sizeof(double));
        pnew = (double*)malloc(N*sizeof(double));
        cnew = (double*)malloc(N*sizeof(double));
        init = 1;
    }
    for(i=0;i<N;i++) c[i] = in[i];
    t = t1;
    dt = 1e-4;
    dtsuggest = dt;
    factor = 1;
    for(stepno=0;t<t2;stepno++)
    {
        den = cblas_dnrm2(N, c, 1);
        for(i=0;i<N;i++) c[i] /= den;
        for(i=0;i<N;i++) p[i] /= den;
        factor *=den;
        derror = 1;
        while(derror>tolerance)
        {
            dt = dtsuggest;
            derror = errorcontrol(dt, t, c, p, pnew, cnew, z, dz, kx, perm, N);
            dtsuggest = min(0.8*dt*pow(tolerance/derror, 0.5), 100000);
            dtsuggest = min(dtsuggest, fabs(t2-t));
        }
        error += derror;
        t += dt;
        //if(stepno%1000==0)
        //{ 
        //    fprintf(stderr, "t = %le, dt = %le, res = %le\r", t, dt, derror);
        //}
        dtsuggest = min(dtsuggest, fabs(t2-t));
        {double *tmp; tmp = pnew; pnew = p; p = tmp;}
        {double *tmp; tmp = cnew; cnew = c; c = tmp;}
    }
    for(i=0;i<N;i++) out[i] = c[i]*(factor/scaleeig);
    return error;
}

void dnaupd_(int *ido, char *bmat, int *N, char *which, int *nev, double *tol, double *resid, int *ncv, double *v,
             int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);
// resid has length N
// v has size ncv*N
// iparam has length 11
// ipntr has length 14
// workd has length 3N
// workl has length lworlk > 3ncv**2+6ncv

void dneupd_(int *rvec, char *, int *select, double *dr, double *di, double *v, int *ldv, double *sr, double *si, double *workev,
             char *bmat, int *N, char *which, int *nev, double *tol, double *resid, int *ncv, double *, int *, 
             int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);
// workev has size 3*ncv

double Arnoldi(const double s1, const double s2, const double kx, double perm, double *y, double dy, const int N, double *Evec)
{
    static double *resid, *v, *workd, *workl, *dr, *di, *workev;
    //static double *rwork;
    int ido, nev, ncv, ldv=N, iparam[11], ipntr[14], lworkl, info, ierr, N2=N;
    double tol, error;
    int rvec=1, counter;
    static int *select;
    static int init=0;
    int done = 0;
    double sigmar, sigmai, retval;
    int i, nmax=0;
    nev = 5; ncv = 11; tol = 1e-4; ldv = N; lworkl = 2*(3*ncv*ncv+6*ncv); ido = 0;
    iparam[0] = 1; iparam[2] = 5000; iparam[3] = 1; iparam[4] = nev; iparam[6] = 1; 
    info = 0;
    if(!init)
    {
         resid = (double*)malloc(N*sizeof(double));
         v     = (double*)malloc(N*ncv*sizeof(double));
         dr    = (double*)malloc(ncv*sizeof(double));
         di    = (double*)malloc(ncv*sizeof(double));
         workd = (double*)malloc(3*N*sizeof(double));
         workl = (double*)malloc(lworkl*sizeof(double));
         workev= (double*)malloc(3*ncv*sizeof(double));
         select = (int*)malloc(ncv*sizeof(int));
         init = 1;
    }
    {int i; for(i=0;i<N;i++) resid[i] = Evec[i];}
    info = 1;
    counter = 0;
    while(!done)
    {
         counter++;
         dnaupd_(&ido, "I", &N2, "LM", &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
         //if(ido!=1 && ido!=-1) done = 1;
         if(counter%1==0) fprintf(stderr,"iter=%d, ido=%d\r", counter, ido);
         if(ido==99) done = 1;
         else 
         {
             // call av (nx, workd(ipntr(1)), workd(ipntr(2)))
             error = MatVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]), s1, s2, kx, perm, y, dy, N);
         }
    }
    fprintf(stderr, "\ninfo = %d\n", info);
    if(info<0)
    {
            fprintf(stderr," Error in dnaupd_: info = %d ... quiting?\n", info);
            exit(1);
    }
    dneupd_(&rvec, "A", select, dr, di, v, &ldv, &sigmar, &sigmai, workev, "I", &N2, "LM", &nev, &tol, resid, &ncv, v, &ldv, 
            iparam, ipntr, workd, workl, &lworkl, &ierr);
    fprintf(stderr,"out of dneupd: ierr=%d\n", ierr);
    //for(i=0;i<nev;i++) fprintf(stderr, "%le+i%le\n", dr[i], di[i]);
    retval = 0;
    for(i=0;i<nev;i++) 
    {
        double val = sqrt(dr[i]*dr[i]+di[i]*di[i]);
        if(val>retval)
        {
            nmax = i;
            retval = val;
        }
    }
    for(i=0;i<N;i++) Evec[i] = v[i+N*nmax];
    return retval;
}
#define EIGROUT Arnoldi

void NewtonEvalue(double s1, double s2, double kx1, double kx2, double perm, double *y, double dy, int N, double *kxmax, double *Evaluemax, double *Evec)
{
    double kx3, kxnew;
    double eval1, eval2, eval3;
    double f1, f2;
    eval1 = EIGROUT(s1, s2, kx1, perm, y, dy, N, Evec);
    eval2 = EIGROUT(s1, s2, kx2, perm, y, dy, N, Evec);
    kx3 = 0.5*(kx1+kx2);                                  
    eval3 = EIGROUT(s1, s2, kx3, perm, y, dy, N, Evec);
    fprintf(stderr, "k = (%le, %le, %le)\n", kx1, kx3, kx2);
    fprintf(stderr, "f = (%le, %le, %le)\n", eval1, eval3, eval2);
    while(fabs(kx2-kx1)>1e-5)
    {
		f1 = (eval2-eval1)/(kx2-kx1);
		f2 = (eval2 + eval1 - 2*eval3)/pow(0.5*(kx2-kx1),2);
		kxnew = kx3 - f1/f2;
		kx1 = kx3; eval1 = eval3;
		kx2 = kx1 + 2*(kxnew-kx1); eval2 = EIGROUT(s1, s2, kx2, perm, y, dy, N, Evec);
        kx3 = kxnew;               eval3 = EIGROUT(s1, s2, kx3, perm, y, dy, N, Evec);	
		fprintf(stderr, "k = (%le, %le, %le)\n", kx1, kx3, kx2);
		fprintf(stderr, "f = (%le, %le, %le)\n", eval1, eval3, eval2);
		*kxmax = kx3; *Evaluemax = eval3;
    }
    return;
}


void MaximizeEvalue(double s1, double s2, double kx1, double kx2, double perm, double *y, double dy, int N, double *kxmax, double *Evaluemax, double *Evec)
{
    //static complex double *Evec;
    static int init=0;
    int iter;
    double kxm1, kxm2;
    double eval1=0, eval2=0, evalm1=0, evalm2=0;
    double golden_ratio = (3-sqrt(5))/2;
    int done;
    if(!init)
    {
	init = 1;
    }
    //for(i=0;i<N;i++) Evalues[i] = 0;
    kxm1 = kx1 + (kx2-kx1)*golden_ratio;
    kxm2 = kx2 - (kx2-kx1)*golden_ratio;
    eval1 = EIGROUT(s1, s2, kx1, perm, y, dy, N, Evec);
    fprintf(stderr, "kx1=%le, kxm1=%le, kxm2=%le, kx2=%le\n", kx1, kxm1, kxm2, kx2);
    fprintf(stderr, "ex1=%le, exm1=%le, exm2=%le, ex2=%le\n", (eval1), (evalm1), (evalm2), (eval2));
    eval2 = EIGROUT(s1, s2, kx2, perm, y, dy, N, Evec);
    fprintf(stderr, "kx1=%le, kxm1=%le, kxm2=%le, kx2=%le\n", kx1, kxm1, kxm2, kx2);
    fprintf(stderr, "ex1=%le, exm1=%le, exm2=%le, ex2=%le\n", (eval1), (evalm1), (evalm2), (eval2));
    evalm1 = EIGROUT(s1, s2, kxm1, perm, y, dy, N, Evec);
    fprintf(stderr, "kx1=%le, kxm1=%le, kxm2=%le, kx2=%le\n", kx1, kxm1, kxm2, kx2);
    fprintf(stderr, "ex1=%le, exm1=%le, exm2=%le, ex2=%le\n", (eval1), (evalm1), (evalm2), (eval2));
    evalm2 = EIGROUT(s1, s2, kxm2, perm, y, dy, N, Evec);
    fprintf(stderr, "kx1=%le, kxm1=%le, kxm2=%le, kx2=%le\n", kx1, kxm1, kxm2, kx2);
    fprintf(stderr, "ex1=%le, exm1=%le, exm2=%le, ex2=%le\n", (eval1), (evalm1), (evalm2), (eval2));
    done = 0;
    for(iter=0;!done;iter++)
    {
	//fprintf(stderr, "kx1 = %le, eval1=%le\n", kx1, eval1);
	//fprintf(stderr, "kx2 = %le, eval2=%le\n", kx2, eval2);
	//fprintf(stderr, "kxm1 = %le, evalm1=%le\n", kxm1, evalm1);
	//fprintf(stderr, "kxm2 = %le, evalm2=%le\n", kxm2, evalm2);
	//fprintf(stderr, "... done iteration %d.\n", iter);
	if( ((evalm1)>(evalm2) && (evalm1)>(eval1)))
	{
	    //fprintf(stderr, "Maximum between kx1=%le and kxm2=%le\n", kx1, kxm2);
	    kx2 = kxm2;
	    eval2 = evalm2;
	    kxm2 = kxm1;
	    evalm2 = evalm1;
	    kxm1 = kx1 + (kx2-kx1)*golden_ratio;
        evalm1 = EIGROUT(s1, s2, kxm1, perm, y, dy, N, Evec);
	}
    else if((evalm2)>(evalm1) && (evalm2)>(eval2))
	{
            //fprintf(stderr, "Maximum between kxm1=%le and kx2=%le.\n", kxm1, kx2);
	    kx1 = kxm1;
	    eval1 = evalm1;
	    kxm1 = kxm2;
	    evalm1 = evalm2;
        kxm2 = kx2 - (kx2-kx1)*golden_ratio;
        evalm2 = EIGROUT(s1, s2, kxm2, perm, y, dy, N, Evec);
	}
	else if((eval1)>(evalm1))
	{
	    //fprintf(stderr, "Maximum at kx1=%le.\n", kx1);
	    kx2 = kxm2;
	    eval2 = evalm2;
	    kxm2 = kxm1;
	    evalm2 = evalm1;
	    kxm1 = kx1 + (kx2-kx1)*golden_ratio;
        evalm1 = EIGROUT(s1, s2, kxm1, perm, y, dy, N, Evec);
	}
    else
    {
        //fprintf(stderr, "Maximum at kx2=%le.\n", kx2);
	    kx1 = kxm1;
	    eval1 = evalm1;
	    kxm1 = kxm2;
	    evalm1 = evalm2;
        kxm2 = kx2 - (kx2-kx1)*golden_ratio;
        evalm2 = EIGROUT(s1, s2, kxm2, perm, y, dy, N, Evec);
	}

        if(fabs(kx2-kx1)/max(0.05, fabs(kx1+kx2))<5e-2) done = 1;
        fprintf(stderr, "kx1=%le, kxm1=%le, kxm2=%le, kx2=%le\n", kx1, kxm1, kxm2, kx2);
        fprintf(stderr, "ex1=%le, exm1=%le, exm2=%le, ex2=%le\n", (eval1), (evalm1), (evalm2), (eval2));
    }
    {
        double f2, f1, xx1, xx2, xx3, ff1, ff2, ff3;
        xx1 = 0.5*(kx1+kxm1); xx2 = 0.5*(kxm1+kxm2); xx3 = 0.5*(kxm2+kx2);
        ff1 = 0.5*(eval1+evalm1); ff2 = 0.5*(evalm1+evalm2); ff3 = 0.5*(evalm2+eval2);

        f2 = (ff3+ff1-2*ff2)/pow(xx3-xx2,2);
        f1 = (ff3-ff1)/(xx3-xx1);
        *kxmax = (xx1+xx3)/2 - f1/f2;
        *Evaluemax = ff2;
    }
    // NewtonEvalue(s1, s2, kx1, kx2, perm, y, dy, N, kxmax, Evaluemax, Evec);
	// kxmax = 0.5*(kx1+kx2); 
    // Evaluemax = 
}
#define NS 19
int main(int argc, char **argv)
{
    int N=4801;
    int i,j;
    double L=600, dz, kx=2e-3, t1, t2, kx1, kx2;
    double res;
    double *z;
    double *p, *c;
    double perm=0.5;
    double T[NS];
    FILE *growth, *vec;
    char vecname[100];
    double tstart = 1, tend = 100;
    double ntfactor = round(log(tend/tstart)/log(10.0));
    if(argc>1) sscanf(argv[1], "%le", &perm);
    if(argc>2) sscanf(argv[2], "%le", &tstart);
    if(argc>3) sscanf(argv[3], "%le", &tend);
    ntfactor = round(log(tend/tstart)/log(10.0));
    fprintf(stderr, "perm = %le, tstart = %le, tend = %le\n", perm, tstart, tend);
    ntfactor = (log(tend/tstart)/log(10.0));
    z = (double*)malloc(N*sizeof(double));
    p = (double*)malloc((N-1)*sizeof(double));
    c = (double*)malloc(N*sizeof(double));
    dz = L*1.0/(N-1);
    for(j=0;j<N;j++) z[j] = -j*dz;
    for(j=0;j<N;j++) c[j] = 0;
    c[1] = 1;
    
    //for(j=0;j<NS;j++) T[j] = 1e6*pow((j+1)*1.0/(NS),2);
    //for(j=0;j<NS;j++) T[j] = tstart*pow(10, (j)*ntfactor/(NS-1));
    for(j=0;j<NS;j++) T[j] = tstart+(tend-tstart)*j*1.0/(NS-1);
    growth = fopen("growth.dat", "w");
    for(i=0;i<NS;i++)
    {
        kx = 2e-3;
        t1 = T[i];
        for(j=0;j<=i;j++)
        {
			int k;
            t2 = T[j];
			for(k=0;k<21; k++)
			{
				 kx = 1e-4 + k*(0.1-1e-4)/20; 
				 fprintf(growth, "%le\t%le\t%le\t%le\n", t1, t2, kx, 0.0); fflush(growth);
			}
			fprintf(growth, "\n");
        }
        kx1 = -0.9e-1; kx2 = 1e-1;
        // kx1 = kx-0.0; kx2 = kx+0.05;
        scaleeig = 1;
        for(j=i+1;j<NS;j++)
        {
            int k;
            t2 = T[j];
            fprintf(stderr,"\ni=%d, j=%d, t1=%le, t2=%le\n",i, j, t1, t2);
            for(k=0;k<N;k++) c[k] = 0;
            c[1] = 1;

			/*
            //res = Arnoldi(t1, t2, kx, perm, z, dz, N, c);
            MaximizeEvalue(t1, t2, kx1, kx2, perm, z, dz, N, &kx, &res, c);
            //NewtonEvalue(t1, t2, kx1, kx2, perm, z, dz, N, &kx, &res, c);
            C2P(N-1, p, c, dz, kx, perm);
            scaleeig *= res;
            kx1 = -9e-2; kx2 = +1e-1;
            //kx1 = -1e-1; kx2 = 1e-1;

            fprintf(stderr, "\nAmplification = %Le at k=%le\n", scaleeig, kx);
            fprintf(growth, "%le\t%le\t%le\t%Le\n", t1, t2, kx, logl(scaleeig)); fflush(growth);
            
            sprintf(vecname, "vecs/vec%d_%d.dat", i, j);
            vec = fopen(vecname, "w");
            for(k=0;k<N;k++) fprintf(vec, "%le\t%le\t%le\n", z[k], c[k], p[k]);
            fclose(vec);
			*/
			for(k=0;k<21; k++)
			{
				int k2;
				kx = 1e-4 + k*(0.1-1e-4)/20; 
				res = EIGROUT(t1, t2, kx, perm, z, dz, N, c);
				C2P(N, p, c, dz, kx, perm);
					// scaleeig *= res;

				 fprintf(stderr, "\nAmplification = %le at k=%le\n", res, kx);
				 fprintf(growth, "%le\t%le\t%le\t%le\n", t1, t2, kx, log(res)); fflush(growth);
				 sprintf(vecname, "vecs/vec%d_%d_%d.dat", i, j, k);
				 vec = fopen(vecname, "w");
				 fprintf(vec, "# perm = %le, tstart = %le, tend = %le\n", perm, t1, t2);
				 for(k2=0;k2<N;k2++) fprintf(vec, "%le\t%le\t%le\n", z[k2], c[k2], p[k2]);
				 fclose(vec);
			}
			fprintf(growth, "\n");
        }
        fprintf(growth, "\n");
    }
    //for(j=0;j<N;j++) printf("%le\t%le\n", z[j], c[j]);
    return 0;
}
