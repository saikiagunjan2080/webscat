package coreservlets;
public class nonSphericalClass {
public static final int npn1= 100;
public static final int npng1= 300;
public static final int npng2 =2*npng1;
public static final int npn2=2*npn1;
public static final int npl =npn2+1;
public static final int npn3 =npn1+1;
public static final int npn4= 80;
public static final int npn5= 2*npn4;
public static final int npn6= npn4+1;
public static final int npl1= npn5+1;
static int ichoice;
static double ss11[]=new double[181];
static double ss12[]=new double[181];
static double ss33[]=new double[181];
static double ss34[]=new double[181];
static double trgqr[][] = new double[npn2][npn2];
static double trgqi[][] = new double[npn2][npn2];
static double qr[][] = new double[npn2][npn2];
static double qi[][] = new double[npn2][npn2];
static double rgqr[][] = new double[npn2][npn2];
static double rgqi[][] = new double[npn2][npn2];
static double ssign[]=new double[900];
static double aa,bb,a,b;
static double ppi;
static double pir;
static double pii;
static int lmax;
static double rat;
static double reff;
static double veff;
static double cond;
static double csca;
static double cextin;
static double cscat;
static double asymm;
static double cabsin;
static double walb;
public static String nonsp(double imrr,double imri,int indistr,double iaxi,
double ir1,double ir2 ,double ib,double ieps,double iddelt,int nkmax,int np,double lam) {
String  ret="";
int i;
int l;
int m;
int n;
int l1;
int m1;
int n1;
int n2;
int n11;
int ii;
int n22;
int nk;
int nm;
int l1m;
int nm1;
int nn1;
int nn2;
int nma;
int iax;
int ink;
int nnm;
int inm1;
int nggg;
int ndgs;
int npna;
int nmin;
int mmax;
int nmax=0;
int ixxx;
int l1max;
int nmax1=0;
int itime;
int ngaus;
int npnax;
int ncheck;
int nnnggg;
int ngauss=0;
int ndistr;
double p;
double r1;
double r2;
double z1;
double z2;
double z3;
double zz1;
double zz2;
double zz3;
double zz4;
double zz5;
double zz6;
double zz7;
double zz8;
double gam;
double dax;
double axi;
double mri;
double eps;
double mrr;
double xev;
double qsc;
double wgi;
double qxt;
double dsca;
double dn1;
double qsca;
double time;
double wgii;
double dext;
double cext;
double qext;
double qsca1;
double ti1nn;
double qext1;
double tr1nn;
double ddelt;
double dqsca;
double ti1nn1;
double axmax;
double tr1nn1;
double dqext;
double coeff1;
double wg[]=new double [1000];
double xg[]=new double [1000];
double wg1[]=new double [2000];
double xg1[]=new double [2000];
double r[]=new double[npng2] ;
double s[]=new double[npng2];
double w[]=new double[npng2];
double x[]=new double[npng2];
double an[]=new double[npn1];
double dr[]=new double[npng2] ;
double ss[]=new double[npng2] ;
double be1[]=new double[npl] ;
double be2[]=new double[npl] ;
double al1[]=new double[npl] ;
double al2[]=new double[npl] ;
double al3[]=new double[npl] ;
double al4[]=new double[npl] ;
double ddr[]=new double[npng2] ;
double dri[]=new double[npng2] ;
double drr[]=new double[npng2] ;
double bet1[]=new double[npl] ;
double bet2[]=new double[npl] ;
double alph1[]=new double[npl] ;
double alph2[]=new double[npl] ;
double alph3[]=new double[npl] ;
double alph4[]=new double[npl] ;
double ann[][]=new double[npn1][npn1];

nonSphericalClass nonSphericalClass=new nonSphericalClass();
mrr = imrr;
mri =imri;
ndistr=indistr;
axi = iaxi;
r1=ir1;
r2=ir2;
b = ib;
eps =ieps;
p=Math.acos(-1);
rat=0.5;
npnax = 1;
gam = .5;
ddelt = iddelt;
npna = 181;
ndgs = 2;
ichoice = 2;
ncheck = 0;
if (np == -1 || np == -2) 
{
ncheck = 1;
}
if (np > 0 && Math.pow((-1),np) == 1) 
{
ncheck = 1;
}
if ((Math.abs(rat - 1.0)) > 0.00000001 && np == -1) 
{
nonSphericalClass.sarea(eps);
}
if ((Math.abs(rat - 1.0)) > 0.00000001 && np == -2) 
{
nonSphericalClass.sareac(eps);
}
ddelt = 0.1*ddelt;
tmatBlock tmatblock=new tmatBlock(npng2,npn1);
trtiBlock trtiblock=new trtiBlock(npn2);
rirgigBlock rirgigblock=new rirgigBlock(npn1);
cbassBlock cbassblock=new cbassBlock(npng2,npn1);
nk = (int) (nkmax + 2);
if (nk > 1000) 
{
ret="NK , I.E. IS GREATER THAN 1000 : EXECUTION TERMINATED";
return ret;
}
nonSphericalClass.gauss(nk,0,0, xg, wg);
z1 = (r2 - r1) * 0.5;
z2 = (r1 + r2) * 0.5;
z3 = r1 * 0.5;
for (i= 1; i<=nk; i++) 
{
xg1[i-1] = z1 * xg[i-1] + z2;
wg1[i-1] = wg[i-1] * z1;
}
nonSphericalClass.distrb(nk, xg1, wg1, ndistr, axi, b, gam, r1, r2, p);
for (i=1;i<=npl;i++) 
{
alph1[i - 1] = 0.;
alph2[i - 1] = 0.;
alph3[i - 1] = 0.;
alph4[i - 1] = 0.;
bet1[i - 1] = 0.;
bet2[i - 1] = 0.;
}
cscat = 0.;
cextin = 0.;
l1max = 0;
for (ink = 1; ink <= nk; ink++) {
i = nk - ink + 1;
a = rat * xg1[i - 1];
xev = p * 2. * a / lam;
ixxx =(int)(xev + Math.pow(xev,0.333333)*4.05);
inm1 = Math.max(4,ixxx);
if (inm1 >= npn1) 
{
ret="CONVERGENCE IS NOT OBTAINED FOR NPN1  : EXECUTION TERMINATED"+npn1;
return ret;
}
qext1 = 0.0;
qsca1 = 0.0;
for (nma = inm1; nma <= npn1; nma++) 
{
nmax = nma;
mmax = 1;
ngauss = nmax * ndgs;
if (ngauss > npng1) 
{
ret="NGAUSS =  I.E. IS GREATER THAN NPNG1 : EXECUTION TERMINATED";
return ret;
}
nonSphericalClass.constt(ngauss, nmax, mmax, p, x, w, an, ann, s, ss, np, eps);
ret=nonSphericalClass.vary(lam, mrr, mri, a, eps, np, ngauss, x, p, r, dr, ddr, drr, dri, nmax,cbassblock);
if(null!=ret&&!ret.equalsIgnoreCase(""))
{
return ret;
}
nonSphericalClass.tmatr0(ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, ddr, drr, dri, nmax, ncheck,tmatblock,trtiblock,rirgigblock,cbassblock);
qext = 0.0;
qsca = 0.0;
for (n = 1; n <= nmax; n++) 
{
n1 = n + nmax;
tr1nn =trtiblock.tr1[n-1][n-1];
ti1nn = trtiblock.ti1[n-1][n-1];
tr1nn1 = trtiblock.tr1[n1-1][n1-1];
ti1nn1 = trtiblock.ti1[n1-1][n1-1];
dn1 = (double) (2*n + 1);
qsca = qsca + dn1 * (tr1nn * tr1nn + ti1nn * ti1nn + tr1nn1 * tr1nn1 + ti1nn1 * ti1nn1);
qext = qext + (tr1nn + tr1nn1) * dn1;
}
dsca = (Math.abs((qsca1 - qsca) / qsca));
dext = (Math.abs((qext1 - qext) / qext));
qext1 = qext;
qsca1 = qsca;
nmin = (int) ((double) nmax / 2.0 + 1.);
for (n = nmin; n <= nmax; ++n) 
{
n1 = n + nmax;
tr1nn = trtiblock.tr1[n-1][n-1];
ti1nn = trtiblock.ti1[n-1][n-1];
tr1nn1 = trtiblock.tr1[n1-1][n1-1];
ti1nn1 = trtiblock.ti1[n1-1][n1-1];
dn1 = (double) (2*n + 1);
dqsca = dn1 * (tr1nn * tr1nn + ti1nn * ti1nn + tr1nn1 * tr1nn1 + ti1nn1 * ti1nn1);
dqext = (tr1nn + tr1nn1) * dn1;
dqsca = (Math.abs(dqsca / qsca));
dqext = (Math.abs(dqext / qext));
nmax1 = n;
if (dqsca <= ddelt && dqext <= ddelt) 
{
break;
}
}
if (dsca <= ddelt && dext <= ddelt) 
{
break;
}
if (nma == 100) 
{
ret="CONVERGENCE IS NOT OBTAINED FOR NPN1 =  : EXECUTION TERMINATED \n"+npn1;
return ret;
}
}
nnnggg = ngauss + 1;
if (ngauss == npng1) 
{
System.out.println("WARNING:NGAUSS=NPNG1");
}
mmax = nmax1;
for (ngaus = nnnggg; ngaus <= npng1; ngaus++)
{
ngauss = ngaus;
nggg = 2*ngauss;
nonSphericalClass.constt(ngauss, nmax, mmax, p, x, w, an, ann, s, ss, np, eps);
ret=nonSphericalClass.vary(lam, mrr, mri, a, eps, np, ngauss, x, p, r, dr, ddr, drr, dri, nmax,cbassblock);
if(null!=ret&&!ret.equalsIgnoreCase(""))
{
return ret;
}
nonSphericalClass.tmatr0(ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r,dr, ddr, drr, dri, nmax, ncheck,tmatblock,trtiblock,rirgigblock,cbassblock);
qext = 0.0;
qsca = 0.0;
for (n = 1; n <= nmax; ++n) 
{
n1 = n + nmax;
tr1nn =trtiblock. tr1[n-1][n-1];
ti1nn = trtiblock.ti1[n-1][n-1];
tr1nn1 = trtiblock.tr1[n1-1][n1-1];
ti1nn1 = trtiblock.ti1[n1-1][n1-1];
dn1 = (double) (2*n + 1);
qsca = qsca + dn1 * (tr1nn * tr1nn + ti1nn * ti1nn + tr1nn1 * tr1nn1 + ti1nn1 * ti1nn1);
qext = qext + (tr1nn + tr1nn1) * dn1;
}
dsca = (Math.abs((qsca1 - qsca) / qsca));
dext = (Math.abs((qext1 - qext) / qext));
qext1 = qext;
qsca1 = qsca;
if (dsca <= ddelt && dext <= ddelt) 
{
break;
}
}
qsca = 0.0;
qext = 0.0;
nnm = 2*nmax;
for (n = 1; n <=nnm; ++n) 
{
qext = qext + trtiblock.tr1[n-1][n-1];
}
if (nmax1 > npn4) 
{
ret="NMAX IS GREATER THAN NPN4: EXECUTION TERMINATED";
return ret;
}
gspBlock gspblock=new gspBlock(npn4,npn6);
for (n2 = 1; n2 <= nmax1; n2++) 
{
nn2 = n2 + nmax;
for (n1 = 1; n1 <= nmax1; n1++) 
{
nn1 = n1 + nmax;
zz1 = trtiblock.tr1[n1-1][n2-1];
gspblock.rt11[0][n1-1][n2-1] = zz1;
zz2 = trtiblock.ti1[n1-1][n2-1];
gspblock.it11[0][n1-1][n2-1]= zz2;
zz3 = trtiblock.tr1[n1-1][nn2-1];
gspblock.rt12[0][n1-1][n2-1] = zz3;
zz4 = trtiblock.ti1[n1-1][nn2-1];
gspblock.it12[0][n1-1][n2-1] = zz4;
zz5 =trtiblock. tr1[nn1-1][n2-1];
gspblock.rt21[0][n1-1][n2-1] = zz5;
zz6 = trtiblock.ti1[nn1-1][n2-1];
gspblock.it21[0][n1-1][n2-1] = zz6;
zz7 = trtiblock.tr1[nn1-1][nn2-1];
gspblock.rt22[0][n1-1][n2-1] = zz7;
zz8 = trtiblock.ti1[nn1-1][nn2-1];
gspblock.it22[0][n1-1][n2-1] = zz8;
qsca = qsca + zz1 * zz1 + zz2 * zz2 + zz3 * zz3 + zz4 * zz4 + zz5 * zz5 + zz6 * zz6 + zz7 * zz7 + zz8 * zz8;
}
}
for (m = 1; m <= nmax1; ++m) 
{
nonSphericalClass.tmatr(m, ngauss, x, w, an, ann, s, ss, ppi, pir, pii, r, dr, ddr, drr, dri, nmax, ncheck,tmatblock,trtiblock,rirgigblock,cbassblock);
nm = nmax - m + 1;
nm1 = nmax1 - m + 1;
m1 = m + 1;
qsc = 0.0;
for (n2 = 1; n2 <= nm1; n2++) 
{
nn2 = n2 + m - 1;
n22 = n2 + nm;
for (n1 = 1; n1 <= nm1; ++n1) 
{
nn1 = n1 + m - 1;
n11 = n1 + nm;
zz1 = trtiblock.tr1[n1-1][n2-1];
gspblock.rt11[m1-1][nn1-1][nn2-1] = zz1;
zz2 = trtiblock.ti1[n1-1][n2-1];
gspblock.it11[m1-1][nn1-1][nn2-1] = zz2;
zz3 = trtiblock.tr1[n1-1][n22-1];
gspblock.rt12[m1-1][nn1-1][nn2-1] = zz3;
zz4 =trtiblock. ti1[n1-1][n22-1];
gspblock.it12[m1-1][nn1-1][nn2-1] = zz4;
zz5 =trtiblock. tr1[n11-1][n2-1];
gspblock.rt21[m1-1][nn1-1][nn2-1] = zz5;
zz6 = trtiblock.ti1[n11-1][n2-1];
gspblock.it21[m1-1][nn1-1][nn2-1] = zz6;
zz7 =trtiblock. tr1[n11-1][n22-1];
gspblock.rt22[m1-1][nn1-1][nn2-1] = zz7;
zz8 =trtiblock. ti1[n11-1][n22-1];
gspblock.it22[m1-1][nn1-1][nn2-1] = zz8;
qsc = qsc+ (zz1 * zz1 + zz2 * zz2 + zz3 * zz3 + zz4 * zz4 + zz5 * zz5 + zz6 * zz6 + zz7 * zz7 + zz8 *
zz8) * 2.;
}
}
nnm = 2*nm;
qxt = 0.0;
for (n = 1; n <= nnm; n++) 
{
qxt = qxt+trtiblock.tr1[n-1][n-1] * 2.0;
}
qsca = qsca+qsc;
qext = qext+qxt;
}
coeff1 = lam * lam * 0.5 / p;
csca = qsca * coeff1;
cext = -qext * coeff1;
nonSphericalClass.gsp(nmax1, csca, lam, al1, al2, al3, al4, be1, be2,gspblock);
l1m = lmax + 1;
l1max = Math.max(l1max,l1m);
wgii = wg1[i - 1];
wgi = wgii * csca;
for (l1 = 1; l1 <= l1m; ++l1) 
{
alph1[l1 - 1] = alph1[l1 - 1]+al1[l1 - 1] * wgi;
alph2[l1 - 1] = alph2[l1 - 1]+al2[l1 - 1] * wgi;
alph3[l1 - 1] = alph3[l1 - 1]+al3[l1 - 1] * wgi;
alph4[l1 - 1] = alph4[l1 - 1]+al4[l1 - 1] * wgi;
bet1[l1 - 1] = bet1[l1 - 1]+be1[l1 - 1] * wgi;
bet2[l1 - 1] = bet2[l1 - 1]+be2[l1 - 1] * wgi;
}
cscat = cscat+wgi;
cextin = cextin+cext * wgii;
}
for (l1 = 1; l1 <= l1max; ++l1) 
{
alph1[l1 - 1] = alph1[l1 - 1]/cscat;
alph2[l1 - 1] = alph2[l1 - 1]/cscat;
alph3[l1 - 1] = alph3[l1 - 1]/cscat;
alph4[l1 - 1] = alph4[l1 - 1]/cscat;
bet1[l1 - 1] = bet1[l1 - 1]/cscat;
bet2[l1 - 1] = bet2[l1 - 1]/cscat;
}
walb = cscat / cextin;
cabsin=cextin-cscat;
nonSphericalClass.hovenr(l1max, alph1, alph2, alph3, alph4, bet1, bet2);
asymm = alph1[1] / 3.0;
lmax = l1max - 1;
nonSphericalClass.matr(alph1, alph2, alph3, alph4, bet1, bet2, lmax, npna);
return ret;
}

public int constt(int ngauss, int nmax, int mmax,double p, double x[], double w[], double an[],double ann[][], double s[],
  double ss[], int np,double eps)
{
int i;
int n;
int n1;
int ng;
int nn;
int ng1;
int ng2;
double d;
double y;
double dd[]=new double[100];
double xx;
double ddd;
double w1[]=new double[npng2];
double x1[]=new double[npng2];
double x2[]=new double[npng2];
double w2[]=new double[npng2];
boolean flagIn=false;
for (n = 1; n <= nmax; ++n) 
{
nn = n * (n + 1);
an[n-1] = (double) nn;
d = Math.sqrt((double) (2*n + 1) / (double) nn);
dd[n - 1] = d;
for (n1 = 1; n1 <= n; ++n1) 
{
ddd = d * dd[n1 - 1] * 0.5;
ann[n-1][n1-1] = ddd;
ann[n1-1][n-1] = ddd;
}
}
ng = 2*(ngauss);
if (np == -2) 
{
flagIn=true;
}

if(flagIn)
{
ng1 = (int) ((double) (ngauss) / 2.0);
ng2 = ngauss - ng1;
xx = -Math.cos(Math.atan(eps));
gauss(ng1, 0, 0, x1, w1);
gauss(ng2, 0, 0, x2, w2);
for (i = 1; i <= ng1; ++i) 
{
w[i-1] = (xx + 1.0) * 0.5 * w1[i - 1];
x[i-1] = (xx + 1.0) * 0.5 * x1[i - 1] + (xx - 1.0) * 0.5;
}
for (i= 1; i<=ng2; ++i) 
{
w[i+ ng1-1] = xx * (-0.5) * w2[i - 1];
x[i+ ng1-1] = xx * (-0.5) * x2[i - 1] + xx * 0.5;
}
for (i = 1; i <= ngauss; ++i) 
{
w[ng - i] = w[i-1];
x[ng - i] = -x[i-1];
}
for (i = 1; i <= ngauss; ++i) 
{
y = x[i-1];
y = 1.0 / (1.0 - y * y);
ss[i-1] = y;
ss[ng - i] = y;
y = Math.sqrt(y);
s[i-1] = y;
s[ng - i] = y;
}
}
else
{
gauss(ng, 0, 0, x, w);//////check/////
for (i = 1; i <= ngauss; ++i) 
{
y = x[i-1];
y = 1.0 / (1.0 - y * y);
ss[i-1] = y;
ss[ng - i] = y;
y = Math.sqrt(y);
s[i-1] = y;
s[ng - i] = y;
}
}
return 0;
} 

public String vary(double lam, double mrr, double mri, double a, double eps, int np, int ngauss, double x[], double p, 
double r[], double dr[], double ddr[], double drr[], double dri[], int nmax,cbassBlock cbassblock)
{
int i;
int ng;
int nnmax1;
int  nnmax2;
double v;
double  v1;
double  v2;
double ta;
double  tb;
double  pi;
double  vv;
double  pri;
double  prr;
double z[]=new double[npng2] ;
double zi[]=new double[npng2] ;
double zr[]=new double[npng2] ;
String strret="";
ng = 2*(ngauss);
if (np == -1) 
{
rsp1(x, ng, ngauss, a, eps, np, r, dr);
}
if (np == -2) 
{
rsp3(x, ng, ngauss, a, eps, r, dr);
}
pi = p * 2. / lam;
ppi = Math.PI * Math.PI;
pir = ppi * mrr;
pii = ppi * mri;
v = 1.0 / (mrr * mrr + mri * mri);
prr = mrr * v;
pri = -(mri) * v;
ta = 0.;
for (i = 1; i <= ng; ++i) 
{
vv = Math.sqrt(r[i-1]);
v = vv * pi;
ta = Math.max(ta,v);
vv = 1.0 / v;
ddr[i-1] = vv;
drr[i-1] = prr * vv;
dri[i-1] = pri * vv;
v1 = v * mrr;
v2 = v * mri;
z[i - 1] = v;
zr[i - 1] = v1;
zi[i - 1] = v2;
}
if (nmax > 100) 
{
strret="NMAX = "+nmax+", I.E. GREATER THAN NPN1 = "+npn1;
return strret;
}
tb = ta * Math.sqrt(mrr * mrr + mri * mri);
tb = Math.max(tb,(double) (nmax));
nnmax1 = (int) (Math.sqrt((Math.max(ta,(double) (nmax)))) * 1.2 + 3.0);
nnmax2 = (int) (tb +Math.pow(tb, 0.33333) * 4.0 + Math.sqrt(tb) * 1.2);
nnmax2 = nnmax2 - nmax + 5;
bess(z, zr, zi, ng, nmax, nnmax1, nnmax2,cbassblock);
return strret;
}

public int rsp1(double x[], int ng, int ngauss,double rev, double eps, int np, double r[],double dr[])
{
int i;
double a;
double c;
double s;
double aa;
double cc;
double ee;
double rr;
double ss;
double ee1;

a = rev * Math.pow(eps, 1.0/3.0);
aa = a * a;
ee = eps * eps;
ee1 = ee - 1.0;
for (i = 1; i <= ngauss; ++i) 
{
c = x[i-1];
cc = c * c;
ss = 1.0 - cc;
s = Math.sqrt(ss);
rr = 1.0 / (ss + ee * cc);
r[i-1] = aa * rr;
r[ng - i] = r[i-1];
dr[i-1] = rr * c * s * ee1;
dr[ng - i] = -dr[i-1];
}
return 0;
}

public int rsp3(double x[], int ng, int ngauss,double rev, double eps, double r[], double dr[])
{
int i;
double a;
double h;
double co;
double si;
double rad;
double rthet;

h = rev * Math.pow((2.0 / (eps * 3.0 * eps)), 1.0/3.0);
a = h * eps;
for (i = 1; i <= ngauss; ++i) 
{
co = -x[i-1];
si = Math.sqrt(1.0 - co * co);
if (si / co > a / h) 
{
rad = a / si;
rthet = -a * co / (si * si);
r[i-1] = rad * rad;
r[ng - i] = r[i-1];
dr[i-1] = -rthet / rad;
dr[ng - i] = -dr[i-1];
}
else
{
rad = h / co;
rthet = h * si / (co * co);
r[i-1] = rad * rad;
r[ng - i] = r[i-1];
dr[i-1] = -rthet / rad;
dr[ng - i] = -dr[i-1];
}}
return 0;
}

public int bess(double x[], double xr[], double xi[],int ng, int nmax, int nnmax1, int nnmax2,cbassBlock cbassblock)
{
int i;
int l;
int n;
double yi;
double yr;
double xx;
double aj[]=new double[npn1];
double ay[]=new double[npn1];
double adj[]=new double[npn1];
double aji[]=new double[npn1];
double ajr[]=new double[npn1];
double ady[]=new double[npn1];
double adji[]=new double[npn1];
double adjr[]=new double[npn1];

for (i = 1; i <= ng; ++i) 
{
xx = x[i-1];
rjb(xx, aj, adj, nmax, nnmax1);
ryb(xx, ay, ady, nmax);
yr = xr[i-1];
yi = xi[i-1];
cjb(yr, yi, ajr, aji, adjr, adji, nmax, nnmax2);
for (n = 1; n <= nmax; ++n) 
{
cbassblock.j[i-1][n-1] = aj[n - 1];
cbassblock.y[i-1][n-1] = ay[n - 1];
cbassblock.jr[i-1][n-1] = ajr[n - 1];
cbassblock.ji[i-1][n-1] = aji[n - 1];

cbassblock.dj[i-1][n-1] = adj[n - 1];
cbassblock.dy[i-1][n-1] = ady[n - 1];
cbassblock.djr[i-1][n-1] = adjr[n - 1];
cbassblock.dji[i-1][n-1] = adji[n - 1];
}
}

return 0;
}

public int rjb(double x, double y[], double u[],int nmax, int nnmax)
{
int i;
int l;
int i1;
int l1;
double z[]=new double [800];
double y0;
double z0;
double y1;
double yi;
double xx;
double yi1;

l = nmax + nnmax;
xx = 1.0 / x;
z[l - 1] = 1.0 / ((double) (2*l + 1) * xx);
l1 = l - 1;
for (i = 1; i<= l1; ++i) 
{
i1 = l - i;
z[i1 - 1] = 1.0 / ((double) (2*i1 + 1) * xx - z[i1]);
}
z0 = 1.0 / (xx - z[0]);
y0 = z0 * Math.cos(x) * xx;
y1 = y0 * z[0];
u[0] = y0 - y1 * xx;
y[0] = y1;
for (i = 2; i <= nmax; ++i) 
{
yi1 = y[i - 2];
yi = yi1 * z[i - 1];
u[i-1] = yi1 - (double) i * yi * xx;
y[i-1] = yi;
}
return 0;
}

public int ryb(double x, double y[], double v[], int nmax)
{
int i;
int nmax1;
double c;
double s;
double x1;
double x2;
double x3;
double y1;

c = Math.cos(x);
s = Math.sin(x);
x1 = 1.0 / x;
x2 = x1 * x1;
x3 = x2 * x1;
y1 = -c * x2 - s * x1;
y[0] = y1;
y[1] = (x3 * -3.0 + x1) * c - x2 * 3.0 * s;
nmax1 = nmax - 1;
for (i= 2; i<= nmax1; ++i) 
{
y[i] = (double) (2*i + 1) * x1 * y[i-1] - y[i - 2];
}
v[0] = -x1 * (c + y1);
for (i= 2; i <= nmax; ++i) 
{
v[i-1] = y[i - 2] - (double) i * x1 * y[i-1];
}
return 0;
}

public int cjb(double xr, double xi, double yr[], double yi[], double ur[], double ui[], int nmax, int nnmax)
{
int i;
int l;
int i1;
int l1;
double ai;
double  ci;
double  ar;
double  cr;
double  qf;
double  qi;
double  ari;
double  cui[]=new double [npn1];
double  cyi[]=new double [npn1];
double  czi[]=new double [1200];
double  cur[]=new double [npn1];
double  cyr[]=new double [npn1];
double  czr[]=new double [1200];
double  cu1i;
double  cy0i;
double  cz0i;
double  cy1i;
double  cu1r;
double  cz0r;
double  cy0r;
double  cy1r;
double  cuii;
double  cyii;
double  cuir;
double  cyir;
double  cxxi;
double  cxxr;
double  xrxi;
double cyi1i;
double  cyi1r;

l = nmax + nnmax;
xrxi = 1.0 / (xr * xr + xi * xi);
cxxr = xr * xrxi;
cxxi = -(xi) * xrxi;
qf = 1.0 / (double) (2*l + 1);
czr[l - 1] = xr * qf;
czi[l - 1] = xi * qf;
l1 = l - 1;
for (i = 1; i <=l1; ++i) 
{
i1 = l - i;
qf = (double) (2*i1 + 1);
ar = qf * cxxr - czr[i1];
ai = qf * cxxi - czi[i1];
ari = 1.0 / (ar * ar + ai * ai);
czr[i1 - 1] = ar * ari;
czi[i1 - 1] = -ai * ari;
}
ar = cxxr - czr[0];
ai = cxxi - czi[0];
ari = 1.0 / (ar * ar + ai * ai);
cz0r = ar * ari;
cz0i = -ai * ari;
cr = Math.cos(xr) * Math.cosh(xi);
ci = -Math.sin(xr) * Math.sinh(xi);
ar = cz0r * cr - cz0i * ci;
ai = cz0i * cr + cz0r * ci;
cy0r = ar * cxxr - ai * cxxi;
cy0i = ai * cxxr + ar * cxxi;
cy1r = cy0r * czr[0] - cy0i * czi[0];
cy1i = cy0i * czr[0] + cy0r * czi[0];
ar = cy1r * cxxr - cy1i * cxxi;
ai = cy1i * cxxr + cy1r * cxxi;
cu1r = cy0r - ar;
cu1i = cy0i - ai;
cyr[0] = cy1r;
cyi[0] = cy1i;
cur[0] = cu1r;
cui[0] = cu1i;
yr[0] = cy1r;
yi[0] = cy1i;
ur[0] = cu1r;
ui[0] = cu1i;
for (i = 2; i <= nmax; ++i) 
{
qi = (double) i;
cyi1r = cyr[i - 2];
cyi1i = cyi[i - 2];
cyir = cyi1r * czr[i - 1] - cyi1i * czi[i - 1];
cyii = cyi1i * czr[i - 1] + cyi1r * czi[i - 1];
ar = cyir * cxxr - cyii * cxxi;
ai = cyii * cxxr + cyir * cxxi;
cuir = cyi1r - qi * ar;
cuii = cyi1i - qi * ai;
cyr[i - 1] = cyir;
cyi[i - 1] = cyii;
cur[i - 1] = cuir;
cui[i - 1] = cuii;
yr[i-1] = cyir;
yi[i-1] = cyii;
ur[i-1] = cuir;
ui[i-1] = cuii;
}
return 0;
}

public int tmatr0(int ngauss, double x[], double w[],double an[], double ann[][], double s[], double ss[],double ppi, double pir,
  double pii, double r[],double dr[], double ddr[], double drr[], double dri[], int nmax, int ncheck,tmatBlock tmatblock,trtiBlock trtiblock,rirgigBlock rirgigblock,
  cbassBlock cbassblock)
{
int i;
int  n;
int i1;
int  i2;
int  k1;
int  k2;
int  n1;
int  n2;
int ng;
int nm;
int mm1;
int kk1;
int  kk2;
int ngss;
int nnmax;
double f1;
double  f2;
double a12;
double  a21;
double  a22;
double si;
double rr[]=new double[npng2];
double aa1;
double  dd1;
double  dd2;
double  b1i;
double  c1i;
double  c2i;
double  b2i;
double  an1;
double  an2;
double  c3i;
double  b3i;
double c4i;
double  b1r;
double  c1r;
double  c2r;
double  b2r;
double  c3r;
double  b3r;
double  dv1[]=new double[npn1];
double  dv2[]=new double[npn2];
double qj1;
double  c4r;
double  b4r;
double  b4i;
double  c5r;
double  c5i;
double  b5r;
double  b5i;
double qy1;
double  ai12;
double  ai21;
double  ar12;
double  ar21;
double  gi12;
double  sig[]=new double[npn2];
double  gr12;
double  gr21;
double  gi21;
double  d1n1;
double  d2n1;
double  d1n2;
double  d2n2;
double  uri;
double  rri;
double  an12;
double  qdj1;
double  qji2;
double  qjr2;
double  qdy1;
double  tai12;
double  tai21;
double  ddri;
double  tgi12;
double  drii;
double  tgi21;
double  tar12;
double  tar21;
double  tgr12;
double  drri;
double  tgr21;
double  tpii;
double tppi;
double  tpir;
double  qdji2;
double  qdjr2;
double  factor;

mm1 = 1;
nnmax = nmax + nmax;
ng = 2*(ngauss);
ngss = ng;
factor = 1.0;
if (ncheck == 1) 
{
ngss = ngauss;
factor = 2.0;
} 
else 
{
//CONTINUE
}
si = 1.0;
for (n = 1; n <= nnmax; ++n) 
{
si = -si;
sig[n - 1] = si;
}
for (i = 1; i<= ngauss; ++i) 
{
i1 = ngauss + i;
i2 = ngauss - i + 1;
vig(x[i1-1], nmax, 0, dv1, dv2);
for (n = 1; n <=nmax; ++n) 
{
si = sig[n - 1];
dd1 = dv1[n - 1];
dd2 = dv2[n - 1];
tmatblock.d1[i1-1][n-1] = dd1;
tmatblock.d2[i1-1][n-1] = dd2;
tmatblock.d1[i2-1][n-1] = dd1 * si;
tmatblock.d2[i2-1][n-1] = -dd2 * si;
}
}
for (i = 1; i <= ngss; ++i) 
{
rr[i - 1] = w[i-1] * r[i-1];
}
for (n1 = mm1; n1 <=nmax; ++n1) 
{
an1 = an[n1-1];
for (n2 = mm1; n2 <= nmax; ++n2) 
{
an2 = an[n2-1];
ar12 = 0.0;
ar21 = 0.0;
ai12 = 0.0;
ai21 = 0.0;
gr12 = 0.0;
gr21 = 0.0;
gi12 = 0.0;
gi21 = 0.0;
if (ncheck == 1 && sig[n1 + n2 - 1] < 0.0) 
{
an12 = ann[n1-1][n2-1] * factor;
rirgigblock.r12[n1-1][n2-1] = ar12 * an12;
rirgigblock.r21[n1-1][n2-1] = ar21 * an12;
rirgigblock.i12[n1-1][n2-1] = ai12 * an12;
rirgigblock.i21[n1-1][n2-1] = ai21 * an12;
rirgigblock.rg12[n1-1][n2-1] = gr12 * an12;
rirgigblock.rg21[n1-1][n2-1] = gr21 * an12;
rirgigblock.ig12[n1-1][n2-1] = gi12 * an12;
rirgigblock.ig21[n1-1][n2-1] = gi21 * an12;
}
else
{
for (i = 1; i <= ngss; ++i) 
{
d1n1 = tmatblock.d1[i-1][n1-1];
d2n1 = tmatblock.d2[i-1][n1-1];
d1n2 = tmatblock.d1[i-1][n2-1];
d2n2 = tmatblock.d2[i-1][n2-1];
a12 = d1n1 * d2n2;
a21 = d2n1 * d1n2;
a22 = d2n1 * d2n2;
aa1 = a12 + a21;
qj1 = cbassblock.j[i-1][n1-1];
qy1 = cbassblock.y[i-1][n1-1];
qjr2 = cbassblock.jr[i-1][n2-1];
qji2 = cbassblock.ji[i-1][n2-1];
qdjr2 = cbassblock.djr[i-1][n2-1];
qdji2 = cbassblock.dji[i-1][n2-1];
qdj1 = cbassblock.dj[i-1][n1-1];
qdy1 = cbassblock.dy[i-1][n1-1];
c1r = qjr2 * qj1;
c1i = qji2 * qj1;
b1r = c1r - qji2 * qy1;
b1i = c1i + qjr2 * qy1;
c2r = qjr2 * qdj1;
c2i = qji2 * qdj1;
b2r = c2r - qji2 * qdy1;
b2i = c2i + qjr2 * qdy1;
ddri = ddr[i-1];
c3r = ddri * c1r;
c3i = ddri * c1i;
b3r = ddri * b1r;
b3i = ddri * b1i;
c4r = qdjr2 * qj1;
c4i = qdji2 * qj1;
b4r = c4r - qdji2 * qy1;
b4i = c4i + qdjr2 * qy1;
drri = drr[i-1];
drii = dri[i-1];
c5r = c1r * drri - c1i * drii;
c5i = c1i * drri + c1r * drii;
b5r = b1r * drri - b1i * drii;
b5i = b1i * drri + b1r * drii;
uri = dr[i-1];
rri = rr[i - 1];
f1 = rri * a22;
f2 = rri * uri * an1 * a12;
ar12 = ar12 + f1 * b2r + f2 * b3r;
ai12 = ai12 + f1 * b2i + f2 * b3i;
gr12 = gr12 + f1 * c2r + f2 * c3r;
gi12 = gi12 + f1 * c2i + f2 * c3i;
f2 = rri * uri * an2 * a21;
ar21 = ar21 + f1 * b4r + f2 * b5r;
ai21 = ai21 + f1 * b4i + f2 * b5i;
gr21 = gr21 + f1 * c4r + f2 * c5r;
gi21 = gi21 + f1 * c4i + f2 * c5i;
}
an12 = ann[n1-1][n2-1] * factor;
rirgigblock.r12[n1-1][n2-1] = ar12 * an12;
rirgigblock.r21[n1-1][n2-1] = ar21 * an12;
rirgigblock.i12[n1-1][n2-1] = ai12 * an12;
rirgigblock.i21[n1-1][n2-1] = ai21 * an12;
rirgigblock.rg12[n1-1][n2-1] = gr12 * an12;
rirgigblock.rg21[n1-1][n2-1] = gr21 * an12;
rirgigblock.ig12[n1-1][n2-1] = gi12 * an12;
rirgigblock.ig21[n1-1][n2-1] = gi21 * an12;
}
}
}
tpir = pir;
tpii = pii;
tppi = ppi;
nm = nmax;
for (n1 = mm1; n1 <= nmax; ++n1) 
{
k1 = n1 - mm1+1;
kk1 = k1 + nm;
for (n2 = mm1; n2 <= nmax; ++n2) 
{
k2 = n2 - mm1+1;
kk2 = k2 + nm;
tar12 = rirgigblock.i12[n1-1][n2-1];
tai12 = -rirgigblock.r12[n1-1][n2-1];
tgr12 = rirgigblock.ig12[n1-1][n2-1];
tgi12 = -rirgigblock.rg12[n1-1][n2-1];
tar21 = -rirgigblock.i21[n1-1][n2-1];
tai21 = rirgigblock.r21[n1-1][n2-1];
tgr21 = -rirgigblock.ig21[n1-1][n2-1];
tgi21 = rirgigblock.rg21[n1-1][n2-1];
trtiblock.tqr[k1-1][k2-1] = tpir * tar21 - tpii * tai21 + tppi * tar12;
trtiblock.tqi[k1-1][k2-1] = tpir * tai21 + tpii * tar21 + tppi * tai12;
trgqr[k1-1][k2-1] = tpir * tgr21 - tpii * tgi21 + tppi * tgr12;
trgqi[k1-1][k2-1] = tpir * tgi21 + tpii * tgr21 + tppi * tgi12;
trtiblock.tqr[k1-1][kk2-1]= 0.0;
trtiblock.tqi[k1-1][kk2-1]= 0.0;
trgqr[k1-1][kk2-1]= 0.0;
trgqi[k1-1][kk2-1]= 0.0;
trtiblock.tqr[kk1-1][k2-1] = 0.0;
trtiblock.tqi[kk1-1][k2-1] = 0.0;
trgqr[kk1-1][k2-1] = 0.0;
trgqi[kk1-1][k2-1] = 0.0;
trtiblock.tqr[kk1-1][kk2-1] = tpir * tar12 - tpii * tai12 + tppi * tar21;
trtiblock.tqi[kk1-1][kk2-1] = tpir * tai12 + tpii * tar12 + tppi * tai21;
trgqr[kk1-1][kk2-1] = tpir * tgr12 - tpii * tgi12 + tppi * tgr21;
trgqi[kk1-1][kk2-1] = tpir * tgi12 + tpii * tgr12 + tppi * tgi21;
}
}
nnmax = 2*nm;
for (n1 = 1; n1 <= nnmax; ++n1) 
{
for (n2 = 1; n2 <= nnmax; ++n2) 
{
qr[n1-1][n2-1] = trtiblock.tqr[n1-1][n2-1];
qi[n1-1][n2-1] = trtiblock.tqi[n1-1][n2-1];
rgqr[n1-1][n2-1] = trgqr[n1-1][n2-1];
rgqi[n1-1][n2-1] = trgqi[n1-1][n2-1];
}
}
tt(nmax, ncheck,trtiblock);
return 0;
}

public int tmatr(int m, int ngauss, double x[],double w[], double an[], double ann[][], double s[], double ss[], double ppi,
 double pir, double pii,double r[], double dr[], double ddr[], double drr[], double dri[], int nmax, int ncheck,tmatBlock tmatblock,trtiBlock trtiblock,rirgigBlock rirgigblock,cbassBlock cbassblock)
{
int i;
int  n;
int i1;
int  i2;
int k1;
int  n1;
int  n2;
int  k2;
int ng;
int nm;
int mm1;
int kk1;
int  kk2;
int ngss;
int nnmax;
double e1;
double  f1;
double  f2;
double e2;
double  e3;
double a11;
double  a12;
double  a21;
double  a22;
double  ds[]=new double[npng2];
double qm;
double  si;
double  rr[]=new double[npng2];
double  wr;
double  aa1;
double  aa2;
double  dd1;
double  dd2;
double  b1i;
double  c1i;
double  c2i;
double  b2i;
double  an1;
double  an2;
double  c3i;
double  b3i;
double  c4i;
double  b1r;
double  c1r;
double  c2r;
double  b2r;
double  c3r;
double  b3r;
double  dv1[]=new double[npn1];
double  dv2[]=new double[npn1];
double  qj1;
double  c4r;
double  b4r;
double  b4i;
double  c5r;
double c5i;
double  b5r;
double  b5i;
double  c6r;
double  c6i;
double  b6r;
double  qy1;
double  b6i;
double  c7r;
double  c7i;
double  b7r;
double  b7i;
double  c8r;
double  c8i;
double  b8r;
double  b8i;
double ai11;
double  ai12;
double  ai21;
double  ai22;
double  ar11;
double  ar12;
double  ar21;
double sig[]=new double[npn2];
double ar22;
double  gr11;
double gr12;
double  dss[]=new double[npng2];
double  qmm;
double  gr21;
double  gr22;
double  gi11;
double  gi12;
double  gi21;
double  gi22;
double  d1n1;
double  d2n1;
double  d1n2;
double  d2n2;
double  uri;
double  dsi;
double  rri;
double  an12;
double  qdj1;
double  qji2;
double  qjr2;
double  qdy1;
double  tai11;
double  tai12;
double  tai21;
double  ddri;
double  tai22;
double  tgi11;
double  tgi12;
double  drii;
double  tar11;
double  tar12;
double  tar21;
double  tgi21;
double  tar22;
double  tgi22;
double  tgr11;
double  tgr12;
double  drri;
double  tgr21;
double  dssi;
double  tpii;
double  tgr22;
double tppi;
double  tpir;
double  qdji2;
double  qdjr2;
double factor;

mm1 = m;
qm = (double) (m);
qmm = qm * qm;
ng = 2*(ngauss);
ngss = ng;
factor = 1.;
if (ncheck == 1) 
{
ngss = ngauss;
factor = 2.;
   } 
else 
{
//CONTINUE
}
si = 1.;
nm = nmax + nmax;
for (n = 1; n <= nm; ++n) 
{
si = -si;
sig[n - 1] = si;
}
for (i = 1; i <= ngauss; ++i) 
{
i1 = ngauss + i;
i2 = ngauss - i + 1;
vig(x[i1-1], nmax, m, dv1, dv2);
for (n = 1; n <= nmax; ++n) 
{
si = sig[n - 1];
dd1 = dv1[n - 1];
dd2 = dv2[n - 1];
tmatblock.d1[i1-1][n-1] = dd1;
tmatblock.d2[i1-1][n-1] = dd2;
tmatblock.d1[i2-1][n-1] = dd1 * si;
tmatblock.d2[i2-1][n-1] = -dd2 * si;
}
}
for (i = 1; i <= ngss; ++i) 
{
wr = w[i-1] * r[i-1];
ds[i - 1] = s[i-1] * qm * wr;
dss[i - 1] = ss[i-1] * qmm;
rr[i - 1] = wr;
}
for (n1 = mm1; n1 <= nmax; ++n1) 
{
an1 = an[n1-1];
for (n2 = mm1; n2 <= nmax; ++n2) 
{
an2 = an[n2-1];
ar11 = 0.0;
ar12 = 0.0;
ar21 = 0.0;
ar22 = 0.0;
ai11 = 0.0;
ai12 = 0.0;
ai21 = 0.0;
ai22 = 0.0;
gr11 = 0.0;
gr12 = 0.0;
gr21 = 0.0;
gr22 = 0.0;
gi11 = 0.0;
gi12 = 0.0;
gi21 = 0.0;
gi22 = 0.0;
si = sig[n1 + n2 - 1];
for (i = 1; i<= ngss; ++i) 
{
d1n1 = tmatblock.d1[i-1][n1-1];
d2n1 = tmatblock.d2[i-1][n1-1];
d1n2 = tmatblock.d1[i-1][n2-1];
d2n2 = tmatblock.d2[i-1][n2-1];
a11 = d1n1 * d1n2;
a12 = d1n1 * d2n2;
a21 = d2n1 * d1n2;
a22 = d2n1 * d2n2;
aa1 = a12 + a21;
aa2 = a11 * dss[i - 1] + a22;
qj1 = cbassblock.j[i-1][n1-1];
qy1 =cbassblock.y[i-1][n1-1];
qjr2 = cbassblock.jr[i-1][n2-1];
qji2 =cbassblock.ji[i-1][n2-1];
qdjr2 = cbassblock.djr[i-1][n2-1];
qdji2 = cbassblock.dji[i-1][n2-1];
qdj1 = cbassblock.dj[i-1][n1-1];
qdy1 = cbassblock.dy[i-1][n1-1];
c1r = qjr2 * qj1;
c1i = qji2 * qj1;
b1r = c1r - qji2 * qy1;
b1i = c1i + qjr2 * qy1;
c2r = qjr2 * qdj1;
c2i = qji2 * qdj1;
b2r = c2r - qji2 * qdy1;
b2i = c2i + qjr2 * qdy1;
ddri = ddr[i-1];
c3r = ddri * c1r;
c3i = ddri * c1i;
b3r = ddri * b1r;
b3i = ddri * b1i;
c4r = qdjr2 * qj1;
c4i = qdji2 * qj1;
b4r = c4r - qdji2 * qy1;
b4i = c4i + qdjr2 * qy1;
drri = drr[i-1];
drii = dri[i-1];
c5r = c1r * drri - c1i * drii;
c5i = c1i * drri + c1r * drii;
b5r = b1r * drri - b1i * drii;
b5i = b1i * drri + b1r * drii;
c6r = qdjr2 * qdj1;
c6i = qdji2 * qdj1;
b6r = c6r - qdji2 * qdy1;
b6i = c6i + qdjr2 * qdy1;
c7r = c4r * ddri;
c7i = c4i * ddri;
b7r = b4r * ddri;
b7i = b4i * ddri;
c8r = c2r * drri - c2i * drii;
c8i = c2i * drri + c2r * drii;
b8r = b2r * drri - b2i * drii;
b8i = b2i * drri + b2r * drii;
uri = dr[i-1];
dsi = ds[i - 1];
dssi = dss[i - 1];
rri = rr[i - 1];
if (ncheck == 1) 
{
if(si > 0.0)
{
f1 = rri * aa2;
f2 = rri * uri * an1 * a12;
ar12 = ar12 + f1 * b2r + f2 * b3r;
ai12 = ai12 + f1 * b2i + f2 * b3i;
gr12 = gr12 + f1 * c2r + f2 * c3r;
gi12 = gi12 + f1 * c2i + f2 * c3i;
f2 = rri * uri * an2 * a21;
ar21 = ar21 + f1 * b4r + f2 * b5r;
ai21 = ai21 + f1 * b4i + f2 * b5i;
gr21 = gr21 + f1 * c4r + f2 * c5r;
gi21 = gi21 + f1 * c4i + f2 * c5i;
}
else
{
e1 = dsi * aa1;
ar11 = ar11+e1 * b1r;
ai11 = ai11+e1 * b1i;
gr11 = gr11+e1 * c1r;
gi11 = gi11+e1 * c1i;
e2 = dsi * uri * a11;
e3 = e2 * an2;
e2 = e2 * an1;
ar22 = ar22 + e1 * b6r + e2 * b7r + e3 * b8r;
ai22 = ai22 + e1 * b6i + e2 * b7i + e3 * b8i;
gr22 = gr22 + e1 * c6r + e2 * c7r + e3 * c8r;
gi22 = gi22 + e1 * c6i + e2 * c7i + e3 * c8i;
}
}
else
{
e1 = dsi * aa1;
ar11 = ar11+e1 * b1r;
ai11 = ai11+e1 * b1i;
gr11 = gr11+e1 * c1r;
gi11 = gi11+e1 * c1i;
f1 = rri * aa2;
f2 = rri * uri * an1 * a12;
ar12 = ar12 + f1 * b2r + f2 * b3r;
ai12 = ai12 + f1 * b2i + f2 * b3i;
gr12 = gr12 + f1 * c2r + f2 * c3r;
gi12 = gi12 + f1 * c2i + f2 * c3i;
f2 = rri * uri * an2 * a21;
ar21 = ar21 + f1 * b4r + f2 * b5r;
ai21 = ai21 + f1 * b4i + f2 * b5i;
gr21 = gr21 + f1 * c4r + f2 * c5r;
gi21 = gi21 + f1 * c4i + f2 * c5i;
e2 = dsi * uri * a11;
e3 = e2 * an2;
e2 = e2 * an1;
ar22 = ar22 + e1 * b6r + e2 * b7r + e3 * b8r;
ai22 = ai22 + e1 * b6i + e2 * b7i + e3 * b8i;
gr22 = gr22 + e1 * c6r + e2 * c7r + e3 * c8r;
gi22 = gi22 + e1 * c6i + e2 * c7i + e3 * c8i;
}
}
an12 = ann[n1-1][n2-1] * factor;
rirgigblock.r11[n1-1][n2-1] = ar11 * an12;
rirgigblock.r12[n1-1][n2-1] = ar12 * an12;
rirgigblock.r21[n1-1][n2-1] = ar21 * an12;
rirgigblock.r22[n1-1][n2-1] = ar22 * an12;
rirgigblock.i11[n1-1][n2-1] = ai11 * an12;
rirgigblock.i12[n1-1][n2-1] = ai12 * an12;
rirgigblock.i21[n1-1][n2-1] = ai21 * an12;
rirgigblock.i22[n1-1][n2-1] = ai22 * an12;
rirgigblock.rg11[n1-1][n2-1] = gr11 * an12;
rirgigblock.rg12[n1-1][n2-1] = gr12 * an12;
rirgigblock.rg21[n1-1][n2-1] = gr21 * an12;
rirgigblock.rg22[n1-1][n2-1] = gr22 * an12;
rirgigblock.ig11[n1-1][n2-1] = gi11 * an12;
rirgigblock.ig12[n1-1][n2-1] = gi12 * an12;
rirgigblock.ig21[n1-1][n2-1] = gi21 * an12;
rirgigblock.ig22[n1-1][n2-1] = gi22 * an12;
}
}
tpir = pir;
tpii = pii;
tppi = ppi;
nm = nmax - mm1 + 1;
for (n1 = mm1; n1 <=nmax; ++n1) 
{
k1 = n1 - mm1 + 1;
kk1 = k1 + nm;
for (n2 = mm1; n2 <= nmax; ++n2) 
{
k2 = n2 - mm1 + 1;
kk2 = k2 + nm;
tar11 = -rirgigblock.r11[n1-1][n2-1];
tai11 = -rirgigblock.i11[n1-1][n2-1];
tgr11 = -rirgigblock.rg11[n1-1][n2-1];
tgi11 = -rirgigblock.ig11[n1-1][n2-1];
tar12 = rirgigblock.i12[n1-1][n2-1];
tai12 = -rirgigblock.r12[n1-1][n2-1];
tgr12 = rirgigblock.ig12[n1-1][n2-1];
tgi12 = -rirgigblock.rg12[n1-1][n2-1];
tar21 = -rirgigblock.i21[n1-1][n2-1];
tai21 = rirgigblock.r21[n1-1][n2-1];
tgr21 = -rirgigblock.ig21[n1-1][n2-1];
tgi21 = rirgigblock.rg21[n1-1][n2-1];
tar22 = -rirgigblock.r22[n1-1][n2-1];
tai22 = -rirgigblock.i22[n1-1][n2-1];
tgr22 = -rirgigblock.rg22[n1-1][n2-1];
tgi22 = -rirgigblock.ig22[n1-1][n2-1];
trtiblock.tqr[k1-1][k2-1] = tpir * tar21 - tpii * tai21 + tppi * tar12;
trtiblock.tqi[k1-1][k2-1] = tpir * tai21 + tpii * tar21 + tppi * tai12;
trgqr[k1-1][k2-1] = tpir * tgr21 - tpii * tgi21 + tppi * tgr12;
trgqi[k1-1][k2-1] = tpir * tgi21 + tpii * tgr21 + tppi * tgi12;
trtiblock.tqr[k1-1][kk2-1] = tpir * tar11 - tpii * tai11 + tppi * tar22;
trtiblock.tqi[k1-1][kk2-1] = tpir * tai11 + tpii * tar11 + tppi * tai22;
trgqr[k1-1][kk2-1] = tpir * tgr11 - tpii * tgi11 + tppi * tgr22;
trgqi[k1-1][kk2-1] = tpir * tgi11 + tpii * tgr11 + tppi * tgi22;
trtiblock.tqr[kk1-1][k2-1] = tpir * tar22 - tpii * tai22 + tppi * tar11;
trtiblock.tqi[kk1-1][k2-1] = tpir * tai22 + tpii * tar22 + tppi * tai11;
trgqr[kk1-1][k2-1] = tpir * tgr22 - tpii * tgi22 + tppi * tgr11;
trgqi[kk1-1][k2-1] = tpir * tgi22 + tpii * tgr22 + tppi * tgi11;
trtiblock.tqr[kk1-1][kk2-1] = tpir * tar12 - tpii * tai12 + tppi * tar21;
trtiblock.tqi[kk1-1][kk2-1] = tpir * tai12 + tpii * tar12 + tppi * tai21;
trgqr[kk1-1][kk2-1] = tpir * tgr12 - tpii * tgi12 + tppi * tgr21;
trgqi[kk1-1][kk2-1] = tpir * tgi12 + tpii * tgr12 + tppi * tgi21;
}
}
nnmax = 2*nm;
for (n1 = 1; n1 <= nnmax; ++n1) 
{
for (n2 = 1; n2 <= nnmax; ++n2) 
{
qr[n1-1][n2-1] = trtiblock.tqr[n1-1][n2-1];
qi[n1-1][n2-1] = trtiblock.tqi[n1-1][n2-1];
rgqr[n1-1][n2-1] = trgqr[n1-1][n2-1];
rgqi[n1-1][n2-1] = trgqi[n1-1][n2-1];
}
}
tt(nm, ncheck,trtiblock);
return 0;
}

public int vig(double x, int nmax, int m,double dv1[], double dv2[])
{
int i;
int n;
int i2;
double a;
double d1;
double d2;
double d3;
double qn;
double qs;
double qn1;
double qn2;
double qs1;
double der;
double qmm;
double qnm;
double qnm1;

a = 1.0;
qs = Math.sqrt(1.0 - x * x);
qs1 = 1.0 / qs;
for (n = 1; n <= nmax; ++n) 
{
dv1[n-1] = 0.0;
dv2[n-1] = 0.0;
}
if (m == 0) 
{
d1 = 1.0;
d2 = x;
for (n = 1; n <= nmax; ++n) 
{
qn = (double) n;
qn1 = (double) (n + 1);
qn2 = (double) (2*n + 1);
d3 = (qn2 * x * d2 - qn * d1) / qn1;
der = qs1 * (qn1 * qn / qn2) * (-d1 + d3);
dv1[n-1] = d2;
dv2[n-1] = der;
d1 = d2;
d2 = d3;
}
return 0;
}
qmm = (double) (m * m);
for (i = 1; i <= m; ++i) 
{
i2 = 2*i;
a = a * Math.sqrt((double) (i2 - 1) / (double) i2) * qs;
}
d1 = 0.0;
d2 = a;
for (n = m; n <= nmax; ++n) 
{
qn = (double) n;
qn2 = (double) (2*n + 1);
qn1 = (double) (n + 1);
qnm = Math.sqrt(qn * qn - qmm);
qnm1 = Math.sqrt(qn1 * qn1 - qmm);
d3 = (qn2 * x * d2 - qnm * d1) / qnm1;
der = qs1 * (-qn1 * qnm * d1 + qn * qnm1 * d3) / qn2;
dv1[n-1] = d2;
dv2[n-1] = der;
d1 = d2;
d2 = d3;
}
return 0;
} 

public int tt(int nmax, int ncheck,trtiBlock trtiblock)
{
int i;
int l;
int n1;
int n2;
int ndim;
int ipvt[]=new int[npn2];
int nnmax;
double b[]=new double[npn2];
double cond;
double work[]=new double[npn2];
double a[][]=new double[npn2][npn2];
double f[][]=new double[npn2][npn2];
double c[][]=new double[npn2][npn2];
double d[][]=new double[npn2][npn2];
double e[][]=new double[npn2][npn2];

ndim = 200;
nnmax = 2*(nmax);
for (n1 = 1; n1 <= nnmax; ++n1) 
{
for (n2 = 1; n2 <= nnmax; ++n2) 
{
f[n1-1][n2-1] = qi[n1-1][n2-1];
}
}
if (ncheck == 1) 
{
inv1(nmax, f, a);
} 
else 
{
invert(ndim, nnmax, f, a, ipvt, work, b);
}
prod(qr, a, c, ndim, nnmax);
prod(c, qr, d, ndim, nnmax);
for (n1 = 1; n1 <= nnmax; ++n1) 
{
for (n2 = 1; n2 <= nnmax; ++n2) 
{
c[n1-1][n2-1] = d[n1-1][n2-1] + qi[n1-1][n2-1];
}
}
if (ncheck == 1) 
{
inv1(nmax, c, qi);
} 
else 
{
invert(ndim, nnmax, c, qi, ipvt, work, b);
}
prod(a, qr, d, ndim, nnmax);
prod(d, qi, qr, ndim, nnmax);
prod(rgqr, qr, a, ndim, nnmax);
prod(rgqi, qi, c, ndim, nnmax);
prod(rgqr, qi, d, ndim, nnmax);
prod(rgqi, qr, e, ndim, nnmax);
for (n1 = 1; n1 <= nnmax; ++n1) 
{
for (n2 = 1; n2 <= nnmax; ++n2) 
{
trtiblock.tr1[n1-1][n2-1] = -a[n1-1][n2-1] - c[n1-1][n2-1];
trtiblock.ti1[n1-1][n2-1] = d[n1-1][n2-1] - e[n1-1][n2-1];
}
}
return 0;
}

public int prod(double ax[][], double bx[][], double c[][],int ndim, int n)
{
int i;
int j;
int k;
int l;
double cij;

for (i = 1; i <= n; ++i) 
{
for (j = 1; j <= n; ++j) 
{
cij = 0.;
for (k = 1; k <= n; ++k) 
{
cij = cij + ax[i-1][k-1] * bx[k-1][j-1];
}
c[i-1][j-1] = cij;
}
}
return 0;
} 

public int inv1(int nmax, double f[][], double a[][])
{
int i;
int  l;
int j;
int  i1;
int  i2;
int  j1;
int  j2;
int  nn1;
int  nn2;
int  ind1[]= new int[npn1];
int  ind2[]= new int[npn1];
int  ndim;
int  ipvt[]= new int[npn1];
int  nnmax;
double b[]= new double[npn1];
double cond;
double work[]=new double[npn1];
double p1[][]= new double[npn1][npn1];
double q1[][]= new double[npn1][npn1];
double p2[][]= new double[npn1][npn1];
double q2[][]= new double[npn1][npn1];

ndim = npn1;
nn1 = (int) (((double) (nmax) - 0.1) * 0.5 + 1.0);
nn2 = nmax - nn1;
for (i = 1; i <= nmax; ++i) 
{
ind1[i - 1] = 2*i - 1;
if (i > nn1) 
{
ind1[i - 1] = nmax + 2*(i - nn1);
}
ind2[i - 1] = 2*i;
if (i > nn2) 
{
ind2[i - 1] = nmax + 2*(i - nn2)-1;
}
}
nnmax = 2*(nmax);
for (i = 1; i <= nmax; ++i) 
{
i1 = ind1[i - 1];
i2 = ind2[i - 1];
for (j = 1; j <= nmax; ++j) 
{
j1 = ind1[j - 1];
j2 = ind2[j - 1];
q1[j-1][i-1] = f[j1-1][i1-1];
q2[j-1][i-1] = f[j2-1][i2-1];
}
}
invert(ndim, nmax, q1, p1,  ipvt, work, b);
invert(ndim, nmax, q2, p2,  ipvt, work, b);
for (i = 1; i <= nnmax; ++i) 
{
for (j = 1; j <= nnmax; ++j) 
{
a[j-1][i-1] = 0.0;
}
}
for (i = 1; i <= nmax; ++i) 
{
i1 = ind1[i - 1];
i2 = ind2[i - 1];
for (j = 1; j <= nmax; ++j) 
{
j1 = ind1[j - 1];
j2 = ind2[j - 1];
a[j1-1][i1-1] = p1[j-1][i-1];
a[j2-1][i2-1] = p2[j-1][i-1];
}
}
return 0;
}

public int invert(int ndim, int n, double a[][],double xa[][], int ipvt[], double work[], double b[])
{
int i, j, k, l,m,jo, kb, km1, nm1, kp1;
double t, ek, anorm, ynorm, znorm;

ipvt[n-1] = 1;
if (n != 1) 
{
nm1 = n - 1;
anorm = 0.0;
for (j = 1; j <= n; ++j) 
{
t = 0.0;
for (i = 1; i <= n; ++i) 
{
t = t+ Math.abs(a[i-1][j-1]);
}
if (t > anorm) 
{
anorm = t;
}
}
for (k = 1; k <= nm1; ++k)
{
kp1 = k + 1;
m = k;
jo=k;
for (i = kp1; i <= n; ++i) 
{
if (Math.abs(a[i-1][k-1]) > Math.abs(a[m-1][k-1]))
{
m = i;
}
if(Math.abs(a[i-1][k-1]) < Math.abs(a[m-1][k-1]))
{
m=m;
}
}
ipvt[k-1] = m;
if (m != k) 
{
ipvt[n-1] = -ipvt[n-1];
}
t = a[m-1][k-1];
a[m-1][k-1] = a[k-1][k-1];
a[k-1][k-1] = t;
if (t != 0.) 
{
for (i = kp1; i <= n; ++i)
{
a[i-1][k-1] = -a[i-1][k-1]/ t;
}
for (j = kp1; j <= n; ++j) 
{
t = a[m-1][j-1];
a[m-1][j-1] = a[k-1][j-1];
a[k-1][j-1] = t;
if (t != 0.0) 
{
for (i= kp1; i <=n; ++i) 
{
a[i-1][j-1] = a[i-1][j-1]+ a[i-1][k-1] * t;
}
}
}
}
}
for (k = 1; k <= n; ++k) 
{
t = 0.0;
if (k != 1) 
{
km1 = k - 1;
for (i = 1; i <= km1; ++i) 
{
t = t + a[i-1][k-1] * work[i-1];
}
}
ek = 1.0;
if (t < 0.0) 
{
ek = -1.0;
}
if (a[k-1][k-1]== 0.0) 
{
cond = 1e+52;

if (cond + 1.0 == cond) 
{
System.out.println("THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY COND = "+cond);
}
for (i= 1; i<= n; ++i) 
{
for (j = 1; j <= n; ++j) 
{
b[j-1] = 0.0;
if (j == i) 
{
b[j-1] = 1.0;
}
}
solve(ndim, n, a, b, ipvt);
for (j = 1; j <= n; ++j) 
{
xa[j-1][i-1] = b[j-1];
}
}
return 0;
}
work[k-1] = -(ek + t) / a[k-1][k-1];
}
for (kb = 1; kb <= nm1; ++kb) 
{
k = n - kb;
t = 0.0;
kp1 = k + 1;
for (i = kp1; i <= n; ++i) 
{
t = t + a[i-1][k-1]* work[k-1];
}
work[k-1] = t;
m = ipvt[k-1];
if (m != k) 
{
t = work[m-1];
work[m-1] = work[k-1];
work[k-1] = t;
}
}
ynorm = 0.0;
for (i = 1; i <= n; ++i) 
{
ynorm = ynorm+(Math.abs(work[i-1]));
}
solve(ndim, n, a, work, ipvt);
znorm = 0.0;
for (i = 1; i <= n; ++i) 
{
znorm = znorm+(Math.abs(work[i-1]));
}
cond = anorm * znorm / ynorm;
if (cond < 1.0) 
{
cond = 1.0;
}
}
cond = 1.0;
if (a[0][0] == 0.0) 
{
cond = 1e+52;
}
if (cond + 1.0 == cond) 
{
System.out.println("THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY COND = "+cond);
}
for (i= 1; i<= n; ++i) 
{
for (j = 1; j <= n; ++j) 
{
b[j-1] = 0.0;
if (j == i) 
{
b[j-1] = 1.0;
}
}
solve(ndim, n, a, b, ipvt);
for (j = 1; j <= n; ++j) 
{
xa[j-1][i-1] = b[j-1];
/* L30: */
}
}
return 0;
}

public int solve(int ndim, int n, double a[][],double b[], int ipvt[])
{
int i;
int  k;
int  m;
int  kb;
int  km1;
int  nm1;
int  kp1;
double t;

if (n != 1) 
{
nm1 = n - 1;
for (k = 1; k <= nm1; ++k) 
{
kp1 = k + 1;
m = ipvt[k-1];
t = b[m-1];
b[m-1] = b[k-1];
b[k-1] = t;
for (i = kp1; i <= n; ++i) 
{
b[i-1] = b[i-1] + a[i-1][k-1] * t;
}
}
for (kb = 1; kb <= nm1; ++kb) 
{
km1 = n - kb;
k = km1 + 1;
b[k-1] = b[k-1] / a[k-1][k-1];
t = -b[k-1];
for (i = 1; i <= km1; ++i) 
{
b[i-1] = b[i-1] + a[i-1][k-1] * t;
}
}
}
b[0] = b[0] / a[0][0];
return 0;
}

public int gsp(int nmax, double csca, double lam, double alf1[], double alf2[], double alf3[], double alf4[], double bet1[],
double bet2[],gspBlock gspblock)//, int *lmax
{
int i;
int  j;
int  l;
int  m;
int  n;
int i1;
int  k1;
int  k2;
int  k3;
int  k4;
int  k5;
int  k6;
int  m1;
int  l1;
int  n1;
int  m2;
int  kn;
int  nl;
int  nn;
int  nn1;
int  nnn;
int  mmin;
int  mmax;
int  m1min;
int  l1max;
int  nmax1;
int  m1max;
int nnmin;
int  nnmax;
int  nn1min;
int  nn1max;
double  t1;
double t2;
double t3;
double t4;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double ff[][]=new double[npn4][npn4];
double si;
double sj;
double sl;
double xi;
double xr;
double xx;
double dd1;
double tr1[][]=new double[npl1][npn4];
double tr2[][]=new double[npl1][npn4];
double dd2;
double ai1[]=new double[npn4];
double ai2[]=new double[npn4];
double dd3;
double dd4;
double dm1;
double ar1[]=new double[npn4];
double ar2[]=new double[npn4];
double g1l;
double g2l;
double g3l;
double ti1[][]=new double[npl1][npn4];
double ti2[][]=new double[npl1][npn4];
double ri2;
double g4l;
double dm2;
double dm3;
double dm4;
double rr1;
double rr2;
double tt1;
double tt2;
double ffn;
double sig;
double tt3;
double tt4;
double tt5;
double tt6;
double tt7;
double tt8;
double ssi[]=new double[npl];
double ssj[]=new double[npn1];
double sss;
double aai1;
double aai2;
double bbi1;
double bbi2;
double aar1;
double aar2;
double bbr1;
double bbr2;
double dd5i;
double dd5r;
double dm5i;
double g5li;
double dm5r;
double g5lr;
double sss1=0.0;
double fi[][]=new double[npn4][npn4];
double dk;
double fr[][]=new double[npn4][npn4];
double ri1;
double d1[][][]=new double[npl1][npn4][npn4];
double d2[][][]=new double[npl1][npn4][npn4];
double d3[][][]=new double[npl1][npn4][npn4];
double d4[][][]=new double[npl1][npn4][npn4];
double d5r[][][]=new double[npl1][npn4][npn4]; 
double d5i[][][]=new double[npl1][npn4][npn4];
double g1[][]=new double[npl1][npn6];
double g2[][]=new double[npl1][npn6];
double m1r[][][]=new double[npl1][npl1][npn4];
double m1i[][][]=new double[npl1][npl1][npn4];
double m2r[][][]=new double[npl1][npl1][npn4];
double m2i[][][]=new double[npl1][npl1][npn4];
Complex ci;
Complex cci;
Complex ccj;
Complex cim[]=new Complex[100];
Complex c1;

signum();
lmax = 2*(nmax);
l1max = lmax + 1;
ci = new Complex(0.0,1.0);
cim[0]=new Complex(ci.real(),ci.imag());
for (i = 2; i <= nmax; ++i) 
{
cim[i-1] = cim[i-2].times(ci);
}
ssi[0] = 1.0;
for (i = 1; i<= lmax; ++i) 
{
i1 = i + 1;
si = (double) (2*i + 1);
ssi[i1 - 1] = si;
if (i <= nmax) 
{
ssj[i - 1] = Math.pow(si,0.5);
}
}
ci=new Complex(-ci.real(),-ci.imag());
for (i = 1; i<= nmax; ++i) 
{
si = ssj[i - 1];
cci = new Complex(cim[i-1].real(),cim[i-1].imag());
for (j = 1; j <= nmax; ++j) 
{
sj = 1.0 / ssj[j - 1];
c1=new Complex(sj*cim[j-1].real(),sj*cim[j-1].imag());
ccj=c1.div(cci);
fr[j-1][i-1]=ccj.real();
fi[j-1][i-1]=(ccj.times(ci)).real();
ff[j-1][i-1] = si * sj;
}
}
nmax1 = nmax + 1;
k1 = 1;
k2 = 0;
k3 = 0;
k4 = 1;
k5 = 1;
k6 = 2;
for (n = 1; n <= nmax; ++n) 
{
for (nn = 1; nn <= nmax; ++nn) 
{
m1max = Math.min(n,nn) + 1;
for (m1 = 1; m1 <= m1max; ++m1) 
{
m = m1 - 1;
l1 = npn6+m;
tt1 = gspblock.rt11[m1-1][n-1][nn-1];
tt2 = gspblock.rt12[m1-1][n-1][nn-1];
tt3 = gspblock.rt21[m1-1][n-1][nn-1];
tt4 =gspblock.rt22[m1-1][n-1][nn-1];
tt5 = gspblock.it11[m1-1][n-1][nn-1];
tt6 = gspblock.it12[m1-1][n-1][nn-1];
tt7 = gspblock.it21[m1-1][n-1][nn-1];
tt8 = gspblock.it22[m1-1][n-1][nn-1];
t1 = tt1 + tt2;
t2 = tt3 + tt4;
t3 = tt5 + tt6;
t4 = tt7 + tt8;
tr1[l1-1][nn-1] = t1 + t2;
tr2[l1-1][nn-1] = t1 - t2;
ti1[l1-1][nn-1] = t3 + t4;
ti2[l1-1][nn-1] = t3 - t4;
if (m != 0) 
{
l1 = 81 - m;
t1 = tt1 - tt2;
t2 = tt3 - tt4;
t3 = tt5 - tt6;
t4 = tt7 - tt8;
tr1[l1-1][nn-1] = t1 - t2;
tr2[l1-1][nn-1] = t1 + t2;
ti1[l1-1][nn-1] = t3 - t4;
ti2[l1-1][nn-1] = t3 + t4;
}
}
}
nn1max = nmax1 + n;
for (nn1 = 1; nn1 <= nn1max; ++nn1) 
{
n1 = nn1 - 1;
ccg(n, n1, nmax, k1, k2, g1);
nnmax = Math.min(nmax,n1 + n);
nnmin = Math.max(1,Math.abs(n - n1));
kn = n + nn1;
for (nn = nnmin; nn <= nnmax; ++nn) 
{
nnn = nn + 1;
sig = ssign[kn + nn - 1];
m1max = npn6+Math.min(n,nn);
aar1 = 0.0;
aar2 = 0.0;
aai1 = 0.0;
aai2 = 0.0;
for (m1 = npn6; m1 <= m1max; ++m1) 
{
m = -(npn6-m1);
sss = g1[m1-1][nnn-1];
rr1 = tr1[m1-1][nn-1];
ri1 = ti1[m1-1][nn-1];
rr2 = tr2[m1-1][nn-1];
ri2 = ti2[m1-1][nn-1];
if (m != 0) 
{
m2 = npn6 - m;
rr1 = rr1 + tr1[m2-1][nn-1] * sig;
ri1 = ri1 + ti1[m2-1][nn-1] * sig;
rr2 = rr2 + tr2[m2-1][nn-1] * sig;
ri2 = ri2 + ti2[m2-1][nn-1] * sig;
}
aar1 += sss * rr1;
aai1 += sss * ri1;
aar2 += sss * rr2;
aai2 += sss * ri2;
}
xr = fr[nn-1][n-1];
xi = fi[nn-1][n-1];
ar1[nn - 1] = aar1 * xr - aai1 * xi;
ai1[nn - 1] = aar1 * xi + aai1 * xr;
ar2[nn - 1] = aar2 * xr - aai2 * xi;
ai2[nn - 1] = aar2 * xi + aai2 * xr;
}
ccg(n, n1, nmax, k3, k4, g2);
m1 = Math.max(-n1 + 1,-n);
m2 = Math.min(n1 + 1,n);
m1max = npn6+m2;
m1min = npn6+m1;
for (m1 = m1min; m1 <= m1max; ++m1) //m1min
{
bbr1 = 0.0;
bbi1 = 0.0;
bbr2 = 0.0;
bbi2 = 0.0;
for (nn = nnmin; nn <= nnmax; ++nn) 
{
nnn = nn + 1;
sss = g2[m1-1][nnn-1];
bbr1 = bbr1+sss * ar1[nn - 1];
bbi1 = bbi1+sss * ai1[nn - 1];
bbr2 = bbr2+sss * ar2[nn - 1];
bbi2 = bbi2+sss * ai2[nn - 1];
}
m1r[nn1-1][m1-1][n-1] = bbr1;
m1i[nn1-1][m1-1][n-1] = bbi1;
m2r[nn1-1][m1-1][n-1] = bbr2;
m2i[nn1-1][m1-1][n-1] = bbi2;
}
}
}
for (n = 1; n <= nmax; ++n) 
{
for (nn = 1; nn <= nmax; ++nn) 
{
m1 = Math.min(n,nn);
m1max = npn6+m1;
m1min = npn6 - m1;
nn1max = nmax1 + Math.min(n,nn);
for (m1 = m1min; m1 <= m1max; ++m1) 
{
m = -(npn6-m1);
nn1min = (Math.abs(m - 1)) + 1;
dd1 = 0.0;
dd2 = 0.0;
for (nn1 = nn1min; nn1 <= nn1max; ++nn1) 
{
xx = ssi[nn1 - 1];
x1 = m1r[nn1-1][m1-1][n-1];
x2 = m1i[nn1-1][m1-1][n-1];
x3 = m1r[nn1-1][m1-1][nn-1];
x4 = m1i[nn1-1][m1-1][nn-1];
x5 = m2r[nn1-1][m1-1][n-1];
x6 = m2i[nn1-1][m1-1][n-1];
x7 = m2r[nn1-1][m1-1][nn-1];
x8 = m2i[nn1-1][m1-1][nn-1];
dd1 = dd1+xx * (x1 * x3 + x2 * x4);
dd2 = dd2+xx * (x5 * x7 + x6 * x8);
}
d1[m1-1][nn-1][n-1] = dd1;
d2[m1-1][nn-1][n-1] = dd2;
}
mmax = Math.min(n,nn + 2);
mmin = Math.max(-n,-nn + 2);
m1max = npn6+mmax;
m1min = npn6+mmin;
for (m1 = m1min; m1 <= m1max; ++m1) 
{
m = -(npn6-m1);
nn1min = (Math.abs(m - 1)) + 1;
dd3 = 0.0;
dd4 = 0.0;
dd5r = 0.0;
dd5i = 0.0;
m2 = npn6 - m + 2;
for (nn1 = nn1min; nn1 <= nn1max; ++nn1) 
{
xx = ssi[nn1 - 1];
x1 = m1r[nn1-1][m1-1][n-1];
x2 = m1i[nn1-1][m1-1][n-1];
x3 = m2r[nn1-1][m1-1][n-1];
x4 = m2i[nn1-1][m1-1][n-1];
x5 = m1r[nn1-1][m2-1][nn-1];
x6 = m1i[nn1-1][m2-1][nn-1];
x7 = m2r[nn1-1][m2-1][nn-1];
x8 = m2i[nn1-1][m2-1][nn-1];
dd3 = dd3+xx * (x1 * x5 + x2 * x6);
dd4 = dd4+xx * (x3 * x7 + x4 * x8);
dd5r = dd5r+xx * (x3 * x5 + x4 * x6);
dd5i = dd5i+xx * (x4 * x5 - x3 * x6);
}
d3[m1-1][nn-1][n-1] = dd3;
d4[m1-1][nn-1][n-1] = dd4;
d5r[m1-1][nn-1][n-1] = dd5r;
d5i[m1-1][nn-1][n-1] = dd5i;
}
}
}
dk = lam * lam / (csca * 4.0 * Math.acos(-1.));
for (l1 = 1; l1 <= l1max; ++l1) //1
{
g1l = 0.0;
g2l = 0.0;
g3l = 0.0;
g4l = 0.0;
g5lr = 0.0;
g5li = 0.0;
l = l1 - 1;
sl = ssi[l1 - 1] * dk;
for (n = 1; n <= nmax; ++n) 
{
nnmin = Math.max(1,Math.abs(n-l));
nnmax = Math.min(nmax,n + l);
if (nnmax >= nnmin) 
{
ccg(n, l, nmax, k1, k2, g1);
if (l >= 2) 
{
ccg(n, l, nmax, k5, k6, g2);
}
nl = n + l;
for (nn = nnmin; nn <= nnmax; ++nn) 
{
nnn = nn + 1;
mmax = Math.min(n,nn);
m1min = npn6 - mmax;
m1max = npn6+mmax;
si = ssign[nl + nnn - 1];
dm1 = 0.0;
dm2 = 0.0;
for (m1 = m1min; m1 <= m1max; ++m1) 
{
m = -(npn6-m1);
if (m >= 0) 
{
sss1 = g1[m1-1][nnn-1];
}
if (m < 0) 
{
sss1 = g1[npn6-m-1][nnn-1] * si;
}
dm1 = dm1+sss1 * d1[m1-1][nn-1][n-1];
dm2 = dm2+sss1 * d2[m1-1][nn-1][n-1];
}
ffn = ff[nn-1][n-1];
sss = g1[npn6][nnn-1] * ffn;
g1l = g1l+sss * dm1;
g2l = g2l+sss * dm2 * si;
if (l >= 2) 
{
dm3 = 0.;
dm4 = 0.;
dm5r = 0.;
dm5i = 0.;
mmax = Math.min(n,nn+2);
mmin = Math.max(-n,-nn+2);
m1max = npn6+mmax;
m1min = npn6+mmin;
for (m1 = m1min; m1 <= m1max; ++m1) 
{
m = -(npn6-m1);
sss1 = g2[npn6-m-1][nnn-1];
dm3 = dm3+sss1 * d3[m1-1][nn-1][n-1];
dm4 = dm4+sss1 * d4[m1-1][nn-1][n-1];
dm5r = dm5r+sss1 * d5r[m1-1][nn-1][n-1];
dm5i = dm5i+sss1 * d5i[m1-1][nn-1][n-1];
}
g5lr = g5lr-sss * dm5r;
g5li = g5li-sss * dm5i;
sss = g2[npn4-1][nnn-1] * ffn;
g3l = g3l+sss * dm3;
g4l = g4l+sss * dm4 * si;
}
}
}
}
g1l = g1l*sl;
g2l = g2l*sl;
g3l = g3l*sl;
g4l = g4l*sl;
g5lr = g5lr*sl;
g5li = g5li*sl;
alf1[l1-1] = g1l + g2l;
alf2[l1-1] = g3l + g4l;
alf3[l1-1] = g3l - g4l;
alf4[l1-1] = g1l - g2l;
bet1[l1-1] = g5lr * 2.;
bet2[l1-1] = g5li * 2.;
lmax = l;
if (Math.abs(g1l) < 1e-6) 
{
break;
}
}
return 0;
} 

public int signum()
{
int n;
ssign[0] = 1.0;
for (n = 2; n <= 899; ++n) 
{
ssign[n - 1] = -ssign[n - 2];
}
return 0;
}

public int ccg(int n, int n1, int nmax, int k1, int k2, double gg[][])
{
int l;
int m;
int m1;
int mf;
int mm;
int nn;
int nnf;
int min;
int nnl;
int nnm;
int nnu;
int mind;
int r;
int s;
double a;
double b;
double c=0.0;
double d;
double c1;
double c2;
double cd[]=new double[npn5];
double cu[]=new double[npn5];;

if (nmax <= npn4 && 0 <= n1 && n1 <= nmax + n && n >= 1 && n <= nmax) 
{
nnf = Math.min(n + n1,nmax);
min = npn6 - n;
mf = npn6+n;
if (k1 == 1 && k2 == 0) 
{
min = npn6;
}
m=0;
for (mind = min; mind <= mf; ++mind) 
{
m=-(npn6-mind);
mm = m * k1 + k2;
m1 = mm - m;
if (Math.abs(m1) <= n1)
{
nnl = Math.max(Math.abs(mm),Math.abs(n - n1));
if (nnl <=nnf) 
{
nnu = n + n1;
nnm = (int) ((nnu + nnl) * .5);
if (nnu == nnl) 
{
nnm = nnl;
}
c=ccgin(n, n1, m, mm, c);
if(nnl-1<0)
{
cu[0] = c;
}
else
{
cu[nnl-1] = c;
}
if (nnl != nnf) 
{
c2 = 0.0;
c1 = c;
for (nn = nnl + 1; nn <= Math.min(nnm,nnf); ++nn) 
{
a = (double) ((nn + mm) * (nn - mm) * (n1 - n + nn));
a = a*(double) ((n - n1 + nn) * (n + n1 - nn + 1) * (n + n1 + nn + 1));
a = (double) (4* nn * nn) / a;
a = a*(double) ((2*nn + 1) * (2*nn - 1));
a = Math.sqrt(a);
b = (double) (m - m1) * .5;
d = 0.;
if (nn != 1) 
{
b = (double) (2*nn * (nn - 1));
b = (double) ((2*m - mm) * nn * (nn - 1) - mm * n * (n + 1) + mm * n1 * (n1 + 1)) / b;
d = (double) (4*(nn - 1) * (nn - 1));
d = d*(double) ((2*nn - 3) * (2*nn - 1));
d = (double) ((nn - mm - 1) * (nn + mm - 1) * (n1 - n + nn - 1)) / d;
d = d*(double) ((n - n1 + nn - 1) * (n + n1 - nn + 2) * (n + n1 + nn));
d = Math.sqrt(d);
}
c = a * (b * c1 - d * c2);
c2 = c1;
c1 = c;
cu[nn-1] = c;
}
if (nnf > nnm) 
{
c=direct(n, m, n1, m1, nnu, mm, c);
cd[nnu-1] = c;
if (nnu != nnm + 1) 
{
c2 = 0.;
c1 = c;
for (nn = nnu - 1; nn >= nnm + 1; --nn) 
{
a = (double) ((nn - mm + 1) * (nn + mm + 1) * (n1 - n + nn + 1));
a = a*(double) ((n - n1 + nn + 1) * (n + n1 - nn) * (n + n1 + nn + 2));
a = (double) (4*(nn + 1) * (nn + 1)) / a;
a = a*(double) ((2*nn + 1) * (2*nn + 3));
a = Math.sqrt(a);
b = (double) (2*(nn + 2) * (nn + 1));
b = (double) ((2*m - mm) * (nn + 2) * (nn + 1) - mm * n * (n + 1) + mm * n1 * (n1 + 1)) / b;
d = (double) (4*(nn + 2) * (nn + 2));
d = d*(double) ((2*nn + 5) * (2*nn + 3));
d = (double) ((nn + mm + 2) * (nn - mm + 2) * (n1 - n + nn + 2)) / d;
d = d*(double) ((n - n1 + nn + 2) * (n + n1 - nn - 1) * (n + n1 + nn + 3));
d = Math.sqrt(d);
c = a * (b * c1 - d * c2);
c2 = c1;
c1 = c;
cd[nn-1] = c;
}
}
}
}
for (nn = nnl; nn <= nnf; ++nn) 
{
if (nn <= nnm) 
{
if(nn-1<0)
{
gg[mind-1][nn] = cu[0];
}
else
{
gg[mind-1][nn] = cu[nn-1];
}
}
if (nn > nnm) 
{
gg[mind-1][nn] = cd[nn-1];
}
}
}
}
}
}
else
{
System.out.println("ERROR IN SUBROUTINE CCG");
System.exit(0);
}
return 0;
}

public double direct(int n, int m, int n1, int m1, int nn, int mm, double c)
{
int i;
int i1;
double f[]=new double[900];

f[0] = 0.0;
f[1] = 0.0;
for (i = 3; i<= 900; ++i) 
{
i1 = i - 1;
f[i - 1] = f[i1 - 1] + Math.log((double) i1) * 0.5;
}
c = f[n * 2] + f[n1 * 2] + f[n + n1 + m + m1] + f[n + n1 - m - m1];
c = c - f[(n + n1) * 2] - f[n + m] - f[n - m] - f[n1 + m1] - f[n1 - m1];
c = Math.exp(c);
return c;
}

public double ccgin(int n, int n1, int m, int mm, double g)
{
int i;
int i1;
int k;
int  l1;
int  m1;
int  l2;
int  l3;
int  n2;
int  m2;
int  m12;
int  n12;
double a;
double f[]=new double[900];

f[0] = 0.0;
f[1] = 0.0;
for (i = 3; i<= 900; ++i) 
{
i1 = i - 1;
f[i - 1] = f[i1 - 1] + Math.log((double) i1) * 0.5;
}
m1 = mm - m;
if (n >= Math.abs(m) && n1 >= Math.abs(m1) && Math.abs(mm) <= n + n1) 
{
if (Math.abs(mm) <= Math.abs(n - n1)) 
{
l1 = n;
l2 = n1;
l3 = m;
if (n1 > n) 
{
k = n;
n = n1;
n1 = k;
k = m;
m = m1;
m1 = k;
}
n2 = 2*n;
m2 = 2*m;
n12 = 2*n1;
m12 = 2*m1;
g = ssign[n1 + m1] * Math.exp(f[n + m] + f[n - m] + f[n12] + f[n2 - n12 + 1] - f[n2 + 1] - f[n1 + m1] - f[n1 - m1] -
f[n - n1 + mm] - f[n - n1 - mm]);
n = l1;
n1 = l2;
m = l3;
return g;
}
a = 1.;
l1 = m;
l2 = mm;
if (mm < 0) 
{
mm = -(mm);
m = -(m);
m1 = -m1;
a = ssign[mm + n + n1];
}
g = a * ssign[n + m] * Math.exp(f[2*(mm) + 1] + f[n +n1 - mm] + f[n + m] + f[n1 + m1] - f[n + n1 + mm + 1] -
f[n - n1 + mm] - f[-(n) + n1 + mm] - f[n - m] - f[n1 - m1]);
m = l1;
mm = l2;
return g;
}
else
{
System.out.println("ERROR IN SUBROUTINE CCGIN");
System.exit(0);
}
return 0.0;
}

public int sarea(double d)
{
double e;
double r;
if (d < 1) 
{
e = Math.sqrt(1. - d * d);
r = (Math.pow(d, 2.0/3.0) + Math.pow(d,1.0/3.0) * Math.asin(e) / e) * 0.5;
r = Math.sqrt(r);
rat = 1.0 / r;
return 0;
}
e = Math.sqrt(1.0 - 1.0 / (d * d));
r = (Math.pow(d, 2.0/3.0) * 2.0 + Math.pow(d, (-4.0/3.0)) * Math.log((e + 1.0) / (1.0 - e)) / e) * 0.25;
r = Math.sqrt(r);
rat = 1.0 / r;
return 0;
}

public int sareac(double eps)
{
rat = Math.pow(1.5 / eps, 1.0/3.0);
rat = rat/Math.sqrt((eps + 2.0) / (eps * 2.0));
return 0;
}

public int gauss(int n, int ind1, int ind2,double z[], double w[])
{
int i;
int  j;
int  k;
int  m;
int  ind;
int  niter;
int boolFlag=1;
double a;
double  b;
double  c;
double  f;
double  x=0.0;
double  dj;
double  pa;
double  pb;
double  pc;
double  zz;
double check;

a = 1.;
b = 2.;
c = 3.;
ind = n % 2;
k = n / 2 + ind;
f = (double) (n);
for (i= 1; i <= k; ++i) 
{
m = n + 1 - i;
if (i == 1) 
{
x = a - b / ((f + a) * f);
}
if (i == 2) 
{
x = (z[n-1] - a) * 4. + z[n-1];
}
if (i == 3) 
{
x = (z[n - 2] - z[n-1]) * 1.6 + z[n - 2];
}
if (i > 3) 
{
x = (z[m] - z[m + 1]) * c + z[m + 2];
}
if (i == k && ind == 1) 
{
x = 0.;
}
niter = 0;
check = 1e-16;
do{
pb = 1.;
niter=niter+1;
if (niter > 100) 
{
check = check*10.;
}
pc = x;
dj = a;
for (j = 2; j <= n; ++j) 
{
dj = dj+a;
pa = pb;
pb = pc;
pc = x * pb + (x * pb - pa) * (dj - a) / dj;
}
pa = a / ((pb - x * pc) * f);
pb = pa * pc * (a - x * x);
x = x - pb;
if (Math.abs(pb) > check * Math.abs(x)) 
{
boolFlag=1;
}
else
{
boolFlag=0;
}
}while(boolFlag==1);
z[m-1] = x;
w[m-1] = pa * pa * (a - x * x);
if (ind1 == 0) 
{
w[m-1] = b * w[m-1];
}
if (!(i == k && ind == 1)) 
{
z[i-1] = -z[m-1];
w[i-1] = w[m-1];
}
}
if (ind2 == 1) 
{
for (i = 1; i <= k; ++i) 
{
zz = -z[i-1];
}
}
if (ind1 != 0) 
{
for (i = 1; i<= n; ++i) 
{
z[i-1] = (a + z[i-1]) / b;
}
}
return 0;
}

public int distrb(int nnk, double yy[], double wy[], int ndistr, double aa, double bb, double gam, double r1, double r2,double pi)
{
int i;
double g, x, y, b2, da, xi, dab, sum;

if (ndistr == 1) 
{
b2 = (1. - bb * 3.) / bb;
dab = 1. / (aa * bb);

for (i = 1; i <= nnk; ++i) 
{
x = yy[i-1];
x = Math.pow(x, b2) * Math.exp(-x * dab);
wy[i-1] *= x;

}
}


if (ndistr == 2) 
{
for (i = 1; i <= nnk; ++i) 
{
x = yy[i-1];
y = x - aa;
y = Math.exp(-y * y * .5 / bb);
wy[i-1] *= y;
}
}


if (ndistr == 3) 
{
da = 1. / aa;
for (i = 1; i <= nnk; ++i) 
{
x = yy[i-1];
y = Math.log(x * da);
y = Math.exp(-y * y * .5 / bb) / x;
wy[i-1] *= y;
}
}
sum = 0.;
for (i = 1; i <= nnk; ++i) 
{
sum = sum+wy[i-1];
}
sum = 1. / sum;
for (i = 1; i <= nnk; ++i) 
{
wy[i-1] = wy[i-1]*sum;
}
g = 0.;
for (i = 1; i<= nnk; ++i) 
{
x = yy[i-1];
g = g+x * x * wy[i-1];
}
reff = 0.0;
for (i = 1; i <= nnk; ++i) 
{
x = yy[i-1];
reff = reff+x * x * x * wy[i-1];
}
reff /= g;
veff = 0.;
for (i = 1; i <= nnk; ++i) 
{
x = yy[i-1];
xi = x - reff;
veff = veff + xi * xi * x * x * wy[i-1];
}
veff = veff / (g * reff * reff);
return 0;
}

public int hovenr(int l1, double a1[], double a2[], double a3[], double a4[], double b1[], double b2[])
{
int i;
int  l;
int ll;
int kontr=0;

double c;
double c1;
double  c2;
double  c3;
double  cc;
double  dl;
double aa1;
double  aa2;
double  aa3;
double  aa4;
double  bb1;
double  bb2;
double  ddl;

for (l = 1; l <= l1; ++l) 
{
kontr = 1;
ll = l - 1;
dl = (double) ll * 2. + 1.;
ddl = dl * .48;
aa1 = a1[l-1];
aa2 = a2[l-1];
aa3 = a3[l-1];
aa4 = a4[l-1];
bb1 = b1[l-1];
bb2 = b2[l-1];
if (ll >= 1 && Math.abs(aa1) >= dl) 
{
kontr = 2;
}
if (Math.abs(aa2) >= dl) 
{
kontr = 2;
}
if (Math.abs(aa3) >= dl) 
{
kontr = 2;
}
if (Math.abs(aa4) >= dl) 
{
kontr = 2;
}
if (Math.abs(bb1) >= ddl) 
{
kontr = 2;
}
if (Math.abs(bb2) >= ddl) 
{
kontr = 2;
}
if (kontr == 2) 
{
System.out.println("TEST FOR VAN DER MEE & HOVENIER IS NOT SATISFIED, L="+ll);
}
c = -0.1;
for (i = 1; i <= 11; ++i) 
{
c=c+ 0.1;
cc = c * c;
c1 = cc * bb2 * bb2;
c2 = c * aa4;
c3 = c * aa3;
if ((dl - c * aa1) * (dl - c * aa2) - cc * bb1 * bb1 <= -1e-4)
 {
kontr = 2;
}
if ((dl - c2) * (dl - c3) + c1 <= -1e-4) 
{
kontr = 2;
}
if ((dl + c2) * (dl - c3) - c1 <= -1e-4) 
{
kontr = 2;
}
if ((dl - c2) * (dl + c3) - c1 <= -1e-4) 
{
kontr = 2;
}
if (kontr == 2) 
{
System.out.println("TEST FOR VAN DER MEE & HOVENIER IS NOT SATISFIED, L & A = \n"+ll+","+c);
}
}
}
if (kontr == 1) 
{
//System.out.println("TEST FOR VAN DER MEE & HOVENIER IS SATISFIED");
}
return 0;
}

public int matr(double a1[], double a2[], double a3[],double a4[], double b1[], double b2[], int lmax, int npna)
{
int l;
int n;
int i1;
int l1;
int l1max;
double p;
double u;
double f2;
double f3;
double d6;
double p1;
double p2;
double p3;
double p4;
double da;
double db;
double f11;
double f12;
double f22;
double f33;
double f34;
double f44;
double dl;
double dn;
double tb;
double dl1;
double pl1=0.0;
double pl2=0.0;
double pl3=0.0;
double pl4=0.0;
double pp1;
double pp2;
double pp3;
double pp4;
double taa;
int count=0;
double temp=0.0;

n = npna;
dn = 1. / (double) (n - 1);
da = Math.acos(-1.) * dn;
db = dn * 180.;
l1max = lmax + 1;
for (l1 = 1; l1 <= l1max; ++l1) 
{
l = l1 - 1;
}
tb = -db;
taa = -da;
d6 = Math.sqrt(6.) * .25;
for (i1 = 1; i1 <= n; ++i1) 
{
taa += da;
tb += db;
u = Math.cos(taa);
f11 = 0.;
f2 = 0.;
f3 = 0.;
f44 = 0.;
f12 = 0.;
f34 = 0.;
p1 = 0.;
p2 = 0.;
p3 = 0.;
p4 = 0.;
pp1 = 1.;
pp2 = (u + 1.) * .25 * (u + 1.);
pp3 = (1. - u) * .25 * (1. - u);
pp4 = d6 * (u * u - 1.);
L400:
for (l1 = 1; l1 <= l1max; ++l1) 
{
l = l1 - 1;
dl = (double) l;
dl1 = (double) l1;
f11 = f11 + a1[l1-1] * pp1;
f44 = f44 + a4[l1-1] * pp1;
if (l != lmax) 
{
pl1 = (double) (2*l + 1);
p = (pl1 * u * pp1 - dl * p1) / dl1;
p1 = pp1;
pp1 = p;
}
if (l < 2) 
{
continue L400;
}
f2 = f2 + (a2[l1-1] + a3[l1-1]) * pp2;
f3 = f3 + (a2[l1-1] - a3[l1-1]) * pp3;
f12 = f12 + b1[l1-1] * pp4;
f34 = f34 + b2[l1-1] * pp4;
if (l == lmax) 
{
continue L400;
}
pl2 = (double) (l * l1) * u;
pl3 = (double) (l1 * (l * l - 4));
pl4 = 1. / (double) (l * (l1 * l1 - 4));
p = (pl1 * (pl2 - 4.) * pp2 - pl3 * p2) * pl4;
p2 = pp2;
pp2 = p;
p = (pl1 * (pl2 + 4.) * pp3 - pl3 * p3) * pl4;
p3 = pp3;
pp3 = p;
p = (pl1 * u * pp4 - Math.sqrt((double) (l * l - 4)) * p4) / Math.sqrt((double) (l1 * l1 - 4));
p4 = pp4;
pp4 = p;
}
if(i1==1)
{
temp=f11;
}
f22 = (f2 + f3) * .5;
f33 = (f2 - f3) * .5;
f22 = f22 / f11;
f33 = f33 / f11;
f44 = f44 / f11;
f12 = -f12 / f11;
f34 = f34 / f11;
ss11[count]=(f11/temp);
ss12[count]=f12;
ss33[count]=f33;
ss34[count]=f34;
count++;
}
return 0;
}
}
class gspBlock{
double rt11[][][];
double rt12[][][];
double rt21[][][];
double rt22[][][];
double it11[][][];
double it12[][][];
double it21[][][];
double it22[][][];
gspBlock(int npn4,int npn6)
{
rt11=new double[npn6][npn4][npn4];
rt12=new double[npn6][npn4][npn4];
rt21=new double[npn6][npn4][npn4];
rt22=new double[npn6][npn4][npn4];
it11=new double[npn6][npn4][npn4];
it12=new double[npn6][npn4][npn4];
it21=new double[npn6][npn4][npn4];
it22=new double[npn6][npn4][npn4];
}
protected void finalize() throws Throwable
{
super.finalize();
}
}
class tmatBlock
{
double d1[][] ;
double d2[][] ;
tmatBlock(int npng2,int npn1)
{
d1 = new double[npng2][npn1];
d2= new double[npng2][npn1];
}
protected void finalize() throws Throwable
{
super.finalize();
}
}
class trtiBlock
{
static double tr1[][];
static double ti1[][];
static double tqi[][] ;
static double tqr[][] ;
trtiBlock(int npn2)
{
tr1=new double[npn2][npn2];
ti1=new double[npn2][npn2];
tqi = new double[npn2][npn2];
tqr = new double[npn2][npn2];
}
protected void finalize() throws Throwable
{
super.finalize();
}
}
class rirgigBlock
{
static double r11[][] ;
static double r12[][] ;
static double r21[][] ;
static double r22[][] ;
static double i11[][] ;
static double i12[][] ;
static double i21[][] ;
static double i22[][] ;
static double rg11[][] ;
static double rg12[][] ;
static double rg21[][] ;
static double rg22[][] ;
static double ig11[][] ;
static double ig12[][] ;
static double ig21[][] ;
static double ig22[][] ;
rirgigBlock(int npn1)
{
r11 = new double[npn1][npn1];
r12 = new double[npn1][npn1];
r21 = new double[npn1][npn1];
r22 = new double[npn1][npn1];
i11 = new double[npn1][npn1];
i12 = new double[npn1][npn1];
i21 = new double[npn1][npn1];
i22 = new double[npn1][npn1];
rg11 = new double[npn1][npn1];
rg12 = new double[npn1][npn1];
rg21 = new double[npn1][npn1];
rg22 = new double[npn1][npn1];
ig11 = new double[npn1][npn1];
ig12 = new double[npn1][npn1];
ig21 = new double[npn1][npn1];
ig22 = new double[npn1][npn1];
}
protected void finalize() throws Throwable
{
super.finalize();
}
}
class cbassBlock
{
static double j[][];
static double y[][];
static double jr[][];
static double ji[][];
static double dj[][];
static double dy[][];
static double djr[][];
static double dji[][];
cbassBlock(int npng2,int npn1)
{
j=new double[npng2][npn1];
y=new double[npng2][npn1];
jr=new double[npng2][npn1];
ji=new double[npng2][npn1];
dj=new double[npng2][npn1];
dy=new double[npng2][npn1];
djr=new double[npng2][npn1];
dji=new double[npng2][npn1];
}
}