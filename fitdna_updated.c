// Fitting supercoiled beads to plectonemic model, Asli Yildirim 2016

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>

#include "vector.h"
#include "pdb.h"

double radius=5.0;
double openingangle=1.073377; // default=61.5 degrees
double PI=3.14159265;
double bead=5.0;
double rise=0.340;
int bpbead=(int)(bead/rise+0.5);
double kT=0.00005;
int maxtrials=20000;
double rfac;
double rfac2;
int natom;
int deltabead;
int target;
int final;
int action;
int ndomain;
int nbranch;
int extrabead=0;

int mbps;                      // number of basepairs

double bpssingle=1.0/rise;      // bps per nm
double bpsscoil=2.0/rise/0.879;   // bps per nm

double setringlen=20;       // nm
double setchainlen;    // nm
double halfchainlen;

int cgbeadperchain;
int cgbeadperring;
int cgbeadforhalfchain;

class Bond;

class Bead {
  public:
    int inx;
    double radius;
    double xcoor,ycoor,zcoor;
    Bead *link[3];
    int nlink;
    int chlinkinx;
    Bond *bond[3];
    int nbond;
    int type; // 1: ring, 2: branch
  public:
    Bead(Vector p, int t, int i, double r) :
      xcoor(p.x()), ycoor(p.y()), zcoor(p.z()), type(t), inx(i),radius(r),nlink(0),chlinkinx(0),nbond(0) {}
};

class Bond {
   public:
     Bead *b1;
     Bead *b2;
     double dx,dy,dz;
     double d;
     double *xfit,*yfit,*zfit;
     int nfit;
     int maxfit;
     int assignedstrand;
     int annealforwarddone;
     int annealbackwarddone;
     double phase;
     int sfrom[2];
     int sto[2];
     int done[2];
     Bond() : nfit(0),maxfit(1000),phase(0.0),assignedstrand(-1),annealforwarddone(0),annealbackwarddone(0) { 
       xfit=new double[maxfit]; 
       yfit=new double[maxfit]; 
       zfit=new double[maxfit]; 
       done[0]=done[1]=0;
     }
     Bond(Bead *bead1, Bead *bead2) : nfit(0),maxfit(1000),phase(0.0),assignedstrand(-1),annealforwarddone(0),annealbackwarddone(0),b1(bead1),b2(bead2) {
       xfit=new double[maxfit]; 
       yfit=new double[maxfit]; 
       zfit=new double[maxfit]; 
       dx=b2->xcoor-b1->xcoor;
       dy=b2->ycoor-b1->ycoor;
       dz=b2->zcoor-b1->zcoor;
       if (dx==0 && dy==0) {
         fprintf(stderr,"Warning! Had to play with the x/y coordinates of one plectonemic bond.\n");
         dx=0.01;
         dy=0.01;
       }
       d=sqrt(dx*dx+dy*dy+dz*dz);
       done[0]=done[1]=0;
     }
     ~Bond() {
       delete xfit;
       delete yfit;
       delete zfit;
     }
};

class Bondconnection {
   public:
     Bond *b1;
     Bond *b2;
     int fromstrand;
     int tostrand;
     int fromdirection;
     int todirection;
     int fromfitinx;
     int tofitinx;
   public:
     Bondconnection(Bond *bond1, Bond *bond2) : fromstrand(-1),tostrand(-1),fromdirection(-1),todirection(-1),fromfitinx(-999),tofitinx(-999),b1(bond1),b2(bond2) {}
}; 

class Connection {
  public:
    int inx;
    int jnx;
};

int getnum(char *cptr, int n) {
  char buf[512];
  strncpy(buf,cptr,n);
  buf[n]=0;
  return atoi(buf);
}

void readPDB(char *pdbfilename, PDBEntry *pdb, int& natom, Connection *connect, int& nconnect, int& ndomain) {
     FILE *fptr;
     int lastnum=-1;
     natom=0;
     ndomain=0;
     fptr=fopen(pdbfilename,"r");
     if (fptr==0) {
       fprintf(stderr,"Cannot open PDB file %s\n",pdbfilename);
       exit(1);
     }
     while (!feof(fptr)) {
       if (pdb[natom].read(fptr)>0 && pdb[natom].residueNumber()!=lastnum) {
           lastnum=pdb[natom].residueNumber();
           if (!strcmp(pdb[natom].atomName(),"C")) {
             ndomain++;
           }
           natom++;
       }
     }
     fclose(fptr);
     nconnect=0;
     char line[1024];
     fptr=fopen(pdbfilename,"r");
     while (!feof(fptr)) {
       if (fgets(line,1024,fptr)) {
         if (!strncmp(line,"CONECT",6)) {
           int inx=atoi(PDBEntry::substr(line,6,5));
           int jnx=atoi(PDBEntry::substr(line,11,5));
           connect[nconnect].inx=inx;
           connect[nconnect].jnx=jnx;
           nconnect++;
         }
       }
     }
     fclose(fptr);
}

void transform(double &x, double &y, double &z, double dx, double dy, double dz) {
  double m11=dz*dx/sqrt(1-dz*dz);
  double m12=-dy/sqrt(1-dz*dz);
  double m13=dx;
  double m21=dz*dy/sqrt(1-dz*dz);
  double m22=dx/sqrt(1-dz*dz); 
  double m23=dy;
  double m31=-sqrt(1-dz*dz);
  double m32=0;
  double m33=dz;
  double nx=(m11*x+m12*y+m13*z);
  double ny=(m21*x+m22*y+m23*z);
  double nz=(m31*x+m32*y+m33*z);
//  fprintf(stderr,"m11: %lf m12: %lf m13: %lf m21: %lf m22: %lf m23: %lf m31: %lf m32: %lf m33: %lf\n",m11,m12,m13,m21,m22,m23,m31,m32,m33);
  x=nx; y=ny; z=nz;  
}

void writePDB(int n, int nres, double x, double y, double z, char *atomname,char c) {
  if (n<=9999) {
    fprintf(stdout,"ATOM %6d  %-3s%c%-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",n,atomname,' ',"BD",c,n,x/10,y/10,z/10,0.0,0.0,"");
  } else {
    if (n<=99999) {
      fprintf(stdout,"ATOM %6d  %-3s%c%-4s%c%5d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",n,atomname,' ',"BD",c,n,x/10,y/10,z/10,0.0,0.0,"");
    } else {
      if (n<=999999) {
        fprintf(stdout,"ATOM %6d  %-3s%c%-4s%c%6d  %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",n,atomname,' ',"BD",c,n,x/10,y/10,z/10,0.0,0.0,"");
      } else {
          fprintf(stdout,"ATOM %6d  %-3s%c%-4s%c%7d %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",n-1000000,atomname,' ',"BD",c,n,x/10,y/10,z/10,0.0,0.0,"");
        }
      }
    }
}

void fitDNA(Bond *b) {
  double angle,x,y,z;
  double dd;
  b->nfit=0;
  b->sfrom[0]=b->nfit;
  double fitbead;
  double fitradius;
  fitbead=bead;
  if (b->b2->inx!=ndomain+nbranch-1) {
    bpsscoil=int(cgbeadperchain-0.5)*(int)(fitbead*bpssingle+0.5)/b->d;
  } else {
    bpsscoil=(cgbeadforhalfchain+int(cgbeadperchain-0.5))*(int)(fitbead*bpssingle+0.5)/b->d;
  }
  rfac=2/rise/bpsscoil;
  if (rfac>1) {
    do {
      fitbead=fitbead+0.5;
      if (b->b2->inx!=ndomain+nbranch-1) {
        bpsscoil=cgbeadperchain*(int)(fitbead*bpssingle+0.5)/b->d;
      } else {
        bpsscoil=(cgbeadforhalfchain+int(cgbeadperchain-0.5))*(int)(fitbead*bpssingle+0.5)/b->d;
      }
      rfac=2/rise/bpsscoil;
    } while (rfac>1); 
  }
  openingangle=asin(rfac);
  rfac2=cos(openingangle)/sin(openingangle);
//  fprintf(stderr,"%d:%d bead %.1f rfac %.2f rfac2 %.2f openingangle %.2f\n",b->b1->inx,b->b2->inx,fitbead,rfac,rfac2,openingangle*180/PI);
//  fprintf(stderr,"%d :: %.3f %.3f %.3f %d :: %.3f %.3f %.3f\n",b->b1->inx,b->b1->xcoor,b->b1->ycoor,b->b1->zcoor,b->b2->inx,b->b2->xcoor,b->b2->ycoor,b->b2->zcoor);
  double ddelta=0.0;
  fitradius=fitbead;
  if (b->b2->nlink==1) {
    ddelta=fitradius/2;
  }
  for (dd=0; dd<b->d-ddelta; dd+=fitbead*rfac) {
    angle=rfac2*dd/fitradius+b->phase;
    x=fitradius*cos(angle);
    y=fitradius*sin(angle);
    z=dd;
//    fprintf(stderr,"Before transform %lf %lf %lf %lf %lf %lf %lf\n",x,y,z,b->dx,b->dy,b->dz,b->d);
    transform(x,y,z,b->dx/b->d,b->dy/b->d,b->dz/b->d); 
    b->xfit[b->nfit]=x+b->b1->xcoor;
    b->yfit[b->nfit]=y+b->b1->ycoor;
    b->zfit[b->nfit]=z+b->b1->zcoor;
//    fprintf(stderr,"Case1 %d %lf %lf %lf\n",b->nfit,b->xfit[b->nfit],b->yfit[b->nfit],b->zfit[b->nfit]);
    b->nfit++;
    if (b->b2->nlink!=1) {
      if (b->nfit==cgbeadperchain/2 && b->b2->inx!=ndomain+nbranch-1) { break; }
    } else {
      if (b->b2->inx!=ndomain+nbranch-1) {
        if (b->nfit==(cgbeadperchain/2)-2) { break; }
      } else {
        if (b->nfit==(cgbeadperchain/2)+(int)(cgbeadforhalfchain/2)+1-2) { break; }
      }
    } 
  }
//  fprintf(stderr,"1.fitdna %d:%d %d\n",b->b1->inx,b->b2->inx,b->nfit);
  if (b->b2->nlink==1) {
    double start;
    double add;
/*  if (difference==0) {
      start=0.3;
      add=1.27;
    } else if (difference==1) {
      start=0.785;
      add=1.57;
    } else if (difference==2) { 
      start=1.57;
      add=3.14;
    } else if (difference==-1) { 
      start=0.2;
      add=0.913;
    } else if (difference==-2) { 
      start=0.1;
      add=0.735;
    } else if (difference==-3) {
      start=0.1;
      add=0.588;
    }*/
    start=0.2;
    add=0.913;
    double lastz=dd;
    angle=rfac2*lastz/(fitradius)+b->phase;
    double theta;
    for (theta=start; theta<=PI/2.0; theta+=add) {
      double x=fitradius*cos(angle)*sin(PI/2.0-theta);
      double y=fitradius*sin(angle)*sin(PI/2.0-theta);
      double z=fitradius*cos(PI/2.0-theta)+lastz;
//      fprintf(stderr,"Before transform %lf %lf %lf %lf %lf %lf %lf\n",x,y,z,b->dx,b->dy,b->dz,b->d);
      transform(x,y,z,b->dx/b->d,b->dy/b->d,b->dz/b->d); 
      b->xfit[b->nfit]=x+b->b1->xcoor;
      b->yfit[b->nfit]=y+b->b1->ycoor;
      b->zfit[b->nfit]=z+b->b1->zcoor;
//      fprintf(stderr,"Case2 %d %lf %lf %lf\n",b->nfit,b->xfit[b->nfit],b->yfit[b->nfit],b->zfit[b->nfit]);
      b->nfit++;
    }
    for (; theta<=PI-start; theta+=add) {
      double x=fitradius*cos(angle+PI)*sin(theta-PI/2.0);
      double y=fitradius*sin(angle+PI)*sin(theta-PI/2.0);
      double z=fitradius*cos(theta-PI/2.0)+lastz;
//      fprintf(stderr,"Before transform %lf %lf %lf %lf %lf %lf %lf\n",x,y,z,b->dx,b->dy,b->dz,b->d);
      transform(x,y,z,b->dx/b->d,b->dy/b->d,b->dz/b->d); 
      b->xfit[b->nfit]=x+b->b1->xcoor;
      b->yfit[b->nfit]=y+b->b1->ycoor;
      b->zfit[b->nfit]=z+b->b1->zcoor;
//      fprintf(stderr,"Case3 %d %lf %lf %lf\n",b->nfit,b->xfit[b->nfit],b->yfit[b->nfit],b->zfit[b->nfit]);
      b->nfit++;
    }
  }
//  fprintf(stderr,"2.fitdna %d:%d %d\n",b->b1->inx,b->b2->inx,b->nfit);
  b->sto[0]=b->nfit-1;
  b->sfrom[1]=b->nfit;
  for (; dd>=-0.0001; dd-=fitbead*rfac) {
    angle=rfac2*dd/(fitradius)+b->phase;
    x=fitradius*cos(angle+PI);
    y=fitradius*sin(angle+PI);
    z=dd;
//    fprintf(stderr,"Before transform %lf %lf %lf %lf %lf %lf %lf\n",x,y,z,b->dx,b->dy,b->dz,b->d);
    transform(x,y,z,b->dx/b->d,b->dy/b->d,b->dz/b->d);
    b->xfit[b->nfit]=x+b->b1->xcoor;
    b->yfit[b->nfit]=y+b->b1->ycoor;
    b->zfit[b->nfit]=z+b->b1->zcoor;
//    fprintf(stderr,"Case4 %d %lf %lf %lf\n",b->nfit,b->xfit[b->nfit],b->yfit[b->nfit],b->zfit[b->nfit]);
    b->nfit++;
    if (b->nfit==cgbeadperchain+cgbeadforhalfchain && b->b2->inx==ndomain+nbranch-1) { break; }
    if (b->nfit==cgbeadperchain && b->b2->inx!=ndomain+nbranch-1) { break; }
  }
  b->sto[1]=b->nfit-1;
//  fprintf(stderr,"3.fitdna %d:%d %d\n",b->b1->inx,b->b2->inx,b->nfit);
}

void fitLinearDNA(Bond *b) {
  double x,y,z;
  double dd;
  b->nfit=0;
  b->sfrom[0]=b->nfit;
  double distincrement=(b->d)/3;
  double incr;
  for (dd=0; dd<cgbeadperring; dd++) {
    incr=dd*distincrement;
    x=0;
    y=0;
    z=incr;
    transform(x,y,z,b->dx/b->d,b->dy/b->d,b->dz/b->d); 
    b->xfit[b->nfit]=x+b->b1->xcoor;
    b->yfit[b->nfit]=y+b->b1->ycoor;
    b->zfit[b->nfit]=z+b->b1->zcoor;
    b->nfit++;
  }
  b->sto[0]=b->nfit-1;
}

double koff=0.0;
double kcontact=1.0;
double getEnergy(int nlookup, int *lookup, Bond **b) {
  double energy=0.0;
  for (int i=0; i<nlookup; i+=4) {
    int bi=lookup[i];
    int ii=lookup[i+1];
    int bj=lookup[i+2];
    int ij=lookup[i+3];
    double dx=b[bi]->xfit[ii]-b[bj]->xfit[ij];
    double dy=b[bi]->yfit[ii]-b[bj]->yfit[ij];
    double dz=b[bi]->zfit[ii]-b[bj]->zfit[ij];
    double d=dx*dx+dy*dy+dz*dz;
    energy+=kcontact*exp(-(sqrt(d)+koff));
  } 
  return energy;
}

void findClosest(Bond *b1, Bond *b2, int s0, int s1, double &mind, int &mini, int &minj) {
  mind=999999;
  mini=-1;
  minj=-1;
  Bead *b1bead;
  Bead *b2bead;
  double ax,ay,az,a;
  double dx,dy,dz,d,tx,ty,tz,t; 
  double bx,by,bz,b,px,py,pz,p;
  if (b1->b2->inx == b2->b1->inx) {
    b1bead=b1->b2;
    b2bead=b2->b1;
  } else if (b1->b2->inx == b2->b2->inx) {
    b1bead=b1->b2;
    b2bead= b2->b2;
  } else if (b1->b1->inx == b2->b1->inx) {
    b1bead=b1->b1;
    b2bead= b2->b1;
  } else {
    b1bead=b1->b1;
    b2bead=b2->b2;
  }
  dx=b1bead->xcoor-b1->xfit[b1->sfrom[s0]];
  dy=b1bead->ycoor-b1->yfit[b1->sfrom[s0]];
  dz=b1bead->zcoor-b1->zfit[b1->sfrom[s0]];
  tx=b1bead->xcoor-b1->xfit[b1->sto[s0]];
  ty=b1bead->ycoor-b1->yfit[b1->sto[s0]];
  tz=b1bead->zcoor-b1->zfit[b1->sto[s0]];
  d=dx*dx+dy*dy+dz*dz;
  t=tx*tx+ty*ty+tz*tz;
  bx=b2bead->xcoor-b2->xfit[b2->sfrom[s1]];
  by=b2bead->ycoor-b2->yfit[b2->sfrom[s1]];
  bz=b2bead->zcoor-b2->zfit[b2->sfrom[s1]];
  px=b2bead->xcoor-b2->xfit[b2->sto[s1]];
  py=b2bead->ycoor-b2->yfit[b2->sto[s1]];
  pz=b2bead->zcoor-b2->zfit[b2->sto[s1]];
  b=bx*bx+by*by+bz*bz;
  p=px*px+py*py+pz*pz; 
  if (d<t && b<p) {
    for (int i=b1->sfrom[s0]; i<=b1->sfrom[s0]+0; i++) {
      for (int j=b2->sfrom[s1]; j<=b2->sfrom[s1]+0; j++) {
        ax=b1->xfit[i]-b2->xfit[j];
        ay=b1->yfit[i]-b2->yfit[j];
        az=b1->zfit[i]-b2->zfit[j];
        a=ax*ax+ay*ay+az*az;
        if (a<mind) {
          mind=a;
          mini=i;
          minj=j;
        }
      }
    }
  } else if (d<t && b>p) {
     for (int i=b1->sfrom[s0]; i<=b1->sfrom[s0]+0; i++) {
        for (int j=b2->sto[s1]-0; j<=b2->sto[s1]; j++) {
          ax=b1->xfit[i]-b2->xfit[j];
          ay=b1->yfit[i]-b2->yfit[j];
          az=b1->zfit[i]-b2->zfit[j];
          a=ax*ax+ay*ay+az*az;
          if (a<mind) {
          mind=a;
          mini=i;
          minj=j;
          }
        }
      }
  } else if (d>t && b<p) {
     for (int i=b1->sto[s0]-0; i<=b1->sto[s0]; i++) {
       for (int j=b2->sfrom[s1]; j<=b2->sfrom[s1]+0; j++) {
         ax=b1->xfit[i]-b2->xfit[j];
         ay=b1->yfit[i]-b2->yfit[j];
         az=b1->zfit[i]-b2->zfit[j];
         a=ax*ax+ay*ay+az*az;
         if (a<mind) {
           mind=a;
           mini=i;
           minj=j;
         }
       }
     }
   } else {
     for (int i=b1->sto[s0]-0; i<=b1->sto[s0]; i++) {
       for (int j=b2->sto[s1]-0; j<=b2->sto[s1]; j++) {
         ax=b1->xfit[i]-b2->xfit[j];
         ay=b1->yfit[i]-b2->yfit[j];
         az=b1->zfit[i]-b2->zfit[j];
         a=ax*ax+ay*ay+az*az;
         if (a<mind) {
           mind=a;
           mini=i;
           minj=j;
         }
       }
     }
  }
}

void enterBond(Bond *b1, Bond *b2, int sf, int st, int df, int dt, int mf, int mt, int prevmt, int &natom) {
  int nfitatom=0;
  if (df==1) {
    for (int i=prevmt; i<=mf; i++) {
      if (final==1) {
        writePDB(natom++,1,b1->xfit[i],b1->yfit[i],b1->zfit[i],"CG",'A');
      } else {
        natom++;
        nfitatom++;
      }
    }
  } else {
    for (int i=prevmt; i>=mf; i--) {
      if (final==1) {
        writePDB(natom++,1,b1->xfit[i],b1->yfit[i],b1->zfit[i],"CG",'A');    
      } else {
        natom++;
        nfitatom++;
      }
    }
  }
//  fprintf(stderr,"enterbond %d:%d - %d:%d nfitatom %d\n",b1->b1->inx,b1->b2->inx,b2->b1->inx,b2->b2->inx,nfitatom);
//  fprintf(stderr,"enterbond sf %d st %d df %d dt %d mf %d mt %d prevmt %d\n",sf,st,df,dt,mf,mt,prevmt); 
}

void annealforwardforward(Bondconnection *c1) {
//  fprintf(stderr,"annealforwardforward %d:%d - %d:%d\n",c1->b1->b1->inx,c1->b1->b2->inx,c1->b2->b1->inx,c1->b2->b2->inx);
  double mind11,mind12;
  int mini11,mini12;
  int minj11,minj12;
  if (c1->b1->assignedstrand==0) {
    findClosest(c1->b1,c1->b2,0,0,mind11,mini11,minj11);
    findClosest(c1->b1,c1->b2,0,1,mind12,mini12,minj12);
    if (mind11<mind12) {
      c1->b2->assignedstrand=0;
      c1->fromstrand=0;
      c1->tostrand=0;
      c1->fromdirection=1;
      c1->todirection=1;
      c1->fromfitinx=mini11;
      c1->tofitinx=minj11;
    } else {
      c1->b2->assignedstrand=1;
      c1->fromstrand=0;
      c1->tostrand=1;
      c1->fromdirection=1;
      c1->todirection=0;
      c1->fromfitinx=mini12;
      c1->tofitinx=minj12;
    }
  } else {
    findClosest(c1->b1,c1->b2,1,0,mind11,mini11,minj11);
    findClosest(c1->b1,c1->b2,1,1,mind12,mini12,minj12);    
    if (mind11<mind12) {
      c1->b2->assignedstrand=0;
      c1->fromstrand=1;
      c1->tostrand=0;
      c1->fromdirection=0;
      c1->todirection=1;
      c1->fromfitinx=mini11;
      c1->tofitinx=minj11;
    } else {
      c1->b2->assignedstrand=1;
      c1->fromstrand=1;
      c1->tostrand=1;
      c1->fromdirection=0;
      c1->todirection=0;
      c1->fromfitinx=mini12;
      c1->tofitinx=minj12; 
    }
  }
}

void annealbackwardforward(Bondconnection *c1) {
//  fprintf(stderr,"annealbackwardforward %d:%d - %d:%d\n",c1->b1->b1->inx,c1->b1->b2->inx,c1->b2->b1->inx,c1->b2->b2->inx);
  double mind11,mind12;
  int mini11,mini12;
  int minj11,minj12;
  if (c1->b1->assignedstrand==0) {
    findClosest(c1->b1,c1->b2,1,0,mind11,mini11,minj11);
    findClosest(c1->b1,c1->b2,1,1,mind12,mini12,minj12);
    if (mind11<mind12) {
      c1->b2->assignedstrand=0;
      c1->fromstrand=1;
      c1->tostrand=0;
      c1->fromdirection=1;
      c1->todirection=1;
      c1->fromfitinx=mini11;
      c1->tofitinx=minj11;
    } else {
      c1->b2->assignedstrand=1;
      c1->fromstrand=1;
      c1->tostrand=1;
      c1->fromdirection=1;
      c1->todirection=0;
      c1->fromfitinx=mini12;
      c1->tofitinx=minj12;
    }
  } else {
    findClosest(c1->b1,c1->b2,0,0,mind11,mini11,minj11);
    findClosest(c1->b1,c1->b2,0,1,mind12,mini12,minj12);
    if (mind11<mind12) {
      c1->b2->assignedstrand=0;
      c1->fromstrand=0;
      c1->tostrand=0;
      c1->fromdirection=0;
      c1->todirection=1;
      c1->fromfitinx=mini11;
      c1->tofitinx=minj11;
    } else {
      c1->b2->assignedstrand=1;
      c1->fromstrand=0;
      c1->tostrand=1;
      c1->fromdirection=0;
      c1->todirection=0;
      c1->fromfitinx=mini12;
      c1->tofitinx=minj12;
    }
  }
}

void annealcap(Bondconnection *c1) {
//  fprintf(stderr,"annealcap %d:%d - %d:%d\n",c1->b1->b1->inx,c1->b1->b2->inx,c1->b2->b1->inx,c1->b2->b2->inx);
  if (c1->b1->assignedstrand==0) {
      c1->fromstrand=0;
      c1->tostrand=1;
      c1->fromdirection=1;
      c1->todirection=1;
      c1->fromfitinx=c1->b1->sto[0];
      c1->tofitinx=c1->b1->sfrom[1];
  } else {
      c1->fromstrand=1;
      c1->tostrand=0;
      c1->fromdirection=0;
      c1->todirection=0;
      c1->fromfitinx=c1->b1->sfrom[1];
      c1->tofitinx=c1->b1->sto[0];
  }
}

void annealbackward(Bondconnection *c1) {
//  fprintf(stderr,"annealbackward %d:%d - %d:%d\n",c1->b1->b1->inx,c1->b1->b2->inx,c1->b2->b1->inx,c1->b2->b2->inx);
  double mind;
  int mini,minj;
  if ( c1->b1->assignedstrand==0 && c1->b2->assignedstrand==0) {
    findClosest(c1->b1,c1->b2,1,1,mind,mini,minj);
    c1->fromstrand=1;
    c1->tostrand=1;
    c1->fromdirection=1;
    c1->todirection=1;
    c1->fromfitinx=mini;
    c1->tofitinx=minj;
  } else if (c1->b1->assignedstrand==0 && c1->b2->assignedstrand==1) {
    findClosest(c1->b1,c1->b2,1,0,mind,mini,minj);
    c1->fromstrand=1;
    c1->tostrand=0;
    c1->fromdirection=1;
    c1->todirection=0;
    c1->fromfitinx=mini;
    c1->tofitinx=minj;
  } else if (c1->b1->assignedstrand==1 && c1->b2->assignedstrand==0) { 
    findClosest(c1->b1,c1->b2,0,1,mind,mini,minj);
    c1->fromstrand=0;
    c1->tostrand=1;
    c1->fromdirection=0;
    c1->todirection=1;
    c1->fromfitinx=mini;
    c1->tofitinx=minj;
  } else {
    findClosest(c1->b1,c1->b2,0,0,mind,mini,minj);
    c1->fromstrand=0;
    c1->tostrand=0;
    c1->fromdirection=0;
    c1->todirection=0;
    c1->fromfitinx=mini;
    c1->tofitinx=minj;
  }
}

void annealRingforward(Bondconnection *c1) {
//  fprintf(stderr,"annealringforward %d:%d - %d:%d\n",c1->b1->b1->inx,c1->b1->b2->inx,c1->b2->b1->inx,c1->b2->b2->inx);
  int mini00,minj00;
  double mind00;
  findClosest(c1->b1,c1->b2,0,0,mind00,mini00,minj00);
  int mini01,minj01;
  double mind01;
  findClosest(c1->b1,c1->b2,0,1,mind01,mini01,minj01);
  if (mind00<mind01) {
    c1->fromstrand=0;
    c1->tostrand=0;
    c1->fromdirection=1;
    c1->todirection=1;
    c1->fromfitinx=mini00;
    c1->tofitinx=minj00;
    c1->b2->assignedstrand=0;
  } else {
    c1->fromstrand=0;
    c1->tostrand=1;
    c1->fromdirection=1;
    c1->todirection=0;
    c1->fromfitinx=mini01;
    c1->tofitinx=minj01;
    c1->b2->assignedstrand=1;
  }
} 

void annealRingbackward(Bondconnection *c1) {
//  fprintf(stderr,"annealringbackward %d:%d - %d:%d\n",c1->b1->b1->inx,c1->b1->b2->inx,c1->b2->b1->inx,c1->b2->b2->inx);
  if (c1->b1->assignedstrand==0) {
    int mini10,minj10;
    double mind10;
    findClosest(c1->b1,c1->b2,1,0,mind10,mini10,minj10);
    c1->fromstrand=1;
    c1->tostrand=0;
    c1->fromdirection=1;
    c1->todirection=1;
    c1->fromfitinx=mini10;
    c1->tofitinx=minj10;
  } else {
    int mini00,minj00;
    double mind00;
    findClosest(c1->b1,c1->b2,0,0,mind00,mini00,minj00);
    c1->fromstrand=0;
    c1->tostrand=0;
    c1->fromdirection=0;
    c1->todirection=1;
    c1->fromfitinx=mini00;
    c1->tofitinx=minj00;
  }
}

int main(int argc, char **argv) {

  srandom(time(NULL));
  srandom(0);
  PDBEntry *dnapdb;
  dnapdb = new PDBEntry[100000];
  int natompdb;
  Connection *connectpdb;
  connectpdb=new Connection[100000];
  int nconnect;
  ndomain=0;
  readPDB(argv[1], dnapdb, natompdb, connectpdb, nconnect, ndomain);
  setchainlen=atoi(argv[2]);

  mbps=atoi(argv[3]);

  nbranch=natompdb-ndomain;
  Bead **ringbead=new Bead*[ndomain];
  Bead **scbead=new Bead*[nbranch];
  for (int i=0; i<ndomain; i++) {
    ringbead[i]=new Bead(dnapdb[i].coordinates()*10.0,1,i,radius);
  }
  for (int i=0; i<nbranch; i++) {
    scbead[i]=new Bead(dnapdb[i+ndomain].coordinates()*10.0,2,i+ndomain,radius);
  }
  Bond **ringbond=new Bond*[ndomain];
  Bond **scbond=new Bond*[nbranch];
  for (int i=0; i<ndomain; i++) {
    Bead *next=(i<ndomain-1)?ringbead[i+1]:ringbead[0];
    ringbond[i]=new Bond(ringbead[i],next);
    ringbead[i]->bond[ringbead[i]->nlink]=ringbond[i];
    ringbead[i]->link[0]=next;
    ringbead[i]->nlink++;
    next->bond[next->nlink]=ringbond[i];
    next->link[1]=ringbead[i];
    next->nlink++;
  }
  int nscbond=0;
  for (int i=0; i<nconnect; i++) {
    if (connectpdb[i].jnx>ndomain) {
      int inx=connectpdb[i].inx;
      int jnx=connectpdb[i].jnx;
      Bead *prev=(inx<=ndomain)?ringbead[inx-1]:scbead[inx-1-ndomain];
      Bead *next=scbead[jnx-1-ndomain];
      scbond[nscbond]=new Bond(prev,next);
      prev->bond[prev->nlink]=scbond[nscbond]; 
      prev->link[prev->nlink++]=next; 
      next->bond[next->nlink]=scbond[nscbond];
      next->link[next->nlink++]=prev; 
      nscbond++;
    }
  }
  delete connectpdb;
  delete dnapdb;

  maxtrials=10;
  int headinx=0; 
  for (int i=0; i<nscbond; i++) {
    if (scbond[i]->b2->nlink==1 || scbond[i]->b1->nlink==1) {
      headinx++;
    }
  }
  int head[headinx];
  int totalfits;

  int totalcgbeads=(int)(mbps/(int)(bead*bpssingle+0.5)); // totalcgbeads: 269529
  int totalringcgbeads=ndomain*setringlen/bead; // 200*20/5
  int totalchaincgbeads=totalcgbeads-totalringcgbeads; //269529-200*20/5
  
//  fprintf(stderr,"ndomain %d nbranch %d setchainlen %.1f\n",ndomain,nbranch,setchainlen); 
//  fprintf(stderr,"totalcgbeads %d totalringcgbeads %d totalchaincgbeads %d\n",totalcgbeads,totalringcgbeads,totalchaincgbeads);

  cgbeadperring=totalringcgbeads/ndomain; // 20/5
  cgbeadperchain=(int)((setchainlen*bpsscoil/(int)(bead*bpssingle+0.5))+0.5); //40*6.69/15
  if (cgbeadperchain%2!=0) {
   cgbeadperchain++;
  }
  
  cgbeadforhalfchain=totalcgbeads-(nbranch*cgbeadperchain+ndomain*cgbeadperring);
  halfchainlen=(double)(cgbeadforhalfchain*setchainlen/cgbeadperchain);
  
//  fprintf(stderr,"cgbeadperring %d cgbeadperchain %d beadforhalfchain %d halfchainlen %.2f\n",cgbeadperring,cgbeadperchain,cgbeadforhalfchain,halfchainlen);

  totalfits=0;
  int hinx=0;
  for (int i=0; i<nscbond; i++) {
    fitDNA(scbond[i]);
    if (scbond[i]->b2->nlink==1 || scbond[i]->b1->nlink==1) {
      head[hinx]=i;
      hinx++;
    }
    totalfits+=scbond[i]->nfit;
  }
  int ringtotalfits=0;
  for (int i=0; i<ndomain; i++) {
    fitLinearDNA(ringbond[i]);
    ringtotalfits+=ringbond[i]->nfit;
  }

  int natom=1;
  std::vector<int> alltheway;
  Bond **bondsall=new Bond*[nscbond+ndomain];
  Bondconnection **connectedbonds=new Bondconnection*[2*nscbond+ndomain];
  int connectioninx=0;
  for (int j=0; j<ndomain; j++) {
    bondsall[j]=ringbond[j];
  }
  for (int j=0; j<nscbond; j++) {
    bondsall[j+ndomain]=scbond[j];
  }
  for (int i=0; i<ndomain; i++) {
    Bead *rb=ringbead[i];
    Bead *rbplus;
    if (i==ndomain-1) {
      rbplus=ringbead[0];
    } else {
      rbplus=ringbead[i+1];
    }
    int scinx=rbplus->bond[2]->b2->inx;
    Bead *b=scbead[scinx-ndomain];
    int last=0;
    alltheway.push_back(rb->inx);
    alltheway.push_back(rbplus->inx);
    alltheway.push_back(b->inx);
    if (i==0) {
      connectedbonds[connectioninx]=new Bondconnection(rb->bond[0],rbplus->bond[2]);
      annealRingforward(connectedbonds[connectioninx]);
      connectioninx++;
    } else {
      connectedbonds[connectioninx]=new Bondconnection(rb->bond[1],rbplus->bond[2]);
      annealRingforward(connectedbonds[connectioninx]);
      connectioninx++;        
    }
    Bead *whosnext=b;
    Bead *prev=rbplus;
    Bond *prevbond=rbplus->bond[2];
    prevbond->annealforwarddone=1;
    Bond *nextbond;
    do {
      int nextinx;
      if (whosnext->nlink==3) {
        int indices[6];
        indices[0]=whosnext->bond[0]->b1->inx;
        indices[1]=whosnext->bond[0]->b2->inx;
        indices[2]=whosnext->bond[1]->b1->inx;
        indices[3]=whosnext->bond[1]->b2->inx;
        indices[4]=whosnext->bond[2]->b1->inx;
        indices[5]=whosnext->bond[2]->b2->inx;
        int cand1;
        int cand2;
        int cand1nlink;
        int cand2nlink;
        int grade1;
        int grade2;
        int candcount=0;
        for (int i=0; i<6; i++) {
          if (indices[i]!=whosnext->inx && indices[i]!=prev->inx && candcount==0) {
            cand1=indices[i];
            grade1 = std::count( alltheway.begin(), alltheway.end(), cand1);
            candcount++;
          } else if (indices[i]!=whosnext->inx && indices[i]!=prev->inx && candcount==1) {
            cand2=indices[i];
            grade2 = std::count( alltheway.begin(), alltheway.end(), cand2);
          }
        }
        if (cand1>=ndomain) {
          cand1nlink=scbead[cand1-ndomain]->nlink;
        } else {
          cand1nlink=ringbead[cand1]->nlink;
        }
        if (cand2>=ndomain) {
          cand2nlink=scbead[cand2-ndomain]->nlink;
        } else {
          cand2nlink=ringbead[cand2]->nlink;
        }
        if (grade1==grade2) {
          if (grade1==0) {
            if (cand1>cand2) {
              nextinx=cand1;
            } else {
              nextinx=cand2;
            }
          } else {
            if (cand1nlink > cand2nlink) {
              nextinx=cand1;
            } else {
              nextinx=cand2;
            }
          }
        } else {
          if (grade1>grade2) {
            nextinx=cand2;
          } else {
            nextinx=cand1;
          }
        }
      } else if (whosnext->nlink==2) {
        int indices[4];
        indices[0]=whosnext->bond[0]->b1->inx;
        indices[1]=whosnext->bond[0]->b2->inx;
        indices[2]=whosnext->bond[1]->b1->inx;
        indices[3]=whosnext->bond[1]->b2->inx;
        for (int i=0; i<4; i++) {
          if (indices[i]!=whosnext->inx && indices[i]!=prev->inx) {
            nextinx=indices[i];
          }
        }
      } else {
        int indices[2];
        indices[0]=whosnext->bond[0]->b1->inx;
        indices[1]=whosnext->bond[0]->b2->inx;
        for (int i=0; i<2; i++) {
          if (indices[i]!=whosnext->inx) {
            nextinx=indices[i];
          } else  {
            alltheway.push_back(indices[i]);
          }
        }
      }
      prev=whosnext;
      alltheway.push_back(nextinx);
      if (nextinx>=ndomain) {
        whosnext=scbead[nextinx-ndomain];
      } else {
        whosnext=ringbead[nextinx];
      }
      for (int i=0; i<ndomain+nscbond; i++) {
        Bond *bb=bondsall[i];
        if ( (bb->b1->inx==prev->inx && bb->b2->inx==whosnext->inx) || (bb->b1->inx==whosnext->inx && bb->b2->inx==prev->inx)) {
          nextbond=bb;
        }
      }
      if (whosnext->inx==rbplus->inx) {
        connectedbonds[connectioninx]=new Bondconnection(prevbond,nextbond);
        if (prevbond!=nextbond) {
          annealbackward(connectedbonds[connectioninx]);
        } else {
          annealcap(connectedbonds[connectioninx]);
          prevbond->annealforwarddone=1;
        }
        nextbond->annealbackwarddone=1;
        connectioninx++;
        if (whosnext->inx==0) { 
          connectedbonds[connectioninx]=new Bondconnection(nextbond,rbplus->bond[0]);
          annealRingbackward(connectedbonds[connectioninx]);
          connectioninx++;
        } else {
          connectedbonds[connectioninx]=new Bondconnection(nextbond,rbplus->bond[1]);
          annealRingbackward(connectedbonds[connectioninx]);
          connectioninx++;
        }
        last=1;
      } else {
        if (prevbond!=nextbond) { 
          connectedbonds[connectioninx]=new Bondconnection(prevbond,nextbond);
          if (nextbond->annealforwarddone==0) {
             if (prevbond->annealbackwarddone==1) {
               annealbackwardforward(connectedbonds[connectioninx]);
             } else { 
               annealforwardforward(connectedbonds[connectioninx]);
             }
            nextbond->annealforwarddone=1;
          } else if (nextbond->annealforwarddone==1) {
             annealbackward(connectedbonds[connectioninx]);
             nextbond->annealbackwarddone=1;
          }
          connectioninx++;
        } else { //annealcap
          connectedbonds[connectioninx]=new Bondconnection(prevbond,nextbond);
          annealcap(connectedbonds[connectioninx]);
          prevbond->annealforwarddone=1;
          nextbond->annealbackwarddone=1;
          connectioninx++;
        }
      }
      prevbond=nextbond;
    } while (last==0);
    alltheway.clear();
  }
  int prevmt;
  for (int i=0; i<2*nscbond+ndomain; i++) {
    if (i==0) {
      prevmt=connectedbonds[2*nscbond+ndomain-1]->tofitinx;
    } else {
      prevmt=connectedbonds[i-1]->tofitinx; 
    }
    enterBond(connectedbonds[i]->b1,connectedbonds[i]->b2,connectedbonds[i]->fromstrand,connectedbonds[i]->tostrand,connectedbonds[i]->fromdirection,connectedbonds[i]->todirection,connectedbonds[i]->fromfitinx,connectedbonds[i]->tofitinx,prevmt,natom);
  } 
//  fprintf(stderr,"natom %d\n",natom-1);
//  exit(1);


  for (int i=0; i<ndomain+nscbond; i++) {
    Bond *bb=bondsall[i];
    bb->assignedstrand=-1;
    bb->annealforwarddone=0;
    bb->annealbackwarddone=0;
  }

// fprintf(stderr,"optimizing ...\n");

  long maxlookup=(long)totalfits*(long)totalfits/(long)1000;
  int *lookup=new int[maxlookup];
  int nlookup=0;
  for (int i=0; i<nscbond-1; i++) {
//    fprintf(stderr,"checking distances for bond %d/%d\n",i,nscbond);
    Bond *bi=scbond[i];
    for (int j=i+1; j<nscbond; j++) {
      Bond *bj=scbond[j];
      for (int ii=0; ii<bi->nfit; ii++) {
        for (int ij=0; ij<bj->nfit; ij++) {
          double dx=bi->xfit[ii]-bj->xfit[ij];
          double dy=bi->yfit[ii]-bj->yfit[ij];
          double dz=bi->zfit[ii]-bj->zfit[ij];
          double d=dx*dx+dy*dy+dz*dz;
          if (d<10*10) {
            if (nlookup+4>maxlookup) {
              fprintf(stderr,"reached maxlookup\n");
              exit(1);
            }
            lookup[nlookup++]=i;
            lookup[nlookup++]=ii;
            lookup[nlookup++]=j;
            lookup[nlookup++]=ij;
          }
        }
      }
    }
  }

//  fprintf(stderr,"lookup table has %d entries\n",nlookup);
  int n=0;
  double lastenergy=getEnergy(nlookup,lookup,scbond);
//  fprintf(stderr,"initial energy: %lf\n",lastenergy);
  do {
    int ib=random()%nscbond;
    double delta=(double((random()%9999999))/9999999.0-0.5)*1.0-0.5;
    scbond[ib]->phase+=delta;
    fitDNA(scbond[ib]);
    double newenergy=getEnergy(nlookup,lookup,scbond);
//    fprintf(stderr,"new energy: %lf\n",newenergy);
    double rnum=double((random()%9999999))/9999999.0;
    if (newenergy<lastenergy  || exp(-(newenergy-lastenergy)/kT)>=rnum) {
      lastenergy=newenergy;
//      fprintf(stderr,"%d energy: %lf\n",n,lastenergy);
    } else {
      scbond[ib]->phase-=delta;
      fitDNA(scbond[ib]);
    }
  } while (++n<maxtrials) ;

//  Bond **bondsall=new Bond*[nscbond+ndomain];
//  Bondconnection **connectedbonds=new Bondconnection*[2*nscbond+ndomain];
  connectioninx=0;
// for (int j=0; j<ndomain; j++) {
//    bondsall[j]=ringbond[j];
//  }
//  for (int j=0; j<nscbond; j++) {
//    bondsall[j+ndomain]=scbond[j];
//  }
  for (int i=0; i<ndomain; i++) {
    Bead *rb=ringbead[i];
    Bead *rbplus;
    if (i==ndomain-1) {
      rbplus=ringbead[0];
    } else {
      rbplus=ringbead[i+1];
    }
    int scinx=rbplus->bond[2]->b2->inx;
    Bead *b=scbead[scinx-ndomain];
    int last=0;
    alltheway.push_back(rb->inx);
    alltheway.push_back(rbplus->inx);
    alltheway.push_back(b->inx);
    if (i==0) {
      connectedbonds[connectioninx]=new Bondconnection(rb->bond[0],rbplus->bond[2]);
      annealRingforward(connectedbonds[connectioninx]);
      connectioninx++;
    } else {
      connectedbonds[connectioninx]=new Bondconnection(rb->bond[1],rbplus->bond[2]);
      annealRingforward(connectedbonds[connectioninx]);
      connectioninx++;      
    }
    Bead *whosnext=b;
    Bead *prev=rbplus;
    Bond *prevbond=rbplus->bond[2];
    prevbond->annealforwarddone=1;
    Bond *nextbond;
    do {
      int nextinx;
      if (whosnext->nlink==3) {
        int indices[6];
        indices[0]=whosnext->bond[0]->b1->inx;
        indices[1]=whosnext->bond[0]->b2->inx;
        indices[2]=whosnext->bond[1]->b1->inx;
        indices[3]=whosnext->bond[1]->b2->inx;
        indices[4]=whosnext->bond[2]->b1->inx;
        indices[5]=whosnext->bond[2]->b2->inx;
        int cand1;
        int cand2;
        int cand1nlink;
        int cand2nlink;
        int grade1;
        int grade2;
        int candcount=0;
        for (int i=0; i<6; i++) {
          if (indices[i]!=whosnext->inx && indices[i]!=prev->inx && candcount==0) {
            cand1=indices[i];
            grade1 = std::count( alltheway.begin(), alltheway.end(), cand1);
            candcount++;
          } else if (indices[i]!=whosnext->inx && indices[i]!=prev->inx && candcount==1) {
            cand2=indices[i];
            grade2 = std::count( alltheway.begin(), alltheway.end(), cand2);
          }
        }
        if (cand1>=ndomain) {
          cand1nlink=scbead[cand1-ndomain]->nlink;
        } else {
          cand1nlink=ringbead[cand1]->nlink;
        }
        if (cand2>=ndomain) {
          cand2nlink=scbead[cand2-ndomain]->nlink;
        } else {
          cand2nlink=ringbead[cand2]->nlink;
        }
        if (grade1==grade2) {
          if (grade1==0) {
            if (cand1>cand2) {
              nextinx=cand1;
            } else {
              nextinx=cand2;
            }
          } else {
            if (cand1nlink > cand2nlink) {
              nextinx=cand1;
            } else {
              nextinx=cand2;
            }
          }
        } else {
          if (grade1>grade2) {
            nextinx=cand2;
          } else {
            nextinx=cand1;
          }
        }
      } else if (whosnext->nlink==2) {
        int indices[4];
        indices[0]=whosnext->bond[0]->b1->inx;
        indices[1]=whosnext->bond[0]->b2->inx;
        indices[2]=whosnext->bond[1]->b1->inx;
        indices[3]=whosnext->bond[1]->b2->inx;
        for (int i=0; i<4; i++) {
          if (indices[i]!=whosnext->inx && indices[i]!=prev->inx) {
            nextinx=indices[i];
          }
        }
      } else {
        int indices[2];
        indices[0]=whosnext->bond[0]->b1->inx;
        indices[1]=whosnext->bond[0]->b2->inx;
        for (int i=0; i<2; i++) {
          if (indices[i]!=whosnext->inx) {
            nextinx=indices[i];
          } else  {
            alltheway.push_back(indices[i]);
          }
        }
      }
      prev=whosnext;
      alltheway.push_back(nextinx);
      if (nextinx>=ndomain) {
        whosnext=scbead[nextinx-ndomain];
      } else {
        whosnext=ringbead[nextinx];
      }
      for (int i=0; i<ndomain+nscbond; i++) {
        Bond *bb=bondsall[i];
        if ( (bb->b1->inx==prev->inx && bb->b2->inx==whosnext->inx) || (bb->b1->inx==whosnext->inx && bb->b2->inx==prev->inx)) {
          nextbond=bb;
        }
      }
      if (whosnext->inx==rbplus->inx) {
        connectedbonds[connectioninx]=new Bondconnection(prevbond,nextbond);
        if (prevbond!=nextbond) {
          annealbackward(connectedbonds[connectioninx]);
        } else {
          annealcap(connectedbonds[connectioninx]);
          prevbond->annealforwarddone=1;
        }
        nextbond->annealbackwarddone=1;
        connectioninx++;
        if (whosnext->inx==0) { 
          connectedbonds[connectioninx]=new Bondconnection(nextbond,rbplus->bond[0]);
          annealRingbackward(connectedbonds[connectioninx]);
          connectioninx++;
        } else {
          connectedbonds[connectioninx]=new Bondconnection(nextbond,rbplus->bond[1]);
          annealRingbackward(connectedbonds[connectioninx]);
          connectioninx++;
        }
        last=1;
      } else {
        if (prevbond!=nextbond) { 
          connectedbonds[connectioninx]=new Bondconnection(prevbond,nextbond);
          if (nextbond->annealforwarddone==0) {
            if (prevbond->annealbackwarddone==1) {
              annealbackwardforward(connectedbonds[connectioninx]);
            } else { 
              annealforwardforward(connectedbonds[connectioninx]);
            }
            nextbond->annealforwarddone=1;
          } else if (nextbond->annealforwarddone==1) {
             annealbackward(connectedbonds[connectioninx]);
             nextbond->annealbackwarddone=1;
          }
          connectioninx++;
        } else { //annealcap
          connectedbonds[connectioninx]=new Bondconnection(prevbond,nextbond);
          annealcap(connectedbonds[connectioninx]);
          prevbond->annealforwarddone=1;
          nextbond->annealbackwarddone=1;
          connectioninx++;
        }
      }
      prevbond=nextbond;
    } while (last==0);
    alltheway.clear();
  }
  final=1; 
  natom=1;
  for (int i=0; i<2*nscbond+ndomain; i++) {
    if (i==0) {
      prevmt=connectedbonds[2*nscbond+ndomain-1]->tofitinx;
    } else {
      prevmt=connectedbonds[i-1]->tofitinx; 
    }
    enterBond(connectedbonds[i]->b1,connectedbonds[i]->b2,connectedbonds[i]->fromstrand,connectedbonds[i]->tostrand,connectedbonds[i]->fromdirection,connectedbonds[i]->todirection,connectedbonds[i]->fromfitinx,connectedbonds[i]->tofitinx,prevmt,natom);
  }

  for (int i=0; i<ndomain+nscbond; i++) {
    Bond *bb=bondsall[i];
    bb->assignedstrand=-1;
    bb->annealforwarddone=0;
    bb->annealbackwarddone=0;
  }

  delete scbead; 
  delete lookup;
}
