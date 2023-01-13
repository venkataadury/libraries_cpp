#include "maths/maths.h"

double fHD=3,fHr=3,fh=0,PROPK=0.5;
double K=1e9;
Matrix propM(const Vector& pop)
{
  double tP=sum(pop);
  double pHD=PROPK*pop[0][0]/tP,ph=PROPK*pop[1][0]/tP,pHr=PROPK*pop[2][0]/tP;
  Matrix r(3,3);
  r[0][0]=fHD*(pHD+ph/2); r[0][1]=fh*(ph/4+pHD/2); r[0][2]=0;
  r[1][0]=fHD*(ph/2+pHr); r[1][1]=fh/2; r[1][2]=fHr*(pHD+ph/2);
  r[2][0]=0; r[2][1]=fh*(ph/4+pHr/2); r[2][2]=fHr*(ph/2+pHr);
  return r;
}
int main()
{
  Vector pop(250,500,250); //Hardy weinberg Equilibrium
  int nmoves=250;
  double tP=1000;
  for(int i=0;i<nmoves;i++)
  {
    cout <<i<<" "<< pop[0][0]<<" "<<pop[1][0]<<" "<<pop[2][0]<<" "<<tP<<"\t"<<(pop[0][0]/tP)<<" "<<(pop[1][0]/tP)<<" "<<(pop[2][0]/tP)<<"\n";
    pop=pop+(propM(pop)*pop);
    pop=(1-tP/K)*pop;
    tP=sum(pop);
  }
}
