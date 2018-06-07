// Adjusted version of code provided by Mikhail Bashkanov

#include "TMinuit.h"

Double_t z[359],E[359],C[359],errorz[359];

//______________________________________________________________________________
Double_t func(double E,double C,Double_t *par)
{

Double_t P0=par[0]*exp(-1*(E-420)*(E-420)/2/100/100)+par[1]*exp(-1*(E-520)*(E-520)/2/100/100)+par[2]*exp(-1*(E-620)*(E-620)/2/100/100)+par[3]*70*70/4/((E-570)*(E-570)+70*70/4);
Double_t P1=par[4]*exp(-1*(E-420)*(E-420)/2/100/100)+par[5]*exp(-1*(E-520)*(E-520)/2/100/100)+par[6]*exp(-1*(E-620)*(E-620)/2/100/100)+par[7]*70*70/4/((E-570)*(E-570)+70*70/4);
Double_t P2=par[8]*exp(-1*(E-420)*(E-420)/2/100/100)+par[9]*exp(-1*(E-520)*(E-520)/2/100/100)+par[10]*exp(-1*(E-620)*(E-620)/2/100/100)+par[11]*70*70/4/((E-570)*(E-570)+70*70/4);
Double_t P3=par[12]*exp(-1*(E-420)*(E-420)/2/100/100)+par[13]*exp(-1*(E-520)*(E-520)/2/100/100)+par[14]*exp(-1*(E-620)*(E-620)/2/100/100)+par[15]*70*70/4/((E-570)*(E-570)+70*70/4);
Double_t P4=par[16]*exp(-1*(E-420)*(E-420)/2/100/100)+par[17]*exp(-1*(E-520)*(E-520)/2/100/100)+par[18]*exp(-1*(E-620)*(E-620)/2/100/100)+par[19]*70*70/4/((E-570)*(E-570)+70*70/4);
Double_t P5=par[20]*exp(-1*(E-420)*(E-420)/2/100/100)+par[21]*exp(-1*(E-520)*(E-520)/2/100/100)+par[22]*exp(-1*(E-620)*(E-620)/2/100/100)+par[23]*70*70/4/((E-570)*(E-570)+70*70/4);
Double_t P6=par[24]*exp(-1*(E-420)*(E-420)/2/100/100)+par[25]*exp(-1*(E-520)*(E-520)/2/100/100)+par[26]*exp(-1*(E-620)*(E-620)/2/100/100)+par[27]*70*70/4/((E-570)*(E-570)+70*70/4);
Double_t P7=par[28]*exp(-1*(E-420)*(E-420)/2/100/100)+par[29]*exp(-1*(E-520)*(E-520)/2/100/100)+par[30]*exp(-1*(E-620)*(E-620)/2/100/100)+par[31]*70*70/4/((E-570)*(E-570)+70*70/4);

 Double_t value=(1-C*C)*(P0*3+P1*15*C+P2*15.0/2*(7*C*C-1)+P3*105.0/2*C*(3*C*C-1)+P4*105.0/8*(33*C*C*C*C-18*C*C+1)+P5*63.0/8*C*(143*C*C*C*C-110*C*C+15)+P6*315.0/16*(143*C*C*C*C*C*C-143*C*C*C*C+33*C*C-1)+P7*495.0/16*(221*C*C*C*C*C*C*C-273*C*C*C*C*C+91*C*C*C-7*C));

 return value;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   const Int_t nbins = 359;
   Int_t i;

    //calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (i=0;i<nbins; i++) {
     delta  = (z[i]-func(E[i],C[i],par))/errorz[i];
     chisq += delta*delta;
   }
   f = chisq;
}

//______________________________________________________________________________
void Fit_SigmaEDep_3GausBW()
{
    Int_t inum=0;
    Double_t Et,Ct,Zt,Zet;
    Int_t itt;
    ifstream infile1("/scratch/Mainz_Software/a2GoAT/Sigma_NS18.dat");//359 points
    for(Int_t i=0;i<359;i++){
        infile1 >> itt >> Et >> Ct >> Zt >> Zet;
        E[inum]=Et;
        C[inum]=Ct;
        z[inum]=Zt;
        errorz[inum]=Zet;
        inum+=1;
    }

    TMinuit *gMinuit = new TMinuit(32);  //initialize TMinuit with a maximum of 5 params
    gMinuit->SetFCN(fcn);

    Double_t arglist[64];
    Int_t ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

   static Double_t vstart[32] = {-0.108754,-0.00356053,0.0745446,-0.0212027,
                                0.0238134,-0.00919569,0.0305643,-0.0102451,
                                0.00121507,-0.00177259,0.0044306,-0.00475281,
                                0.00533593,-0.000667558,0.00358099,-0.00188308,
                                0.0029606,-0.00206478,0.00386589,-0.00364641,
                                0.00146578,-0.000521466,0.0013484,-0.00115159,
                                0,0,0,0,};

      static Double_t step[32] = {0.0051 , 0.0051 , 0.0051 , 0.0051,
			       0.0051 , 0.0051 , 0.0051 , 0.0051,
			       0.0051 , 0.0051 , 0.0051 , 0.0051,
			       0.0051 , 0.0051 , 0.0051 , 0.0051,
			       0.0051 , 0.0051 , 0.0051 , 0.0051,
			       0.0051 , 0.0051 , 0.0051 , 0.0051,
			       0.0051 , 0.0051 , 0.0051 , 0.0051,
			       0.0051 , 0.0051 , 0.0051 , 0.0051

   };

    for(Int_t i=0;i<32;i++){
        gMinuit->mnparm(i, Form("a%i",i), vstart[i], step[i], 0,0,ierflg);
    }

    gMinuit->FixParameter(24);
    gMinuit->FixParameter(25);
    gMinuit->FixParameter(26);
    gMinuit->FixParameter(27);
    gMinuit->FixParameter(28);
    gMinuit->FixParameter(29);
    gMinuit->FixParameter(30);
    gMinuit->FixParameter(31);

    // Now ready for minimization step
    arglist[0] = 15000;
    arglist[1] = 1.;
    gMinuit->SetMaxIterations(15000);

    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    Double_t ParF[32],eParF[32];

    Double_t PPM[3];
    for(Int_t i=0;i<32;i++){gMinuit->GetParameter(i,ParF[i],eParF[i]);}

    cout<<"Values"<<endl;
    for(Int_t i=0;i<7;i++){cout<<ParF[0+i*4]<<","<<ParF[1+i*4]<<","<<ParF[2+i*4]<<","<<ParF[3+i*4]<<","<<endl;}
    cout<<"Errors"<<endl;
    for(Int_t i=0;i<7;i++){cout<<eParF[0+i*4]<<","<<eParF[1+i*4]<<","<<eParF[2+i*4]<<","<<eParF[3+i*4]<<","<<endl;}

}
