// Adjusted version of code provided by Mikhail Bashkanov

#include "TMinuit.h"

Double_t z[359],E[359],C[359],errorz[359];

//______________________________________________________________________________
Double_t func(double E,double C,Double_t *par)
{
    Double_t P0=par[0]*exp(-1*(E-420)*(E-420)/2/60/60)+par[1]*exp(-1*(E-450)*(E-450)/2/60/60)+par[2]*exp(-1*(E-480)*(E-480)/2/60/60)+par[3]*exp(-1*(E-510)*(E-510)/2/60/60)+par[4]*exp(-1*(E-540)*(E-540)/2/60/60)+par[5]*exp(-1*(E-570)*(E-570)/2/60/60)+par[6]*exp(-1*(E-600)*(E-600)/2/60/60);
    Double_t P1=par[7]*exp(-1*(E-420)*(E-420)/2/60/60)+par[8]*exp(-1*(E-450)*(E-450)/2/60/60)+par[9]*exp(-1*(E-480)*(E-480)/2/60/60)+par[10]*exp(-1*(E-510)*(E-510)/2/60/60)+par[11]*exp(-1*(E-540)*(E-540)/2/60/60)+par[12]*exp(-1*(E-570)*(E-570)/2/60/60)+par[13]*exp(-1*(E-600)*(E-600)/2/60/60);
    Double_t P2=par[14]*exp(-1*(E-420)*(E-420)/2/60/60)+par[15]*exp(-1*(E-450)*(E-450)/2/60/60)+par[16]*exp(-1*(E-480)*(E-480)/2/60/60)+par[17]*exp(-1*(E-510)*(E-510)/2/60/60)+par[18]*exp(-1*(E-540)*(E-540)/2/60/60)+par[19]*exp(-1*(E-570)*(E-570)/2/60/60)+par[20]*exp(-1*(E-600)*(E-600)/2/60/60);
    Double_t P3=par[21]*exp(-1*(E-420)*(E-420)/2/60/60)+par[22]*exp(-1*(E-450)*(E-450)/2/60/60)+par[23]*exp(-1*(E-480)*(E-480)/2/60/60)+par[24]*exp(-1*(E-510)*(E-510)/2/60/60)+par[25]*exp(-1*(E-540)*(E-540)/2/60/60)+par[26]*exp(-1*(E-570)*(E-570)/2/60/60)+par[27]*exp(-1*(E-600)*(E-600)/2/60/60);
    Double_t P4=par[28]*exp(-1*(E-420)*(E-420)/2/60/60)+par[29]*exp(-1*(E-450)*(E-450)/2/60/60)+par[30]*exp(-1*(E-480)*(E-480)/2/60/60)+par[31]*exp(-1*(E-510)*(E-510)/2/60/60)+par[32]*exp(-1*(E-540)*(E-540)/2/60/60)+par[33]*exp(-1*(E-570)*(E-570)/2/60/60)+par[34]*exp(-1*(E-600)*(E-600)/2/60/60);
    Double_t P5=par[35]*exp(-1*(E-420)*(E-420)/2/60/60)+par[36]*exp(-1*(E-450)*(E-450)/2/60/60)+par[37]*exp(-1*(E-480)*(E-480)/2/60/60)+par[38]*exp(-1*(E-510)*(E-510)/2/60/60)+par[39]*exp(-1*(E-540)*(E-540)/2/60/60)+par[40]*exp(-1*(E-570)*(E-570)/2/60/60)+par[41]*exp(-1*(E-600)*(E-600)/2/60/60);
    Double_t P6=par[42]*exp(-1*(E-420)*(E-420)/2/60/60)+par[43]*exp(-1*(E-450)*(E-450)/2/60/60)+par[44]*exp(-1*(E-480)*(E-480)/2/60/60)+par[45]*exp(-1*(E-510)*(E-510)/2/60/60)+par[46]*exp(-1*(E-540)*(E-540)/2/60/60)+par[47]*exp(-1*(E-570)*(E-570)/2/60/60)+par[48]*exp(-1*(E-600)*(E-600)/2/60/60);
    Double_t P7=par[49]*exp(-1*(E-420)*(E-420)/2/60/60)+par[50]*exp(-1*(E-450)*(E-450)/2/60/60)+par[51]*exp(-1*(E-480)*(E-480)/2/60/60)+par[52]*exp(-1*(E-510)*(E-510)/2/60/60)+par[53]*exp(-1*(E-540)*(E-540)/2/60/60)+par[54]*exp(-1*(E-570)*(E-570)/2/60/60)+par[55]*exp(-1*(E-600)*(E-600)/2/60/60);

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
void Fit_SigmaEDep_7Gaus()
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

    TMinuit *gMinuit = new TMinuit(56);  //initialize TMinuit with a maximum of 5 params
    gMinuit->SetFCN(fcn);

    Double_t arglist[56];
    Int_t ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    static Double_t vstart[56] = {2.99203,-11.6753,22.4077,-27.858,23.5892,-12.9431,3.63412,
    -1.67285,5.72189,-9.9412,11.3262,-8.83277,4.43907,-1.10018,
    0.648495,-1.96551,3.04773,-3.14691,2.31897,-1.16499,0.309782,
    -0.884155,3.09574,-5.41028,6.06124,-4.58708,2.23079,-0.544889,
    -0.134468,0.560091,-1.0877,1.2847,-0.972672,0.445549,-0.0954476,
    -0.460013,1.41646,-2.17844,2.15514,-1.45138,0.633165,-0.139968,
    0,0,0,0,0,0,0,
    };

    static Double_t step[56] = {0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,
    0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,
    0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,
    0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,
    0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,
    0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,
    0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,
    0.51 , 0.51 , 0.51 , 0.51, 0.51, 0.51, 0.51,

    };

    gMinuit->mnparm(0, "a0", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "a1", vstart[1], step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "a2", vstart[2], step[2], 0,0,ierflg);
    gMinuit->mnparm(3, "a3", vstart[3], step[3], 0,0,ierflg);
    gMinuit->mnparm(4, "a4", vstart[4], step[4], 0,0,ierflg);
    gMinuit->mnparm(5, "a5", vstart[5], step[5], 0,0,ierflg);
    gMinuit->mnparm(6, "a6", vstart[6], step[6], 0,0,ierflg);

    gMinuit->mnparm(7, "a7", vstart[7], step[7], 0,0,ierflg);
    gMinuit->mnparm(8, "a8", vstart[8], step[8], 0,0,ierflg);
    gMinuit->mnparm(9, "a9", vstart[9], step[9], 0,0,ierflg);
    gMinuit->mnparm(10, "a10", vstart[10], step[10], 0,0,ierflg);
    gMinuit->mnparm(11, "a11", vstart[11], step[11], 0,0,ierflg);
    gMinuit->mnparm(12, "a12", vstart[12], step[12], 0,0,ierflg);
    gMinuit->mnparm(13, "a13", vstart[13], step[13], 0,0,ierflg);

    gMinuit->mnparm(14, "a14", vstart[14], step[14], 0,0,ierflg);
    gMinuit->mnparm(15, "a15", vstart[15], step[15], 0,0,ierflg);
    gMinuit->mnparm(16, "a16", vstart[16], step[16], 0,0,ierflg);
    gMinuit->mnparm(17, "a17", vstart[17], step[17], 0,0,ierflg);
    gMinuit->mnparm(18, "a18", vstart[18], step[18], 0,0,ierflg);
    gMinuit->mnparm(19, "a19", vstart[19], step[19], 0,0,ierflg);
    gMinuit->mnparm(20, "a20", vstart[20], step[20], 0,0,ierflg);

    gMinuit->mnparm(21, "a21", vstart[21], step[21], 0,0,ierflg);
    gMinuit->mnparm(22, "a22", vstart[22], step[22], 0,0,ierflg);
    gMinuit->mnparm(23, "a23", vstart[23], step[23], 0,0,ierflg);
    gMinuit->mnparm(24, "a24", vstart[24], step[24], 0,0,ierflg);
    gMinuit->mnparm(25, "a25", vstart[25], step[25], 0,0,ierflg);
    gMinuit->mnparm(26, "a26", vstart[26], step[26], 0,0,ierflg);
    gMinuit->mnparm(27, "a27", vstart[27], step[27], 0,0,ierflg);

    gMinuit->mnparm(28, "a28", vstart[28], step[28], 0,0,ierflg);
    gMinuit->mnparm(29, "a29", vstart[29], step[29], 0,0,ierflg);
    gMinuit->mnparm(30, "a30", vstart[30], step[30], 0,0,ierflg);
    gMinuit->mnparm(31, "a31", vstart[31], step[31], 0,0,ierflg);
    gMinuit->mnparm(32, "a32", vstart[32], step[32], 0,0,ierflg);
    gMinuit->mnparm(33, "a33", vstart[33], step[33], 0,0,ierflg);
    gMinuit->mnparm(34, "a34", vstart[34], step[34], 0,0,ierflg);

    gMinuit->mnparm(35, "a35", vstart[35], step[35], 0,0,ierflg);
    gMinuit->mnparm(36, "a36", vstart[36], step[36], 0,0,ierflg);
    gMinuit->mnparm(37, "a37", vstart[37], step[37], 0,0,ierflg);
    gMinuit->mnparm(38, "a38", vstart[38], step[38], 0,0,ierflg);
    gMinuit->mnparm(39, "a39", vstart[39], step[39], 0,0,ierflg);
    gMinuit->mnparm(40, "a40", vstart[40], step[40], 0,0,ierflg);
    gMinuit->mnparm(41, "a41", vstart[41], step[41], 0,0,ierflg);

    gMinuit->mnparm(42, "a42", vstart[42], step[42], 0,0,ierflg);
    gMinuit->mnparm(43, "a43", vstart[43], step[43], 0,0,ierflg);
    gMinuit->mnparm(44, "a44", vstart[44], step[44], 0,0,ierflg);
    gMinuit->mnparm(45, "a45", vstart[45], step[45], 0,0,ierflg);
    gMinuit->mnparm(46, "a46", vstart[46], step[46], 0,0,ierflg);
    gMinuit->mnparm(47, "a47", vstart[47], step[47], 0,0,ierflg);
    gMinuit->mnparm(48, "a48", vstart[48], step[48], 0,0,ierflg);

    gMinuit->mnparm(49, "a49", vstart[49], step[49], 0,0,ierflg);
    gMinuit->mnparm(50, "a50", vstart[50], step[50], 0,0,ierflg);
    gMinuit->mnparm(51, "a51", vstart[51], step[51], 0,0,ierflg);
    gMinuit->mnparm(52, "a52", vstart[52], step[52], 0,0,ierflg);
    gMinuit->mnparm(53, "a53", vstart[53], step[53], 0,0,ierflg);
    gMinuit->mnparm(54, "a54", vstart[54], step[54], 0,0,ierflg);
    gMinuit->mnparm(55, "a55", vstart[55], step[55], 0,0,ierflg);

    gMinuit->FixParameter(42);
    gMinuit->FixParameter(43);
    gMinuit->FixParameter(44);
    gMinuit->FixParameter(45);
    gMinuit->FixParameter(46);
    gMinuit->FixParameter(47);
    gMinuit->FixParameter(48);

    gMinuit->FixParameter(49);
    gMinuit->FixParameter(50);
    gMinuit->FixParameter(51);
    gMinuit->FixParameter(52);
    gMinuit->FixParameter(53);
    gMinuit->FixParameter(54);
    gMinuit->FixParameter(55);

    // Now ready for minimization step
    arglist[0] = 15000;
    arglist[1] = 1.;
    gMinuit->SetMaxIterations(15000);

    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    //gMinuit->mnprin(3,amin);
    Double_t ParF[56],eParF[56];

    Double_t PPM[3];
    for(Int_t i=0;i<56;i++){gMinuit->GetParameter(i,ParF[i],eParF[i]);}

    cout<<"Values"<<endl;
    for(Int_t i=0;i<7;i++){cout<<ParF[0+i*7]<<","<<ParF[1+i*7]<<","<<ParF[2+i*7]<<","<<ParF[3+i*7]<<","<<ParF[4+i*7]<<","<<ParF[5+i*7]<<","<<ParF[6+i*7]<<","<<endl;}
    cout<<"Errors"<<endl;
    for(Int_t i=0;i<7;i++){cout<<eParF[0+i*7]<<","<<eParF[1+i*7]<<","<<eParF[2+i*7]<<","<<eParF[3+i*7]<<","<<eParF[4+i*7]<<","<<eParF[5+i*7]<<","<<eParF[6+i*7]<<","<<endl;}

}
