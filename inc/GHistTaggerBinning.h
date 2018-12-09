#ifndef __GHistTaggerBinning_h__
#define __GHistTaggerBinning_h__


#include <TROOT.h>
#include <TDirectory.h>

#include "GHistBGSub2.h"


class   GTreeTagger;

class  GHistTaggerBinning  : public GHistLinked
{
private:

protected:
    GHistBGSub* sum;
    GHistBGSub* array;


    GHistTaggerBinning(Bool_t linkHistogram);

public:
    GHistTaggerBinning();
    GHistTaggerBinning(const char* name, const char* title,
                       const Int_t nbinsx, const Double_t xlow, const Double_t xup,
                       const Int_t nTaggerBins = 352,
                       const Bool_t linkHistogram = kTRUE);
    virtual ~GHistTaggerBinning();

    virtual Bool_t          Add(const GHistTaggerBinning* h, Double_t c = 1);
    virtual void            CalcResult();
    const   GHistBGSub*     GetArray()  const   {return array;}
    const   GHistBGSub*     GetSum()    const   {return sum;}
    virtual Int_t           Fill(const Double_t value)                                                          {return sum->Fill(value);}
    virtual Int_t           Fill(const Double_t value, const Double_t taggerTime)                                  {return sum->Fill(value, taggerTime);}
    virtual Int_t           FillWeighted(const Double_t value, const Double_t taggerTime, const Double_t Weight)                                  {return sum->FillWeighted(value, taggerTime, Weight);}
    virtual Int_t           Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel)       {return ((GHistBGSub2*)array)->Fill(value, taggerChannel, taggerTime);}
    virtual Int_t           Fill(const Double_t value, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void            PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void            Reset(Option_t* option = "")                    {sum->Reset(option); array->Reset(option);}
    virtual void        	Scale(Double_t c1 = 1, Option_t* option = "")   {sum->Scale(c1, option); array->Scale(c1, option);}
    virtual void            ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t           WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};




#endif
