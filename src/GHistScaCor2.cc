#include "GHistScaCor2.h"

#include <iostream>


GHistScaCor2::GHistScaCor2(const Bool_t linkHistogram) :
    GHistScaCor(linkHistogram)
{
}

GHistScaCor2::GHistScaCor2() :
    GHistScaCor(kTRUE)
{
    buffer                  = new TH2D();
    accumulated             = new TH2D();
    accumulatedCorrected    = new TH2D();

    buffer->SetDirectory(0);
    accumulated->SetDirectory(0);
    accumulatedCorrected->SetDirectory(0);
}

GHistScaCor2::GHistScaCor2(const char* name, const char* title, const Int_t nbinsx, const Double_t xlow, const Double_t xup, const Int_t nbinsy, const Double_t ylow, const Double_t yup, const Bool_t linkHistogram) :
    GHistScaCor(linkHistogram)
{
    buffer                  = new TH2D(TString(name).Append(GHSC_bufferNameSuffix).Data(),
                                       TString(title).Append(GHSC_bufferTitleSuffix).Data(),
                                       nbinsx, xlow, xup, nbinsy, ylow, yup);
    accumulated             = new TH2D(TString(name).Append(GHSC_accumulatedNameSuffix).Data(),
                                       TString(title).Append(GHSC_accumulatedTitleSuffix).Data(),
                                       nbinsx, xlow, xup, nbinsy, ylow, yup);
    accumulatedCorrected    = new TH2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);

    buffer->SetDirectory(0);
    accumulated->SetDirectory(0);
    accumulatedCorrected->SetDirectory(0);
}

GHistScaCor2::~GHistScaCor2()
{
}

void    GHistScaCor2::CreateSingleScalerRead()
{
    TH2D*    uncor   = new TH2D(TString(accumulatedCorrected->GetName()).Append(GHSC_singleScalerReadNameSuffix).Append(TString::Itoa(singleScalerReads.GetEntriesFast(), 10)),
                                TString(accumulatedCorrected->GetTitle()).Append(GHSC_singleScalerReadTitleSuffix).Append(TString::Itoa(singleScalerReads.GetEntriesFast(), 10)),
                                accumulatedCorrected->GetNbinsX(),
                                accumulatedCorrected->GetXaxis()->GetXmin(),
                                accumulatedCorrected->GetXaxis()->GetXmax(),
                                accumulatedCorrected->GetNbinsY(),
                                accumulatedCorrected->GetYaxis()->GetXmin(),
                                accumulatedCorrected->GetYaxis()->GetXmax());
    uncor->SetDirectory(0);
    singleScalerReads.Add(uncor);

    TH2D*    cor     = new TH2D(TString(accumulatedCorrected->GetName()).Append(GHSC_singleScalerReadCorrectedNameSuffix).Append(TString::Itoa(singleScalerReadsCorrected.GetEntriesFast(), 10)),
                                TString(accumulatedCorrected->GetTitle()).Append(GHSC_singleScalerReadCorrectedTitleSuffix).Append(TString::Itoa(singleScalerReadsCorrected.GetEntriesFast(), 10)),
                                accumulatedCorrected->GetNbinsX(),
                                accumulatedCorrected->GetXaxis()->GetXmin(),
                                accumulatedCorrected->GetXaxis()->GetXmax(),
                                accumulatedCorrected->GetNbinsY(),
                                accumulatedCorrected->GetYaxis()->GetXmin(),
                                accumulatedCorrected->GetYaxis()->GetXmax());
    cor->SetDirectory(0);
    singleScalerReadsCorrected.Add(cor);
}

Int_t	GHistScaCor2::Fill(Double_t x)
{
    std::cout << "ERROR: You tried to fill a 2 dim. Hist with only 1 value." << std::endl;
    return 0;
}

void	GHistScaCor2::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax)
{
    buffer->SetBins(nx, xmin, xmax, ny, ymin, ymax);
    accumulated->SetBins(nx, xmin, xmax, ny, ymin, ymax);
    accumulatedCorrected->SetBins(nx, xmin, xmax, ny, ymin, ymax);
    for(Int_t i=0; i<singleScalerReads.GetEntriesFast(); i++)
    {
        ((TH1D*)singleScalerReads.At(i))->SetBins(nx, xmin, xmax, ny, ymin, ymax);
        ((TH1D*)singleScalerReadsCorrected.At(i))->SetBins(nx, xmin, xmax, ny, ymin, ymax);
    }
}


GHistScaCor*    GHistScaCor2::ProjectionX(const char* name, Int_t firstybin, Int_t lastybin, Option_t* option)
{
    GHistScaCor*    ret = new GHistScaCor(name, name, buffer->GetNbinsX(), buffer->GetXaxis()->GetXmin(), buffer->GetXaxis()->GetXmax(), kFALSE);
    TH1D*   help1   = ((TH2D*)buffer)->ProjectionX("_pbuf", firstybin, lastybin, option);
    TH1D*   help2   = ((TH2D*)accumulated)->ProjectionX("_pacc", firstybin, lastybin, option);
    TH1D*   help3   = ((TH2D*)accumulatedCorrected)->ProjectionX("_pacccor", firstybin, lastybin, option);
    ret->Add(help1, help2, help3, IsCorrected());
    if(help1)   delete help1;
    if(help2)   delete help2;
    if(help3)   delete help3;
    return ret;
}

GHistScaCor*   GHistScaCor2::ProjectionY(const char* name, Int_t firstxbin, Int_t lastxbin, Option_t* option)
{
    GHistScaCor*    ret = new GHistScaCor(name, name, buffer->GetNbinsY(), buffer->GetYaxis()->GetXmin(), buffer->GetYaxis()->GetXmax(), kFALSE);
    TH1D*   help1   = ((TH2D*)buffer)->ProjectionY("_pbuf", firstxbin, lastxbin, option);
    TH1D*   help2   = ((TH2D*)accumulated)->ProjectionY("_pacc", firstxbin, lastxbin, option);
    TH1D*   help3   = ((TH2D*)accumulatedCorrected)->ProjectionY("_pacccor", firstxbin, lastxbin, option);
    ret->Add(help1, help2, help3, IsCorrected());
    if(help1)   delete help1;
    if(help2)   delete help2;
    if(help3)   delete help3;
    return ret;
}
