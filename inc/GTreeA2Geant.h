#ifndef __GTreeA2Geant_h__
#define __GTreeA2Geant_h__

#include "GTree.h"
#include "Rtypes.h"

#include "TVector3.h"
#include "TLorentzVector.h"

#include <vector>
#include <stdexcept>

#define GEANT_TREE_NAME     "h12"

//keep in syn with A2CBoutput.h in a2geant
#define GEANT_MAX_TAPSHITS  438
#define GEANT_MAX_CBHITS    720
#define GEANT_MAX_PIDHITS    24
#define GEANT_MAX_MWPCHITS  400
#define GEANT_MAX_PART      100

/**
 * @brief A GTree for a2geant output files.
 *
 *   Reads the TTree named "h12". Supplies Geant4 based detector hit inforamtion and also
 *   copies the MC-True inforamtion.
 *
 * Use the "hadd" program from ROOT to copy the a2geant tree into the acqu output file:
 * $ hadd target.root Acqu.root geant.root
 * And then run goat on target.root.
 */
class  GTreeA2Geant : public GTree
{
public:
       typedef std::vector<Double_t> hitvector;

protected:

        // Brach memories
       Int_t           fnhits;
       Int_t           fnpart;
       Int_t           fntaps;
       Int_t           fnvtaps;
       Int_t           fvhits;
       Float_t         plab[GEANT_MAX_PART];
       Float_t         tctaps[GEANT_MAX_TAPSHITS];
       Float_t         fvertex[3];
       Float_t         fbeam[5];
       Float_t         dircos[GEANT_MAX_PART][3];
       Float_t         ecryst[GEANT_MAX_CBHITS];
       Float_t         tcryst[GEANT_MAX_CBHITS];
       Float_t         ectapfs[GEANT_MAX_TAPSHITS];
       Float_t         ectapsl[GEANT_MAX_TAPSHITS];
       Float_t         elab[GEANT_MAX_PART];
       Float_t         feleak;
       Float_t         fenai;
       Float_t         fetot;
       Float_t         eveto[GEANT_MAX_PIDHITS];
       Float_t         tveto[GEANT_MAX_PIDHITS];
       Float_t         evtaps[GEANT_MAX_TAPSHITS];
       Int_t           icryst[GEANT_MAX_CBHITS];
       Int_t           ictaps[GEANT_MAX_TAPSHITS];
       Int_t           ivtaps[GEANT_MAX_TAPSHITS];
       Int_t           idpart[GEANT_MAX_PART];
       Int_t           iveto[GEANT_MAX_PIDHITS];
       Int_t           fnmwpc;
       Int_t           imwpc[GEANT_MAX_MWPCHITS];
       Float_t         mposx[GEANT_MAX_MWPCHITS];
       Float_t         mposy[GEANT_MAX_MWPCHITS];
       Float_t         mposz[GEANT_MAX_MWPCHITS];
       Float_t         emwpc[GEANT_MAX_MWPCHITS];
       Long64_t        mc_evt_id;
       Long64_t        mc_rnd_id;


    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

    static Double_t sumArray( const Float_t* const data, const Int_t size);
    static void buildPattern( hitvector& pattern, const Int_t* const indices, const Float_t* const data, const Int_t nhits, const Int_t patternsize);

public:

    GTreeA2Geant(GTreeManager *Manager);
    virtual ~GTreeA2Geant();

    virtual TLorentzVector GetBeam() const;
    virtual TVector3    GetVertex() const;

    virtual Double_t    GetCBESum() const;
    virtual Double_t    GetPIDESum() const;
    virtual Double_t    GetTAPSESum() const;
    virtual Double_t    GetTAPSVetoESum() const;

    virtual Int_t                   GetNCBHits() const;
    virtual Int_t                   GetCBHitIndex(const UInt_t n) const throw (std::out_of_range);
    virtual const Int_t * const     GetCBHitIndices() const;
    virtual Float_t                 GetCBHitEnergy(const UInt_t n) const throw (std::out_of_range);
    virtual const Float_t * const   GetCBHitEnergies() const;

    virtual Int_t                   GetNPIDHits() const;
    virtual Int_t                   GetPIDHitIndex(const UInt_t n) const throw (std::out_of_range);
    virtual const Int_t * const     GetPIDHitIndices() const;
    virtual Float_t                 GetPIDHitEnergy(const UInt_t n) const throw (std::out_of_range);
    virtual const Float_t * const   GetPIDHitEnergies() const;
    virtual Float_t                 GetPIDHitTime(const UInt_t n) const throw (std::out_of_range);
    virtual const Float_t * const   GetPIDHitTimes() const;

    virtual Int_t                   GetNTAPSHits() const;
    virtual Int_t                   GetTAPSHitIndex(const UInt_t n) const throw (std::out_of_range);
    virtual const Int_t * const     GetTAPSHitIndices() const;
    virtual Float_t                 GetTAPSHitEnergyLong(const UInt_t n) const throw (std::out_of_range);
    virtual const Float_t * const   GetTAPSHitEnergiesLong() const;
    virtual Float_t                 GetTAPSHitEnergyShort(const UInt_t n) const throw (std::out_of_range);
    virtual const Float_t * const   GetTAPSHitEnergiesShort() const;
    virtual Float_t                 GetTAPSHitTime(const UInt_t n) const throw (std::out_of_range);
    virtual const Float_t * const   GetTAPSHitTimes() const;

    virtual Int_t                   GetNTAPSVetoHits() const;
    virtual Int_t                   GetTAPSVetoHitIndex(const UInt_t n) const throw (std::out_of_range);
    virtual const Int_t * const     GetTAPSVetoHitIndices() const;
    virtual Float_t                 GetTAPSVetoHitEnergy(const UInt_t n) const throw (std::out_of_range);
    virtual const Float_t * const   GetTAPSVetoHitEnergies() const;

    virtual Int_t                   GetNMWPCHits() const;
    virtual Int_t                   GetMWPCIndex( const UInt_t n )  const throw (std::out_of_range);
    virtual const Int_t * const     GetMWPCHitIndices() const;
    virtual Float_t                 GetMWPCEnergy( const UInt_t n ) const throw (std::out_of_range);
    virtual const Float_t * const   GetMWPCHitEnergies() const;
    virtual const Float_t * const   GetMWPCHitPosX() const;
    virtual const Float_t * const   GetMWPCHitPosY() const;
    virtual const Float_t * const   GetMWPCHitPosZ() const;
    virtual TVector3                GetMWPCVector( const UInt_t n ) const throw (std::out_of_range);

    virtual void        BuildCBHitPattern( hitvector& pattern) const;
    virtual void        BuildTAPSHitPattern( hitvector& pattern) const;

    virtual UInt_t      GetNTrueParticles() const;
    virtual UInt_t      GetTrueID( const UInt_t n ) const throw (std::out_of_range);
    virtual TLorentzVector GetTrueVector( const UInt_t n ) const throw (std::out_of_range);

    virtual void        Clear();

};


#endif
