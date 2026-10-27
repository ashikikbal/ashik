/* ============================================================================
 * Project: nOOn Monte Carlo Generator - GDR Physics Extension
 * 
 * ORIGINAL nOOn FRAMEWORK:
 * Authors: M. Broz, J. G. Contreras, J. D. Tapia Takaki
 * 
 * GDR PHYSICS MODIFICATIONS:
 * New Implementations: Giant Dipole Resonance (GDR) physics, neutron azimuthal correlations
 *                 and integration to the default nOOn generator
 *                 
 * Authors: Ashik Ikbal Sheikh and Prithwish Tribedy
 * Copyright (c) 2026 Ashik Ikbal Sheikh & Prithwish Tribedy (for new implementations)
 * 
 * LICENSE & CITATION: 
 * The original nOOn framework is subject to its respective license. 
 * Use of the GDR physics modifications requires citation of the 
 * corresponding publications by A. I. Sheikh, P. Tribedy, et al. detailing 
 * this specific implementation.
 * ============================================================================ */

#ifndef NEUTRONGENERATOR_H
#define NEUTRONGENERATOR_H

#include <fstream>

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

class NeutronGenerator : public TObject
{
public:
  enum RunMode_t {
    kMassRapidity,
    kStarlightAscii,
    kFlatMultiplicity,
    kInterface,
    k1n1n
  };
  enum HadronicIntModel_t {
    kGlauber,
    kHardSphere
  };
  enum ProductionMode_t{
    kPhotonPomeron,
    kTwoPhoton
  };
  enum Nucleus_t {
    kPb208,
    kAu197,
    kU238,
    kXe129
  };

  NeutronGenerator();
  virtual ~NeutronGenerator();
  
  void		GenerateEvent(const Double_t photonK1, const Double_t photonK2 = -1);
  TClonesArray* ImportParticles(){return fParticles;}
  
  void Initialize();
  void SetBeamParameters(Nucleus_t nucleus, Double_t gamma);
  void SetBeamParameters(Int_t nuclZ, Int_t nuclA, Double_t gamma);
  void Setup();
  void Run(const UInt_t nEvents);
  
  void SetDataPath(const char* path){fDataPath = path;}
  void SetRunMode(RunMode_t mode,const char *filename = " ", const char *name1 = " ", const char *name2 = " ");
  void SetRapidityCut(Double_t minY, Double_t maxY){fRapMin = minY; fRapMax = maxY;}
  void SetMassCut(Double_t minM, Double_t maxM){fMassMin = minM; fMassMax = maxM;}
  void SetStoreQA(){kStoreQA = kTRUE;}
  void SetStoreGeneratorFunctions(){kStoreGen = kTRUE;}
  void SetHadronicInteractionModel(HadronicIntModel_t model){fHadronicIntModel = model;}
  void SetProductionMode(ProductionMode_t mode){fProductionMode = mode;}
  ProductionMode_t GetProductionMode(){return fProductionMode;}
  
  void LoadENDF(const char *filename = "hENDF.root",const char *histname = "hENDF_2D"); // Edited by Ashik (moved from protected to public)
  void 	FinishEvent();
  void 	FinishProduction();
  
  //–––––––––––––––––––––––––––––––––––––––––––––  
  /// Public helper: photon density d^3n/(dk·d^2b) at (b, photonK).
  /// Internally calls the protected PhotonDensity(b, photonK).
  Double_t GetPhotonDensity(Double_t b, Double_t photonK) { // b in fm and photonK in GeV
    return PhotonDensity(b, photonK);
  }

  void SetImpactParameterPhi(Double_t bAngle) { fImpactParameterPhi = bAngle; } // Added by Ashik   
  void SetNRun(Int_t nRun) { fNRun = nRun; } // Added by Ashik   
  void SetImpactParameterMag(Double_t bMag) { fImpactParameterMag = bMag; } // Added by Ashik   
  void SetEMkickMode(Bool_t emkickOn = kTRUE) { fEMkickIsOn = emkickOn; } // Added by Ashik
  void SetClassicalGDRMode(Bool_t gdrOn = kTRUE) { fGDRIsOn = gdrOn; } // Added by Ashik
  void SetEntangledGDRMode(Bool_t entanglementOn = kTRUE) { fEntanglementIsOn = entanglementOn; } // Added by Ashik

  Double_t GetEMDCrossSection();
  Double_t GetBreakupProbability(const Double_t photonK, const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2);
  Double_t GetTotalFlux(const Double_t photonK1, const Double_t photonK2 = -1);
  Double_t GetBreakupProbability(const Double_t photonK1, const Double_t photonK2, const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2, const Bool_t interpolate = kTRUE);
  
  Double_t *PhotonFlux(const Double_t energyPhoton, const Bool_t printout);
  Double_t HadronicInteractionProbability(const Double_t impactPar); //modified by Prithwish Jun 3,2025
  Double_t SampleGDRTime(double Gamma);  

  Double_t GetProbability(Double_t dphi, Double_t t, Double_t delta_phase, Double_t omega, Double_t Gamma, Double_t alpha_target);
  Double_t GetMaxEnvelope(Double_t t, Double_t Gamma, Double_t alpha_target);
  void GenerateOneNeutronPhiFixedTime(TRandom3* rnd, Double_t phi_gamma, Double_t delta_phase, Double_t t_fixed, Double_t omega, Double_t Gamma, Double_t alpha_target, Double_t &phi_out);

/*
 /// Accessor for hadronic interaction probability (read-only use)
  Double_t EvalHadronicInteractionProbability(Double_t b) const {
    return HadronicInteractionProbability(b);
  }

  void BuildHadronicInteractionProbabilityTable();  // add this line here
*/

// --------------------------------------------------------------------------------
// GenerateEventAtFixedB: use a sampled b and photon k to produce neutrons exclusively
//   at that b (instead of integrating over b).
//
//   b_sampled : the sampled impact parameter (fm)
//   photonK   : the photon energy (GeV) to use in flux/ENDF lookups
//

 inline Double_t FmToInvMeV(Double_t x_fm);

 inline Double_t WrapToPi(Double_t x);

 inline Double_t ComputeSigmaDeltaPhysics(Double_t omega1, Double_t omega2, Double_t Gamma1, Double_t Gamma2);


 const Bool_t use_phase = kTRUE;
// const Bool_t delta_shared_between_nuclei = kTRUE; // this should give 1+a1cos(phi)
 const Bool_t delta_shared_between_nuclei = kFALSE; // this should give 1+a2cos(2phi)

  void GenerateEventAtFixedB(Double_t b_sampled, Double_t photonK);

  /// Return vector {P_0(b), P_1(b),...., P_{maxNeutrons}(b)} by calling the protected helper

  std::vector<Double_t> GetBreakupProbabilities(Double_t b){
    // 1) Call the protected method, which returns a raw array
    Double_t* raw = NucleusBreakupProbability(b);
    // 2) Copy into an std::vector
  //  std::vector<Double_t> v(raw, raw + (maxNeutrons + 1));
    std::vector<Double_t> v(raw, raw + (maxNeutrons)); // edited by Ashik from above line to this one
    // 3) Free the raw array (NucleusBreakupProbability used new[] internally)
    delete[] raw;
    return v;
  }



protected:

  Double_t ComputeXnProbability(const Double_t impactPar);
  void BuildNucleusBreakupProbabilityTable();
  Double_t *NucleusBreakupProbability(const Double_t impactPar);
  Double_t GetBR(const Double_t energyPhoton, const Int_t nNeutrons); 
  void 	CreateNeutrons(const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2);
  void 	CreateGDRNeutrons(const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2); // Added by Ashik 
  void 	CreateGDREntangledNeutrons(const Int_t nNeutronsBeam1, const Int_t nNeutronsBeam2); // Added by Ashik 


  Double_t *ReadXSFile(const char *filename, Double_t scale = 1.0);
  Double_t *MakeReggeParametrization(Double_t startEnergy, Int_t nPoints);
  void InsertDataset(TGraph& graph, Double_t *dataset);
  void InsertPoint(TGraph& graph, Double_t x, Double_t y);
  void InitQAhistograms();
  void ReadENDF(const char *filename, Bool_t saveToFile = kFALSE);
 // void LoadENDF(const char *filename = "hENDF.root",const char *histname = "hENDF_2D");

  
  void BuildPhotonFluxModulationTables();
  void BuildTwoPhotonFluxModulationTables();
 
  Double_t *TwoPhotonFlux(const Double_t energyPhoton1, const Double_t energyPhoton2);
  
  void BuildHadronicInteractionProbabilityTable(); //Modified  by PT June 3,2025
//  Double_t HadronicInteractionProbability(const Double_t impactPar) const; //Modifed by Prithwish June 3,2025
  
  Double_t PhotonDensity(const Double_t distance, const Double_t energyPhoton);
  Double_t PhotonDensity2(const Double_t distance, const Double_t energyPhoton);
  TH1D* CreateHist1D(const char* name, const char* title,Int_t nBins, Double_t xMin, Double_t xMax, const char* xLabel, const char* yLabel);
  TH1D* CreateHist1D(const char* name, const char* title,Int_t nBins, Double_t* xBins, const char* xLabel, const char* yLabel);
  TH2D* CreateHist2D(const char* name, const char* title,Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t yMin, Double_t yMax, const char* xLabel, const char* yLabel, const char* zLabel);
  Int_t FromMatrixToVector(Int_t i, Int_t j);
  void  FromVectorToMatrix(Int_t index, Int_t &row, Int_t &col); 

private:

  NeutronGenerator(const NeutronGenerator &o);
  NeutronGenerator &operator=(const NeutronGenerator &o);
  
  RunMode_t		fRunMode;
  HadronicIntModel_t	fHadronicIntModel;
  ProductionMode_t	fProductionMode;
  Nucleus_t		fNucleus;
  UInt_t 		iEvent;
  TH1			*hInputRapidity;
  TH1			*hInputMass;
  Double_t		fRapMin;
  Double_t		fRapMax;
  Double_t		fMassMin;
  Double_t		fMassMax;
  //ifstream 		fInputStarlightAscii;
  //ofstream 		fOutputStarlightAscii;
  TString 		lineString;
  TString		fDataPath;
  
  Int_t fNRun; // Added by Ashik 
  Double_t fImpactParameterPhi; // Added by Ashik 
  Double_t fImpactParameterMag; // Added by Ashik 
  Bool_t   fEMkickIsOn; // Added by Ashik 
  Bool_t   fGDRIsOn; // Added by Ashik 
  Bool_t   fEntanglementIsOn; // Added by Ashik 

  const Double_t 	neutron_M = 939.565;// MeV
  const Double_t 	hbarc = 0.1973269718;
  const Double_t 	hbarcmev = 197.3269718;
  const Double_t 	pi = 3.14159;
  const Double_t 	alpha = 1/137.035999074;
  const Int_t 		nSteps_impactPar = 120;
  const Int_t 		nSteps_GG = 10;
  const Int_t 		nSteps_R = 60;  
  const Int_t 		nSteps_Phi = 40;
  const Int_t		nSteps_Energy = 100;
  const Int_t		maxNeutrons = 50;
  Int_t 		nFluxes;
  
  const Double_t DeltaOmega      = 3.0;   // MeV  (try 2–5 MeV)
 const Double_t sigma_t_extra_fm= 0.0;   // fm   (optional extra timing jitter beyond b/γ, e.g. 1–5 fm)
 const Double_t sigma_cascade   = 0.30;  // rad  (increase for XnXn if desired)
   

  Int_t 	nucleus_Z; 
  Int_t 	nucleus_A;
  Double_t 	nucleus_R;
  Double_t 	beamGamma;
  Double_t 	gammaTarget;
  Double_t 	neutronSepThr;
  Double_t	saturationEnergy;

  TGraph	*gXsection;
  TH1D		*hXsection;
  TGraph 	*gPhotonFluxTable;
  TGraph	*gNucleusBreakupTable;
  TGraph2D	*gTwoPhotonFluxTable;
  TGraph	*gNuclearThickness;
  TGraph 	*gHadronicProbabilityTable;
  TGraph 	*gMeanNeutronN;
  TGraph 	*gWidthNeutronN;
  TF1		*fitMeanNeutronN;
  TF1		*fitWidthNeutronN;
  TH2D		*hENDF_2D;
  TH1D		*hENDF_1D;
  TFile		*fENDFFile; 
  TH2D		*hBranchingRatioMap;
  TH1D		*hEventBreakupMap;
  TGraph	*gUnitaryLeak;
  
  TTree 	*fEventTree;
  TClonesArray	*fParticles;
  TFile 	*fOutputFile;
  
  TGraph	*fImpProfile;
  TList		*fQAhistList;
  TH2D		*hNeutronMultiplicity;
  TH1D		*hEnergyGen;
  TH1D		*hKinEnergyGen;
  TH1D 		*hMomGen;
  TH1D 		*hPhiGen;
  TH1D 		*hThetaGen;
  TH1D 		*hEnergyBoosted;
  TH1D 		*hMomBoosted;
  TH1D 		*hPhiBoosted;
  TH1D 		*hThetaBoosted;
  TH1D 		*hNeutronRapidity;
  TH1D 		*hNeutronEta;
  TH1D 		*hEnergyBin;
  TH1D 		*hEnergyForNeutronMulti;
  TH1D		*hProbabilityXn;
  TH1D		*hPhotonK;
  TH1D		*hRapidityVM;
  TH1D		*hMassVM;
  Bool_t	kStoreQA;
  Bool_t	kStoreGen;
  
 // store the GDR neutrons
  TH1D* h_all;
  TH1D* h_east;
  TH1D* h_west;
  TH1D* h_prompt;
  TH1D* h_delayed;
  TH1D* h_emkick;

  TH1D* hPhotonEnergy_1n1n;
  TH1D* hPhotonEnergy_1nXn;
  TH1D* hPhotonEnergy_XnXn;


  TH2D* h2_dphi_nE_all;
  TH2D* h2_dphi_nE_prompt;
  TH2D* h2_dphi_nE_delayed;
  TH2D* h2_dphi_photonE_prompt;

  TH2D* h2_dphi_n_all;
  TH2D* h2_dphi_n_prompt;
  TH2D* h2_dphi_n_delayed;

  ClassDef(NeutronGenerator,1);
};

#endif
