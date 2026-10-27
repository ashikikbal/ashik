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

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TH1.h"
#include "TDirectory.h"

#include "NeutronGenerator.h"
#include <iostream>


void run_gdrn(int nRun = 2) {
 
 TH1::AddDirectory(kFALSE);
  gErrorIgnoreLevel = kError;
 
  cout<<"nRun = "<<nRun<<endl;
  
   UInt_t nEvents = 1000000;  
 
  //  UInt_t nEvents = 100000;  
  //   int nRun = 10;
  //  double VMma0s = 3.09;// in GeV for J/psi
    double VMmass = 0.775;// in GeV for rho^0
    double fRapMin = -1.0;
    double fRapMax = 1.0;
    double VMrapidity = 0.0;
    double MeVtoGeV = 0.001;
 
  TFile *fOutputFile = new TFile(Form("SLoutput_%d.root",nRun),"RECREATE");
 // fOutputFile->cd();
  TClonesArray *fParticles = new TClonesArray("TParticle", 200);
  TTree *fEventTree = new TTree("fEventTree", "fEventTree");
  fEventTree ->Branch("fParticles", &fParticles);


     gROOT->LoadMacro("NeutronGenerator.cxx+"); // come back here
   //  NeutronGenerator* gen = new NeutronGenerator();
     std::cout << " ---> Got the NeutronGenerator class pointer...\n";
   //  gen->SetRunMode(NeutronGenerator::kFlatMultiplicity);
   //  gen->SetRunMode(NeutronGenerator::kInterface);
   //  gen->SetProductionMode(NeutronGenerator::kPhotonPomeron);
   //  gen->SetBeamParameters(NeutronGenerator::kAu197, 107.0);
   //  gen->Initialize();
    // gen->LoadENDF();


    TClonesArray *fNGparticles = NULL;
    NeutronGenerator *gen = new NeutronGenerator();
    //gen->SetHadronicInteractionModel(NeutronGenerator::kHardSphere);
     gen->SetRunMode(NeutronGenerator::kInterface);
    gen->SetBeamParameters(NeutronGenerator::kAu197, 107.3); // for Au+Au 200 GeV
   // gen->SetBeamParameters(NeutronGenerator::kPb208, 2694); // for Pb+Pb 5.02 TeV
    //gen->SetProductionMode(NeutronGenerator::kPhotonPomeron);
    gen->Initialize();
    gen->LoadENDF("hENDF.root", "hENDF_2D"); 
    gen->SetNRun(nRun);
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->Setup();
    
    gen->SetEMkickMode(true); // Set EMkickMode on
    gen->SetClassicalGDRMode(false); // Set GDRMode off
    gen->SetEntangledGDRMode(true); // Set EntanglementMode on

    TFile *f = TFile::Open("noon_flux_au.root", "READ"); // File for flux: d3N/d2bdk * 2pib * Exp(-P_H) 
   //  TFile *f = TFile::Open("noon_flux_pb.root", "READ"); // File for flux: d3N/d2bdk * 2pib * Exp(-P_H) 
    if (!f || f->IsZombie()) return;

    TH2D* h2 = dynamic_cast<TH2D*>(f->Get("h2fluxNoon_wPnohad2"));
    if (!h2) return;

 //   TH1D* h1 = h2->ProjectionX("projX");
 //   h1->Scale(1.0 / h1->Integral());

   TRandom3* myRand = new TRandom3(0);   
 //  TRandom3* myRand = new TRandom3(12345 + nRun);
   gRandom = myRand;

    double impact_parameter_b = 0.0;   // Magnitude |b|
    double impact_parameter_phi = 0.0; // Azimuthal angle of b vector
    double photonK = 0.0; // photon energy
     
   cout<<"Running production"<<endl;
for(UInt_t iEvent = 0; iEvent < nEvents; ++iEvent){
    if(iEvent%(nEvents/10) == 0){
      if(iEvent != 0){ printf("\033[1A"); printf("\033[K");}
      cout<<100*iEvent/nEvents<<"% "<<endl;
      }
     int nTotalPart = 0;

    // photonK = h1->GetRandom(); // in MeV
   // std::cout << "Before: K = " << photonK << std::endl;

     h2->GetRandom2(photonK, impact_parameter_b); // photonK in MeV and impact_parameter_b in fm     
     impact_parameter_phi = myRand->Uniform(0, 2 * TMath::Pi());
     gen->SetImpactParameterPhi(impact_parameter_phi);
     gen->SetImpactParameterMag(impact_parameter_b);
   //  cout << " b = "<<impact_parameter_b<< " phi = "<< impact_parameter_phi<<" K = " << photonK << std::endl;

     VMrapidity = TMath::Abs(TMath::Log(2*photonK*MeVtoGeV/VMmass));
   //  cout<<" y = "<<VMrapidity<<endl;
     if(VMrapidity>fRapMax || VMrapidity<fRapMin) continue;
    // cout<<" y = "<<VMrapidity<<endl;
 
   // std::cout << "After: K = " << photonK << std::endl;
 
 //---------------------- Get neutron distributions ------
  
     gen->GenerateEvent(photonK*MeVtoGeV);

// cout<<" Event generated!"<<endl;
// turn the following on if you want to write the TTree info in SLOutput.root
/*
 fNGparticles = gen->ImportParticles();
    for(Int_t i = 0; i<fNGparticles->GetEntriesFast(); i++){
      TParticle *part(dynamic_cast<TParticle*>(fNGparticles->At(i)));
      new((*fParticles)[nTotalPart++]) TParticle(*part);
      }

*/
//    fEventTree->Fill();
    fParticles->Clear("C");  
   gen->FinishEvent();

   } // end of event loop 
 

  fOutputFile->cd();
//  fEventTree->Write(); // turn the following on if you want to write the TTree info in SLOutput.root
 
  gen->FinishProduction();
  fOutputFile->Close();
  f->Close();

}

