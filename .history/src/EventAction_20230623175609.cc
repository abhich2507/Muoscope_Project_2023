//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4/B4a/src/EventAction.cc
/// \brief Implementation of the B4a::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "Randomize.hh"
#include <iomanip>
#include "SensitiveDetector.hh"
#include "G4PrimaryVertex.hh"
#include "ROOTManager.hh"
#include <vector>
#include <math.h>
using std::array;
using std::vector;


namespace B1
{
  //G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  // auto hce = event->GetHCofThisEvent();
  // auto hc = hce->GetHC(collId);
  // return hc;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
   auto myrootManager = ROOTManager::Instance();
  
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 void EventAction::EndOfEventAction(const G4Event* event)
{
  auto analysisManager = G4AnalysisManager::Instance();
   auto myrootManager = ROOTManager::Instance();
  
   
   //myrootManager->ROOTTreeStruct.NScintHit= nScintHit;

  if(hcID<0) hcID = G4SDManager::GetSDMpointer()->GetCollectionID("SensitiveDetectorHitsCollection");
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  SensitiveDetectorHitsCollection* fHitsCollection = 0;
  if(hce) fHitsCollection = (SensitiveDetectorHitsCollection*)(hce->GetHC(hcID));
  if(fHitsCollection) {
   G4int nHits = fHitsCollection->entries();
   // G4cout<<" numbER of enTRIES:   "<<nHits<<G4endl;
 
   G4float eDep=0.;
   G4int nScintHit=0;
 
  for (G4int iHit=0; iHit<nHits; ++iHit){
    SensitiveDetectorHit*  newHit = (*fHitsCollection)[iHit];
    //   myrootManager->ROOTTreeStruct.ScintHitE[myrootManager->ROOTTreeStruct.NScintHit] =
    //      (float_t)newHit->GetEdep();
       nScintHit+=1;
       
    G4float edep = newHit->GetEdep();
    G4int eventid= newHit->GetEventId();
    
   
    eDep+=edep;
    
    G4ThreeVector pos = newHit->GetPos();
    G4float Pz= newHit->GetPz();
    G4int PlaneNb =   newHit->GetPlaneNb();
    G4int pdgid = newHit->GetPDGID();
    G4int trackid = newHit->GetTrackID();
    G4float tracklength = newHit->GetTrackLength();
   // G4float tracklength_test= tracklength*pow(10,-10);
    //  G4cout<<PlaneNb<<G4endl;
    
    analysisManager->FillNtupleDColumn(0,0, edep);
    analysisManager->FillNtupleDColumn(0,1, pos[0]);
    analysisManager->FillNtupleDColumn(0,2, pos[1]);
    analysisManager->FillNtupleDColumn(0,3, pos[2]);
    analysisManager->FillNtupleDColumn(0,4,  Pz ); 
    analysisManager->FillNtupleDColumn(0,5,  PlaneNb );  
    analysisManager->FillNtupleDColumn(0,6,  pdgid );
    analysisManager->FillNtupleDColumn(0,7, eventid);
    analysisManager->AddNtupleRow(0);
     
     
      myrootManager->ROOTTreeStruct.NScintHit= nScintHit-1;
      myrootManager->ROOTTreeStruct.ScintHitE[myrootManager->ROOTTreeStruct.NScintHit]=edep;
      myrootManager->ROOTTreeStruct.ScintHitPosX[myrootManager->ROOTTreeStruct.NScintHit]= (G4float) pos[0];
      myrootManager->ROOTTreeStruct.ScintHitPosY[myrootManager->ROOTTreeStruct.NScintHit]= (G4float) pos[1];
      myrootManager->ROOTTreeStruct.ScintHitPosZ[myrootManager->ROOTTreeStruct.NScintHit]= (G4float) pos[2];
      myrootManager->ROOTTreeStruct.ScintHitPz[myrootManager->ROOTTreeStruct.NScintHit]=Pz;
      myrootManager->ROOTTreeStruct.ScintHitStation[myrootManager->ROOTTreeStruct.NScintHit]= (G4int) PlaneNb;
      myrootManager->ROOTTreeStruct.HitPDG[myrootManager->ROOTTreeStruct.NScintHit]= (G4int) pdgid ;
      myrootManager->ROOTTreeStruct.ScintHitPrimaryID[myrootManager->ROOTTreeStruct.NScintHit]= (G4int) eventid ;
      myrootManager->ROOTTreeStruct.ScintHitTrackID[myrootManager->ROOTTreeStruct.NScintHit]= (G4int) trackid;
      myrootManager->ROOTTreeStruct.ScintHitTrackLength[myrootManager->ROOTTreeStruct.NScintHit]= (G4float) tracklength;
      
        
  }

  // if (myrootManager != nullptr) {
    

   
  //    myrootManager->Fill();
    // } 
  EDep.clear();
 
  
  if(eDep!=0){
  EDep.push_back(eDep);
  }

  


   auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo == 0 ) || ( eventID % printModulo != 0 ) ) 
    return;
  
      G4float  KE=0;
      G4float Theta =0;
    G4int nGenPart =0;
    

     for(G4int i =0; i<event->GetNumberOfPrimaryVertex(); i++){
       nGenPart=event->GetPrimaryVertex(i)->GetNumberOfParticle();
        myrootManager->ROOTTreeStruct.NGenPart= nGenPart;
      for(G4int q=0; q<event->GetPrimaryVertex(i)->GetNumberOfParticle(); q++){
	
	auto  primary= event->GetPrimaryVertex(i)->GetPrimary(q);
      G4float  ke  = primary->GetKineticEnergy();
      G4float theta = primary->GetMomentum().theta();
       G4double phi = primary->GetMomentum().phi();
       G4int pdgid= primary->GetPDGcode();
       analysisManager->FillNtupleDColumn(1,0,ke);
       analysisManager->FillNtupleDColumn(1,1,theta);
       analysisManager->FillNtupleDColumn(1,2,phi);
       analysisManager->FillNtupleDColumn(1,3,pdgid);
      analysisManager->AddNtupleRow(1);
      nGenPart+=1;
      KE+=ke;
      Theta+=theta;
      //  G4cout<<"spor me   "<< nGenPart<<G4endl;
	// myrootManager->ROOTTreeStruct.GenPartE[ myrootManager->ROOTTreeStruct.NGenPart]= KE;
	myrootManager->ROOTTreeStruct.GenPartE= KE;
       myrootManager->ROOTTreeStruct.GenPartTheta=Theta;
                 }
      //   myrootManager->Fill();
     
     
  }//main for loop end
  
       
    myrootManager->Fill();
 
  
  }

  //  myrootManager->Save();
  
  
  }
  
  

// analysisManager->AddNtupleRow();
  // fHitsCollection = GetHitCollection(hce, "SensitiveDetectorHitsCollection");
  // if(fHitsCollection) {
  //const G4int nHits = fHitsCollection->entries();
  //for (G4int iHit=0; iHit<nHits; ++iHit){
  //    analysisManager->FillNtupleDColumn(0, energyDeposition);
      
  
  // Accumulate statistics
     
  
  //G4VSensitiveDetector* sensitiveDetector = sdManager->FindSensitiveDetector("SensitiveDetectorName");

  // get analysis manager
 
  
  

  // fill histograms
  // analysisManager->FillH1(0, fEnergyAbs);
  //analysisManager->FillH1(1, fEnergyGap);
  // analysisManager->FillH1(2, fTrackLAbs);
  //analysisManager->FillH1(3, fTrackLGap);

  // fill ntuple
  // analysisManager->FillNtupleDColumn(1,pz );
  //analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  //analysisManager->FillNtupleDColumn(1, fEnergyGap);
  //analysisManager->FillNtupleDColumn(2, fTrackLAbs);
  //analysisManager->FillNtupleDColumn(3, fTrackLGap);
  //analysisManager->AddNtupleRow();

  // Print per event (modulo n)
  //
  /* auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;

    G4cout
       << "   Absorber: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyAbs,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(fTrackLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(fTrackLGap,"Length")
       << G4endl;
       }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


