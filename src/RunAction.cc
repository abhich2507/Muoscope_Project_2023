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
/// \file B4/B4a/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "ROOTManager.hh"
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction)
  : fEventAction(eventAction)
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();
  auto myrootManager = ROOTManager::Instance();
  if (myrootManager != nullptr) {
    // Access methods or variables of myrootManager
    myrootManager->Init();
  }
  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  analysisManager->OpenFile("RPC_project.root");
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  // analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 5*MeV);
  //analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
  //analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
  // analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);

  // Creating ntuple
  //
  std::vector<G4double> EDep;
  if (fEventAction) {
    analysisManager->CreateNtuple("Hit", "Hit");
    analysisManager->CreateNtupleDColumn("edep");
    analysisManager->CreateNtupleDColumn("x");
    analysisManager->CreateNtupleDColumn("y");
    analysisManager->CreateNtupleDColumn("z");
    analysisManager->CreateNtupleDColumn("Pz");
    analysisManager->CreateNtupleDColumn("Plane");
    analysisManager->CreateNtupleDColumn("PDG");
    analysisManager->CreateNtupleDColumn("EventId");
    analysisManager->CreateNtupleDColumn("EDep", fEventAction->GetEDep());
    //  analysisManager->CreateNtupleDColumn("GenX");
    // analysisManager->CreateNtupleDColumn("GenY");
    // analysisManager->CreateNtupleDColumn("GenZ");
    
    
  
    /*
   analysisManager->CreateNtupleDColumn("x0");
  analysisManager->CreateNtupleDColumn("y0");
  analysisManager->CreateNtupleDColumn("z0");
  
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("y1");
  analysisManager->CreateNtupleDColumn("z1");

   analysisManager->CreateNtupleDColumn("x2");
  analysisManager->CreateNtupleDColumn("y2");
  analysisManager->CreateNtupleDColumn("z2");
  
   analysisManager->CreateNtupleDColumn("x3");
  analysisManager->CreateNtupleDColumn("y3");
  analysisManager->CreateNtupleDColumn("z3");
  analysisManager->CreateNtupleDColumn("edepvector", fEventAction->GetEdep());
  
 analysisManager->CreateNtupleDColumn("eventID");
  analysisManager->CreateNtupleDColumn("Pz");
  analysisManager->CreateNtupleDColumn("Phi");
  analysisManager->CreateNtupleDColumn("edep");
  
  analysisManager->CreateNtupleDColumn("theta");
  // analysisManager->CreateNtupleDColumn("Plane");
  analysisManager->CreateNtupleDColumn("x0");
  analysisManager->CreateNtupleDColumn("y0");
  analysisManager->CreateNtupleDColumn("z0");
  
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("y1");
  analysisManager->CreateNtupleDColumn("z1");

   analysisManager->CreateNtupleDColumn("x2");
  analysisManager->CreateNtupleDColumn("y2");
  analysisManager->CreateNtupleDColumn("z2");
  
   analysisManager->CreateNtupleDColumn("x3");
  analysisManager->CreateNtupleDColumn("y3");
  analysisManager->CreateNtupleDColumn("z3");
  */
  
  analysisManager->FinishNtuple();


   analysisManager->CreateNtuple("Gen", "Gen");
    analysisManager->CreateNtupleDColumn("Ke");
  analysisManager->CreateNtupleDColumn("theta");
  analysisManager->CreateNtupleDColumn("phi");
 analysisManager->CreateNtupleDColumn("PDGID");
 analysisManager->CreateNtupleDColumn("x-position");
 analysisManager->CreateNtupleDColumn("y-position");
  analysisManager->CreateNtupleDColumn("GenM");

 analysisManager->FinishNtuple();


  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
   auto myrootManager = ROOTManager::Instance();
   // myrootManager->Init();
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  // auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  // G4String fileName = "B4.root";
  // Other supported output types:
  // G4String fileName = "B4.csv";
  // G4String fileName = "B4.hdf5";
  // G4String fileName = "B4.xml";
  // analysisManager->OpenFile(fileName);
  //  G4cout << "Using " << analysisManager->GetType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  auto analysisManager = G4AnalysisManager::Instance();
    auto myrootManager = ROOTManager::Instance();
  
  /* 
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl;
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl;
    }

    G4cout << " EAbs : mean = "
       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

    G4cout << " EGap : mean = "
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

    G4cout << " LAbs : mean = "
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;

    G4cout << " LGap : mean = "
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
      } */

  // save histograms & ntuple
  //
    // if (myrootManager != nullptr) {
   
    myrootManager->Save();
    // } 
     analysisManager->Write();
    analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
