
//Marwa
#include "ROOTManager.hh"

#include <TROOT.h>
#include <TFile.h>
#include "TBranch.h"
#include "TTree.h"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4UnitsTable.hh"
//#include"MuravesMessenger.hh"

namespace B1{
  // ROOTManager* ROOTManager::fgInstance = 0;
  ROOTManager* ROOTManager::fgInstance = nullptr;
  std::mutex ROOTManager::fgMutex;
  
  ROOTManager* ROOTManager::Instance()
{std::lock_guard<std::mutex> lock(fgMutex); // Lock the mutex

  if (fgInstance == nullptr) {
        fgInstance = new ROOTManager();
    }
  return fgInstance;
}      
  
ROOTManager::ROOTManager()
{ 
  fgInstance = this;
}

ROOTManager::~ROOTManager()
{
  // if ( ROOTFile ) delete ROOTFile;
  fgInstance = 0;  
}

void ROOTManager::Init()
{ 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
  // create ROOT file (dove alla fine salvare dati e grafici)
  //auto Label = MuravesMessenger::Instance()->GetRunLabel();
  //G4String input = "Muraves"+ Label + ".root";
  ROOTFile = new TFile("RPC.root", "RECREATE");
  if (!ROOTFile) {
    G4cout << " problem creating the ROOT TFile" << G4endl;
    return;
  }

  // creazione ttree p_incident_TOT
  ROOTTree = new TTree("Muraves_Generator", "Muraves_Generator");
   ROOTTreeG4Particle = new TTree("Muraves_G4Particle", "Muraves_G4Particle");
  ROOTTreeDefault = new TTree("Default", "Default");
  //tree for separetad station
  /*ROOTTreeStn_1 = new TTree("Station1", "Station1");
  ROOTTreeStn_2 = new TTree("Station2", "Station2");
  ROOTTreeStn_3 = new TTree("Station3", "Station3");
  ROOTTreeStn_4 = new TTree("Station4", "Station4");
  */




//Generator
  ROOTTree->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
   ROOTTree->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE[NGenPart]/F");
  ROOTTree->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta[NGenPart]/F");
 
  /* ROOTTree->Branch("Event", &ROOTTreeStruct.Event, "Event/I");
  ROOTTree->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
  ROOTTree->Branch("GenPartID", &ROOTTreeStruct.GenPartID, "GenPartID[NGenPart]/I");
 ROOTTree->Branch("GenPartPDG", &ROOTTreeStruct.GenPartPDG, "GenPartPDG[NGenPart]/I");
  ROOTTree->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE[NGenPart]/F");
  ROOTTree->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta[NGenPart]/F");
  ROOTTree->Branch("GenPartPhi", &ROOTTreeStruct.GenPartPhi, "GenPartPhi[NGenPart]/F");
   ROOTTree->Branch("NScintHit", &ROOTTreeStruct.NScintHit, "NScintHit/I");
  */

//G4Particule
   ROOTTreeG4Particle->Branch("NScintHit", &ROOTTreeStruct.NScintHit, "NScintHit/I");
   ROOTTreeG4Particle->Branch("ScintHitE", &ROOTTreeStruct.ScintHitE, "ScintHitE[NScintHit]/F");
    ROOTTreeG4Particle->Branch("ScintHitPosX", &ROOTTreeStruct.ScintHitPosX, "ScintHitPosX[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitPosY", &ROOTTreeStruct.ScintHitPosY, "ScintHitPosY[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitPosZ", &ROOTTreeStruct.ScintHitPosZ, "ScintHitPosZ[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitE", &ROOTTreeStruct.ScintHitPz, "ScintHitPz[NScintHit]/F");
   ROOTTreeG4Particle->Branch("ScintHitStation", &ROOTTreeStruct.ScintHitStation, "ScintHitStation[NScintHit]/I");
    ROOTTreeG4Particle->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NScintHit]/I");
   /*
  ROOTTreeG4Particle->Branch("ScintHitPrimaryID", &ROOTTreeStruct.ScintHitPrimaryID, "ScintHitPrimaryID[NScintHit]/I");
  ROOTTreeG4Particle->Branch("ScintHitE", &ROOTTreeStruct.ScintHitE, "ScintHitE[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitPosX", &ROOTTreeStruct.ScintHitPosX, "ScintHitPosX[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitPosY", &ROOTTreeStruct.ScintHitPosY, "ScintHitPosY[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitPosZ", &ROOTTreeStruct.ScintHitPosZ, "ScintHitPosZ[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitStation", &ROOTTreeStruct.ScintHitStation, "ScintHitStation[NScintHit]/I");
  ROOTTreeG4Particle->Branch("ScintHitModule", &ROOTTreeStruct.ScintHitModule, "ScintHitModule[NScintHit]/I");
  ROOTTreeG4Particle->Branch("ScintHitBar", &ROOTTreeStruct.ScintHitBar, "ScintHitBar[NScintHit]/I");
  ROOTTreeG4Particle->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NScintHit]/I");
  ROOTTreeG4Particle->Branch("ScintHitTime", &ROOTTreeStruct.ScintHitTime, "ScintHitTime[NScintHit]/F");

*/

  
 // ROOTTreeG4Particle->Branch("ScintHitKe", &ROOTTreeStruct.ScintHitKe, "ScintHitKe[NScintHit]/F");  
 // ROOTTreeG4Particle->Branch("ScintHitPhi",&ROOTTreeStruct.ScintHitPhi,"ScintHitPhi[NScintHit]/F");
  //ROOTTreeG4Particle->Branch("ScintHitPz",&ROOTTreeStruct.ScintHitPz,"ScintHitPz[NScintHit]/F");


   /*
  ROOTTreeG4Particle->Branch("ScintHitPosXWo", &ROOTTreeStruct.ScintHitPosXWo, "ScintHitPosXWo[NScintHit]/F");
  ROOTTreeG4Particle->Branch("ScintHitPosYWo", &ROOTTreeStruct.ScintHitPosYWo, "ScintHitPosYWo[NScintHit]/F");
   */

  
 // ROOTTreeG4Particle->Branch("ScintHitSmear", &ROOTTreeStruct.ScintHitSmear, "ScintHitSmear[NScintHit]/F");  
//Default
   // ROOTTreeDefault->Branch("Event", &ROOTTreeStruct.Event, "Event/I");
     ROOTTreeDefault->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
     ROOTTreeDefault->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE/F");
     ROOTTreeDefault->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta/F");
     
      ROOTTreeDefault->Branch("NScintHit", &ROOTTreeStruct.NScintHit, "NScintHit/I");
      ROOTTreeDefault->Branch("ScintHitPrimaryID", &ROOTTreeStruct.ScintHitPrimaryID, "ScintHitPrimaryID[NScintHit]/I");
      ROOTTreeDefault->Branch("ScintHitTrackID", &ROOTTreeStruct.ScintHitTrackID, "ScintHitTrackID[NScintHit]/I");
      ROOTTreeDefault->Branch("ScintHitTrackLength", &ROOTTreeStruct.ScintHitTrackLength, "ScintHitTrackLength[NScintHit]/I");
      ROOTTreeDefault->Branch("ScintHitE", &ROOTTreeStruct.ScintHitE, "ScintHitE[NScintHit]/F");
      ROOTTreeDefault->Branch("ScintHitPz", &ROOTTreeStruct.ScintHitPz, "ScintHitPz[NScintHit]/F");
      ROOTTreeDefault->Branch("ScintHitPosX", &ROOTTreeStruct.ScintHitPosX, "ScintHitPosX[NScintHit]/F");
      ROOTTreeDefault->Branch("ScintHitPosY", &ROOTTreeStruct.ScintHitPosY, "ScintHitPosY[NScintHit]/F");
      ROOTTreeDefault->Branch("ScintHitPosZ", &ROOTTreeStruct.ScintHitPosZ, "ScintHitPosZ[NScintHit]/F");
      ROOTTreeDefault->Branch("ScintHitStation", &ROOTTreeStruct.ScintHitStation, "ScintHitStation[NScintHit]/I");
      ROOTTreeDefault->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NScintHit]/I");
      
   /*
  ROOTTreeDefault->Branch("Event", &ROOTTreeStruct.Event, "Event/I");
  ROOTTreeDefault->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
  ROOTTreeDefault->Branch("GenPartID", &ROOTTreeStruct.GenPartID, "GenPartID[NGenPart]/I");
  ROOTTreeDefault->Branch("GenPartPDG", &ROOTTreeStruct.GenPartPDG, "GenPartPDG[NGenPart]/I");
  ROOTTreeDefault->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE[NGenPart]/F");
  ROOTTreeDefault->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta[NGenPart]/F");
  ROOTTreeDefault->Branch("GenPartPhi", &ROOTTreeStruct.GenPartPhi, "GenPartPhi[NGenPart]/F");


  ROOTTreeDefault->Branch("NScintHit", &ROOTTreeStruct.NScintHit, "NScintHit/I");
  ROOTTreeDefault->Branch("ScintHitPrimaryID", &ROOTTreeStruct.ScintHitPrimaryID, "ScintHitPrimaryID[NScintHit]/I");
  ROOTTreeDefault->Branch("ScintHitE", &ROOTTreeStruct.ScintHitE, "ScintHitE[NScintHit]/F");
  ROOTTreeDefault->Branch("ScintHitPosX", &ROOTTreeStruct.ScintHitPosX, "ScintHitPosX[NScintHit]/F");
  ROOTTreeDefault->Branch("ScintHitPosY", &ROOTTreeStruct.ScintHitPosY, "ScintHitPosY[NScintHit]/F");
  ROOTTreeDefault->Branch("ScintHitPosZ", &ROOTTreeStruct.ScintHitPosZ, "ScintHitPosZ[NScintHit]/F");
  ROOTTreeDefault->Branch("ScintHitStation", &ROOTTreeStruct.ScintHitStation, "ScintHitStation[NScintHit]/I");
  ROOTTreeDefault->Branch("ScintHitModule", &ROOTTreeStruct.ScintHitModule, "ScintHitModule[NScintHit]/I");
  ROOTTreeDefault->Branch("ScintHitBar", &ROOTTreeStruct.ScintHitBar, "ScintHitBar[NScintHit]/I");

   

  
  ROOTTreeDefault->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NScintHit]/I");
  ROOTTreeDefault->Branch("ScintHitPosXWo", &ROOTTreeStruct.ScintHitPosXWo, "ScintHitPosXWo[NScintHit]/F");
  ROOTTreeDefault->Branch("ScintHitPosYWo", &ROOTTreeStruct.ScintHitPosYWo, "ScintHitPosYWo[NScintHit]/F");  
  ROOTTreeDefault->Branch("ScintHitTime", &ROOTTreeStruct.ScintHitTime, "ScintHitTime[NScintHit]/F");

*/
  
 // ROOTTreeDefault->Branch("ScintHitKe", &ROOTTreeStruct.ScintHitKe, "ScintHitKe[NScintHit]/F");
 
//if you want to separet the station in the root output
// first station(1||9)
 /* ROOTTreeStn_1->Branch("NScintHitX1", &ROOTTreeStruct.NScintHitX1, "NScintHitX1/I");
  ROOTTreeStn_1->Branch("ScintHitPosXStationX1", &ROOTTreeStruct.ScintHitPosXStationX1, "ScintHitPosXStationX1[NScintHitX1]/F");
  ROOTTreeStn_1->Branch("ScintHitPosYStationX1", &ROOTTreeStruct.ScintHitPosYStationX1, "ScintHitPosYStationX1[NScintHitX1]/F");
  ROOTTreeStn_1->Branch("ScintHitPosZStationX1", &ROOTTreeStruct.ScintHitPosZStationX1, "ScintHitPosZStationX1[NScintHitX1]/F"); 
  ROOTTreeStn_1->Branch("ScintHitBarStationX1", &ROOTTreeStruct.ScintHitBarStationX1, "ScintHitBarStationX1[NScintHitX1]/I");
  ROOTTreeStn_1->Branch("ScintHitModuleStationX1", &ROOTTreeStruct.ScintHitModuleStationX1, "ScintHitModuleStationX1[NScintHitX1]/I");
  ROOTTreeStn_1->Branch("ScintHitEStationX1", &ROOTTreeStruct.ScintHitEStationX1, "ScintHitEStationX1[NScintHitX1]/F");
  ROOTTreeStn_1->Branch("ScintHitKeStationX1", &ROOTTreeStruct.ScintHitKeStationX1, "ScintHitKeStationX1[NScintHitX1]/F");

  ROOTTreeStn_1->Branch("NScintHitY1", &ROOTTreeStruct.NScintHitY1, "NScintHitY1/I");
  ROOTTreeStn_1->Branch("ScintHitPosXStationY1", &ROOTTreeStruct.ScintHitPosXStationY1, "ScintHitPosXStationY1[NScintHitY1]/F");
  ROOTTreeStn_1->Branch("ScintHitPosYStationY1", &ROOTTreeStruct.ScintHitPosYStationY1, "ScintHitPosYStationY1[NScintHitY1]/F");
  ROOTTreeStn_1->Branch("ScintHitPosZStationY1", &ROOTTreeStruct.ScintHitPosZStationY1, "ScintHitPosZStationY1[NScintHitY1]/F"); 
  ROOTTreeStn_1->Branch("ScintHitBarStationY1", &ROOTTreeStruct.ScintHitBarStationY1, "ScintHitBarStationY1[NScintHitY1]/I");
  ROOTTreeStn_1->Branch("ScintHitModuleStationY1", &ROOTTreeStruct.ScintHitModuleStationY1, "ScintHitModuleStationY1[NScintHitY1]/I");
  ROOTTreeStn_1->Branch("ScintHitEStationY1", &ROOTTreeStruct.ScintHitEStationY1, "ScintHitEStationY1[NScintHitY1]/F");
  ROOTTreeStn_1->Branch("ScintHitKeStationY1", &ROOTTreeStruct.ScintHitKeStationY1, "ScintHitKeStationY1[NScintHitY1]/F");

// station3
  ROOTTreeStn_3->Branch("NScintHitX3", &ROOTTreeStruct.NScintHitX3, "NScintHitX3/I");
  ROOTTreeStn_3->Branch("ScintHitPosXStationX3", &ROOTTreeStruct.ScintHitPosXStationX3, "ScintHitPosXStationX3[NScintHitX3]/F");
  ROOTTreeStn_3->Branch("ScintHitPosYStationX3", &ROOTTreeStruct.ScintHitPosYStationX3, "ScintHitPosYStationX3[NScintHitX3]/F");
  ROOTTreeStn_3->Branch("ScintHitPosZStationX3", &ROOTTreeStruct.ScintHitPosZStationX3, "ScintHitPosZStationX3[NScintHitX3]/F"); 
  ROOTTreeStn_3->Branch("ScintHitBarStationX3", &ROOTTreeStruct.ScintHitBarStationX3, "ScintHitBarStationX3[NScintHitX3]/I");
  ROOTTreeStn_3->Branch("ScintHitModuleStationX3", &ROOTTreeStruct.ScintHitModuleStationX3, "ScintHitModuleStationX3[NScintHitX3]/I");
  ROOTTreeStn_3->Branch("ScintHitEStationX3", &ROOTTreeStruct.ScintHitEStationX3, "ScintHitEStationX3[NScintHitX3]/F");
  ROOTTreeStn_3->Branch("ScintHitKeStationX3", &ROOTTreeStruct.ScintHitKeStationX3, "ScintHitKeStationX3[NScintHitX3]/F");

  ROOTTreeStn_3->Branch("NScintHitY3", &ROOTTreeStruct.NScintHitY3, "NScintHitY3/I");
  ROOTTreeStn_3->Branch("ScintHitPosXStationY3", &ROOTTreeStruct.ScintHitPosXStationY3, "ScintHitPosXStationY3[NScintHitY3]/F");
  ROOTTreeStn_3->Branch("ScintHitPosYStationY3", &ROOTTreeStruct.ScintHitPosYStationY3, "ScintHitPosYStationY3[NScintHitY3]/F");
  ROOTTreeStn_3->Branch("ScintHitPosZStationY3", &ROOTTreeStruct.ScintHitPosZStationY3, "ScintHitPosZStationY3[NScintHitY3]/F"); 
  ROOTTreeStn_3->Branch("ScintHitBarStationY3", &ROOTTreeStruct.ScintHitBarStationY3, "ScintHitBarStationY3[NScintHitY3]/I");
  ROOTTreeStn_3->Branch("ScintHitModuleStationY3", &ROOTTreeStruct.ScintHitModuleStationY3, "ScintHitModuleStationY3[NScintHitY3]/I");
  ROOTTreeStn_3->Branch("ScintHitEStationY3", &ROOTTreeStruct.ScintHitEStationY3, "ScintHitEStationY3[NScintHitY3]/F");
  ROOTTreeStn_3->Branch("ScintHitKeStationY3", &ROOTTreeStruct.ScintHitKeStationY3, "ScintHitKeStationY3[NScintHitY3]/F");

//station4
  ROOTTreeStn_4->Branch("NScintHitX4", &ROOTTreeStruct.NScintHitX4, "NScintHitX4/I");
  ROOTTreeStn_4->Branch("ScintHitPosXStationX4", &ROOTTreeStruct.ScintHitPosXStationX4, "ScintHitPosXStationX4[NScintHitX4]/F");
  ROOTTreeStn_4->Branch("ScintHitPosYStationX4", &ROOTTreeStruct.ScintHitPosYStationX4, "ScintHitPosYStationX4[NScintHitX4]/F");
  ROOTTreeStn_4->Branch("ScintHitPosZStationX4", &ROOTTreeStruct.ScintHitPosZStationX4, "ScintHitPosZStationX4[NScintHitX4]/F"); 
  ROOTTreeStn_4->Branch("ScintHitBarStationX4", &ROOTTreeStruct.ScintHitBarStationX4, "ScintHitBarStationX4[NScintHitX4]/I");
  ROOTTreeStn_4->Branch("ScintHitModuleStationX4", &ROOTTreeStruct.ScintHitModuleStationX4, "ScintHitModuleStationX4[NScintHitX4]/I");
  ROOTTreeStn_4->Branch("ScintHitEStationX4", &ROOTTreeStruct.ScintHitEStationX4, "ScintHitEStationX4[NScintHitX4]/F");
  ROOTTreeStn_4->Branch("ScintHitKeStationX4", &ROOTTreeStruct.ScintHitKeStationX4, "ScintHitKeStationX4[NScintHitX4]/F");

  ROOTTreeStn_4->Branch("NScintHitY4", &ROOTTreeStruct.NScintHitY4, "NScintHitY4/I");
  ROOTTreeStn_4->Branch("ScintHitPosXStationY4", &ROOTTreeStruct.ScintHitPosXStationY4, "ScintHitPosXStationY4[NScintHitY4]/F");
  ROOTTreeStn_4->Branch("ScintHitPosYStationY4", &ROOTTreeStruct.ScintHitPosYStationY4, "ScintHitPosYStationY4[NScintHitY4]/F");
  ROOTTreeStn_4->Branch("ScintHitPosZStationY4", &ROOTTreeStruct.ScintHitPosZStationY4, "ScintHitPosZStationY4[NScintHitY4]/F"); 
  ROOTTreeStn_4->Branch("ScintHitBarStationY4", &ROOTTreeStruct.ScintHitBarStationY4, "ScintHitBarStationY4[NScintHitY4]/I");
  ROOTTreeStn_4->Branch("ScintHitModuleStationY4", &ROOTTreeStruct.ScintHitModuleStationY4, "ScintHitModuleStationY4[NScintHitY4]/I");
  ROOTTreeStn_4->Branch("ScintHitEStationY4", &ROOTTreeStruct.ScintHitEStationY4, "ScintHitEStationY4[NScintHitY4]/F");
  ROOTTreeStn_4->Branch("ScintHitKeStationY4", &ROOTTreeStruct.ScintHitKeStationY4, "ScintHitKeStationY4[NScintHitY4]/F");

//station2(0||8)
  ROOTTreeStn_2->Branch("NScintHitX2", &ROOTTreeStruct.NScintHitX2, "NScintHitX2/I");
  ROOTTreeStn_2->Branch("ScintHitPosXStationX2", &ROOTTreeStruct.ScintHitPosXStationX2, "ScintHitPosXStationX2[NScintHitX2]/F");
  ROOTTreeStn_2->Branch("ScintHitPosYStationX2", &ROOTTreeStruct.ScintHitPosYStationX2, "ScintHitPosYStationX2[NScintHitX2]/F");
  ROOTTreeStn_2->Branch("ScintHitPosZStationX2", &ROOTTreeStruct.ScintHitPosZStationX2, "ScintHitPosZStationX2[NScintHitX2]/F");
  ROOTTreeStn_2->Branch("ScintHitBarStationX2", &ROOTTreeStruct.ScintHitBarStationX2, "ScintHitBarStationX2[NScintHitX2]/I");
  ROOTTreeStn_2->Branch("ScintHitModuleStationX2", &ROOTTreeStruct.ScintHitModuleStationX2, "ScintHitModuleStationX2[NScintHitX2]/I");
  ROOTTreeStn_2->Branch("ScintHitEStationX2", &ROOTTreeStruct.ScintHitEStationX2, "ScintHitEStationX2[NScintHitX2]/F");
  ROOTTreeStn_2->Branch("ScintHitKeStationX2", &ROOTTreeStruct.ScintHitKeStationX2, "ScintHitKeStationX2[NScintHitX2]/F");

  ROOTTreeStn_2->Branch("NScintHitY2", &ROOTTreeStruct.NScintHitY2, "NScintHitY2/I");
  ROOTTreeStn_2->Branch("ScintHitPosXStationY2", &ROOTTreeStruct.ScintHitPosXStationY2, "ScintHitPosXStationY2[NScintHitY2]/F");
  ROOTTreeStn_2->Branch("ScintHitPosYStationY2", &ROOTTreeStruct.ScintHitPosYStationY2, "ScintHitPosYStationY2[NScintHitY2]/F");
  ROOTTreeStn_2->Branch("ScintHitPosZStationY2", &ROOTTreeStruct.ScintHitPosZStationY2, "ScintHitPosZStationY2[NScintHitY2]/F");
  ROOTTreeStn_2->Branch("ScintHitBarStationY2", &ROOTTreeStruct.ScintHitBarStationY2, "ScintHitBarStationY2[NScintHitY2]/I");
  ROOTTreeStn_2->Branch("ScintHitModuleStationY2", &ROOTTreeStruct.ScintHitModuleStationY2, "ScintHitModuleStationY2[NScintHitY2]/I");
  ROOTTreeStn_2->Branch("ScintHitEStationY2", &ROOTTreeStruct.ScintHitEStationY2, "ScintHitEStationY2[NScintHitY2]/F");
  ROOTTreeStn_2->Branch("ScintHitKeStationY2", &ROOTTreeStruct.ScintHitKeStationY2, "ScintHitKeStationY2[NScintHitY2]/F");*/


}

void ROOTManager::Save()
{ 
  if (ROOTFile) {
    ROOTFile->Write();       
    ROOTFile->Close();       
    G4cout << "ROOT Tree closed" << G4endl;
  }
}

void ROOTManager::Fill()
{
  //ROOTTree->Fill();
  //ROOTTreeG4Particle->Fill();
   ROOTTreeDefault->Fill();
  //separeted station
 /* ROOTTreeStn_1->Fill();
  ROOTTreeStn_2->Fill();
  ROOTTreeStn_3->Fill();
  ROOTTreeStn_4->Fill();*/
}


}
