#include "SensitiveDetectorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include <iomanip>

namespace B1
{

  G4ThreadLocal G4Allocator <SensitiveDetectorHit>*   SensitiveDetectorHitAllocator = nullptr;

  // SensitiveDetec  G4double fKe;torHit::SensitiveDetectorHit()
  //{}

  /////////////

  // SensitiveDetectorHit::~  SensitiveDetectorHit() {}

  ///////////////

  G4bool   SensitiveDetectorHit::operator==(const   SensitiveDetectorHit& right) const
    {
      return(this ==&right) ? true : false;
    }

  void SensitiveDetectorHit::Draw()
  { G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
      {
	G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
      }

  }

  void   SensitiveDetectorHit::Print()

  {
    //analysisManager->FillNtupleDColumn(0,fEdep);
    // analysisManager->FillNtupleDColumn(0,fEdep);
    // analysisManager->AddNtupleRow();
    
    /*   G4cout
      << "trackID: " << fTrackID << "      PlaneNb: "<< fPlaneNb << "     Edep: "<< G4BestUnit(fEdep, "Energy")
      << "     Position: " << G4BestUnit(fPos, "Length")<<"      Pz: "<<G4BestUnit(fPz, "Momentum")<<"      Kinetic Energy: "<< G4BestUnit(fKe, "Energy")<<"      Theta: "<<G4BestUnit(fTheta, "Angle")<<"    Phi: "<<G4BestUnit(fPhi,"Angle") <<G4endl; */
  }

}
