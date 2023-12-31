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
/// \file B1/include/DetectorConstruction.hh
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "tls.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UserLimits;

/// Detector construction class to define materials and geometry.

namespace B1
{

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() ;
    ~DetectorConstruction() override ;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  //Set Methods
   void SetMaxStep (G4double );
   void SetCheckOverlaps(G4bool );


    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  protected:
    G4LogicalVolume** fLogicPlane = nullptr;
    G4LogicalVolume** fLogicPlaneG1 = nullptr;
    G4LogicalVolume** fLogicPlaneG2= nullptr;
    G4LogicalVolume** fLogicPlaneAc= nullptr;
    G4LogicalVolume** fLogicPlaneSheet= nullptr;
    G4LogicalVolume** fLogicPlanePCB= nullptr;
     G4LogicalVolume** fLogicPlanestrips= nullptr;
     G4LogicalVolume** fLogicPlanePlates= nullptr;

    G4LogicalVolume* fScoringVolume = nullptr;
    G4UserLimits* fStepLimit = nullptr; // pointer to user step limits
    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
