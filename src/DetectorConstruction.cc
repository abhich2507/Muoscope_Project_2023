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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4GDMLParser.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 150*cm, env_sizeZ = 150*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");


 
  //creating CO2
  
G4double a = 12*g/mole;  //carbon
G4int z;
G4Element* elC = new G4Element("Carbon", "C", z=6, a);
 
 a = 16.00*g/mole;     //
G4Element* elO = new G4Element("Oxygen", "O", z=8, a);
G4Element* elSi = new G4Element("Silicon", "Si", z=14., a= 28.0855*g/mole);
  double density = 2.64 *g/cm3;
  int nel;
  int natoms;
    G4Material* Quartz = new G4Material("Quartz", density,  nel= 2);
     Quartz-> AddElement(elSi,  natoms=1);
     Quartz-> AddElement(elO,  natoms=2);

density;
G4int n;
G4Material* CO2 = new G4Material("Carbon Dioxide", density=1.87*g/cm3 ,n=2);
G4int nAtoms ;
CO2->AddElement(elC, nAtoms=1); 
CO2->AddElement(elO, nAtoms=2);

 
a = 39.948*g/mole;
 G4Element* elAr = new G4Element("Argon","Ar", z=18, a);

G4double ncomponents = 2;
G4Material* Aerog = new G4Material("Aerogel",density=1.83*kg/m3,ncomponents);
G4double fractionmassAr=70*perCent;
G4double fractionmassCO2=30*perCent;
Aerog->AddElement(elAr, fractionmassAr);
Aerog->AddMaterial(CO2, fractionmassCO2);


 a= 102.032*g/mole;
 density = 0.00425 * g/cm3; // Set the desired density for R134a
  // Create the gas mixture material
 G4Material* r134a = new G4Material("r134a", density, 3);

   // Add the component to the gas mixture
  r134a->AddElement(nist->FindOrBuildElement("C"), 2); // C2
  r134a->AddElement(nist->FindOrBuildElement("H"), 2); // H2
  r134a->AddElement(nist->FindOrBuildElement("F"), 4); // F4

  //butane
  a=58.124*g/mole;
   G4Material* butane = new G4Material("butane", 2.489 *mg/cm3, 2);

    // Add the components to the butane material
    butane->AddElement(nist->FindOrBuildElement("C"), 4); // Carbon (C)
    butane->AddElement(nist->FindOrBuildElement("H"), 10); // Hydrogen (H)


    //SF6
     density = 6.43*mg/cm3;
    G4Material* hexafluoride = new G4Material("hexafluoride",density,2);
    hexafluoride->AddElement(nist->FindOrBuildElement("S"),1);
    hexafluoride->AddElement(nist->FindOrBuildElement("F"),6);
    
    //mixture of r134a(95%) + Butane(5%)
    G4Material* Aerog1 = new G4Material("Aerogel1",4.18364 *mg/cm3, 3);
    Aerog1->AddMaterial(r134a, 95.2*perCent);
    Aerog1->AddMaterial(butane,4.5*perCent);
    Aerog1->AddMaterial(hexafluoride,0.3*perCent); 

 G4Material*  sc = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");///scintillator
   G4Material* sc1 = new G4Material("sc1",1.05 *g/cm3, 1);
    sc1->AddMaterial(sc, 100*perCent);


  G4NistManager* nistManager = G4NistManager::Instance();

    // Define the galactic vacuum material
    G4String materialName = "GalacticVacuum";
     density = 1.e-25 * g/cm3; // Very low density, effectively zero
    G4double pressure = 1.e-19 * pascal; // Very low pressure, effectively zero
    G4double temperature = 2.7 * kelvin; // Cosmic microwave background temperature

    // Create the galactic vacuum material
    G4Material* galacticVacuumMaterial = new G4Material(materialName, density, 1, kStateGas, temperature, pressure);
    galacticVacuumMaterial->AddMaterial(nistManager->FindOrBuildMaterial("G4_Galactic"), 1.0); 

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1*env_sizeXY;
  G4double world_sizeZ  = 1*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");


  G4Box* solidWorld =
    new G4Box("World",                       //its name
      .1*world_sizeXY, 0.1*world_sizeXY, 0.3*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,galacticVacuumMaterial,  "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //
  // Envelope
  
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        world_mat,             //its material
                        "Envelope");         //its name

     new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  

  G4int fNbOfPlanes = 4;
  G4double G10= 0.25*cm;
  G4double z0=0.5*cm;
  G4double G20= 0.75*cm;

  G4double G11= 27.4*cm;
  G4double z1=27.65*cm;
  G4double G12= 27.90*cm;

  /*G4double G11= 10.4*cm;
  G4double z1=10.5*cm;
  G4double G12= 10.6*cm;*/

  G4double G1_pos;
  G4double z_pos;
  G4double G2_pos;

 fLogicPlaneG1 = new G4LogicalVolume*[fNbOfPlanes] ;
 fLogicPlane = new G4LogicalVolume*[fNbOfPlanes];
 fLogicPlaneG2 = new G4LogicalVolume*[fNbOfPlanes] ;

  //Building four parrallel planes
  for (G4int copyNo=0; copyNo<fNbOfPlanes; copyNo++) {
   if(copyNo<2){ 
       z_pos=  z0 + (10*cm)*copyNo;
       G1_pos=  G10 + (10*cm)*copyNo;
       G2_pos=  G20 + (10*cm)*copyNo;
       }

    if(copyNo>=2){ 
        z_pos=  z1 + (10*cm)*copyNo;
        G1_pos=  G11 + (10*cm)*copyNo;
        G2_pos=  G12 + (10*cm)*copyNo;
        }
    

   
       /*z_pos=  z0 + (2.5*cm)*copyNo;
       G1_pos=  G10 + (2.5*cm)*copyNo;
       G2_pos=  G20 + (2.5*cm)*copyNo;*/
       

    

        
  G4ThreeVector pos1 = G4ThreeVector(0, 0*cm, z_pos);
  G4ThreeVector posG1 = G4ThreeVector(0, 0*cm, G1_pos);
  G4ThreeVector posG2 = G4ThreeVector(0, 0*cm, G2_pos);



  //a 2D plane
  
   G4Box * solidShape1  =
     new G4Box("Shape1",
	       (16/2)*cm,
	       (16/2)*cm,
	       (0.2/2)*cm);

  G4Box * glass1  =
     new G4Box("glass1",
	       (16/2)*cm,
	       (16/2)*cm,
	       (0.3/2)*cm);

  G4Box * glass2  =
     new G4Box("glass2",
	       (16/2)*cm,
	       (16/2)*cm,
	       (0.3/2)*cm);
  


 
 
   fLogicPlane[copyNo]=
    new G4LogicalVolume(solidShape1,         //its solid
                        Aerog1,          //its material
                        "Shape1");           //its name

   new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
		    // logicShape1,
		    fLogicPlane[copyNo],    //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    copyNo,                       //copy number
                    checkOverlaps);          //overlaps checking

  fLogicPlaneG1[copyNo]=
    new G4LogicalVolume(glass1,         //its solid
                        Quartz,          //its material
                        "glass1");           //its name

  
  new G4PVPlacement(0,                       //no rotation
                    posG1,                    //at position
		    // logicShape1,
		    fLogicPlaneG1[copyNo],    //its logical volume
                    "glass1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    copyNo,                       //copy number
                    checkOverlaps);          //overlaps checking

  
  fLogicPlaneG2[copyNo]=
    new G4LogicalVolume(glass2,         //its solid
                        Quartz,          //its material
                        "glass2");           //its name

 

  new G4PVPlacement(0,                       //no rotation
                    posG2,                    //at position
		    // logicShape1,
		    fLogicPlaneG2[copyNo],    //its logical volume
                    "glass2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    copyNo,                       //copy number
                    checkOverlaps);          //overlaps checking
  

 G4VisAttributes *red = new G4VisAttributes(G4Colour(1.,0.,0.));

  fLogicPlane[copyNo]->SetVisAttributes(red);
  fScoringVolume = fLogicPlane[copyNo]; 
}

 // G4ThreeVector pos4 = G4ThreeVector(10*cm,10*cm, 10*cm);
 G4ThreeVector pos4 = G4ThreeVector(10*cm,10*cm, 10*cm);
  G4Box * solidShape4  = new G4Box("Shape4", (10/2)*cm, (10/2)*cm, (8/2)*cm);

  G4ThreeVector pos6 = G4ThreeVector(-10*cm,-10*cm, 10*cm);
  G4Box * solidShape6  = new G4Box("Shape6", (10/2)*cm, (10/2)*cm, (6/2)*cm);

G4VSolid* box = new G4Box("Box",10*cm,10*cm,5*cm);
G4VSolid* cylinder = new G4Tubs("Cylinder",0.,5.*cm,2.*cm,0.,2*M_PI*rad);
G4VSolid* myunion = new G4UnionSolid("Box+Cylinder", box, cylinder); 
  
  G4double A = 207.2 *g/mole;//lead
  density = 11.35 *g/cm3;//lead
 G4Material* Lead =  new G4Material("Lead", z=82, A, density);//lead
 

 A = 27 *g/mole;
density = 2.7 *g/cm3;
G4Material* Aluminium =  new G4Material("Aluminium", z=13, A, density);


 A = 238.03 *g/mole;
  density = 18.95 *g/cm3;
 G4Material* Uranium =  new G4Material("Uranium", z=92, A, density);

  A = 55.85 *g/mole;
  density = 7.9 *g/cm3;
 G4Material* Iron =  new G4Material("Iron", z=26, A, density);




  /*G4LogicalVolume* logicShape4 = new G4LogicalVolume(solidShape4, Uranium, "Shape4" );
   new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logicShape4,             //its logical volume
                    "Shape4",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

 G4LogicalVolume* logicShape6 = new G4LogicalVolume(solidShape6,Aluminium, "Shape6" );
   new G4PVPlacement(0,                       //no rotation
                    pos6,                    //at position
                    logicShape6,             //its logical volume
                    "Shape6",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  */
                    

/*G4LogicalVolume* logicShape7  = new G4LogicalVolume(myunion,Uranium, "Box+Cylinder" );
   new G4PVPlacement(0,                       //no rotation
                    pos6,                    //at position
                    logicShape7,             //its logical volume
                    "Box+Cylinder",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);  */
   
  
  //always return the physical World

  //////////////////////////////////////////////////////////////////////////////////////////////
 


    ////////////////////////////////////////////////////////////////////
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  void DetectorConstruction::ConstructSDandField()
  {
     G4String SensitiveDetectorname = "SensitiveDetectorname";
   auto aSensitiveDetector = new SensitiveDetector(SensitiveDetectorname, "SensitiveDetectorHitsCollection");
   
   G4SDManager::GetSDMpointer()->AddNewDetector(aSensitiveDetector);
   SetSensitiveDetector("Shape1",aSensitiveDetector, true); 
  }

  
  void DetectorConstruction::SetMaxStep(G4double maxStep)
  {
    if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
  }


  
  void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}

}
