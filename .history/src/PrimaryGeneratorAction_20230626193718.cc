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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
   double J(double p, double theta) {
double A = 0.14*pow(p, -2.7);
double B = 1. / (1. + 1.1*p*cos(theta)/115.);
double C = 0.054 / (1. + 1.1*p*cos(theta)/850.);
return A*(B+C);
    }

  PrimaryGeneratorAction::PrimaryGeneratorAction():G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), mu_plus(0), mu_minus(0)


						   
  {

    
    /* fMuonGen.SetUseHSphere();
    fMuonGen.SetHSphereRadius(2 fMuonGen.SetUseHSphere();
    fMuonGen.SetHSphereRadius(25.*cm);
    fMuonGen.SetHSphereCenterPosition({{0., 0., 20.*cm}});5.*cm);
    fMuonGen.SetHSphereCenterPosition({{0., 0., 20.*cm}});
    */
  fMuonGen.SetUseSky();
  fMuonGen.SetSkySize({{15.*cm, 15.*cm}});
  fMuonGen.SetSkyCenterPosition({{0., 0.,-.5*cm}});   
  fMuonGen.SetDifferentialFlux(&J);
  //fMuonGen.SetMaximumMomentum(50.*MeV);
   //fMuonGen.SetMaximumMomentum(1.*MeV);
  
  //fMuonGen.SetMinimumTheta(0.001);
  //fMuonGen.SetMaximumTheta(0.002);
  //fMuonGen.SetMinimumPhi(0.001);
  //fMuonGen.SetMaximumPhi(0.002);

  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  mu_minus = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  //mu_plus = G4ParticleTable::GetParticleTable()->FindParticle("mu+");


  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  // delete mu_plus;
  //delete mu_minus; 
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
 
  
   fMuonGen.Generate();
  std:: array<double,3> muon_pos = fMuonGen.GetGenerationPosition();
  double muon_ptot = -1.*fMuonGen.GetGenerationMomentum();
  double muon_theta = fMuonGen.GetGenerationTheta();
  double muon_phi = fMuonGen.GetGenerationPhi();

  fParticleGun->SetParticlePosition(G4ThreeVector(
						  muon_pos[0]*mm,
						  muon_pos[1]*mm,
						  muon_pos[2]*mm));

  fParticleGun->SetParticleMomentum( G4ParticleMomentum(muon_ptot*sin(muon_theta)*cos(muon_phi)*GeV,
							muon_ptot*sin(muon_theta)*sin(muon_phi)*GeV,
							muon_ptot*cos(muon_theta)*GeV) );


  
    
  /*  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    // G4double  envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    // G4double  envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  G4double size = 25*cm;
  G4double x0 = -size +2 * size*G4UniformRand();
  G4double y0 =-size +2 * size*G4UniformRand();
  G4double z0 = -1.*cm; //* envSizeZ;

  //  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  // fParticleGun->GeneratePrimaryVertex(anEvent);
  */
  if (fMuonGen.GetCharge() <0 ) {
    fParticleGun->SetParticleDefinition(mu_minus);
  } 
 // else { fParticleGun->SetParticleDefinition(mu_plus);}
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


