//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 27 17:21:01 2023 by ROOT version 6.26/06
// from TTree Default/Default
// found on file: RPC.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
struct Vector3D {
    double x;
    double y;
    double z;

    Vector3D(double _x=0., double _y=0., double _z=0.) : x(_x), y(_y), z(_z) {}
};

Vector3D crossProduct(const Vector3D& a, const Vector3D& b) {
    Vector3D result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

double dotProduct(const Vector3D& v1, const Vector3D& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double magnitude(const Vector3D& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double calculateAngle(const Vector3D& vA, const Vector3D& vB) {
    double dot = dotProduct(vA, vB);
    double magA = magnitude(vA);
    double magB = magnitude(vB);
    double theta = std::acos(dot / (magA * magB));
    
   Vector3D crossProductAB = crossProduct(vA, vB);
    double zComponent = crossProductAB.z;
    //int sign = (zComponent >= 0) ? 1 : -1;
    int sign;
/*if (zComponent >= 0) {
    sign = 1;
} else {
    sign = -1;
}*/

    
    // Calculate the scattering angle with sign
    //double scatteringAngleWithSign = sign * theta;
    double scatteringAngleWithSign;

    /*if (theta>=0 && theta<M_PI/2.){
      scatteringAngleWithSign=theta;
    }
    if (theta>M_PI/2.  && theta<=M_PI )
    {  scatteringAngleWithSign = M_PI-theta;}
*/
    //return scatteringAngleWithSign;
    return theta;
}


class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           NGenPart;
   Float_t         GenPartE;
   Float_t         GenPartTheta;
   Float_t         GenPartPhi;
   Int_t           NScintHit;
   Int_t           ScintHitPrimaryID[54];   //[NScintHit]
   Int_t           ScintHitTrackID[54];   //[NScintHit]
   Int_t           ScintHitTrackLength[54];   //[NScintHit]
   Float_t         ScintHitE[54];   //[NScintHit]
   Float_t         ScintHitPz[54];   //[NScintHit]
   Float_t         ScintHitPosX[54];   //[NScintHit]
   Float_t         ScintHitPosY[54];   //[NScintHit]
   Float_t         ScintHitPosZ[54];   //[NScintHit]
   Int_t           ScintHitStation[54];   //[NScintHit]
   Int_t           HitPDG[54];   //[NScintHit]

   // List of branches
   TBranch        *b_NGenPart;   //!
   TBranch        *b_GenPartE;   //!
   TBranch        *b_GenPartTheta;   //!
   TBranch        *b_GenPartPhi;   //!
   TBranch        *b_NScintHit;   //!
   TBranch        *b_ScintHitPrimaryID;   //!
   TBranch        *b_ScintHitTrackID;   //!
   TBranch        *b_ScintHitTrackLength;   //!
   TBranch        *b_ScintHitE;   //!
   TBranch        *b_ScintHitPz;   //!
   TBranch        *b_ScintHitPosX;   //!
   TBranch        *b_ScintHitPosY;   //!
   TBranch        *b_ScintHitPosZ;   //!
   TBranch        *b_ScintHitStation;   //!
   TBranch        *b_HitPDG;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RPC.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RPC.root");
      }
      f->GetObject("Default",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NGenPart", &NGenPart, &b_NGenPart);
   fChain->SetBranchAddress("GenPartE", &GenPartE, &b_GenPartE);
   fChain->SetBranchAddress("GenPartTheta", &GenPartTheta, &b_GenPartTheta);
   fChain->SetBranchAddress("GenPartPhi", &GenPartPhi, &b_GenPartPhi);
   fChain->SetBranchAddress("NScintHit", &NScintHit, &b_NScintHit);
   fChain->SetBranchAddress("ScintHitPrimaryID", ScintHitPrimaryID, &b_ScintHitPrimaryID);
   fChain->SetBranchAddress("ScintHitTrackID", ScintHitTrackID, &b_ScintHitTrackID);
   fChain->SetBranchAddress("ScintHitTrackLength", ScintHitTrackLength, &b_ScintHitTrackLength);
   fChain->SetBranchAddress("ScintHitE", ScintHitE, &b_ScintHitE);
   fChain->SetBranchAddress("ScintHitPz", ScintHitPz, &b_ScintHitPz);
   fChain->SetBranchAddress("ScintHitPosX", ScintHitPosX, &b_ScintHitPosX);
   fChain->SetBranchAddress("ScintHitPosY", ScintHitPosY, &b_ScintHitPosY);
   fChain->SetBranchAddress("ScintHitPosZ", ScintHitPosZ, &b_ScintHitPosZ);
   fChain->SetBranchAddress("ScintHitStation", ScintHitStation, &b_ScintHitStation);
   fChain->SetBranchAddress("HitPDG", HitPDG, &b_HitPDG);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
