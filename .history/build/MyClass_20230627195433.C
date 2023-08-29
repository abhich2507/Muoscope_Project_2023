#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include  "TVector.h"

Double_t MyFitFunction(Double_t* x, Double_t* par) {
       Double_t value =par[0]*TMath::Power(TMath::Cos(x[0]), par[1]);  // cos^2(x)
       
    
    return value;
    }

void MyClass::Loop()
{ 
  vector<float>           S0;
vector<float>          S1;
vector<float>           S2;
vector<float>         S3;
vector<float> Sx[4]= {S0,S1,S2,S3};  
vector<float> Sy[4]= {S0,S1,S2,S3};  
vector<float> Sz[4]= {S0,S1,S2,S3};  

   
   
   
vector<float> Planenb;

fstream outFile("test11.csv",ios::out);
vector<float>           HE0;
vector<float>          HE1;
vector<float>           HE2;
vector<float>         HE3;
vector<float> HE[4]= {HE0,HE1,HE2,HE3};

 TH1F* histAngle = new TH1F("histAngle","scattering angle",500,-.10,.10);

   TH1F* histGenE = new TH1F("histGenE","KE",100,-100,14000);
   histGenE->GetXaxis()->SetTitle("Kinetic Energy[MeV]");
   histGenE->GetYaxis()->SetTitle("N");
    TH1F* histGenTheta = new TH1F("histGenTheta","theta",100,-1,1.6);
   histGenTheta->GetXaxis()->SetTitle("Radian");
   histGenTheta->GetYaxis()->SetTitle("N");


    char nameE[15], titleE[20];
   TH1F *HlistE[4];      // create an array of Histograms
   TH1F* hStationNumberE;                 // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++) {
      sprintf(nameE,"hStationNumberE%d",i);
      sprintf(titleE,"station:%d",i);
      hStationNumberE = new TH1F(nameE,titleE,100,-0.0001,0.0030);
      hStationNumberE->GetXaxis()->SetTitle("Deposited Energy[MeV]");
      hStationNumberE->GetYaxis()->SetTitle("numHits");
      HlistE[i]=hStationNumberE;
   }

 char nameX[15], titleX[20];
   TH1F *HlistX[4];      // create an array of Histograms
   TH1F* hStationNumberX;                 // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++) {
      sprintf(nameX,"hStationNumberX%d",i);
      sprintf(titleX,"station:%d",i);
      hStationNumberX = new TH1F(nameX,titleX,520,-260,260);
      hStationNumberX->GetXaxis()->SetTitle("mm");
      hStationNumberX->GetYaxis()->SetTitle("numHits");
      HlistX[i]=hStationNumberX;
   }

    char nameY[15], titleY[20];
   TH1F *HlistY[10];      // create an array of Histograms
   TH1F* hStationNumberY;                 // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++) {
      sprintf(nameY,"hStationNumberY%d",i);
      sprintf(titleY,"station:%d",i);
      hStationNumberY = new TH1F(nameY,titleY,100,-260,260);
      hStationNumberY->GetXaxis()->SetTitle("mm");
      hStationNumberY->GetYaxis()->SetTitle("numHits");
      HlistY[i]=hStationNumberY;
   }

   
 
 char nameZ[15], titleZ[20];
   TH1F *HlistZ[10];      // create an array of Histograms
   TH1F* hStationNumberZ;                 // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++) {
      sprintf(nameZ,"hStationNumberZ%d",i);
      sprintf(titleZ,"station:%d",i);
      hStationNumberZ = new TH1F(nameZ,titleZ,100,-10,350);
      hStationNumberZ->GetXaxis()->SetTitle("mm");
      hStationNumberZ->GetYaxis()->SetTitle("numHits");
      HlistZ[i]=hStationNumberZ;
   }

       TCanvas *c1 = new TCanvas("c1","Histograms related Planes edep ");
        c1->Divide(2,2);
       TCanvas *c2 = new TCanvas("c2", "Histograms related to the Planes posZ ");
       c2->Divide(2,2);
       TCanvas *c3 = new TCanvas("c3", "Histograms related to  the Planes posY ");
       c3->Divide(2,2);
       TCanvas *c4 = new TCanvas("c4", "Histograms related to  the Planes posX");
       c4->Divide(2,2);

       TCanvas *c5 = new TCanvas("c5","Histograms related Generator");
      c5->Divide(1,2);
      
      TCanvas *c6 = new TCanvas("c6","Histograms related to Scattering angle");
     TCanvas *c8 = new TCanvas("c8","Histograms related fit");


   int t=0;
  
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Planenb.clear();
      if (NScintHit!=0  ){
         for(int i=0; i< NScintHit; i++){
            if(abs(HitPDG[i])==13){
             //  cout<<"first loop event "<<t<<endl;
               Planenb.push_back(ScintHitStation[i]);
                   }
            }
      }
       auto it = std::find(Planenb.begin(), Planenb.end(), 3);
       if (it!= Planenb.end()){
      	for (int j=0; j<4; j++) {
	   double E=0.;
	   double Z=0.;
      double Y=0.;
      double X=0.;
      double t=0;
      double H=0.;
	  
      if (NScintHit!=0){
 
       
      for(int i=0; i< NScintHit; i++){
	  if(abs(HitPDG[i])==13){
     
     //  cout<<"second loop event "<<t<<endl;
     // Planenb.push_back(ScintHitStation[i]);
       
         if (  ScintHitStation[i]==j  ){
            
	        t+=1;

	         Z+=(ScintHitPosZ[i]);
	         Y+=(ScintHitPosY[i]);
	         X+=(ScintHitPosX[i]);
       	   H+=HitPDG[i];
	         E+=ScintHitE[i];
		 //  cout<<"EventID "<<ScintHitPrimaryID[i]<<endl;
                   
	            }
    }
	       }   
      }      
            
     
    
      if(t!=0){
      double Zavg= Z/t;
      double Yavg= Y/t;
      double Xavg= X/t;
      double Havg= H/t;//caution: loosly defined
      
      //  cout<<"X ="<<X<<endl;
      //	cout<<"t ="<<t<<endl;
     	//cout<<Xavg<<endl;
      HlistZ[j]->Fill(Zavg);
     
      HlistX[j]->Fill(Xavg);
        
      HlistY[j]->Fill(Yavg);

      Sx[j].push_back(Xavg);
      Sy[j].push_back(Yavg);
      Sz[j].push_back(Zavg);
     // cout<<"sizex = "<<Sx[j].size()<<"sizey = "<<Sy[j].size()<<"sizez = "<<Sz[j].size() <<endl; 
      // HlistPDG[j]->Fill(Havg);
      HlistE[j]->Fill(E);
      
      HE[j].push_back(E);
              }
        }
       }

      
      histGenE->Fill(GenPartE);
      histGenTheta->Fill(GenPartTheta);

      

   }


for(int q =0; q<Sx[1].size(); q++){ 
 //Vector3D v1(Sx[1][q]-Sx[0][q],Sy[1][q]-Sy[0][q],Sz[1][q]-Sz[0][q]);
 //Vector3D v2(Sx[3][q]-Sx[2][q],Sy[3][q]-Sy[2][q],Sz[3][q]-Sz[2][q]);

TVector3 v1(Sx[1][q]-Sx[0][q],Sy[1][q]-Sy[0][q],Sz[1][q]-Sz[0][q]);
TVector3 v2(Sx[3][q]-Sx[2][q],Sy[3][q]-Sy[2][q],Sz[3][q]-Sz[2][q]);
   double TAngle= v1.Angle(v2);
   TVector3 CP = v1.Cross(v2);
   //xz plane
TVector3 v1xz(Sx[1][q]-Sx[0][q],0.,Sz[1][q]-Sz[0][q]);
TVector3 v2xz(Sx[3][q]-Sx[2][q],0.,Sz[3][q]-Sz[2][q]);
//yz plane
TVector3 v1yz(0,Sy[1][q]-Sy[0][q],Sz[1][q]-Sz[0][q]);
TVector3 v2yz(0,Sy[3][q]-Sy[2][q],Sz[3][q]-Sz[2][q]);

TVector3 CPxz = v1xz.Cross(v2xz);
TVector3 CPyz = v1yz.Cross(v2yz);

   double TAnglexz = v1xz.Angle(v2xz);
   double TAngleyz = v1yz.Angle(v2yz);

   int sign;
    if(CPxz.Y()<0){sign=-1;}
    if(CPxz.Y()>0){sign=1;}
   
 
   //double angle = calculateAngle(v1, v2);
   //double angleDegrees = angle * 180.0 / M_PI;
   
   histAngle->Fill(sign*TAnglexz);
   }



TF1* fitFunc = new TF1("fitFunc", MyFitFunction, -2, 4,2);
      Double_t max = histGenTheta->GetMaximum();
      std::cout<<max<<endl;
     fitFunc->SetParameters(0,1000);
     fitFunc->SetParameters(1,2.5);
     c8->cd();
    histGenTheta->Draw();
    histGenTheta->Fit("fitFunc","R");


for(int j =0; j<4; j++){ 
   for(int i=0; i<HE[j].size(); i++){
      outFile<<HE[j][i]<<",";
   }
   outFile<<endl;
}


      outFile.close();


      for (int i=0; i<4; i++) {
	  
       c1->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistE[i]->Draw();

        c2->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistZ[i]->Draw();

        c3->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistY[i]->Draw();

        c4->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistX[i]->Draw();
       }

       c5->cd(1);
       histGenE->Draw();
       c5->cd(2);
      // histGenTheta->Draw();

       c6->cd();
       histAngle->Draw();
       
}
