#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include "TVector.h"
#include "TMath.h"

/*struct TVector3
{
   double x, y, z;
};
*/

// Function to find the mutual closest point between two skew vectors
TVector3 findClosestPoint(const TVector3 &p1, const TVector3 &v1, const TVector3 &p2, const TVector3 &v2)
{
   TVector3 closestPoint;
   TVector3 closestPoint1;
   TVector3 closestPoint2;

   // Calculate the direction vectors of the two lines
   TVector3 dir1(v1.X() - p1.X(), v1.Y() - p1.Y(), v1.Z() - p1.Z());
   TVector3 dir2(v2.X() - p2.X(), v2.Y() - p2.Y(), v2.Z() - p2.Z());

   // Calculate the dot products
   /* double dot1 = dir1.x * dir1.x + dir1.y * dir1.y + dir1.z * dir1.z;
    double dot2 = dir2.x * dir2.x + dir2.y * dir2.y + dir2.z * dir2.z;
    double dotProduct = dir1.x * dir2.x + dir1.y * dir2.y + dir1.z * dir2.z;
 */

   double dot1 = dir1.Dot(dir1);
   double dot2 = dir2.Dot(dir2);
   double dotProduct = dir1.Dot(dir2);

   // Calculate the parameters t1 and t2
   double t1 = (dotProduct * dot2 - dot1) / (dot1 * dot2 - dotProduct * dotProduct);
   double t2 = (dotProduct * dot1 - dot2) / (dot1 * dot2 - dotProduct * dotProduct);

   // Calculate the closest points on the two lines
   closestPoint1.SetX(p1.X() + t1 * dir1.X());
   closestPoint1.SetY(p1.Y() + t1 * dir1.Y());
   closestPoint1.SetZ(p1.Z() + t1 * dir1.Z());

   closestPoint2.SetX(p2.X() + t2 * dir2.X());
   closestPoint2.SetY(p2.Y() + t2 * dir2.Y());
   closestPoint2.SetZ(p2.Z() + t2 * dir2.Z());

   TVector3 cpVector(closestPoint1.X() - closestPoint2.X(), closestPoint1.Y() - closestPoint2.Y(), closestPoint1.Z() - closestPoint2.Z());

   double distance = cpVector.Mag();

   if (distance > 0.)
   {
      closestPoint.SetX((closestPoint1.X() - closestPoint2.X()) / 2.);
      closestPoint.SetY((closestPoint1.Y() - closestPoint2.Y()) / 2.);
      closestPoint.SetZ((closestPoint1.Z() - closestPoint2.Z()) / 2.);
   }
   return closestPoint;
}

Double_t MyFitFunction(Double_t *x, Double_t *par)
{
   Double_t value = par[0] * TMath::Power(TMath::Cos(x[0]), par[1]); // cos^2(x)

   return value;
}

void MyClass::Loop()
{
   vector<float> S0;
   vector<float> S1;
   vector<float> S2;
   vector<float> S3;
   vector<float> Sx[4] = {S0, S1, S2, S3};
   vector<float> Sy[4] = {S0, S1, S2, S3};
   vector<float> Sz[4] = {S0, S1, S2, S3};

   vector<float> Planenb;

   fstream outFile("test11.csv", ios::out);
   vector<float> HE0;
   vector<float> HE1;
   vector<float> HE2;
   vector<float> HE3;
   vector<float> HE[4] = {HE0, HE1, HE2, HE3};

   TH1F *histAngle = new TH1F("histAngle", "scattering angle", 200, 0., 0.5);

   TH1F *histGenE = new TH1F("histGenE", "KE", 100, -100, 14000);
   histGenE->GetXaxis()->SetTitle("Kinetic Energy[MeV]");
   histGenE->GetYaxis()->SetTitle("N");
   TH1F *histGenTheta = new TH1F("histGenTheta", "theta", 100, -1, 1.6);
   histGenTheta->GetXaxis()->SetTitle("Radian");
   histGenTheta->GetYaxis()->SetTitle("N");

   TH1F *histGenPhi = new TH1F("histGenPhi", "Phi", 100, -3.2, 3.2);
   histGenPhi->GetXaxis()->SetTitle("Radian");
   histGenPhi->GetYaxis()->SetTitle("N");

   char nameE[15], titleE[20];
   TH1F *HlistE[4];       // create an array of Histograms
   TH1F *hStationNumberE; // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++)
   {
      sprintf(nameE, "hStationNumberE%d", i);
      sprintf(titleE, "station:%d", i);
      hStationNumberE = new TH1F(nameE, titleE, 100, -0.0001, 0.0030);
      hStationNumberE->GetXaxis()->SetTitle("Deposited Energy[MeV]");
      hStationNumberE->GetYaxis()->SetTitle("numHits");
      HlistE[i] = hStationNumberE;
   }

   char nameX[15], titleX[20];
   TH1F *HlistX[4];       // create an array of Histograms
   TH1F *hStationNumberX; // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++)
   {
      sprintf(nameX, "hStationNumberX%d", i);
      sprintf(titleX, "station:%d", i);
      hStationNumberX = new TH1F(nameX, titleX, 100, -260, 260);
      hStationNumberX->GetXaxis()->SetTitle("mm");
      hStationNumberX->GetYaxis()->SetTitle("numHits");
      HlistX[i] = hStationNumberX;
   }

   char nameY[15], titleY[20];
   TH1F *HlistY[10];      // create an array of Histograms
   TH1F *hStationNumberY; // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++)
   {
      sprintf(nameY, "hStationNumberY%d", i);
      sprintf(titleY, "station:%d", i);
      hStationNumberY = new TH1F(nameY, titleY, 100, -260, 260);
      hStationNumberY->GetXaxis()->SetTitle("mm");
      hStationNumberY->GetYaxis()->SetTitle("numHits");
      HlistY[i] = hStationNumberY;
   }

   char nameZ[15], titleZ[20];
   TH1F *HlistZ[10];      // create an array of Histograms
   TH1F *hStationNumberZ; // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++)
   {
      sprintf(nameZ, "hStationNumberZ%d", i);
      sprintf(titleZ, "station:%d", i);
      hStationNumberZ = new TH1F(nameZ, titleZ, 100, -10, 350);
      hStationNumberZ->GetXaxis()->SetTitle("mm");
      hStationNumberZ->GetYaxis()->SetTitle("numHits");
      HlistZ[i] = hStationNumberZ;
   }


 char nameP[15], titleP[20];
   TH2F *HlistP[10];      // create an array of Histograms
   TH2F *hStationNumberP; // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++)
   {
      sprintf(nameP, "hStationNumberP%d", i);
      sprintf(titleP, "station:%d", i);
      hStationNumberP = new TH2F(nameP, titleP, 100 ,-250, 250,100,-250,250);
     // hStationNumberP->GetXaxis()->SetTitle("PDG");
     // hStationNumberP->GetYaxis()->SetTitle("numHits");
      HlistP[i] = hStationNumberP;
   }
   TCanvas *canvas1 = new TCanvas("canvas1", "2D Histogram");
   canvas1->Divide(1, 2);
   TH2F *histogram1 = new TH2F("histogram1", "2D Histogram", 100, -50, 50, 100, -50, 50);
   TH2F *histogram2 = new TH2F("histogram2", "2D Histogram", 100, -250, 250, 100, -250, 250);

   TCanvas *c1 = new TCanvas("c1", "Histograms related Planes edep ");
   c1->Divide(2, 2);
   TCanvas *c2 = new TCanvas("c2", "Histograms related to the Planes posZ ");
   c2->Divide(2, 2);
   TCanvas *c3 = new TCanvas("c3", "Histograms related to  the Planes posY ");
   c3->Divide(2, 2);
   TCanvas *c4 = new TCanvas("c4", "Histograms related to  the Planes posX");
   c4->Divide(2, 2);

   TCanvas *c5 = new TCanvas("c5", "Histograms related Generator");
   c5->Divide(2, 2);

   TCanvas *c7 = new TCanvas("c5", "Histograms related PDG");
   c7->Divide(2, 2);


   TCanvas *c6 = new TCanvas("c6", "Histograms related to Scattering angle");
   TCanvas *c8 = new TCanvas("c8", "Histograms related fit");

   int t = 0;

   if (fChain == 0)
      return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Planenb.clear();
      if (NScintHit != 0)
      {
         for (int i = 0; i < NScintHit; i++)
         {
            if (abs(HitPDG[i]) == 13)
            {
               //  cout<<"first loop event "<<t<<endl;
               Planenb.push_back(ScintHitStation[i]);
            }
         }
      }
      auto it = std::find(Planenb.begin(), Planenb.end(), 3);
      if (it != Planenb.end())
      {
         for (int j = 0; j < 4; j++)
         {
            double E = 0.;
            double Z = 0.;
            double Y = 0.;
            double X = 0.;
            double t = 0;
            double H = 0.;

            if (NScintHit != 0)
            {

               for (int i = 0; i < NScintHit; i++)
               {
                  if (abs(HitPDG[i]) == 13)
                  {

                     //  cout<<"second loop event "<<t<<endl;
                     // Planenb.push_back(ScintHitStation[i]);

                     if (ScintHitStation[i] == j)
                     {

                        t += 1;

                        Z += (ScintHitPosZ[i]);
                        Y += (ScintHitPosY[i]);
                        X += (ScintHitPosX[i]);
                        H += HitPDG[i];
                        E += ScintHitE[i];
                       
                         
                    
                        //  cout<<"EventID "<<ScintHitPrimaryID[i]<<endl;
                     }
                  }
               }
            }
            if (t != 0)
            {
               double Zavg = Z / t;
               double Yavg = Y / t;
               double Xavg = X / t;
               double Havg = H / t; // caution: loosly defined
               double Eavg = E;
              
               //  cout<<"X ="<<X<<endl;
               //	cout<<"t ="<<t<<endl;
               // cout<<Xavg<<endl;
               HlistZ[j]->Fill(Zavg);

               HlistX[j]->Fill(Xavg);

               HlistY[j]->Fill(Yavg);
             // histogram2->Fill(Xavg, Yavg);
               cout<<Xavg<<endl;
               Sx[j].push_back(Xavg);
               Sy[j].push_back(Yavg);
               Sz[j].push_back(Zavg);
               // cout<<"sizex = "<<Sx[j].size()<<"sizey = "<<Sy[j].size()<<"sizez = "<<Sz[j].size() <<endl;
               HlistP[j]->Fill(Xavg,Yavg);
               HlistE[j]->Fill(E);
               HE[j].push_back(E);
            }
         }
      }
      histGenE->Fill(GenPartE);
      histGenTheta->Fill(GenPartTheta);
      histGenPhi->Fill(GenPartPhi);
   }

   for (int q = 0; q < Sx[1].size(); q++)
   {
      // Vector3D v1(Sx[1][q]-Sx[0][q],Sy[1][q]-Sy[0][q],Sz[1][q]-Sz[0][q]);
      // Vector3D v2(Sx[3][q]-Sx[2][q],Sy[3][q]-Sy[2][q],Sz[3][q]-Sz[2][q]);

      TVector3 v1(Sx[1][q] - Sx[0][q], Sy[1][q] - Sy[0][q], Sz[1][q] - Sz[0][q]);
      TVector3 v2(Sx[3][q] - Sx[2][q], Sy[3][q] - Sy[2][q], Sz[3][q] - Sz[2][q]);
      double TAngle = v1.Angle(v2);
      TVector3 CP = v1.Cross(v2);
      // xz plane
      TVector3 v1xz(Sx[1][q] - Sx[0][q], 0., Sz[1][q] - Sz[0][q]);
      TVector3 v2xz(Sx[3][q] - Sx[2][q], 0., Sz[3][q] - Sz[2][q]);
      // yz plane
      TVector3 v1yz(0, Sy[1][q] - Sy[0][q], Sz[1][q] - Sz[0][q]);
      TVector3 v2yz(0, Sy[3][q] - Sy[2][q], Sz[3][q] - Sz[2][q]);

      TVector3 CPxz = v1xz.Cross(v2xz);
      TVector3 CPyz = v1yz.Cross(v2yz);

      // TVector3 test1(0.,0.,1);
      // TVector3 test2(0.,-0.7,-0.7);
      //  cout<<test2.Angle(test1)<<endl;

      double TAnglexz = v1xz.Angle(v2xz);
      double TAngleyz = v1yz.Angle(v2yz);

      int signxz;
      if (CPxz.Y() < 0)
      {
         signxz = -1;
      }
      if (CPxz.Y() > 0)
      {
         signxz = 1;
      }

      int signyz;
      if (CPxz.X() < 0)
      {
         signyz = -1;
      }
      if (CPxz.X() > 0)
      {
         signyz = 1;
      }
      int sign = signxz * signyz;

      // double angle = calculateAngle(v1, v2);
      // double angleDegrees = angle * 180.0 / M_PI;
      double TTAngle = sqrt(TAnglexz * TAnglexz + TAngleyz * TAngleyz);
      histAngle->Fill(TTAngle);
      if (v1.Dot(v2) != 0 && CP.Mag() != 0)
      {
         TVector3 p1(Sx[0][q], Sy[0][q], Sz[0][q]);
         TVector3 V1(Sx[1][q] - Sx[0][q], Sy[1][q] - Sy[0][q], Sz[1][q] - Sz[0][q]);
         TVector3 p2(Sx[2][q], Sy[2][q], Sz[2][q]);
         TVector3 V2(Sx[3][q] - Sx[2][q], Sy[3][q] - Sy[2][q], Sz[3][q] - Sz[2][q]);
         TVector3 poca = findClosestPoint(p1, V1, p2, V2);

        // cout << poca.X() << endl;
         histogram1->Fill(1.*poca.X(),1.*poca.Y());
      }
   }
   TF1 *fitFunc = new TF1("fitFunc", MyFitFunction, -0, 1.5, 2);
   Double_t max = histGenTheta->GetMaximum();
  // std::cout << max << endl;
   fitFunc->SetParameter(0, 1842);
   fitFunc->SetParameter(1, 2.025);
   c8->cd();
   fitFunc->Draw();
   // histGenTheta->Fit("fitFunc","","",0.0,1.4);
   // histGenTheta->Draw("SAME");

   // fitFunc->Draw("SAME");

   for (int j = 0; j < 4; j++)
   {
      for (int i = 0; i < HE[j].size(); i++)
      {
         outFile << HE[j][i] << ",";
      }
      outFile << endl;
   }

   outFile.close();

   for (int i = 0; i < 4; i++)
   {

      c1->cd(i + 1);
      gPad->DrawFrame(0, 0, 1, 1);
      HlistE[i]->Draw();

      c2->cd(i + 1);
      gPad->DrawFrame(0, 0, 1, 1);
      HlistZ[i]->Draw();

      c3->cd(i + 1);
      gPad->DrawFrame(0, 0, 1, 1);
      HlistY[i]->Draw();

      c4->cd(i + 1);
      gPad->DrawFrame(0, 0, 1, 1);
      HlistX[i]->Draw();

      c7->cd(i + 1);
      gPad->DrawFrame(0, 0, 1, 1);
      HlistP[i]->Draw("colz");
   }

   c5->cd(1);
  // histGenE->Draw();
   c5->cd(2);
 //  histGenTheta->Draw();
   c5->cd(3);
   //histGenPhi->Draw();
   c6->cd();
   histAngle->Draw();

   canvas1->cd(1);
   histogram1->Draw("colz");
   canvas1->cd();
  //histogram2->Draw("LEGO2");
}
