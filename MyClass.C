#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MyClass::Loop()
{

    fstream outFile("test11.csv",ios::out);
  /*TCanvas* canvas1 = new TCanvas("canvas1", "2D Histogram");
 canvas1->Divide(1,2);
TH2F* histogram1 = new TH2F("histogram1", "2D Histogram", 100, -250, 250, 100, -250, 250);

//TCanvas* canvas2 = new TCanvas("canvas2", "2D Histogram");
TH2F* histogram2 = new TH2F("histogram2", "2D Histogram", 100, -250, 250, 100, -250, 250); 
  */
  
  char nameE[15], titleE[20];
   TH1F *HlistE[4];      // create an array of Histograms
   TH1F* hStationNumberE;                 // create a pointer to a histogram
   // make  8 histograms and add them to the object array
   for (Int_t i = 0; i < 4; i++) {
      sprintf(nameE,"hStationNumberE%d",i);
      sprintf(titleE,"station:%d",i);
      hStationNumberE = new TH1F(nameE,titleE,100,0,0.02);
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
      hStationNumberX = new TH1F(nameX,titleX,100,-250,250);
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
      hStationNumberY = new TH1F(nameY,titleY,100,-250,250);
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


    TCanvas *c1 = new TCanvas("c1","Histograms related Planes posX");
        c1->Divide(2,2);
      TCanvas *c2 = new TCanvas("c2", "Histograms related to the Planes posY ");
       c2->Divide(2,2);
     TCanvas *c3 = new TCanvas("c3", "Histograms related to  the Planes posZ ");
       c3->Divide(2,2);
       TCanvas *c4 = new TCanvas("c4", "Histograms related to  the Planes dep energy ");
       c4->Divide(2,2);
   
  
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //	cout<<"x-pos  "<< X/t<<endl;
      //	cout<<"eventid   "<< EventId <<endl;
      if(abs(PDG)==13){
	double X =0.;
	double Z =0.;
	double t=0.;
	
	/*	for (int i=0;  i<= EventId; i++ )
	  { X=+x;
	    Z=+z;
	    t=+1.;      }
	//cout<<"t  "<< t<<endl;
		if(t!=0){
		  //	cout<<"x-pos  "<< X/t<<endl;
	cout<<"z-pos  "<< Z/t<<endl;
	}   */
	       
	 for(int i=0; i<4; i++){
		if(Plane==i ){
		
	
		  	HlistX[i]->Fill(x);
		  HlistY[i]->Fill(y);
		  HlistZ[i]->Fill(z);
		  HlistE[i]->Fill(edep);
		  	 
	outFile<< x<<",";
	outFile<< y<<",";
   	outFile<< z<<",";

	//	outFile<<endl;
		}
		
		if(Plane==1){
	outFile<< x<<",";
	outFile<< y<<",";
   	outFile<< z<<",";
	//	 outFile<<endl;
	
		}

	  	if(Plane==2){
	outFile<< x<<",";
	outFile<< y<<",";
   	outFile<< z<<",";
	//	outFile<<endl;
	//	cout<<z<<endl;
		}
      
	
	  	if(Plane==3){
	outFile<< x<<",";
	outFile<< y<<",";
   	outFile<< z<<",";
	//	 cout<<z<<endl;
	outFile<<endl;
          }
		
	
		 }
	
	   } 

	 } //if pdgid 13
     
      
      
      
    outFile.close();

    for (int i=0; i<4; i++) {
	  
       c1->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistX[i]->Draw();

        c2->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistY[i]->Draw();

        c3->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistZ[i]->Draw();

        c4->cd(i+1);
       gPad->DrawFrame(0,0,1,1);
       HlistE[i]->Draw();
       }
   /*  canvas1->cd(1);
   gPad->DrawFrame(0,0,1,1);
   histogram1->Draw("colz");
   
   canvas1->cd(2);
   gPad->DrawFrame(0,0,1,1);
   histogram2->Draw("colz");
   */
   
}
