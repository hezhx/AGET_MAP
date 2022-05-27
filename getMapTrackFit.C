#include "vector_define.h"
#include "TrackFit.h"
#include "getClusterNumber.h"
#include "TTree.h"
#include "TFile.h"

int getMapTrackFit(int id){
// int getMapTrackFit(){
    Double_t slope;
    Double_t deltay;

    TFile *ff = new TFile("getFitResult.root", "RECREATE");
    TTree *TT = new TTree("TT", "Fit Result");

    TT->Branch("slope", &slope, "slope/D");
    TT->Branch("deltay", &deltay, "deltay/D");

    // TH1F *hd = new TH1F("hd", "y position resolution", 1024, -5, 5);
    Double_t *delta_y;
    int k = id;
    // int nentries = 41206; // for N750atm.root
    int nentries = 33917; // for N850atm.root
    // int nentries = 3014; // for N850atm.root
    while (k<nentries){
        // cout << "------\n"
        //      << getClusterNumber(k) 
        //      << " : Next event starts. \n------\n";
        cout << getClusterNumber(k) << endl;
        delta_y = TrackFit(k);
        // cout << (int)delta_y[0] << endl;
        if (delta_y[0]>0){
            for (int i=2; i<=(int)delta_y[0]; i++){
                // cout << delta_y[1] << " " << delta_y[i] << endl;
                slope = delta_y[1];
                deltay = delta_y[i];
                TT->Fill();
                // hd->Fill(delta_y[i]);
            }            
        }
        k = getClusterNumber(k);
    }

    // TCanvas *c2 = new TCanvas("c2", "y position", 600, 600);
    // c2->cd();
    // hd->Draw();
    // c2->SaveAs("deltay.root");
    TT->Write();
    ff->Close();


    return 1;
}
