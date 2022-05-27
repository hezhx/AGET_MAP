#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>

// #include "vector_define.h"

#include "TFile.h"
// #include "TrackFit.h"

using namespace std;

struct Reading{
    int pid;
    int bid;
    int cid;
    int chid;

    double yc;
    double zc;

    int rid;

    double yl;
    double zl;
};

struct ReadingRoot{
    int eid;
    int bid;
    int cid;
    int chid;

    double lev;
    double riv;
    double pep;
    double amp;
};

struct FitData{
    double yc;
    double zc;
    double amp;
    int rid;
};

// int TrackFit(int ev_ID);
int getClusterNumber(int ev_ID){
    // cout << "Please enter map file : ";
    // string mapfilename;
    // cin >> mapfilename;
    // ifstream mapfile{mapfilename};
    ifstream mapfile{"cert.map"};
    // if (!mapfile) error("Can't open map file ", mapfilename);

    vector<Reading> mapvec;
    int pid;
    int bid;
    int cid;
    int chid;

    double yc;
    double zc;

    int rid;

    double yl;
    double zl;

    double amp = 0.;

    while(mapfile >> pid >> bid >> cid >> chid >> yc >> zc >> rid >> yl >> zl){
        mapvec.push_back(Reading{pid, bid, cid, chid, yc, zc, rid, yl, zl});
    }


    int pad_id = 0;
    int board_id = 0;
    int chip_id = 0;
    int channel_id = 0;

    double center_y = 0.;
    double center_z = 0.;

    int row_id = 0;

    double length_y = 0.;
    double length_z = 0.;

    // TH2Poly *TH2_Poly = new TH2Poly();

    for (int i=0; i<mapvec.size(); ++i){
        // cout << mapvec[i].pid << " " << mapvec[i].zl << endl;
        pad_id = mapvec[i].pid;
        board_id = mapvec[i].bid;
        chip_id = mapvec[i].cid;
        channel_id = mapvec[i].chid;

        center_y = mapvec[i].yc;
        center_z = mapvec[i].zc;
        TVector3 center = TVector3(0., center_y, center_z);

        row_id = mapvec[i].rid;

        length_y = mapvec[i].yl;
        length_z = mapvec[i].zl;
        TVector3 length = TVector3(0., length_y, length_z);

        // int pos_num = 20;
        // double y[20] = {0.};
        // double z[20] = {0.};

        // // cout << " test1" << endl;
        // y[0] = center.Y() + length_y*0.5;
        // z[0] = center.Z() - length_z*8./16.;

        // y[1] = center.Y() + length_y*0.;
        // z[1] = center.Z() - length_z*7./16.;

        // y[2] = center.Y() + length_y*1.;
        // z[2] = center.Z() - length_z*5./16.;

        // y[3] = center.Y() + length_y*0.;
        // z[3] = center.Z() - length_z*3./16.;

        // y[4] = center.Y() + length_y*1.;
        // z[4] = center.Z() - length_z*1./16.;

        // y[5] = center.Y() + length_y*0.;
        // z[5] = center.Z() + length_z*1./16.;

        // y[6] = center.Y() + length_y*1.;
        // z[6] = center.Z() + length_z*3./16.;

        // y[7] = center.Y() + length_y*0.;
        // z[7] = center.Z() + length_z*5./16.;

        // y[8] = center.Y() + length_y*1.;
        // z[8] = center.Z() + length_z*7./16.;

        // y[9] = center.Y() + length_y*0.5;
        // z[9] = center.Z() + length_z*8./16.;

        // y[10] = center.Y() - length_y*0.5;
        // z[10] = center.Z() + length_z*8./16.;

        // y[11] = center.Y() - length_y*0.;
        // z[11] = center.Z() + length_z*7./16.;

        // y[12] = center.Y() - length_y*1.;
        // z[12] = center.Z() + length_z*5./16.;

        // y[13] = center.Y() - length_y*0.;
        // z[13] = center.Z() + length_z*3./16.;

        // y[14] = center.Y() - length_y*1.;
        // z[14] = center.Z() + length_z*1./16.;

        // y[15] = center.Y() - length_y*0.;
        // z[15] = center.Z() - length_z*1./16.;

        // y[16] = center.Y() - length_y*1.;
        // z[16] = center.Z() - length_z*3./16.;

        // y[17] = center.Y() - length_y*0.;
        // z[17] = center.Z() - length_z*5./16.;

        // y[18] = center.Y() - length_y*1.;
        // z[18] = center.Z() - length_z*7./16.;

        // y[19] = center.Y() - length_y*0.5;
        // z[19] = center.Z() - length_z*8./16.;
        // cout << " test1" << endl;

        // if (board_id==1 && chip_id==1){
        //     amp = 1.;
        // }else if (board_id==1 && chip_id==2){
        //     amp = 2.;
        // }else if (board_id==1 && chip_id==0){
        //     amp = 3.;
        // }else if (board_id==1 && chip_id==3){
        //     amp = 4.;
        // }else if (board_id==2 && chip_id==0){
        //     amp = 5.;
        // }else if (board_id==2 && chip_id==1){
        //     amp = 6.;
        // }else if (board_id==2 && chip_id==2){
        //     amp = 7.;
        // }else if (board_id==3 && chip_id==0){
        //     amp = 8.;
        // }else if (board_id==3 && chip_id==1){
        //     amp = 9.;
        // }else if (board_id==3 && chip_id==2){
        //     amp = 10.;
        // }else if (board_id==2 && chip_id==3){
        //     amp = 11.;
        // }else if (board_id==4 && chip_id==1){
        //     amp = 12.;
        // }else if (board_id==4 && chip_id==2){
        //     amp = 13.;
        // }else if (board_id==3 && chip_id==3){
        //     amp = 14.;
        // }

        // TH2_Poly->AddBin(pos_num, y, z);

        // TH2_Poly->Fill(center.Y(), center.Z(), gRandom->Gaus(80,500));
        // TH2_Poly->Fill(center.Y(), center.Z(), amp);
    }
    TFile *f1 = new TFile("N850atm.root");

    Int_t eventID;
    Int_t boardID;
    Int_t chipID;
    Int_t channelID;
    Double_t leftVal;
    Double_t rightVal;
    Double_t peakPos;
    Double_t ampVal;

    TTree *T1;
    T1 = (TTree *)f1->Get("T1");
    T1->SetBranchAddress("eventID", &eventID);
    T1->SetBranchAddress("boardID", &boardID);
    T1->SetBranchAddress("chipID", &chipID);
    T1->SetBranchAddress("channelID", &channelID);

    T1->SetBranchAddress("leftVal", &leftVal);
    T1->SetBranchAddress("rightVal", &rightVal);
    T1->SetBranchAddress("peakPos", &peakPos);
    T1->SetBranchAddress("ampVal", &ampVal);

    vector<ReadingRoot> rootvec;
    vector<FitData> fitdata;

    Int_t nentries = T1->GetEntries();
    nentries = 10000;
    for (int i=0; i<nentries; i++){
        T1->GetEntry(i);
        rootvec.push_back(ReadingRoot{eventID, boardID, chipID, channelID, leftVal, rightVal, peakPos, ampVal});
    }

    bool flag = true;
    vector<double> yc_CM;
    vector<double> zc_CM;
    double yc_cm = 0.;
    double zc_cm = 0.;
    double yc_cm_tot = 0.;
    double amp_tot = 0.;
    int r_max = 0;
    int r_min = 0;
    vector<int> row_all; 
    
    // 0, 17, 36, 53, 69
    for (int j=ev_ID; j<nentries && flag==true; j++){
        // tot_energy = 0.;
        // cout << rootvec[j].eid << " " << rootvec[j].amp << endl;
        // nevent = 0;
        if (rootvec[j].eid == rootvec[j+1].eid){
            for (int k=0; k<mapvec.size(); ++k){
                    if (mapvec[k].bid==rootvec[j].bid && mapvec[k].cid==rootvec[j].cid && mapvec[k].chid==rootvec[j].chid){
                        // TH2_Poly->Fill(mapvec[k].yc, mapvec[k].zc, rootvec[j].amp);
                        // cout << rootvec[j].eid << " " 
                        //      << mapvec[k].yc << " " << mapvec[k].zc << " "
                        //      << mapvec[k].rid << " "
                        //      << " " << rootvec[j].amp << "\n";
                        fitdata.push_back(FitData{mapvec[k].yc, mapvec[k].zc, rootvec[j].amp, mapvec[k].rid});
                        row_all.push_back(mapvec[k].rid);

                    }
            }
            // cout << rootvec[j].eid << " " << rootvec[j].amp << endl;            
        }else{
            for (int k=0; k<mapvec.size(); ++k){
                    if (mapvec[k].bid==rootvec[j].bid && mapvec[k].cid==rootvec[j].cid && mapvec[k].chid==rootvec[j].chid){
                        // TH2_Poly->Fill(mapvec[k].yc, mapvec[k].zc, rootvec[j].amp);
                        // cout << rootvec[j].eid << " " 
                        //      << mapvec[k].yc << " " << mapvec[k].zc << " "
                        //      << mapvec[k].rid << " "
                        //      << " " << rootvec[j].amp << "\n";  
                        fitdata.push_back(FitData{mapvec[k].yc, mapvec[k].zc, rootvec[j].amp, mapvec[k].rid});
                        row_all.push_back(mapvec[k].rid);
                    }
            }
            // cout << rootvec[j].eid << " " << tot_energy << "\n";
            // h1->Fill(tot_energy*4096./50000.);
            // tot_energy = 0.;
            yc_cm = 0.;
            zc_cm = 0.;
            yc_cm_tot = 0.;
            amp_tot = 0.;
            flag = false;
            // cout << fitdata[1].yc << " " << fitdata[5].amp << endl;
            // fitdata is a complete track
            // cout << *min_element(row_all.begin(), row_all.end()) << endl;
            r_max = *max_element(row_all.begin(), row_all.end());
            r_min = *min_element(row_all.begin(), row_all.end());

            for (int r=r_min; r<=r_max; r++){
                for (int g=0; g<fitdata.size(); ++g){
                    if (fitdata[g].rid==r){
                        zc_cm = fitdata[g].zc;
                        yc_cm_tot = fitdata[g].yc * fitdata[g].amp + yc_cm_tot;
                        amp_tot = fitdata[g].amp + amp_tot;
                    }
                }
                yc_cm = yc_cm_tot / amp_tot;
                yc_CM.push_back(yc_cm);
                zc_CM.push_back(zc_cm);
                // cout << yc_cm << " " << zc_cm << endl;
            }
            // cout << j+1 << endl;
            return j+1;
        }
    }
    // TCanvas *c1 = new TCanvas("c1", "testMap", 1000., 600.);
    // c1->Divide(2,1);

    // int n = zc_CM.size();
    // int Position_CM[2][n];
    // for (int i=0; i<n; i++){
    //     Position_CM[0][i] = zc_CM[i];
    //     Position_CM[1][i] = yc_CM[i];
    // }
    // TGraph *gtrack = new TGraph(n, Position_CM[0], Position_CM[1]);
    // TGraph *gtrack2 = new TGraph(n, Position_CM[1], Position_CM[0]);
    // TF1 *linearFit = new TF1("linearFit", "pol1");
    // gtrack->Fit("linearFit","Q");
    // TF1 *fitFunc = gtrack->GetFunction("linearFit");
    // double p0 = fitFunc->GetParameter(0);
    // double p1 = fitFunc->GetParameter(1);
    // cout << "Yc = Zc*(" << p1 << ") + " << p0 << endl;

    // double delta_yc;
    // for (int i=0; i<n; i++){
    //     delta_yc = Position_CM[1][i] - (Position_CM[0][i] * p1 + p0);
    //     cout << delta_yc << endl;
    // }

    // gtrack->SetLineWidth(0);
    // gtrack->SetMarkerSize(5);
    // gtrack->SetMarkerColor(9);
    // gtrack->SetMarkerStyle(22);

    // c1->cd(1);
    // gStyle->SetOptStat(0);
    // TH2_Poly->SetTitle("testMap");

    // TH2_Poly->Draw("colz");
    // c1->cd(2);
    // gtrack2->SetLineWidth(0);
    // gtrack2->SetMarkerSize(5);
    // gtrack2->SetMarkerColor(9);
    // gtrack2->SetMarkerStyle(22);
    // gtrack2->Draw();
    // c1->SaveAs("testMap.png");
    // c1->SaveAs("testMap.root");


    return 0;
}
