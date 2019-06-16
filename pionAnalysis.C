#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"

using namespace std;

void pionAnalysis(char* filename){
    // Read all files
    ifstream input(filename);
    
    auto outName = string(filename).substr(0, 4);
    
    auto chain = new TChain("tree");

    string file;
    while (input >> file) chain->Add(Form("%s", file.c_str()));

    int nevents = chain->GetEntries();

    Float_t t_p3pi[3], t_piKE, t_piGArLen, t_pionExitMom[3];
    Int_t PDG, t_fsPdg[100], t_nFS;
    Float_t t_fsPx[100],t_fsPy[100],t_fsPz[100];
    Float_t t_vtx[3];

    chain->SetBranchAddress("p3pi", &t_p3pi);
    chain->SetBranchAddress("piKE", &t_piKE);
    chain->SetBranchAddress("piPdg", &PDG);
    chain->SetBranchAddress("piGArLen", &t_piGArLen);
    chain->SetBranchAddress("pionExitMom", &t_pionExitMom);

    chain->SetBranchAddress("fsPdg", &t_fsPdg);
    chain->SetBranchAddress("nFS", &t_nFS);
    chain->SetBranchAddress("fsPx", &t_fsPx);
    chain->SetBranchAddress("fsPy", &t_fsPy);
    chain->SetBranchAddress("fsPz", &t_fsPz);
    chain->SetBranchAddress("vtx", &t_vtx);

    auto hist_TrueInformation = new TH2F("hist_TrueInformation", 
					 "Pion True Information; Energy [GeV]; Angle [deg]",
					 60, 0, 20, 180, 0, 180);
    auto hist_1 = new TH2F("hist_1", "Energy within 3DST; Energy [GeV]; Angle",
					 60, 0, 20, 180, 0, 180);
    auto hist_2 = new TH2F("hist_2", "Energy within 3DST + TPC; Energy [GeV]; Angle]",
					 60, 0, 20, 180, 0, 180);
    auto hist_3 = new TH2F("hist_3", "Energy within 3DST + TPC + ECAL; Energy [GeV]; Angle]",
					 60, 0, 20, 180, 0, 180);

    int CC_event_count = 0;
    for (int i = 0; i < nevents; ++i){
	chain->GetEntry(i);

	// TODO: Change this according to the geometry
	if (abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50){

	    for (int index = 0; index < t_nFS; index++){
		// looking for CC events
		if (abs(t_fsPdg[index]) == 11 || abs(t_fsPdg[index]) == 13 ){
		    CC_event_count++;
		    if (abs(PDG) == 211){
			float p2 = t_p3pi[0] * t_p3pi[0] + 
				   t_p3pi[1] * t_p3pi[1] + 
				   t_p3pi[2] * t_p3pi[2];
			
			float px = t_fsPx[index], py = t_fsPy[index], pz = t_fsPz[index];
			
			float theta = TVector3(t_p3pi[0], t_p3pi[1], t_p3pi[2]).Theta()*TMath::RadToDeg();
			float e = (p2 + t_piKE * t_piKE) / (2000 * t_piKE);
			
			hist_TrueInformation->Fill(e, theta);

			// 3DST
			if ((t_pionExitMom[0] == 0) || (t_pionExitMom[1] == 0) || (t_pionExitMom[2] == 0))
			    hist_1->Fill(e, theta);	
			// 3DST + TPC
			else if (t_piGArLen > 20)
			    hist_2->Fill(e, theta);
			else // 3DST + TPC + ECAL
			    hist_3->Fill(e, theta);

			break;
		    }
		}
	    }
	}
    }

    cout << "Nb events: " << nevents << " CC events: " << CC_event_count << " Ratio: " << CC_event_count/float(nevents) << endl;
    cout << "Saving plots" << endl;

    auto hist_4 = (TH2F *)hist_1->Clone("hist_4");
    hist_4->Add(hist_3);
    auto hist_denom = (TH2F *)hist_4->Clone();
    hist_denom->Add(hist_2);
    hist_4->Divide(hist_denom);
    hist_4->SetTitle("Pion coverage within [3DST] + [3DST + TPC]");

    auto _file = new TFile("containmetPion.root", "recreate");
    hist_TrueInformation->Write();
    hist_1->Write();
    hist_2->Write();
    hist_3->Write();
    hist_4->Write();
    _file->Close();

    // Save file{{{
    /*
    auto can1 = new TCanvas("can1", "can", 900, 600);
    can1->SetGrid();
    hist_TrueInformation->Draw("colz");
    can1->SaveAs(Form("pion%s_trueInformation.root", outName.c_str()));

    auto can2 = new TCanvas("can2", "can", 900, 600);
    can2->SetGrid();
    hist_1->Draw("colz");
    can2->SaveAs(Form("pion%s_3DST.root", outName.c_str()));

    auto can3 = new TCanvas("can3", "can", 900, 600);
    can3->SetGrid();
    hist_2->Draw("colz");
    can3->SaveAs(Form("pion%s_3DST_TPC.root", outName.c_str()));

    auto can4 = new TCanvas("can4", "can", 900, 600);
    can4->SetGrid();
    hist_3->Draw("colz");
    can4->SaveAs(Form("pion%s_3DST_TPC_ECAL.root", outName.c_str()));

    auto can5 = new TCanvas("can5", "can", 900, 600);
    can5->SetGrid();
    auto hist_4 = (TH2F *)hist_1->Clone();
    hist_4->Add(hist_3);
    auto hist_denom = (TH2F *)hist_4->Clone();
    hist_denom->Add(hist_2);
    hist_4->Divide(hist_denom);
    hist_4->SetTitle("Pion coverage within [3DST] + [3DST + TPC]");
    hist_4->Draw("colz");
    can5->SaveAs(Form("pion%s_coverage.root", outName.c_str()));

    delete can1, can2, can3, can4, can5, 
	   hist_TrueInformation, hist_1, 
	   hist_2, hist_3, hist_4,
	   chain;
    */
    //}}}
}

int main(int argc, char* argv[]){
    pionAnalysis(argv[1]);

    return 0;
}
