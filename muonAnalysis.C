#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector3.h"

using namespace std;

void muonAnalysis(char* filename){
    // Read all files
    ifstream input(filename);
    
    auto outName = string(filename).substr(0, 4);

    auto chain = new TChain("tree");

    string file;
    while (input >> file) chain->Add(Form("%s", file.c_str()));

    int nevents = chain->GetEntries();

    Float_t t_p3lep[3], t_lepKE, t_muGArLen, t_muonExitMom[3];

    Int_t PDG, t_fsPdg[100], t_nFS;

    Float_t t_fsPx[100],t_fsPy[100],t_fsPz[100],
            t_vtx[3];
    
    chain->SetBranchAddress("vtx", &t_vtx);
    chain->SetBranchAddress("nFS", &t_nFS);
    chain->SetBranchAddress("fsPx", &t_fsPx);
    chain->SetBranchAddress("fsPy", &t_fsPy);
    chain->SetBranchAddress("fsPz", &t_fsPz);
    chain->SetBranchAddress("p3lep", &t_p3lep);
    chain->SetBranchAddress("lepKE", &t_lepKE);
    chain->SetBranchAddress("fsPdg", &t_fsPdg);
    cchain->SetBranchAddress("lepPdg", &PDG);
    hain->SetBranchAddress("muGArLen", &t_muGArLen);
    chain->SetBranchAddress("muonExitMom", &t_muonExitMom);

    auto hist_TrueInformation = new TH2F("hist_TrueInformation", 
					 "Muon True Information; Energy [GeV]; Angle",
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
		    if (abs(PDG) == 13){
			float p2 = t_p3lep[0] * t_p3lep[0] + 
				       t_p3lep[1] * t_p3lep[1] + 
				       t_p3lep[2] * t_p3lep[2];

			float theta = TVector3(t_p3lep[0], t_p3lep[1], t_p3lep[2]).Theta()*TMath::RadToDeg();
			float e = (p2 + t_lepKE * t_lepKE) / (2000 * t_lepKE);
			
			hist_TrueInformation->Fill(e, theta);

			// 3DST
			if ((t_muonExitMom[0] == 0) || (t_muonExitMom[1] == 0) || (t_muonExitMom[2] == 0))
			    hist_1->Fill(e, theta);	
			// 3DST + TPC
			else if (t_muGArLen > 20)
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
    hist_4->Add(hist_2);
    auto hist_denom = (TH2F *)hist_4->Clone();
    hist_denom->Add(hist_3);
    hist_4->Divide(hist_denom);
    hist_4->SetTitle("Muon coverage within [3DST] + [3DST + TPC]");

    auto _file = new TFile("containmetMuon.root", "recreate");
    hist_TrueInformation->Write();
    hist_1->Write();
    hist_2->Write();
    hist_3->Write();
    hist_4->Write();
    _file->Close();

    // Save files {{{
    /*
    auto can1 = new TCanvas("can1", "can", 900, 600);
    can1->SetGrid();
    hist_TrueInformation->Draw("colz");
    can1->SaveAs(Form("muon%s_trueInformation.root", outName.c_str(), hist_TrueInformation->GetTitle()));

    auto can2 = new TCanvas("can2", "can", 900, 600);
    can2->SetGrid();
    hist_1->Draw("colz");
    can2->SaveAs(Form("muon%s_3DST.root", outName.c_str(), hist_1->GetTitle()));

    auto can3 = new TCanvas("can3", "can", 900, 600);
    can3->SetGrid();
    hist_2->Draw("colz");
    can3->SaveAs(Form("muon%s_3DST_TPC.root", outName.c_str(), hist_2->GetTitle()));

    auto can4 = new TCanvas("can4", "can", 900, 600);
    can4->SetGrid();
    hist_3->Draw("colz");
    can4->SaveAs(Form("muon%s_3DST_TPC_ECAL.root", outName.c_str(), hist_3->GetTitle()));
   
    auto can5 = new TCanvas("can5", "can", 900, 600);
    can5->SetGrid();
    auto hist_4 = (TH2F *)hist_1->Clone();
    hist_4->Add(hist_2);
    auto hist_denom = (TH2F *)hist_4->Clone();
    hist_denom->Add(hist_3);
    hist_4->Divide(hist_denom);
    hist_4->SetTitle("Muon coverage within [3DST] + [3DST + TPC]");
    hist_4->Draw("colz");
    can5->SaveAs(Form("muon%s_coverage.root", outName.c_str(), hist_4->GetTitle()));

    delete can1, can2, can3, can4, can5, 
	   hist_TrueInformation, hist_1, 
	   hist_2, hist_3, hist_4,
	   chain;
	   */
    //}}}
}

int main(int argc, char* argv[]){
    muonAnalysis(argv[1]);

    return 0;
}
