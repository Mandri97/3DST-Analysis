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

void hadronAnalysis(char* filename){
    ifstream input(filename);
    auto outName = string(filename).substr(0, 4);
    auto chain = new TChain("tree");

    string file;
    while (input >> file) chain->Add(Form("%s", file.c_str()));

    int nevents = chain->GetEntries();

    Float_t hadTot_3DST, hadTot_TPC, hadTot_ECAL, hadTot_allECAL, hadTot_leak;

    auto hist_3DST = new TH1F("hist_3DST", "Hadronic Energy Deposit within 3DST; Energy [GeV]; Number",
			60, 0, 20);
    auto hist_TPC = new TH1F("hist_TPC", "Hadronic Energy Deposit within TPC; Energy [Gev]; Number",
			60, 0, 20);
    auto hist_ECAL = new TH1F("hist_ECAL", "Hadronic Energy Deposit within ECAL (without radiator); Energy [GeV]; Number", 
			60, 0, 20);
    auto hist_allECAL = new TH1F("hist_allECAL", "Hadronic Energy Deposit within ECAL (with radiator); Energy [GeV]; Number", 
			60, 0, 20);
    auto hist_leak = new TH1F("hist_leak", "Leaked Hadronic Energy Deposit; Energy [GeV]; Number",
			60, 0, 20);

    Float_t t_3DST, t_TPC, t_ECAL, t_allECAL, t_leak;
    Float_t t_vtx[3];
    Int_t t_nFS, t_fsPdg[100];

    chain->SetBranchAddress("hadTot_3DST", &t_3DST);
    chain->SetBranchAddress("hadTot_TPC", &t_TPC);
    chain->SetBranchAddress("hadTot_ECAL", &t_ECAL);
    chain->SetBranchAddress("hadTot_allECAL", &t_allECAL);
    chain->SetBranchAddress("hadTot_leak", &t_leak);
    chain->SetBranchAddress("nFS", &t_nFS);
    chain->SetBranchAddress("fsPdg", &t_fsPdg);
    chain->SetBranchAddress("vtx", &t_vtx);


    int CC_event_count = 0;
    for (int i = 0; i < nevents; ++i){
	chain->GetEntry(i);

	// TODO: Change this according to the geometry
	if (abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50){

	    for (int index = 0; index < t_nFS; index++){
		// looking for CC events
		if (abs(t_fsPdg[index]) == 11 || abs(t_fsPdg[index]) == 13 ){
		    CC_event_count++;
		    hist_3DST->Fill(t_3DST/1000);
		    hist_TPC->Fill(t_TPC/1000);
		    hist_ECAL->Fill(t_ECAL/1000);
		    hist_allECAL->Fill(t_allECAL/1000);
		    hist_leak->Fill(t_leak/1000);

		    break;
		}
	    }
	}
    }

    cout << "Nb events: " << nevents << " CC events: " << CC_event_count << " Ratio: " << CC_event_count/float(nevents) << endl;
    
    auto hist_numerator = (TH1F *)hist_ECAL->Clone();
    auto hist_denom = (TH1F *)hist_numerator->Clone();
    hist_denom->Add(hist_allECAL);
    hist_denom->Add(hist_leak);
    hist_denom->Add(hist_3DST);
    auto hist_ratio = (TH1F *)hist_numerator->Clone("hist_ratio");
    hist_ratio->Divide(hist_denom);
    hist_ratio->SetTitle("Hadronic energy deposit ratio ECAL; Energy [GeV]; Efficiency");

    auto _file = new TFile("containmetHadron.root", "recreate");
    hist_3DST->Write();
    hist_TPC->Write();
    hist_leak->Write();
    hist_ECAL->Write();
    hist_allECAL->Write();
    hist_ratio->Write();
    _file->Close();

    // Save file {{{
    /*
    auto can5 = new TCanvas("can5", "can", 900, 600);
    can5->SetGrid();
    auto can1 = new TCanvas("can1", "can", 900, 600);
    can1->SetGrid();
    hist_3DST->Draw();
    can1->SaveAs(Form("hadron%s_3DST.root", outName.c_str()));
    
    auto can2 = new TCanvas("can2", "can", 900, 600);
    can2->SetGrid();
    hist_TPC->Draw();
    can2->SaveAs(Form("hadron%s_TPC.root", outName.c_str()));

    auto can3 = new TCanvas("can3", "can", 900, 600);
    can3->SetGrid();
    hist_ECAL->Draw();
    can3->SaveAs(Form("hadron%s_ECAL.root", outName.c_str()));
    
    auto can4 = new TCanvas("can4", "can", 900, 600);
    can4->SetGrid();
    hist_allECAL->Draw();
    can4->SaveAs(Form("hadron%s_allECAL.root", outName.c_str()));

    auto can5 = new TCanvas("can5", "can", 900, 600);
    can5->SetGrid();
    hist_leak->Draw();
    can5->SaveAs(Form("hadron%s_leak.root", outName.c_str()));
    
    // Ratio
    auto can6 = new TCanvas("can6", "can", 900, 600);
    can6->SetGrid();
    auto hist_numerator = (TH1F *)hist_3DST->Clone();
    hist_numerator->Add(hist_ECAL);
    hist_numerator->Add(hist_TPC);
    auto hist_denom = (TH1F *)hist_numerator->Clone();
    hist_denom->Add(hist_allECAL);
    hist_denom->Add(hist_leak);
    auto hist_ratio = (TH1F *)hist_numerator->Clone();
    hist_ratio->Divide(hist_denom);
    hist_ratio->SetTitle("Hadronic energy deposit ratio; Energy [GeV]; Efficiency");
    hist_ratio->Draw();
    can6->SaveAs(Form("hadron%s_Efficiency.root", outName.c_str(), hist_ratio->GetTitle()));

    auto can7 = new TCanvas("can7", "can", 900, 600);
    can7->SetGrid();
    auto hist_numerator = (TH1F *)hist_3DST->Clone();
    auto hist_denom = (TH1F *)hist_numerator->Clone();
    hist_denom->Add(hist_allECAL);
    hist_denom->Add(hist_leak);
    hist_denom->Add(hist_ECAL);
    auto hist_ratio = (TH1F *)hist_numerator->Clone();
    hist_ratio->Divide(hist_denom);
    hist_ratio->SetTitle("Hadronic energy deposit ratio 3DST; Energy [GeV]; Efficiency");
    hist_ratio->Draw();
    can7->SaveAs(Form("hadron%s_Efficiency_3DST.root", outName.c_str(), hist_ratio->GetTitle()));

    auto can8 = new TCanvas("can8", "can", 900, 600);
    can8->SetGrid();
    auto hist_numerator = (TH1F *)hist_ECAL->Clone();
    auto hist_denom = (TH1F *)hist_numerator->Clone();
    hist_denom->Add(hist_allECAL);
    hist_denom->Add(hist_leak);
    hist_denom->Add(hist_3DST);
    auto hist_ratio = (TH1F *)hist_numerator->Clone();
    hist_ratio->Divide(hist_denom);
    hist_ratio->SetTitle("Hadronic energy deposit ratio ECAL; Energy [GeV]; Efficiency");
    hist_ratio->Draw();
    can8->SaveAs(Form("hadron%s_Efficiency_ECAL.root", outName.c_str(), hist_ratio->GetTitle()));

    delete can1, can2, can3, can4, can5, can6, can7, can8,
	   hist_allECAL, hist_leak, hist_numerator, 
	   hist_3DST, hist_TPC, hist_ECAL,
	   hist_denom, hist_ratio, chain;
    */
  // }}}
}

int main(int argc, char* argv[]){
    hadronAnalysis(argv[1]);

    return 0;
}
