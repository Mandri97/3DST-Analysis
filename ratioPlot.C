void ratioPlot(){
    auto _f1 = new TFile("containmetMuon.root", "read");
    auto _f2 = new TFile("containmetPion.root", "read");
    auto _f3 = new TFile("containmetHadron.root", "read");

    auto h1 = (TH2F*)_f1->Get("hist_4");
    auto h2 = (TH2F*)_f2->Get("hist_4");

    auto h3_3DST    = (TH1F*)_f3->Get("hist_3DST");
    auto h3_ECAL    = (TH1F*)_f3->Get("hist_ECAL");
    auto h3_allECAL = (TH1F*)_f3->Get("hist_allECAL");
    auto h3_leak    = (TH1F*)_f3->Get("hist_leak");
    
    auto h3_numerator   = (TH1F*)h3_3DST->Clone();
    h3_numerator->Add(h3_ECAL);
    auto h3_denominator = (TH1F*)h3_numerator->Clone();
    h3_denominator->Add(h3_leak);

    auto ratio = (TH1F*)h3_numerator->Clone();
    ratio->Divide(h3_denominator);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);

    auto c = new TCanvas("c", "c", 600, 900);

    /*********** Lever arm **********/
    auto pad1 = new TPad("pad1", "pad1", 0., 0.65, 1., 1.);
    pad1->SetBottomMargin(0);
    pad1->SetGridx();
    pad1->Draw();
    pad1->cd();

    h1->GetYaxis()->SetLabelSize(0.0);
    h1->GetYaxis()->SetTitleFont(43);
    h1->GetYaxis()->SetTitleSize(15);
    h1->GetYaxis()->SetTitleOffset(2);
    h1->GetYaxis()->SetTitle("Angle");

    auto pal = new TPaletteAxis(20.5, 10, 21.5, 170, h1);
    h1->GetListOfFunctions()->Add(pal);
    gStyle->SetNumberContours(255);

    h1->Draw("colz");

    auto axis1 = new TGaxis(0, 20, 0, 180, 20, 180, 8, "");
    axis1->SetLabelFont(43);
    axis1->SetLabelSize(15);
    axis1->Draw();
  
    /*********** Purity *************/
    c->cd();
    auto pad2 = new TPad("pad2", "pad2", 0., 0.3, 1., 0.65);
    pad2->SetBottomMargin(0);
    pad2->SetTopMargin(0);
    pad2->SetGridx();
    pad2->Draw();
    pad2->cd();

    pal = new TPaletteAxis(20.5, 10, 21.5, 170, h2);
    h2->GetListOfFunctions()->Add(pal);
    gStyle->SetNumberContours(255);

    h2->GetYaxis()->SetLabelSize(0.0);
    h2->GetYaxis()->SetLabelSize(0.);
    h2->GetYaxis()->SetTitleFont(43);
    h2->GetYaxis()->SetTitleSize(15);
    h2->GetYaxis()->SetTitleOffset(2);
    h2->GetYaxis()->SetTitle("Angle");

    h2->Draw("colz1");

    auto axis2 = new TGaxis(0, 20, 0, 160, 20, 160, 8, "");
    axis2->SetLabelFont(43);
    axis2->SetLabelSize(15);
    axis2->Draw();

    /************* Resolution ***************/
    c->cd();
    auto pad3 = new TPad("pad3", "pad3", 0., 0., 1., 0.3);
    pad3->SetTopMargin(0);
    pad3->SetBottomMargin(0.2);
    pad3->SetGridx();
    pad3->Draw();
    pad3->cd();

    ratio->GetYaxis()->SetLabelSize(0.);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(15);
    ratio->GetYaxis()->SetTitleOffset(2);
    ratio->GetYaxis()->SetTitle("Hadronic containemnt efficiency");
    //h3->GetYaxis()->SetRange(0, 15);

    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(15);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleSize(15);
    ratio->GetXaxis()->SetTitleOffset(3.67);

    ratio->Draw("hist");

    auto axis3 = new TGaxis(0, 0, 0, 1, 0, 1, 5, "");
    axis3->SetLabelFont(43);
    axis3->SetLabelSize(15);
    axis3->Draw();
}
