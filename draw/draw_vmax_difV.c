#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;
void draw_vmax_difV()
{
    TFile *f_sim1 = new TFile("./dataout20V.root");
    TFile *f_sim2 = new TFile("./dataout40V.root");
    TFile *f_sim3 = new TFile("./dataout60V.root");
    TFile *f_sim4 = new TFile("./dataout80V.root");
    TFile *f_sim5 = new TFile("./dataout100V.root");
    TFile *f_sim6 = new TFile("./dataout120V.root");
    TFile *f_sim7 = new TFile("./dataout140V.root");
    TFile *f_sim8 = new TFile("./dataout160V.root");
    TFile *f_sim9 = new TFile("./dataout180V.root");
    TFile *f_sim10 = new TFile("./dataout200V.root");

    TCanvas *c = new TCanvas("c", "", 100, 100, 800, 600);
    c->Divide(1,1);
    c->cd(1);

	//////RASER///////
    TTree *chain1 = (TTree *)f_sim1->Get("tree");
	TH2F *h1 = new TH2F("RASER1","",4000,-10.5,60.5,1000,0,30);
	chain1->Draw("(-v*1.3e2):(z)>>RASER1");

    TTree *chain2 = (TTree *)f_sim2->Get("tree");
	TH2F *h2 = new TH2F("RASER2","",4000,-10.5,60.5,1000,0,30);
	chain2->Draw("(-v*1.3e2):(z)>>RASER2");

    TTree *chain3 = (TTree *)f_sim3->Get("tree");
	TH2F *h3 = new TH2F("RASER3","",4000,-10.5,60.5,1000,0,30);
	chain3->Draw("(-v*1.3e2):(z)>>RASER3");

    TTree *chain4 = (TTree *)f_sim4->Get("tree");
	TH2F *h4 = new TH2F("RASER4","",4000,-10.5,60.5,1000,0,30);
	chain4->Draw("(-v*1.3e2):(z)>>RASER4");

    TTree *chain5 = (TTree *)f_sim5->Get("tree");
	TH2F *h5 = new TH2F("RASER5","",4000,-10.5,60.5,1000,0,30);
	chain5->Draw("(-v*1.3e2):(z)>>RASER5");

    TTree *chain6 = (TTree *)f_sim6->Get("tree");
	TH2F *h6 = new TH2F("RASER6","",4000,-10.5,60.5,1000,0,30);
	chain6->Draw("(-v*1.3e2):(z)>>RASER6");

    TTree *chain7 = (TTree *)f_sim7->Get("tree");
	TH2F *h7 = new TH2F("RASER7","",4000,-10.5,60.5,1000,0,30);
	chain7->Draw("(-v*1.3e2):(z)>>RASER7");

    TTree *chain8 = (TTree *)f_sim8->Get("tree");
	TH2F *h8 = new TH2F("RASER8","",4000,-10.5,60.5,1000,0,30);
	chain8->Draw("(-v*1.3e2):(z)>>RASER8");

    TTree *chain9 = (TTree *)f_sim9->Get("tree");
	TH2F *h9 = new TH2F("RASER9","",4000,-10.5,60.5,1000,0,30);
	chain9->Draw("(-v*1.3e2):(z)>>RASER9");

    TTree *chain10 = (TTree *)f_sim10->Get("tree");
	TH2F *h10 = new TH2F("RASER10","",4000,-10.5,60.5,1000,0,30);
	chain10->Draw("(-v*1.3e2):(z)>>RASER10");

    h1->GetXaxis()->SetTitle("z or z'[um]");
	h1->GetYaxis()->SetTitle("Vmax[mv]");

	h1->SetMarkerColor(1);    // 1
	h1->SetMarkerStyle(33);
	h1->SetMarkerSize(1);
	h1->SetStats(0);
	h1->Draw();
	h2->SetMarkerColor(2);
	h2->SetMarkerStyle(33);
	h2->SetMarkerSize(1);
	h2->Draw("same");
	h3->SetMarkerColor(3);
	h3->SetMarkerStyle(33);
	h3->SetMarkerSize(1);
	h3->Draw("same");
    h4->SetMarkerColor(4);
	h4->SetMarkerStyle(33);
	h4->SetMarkerSize(1);
	h4->Draw("same");
    h5->SetMarkerColor(5);
	h5->SetMarkerStyle(33);
	h5->SetMarkerSize(1);
	h5->Draw("same");
    h6->SetMarkerColor(6);
	h6->SetMarkerStyle(33);
	h6->SetMarkerSize(1);
	h6->Draw("same");
    h7->SetMarkerColor(7);
	h7->SetMarkerStyle(33);
	h7->SetMarkerSize(1);
	h7->Draw("same");
    h8->SetMarkerColor(8);
	h8->SetMarkerStyle(33);
	h8->SetMarkerSize(1);
	h8->Draw("same");
    h9->SetMarkerColor(9);
	h9->SetMarkerStyle(33);
	h9->SetMarkerSize(1);
	h9->Draw("same");
    h10->SetMarkerColor(11);
	h10->SetMarkerStyle(33);
	h10->SetMarkerSize(1);
	h10->Draw("same");

    TLegend *leg1=new TLegend(0.65,0.6,0.98,0.88,NULL,"brNDC");

	/////RASER/////
	leg1->AddEntry(h1,"20V_RASER");
	leg1->AddEntry(h2,"40V_RASER");
	leg1->AddEntry(h3,"60V_RASER");
    leg1->AddEntry(h4,"80V_RASER");
    leg1->AddEntry(h5,"100V_RASER");
    leg1->AddEntry(h6,"120V_RASER");
    leg1->AddEntry(h7,"140V_RASER");
    leg1->AddEntry(h8,"160V_RASER");
    leg1->AddEntry(h9,"180V_RASER");
    leg1->AddEntry(h10,"200V_RASER");

	leg1->SetBorderSize(1);
	leg1->SetTextFont(62);
	leg1->SetTextSize(0.028);
	leg1->SetLineColor(0);
	leg1->SetLineStyle(1);
	leg1->SetLineWidth(3);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);

	leg1->Draw();
	c->Print("Vmax_pin_raser.pdf");
}