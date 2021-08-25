#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;
void temp_plt_Vmax_pin()
{
	TFile *f_sim1 = new TFile("./dataout120V.root");
	TFile *f_data = new TFile("./data/edge_voltage_2019_10_24_15_12_57_HPK-EPI-W2-200-DS-SE5PINNM-01.root");
	TCanvas *c = new TCanvas("c", "", 100, 100, 800, 600);
	c->Divide(1,1);
	c->cd(1);

	TTree *chain1 = (TTree *)f_sim1->Get("tree");
	TH2F *h1 = new TH2F("RASER","",4000,-10.5,60.5,1000,0,30);
	chain1->Draw("(-v*1.3e2):(z)>>RASER");
		
	TTree *chain2 = (TTree *)f_data->Get("edge");
	TH2F *h2 = new TH2F("data","",4000,-10.5,60.5,1000,0,30);
	chain2->Draw("(-Vmin*1e3):((z*1e3)-11985)>>data","Vbias==-120");

	h1->GetXaxis()->SetTitle("z or z'[um]");
	h1->GetYaxis()->SetTitle("Vmax[mv]");

	h1->SetMarkerColor(1);
	h1->SetMarkerStyle(42);
	h1->SetMarkerSize(1);
	h1->SetStats(0);
	h1->Draw();
	h2->SetMarkerColor(2);
	h2->SetMarkerStyle(42);
	h2->SetMarkerSize(1);
	h2->Draw("same");

	TLegend *leg1=new TLegend(0.14,0.75,0.4,0.85,NULL,"brNDC");

	leg1->AddEntry(h1,"120V_RASER");
	leg1->AddEntry(h2,"120V_TCT_HPK(z'=z-11985)");

	leg1->SetBorderSize(1);
	leg1->SetTextFont(62);
	leg1->SetTextSize(0.028);
	leg1->SetLineColor(0);
	leg1->SetLineStyle(1);
	leg1->SetLineWidth(3);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);

	leg1->Draw();
	c->Print("figure/Vmax_hpk_pin_com_TF_bv120.pdf");


}
