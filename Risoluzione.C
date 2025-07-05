#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtupleD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStopwatch.h>

void Risoluzione(){

    TFile *fin = TFile::Open("fileRoot/output.root");
	TNtupleD *vtx = (TNtupleD*)fin->Get("vtx");

	TFile fout("fileRoot/risoluzione.root","RECREATE");

	double z0, molti, zRic;

	vtx->SetBranchAddress("z0", &z0);
    vtx->SetBranchAddress("molti", &molti);
	vtx->SetBranchAddress("zRic", &zRic);

	const int nMolt = 11, nVer = 12;

	TH1F *hist[nMolt];
	TH1F *hist_v[nVer];
	char titolo[50];
	char nome[30];

	double molteplicita[nMolt] = { 5., 7., 9., 12., 15., 18., 22., 26., 32., 40., 48.};
	double deltaM[nMolt] =       { 0.5, 0.5, 0.5, 0.5 ,1. ,1.  ,2.  ,2.  ,3.  ,3., 4.};

	double vertice[nVer] = {-147., -115., -90., -60., -35., -10., 10., 35., 60., 90., 115., 147.}; 
	double deltaV[nVer]  = {10.,10.,10.,10.,5.,5.,5.,5.,10.,10.,10.,10.};

	double y1[nMolt], y2[nVer], ey1[nMolt], ey2[nVer];

	TF1 *f1 = new TF1("f1", "[0]*TMath::Gaus(x, [1], [2], true)", -500, 500);
	f1->SetParameters(1000,0.,100);

	//********************************************************************************************************

	for(int i = 0; i < nMolt; i++){
		snprintf(nome,15,"histo%i",i+1);
		snprintf(titolo,50,"Istogramma residui #%i ",i+1);
		hist[i] = new TH1F(nome,titolo, 600, -600, 600);
	}

	for(int ev = 0; ev < vtx->GetEntries(); ev++){
		vtx->GetEntry(ev);
		for(int i = 0; i < nMolt; i++){
			if(-159. <= z0 && z0 <= 159.){
				if( (molteplicita[i] - deltaM[i]) <= molti && molti <= (molteplicita[i] + deltaM[i])){
					hist[i]->Fill((zRic - z0)*1000);
				}
			}
		}
	}

	for(int i = 0; i < nMolt; i++){
		//hist[i]->Fit("f1","RQ+N");
		y1[i] = hist[i]->GetStdDev(); 
		ey1[i] = hist[i]->GetStdDevError();
	}

	//********************************************************************************************************

	for(int i = 0; i < nVer; i++){
		snprintf(nome,15,"histo_v%i",i+1);
		snprintf(titolo,50,"Istogramma residui vs Z_{true} #%i ",i+1);
		hist_v[i] = new TH1F(nome,titolo, 150, -700, 700);
	}

	for(int ev = 0; ev < vtx->GetEntries(); ev++){
		vtx->GetEntry(ev);
		for(int i = 0; i < nVer; i++){
			if(-159. <= z0 && z0 <= 159.){
				if( (vertice[i] - deltaV[i]) <= z0 && z0 <= (vertice[i] + deltaV[i])){
					hist_v[i]->Fill((zRic - z0)*1000);
				}
			}
		}
	}

	for(int i = 0; i < nVer; i++){
		//hist_v[i]->Fit("f1","RQ+");
		y2[i] = hist_v[i]->GetStdDev(); 
		ey2[i] = hist_v[i]->GetStdDevError();
	}

	//********************************************************************************************************

	TGraphErrors* ris1 = new TGraphErrors(nMolt, molteplicita, y1, deltaM, ey1);
	TGraphErrors* ris2 = new TGraphErrors(nVer, vertice, y2, deltaV, ey2);

    ris1->SetLineColor(kGreen+3);
    ris1->SetLineWidth(1);
    ris1->SetMarkerColor(kBlack);
    ris1->SetMarkerSize(1);
    ris1->SetMarkerStyle(20);
    ris1->SetTitle("Risoluzione vs Molteplicita'");
    ris1->GetXaxis()->SetTitle("Molteplicita'");
    ris1->GetYaxis()->SetTitle("Risoluzione [micron]");

	ris2->SetLineColor(kGreen+3);
    ris2->SetLineWidth(1);
    ris2->SetMarkerColor(kBlack);
    ris2->SetMarkerSize(1);
    ris2->SetMarkerStyle(20);
    ris2->SetTitle("Risoluzione vs Z_{vert}");
    ris2->GetXaxis()->SetTitle("Z_{vert} [mm]");
    ris2->GetYaxis()->SetTitle("Risoluzione [micron]");

	TGraphErrors* copy1 = (TGraphErrors*)ris1->Clone("copy1");
	TGraphErrors* copy2 = (TGraphErrors*)ris2->Clone("copy2");

    TCanvas *c4 = new TCanvas("c4", "c4", 900, 600);
	copy1->Draw("ACP");
    c4->SetGrid();
    c4->Update();
	c4->SaveAs("Plot/c4.png");

	TCanvas *c5 = new TCanvas("c5", "c5", 900, 600);
	copy2->Draw("ACP");
    c5->SetGrid();
    c5->Update();
	c5->SaveAs("Plot/c5.png");

	ris1->Write();
	ris2->Write();
	fout.Close();

	fin->Close();
	delete fin;
}
