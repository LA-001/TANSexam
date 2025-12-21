#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TNtupleD.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TStopwatch.h"
#include "TEfficiency.h"

#include "Defstruct.h"

using namespace std;

void Efficienza() {

	TFile *fin1 = TFile::Open("fileRoot/residui.root");
    if (!fin1 || fin1->IsZombie()) {
        cout << "Errore: impossibile aprire il file ROOT." << endl;
        return;
    }

    TFile *fin2 = TFile::Open("fileRoot/simulazione.root");
    if (!fin2 || fin2->IsZombie()) {
        cout << "Errore: impossibile aprire il file ROOT." << endl;
        return;
    }

    TTree *T_vtx = (TTree*)fin1->Get("T_vtx");
    TTree *T_vrt = (TTree*)fin2->Get("T_vrt");
    if (!T_vtx || !T_vrt) {
        cout << "Errore: TTree non trovati." << endl;
		fin1->Close();
		fin2->Close();
		delete fin1;
		delete fin2;
        return;
    }

    Ric vertex;
    Vrt MC;
    T_vtx->SetBranchAddress("vtx", &vertex);
    T_vrt->SetBranchAddress("vrt", &MC);

	double edgesMolti[] = {2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
	int nBinsMolti = sizeof(edgesMolti)/sizeof(edgesMolti[0]) - 1;

    const int nBinZ = 16;
    const double minZ = -160;
    const double maxZ = 160;
    int bin_M, bin_V;

    TH1F h_tot_M("h_tot_M","Totale;Molteplicità;Conteggi", nBinsMolti, edgesMolti);
    TH1F h_pass_M("h_pass_M","Passati;Molteplicità;Conteggi", nBinsMolti, edgesMolti);
    TH1F h_tot_V("h_tot_V","Totale;Vertice;Conteggi", nBinZ, minZ, maxZ);
    TH1F h_pass_V("h_pass_V","Passati;Vertice;Conteggi", nBinZ, minZ, maxZ);
 
    for (long ev = 0; ev < T_vrt->GetEntries(); ++ev) {
        T_vrt->GetEntry(ev);
        h_tot_M.Fill(MC.moltiplicita);
        h_tot_V.Fill(MC.z0);
    }

    for (long ev = 0; ev < T_vtx->GetEntries(); ++ev) {
        T_vtx->GetEntry(ev);
        h_pass_M.Fill(vertex.molti);
        h_pass_V.Fill(vertex.z0);
    }

    TEfficiency* eff = new TEfficiency(h_pass_M, h_tot_M);
    eff->SetName("eff");
    eff->SetTitle("Efficienza vs molteplicita'; Molteplicita'; #epsilon");
    TEfficiency* eff_v = new TEfficiency(h_pass_V, h_tot_V);
    eff_v->SetName("eff_v");
    eff_v->SetTitle("Efficienza vs Z_{true}; Z_{true} [mm]; #epsilon");

	TCanvas *c6 = new TCanvas("c6", "Efficienza vs molteplicita'", 1000, 700);
    c6->SetLeftMargin(0.1);
    eff->SetMarkerStyle(20);
    eff->SetMarkerSize(1);
    TEfficiency *eff_clone = (TEfficiency*)eff->Clone("eff_clone");
    eff_clone->SetDirectory(0);
    eff_clone->Draw("APL");
    c6->SetGrid();
    c6->SaveAs("Plot/c6.png");

    TCanvas *c7 = new TCanvas("c7", "Efficienza vs Z_{true}", 1000, 700);
    c7->SetLeftMargin(0.1);
    eff_v->SetMarkerStyle(20);
    eff_v->SetMarkerSize(1);
    TEfficiency *eff_v_clone = (TEfficiency*)eff_v->Clone("eff_v_clone");
    eff_v_clone->SetDirectory(0);
    eff_v_clone->Draw("AP");
    c7->SetGrid();
    c7->SaveAs("Plot/c7.png");

    TFile fout("fileRoot/efficienza.root","RECREATE");
    eff->Write();
    eff_v->Write();
    fout.Close();

    fin1->Close();
    fin2->Close();
}
