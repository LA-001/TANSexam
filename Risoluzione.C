#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "Defstruct.h"

void Risoluzione() {
    TFile *fin = TFile::Open("fileRoot/residui.root");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Errore nell'aprire fileRoot/output.root\n";
        return;
    }

    TTree *T_vtx = (TTree*)fin->Get("T_vtx");
    if (!T_vtx) {
        std::cerr << "TTree 'T_vtx' non trovato\n";
        return;
    }

    Ric vertex;
    
    T_vtx->SetBranchAddress("vtx", &vertex);

    TH2D *hist2D_M = new TH2D("hist2D_M","Istogramma 2D - Residui vs Molteplicita'",60,0.5,60.5,600,-650.,650.);
    hist2D_M->GetXaxis()->SetTitle("Molteplicita'");
    hist2D_M->GetYaxis()->SetTitle("Residui [#mum]");
    TH2D *hist2D_V = new TH2D("hist2D_V","Istogramma 2D - Residui vs Z_{vert}",340,-170.,170.,600,-650.,650.);
    hist2D_V->GetXaxis()->SetTitle("Z_{vert} [mm]");
    hist2D_V->GetYaxis()->SetTitle("Residui [#mum]");

    const int nMolt = 11, nVer = 9;
 
    double molteplicita[nMolt] = { 5., 7., 9., 12., 15., 18., 22., 26., 32., 40., 48. };
    double deltaM[nMolt] = { 0.5, 0.5, 0.5, 0.5 ,1. ,1.  ,2.  ,2.  ,3.  ,3., 4. };

    double vertice[nVer] = {-160., -120., -80., -40., 0., 40., 80., 120., 160.};
    double deltaV[nVer] = {10., 10., 10., 10., 10., 10., 10., 10., 10.};

    double y1[nMolt], y2[nVer], ey1[nMolt], ey2[nVer];

    for (long ev = 0; ev < T_vtx->GetEntries(); ev++) {
        T_vtx->GetEntry(ev);
        hist2D_M->Fill(vertex.molti, (vertex.zRic - vertex.z0) * 1000);
        hist2D_V->Fill(vertex.z0, (vertex.zRic - vertex.z0) * 1000);
    }

    for (int i = 0; i < nMolt; i++) {

        int binM_min = hist2D_M->GetXaxis()->FindBin(molteplicita[i] - deltaM[i]);
        int binM_max = hist2D_M->GetXaxis()->FindBin(molteplicita[i] + deltaM[i]);

        TH1D *proj_M = hist2D_M->ProjectionY(Form("proj_M_%d", i), binM_min, binM_max);
        proj_M->SetDirectory(0);

        y1[i]  = proj_M->GetStdDev();
        ey1[i] = proj_M->GetStdDevError();

        delete proj_M;
    }


    for (int i = 0; i < nVer; i++) {

        int binV_min = hist2D_V->GetXaxis()->FindBin(vertice[i] - deltaV[i]);
        int binV_max = hist2D_V->GetXaxis()->FindBin(vertice[i] + deltaV[i]);

        TH1D *proj_V = hist2D_V->ProjectionY(Form("proj_V_%d", i),binV_min,binV_max);
        proj_V->SetDirectory(0);

        y2[i]  = proj_V->GetStdDev();
        ey2[i] = proj_V->GetStdDevError();

        delete proj_V;
    }


    TGraphErrors* ris1 = new TGraphErrors(nMolt, molteplicita, y1, deltaM, ey1);
    TGraphErrors* ris2 = new TGraphErrors(nVer, vertice, y2, deltaV, ey2);

    ris1->SetLineColor(kGreen+3);
    ris1->SetLineWidth(1);
    ris1->SetMarkerColor(kBlack);
    ris1->SetMarkerSize(1);
    ris1->SetMarkerStyle(20);
    ris1->SetTitle("Risoluzione vs Molteplicita'");
    ris1->GetXaxis()->SetTitle("Molteplicita'");
    ris1->GetYaxis()->SetTitle("Risoluzione [#mum]");

    ris2->SetLineColor(kGreen+3);
    ris2->SetLineWidth(1);
    ris2->SetMarkerColor(kBlack);
    ris2->SetMarkerSize(1);
    ris2->SetMarkerStyle(20);
    ris2->SetTitle("Risoluzione vs Z_{vert}");
    ris2->GetXaxis()->SetTitle("Z_{vert} [mm]");
    ris2->GetYaxis()->SetTitle("Risoluzione [#mum]");

    TCanvas *c4 = new TCanvas("c4", "c4", 900, 600);
    ris1->Draw("ACP");
    c4->SetGrid();
    c4->SaveAs("Plot/c4.png");

    TCanvas *c5 = new TCanvas("c5", "c5", 900, 600);
    ris2->Draw("ALP");
    c5->SetGrid();
    c5->SaveAs("Plot/c5.png");

    TCanvas *c4_2 = new TCanvas("c4_2", "c4_2", 900, 600);
    hist2D_M->SetStats(0);
    hist2D_M->Draw("COLZ");
    c4_2->SaveAs("Plot/c4_2.png");

    TCanvas *c5_2 = new TCanvas("c5_2", "c5_2", 900, 600);
    hist2D_V->SetStats(0);
    hist2D_V->Draw("COLZ");
    c5_2->SaveAs("Plot/c5_2.png");

}
