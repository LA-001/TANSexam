#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TH1D.h"
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

    const int nMolt = 11, nVer = 9;
    int inizio = -160;

    TH1D *hist[nMolt];
    TH1D *hist_v[nVer];
    char titolo[50];
    char nome[30];

    TCanvas *cv[nVer];

    double molteplicita[nMolt] = { 5., 7., 9., 12., 15., 18., 22., 26., 32., 40., 48. };
    double deltaM[nMolt] = { 0.5, 0.5, 0.5, 0.5 ,1. ,1.  ,2.  ,2.  ,3.  ,3., 4. };

    double vertice[nVer];
    double deltaV[nVer];

    for (int i = 0; i < nVer; i++) {
        vertice[i] = inizio;
        inizio += 40;
        deltaV[i] = 10;
    }

    double y1[nMolt], y2[nVer], ey1[nMolt], ey2[nVer];

    for (int i = 0; i < nMolt; i++) {
        snprintf(nome, 15, "histo%i", i+1);
        snprintf(titolo, 50, "Istogramma residui #%i ", i+1);
        hist[i] = new TH1D(nome, titolo, 70, -700, 700);
    }

    for (long ev = 0; ev < T_vtx->GetEntries(); ev++) {
        T_vtx->GetEntry(ev);
        for (int i = 0; i < nMolt; i++) {
            if (-159. <= vertex.z0 && vertex.z0 <= 159.) {
                if ((molteplicita[i] - deltaM[i]) <= vertex.molti &&
                    vertex.molti <= (molteplicita[i] + deltaM[i])) {
                    hist[i]->Fill((vertex.zRic - vertex.z0) * 1000);
                }
            }
        }
    }

    for (int i = 0; i < nMolt; i++) {
        y1[i] = hist[i]->GetStdDev();
        ey1[i] = hist[i]->GetStdDevError();
    }

    for (int i = 0; i < nVer; i++) {
        snprintf(nome, 15, "histo_v%i", i+1);
        snprintf(titolo, 50, "Istogramma residui vs Z_{true} #%i ", i+1);
        hist_v[i] = new TH1D(nome, titolo, 70, -700, 700);
    }

    for (long ev = 0; ev < T_vtx->GetEntries(); ev++) {
        T_vtx->GetEntry(ev);
        for (int i = 0; i < nVer; i++) {
            if (-159. <= vertex.z0 && vertex.z0 <= 159.) {
                if ((vertice[i] - deltaV[i]) <= vertex.z0 &&
                    vertex.z0 <= (vertice[i] + deltaV[i])) {
                    hist_v[i]->Fill((vertex.zRic - vertex.z0) * 1000);
                }
            }
        }
    }

    for (int i = 0; i < nVer; i++) {
        y2[i] = hist_v[i]->GetStdDev();
        ey2[i] = hist_v[i]->GetStdDevError();
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
}
