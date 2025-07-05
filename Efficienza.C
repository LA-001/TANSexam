#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtupleD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TStopwatch.h>
#include <TEfficiency.h>

void Efficienza(){

	TFile *fin1 = TFile::Open("fileRoot/output.root");
    TFile *fin2 = TFile::Open("fileRoot/simulazione.root");
	TNtupleD *vtx = (TNtupleD*)fin1->Get("vtx");
    TNtupleD *vrt = (TNtupleD*)fin2->Get("vrt");

	//TFile fout("fileRoot/efficienza.root","RECREATE");

	double molti, moltiRic, z0, Zric;

    vrt->SetBranchAddress("moltiplicita", &molti);
    vrt->SetBranchAddress("z0", &z0);
    vtx->SetBranchAddress("molti", &moltiRic);
    vtx->SetBranchAddress("z0", &Zric);

	TEfficiency* eff = new TEfficiency("eff", "Efficienza vs molteplicita'; molteplicita' ; #epsilon", 11, 0, 55);
	TEfficiency* eff_v = new TEfficiency("eff_v", "Efficienza vs Z_{true}; Z_{true} [mm] ; #epsilon", 16, -160, 160);

    double conta = 0, contaRic = 0;

//***********************************************************************************************************************

    for(int i = 1; i <= 11; i++){
		for(int ev = 0; ev < vtx->GetEntries(); ev++){
			vtx->GetEntry(ev);
			if((i*5 - 2.5) <= moltiRic && moltiRic <= (i*5 + 2.5)){
				contaRic++;	  
            }                        
        }
        for(int ev = 0; ev < vrt->GetEntries(); ev++){
			vrt->GetEntry(ev);
			if((i*5 - 2.5) <= molti && molti <= (i*5 + 2.5)){
				conta++;	  
            }                       
        }

		eff->SetTotalEvents(i , conta); 
        eff->SetPassedEvents(i , contaRic);

		contaRic = 0;
        conta = 0;	
	}

//***********************************************************************************************************************

	for(int i = 1; i <= 16; i++){
		for(int ev = 0; ev < vtx->GetEntries(); ev++){
			vtx->GetEntry(ev);
			if((-160 + (i - 1)*20) <= Zric && Zric <= (-160 + i*20)){
				contaRic++;	  
            }                        
        }
        for(int ev = 0; ev < vrt->GetEntries(); ev++){
			vrt->GetEntry(ev);
			if((-160 + (i - 1)*20) <= z0 && z0 <= (-160 + i*20)){
				conta++;	  
            }                       
        }

		eff_v->SetTotalEvents(i , conta); 
        eff_v->SetPassedEvents(i, contaRic);

		contaRic = 0;
        conta = 0;
	}

//***********************************************************************************************************************

    TCanvas *c6 = new TCanvas("c6", "c6", 1000, 700);
	c6->SetLeftMargin(0.1);
	eff->SetMarkerStyle(20);
	eff->SetMarkerSize(1);
	eff->Draw("AP");
	c6->SetGrid();
	c6->SaveAs("Plot/c6.png");

	TCanvas *c7 = new TCanvas("c7", "c7", 1000, 700);
	c7->SetLeftMargin(0.1);	
	eff_v->SetMarkerStyle(20);
	eff_v->SetMarkerSize(1);
	eff_v->Draw("AP");
    c7->SetGrid();
    c7->SaveAs("Plot/c7.png");

	// eff->Write();
	// eff_v->Write();
	// fout.Close();

	// fin1->Close();
	// fin2->Close();
	// delete fin1;
	// delete fin2;
}
