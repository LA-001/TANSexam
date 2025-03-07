#include "Simulazione1.h"
#include "Simulazione2.h"
#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TLeaf.h>

using namespace std;

void Simulazione(int numero, unsigned int seed){

    TFile *file = TFile::Open("fileRoot/kinem.root");
    TH1F *hm = (TH1F*)file->Get("hm");
    TH1F *eta = (TH1F*)file->Get("heta2");
	TFile fout("fileRoot/simulazione.root","recreate");

    Simulazione1 *ptr = new Simulazione1(eta,hm,seed);
    Simulazione2 *ptr2 = new Simulazione2(seed);

    TNtuple *nt = new TNtuple("nt","generazione vertice","x0:y0:z0:phi:tetha");
    float xnt[5]; 
    TNtuple *hit = new TNtuple("hit","hit layer","x2:y2:z2:x3:y3:z3:etichetta");
    float xhit[7];
	TNtuple *vrt = new TNtuple("vrt","verità montecarlo","z0:moltiplicita");
    float xvrt[2]; 

    float x0, y0, z0, molti, pseudorap;

    //Generazione dei vertici e delle molteplicità

    for(int tot = 1; tot<= numero; tot++){	
	
        nt->Reset();  

        x0 = ptr->VertexSimXY();  
        y0 = ptr->VertexSimXY();  
        z0 = ptr->VertexSimZ(true);   
        molti = ptr->Multiplicity(true,1); 

        for (int i = 0; i < molti; i++) {
                xnt[0] = x0; 
                xnt[1] = y0; 
                xnt[2] = z0; 
                xnt[3] = ptr->Azimut();

                pseudorap = ptr->Eta();   
                xnt[4] = 2*TMath::ATan(TMath::Exp(-pseudorap));
                nt->Fill(xnt);
                
         }
		
		xvrt[0] = z0;
		xvrt[1] = molti;
		vrt->Fill(xvrt);

        //trasporto delle particelle generate

        for(int ev=0; ev<nt->GetEntries(); ev++){
            nt->GetEntry(ev);

            float xx0 = nt->GetLeaf("x0")->GetValue();
            float yy0 = nt->GetLeaf("y0")->GetValue();
            float zz0 = nt->GetLeaf("z0")->GetValue();
            float pphi = nt->GetLeaf("phi")->GetValue();
            float ttetha = nt->GetLeaf("tetha")->GetValue();

            vector<float> punto1 = ptr2->EquazioneRetta1(xx0, yy0, zz0, pphi, ttetha);  //Beam Pipe

            vector<float> scatter1 = ptr2->Scattering1(ttetha,pphi,true); // scattering Beam Pipe

            vector<float> punto2 = ptr2->EquazioneRetta2(punto1[0], punto1[1], punto1[2], scatter1[0], scatter1[1], scatter1[2], 2);  //L1

            if(-135. <= punto2[2] && punto2[2] <= 135.){

                vector<float> scatter2 = ptr2->Scattering2(scatter1[0], scatter1[1], scatter1[2], true); // scattering L1

                vector<float> punto3 = ptr2->EquazioneRetta2(punto2[0], punto2[1], punto2[2], scatter2[0], scatter2[1], scatter2[2], 3);  //L2

                if(-135. <= punto3[2] && punto3[2] <= 135.){
                  
                    xhit[0] = punto2[0];
                    xhit[1] = punto2[1];
                    xhit[2] = punto2[2];
                    xhit[3] = punto3[0];
                    xhit[4] = punto3[1];
                    xhit[5] = punto3[2];

                    xhit[6] = tot;  //etichetta dello scontro

                    hit->Fill(xhit);
                }else{          //hit nel primo layer ma non nel secondo
                    
                    xhit[0] = punto2[0];
                    xhit[1] = punto2[1];
                    xhit[2] = punto2[2];
                    xhit[3] = punto3[0];
                    xhit[4] = punto3[1];
                    xhit[5] = 10000.;

                    xhit[6] = tot;  //etichetta dello scontro

                    hit->Fill(xhit);
                }
            } 
        }
	}


hit->Write();
vrt->Write();
fout.Close();

delete ptr;
delete ptr2;
}
