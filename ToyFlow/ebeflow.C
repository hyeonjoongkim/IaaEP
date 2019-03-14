#include <TMath.h>

void ebeflow() {
	TFile *fin = TFile::Open("results.root");
	const int NC = 5;
	TH1D *hphi[NC];
	for(int i = 0; i<NC ; i++) {
		TString hname = Form("hPhiEbE_c%02d",i);
		cout << hname << endl;
		hphi[i] = (TH1D*)fin->Get(hname);
	}
	hphi[0]->Draw();

	TCanvas *cFlow =new TCanvas("cFlow","3",1000,1000);
	cFlow-> Divide(2,2,0,0);


	cFlow->cd(1);
	hphi[0]->Draw();

	cFlow->cd(2);
	hphi[1]->Draw();

	cFlow->cd(3);
	hphi[2]->Draw();

	cFlow->cd(4);
	hphi[3]->Draw();			

}