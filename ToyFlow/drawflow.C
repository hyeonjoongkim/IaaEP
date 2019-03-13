
#define NC 6             //Changed for high pT test
static const char *presn[] = {"V0A","V0C","V0P"};

enum RESOLUTION{
	R_V0A,
	R_V0C,
	R_V0P,
	R_COUNT
};

static double CentBins[NC+3] = {0,5,10,20,30,40,50,60,70};

void drawflow() {

	double pi = TMath::Pi();
	string name;
	TFile *fin = new TFile("results.root");                     // Resolution & smearing histos
	TH1D *pah[R_COUNT][NC], *pbh[R_COUNT][NC], *pch[R_COUNT][NC];
	TH1D *resolution[R_COUNT];
	TH2D *evpcorr2d[R_COUNT][NC];

	TH1D *evpdifference[R_COUNT][NC];
	for (int s = 0; s < R_COUNT; s++) {
		resolution[s] = new TH1D(Form("resolution_%s", presn[s]), Form("resolution_%s", presn[s]), NC, CentBins);
		for (int c = 0; c < NC; c++) {
			pah[s][c]=(TH1D*)fin->Get(Form("h_res_%s_a%02u",presn[s],c));
			pbh[s][c]=(TH1D*)fin->Get(Form("h_res_%s_b%02u",presn[s],c));
			pch[s][c]=(TH1D*)fin->Get(Form("h_res_%s_c%02u",presn[s],c));

			evpdifference[s][c] = (TH1D*) fin->Get(Form("h_evpdiff_%s_c%02d", presn[s], c));
			evpdifference[s][c]->Print();
			evpdifference[s][c]->GetXaxis()->SetTitle("Event Plane Difference #Psi_{true}-#Psi{reco}");
			evpdifference[s][c]->GetYaxis()->SetTitleOffset(1.3);
			evpdifference[s][c]->SetTitle("");

			evpcorr2d[s][c] = (TH2D*)fin->Get(Form("h_evpcorr2d_%s_c%02u", presn[s], c));
			resolution[s]->SetBinContent(c+1, TMath::Sqrt(pah[s][c]->GetMean()*pbh[s][c]->GetMean()/pch[s][c]->GetMean()));
			evpcorr2d[s][c]->GetXaxis()->SetRangeUser(-0.5*pi, 0.5*pi);
			evpcorr2d[s][c]->GetYaxis()->SetRangeUser(-0.5*pi, 0.5*pi);
		}
	}

	TCanvas *c1 = new TCanvas();
	/* Drawing and Writing to file */
	for (int s = 0; s < R_COUNT; s++) {
		resolution[s]->Write();

		resolution[s]->Draw();
		c1->SaveAs(Form("figs/Resolution%s.pdf", presn[s]));
		for (int c = 0; c < NC; c++) {
			evpdifference[s][c]->SetStats(0);

			TLatex l2;
			l2.SetTextSize(0.04);
			evpdifference[s][c]->Write();
			evpdifference[s][c]->SetMaximum(evpdifference[s][c]->GetMaximum() * 1.2);
			evpdifference[s][c]->Draw();
			l2.DrawLatexNDC(0.12, 0.85, Form("Toy MC, #sqrt{s_{NN}}=5.02TeV, Centrality %.0f-%.0f%%", CentBins[c], CentBins[c+1]));
			l2.DrawLatexNDC(0.12, 0.8, Form("Event plane reconstructed with V0C Detector"));
			c1->SaveAs(Form("figs/evpdifference%s_%d.pdf", presn[s], c));

			evpcorr2d[s][c]->SetStats(0);
			evpcorr2d[s][c]->Draw("colz");
			c1->SaveAs(Form("figs/%s_evpcorr2d_%d.pdf", presn[s], c));
		}
	}

}
