
#include "include/Filipad.h"

void LoadData();
void compare();
void DrawSignal(int padid, int iPTT, int iPTA);
void DrawPP(int padID, int iPTT, int iPTA);
void DrawIAA(int padID, int iPTT, int iPTA);

double lowx=-1.6;
double highx=3.1;
double ly = 0.15;
double hy = 0.2;

TLatex latexRun;
TString strRun = "pp #sqrt{#it{s}} = 13 TeV";

const int Nsets = 1;
TString infiles[Nsets] = {
	"pythia8240_pp13TeV_Rope_1-2GeV.root"
};
TFile *fin[Nsets];

TString sLeg[Nsets] = {
	"PYTHIA Shoving g=4"
};

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};


// ROOT histograms
const int Nexp = 2;
enum exp { kCMS, kALICE, kEXP };
TString expName[Nexp]= {"CMS,2.0<|#Delta#eta|<4.0","ALICE, 1.5<|#Delta#eta|<1.8"};
const int NMultBins=4;
const double NchBin[NMultBins] = { 0, 35, 80, 105 };
TH1D *hMult[Nsets][Nexp]; // CMS vs ALICE
TH1D *hTriggPtBins[Nsets][Nexp][NMultBins]; // CMS vs ALICE
TH1D *hDeltaPhi[Nsets][Nexp][NMultBins];
TH1D *hDeltaPhi_Ratio[Nsets][Nexp][NMultBins];
TString hname ="";
TString htitle = "";
double ptmin=1.0;
double ptmax=2.0;

int iRef=0;

void Compare();
void DrawDPhi(int padID);

//------------------------------------------------------------------------------------------------
void LoadData() {

	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	for(int iS=0;iS<Nsets;iS++){
		for(int i = 0; i<Nexp; i++) {
			hname = Form("hMult%02d",i);
			htitle = Form("Multiplicity %s",expName[i].Data());
			hMult[iS][i] = (TH1D*)fin[iS]->Get(hname);
			for(int j = 0; j<NMultBins;j++){
				hname = Form("hTriggPtBinE%02dT%02d",i,j);
				htitle = Form("%.1f < p_{T} < %.1f, %.1f < N_{ch} < %.1f",ptmin,ptmax,NchBin[j],NchBin[j+1]);
				hTriggPtBins[iS][i][j] = (TH1D*)fin[iS]->Get(hname);
				hname = Form("hhDeltaPhiE%02dT%02d",i,j);
				hDeltaPhi[iS][i][j] = (TH1D*)fin[iS]->Get(hname);
				double ntrigg = hTriggPtBins[iS][i][j]->GetEntries();
				hDeltaPhi[iS][i][j]->Scale(1./ntrigg,"width");
			}
		}
	}
	// Ratio 
	for(int iS=0;iS<Nsets;iS++){
		for(int i = 0; i<Nexp; i++) {
			for(int j = 0; j<NMultBins;j++){
				hDeltaPhi_Ratio[iS][i][j] = (TH1D*)hDeltaPhi[iS][i][j]->Clone();
				hDeltaPhi_Ratio[iS][i][j]->Divide(hDeltaPhi[iS][kCMS][j]);
			}
		}
	}
}

//------------------------------------------------------------------------------------------------
void Compare(){
	LoadData();

	int ic=0;
	for(int ic=0; ic<NMultBins; ic++){
			DrawDPhi(ic++);
	}

}


//------------------------------------------------------------------------------------------------
void DrawDPhi(int padID) {
	Filipad *fpad[NMultBins];
	for(int ic=0;ic<NMultBins;ic++) {
		fpad[ic] = new Filipad(padID+ic+1, 1.1, 0.3, 100, 100, 0.7,NMultBins);
		fpad[ic]->Draw();
		//==== Upper pad
		TPad *p = fpad[ic]->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaPhi[0][kCMS][ic]->GetMaximum()*1.02;
		//ly = hDeltaPhi[0][kCMS][ic]->GetMinimum()*0.02;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "#Delta#phi", "1/N_{trigg} dN/d#Delta#phi",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.35,0.6,0.85,0.85,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.4, 0.9 ,strRun);
		leg->AddEntry((TObject*)NULL,sLeg[0],"");
		leg->AddEntry((TObject*)NULL,hDeltaPhi[0][kALICE][ic]->GetTitle(),"");

		for(int iS=0;iS<Nsets;iS++) {
			for(int i = 0; i<Nexp; i++) {
				hDeltaPhi[iS][i][ic]->SetMarkerStyle(gMarkers[i]);
				hDeltaPhi[iS][i][ic]->SetMarkerColor(gColors[i]);
				hDeltaPhi[iS][i][ic]->SetLineColor(gColors[i]);
				hDeltaPhi[iS][i][ic]->Draw("p,same");
				leg->AddEntry(hDeltaPhi[iS][i][ic],expName[i],"pl");
			}
		}


		leg->Draw();

		//==== Lower pad
		p = fpad[ic]->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.95, 1.05);
		hset( *hfr1, "#Delta#phi", "ALICE/CMS",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
		hfr1->Draw();
		for(int iS=0;iS<Nsets;iS++){
			hDeltaPhi_Ratio[iS][kALICE][ic]->SetMarkerStyle(20);
			hDeltaPhi_Ratio[iS][kALICE][ic]->Draw("psame");
		}
		
		gPad->GetCanvas()->SaveAs(Form("figs/DeltaPhi_C%02d.pdf",ic));
	}
}
