
#include "include/Filipad.h"

void LoadData();
void compare();
void DrawPP(int padID, int iPTT, int iPTA);
void DrawGF(int padID, int iPTT, int iPTA);

double lowx=-0.1;
double highx=0.6;
double ly = -0.05;
double hy = 0.3;
double lowIAA = -0.2;
double highIAA = 2.2;

TLatex latexRun;
TString strRun = "PYTHIA #sqrt{#it{s}} = 2.76 TeV";

const int Nsets = 4;
TString infiles[Nsets] = {
	"sysErrors/Signal_LHC10h_AOD86_MgFpMgFm_5217_JCIAA_TPCOnly_H0_T0_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF1_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_QF1_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
	//"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	//"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_LHC12f1b_Phojet_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	//"sysErrors/Signal_LHC15o_GlobalSDD_JCIAA_GlobalSDD_pythia8230_pp5.02TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	//"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_LHC17l3b_cent_woSDD_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
};
TFile *fin[Nsets];

TString sLeg[Nsets] = {
	"data",
	"SoftQCD",
	"SoftQCD, Gluon Filtering",
	"SoftQCD, Quark Filtering"
	//"LHC12f1a_Pythia",
	//"LHC12f1b_Phojet",
	//"pythia8230 SoftQCD, #sqrt{#it{s}} = 5.02 TeV",
	//"LHC17l3b, #sqrt{#it{s}} = 5.02 TeV"
};

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
enum dataType { AA, pp };


TH1D *hDeltaEtaSig[Nsets][2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
TH1D *hIAADeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; // background substracted signal IAA

TH1D *hRatiosPP[Nsets];
TH1D *hRatiosFilter[2]; //GF->0 QF->1

TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef=1; // data:

//------------------------------------------------------------------------------------------------
void LoadData() {
	
	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 0;
	TriggPtBorders             = (TVector*) fin[irefD]->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin[irefD]->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin[irefD]->Get("CentBinBorders");
	NumCent[AA]    = CentBinBorders->GetNoElements()-1;
	NumCent[pp]    = 1; 
	NPTT     = TriggPtBorders->GetNoElements()-1;
	NPTA     = AssocPtBorders->GetNoElements()-1;
	cout <<"PbPb"<<endl;
	cout <<"bins:  c="<<  NumCent[AA] <<" ptt="<< NPTT <<" pta="<< NPTA  << endl; 
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	for(int i=0;i<Nsets;i++){
		for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
			for(int ic=0; ic<NumCent[idtyp]; ic++){
				for(int iptt=0; iptt<NPTT; iptt++){
					for(int ipta=0;ipta<NPTA;ipta++) {
						hDeltaEtaSig[i][idtyp][ic][iptt][ipta] = (TH1D *)fin[i]->Get(Form("hDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
						if(idtyp==AA) hIAADeltaEtaSig[i][ic][iptt][ipta] = (TH1D *)fin[i]->Get(Form("hIAADeltaEtaSigC%02dT%02dA%02d",ic,iptt,ipta));
					} // ipta
				} // iptt 
			} // ic
		} // pp or AA
	} // iset

}

//------------------------------------------------------------------------------------------------
void Compare(){
	LoadData();	
	int ic=0;
	for(int iptt=3; iptt<NPTT; iptt++){
		for(int ipta=1;ipta<NPTA;ipta++) {
			//DrawPP(ic++,iptt,ipta);
			DrawGF(ic++,iptt,ipta);
		}
	}

}


//------------------------------------------------------------------------------------------------
void DrawPP(int padID, int iPTT, int iPTA) {
	Filipad *fpad;
	lowx = -0.01;
	fpad = new Filipad(padID+1, 1.1, 0.5, 100, 100, 0.7,NPTA);
	fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEtaSig[iRef][pp][0][iPTT][iPTA]->GetMaximum()*2.0;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hDeltaEtaSig[0][pp][0][iPTT][iPTA]->GetTitle(),"");

		for(int i=0;i<Nsets;i++){
			hDeltaEtaSig[i][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
			hDeltaEtaSig[i][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[i]);
			hDeltaEtaSig[i][pp][0][iPTT][iPTA]->SetLineColor(gColors[i]);
			hDeltaEtaSig[i][pp][0][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hDeltaEtaSig[i][pp][0][iPTT][iPTA],Form("%s",sLeg[i].Data()),"pl");
		}
		
		leg->Draw();
			// Calculation Various Ratios
		for(int i=0;i<Nsets;i++){
			hRatiosPP[i] = (TH1D*)hDeltaEtaSig[i][pp][0][iPTT][iPTA]->Clone();
			hRatiosPP[i]->Divide(hDeltaEtaSig[iRef][pp][0][iPTT][iPTA]);
		}
		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		for(int i=0;i<Nsets;i++){
			hRatiosPP[i]->Draw("p,same");
		}
		gPad->GetCanvas()->SaveAs(Form("figs_iaa/DeltaEta_ppComparison_T%02dA%02d.pdf",iPTT,iPTA));
	}

//------------------------------------------------------------------------------------------------
void DrawGF(int padID, int iPTT, int iPTA) {
	Filipad *fpad;
	lowx = -0.01;
	fpad = new Filipad(padID+1, 1.1, 0.5, 100, 100, 0.7,NPTA);
	fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEtaSig[3][pp][0][iPTT][iPTA]->GetMaximum()*1.3;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.5,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.35, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hDeltaEtaSig[0][pp][0][iPTT][iPTA]->GetTitle(),"");

		int sIndex[3] = {1,2,3};
		for(int i=0;i<3;i++){
			hDeltaEtaSig[sIndex[i]][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
			hDeltaEtaSig[sIndex[i]][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[i]);
			hDeltaEtaSig[sIndex[i]][pp][0][iPTT][iPTA]->SetLineColor(gColors[i]);
			hDeltaEtaSig[sIndex[i]][pp][0][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hDeltaEtaSig[sIndex[i]][pp][0][iPTT][iPTA],Form("%s",sLeg[sIndex[i]].Data()),"pl");
		}
		
		leg->Draw();
			// Calculation Various Ratios
		hRatiosFilter[0] = (TH1D*)hDeltaEtaSig[sIndex[1]][pp][0][iPTT][iPTA]->Clone();
		hRatiosFilter[0]->Divide(hDeltaEtaSig[sIndex[0]][pp][0][iPTT][iPTA]);
		hRatiosFilter[1] = (TH1D*)hDeltaEtaSig[sIndex[2]][pp][0][iPTT][iPTA]->Clone();
		hRatiosFilter[1]->Divide(hDeltaEtaSig[sIndex[0]][pp][0][iPTT][iPTA]);
		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		
		hRatiosFilter[0]->Draw("p,same");
		hRatiosFilter[1]->Draw("p,same");
		gPad->GetCanvas()->SaveAs(Form("figs_iaa/Pythia8_DeltaEta_GF_T%02dA%02d.pdf",iPTT,iPTA));
	}
