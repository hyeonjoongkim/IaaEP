
#include "include/Filipad.h"

void LoadData();
void compare();
void DrawSignal(int padid, int iPTT, int iPTA);
void DrawPP(int padID, int iPTT, int iPTA);
void DrawIAA(int padID, int iPTT, int iPTA);

double lowx=-0.8;
double highx=0.85;
double ly = -0.05;
double hy = 0.3;
double lowIAA = -0.2;
double highIAA = 2.2;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV";

const int Nsets = 5;
TString infiles[Nsets] = {
	"sysErrors/Signal_LHC10h_AOD86_MgFpMgFm_5217_JCIAA_TPCOnly_H0_T0_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3a_JCIAA_EPInclusive_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_5175_JCIAA_TPCOnly_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_5146_JCIAA_GlobalSDD_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
};
TFile *fin[Nsets];

TString sLeg[Nsets] = {
	"LHC10h",
	"AMPT String melting",
	"AMPT String melting w/o hardronic rescattering",
	"15o, TPCOnly",
	"15o, GlobalSDD"
};

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
enum dataType { AA, pp };


TH1D *hDeltaEtaSig[Nsets][2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
TH1D *hIAADeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; // background substracted signal IAA

TH1D *hRatio_DeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; //
TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef=0;

//------------------------------------------------------------------------------------------------
void LoadData() {
	
	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 0;
	TriggPtBorders             = (TVector*) fin[irefD]->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin[irefD]->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin[irefD]->Get("CentBinBorders");
	NumCent[AA]    = CentBinBorders->GetNoElements()-2; // for 5TeV one less
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
	// Calculation Various Ratios
	// 1. Out/In ratio for each data setfor(int i=0;i<Nsets;i++){
	for(int ic=0; ic<NumCent[AA]; ic++){
		for(int iptt=0; iptt<NPTT; iptt++){
			for(int ipta=0;ipta<NPTA;ipta++) {
				for(int i=0;i<Nsets;i++){
					hRatio_DeltaEtaSig[i][ic][iptt][ipta] = (TH1D*)hDeltaEtaSig[i][AA][ic][iptt][ipta]->Clone();
					hRatio_DeltaEtaSig[i][ic][iptt][ipta]->Divide(hDeltaEtaSig[iRef][AA][ic][iptt][ipta]);
				}
			} // ipta
		} // iptt 
	} // ic
}

//------------------------------------------------------------------------------------------------
void Compare(){
	LoadData();
	
	int ic=0;
	for(int iptt=3; iptt<NPTT; iptt++){
		for(int ipta=1;ipta<NPTA;ipta++) {
			DrawSignal(ic++, iptt, ipta);
			//DrawIAA(ic++, iptt, ipta);
	 	}
	}

}


//------------------------------------------------------------------------------------------------
void DrawSignal(int padID, int iPTT, int iPTA) {
	Filipad *fpad[NumCent[AA]];
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++) {
		fpad[ic] = new Filipad(padID+ic+1, 1.1, 0.5, 100, 100, 0.7,NumCent[AA]);
		fpad[ic]->Draw();
		//==== Upper pad
		TPad *p = fpad[ic]->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEtaSig[2][AA][ic][iPTT][iPTA]->GetMaximum()*1.5;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.35,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hDeltaEtaSig[0][AA][ic][iPTT][iPTA]->GetTitle(),"");

		for(int iS=0;iS<Nsets;iS++) {
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->SetMarkerColor(gColors[iS]);
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->SetLineColor(gColors[iS]);
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hDeltaEtaSig[iS][AA][ic][iPTT][iPTA],sLeg[iS],"pl");
		}

		
		leg->Draw();

		//==== Lower pad
		p = fpad[ic]->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", Form("Ratio to %s",sLeg[iRef].Data()),1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		for(int i=0;i<Nsets;i++) {
			hRatio_DeltaEtaSig[i][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
			hRatio_DeltaEtaSig[i][ic][iPTT][iPTA]->SetMarkerColor(gColors[i]);
			hRatio_DeltaEtaSig[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
			hRatio_DeltaEtaSig[i][ic][iPTT][iPTA]->Draw("p,same");
		}
		gPad->GetCanvas()->SaveAs(Form("figs_iaa/DeltaEta_BeamEnergy_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
	//for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];
}

//------------------------------------------------------------------------------------------------
void DrawIAA(int padID, int iPTT, int iPTA) {
	Filipad *fpad[NumCent[AA]];
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++) {
		fpad[ic] = new Filipad(padID+ic+1, 1.1, 0.5, 100, 100, 0.7,NumCent[AA]);
		fpad[ic]->Draw();
		//==== Upper pad
		TPad *p = fpad[ic]->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hIAADeltaEtaSig[2][ic][iPTT][iPTA]->GetMaximum()*1.2;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, lowIAA, highIAA); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "I_{AA}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hIAADeltaEtaSig[0][ic][iPTT][iPTA]->GetTitle(),"");

		for(int iS=0;iS<Nsets;iS++) {
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetLineColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hIAADeltaEtaSig[iS][ic][iPTT][iPTA],sLeg[iS],"pl");
		}

		
		leg->Draw();

		//==== Lower pad
		p = fpad[ic]->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", Form("Ratio to %s",sLeg[iRef].Data()),1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		//for(int i=0;i<Nsets;i++) {
		//	hRatio_IAA[i][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
		//	hRatio_IAA[i][ic][iPTT][iPTA]->SetMarkerColor(gColors[i]);
		//	hRatio_IAA[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
		//	hRatio_IAA[i][ic][iPTT][iPTA]->Draw("p,same");
		//}
		gPad->GetCanvas()->SaveAs(Form("figs_iaa/IAA_BeamEnergy_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
	//for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];
}


//------------------------------------------------------------------------------------------------
void DrawPP(int padID, int iPTT, int iPTA) {
	Filipad *fpad;
	lowx = -0.01;
	fpad = new Filipad(padID+1, 1.1, 0.5, 100, 100, 0.7,5);
	fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEtaSig[iRef][pp][0][iPTT][iPTA]->GetMaximum()*1.2;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hDeltaEtaSig[0][pp][0][iPTT][iPTA]->GetTitle(),"");

		
		hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[0]);
		hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[0]);
		hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetLineColor(gColors[0]);
		hDeltaEtaSig[0][pp][0][iPTT][iPTA]->Draw("p,same");
		leg->AddEntry(hDeltaEtaSig[0][pp][0][iPTT][iPTA],"Data pp","pl");
		hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[1]);
		hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[1]);
		hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetLineColor(gColors[1]);
		hDeltaEtaSig[1][pp][0][iPTT][iPTA]->Draw("p,same");
		leg->AddEntry(hDeltaEtaSig[1][pp][0][iPTT][iPTA],"pythia pp","pl");
		
		leg->Draw();

		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();

		TH1D *hratio = (TH1D*)hDeltaEtaSig[1][pp][0][iPTT][iPTA]->Clone();
		hratio->Divide(hDeltaEtaSig[0][pp][0][iPTT][iPTA]);
		hratio->Draw("p,same");
		//gPad->GetCanvas()->SaveAs(Form("figs/DeltaEta_OUTOIN_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
