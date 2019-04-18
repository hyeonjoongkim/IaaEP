
void setStyle();

void LoadData();
void compare();

double lowx=-1.6;
double highx=3.1;
double ly = 0.15;
double hy = 0.2;

TLatex latexRun;
TString strRun = "pp #sqrt{#it{s}} = 13 TeV";

const int Nsets = 2;
TString infiles[Nsets] = {
	"pythia8240_pp13TeV_RopegAmplitude10.0_1-2GeV.root",
	"pythia8240_pp13TeV_RopegAmplitude10.0_2-3GeV.root"
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
TCanvas *canvas[NMultBins];
TCanvas *c1;
void Compare();
void DrawDPhi(int padID);
void DrawDPhiAll(int iE);

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
	//DrawDPhiAll(kCMS);
	DrawDPhiAll(kALICE);
	//for(int ic=0; ic<NMultBins; ic++){
	//		DrawDPhi(ic++);
	//}

}


//------------------------------------------------------------------------------------------------
void DrawDPhi(int padID) {
	for(int ic=0;ic<NMultBins;ic++) {
		//Legend definition
		canvas[ic] = new TCanvas(Form("c%d",ic+1), "Canvas");
		gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
		//setStyle();
		TLegend *leg = new TLegend(0.35,0.6,0.85,0.85,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
		leg->AddEntry((TObject*)NULL,sLeg[0],"");
		leg->AddEntry((TObject*)NULL,hDeltaPhi[0][kALICE][ic]->GetTitle(),"");
		for(int i = 0; i<Nexp; i++) {
				hDeltaPhi[0][i][ic]->SetMarkerStyle(gMarkers[i]);
				hDeltaPhi[0][i][ic]->SetMarkerColor(gColors[i]);
				hDeltaPhi[0][i][ic]->SetLineColor(gColors[i]);
				hDeltaPhi[0][i][ic]->Draw("p,same");
				leg->AddEntry(hDeltaPhi[0][i][ic],expName[i],"pl");
		}
		hDeltaPhi[0][kALICE][ic]->GetXaxis()->SetTitle("#Delta#phi");
		hDeltaPhi[0][kALICE][ic]->GetXaxis()->CenterTitle();
		auto rp = new TRatioPlot(hDeltaPhi[0][kALICE][ic], hDeltaPhi[0][kCMS][ic]);
		//canvas[ic]->SetTicks(0, 1);
   		rp->Draw();
   		
   		rp->GetLowYaxis()->SetNdivisions(505);
   		rp->GetLowerRefXaxis()->SetTitle("#Delta#phi");
   		rp->GetLowerRefYaxis()->SetTitle("ALICE/CMS");
   		rp->GetUpperRefYaxis()->SetTitle("1/N_{trigg} dN/d#Delta#phi");
   		rp->GetUpperRefXaxis()->CenterTitle();
   		rp->GetUpperRefYaxis()->CenterTitle();
   		rp->GetUpperPad()->cd(); 
   		latexRun.DrawLatexNDC( 0.4, 0.9 ,strRun);
   		leg->Draw();
   		canvas[ic]->Update();
   	}
}

//------------------------------------------------------------------------------------------------
void DrawDPhiAll(int iE) {
	c1 = new TCanvas("c1", "Canvas");
	gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.35,0.6,0.85,0.85,"","brNDC");
	leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
	leg->AddEntry((TObject*)NULL,sLeg[0],"");
	//leg->AddEntry((TObject*)NULL,hDeltaPhi[0][kALICE][0]->GetTitle(),"");
	for(int ic=0;ic<NMultBins;ic++) {
		hDeltaPhi[0][iE][ic]->SetMarkerStyle(gMarkers[ic]);
		hDeltaPhi[0][iE][ic]->SetMarkerColor(gColors[ic]);
		hDeltaPhi[0][iE][ic]->SetLineColor(gColors[ic]);
		hDeltaPhi[0][iE][ic]->GetXaxis()->SetTitle("#Delta#phi");
		hDeltaPhi[0][iE][ic]->GetXaxis()->CenterTitle();
		hDeltaPhi[0][iE][ic]->GetYaxis()->CenterTitle();
		hDeltaPhi[0][iE][ic]->GetYaxis()->SetTitle("1/N_{trigg} dN/d#Delta#phi");
		hDeltaPhi[0][iE][ic]->Draw("p,same");
		leg->AddEntry(hDeltaPhi[0][iE][ic],hDeltaPhi[0][iE][ic]->GetTitle(),"pl");
   	}
   		//rp->GetUpperPad()->cd(); 
   	latexRun.DrawLatexNDC( 0.4, 0.9 ,strRun);
   	leg->Draw();
   	c1->Update();

}

void setStyle(){
        gStyle->SetTitleFont(22,"y");
        gStyle->SetTitleFont(22,"z");
        gStyle->SetLabelFont(22,"x");
        gStyle->SetLabelFont(22,"y");
        gStyle->SetLabelFont(22,"z");
        gStyle->SetTitleSize(0.06,"x");
        gStyle->SetTitleSize(0.06,"y");
        gStyle->SetLabelSize(0.05,"x");
        gStyle->SetLabelSize(0.05,"y");
        gStyle->SetTitleOffset(0.9,"x");
        gStyle->SetTitleOffset(0.8,"y");
        gStyle->SetNdivisions(505);
        gStyle->SetOptStat(0);
        // gStyle->SetOptTitle(0);
        gStyle->SetOptFit(1);
}
