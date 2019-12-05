
#include "include/Filipad.h"

void LoadData();
void compare();
void DrawSystematics(int padid, int iPTT, int iPTA);
void CalculateRatios(int padID, int iPTT, int iPTA);
void DrawSignal(int padid, int iPTT, int iPTA, int iSet);
void DrawPP(int padID, int iPTT, int iPTA);
void DrawIAA(int padID, int iPTT, int iPTA, int iSet);

double lowx = -0.8;
double highx = 0.4;
double ly = -0.05;
double hy = 0.3;
double lowIAA = -0.2;
double highIAA = 2.2;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV";

enum dataType
{
  AA,
  pp
};
const int kMAXD = 20; // maximal number of pT trigger bins
const int kCENT = 10; // maximal number of pT trigger bins
const int kMAXSys = 5;
const int Nsets = 5; // Number of systematic sets
const int Nsysfile[Nsets] = {
    1, 2, 1, 1, 2 };                    // Variable set inside each systematics (usually 2 file)
const int isp2p[Nsets] = {0, 0, 0, 0, 0}; // If it's p2p, put 1, else put 0

TString defaultfile = "sysErrors/Signal_AOD160_TPConly_JDiHadronIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root";
TString systematicfile = "Systematic.root";

TH1D *hSystematic_pbpFull[2][kCENT][kMAXD][kMAXD];
TH1D *hSystematic_scaFull[2][kCENT][kMAXD][kMAXD];
TH1D *hDeltaEtaSig_def[2][kCENT][kMAXD][kMAXD]; // Default background substracted
                                                // signal based on fit
TH1D *hIAADeltaEtaSig_def[kCENT][kMAXD][kMAXD]; // Default background substracted signal IAA

TGraphErrors *gscaFull[kCENT][kMAXD][kMAXD];
TGraphErrors *gpbpFull[kCENT][kMAXD][kMAXD];

TFile *fdefault;
TFile *fsystem;

TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef = 0;
float *ptaborders;
float *pttborders;

void LoadData() {
    fdefault = TFile::Open(defaultfile);
    fsystem = TFile::Open(systematicfile);

  int irefD = 0;
  TriggPtBorders = (TVector *)fdefault->Get("TriggPtBorders");
  AssocPtBorders = (TVector *)fdefault->Get("AssocPtBorders");
  CentBinBorders = (TVector *)fdefault->Get("CentBinBorders");
  NumCent[AA] = CentBinBorders->GetNoElements() - 2; // for 5TeV one less
  NumCent[pp] = 1;
  NPTT = TriggPtBorders->GetNoElements() - 1;
  NPTA = AssocPtBorders->GetNoElements() - 1;
  pttborders = TriggPtBorders->GetMatrixArray();
  ptaborders = AssocPtBorders->GetMatrixArray();
  cout << "PbPb" << endl;
  cout << "bins:  c=" << NumCent[AA] << " ptt=" << NPTT << " pta=" << NPTA
       << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  //------------ R e a d    D a t a ------------
  cout << "Reading data...." << endl;
  for (int idtyp = 0; idtyp < 2;
       idtyp++) { // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
    for (int ic = 0; ic < NumCent[idtyp]; ic++) {
      for (int iptt = 0; iptt < NPTT; iptt++) {
        for (int ipta = 0; ipta < NPTA; ipta++) {
          if (pttborders[iptt] - ptaborders[ipta] < 0.000001)
            continue;
          hDeltaEtaSig_def[idtyp][ic][iptt][ipta] = (TH1D *)fdefault->Get(
              Form("hDeltaEtaSig%02dC%02dT%02dA%02d", idtyp, ic, iptt, ipta));
          if (idtyp == AA)
            hIAADeltaEtaSig_def[ic][iptt][ipta] = (TH1D *)fdefault->Get(
                Form("hIAADeltaEtaSigC%02dT%02dA%02d", ic, iptt, ipta));
          hSystematic_scaFull[idtyp][ic][iptt][ipta] = (TH1D*)fsystem->Get(Form("hSystematic_scaFull%02d%02d%02d%02d", idtyp, ic, iptt, ipta));
          hSystematic_pbpFull[idtyp][ic][iptt][ipta] = (TH1D*)fsystem->Get(Form("hSystematic_pbpFull%02d%02d%02d%02d", idtyp, ic, iptt, ipta));
          if (idtyp == 0) {
              gscaFull[ic][iptt][ipta] = new TGraphErrors();
              gpbpFull[ic][iptt][ipta] = new TGraphErrors();

              gscaFull[ic][iptt][ipta]->SetName(Form("gscaFull%02d%02d%02d", ic, iptt, ipta));
              gpbpFull[ic][iptt][ipta]->SetName(Form("gpbpFull%02d%02d%02d", ic, iptt, ipta));

          }
        }
      }
    }
  }
  for (int ic = 0; ic < NumCent[0]; ic++) {
    for (int iptt = 0; iptt < NPTT; iptt++) {
      for (int ipta = 0; ipta < NPTA; ipta++) {
        if (pttborders[iptt] - ptaborders[ipta] < 0.000001)
          continue;
        double scalingerr = 0;
        { // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
          scalingerr +=
              hSystematic_scaFull[0][ic][iptt][ipta]->GetBinContent(1) *
              hSystematic_scaFull[0][ic][iptt][ipta]->GetBinContent(1); // AA
          scalingerr += hSystematic_scaFull[1][0][iptt][ipta]->GetBinContent(1) *
              hSystematic_scaFull[1][0][iptt][ipta]->GetBinContent(1); // pp
        } // idtyp
        scalingerr = sqrt(scalingerr)*0.01;

        gscaFull[ic][iptt][ipta]->SetPoint(0, 0.2, 1);
        gscaFull[ic][iptt][ipta]->SetPointError(0, 0.5, scalingerr);
        int Nbins = hIAADeltaEtaSig_def[ic][iptt][ipta]->GetNbinsX();
        for (int ibin = 1; ibin < Nbins; ++ibin) {
          double pbperr = 0;

          {
            pbperr +=
                hSystematic_pbpFull[0][ic][iptt][ipta]->GetBinContent(ibin) *
                hSystematic_pbpFull[0][ic][iptt][ipta]->GetBinContent(ibin);
            pbperr +=
                hSystematic_pbpFull[1][0][iptt][ipta]->GetBinContent(ibin) *
                hSystematic_pbpFull[1][0][iptt][ipta]->GetBinContent(ibin);
          }
          pbperr = sqrt(pbperr);
          double xpoint =
              hIAADeltaEtaSig_def[ic][iptt][ipta]->GetBinCenter(ibin);
          double ypoint =
              hIAADeltaEtaSig_def[ic][iptt][ipta]->GetBinContent(ibin);
          double xerr =
              hIAADeltaEtaSig_def[ic][iptt][ipta]->GetBinWidth(ibin) / 2;

          gpbpFull[ic][iptt][ipta]->SetPoint(ibin - 1, xpoint, ypoint);
          gpbpFull[ic][iptt][ipta]->SetPointError(ibin - 1, xerr,
                                                  ypoint * pbperr * 0.01);
        cout << ic << " " << iptt << " " << ipta << " scal : " << pbperr << endl;

                                                  

        } // ibin

        gscaFull[ic][iptt][ipta]->SetFillColorAlpha(kGray, 0.5);
        gpbpFull[ic][iptt][ipta]->SetFillColorAlpha(kRed+2, 0.5);

      }   // ipta
    }
  }

} // Load Data

void DrawSystematics() {
LoadData();
TCanvas *fpad[kCENT][kMAXD][kMAXD];
lowx = -0.01;

  for (int ic = 0; ic < 1; ++ic) {
    for (int iptt = 3; iptt < NPTT; iptt++) {
      for (int ipta = 0; ipta < NPTA; ipta++) {
        if (pttborders[iptt] - ptaborders[ipta] < 0.000001)
          continue;
          fpad[ic][iptt][ipta] = new TCanvas(Form("canvas%02d%02d%02d", ic, iptt, ipta), Form("canvas%02d%02d%02d", ic, iptt, ipta), 800, 600);
          fpad[ic][iptt][ipta]->Draw();
          TPad *p = (TPad*)fpad[ic][iptt][ipta]->GetPad(0);
          p->SetTickx();
          p->cd();
          //    hy = hIAADeltaEtaSig[2][ic][iPTT][iPTA]->GetMaximum() * 1.2;
          TH2F *hfr =
              new TH2F("hfr", " ", 100, lowx, highx, 10, lowIAA,
                       highIAA); // numbers: tics x, low limit x, upper limit
                                 // x, tics y, low limit y, upper limit y
          hset(*hfr, "|#Delta#eta|", "I_{AA}", 1.1, 1.0, 0.09, 0.09, 0.01, 0.01,
               0.04, 0.05, 510,
               505); // settings of the upper pad: x-axis, y-axis
          hfr->Draw();
          hIAADeltaEtaSig_def[ic][iptt][ipta]->SetMarkerStyle(20);
//          hIAADeltaEtaSig_def[ic][iptt][ipta]->Draw("same");
          gscaFull[ic][iptt][ipta]->Draw("2");
          gpbpFull[ic][iptt][ipta]->Draw("2");
          
      }
    }
  }
}
