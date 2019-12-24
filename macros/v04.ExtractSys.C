#include "include/Filipad.h"

void LoadData();
void compare();
void DrawSystematics(int padid, int iPTT, int iPTA);
void CalculateRatios(int padID, int iPTT, int iPTA);
void DrawSignal(int padid, int iPTT, int iPTA, int iSet);
void DrawPP(int padID, int iPTT, int iPTA);
void DrawIAA(int padID, int iPTT, int iPTA, int iSet);

double lowx = -0.8;
double highx = 0.85;
double ly = -0.05;
double hy = 0.3;
double lowIAA = -0.2;
double highIAA = 2.2;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV";

const int kMAXSys = 6;
const int Nsets = 6; // Number of systematic sets
const int Nsysfile[Nsets] = {
    1, 2, 1, 1, 2, 1 };                    // Variable set inside each systematics (usually 2 file)
const int isp2p[Nsets] = {0, 0, 0, 0, 0, 1}; // If it's p2p, put 1, else put 0

TString defaultfile = "sysErrors/Signal_AOD86_TPCOnly_NTM_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_NTM_Iaa_R0.2_1.0_1.60_Near_Wing0.root";
TString outfile = "Systematic.root";

TString infiles[Nsets][3] = { // Nsets / Nsysfile
    { //"sysErrors/Signal_AOD86_RAA_JDiHadronIAA_RAA_H0_T0_LHC11a_AOD113_noSDD_RAA_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
    "sysErrors/Signal_AOD86_RAA_NTM_JDiHadronIAA_RAA_H0_T0_LHC11a_AOD113_noSDD_TPC_NTM_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
    }, // Track Cut
    {
      "sysErrors/Signal_AOD86_TPCOnly_NTM_phi15_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_NTM_phi15_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
"sysErrors/Signal_AOD86_TPCOnly_NTM_phi25_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_NTM_phi25_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
    }, // dphi

    {"sysErrors/Signal_AOD86_TPCOnly_NTM_zvtx8_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_NTM_zvtx8_Iaa_R0.2_1.0_1.60_Near_Wing0.root",

}, // zvtx
{
  "sysErrors/Signal_AOD86_TPCOnly_NTM_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",

}, // Track Merging
{
  "sysErrors/Signal_FR1.4_AOD86_TPCOnly_NTM_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_NTM_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
"sysErrors/Signal_FR1.6_AOD86_TPCOnly_NTM_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_NTM_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
}, // fit range
{
"sysErrors/Signal_gaus_AOD86_TPCOnly_NTM_JCIAA_TPCOnly_H0_T0_LHC11a_AOD113_noSDD_GlobalSDD_NTM_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
} //gaussian

//        "sysErrors/Signal_LHC10h_AOD86_MgFpMgFm_5217_JDiHadronIAA_TPCOnly_H0_T0_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
//        "sysErrors/Signal_LHC10h_AOD86_MgFpMgFm_5217_JDiHadronIAA_TPCOnly_H0_T0_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",

};

TFile *fin[Nsets][kMAXSys];
TFile *fdefault;
TFile *foutfile;

TString sLeg[Nsets][3] = { {"TPCOnly&GlobalSDD", "RAA&TPCOnly"}, {"DPhi 2.0", "DPhi1.5", "DPhi2.5"}, {"ZVtx10", "Zvtx8"}, {"NoTrackMerging", "TrackMerging on pp"}, {"FitRange1.5", "FitRange1.4", "FitRange1.6"},{"Fit GG", "Fit Gaus" } 
};

int gMarkers[] = {20, 24, 21, 25, 23, 27, 29, 30};
int gColors[] = {kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};

const int kMAXD = 20; // maximal number of pT trigger bins
const int kCENT = 10; // maximal number of pT trigger bins
enum dataType
{
  AA,
  pp
};

TH1D *hDeltaEtaSig_def[2][kCENT][kMAXD][kMAXD]; // Default background substracted
                                                // signal based on fit
TH1D *hIAADeltaEtaSig_def[kCENT][kMAXD][kMAXD]; // Default background substracted signal IAA

TH1D *hDeltaEtaSig[Nsets][kMAXSys][2][kCENT][kMAXD]
                  [kMAXD]; // background substracted signal based on fit
TH1D *hIAADeltaEtaSig[Nsets][kMAXSys][kCENT][kMAXD]
                     [kMAXD]; // background substracted signal IAA

TH1D *hRatio_DeltaEtaSig[Nsets][kMAXSys][2][kCENT][kMAXD]
                        [kMAXD]; // Ratio between default - and comparison
TH1D *hRatio_IAADeltaEtaSig[Nsets][kMAXSys][kCENT][kMAXD]
                           [kMAXD]; // Ratio between default - and comparison

TH1D *hSystematic_pbp[Nsets][2][kCENT][kMAXD][kMAXD]; // Systematic point by point
TH1D *hSystematic_sca[Nsets][2][kCENT][kMAXD][kMAXD]; // Systematic scaling
TH1D *hSystematic_pbpFull[2][kCENT][kMAXD][kMAXD];
TH1D *hSystematic_scaFull[2][kCENT][kMAXD][kMAXD];

TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef = 0;
float *ptaborders;
float *pttborders;


//------------------------------------------------------------------------------------------------
void LoadData()
{

  for (int i = 0; i < Nsets; i++)
  {
    for (int isys = 0; isys < Nsysfile[i]; ++isys)
    {
      fin[i][isys] = TFile::Open(infiles[i][isys]);
    }
  }
  fdefault = TFile::Open(defaultfile);

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
  for (int iset = 0; iset < Nsets; iset++)
 {
    for (int isys = 0; isys < Nsysfile[iset]; ++isys)
    {
      for (int idtyp = 0; idtyp < 2;
           idtyp++)
      { // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
        for (int ic = 0; ic < NumCent[idtyp]; ic++)
        {
          for (int iptt = 0; iptt < NPTT; iptt++)
          {
            for (int ipta = 0; ipta < NPTA; ipta++)
            {
              if (pttborders[iptt] - ptaborders[ipta] < 0.000001)
                continue;
              hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta] =
                  (TH1D *)fin[iset][isys]->Get(
                      Form("hDeltaEtaSig%02dC%02dT%02dA%02d", idtyp, ic, iptt,
                           ipta));
              cout << Form("hDeltaEtaSig%02dC%02dT%02dA%02d", idtyp, ic, iptt,
                           ipta)
                   << endl;
              int nbinsx =
                  hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetNbinsX();
              if (isys == 0)
              {

                hSystematic_pbp[iset][idtyp][ic][iptt][ipta] = new TH1D(
                    Form("hSystematic_pbp%02d%02d%02d%02d%02d", iset, idtyp, ic,
                         iptt, ipta),
                    Form("hSystematic_pbp%02d%02d%02d%02d%02d", iset, idtyp, ic,
                         iptt, ipta),
                    nbinsx,
                    hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                        1),
                    hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                        nbinsx + 1));
                hSystematic_sca[iset][idtyp][ic][iptt][ipta] = new TH1D(
                    Form("hSystematic_sca%02d%02d%02d%02d%02d", iset, idtyp, ic,
                         iptt, ipta),
                    Form("hSystematic_sca%02d%02d%02d%02d%02d", iset, idtyp, ic,
                         iptt, ipta),
                    nbinsx,
                    hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                        1),
                    hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                        nbinsx + 1));
                if (iset == 0)
                {
                  hSystematic_scaFull[idtyp][ic][iptt][ipta] = new TH1D(Form("hSystematic_scaFull%02d%02d%02d%02d", idtyp, ic,
                                                                             iptt, ipta),
                                                                        Form("hSystematic_scaFull%02d%02d%02d%02d", idtyp, ic,
                                                                             iptt, ipta),
                                                                        nbinsx,
                                                                        hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                                                                            1),
                                                                        hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                                                                            nbinsx + 1));
                  hSystematic_pbpFull[idtyp][ic][iptt][ipta] = new TH1D(Form("hSystematic_pbpFull%02d%02d%02d%02d", idtyp, ic,
                                                                             iptt, ipta),
                                                                        Form("hSystematic_pbpFull%02d%02d%02d%02d", idtyp, ic,
                                                                             iptt, ipta),
                                                                        nbinsx,
                                                                        hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                                                                            1),
                                                                        hDeltaEtaSig[iset][isys][idtyp][ic][iptt][ipta]->GetBinLowEdge(
                                                                            nbinsx + 1));
                } // iset = 0
              } // isys == 0

          if (idtyp == AA)
            hIAADeltaEtaSig[iset][isys][ic][iptt][ipta] =
                (TH1D *)fin[iset][isys]->Get(
                    Form("hIAADeltaEtaSigC%02dT%02dA%02d", ic, iptt, ipta));
          if (iset == 0 && isys == 0)
          {
            hDeltaEtaSig_def[idtyp][ic][iptt][ipta] = (TH1D *)fdefault->Get(
                Form("hDeltaEtaSig%02dC%02dT%02dA%02d", idtyp, ic, iptt,
                     ipta));
            if (idtyp == AA)
              hIAADeltaEtaSig_def[ic][iptt][ipta] = (TH1D *)fdefault->Get(
                  Form("hIAADeltaEtaSigC%02dT%02dA%02d", ic, iptt, ipta));
              } // Default file if()
            }   // ipta
          }     // iptt
        }       // ic
      }         // pp or AA
    }           // isys
  }             // iset
                // Calculation Various Ratios
                // 1. Out/In ratio for each data set
  for (int iSet = 0; iSet < Nsets; iSet++)
  {
    for (int isys = 0; isys < Nsysfile[iSet]; ++isys)
    {
      for (int idtyp = 0; idtyp < 2; idtyp++)
      {
        for (int ic = 0; ic < NumCent[idtyp]; ic++)
        {
          for (int iptt = 0; iptt < NPTT; iptt++)
          {
            for (int ipta = 0; ipta < NPTA; ipta++)
            {
                            if (pttborders[iptt] - ptaborders[ipta] < 0.000001) continue;

              hRatio_DeltaEtaSig[iSet][isys][idtyp][ic][iptt][ipta] =
                  (TH1D *)hDeltaEtaSig[iSet][isys][idtyp][ic][iptt][ipta]->Clone();
              hRatio_DeltaEtaSig[iSet][isys][idtyp][ic][iptt][ipta]->Divide(
                  hDeltaEtaSig_def[idtyp][ic][iptt][ipta]);
            } // iset
          }   // ipta
        }     // iptt
      }       // ic
    }
  } // isys
} // nsets

//------------------------------------------------------------------------------------------------
void Compare()
{
  LoadData();
  CalculateRatios(1, 3, 3); // no meaning on numbers;;

  int ic = 0;
  for (int iptt = 3; iptt < NPTT; iptt++)
  {
    for (int ipta = 4; ipta < 5; ipta++)
    {
       //DrawSignal(ic++, iptt, ipta, 0);
       for (int iset = 5; iset < 6; iset++){
       DrawIAA(ic, iptt, ipta, iset);

       }
    }
  }
}

//------------------------------------------------------------------------------------------------
void DrawSignal(int padID, int iPTT, int iPTA, int iSet)
{
  Filipad *fpad[NumCent[AA]];
  lowx = -0.01;
  for (int ic = 0; ic < NumCent[AA]; ic++)
  {
    fpad[ic] =
        new Filipad(padID + ic + 1, 1.1, 0.5, 100, 100, 0.7, NumCent[AA]);
    fpad[ic]->Draw();
    //==== Upper pad
    TPad *p = fpad[ic]->GetPad(1); // upper pad
    p->SetTickx();
    p->SetLogx(0);
    p->SetLogy(0);
    p->cd();
    hy = hDeltaEtaSig_def[AA][ic][iPTT][iPTA]->GetMaximum() * 1.5;
    TH2F *hfr = new TH2F("hfr", " ", 100, lowx, highx, 10, ly,
                         hy); // numbers: tics x, low limit x, upper limit x,
                              // tics y, low limit y, upper limit y
    hset(*hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|", 1.1, 1.0, 0.09,
         0.09, 0.01, 0.01, 0.04, 0.05, 510,
         505); // settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    // Legend definition
    TLegend *leg = new TLegend(0.35, 0.4, 0.85, 0.78, "", "brNDC");
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0); // legend settings;

    latexRun.DrawLatexNDC(0.25, 0.85, strRun);

    hDeltaEtaSig_def[AA][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[4]);
    hDeltaEtaSig_def[AA][ic][iPTT][iPTA]->SetMarkerColor(gColors[4]);
    hDeltaEtaSig_def[AA][ic][iPTT][iPTA]->SetLineColor(gColors[4]);
    hDeltaEtaSig_def[AA][ic][iPTT][iPTA]->Draw("p,same");
    leg->AddEntry(hDeltaEtaSig_def[AA][ic][iPTT][iPTA], sLeg[iSet][0], "pl");

    leg->Draw();

    //==== Lower pad
    p = fpad[ic]->GetPad(2);
    p->SetTickx();
    p->SetGridy(1);
    p->SetLogx(0), p->SetLogy(0);
    p->cd();
    TH2F *hfr1 = new TH2F("hfr1", " ", 100, lowx, highx, 10, lowIAA, highIAA);
    hset(*hfr1, "|#Delta#eta|", Form("Ratio to %s", sLeg[iSet][iRef].Data()), 1.1,
         1.0, 0.09, 0.09, 0.01, 0.01, 0.08, 0.08, 510, 505);
    hfr1->Draw();
    for (int isys = 0; isys < Nsysfile[iSet]; isys++)
    {
      hRatio_DeltaEtaSig[iSet][isys][AA][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[isys]);
      hRatio_DeltaEtaSig[iSet][isys][AA][ic][iPTT][iPTA]->SetMarkerColor(gColors[isys]);
      hRatio_DeltaEtaSig[iSet][isys][AA][ic][iPTT][iPTA]->SetLineColor(gColors[isys]);
      hRatio_DeltaEtaSig[iSet][isys][AA][ic][iPTT][iPTA]->Draw("p,same");
    }
    gPad->GetCanvas()->SaveAs(Form(
        "figs_iaa/DeltaEtaSig_AA_Sys%02d_C%02dT%02dA%02d.pdf", iSet, ic, iPTT, iPTA));
  }
  // for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];


  Filipad *fpad2[NumCent[pp]];
  lowx = -0.01;
  for (int ic = 0; ic < NumCent[pp]; ic++)
  {
    fpad2[ic] =
        new Filipad(padID + ic + 10, 1.1, 0.5, 100, 100, 0.7, NumCent[pp]);
    fpad2[ic]->Draw();
    //==== Upper pad
    TPad *p = fpad2[ic]->GetPad(1); // upper pad
    p->SetTickx();
    p->SetLogx(0);
    p->SetLogy(0);
    p->cd();
    hy = hDeltaEtaSig_def[pp][ic][iPTT][iPTA]->GetMaximum() * 1.5;
    TH2F *hfr = new TH2F("hfr", " ", 100, lowx, highx, 10, ly,
                         hy); // numbers: tics x, low limit x, upper limit x,
                              // tics y, low limit y, upper limit y
    hset(*hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|", 1.1, 1.0, 0.09,
         0.09, 0.01, 0.01, 0.04, 0.05, 510,
         505); // settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    // Legend definition
    TLegend *leg = new TLegend(0.35, 0.4, 0.85, 0.78, "", "brNDC");
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0); // legend settings;

    latexRun.DrawLatexNDC(0.25, 0.85, strRun);

    hDeltaEtaSig_def[pp][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iPTA]);
    hDeltaEtaSig_def[pp][ic][iPTT][iPTA]->SetMarkerColor(gColors[iPTA]);
    hDeltaEtaSig_def[pp][ic][iPTT][iPTA]->SetLineColor(gColors[iPTA]);
    hDeltaEtaSig_def[pp][ic][iPTT][iPTA]->Draw("p,same");
    leg->AddEntry(hDeltaEtaSig_def[pp][ic][iPTT][iPTA], sLeg[iSet][0], "pl");

    leg->Draw();

    //==== Lower pad
    p = fpad2[ic]->GetPad(2);
    p->SetTickx();
    p->SetGridy(1);
    p->SetLogx(0), p->SetLogy(0);
    p->cd();
    TH2F *hfr1 = new TH2F("hfr1", " ", 100, lowx, highx, 10, lowIAA, highIAA);
    hset(*hfr1, "|#Delta#eta|", Form("Ratio to %s", sLeg[iSet][iRef].Data()), 1.1,
         1.0, 0.09, 0.09, 0.01, 0.01, 0.08, 0.08, 510, 505);
    hfr1->Draw();
    for (int isys = 0; isys < Nsysfile[iSet]; isys++)
    {
      hRatio_DeltaEtaSig[iSet][isys][pp][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[isys]);
      hRatio_DeltaEtaSig[iSet][isys][pp][ic][iPTT][iPTA]->SetMarkerColor(gColors[isys]);
      hRatio_DeltaEtaSig[iSet][isys][pp][ic][iPTT][iPTA]->SetLineColor(gColors[isys]);
      hRatio_DeltaEtaSig[iSet][isys][pp][ic][iPTT][iPTA]->Draw("p,same");
    }
    gPad->GetCanvas()->SaveAs(Form(
        "figs_iaa/DeltaEtaSig_pp_Sys%02d_C%02dT%02dA%02d.pdf", iSet, ic, iPTT, iPTA));
  }
}

//------------------------------------------------------------------------------------------------
void DrawIAA(int padID, int iPTT, int iPTA, int iSet)
{
  Filipad *fpad[NumCent[AA]];
  lowx = -0.01;
  for (int ic = 0; ic < NumCent[AA]; ic++)
  {
    fpad[ic] =
        new Filipad(padID + ic + 1, 1.1, 0.5, 100, 100, 0.7, NumCent[AA]);
    fpad[ic]->Draw();
    //==== Upper pad
    TPad *p = fpad[ic]->GetPad(1); // upper pad
    p->SetTickx();
    p->SetLogx(0);
    p->SetLogy(0);
    p->cd();
    //    hy = hIAADeltaEtaSig[2][ic][iPTT][iPTA]->GetMaximum() * 1.2;
    TH2F *hfr = new TH2F("hfr", " ", 100, lowx, highx, 10, lowIAA,
                         highIAA); // numbers: tics x, low limit x, upper limit
                                   // x, tics y, low limit y, upper limit y
    hset(*hfr, "|#Delta#eta|", "I_{AA}", 1.1, 1.0, 0.09, 0.09, 0.01, 0.01, 0.04,
         0.05, 510, 505); // settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    // Legend definition
    TLegend *leg = new TLegend(0.45, 0.4, 0.85, 0.78, "", "brNDC");
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0); // legend settings;

    latexRun.DrawLatexNDC(0.25, 0.85, strRun);

    leg->AddEntry(hIAADeltaEtaSig_def[ic][iPTT][iPTA],
                  "Default", "pl");

    hIAADeltaEtaSig_def[ic][iPTT][iPTA]->Draw("p,same");
    for (int isys = 0; isys < Nsysfile[iSet]; isys++)
    {
            hIAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[isys]);
            hIAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->SetMarkerColor(gColors[isys]);
            hIAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->SetLineColor(gColors[isys]);
            hIAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->Draw("p,same");
            leg->AddEntry(hIAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA], sLeg[iSet][isys+1], "pl");
    }

    leg->Draw();

    //==== Lower pad
    p = fpad[ic]->GetPad(2);
    p->SetTickx();
    p->SetGridy(1);
    p->SetLogx(0), p->SetLogy(0);
    p->cd();
    TH2F *hfr1 = new TH2F("hfr1", " ", 100, lowx, highx, 10, lowIAA, highIAA);
    hset(*hfr1, "|#Delta#eta|", Form("Ratio to %s", sLeg[iSet][iRef].Data()), 1.1,
         1.0, 0.09, 0.09, 0.01, 0.01, 0.08, 0.08, 510, 505);
    hfr1->Draw();
     for(int isys=0;isys<Nsysfile[iSet];isys++) {
    	hRatio_IAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[isys]);
    	hRatio_IAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->SetMarkerColor(gColors[isys]);
    	hRatio_IAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->SetLineColor(gColors[isys]);
    	hRatio_IAADeltaEtaSig[iSet][isys][ic][iPTT][iPTA]->Draw("p,same");
    }
    gPad->GetCanvas()->SaveAs(
        Form("figs_iaa/IAA_Set%02d_C%02dT%02dA%02d.pdf", iSet, ic, iPTT, iPTA));
  }
  // for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];
}

/*
//------------------------------------------------------------------------------------------------
void DrawPP(int padID, int iPTT, int iPTA)
{
  Filipad *fpad;
  lowx = -0.01;
  fpad = new Filipad(padID + 1, 1.1, 0.5, 100, 100, 0.7, 5);
  fpad->Draw();
  //==== Upper pad
  TPad *p = fpad->GetPad(1); // upper pad
  p->SetTickx();
  p->SetLogx(0);
  p->SetLogy(0);
  p->cd();
  // hy = hDeltaEtaSig[iRef][pp][0][iPTT][iPTA]->GetMaximum() * 1.2;
  TH2F *hfr = new TH2F("hfr", " ", 100, lowx, highx, 10, ly,
                       hy); // numbers: tics x, low limit x, upper limit x, tics
                            // y, low limit y, upper limit y
  hset(*hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|", 1.1, 1.0, 0.09,
       0.09, 0.01, 0.01, 0.04, 0.05, 510,
       505); // settings of the upper pad: x-axis, y-axis
  hfr->Draw();
  // Legend definition
  TLegend *leg = new TLegend(0.45, 0.4, 0.85, 0.78, "", "brNDC");
  leg->SetTextSize(0.037);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0); // legend settings;

  latexRun.DrawLatexNDC(0.25, 0.85, strRun);

  //  leg->AddEntry((TObject *)NULL, hDeltaEtaSig[0][pp][0][iPTT][iPTA]->GetTitle(),
  //                "");
  //
  //  hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[0]);
  //  hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[0]);
  //  hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetLineColor(gColors[0]);
  //  hDeltaEtaSig[0][pp][0][iPTT][iPTA]->Draw("p,same");
  //  leg->AddEntry(hDeltaEtaSig[0][pp][0][iPTT][iPTA], "Data pp", "pl");
  //  hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[1]);
  //  hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[1]);
  //  hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetLineColor(gColors[1]);
  //  hDeltaEtaSig[1][pp][0][iPTT][iPTA]->Draw("p,same");
  //  leg->AddEntry(hDeltaEtaSig[1][pp][0][iPTT][iPTA], "pythia pp", "pl");
  //
  //  leg->Draw();
  //
  //  //==== Lower pad
  //  p = fpad->GetPad(2);
  //  p->SetTickx();
  //  p->SetGridy(1);
  //  p->SetLogx(0), p->SetLogy(0);
  //  p->cd();
  //  TH2F *hfr1 = new TH2F("hfr1", " ", 100, lowx, highx, 10, lowIAA, highIAA);
  //  hset(*hfr1, "|#Delta#eta|", "Ratio", 1.1, 1.0, 0.09, 0.09, 0.01, 0.01, 0.08,
  //       0.08, 510, 505);
  //  hfr1->Draw();
  //
  //  TH1D *hratio = (TH1D *)hDeltaEtaSig[1][pp][0][iPTT][iPTA]->Clone();
  //  hratio->Divide(hDeltaEtaSig[0][pp][0][iPTT][iPTA]);
  //  hratio->Draw("p,same");
  // gPad->GetCanvas()->SaveAs(Form("figs/DeltaEta_OUTOIN_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
}
*/
//------------------------------------------------------------------------------------------------
void CalculateRatios(int padID, int iPTT, int iPTA)
{
  Filipad *fpad;
  lowx = -0.01;
  fpad = new Filipad(padID + 1, 1.1, 0.5, 100, 100, 0.7, 5);
  fpad->Draw();
  //==== Upper pad
  TPad *p = fpad->GetPad(1); // upper pad
  p->SetTickx();
  p->SetLogx(0);
  p->SetLogy(0);
  p->cd();
  hy = hDeltaEtaSig_def[pp][0][iPTT][iPTA]->GetMaximum() * 1.2;
  TH2F *hfr = new TH2F("hfr", " ", 100, lowx, highx, 10, ly,
                       hy); // numbers: tics x, low limit x, upper limit x, tics
                            // y, low limit y, upper limit y
  hset(*hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|", 1.1, 1.0, 0.09,
       0.09, 0.01, 0.01, 0.04, 0.05, 510,
       505); // settings of the upper pad: x-axis, y-axis
  hfr->Draw();
  // Legend definition
  TLegend *leg = new TLegend(0.45, 0.4, 0.85, 0.78, "", "brNDC");
  leg->SetTextSize(0.037);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0); // legend settings;

  latexRun.DrawLatexNDC(0.25, 0.85, strRun);

  leg->AddEntry((TObject *)NULL, hDeltaEtaSig_def[pp][0][iPTT][iPTA]->GetTitle(),
                "");

  for (int iset = 0; iset < Nsets; ++iset)
  {
    for (int isys = 0; isys < Nsysfile[iset]; ++isys)
    {
      for (int iptt = 0; iptt < NPTT; ++iptt)
      {
        for (int ipta = 0; ipta < NPTA; ++ipta)
        {
                       if (pttborders[iptt] - ptaborders[ipta] < 0.000001) continue;
          hRatio_DeltaEtaSig[iset][isys][pp][0][iptt][ipta] = (TH1D*)
              hDeltaEtaSig[iset][isys][pp][0][iptt][ipta]->Clone(
                  Form("hRatio_DeltaEtaSig%02d%02d%02d%02d%02d%02d", iset,
                       isys, pp, 0, iptt, ipta));
          hRatio_DeltaEtaSig[iset][isys][pp][0][iptt][ipta]->Divide(
              hDeltaEtaSig_def[pp][0][iptt][ipta]);

          hRatio_DeltaEtaSig[iset][isys][pp][0][iptt][ipta]->Fit("pol0", "", "",
                                                                 0, 0.3);

          hRatio_DeltaEtaSig[iset][isys][pp][0][iptt][ipta]->SetMarkerStyle(
              gMarkers[isys]);
          hRatio_DeltaEtaSig[iset][isys][pp][0][iptt][ipta]->SetMarkerColor(
              gColors[isys]);
          hRatio_DeltaEtaSig[iset][isys][pp][0][iptt][ipta]->SetLineColor(
              gColors[isys]);

        } // ipta
      }   // iptt

    } // isys
  }   // iset

  for (int iset = 0; iset < Nsets; ++iset)
  {
    for (int isys = 0; isys < Nsysfile[iset]; ++isys)
    {
      for (int icent = 0; icent < NumCent[AA]; ++icent)
      {
        for (int iptt = 0; iptt < NPTT; ++iptt)
        {
          for (int ipta = 0; ipta < NPTA; ++ipta)
          {
                          if (pttborders[iptt] - ptaborders[ipta] < 0.000001) continue;

            hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta] = (TH1D*)
                hIAADeltaEtaSig[iset][isys][icent][iptt][ipta]->Clone(
                    Form("hRatio_IAADeltaEtaSig%02d%02d%02d%02d%02d", iset,
                         isys, icent, iptt, ipta));
            hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]->Divide(
                hIAADeltaEtaSig_def[icent][iptt][ipta]);

            hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]->Fit(
                "pol0", "", "", 0, 0.3);
            hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]
                ->SetMarkerStyle(gMarkers[isys]);
            hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]
                ->SetMarkerColor(gColors[isys]);
            hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]->SetLineColor(
                gColors[isys]);
          }
        }
      }
    }
  }

  for (int iset = 0; iset < Nsets; ++iset)
  {
    for (int idtyp = 0; idtyp < 2; ++idtyp)
    {
      for (int icent = 0; icent < NumCent[idtyp]; ++icent)
      {
        for (int iptt = 0; iptt < NPTT; ++iptt)
        {
          for (int ipta = 0; ipta < NPTA; ++ipta)
          {
              if (pttborders[iptt] - ptaborders[ipta] < 0.000001) continue;

            for (int isys = 0; isys < Nsysfile[iset]; ++isys)
            {
              if (isp2p[iset])
              {
                int nbins = hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]
                                ->GetNbinsX();
                for (int ibin = 1; ibin <= nbins; ++ibin)
                {
                  double delta =
                      hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]
                          ->GetBinContent(ibin);
                  if (Nsysfile[iset] == 2)
                    delta *= delta / 2;
                  else
                    delta *= delta;
                  delta +=
                      hSystematic_pbp[iset][idtyp][icent][iptt][ipta]->GetBinContent(ibin) *
                      hSystematic_pbp[iset][idtyp][icent][iptt][ipta]->GetBinContent(ibin);

                  double x =
                      hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]
                          ->GetBinCenter(ibin);
                  hSystematic_pbp[iset][idtyp][icent][iptt][ipta]->SetBinContent(
                      ibin, sqrt(delta));
                }
              }
              else
              {
                int nbins = hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]
                                ->GetNbinsX();
                for (int ibin = 1; ibin <= nbins; ++ibin)
                {
                  double delta =
                      hRatio_IAADeltaEtaSig[iset][isys][icent][iptt][ipta]
                          ->GetFunction("pol0")
                          ->GetParameter(0);
                  if (Nsysfile[iset] == 2)
                    delta *= delta / 2;
                  else

                    delta *= delta;
                  delta += hSystematic_pbp[iset][idtyp][icent][iptt][ipta]
                               ->GetBinContent(ibin) *
                           hSystematic_pbp[iset][idtyp][icent][iptt][ipta]
                               ->GetBinContent(ibin);

                  double x = hRatio_IAADeltaEtaSig[iset][isys][icent]
                                                  [iptt][ipta]
                                                      ->GetBinCenter(ibin);
                  hSystematic_sca[iset][idtyp][icent][iptt][ipta]
                      ->SetBinContent(ibin, sqrt(delta));
                } // for ibin
              }   // is p2p or scaling
            }     // isys
          }       // ipta
        }         // iptt
      }           // cent
    }             // idtyp
  }               // iset

  foutfile = new TFile("systematic.root", "recreate");

  for (int idtyp = 0; idtyp < 2; ++idtyp)
  {
    for (int icent = 0; icent < NumCent[idtyp]; ++icent)
    {
      for (int iptt = 0; iptt < NPTT; ++iptt)
      {
        for (int ipta = 0; ipta < NPTA; ++ipta)
        {
          if (pttborders[iptt] - ptaborders[ipta] < 0.000001)
            continue;
          double scasys = 0;
          for (int iset = 0; iset < Nsets; ++iset)
          {
            scasys += hSystematic_sca[iset][idtyp][icent][iptt][ipta]->GetBinContent(1) * hSystematic_sca[iset][idtyp][icent][iptt][ipta]->GetBinContent(1);
          } // iset ipta

          int nbins = hSystematic_pbp[0][idtyp][icent][iptt][ipta]->GetNbinsX();
          for (int ibin = 1; ibin < nbins; ++ibin)
          {
            hSystematic_scaFull[idtyp][icent][iptt][ipta]->SetBinContent(ibin, sqrt(scasys));
          }

          for (int ibin = 1; ibin < nbins; ++ibin)
          {
            double pbpsys = 0;
            for (int iset = 0; iset < Nsets; ++iset)
            {
              pbpsys += hSystematic_pbp[iset][idtyp][icent][iptt][ipta]->GetBinContent(ibin) * hSystematic_pbp[iset][idtyp][icent][iptt][ipta]->GetBinContent(ibin);
            }
            hSystematic_pbpFull[idtyp][icent][iptt][ipta]->SetBinContent(ibin, sqrt(pbpsys));
          }
        } // iptt
      }   // cent
    }     // idtyp
  }       //

  for (int iset = 0; iset < Nsets; ++iset)
  {
    for (int idtyp = 0; idtyp < 2; ++idtyp)
    {
      for (int icent = 0; icent < NumCent[idtyp]; ++icent)
      {
        for (int iptt = 0; iptt < NPTT; ++iptt)
        {
          for (int ipta = 0; ipta < NPTA; ++ipta)
          {
              if (pttborders[iptt] - ptaborders[ipta] < 0.000001) continue;

            for (int isys = 0; isys < Nsysfile[iset]; ++isys)
            {
              hSystematic_pbp[iset][idtyp][icent][iptt][ipta]->Write();
              hSystematic_sca[iset][idtyp][icent][iptt][ipta]->Write();

            }
            if (iset == 0) {
            hSystematic_pbpFull[idtyp][icent][iptt][ipta]->Write();
            hSystematic_scaFull[idtyp][icent][iptt][ipta]->Write();
            }
          }
        }
      }
    }
  }

  

  //          hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[0]);
  //          hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[0]);
  //          hDeltaEtaSig[0][pp][0][iPTT][iPTA]->SetLineColor(gColors[0]);
  //          hDeltaEtaSig[0][pp][0][iPTT][iPTA]->Draw("p,same");
  //          leg->AddEntry(hDeltaEtaSig[0][pp][0][iPTT][iPTA], "Data pp",
  //          "pl");
  //          hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetMarkerStyle(gMarkers[1]);
  //          hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetMarkerColor(gColors[1]);
  //          hDeltaEtaSig[1][pp][0][iPTT][iPTA]->SetLineColor(gColors[1]);
  //          hDeltaEtaSig[1][pp][0][iPTT][iPTA]->Draw("p,same");
  //          leg->AddEntry(hDeltaEtaSig[1][pp][0][iPTT][iPTA], "pythia pp",
  //          "pl");
  //
  //          leg->Draw();
  //
  //          //==== Lower pad
  //          p = fpad->GetPad(2);
  //          p->SetTickx();
  //          p->SetGridy(1);
  //          p->SetLogx(0), p->SetLogy(0);
  //          p->cd();
  //          TH2F *hfr1 =
  //              new TH2F("hfr1", " ", 100, lowx, highx, 10, lowIAA, highIAA);
  //          hset(*hfr1, "|#Delta#eta|", "Ratio", 1.1, 1.0, 0.09, 0.09, 0.01,
  //          0.01,
  //               0.08, 0.08, 510, 505);
  //          hfr1->Draw();
  //
  //          TH1D *hratio = (TH1D
  //          *)hDeltaEtaSig[1][pp][0][iPTT][iPTA]->Clone();
  //          hratio->Divide(hDeltaEtaSig[0][pp][0][iPTT][iPTA]);
  //          hratio->Draw("p,same");
  // gPad->GetCanvas()->SaveAs(Form("figs/DeltaEta_OUTOIN_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
}

void DrawSystematics(int padID, int iPTT, int iPTA) {}

//------------------------------------------------------------------------------------------------

/*
void CalculateSysPP(int padID, int iPTT, int iPTA)
{
  Filipad *fpad[NumCent[AA]];
  lowx = -0.01;
  for (int ic = 0; ic < NumCent[AA]; ic++)
  {
    fpad[ic] =
        new Filipad(padID + ic + 1, 1.1, 0.5, 100, 100, 0.7, NumCent[AA]);
    fpad[ic]->Draw();
    //==== Upper pad
    TPad *p = fpad[ic]->GetPad(1); // upper pad
    p->SetTickx();
    p->SetLogx(0);
    p->SetLogy(0);
    p->cd();
    hy = hDeltaEtaSig_def[pp][0][iPTT][iPTA]->GetMaximum() * 1.2;
    TH2F *hfr = new TH2F("hfr", " ", 100, lowx, highx, 10, lowIAA,
                         highIAA); // numbers: tics x, low limit x, upper limit
                                   // x, tics y, low limit y, upper limit y
    hset(*hfr, "|#Delta#eta|", "I_{AA}", 1.1, 1.0, 0.09, 0.09, 0.01, 0.01, 0.04,
         0.05, 510,
         505); // settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    // Legend definition
    TLegend *leg = new TLegend(0.45, 0.4, 0.85, 0.78, "", "brNDC");
    leg->SetTextSize(0.037);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0); // legend settings;

    latexRun.DrawLatexNDC(0.25, 0.85, strRun);

    leg->AddEntry((TObject *)NULL,
                  hIAADeltaEtaSig[0][ic][iPTT][iPTA]->GetTitle(), "");

    for (int iS = 0; iS < Nsets; iS++)
    {
      hDeltaEtaSig_def[iS][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
      hDeltaEtaSig_def[iS][ic][iPTT][iPTA]->SetMarkerColor(gColors[iS]);
      hDeltaEtaSig_def[iS][ic][iPTT][iPTA]->SetLineColor(gColors[iS]);
      hDeltaEtaSig_def[iS][ic][iPTT][iPTA]->Draw("p,same");
      leg->AddEntry(hDeltaEtaSig_def[iS][ic][iPTT][iPTA], sLeg[iS], "pl");
    }

    leg->Draw();

    //==== Lower pad
    p = fpad[ic]->GetPad(2);
    p->SetTickx();
    p->SetGridy(1);
    p->SetLogx(0), p->SetLogy(0);
    p->cd();
    TH2F *hfr1 = new TH2F("hfr1", " ", 100, lowx, highx, 10, lowIAA, highIAA);
    hset(*hfr1, "|#Delta#eta|", Form("Ratio to %s", sLeg[iRef].Data()), 1.1,
         1.0, 0.09, 0.09, 0.01, 0.01, 0.08, 0.08, 510, 505);
    hfr1->Draw();
    // for(int i=0;i<Nsets;i++) {
    //	hRatio_IAA[i][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
    //	hRatio_IAA[i][ic][iPTT][iPTA]->SetMarkerColor(gColors[i]);
    //	hRatio_IAA[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
    //	hRatio_IAA[i][ic][iPTT][iPTA]->Draw("p,same");
    //}
    gPad->GetCanvas()->SaveAs(
        Form("figs_iaa/IAA_BeamEnergy_C%02dT%02dA%02d.pdf", ic, iPTT, iPTA));
  }
  // for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];
}
*/
