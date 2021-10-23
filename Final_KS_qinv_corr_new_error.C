#include <iostream> 
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector.h"
#include "Final_KS_qinv_corr_new_error.h"
#include <cstdlib>
#include <string>
#include "TString.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include <chrono>

using namespace std;
using namespace TMVA;
using namespace std::chrono;
  
void AssignpTbins(double V0pt, double V0eta, double V0phi,double V0id, double V0mass, double pdau_pt, double pdau_eta, double pdau_phi, double ndau_pt, double ndau_eta, double ndau_phi, int idx, int ctbin);

void AnalyseEvents(std::vector<double> zvtx[NBINS_cent], std::vector<int> evtcount[NBINS_cent], std::vector<float> cntbin[NBINS_cent], std::vector<double> evtno[NBINS_cent], std::vector<double> runno[NBINS_cent], std::vector<double> vect_V0id_trg[NBINS_cent][iiev], std::vector<double> vect_V0id_ass[NBINS_cent][iiev], std:: vector<TLorentzVector> V0_vect_trg[NBINS_cent][iiev], std:: vector<TLorentzVector> daup_vect_trg[NBINS_cent][iiev], std:: vector<TLorentzVector> daun_vect_trg[NBINS_cent][iiev], std:: vector<TLorentzVector> V0_vect_ass[NBINS_cent][iiev], std:: vector<TLorentzVector> daup_vect_ass[NBINS_cent][iiev], std:: vector<TLorentzVector> daun_vect_ass[NBINS_cent][iiev]);

void SignalTwoPartCorr(unsigned int iev, unsigned int ic, float centbin);
void TwoPartCorr(unsigned int iev, unsigned int ic);

void MixedTwoPartCorr(unsigned int iev, unsigned int jev, int ic, float centbin);

double GetQInv(double px1,double py1,double pz1,double E1,double px2,double py2,double pz2,double E2);
double GetInvMass(double e1, double px1, double py1, double pz1, double e2, double px2, double py2, double pz2);
int GetCentbin(float Cent);

int GetCent(float Cent);

int GetPtbin(float pt);
int GetKtbin(float kt);
int Getqinvbin(Float_t qinv);

void Final_KS_qinv_corr_new_error(Int_t kFile = 6339)
{

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  auto start = high_resolution_clock::now();
    
  hV0mass_KS_.resize(NBINS_cent);
  //hV0mass_KS_cbin_qbin_.resize(NBINS_cent);

  
  for(int icent = 0; icent < NBINS_cent; icent++)
    {
      hV0mass_KS_[icent].resize(NBINS_pt);
      //hV0mass_KS_cbin_qbin_[icent].resize(NBINS_qbin);

      hf2_sig_[icent] = new TH1D(Form("hf2_sig_%d",icent),"", 75, 1.2, 1.8);
      hf2_mix_[icent] = new TH1D(Form("hf2_mix_%d",icent),"", 75, 1.2, 1.8);
      
      h_V0kt_[icent] = new TH1D(Form("h_V0kt_%d",icent),"", 95, -0.5, 9.0);
      h_V0mt_[icent] = new TH1D(Form("h_V0mt_%d",icent),"", 95, -0.5, 9.0);
      
      hcent_dist[icent] = new TH1D(Form("hcent_dist%d",icent),"", 200 ,0, 200);
      
      hV0mass_KS_cent_[icent] = new TH1D(Form("hV0mass_KS_cent_%d",icent),"", 260 ,0.43, 0.56);
      
      hV0mass_KS_cent_qinv0p8_[icent] = new TH1D(Form("hV0mass_KS_cent_qinv0p8_%d",icent),"", 260 ,0.43, 0.56);
      hV0mass_KS_cent_qinv1p0_[icent] = new TH1D(Form("hV0mass_KS_cent_qinv1p0_%d",icent),"", 260 ,0.43, 0.56);
      hV0mass_KS_cent_qinv1p5_[icent] = new TH1D(Form("hV0mass_KS_cent_qinv1p5_%d",icent),"", 260 ,0.43, 0.56);
      hV0mass_KS_cent_qinv2p0_[icent] = new TH1D(Form("hV0mass_KS_cent_qinv2p0_%d",icent),"", 260 ,0.43, 0.56);
      
      hsig_KS_[icent] = new TH1D(Form("hsig_KS_%d",icent),"", 500, 0, 10);
      hsig_KS_rotate_[icent] = new TH1D(Form("hsig_KS_rotate_%d",icent),"", 500, 0, 10);
      hsig_KS_invert_[icent] = new TH1D(Form("hsig_KS_invert_%d",icent),"", 500, 0, 10);
      
      hmix_KS_[icent] = new TH1D(Form("hmix_KS_%d",icent),"", 500, 0, 10);
      hmix_KS_rotate_[icent] = new TH1D(Form("hmix_KS_rotate_%d",icent),"", 500, 0, 10);
      hmix_KS_invert_[icent] = new TH1D(Form("hmix_KS_invert_%d",icent),"", 500, 0, 10);

      hsig_KS_AB_[icent] = new TH1D(Form("hsig_KS_AB_%d",icent),"", 500, 0, 10);
      hsig_KS_left_AB_[icent] = new TH1D(Form("hsig_KS_left_AB_%d",icent),"", 500, 0, 10);
      hsig_KS_right_AB_[icent] = new TH1D(Form("hsig_KS_right_AB_%d",icent),"", 500, 0, 10);
      hsig_KS_rotate_AB_[icent] = new TH1D(Form("hsig_KS_rotate_AB_%d",icent),"", 500, 0, 10);
      hsig_KS_invert_AB_[icent] = new TH1D(Form("hsig_KS_invert_AB_%d",icent),"", 500, 0, 10);

      hmix_KS_AB_[icent] = new TH1D(Form("hmix_KS_AB_%d",icent),"", 500, 0, 10);
      hmix_KS_left_AB_[icent] = new TH1D(Form("hmix_KS_left_AB_%d",icent),"", 500, 0, 10);
      hmix_KS_right_AB_[icent] = new TH1D(Form("hmix_KS_right_AB_%d",icent),"", 500, 0, 10);
      hmix_KS_rotate_AB_[icent] = new TH1D(Form("hmix_KS_rotate_AB_%d",icent),"", 500, 0, 10);
      hmix_KS_invert_AB_[icent] = new TH1D(Form("hmix_KS_invert_AB_%d",icent),"", 500, 0, 10);

      hsig_KS_BB_[icent] = new TH1D(Form("hsig_KS_BB_%d",icent),"", 500, 0, 10);
      hsig_KS_BB1_[icent] = new TH1D(Form("hsig_KS_BB1_%d",icent),"", 500, 0, 10);
      hsig_KS_BB2_[icent] = new TH1D(Form("hsig_KS_BB2_%d",icent),"", 500, 0, 10);
      hsig_KS_left_BB_[icent] = new TH1D(Form("hsig_KS_left_BB_%d",icent),"", 500, 0, 10);
      hsig_KS_right_BB_[icent] = new TH1D(Form("hsig_KS_right_BB_%d",icent),"", 500, 0, 10);
      hsig_KS_rotate_BB_[icent] = new TH1D(Form("hsig_KS_rotate_BB_%d",icent),"", 500, 0, 10);
      hsig_KS_invert_BB_[icent] = new TH1D(Form("hsig_KS_invert_BB_%d",icent),"", 500, 0, 10);

      hmix_KS_BB_[icent] = new TH1D(Form("hmix_KS_BB_%d",icent),"", 500, 0, 10);
      hmix_KS_BB1_[icent] = new TH1D(Form("hmix_KS_BB1_%d",icent),"", 500, 0, 10);
      hmix_KS_BB2_[icent] = new TH1D(Form("hmix_KS_BB2_%d",icent),"", 500, 0, 10);
      hmix_KS_left_BB_[icent] = new TH1D(Form("hmix_KS_left_BB_%d",icent),"", 500, 0, 10);
      hmix_KS_right_BB_[icent] = new TH1D(Form("hmix_KS_right_BB_%d",icent),"", 500, 0, 10);
      hmix_KS_rotate_BB_[icent] = new TH1D(Form("hmix_KS_rotate_BB_%d",icent),"", 500, 0, 10);
      hmix_KS_invert_BB_[icent] = new TH1D(Form("hmix_KS_invert_BB_%d",icent),"", 500, 0, 10);
    }
  
  for(int icent = 0; icent < NBINS_cent; ++icent)
    {
      for(int ipt = 0; ipt < NBINS_pt; ++ipt)
        {
          hV0mass_KS_[icent][ipt] = new TH1D(Form("hV0mass_KS_centbin_%d_ptbin_%d",icent,ipt),Form("hV0mass_KS_centbin_%d_ptbin_%d",icent,ipt), 260, 0.43, 0.56);
        }
    }

  /*
  for(int icent = 0; icent < NBINS_cent; ++icent)
    {
      for(int iq = 0; iq < NBINS_qbin; ++iq)
        {
	  hV0mass_KS_cbin_qbin_[icent][iq] = new TH1D(Form("hV0mass_KS_centbin_%d_qbin_%d",icent,iq),Form("hV0mass_KS_centbin_%d_qbin_%d",icent,iq), 260, 0.43, 0.56);
	}
    }
  */

  for(int iq = 0; iq < NBINS_qbin; ++iq)
    {
      hV0mass_KS_cbin_qbin_[iq] = new TH1D(Form("hV0mass_KS_qbin_%d",iq),Form("hV0mass_KS_qbin_%d",iq), 260, 0.43, 0.56);
    }
  
  for(int ikt = 0; ikt < NBINS_kt; ++ikt)
    {
      hsig_KS_kt_[ikt] = new TH1D(Form("hsig_KS_kt_%d",ikt),"", 500, 0, 10);
      hmix_KS_kt_[ikt] = new TH1D(Form("hmix_KS_kt_%d",ikt),"", 500, 0, 10);
    }
  
  hsig_KS_ckt_.resize(NBINS_cent);
  hmix_KS_ckt_.resize(NBINS_cent);
  
   for(int icent = 0; icent < NBINS_cent; ++icent)
     {
       hsig_KS_ckt_[icent].resize(NBINS_kt);
       hmix_KS_ckt_[icent].resize(NBINS_kt);
     }
   
   
   for(int icent = 0; icent < NBINS_cent; ++icent)
     {
       for(int ikt = 0; ikt < NBINS_kt; ++ikt)
	 {
	   hsig_KS_ckt_[icent][ikt] = new TH1D(Form("hsig_KS_ckt_%d%d",icent,ikt),"", 500, 0, 10);
	   hmix_KS_ckt_[icent][ikt] = new TH1D(Form("hmix_KS_ckt_%d%d",icent,ikt),"", 500, 0, 10);	   
	 }
     }
   
   
   for(Int_t z_rsz = 0; z_rsz < NBINS_cent; z_rsz++)
     {
       vzvtx[z_rsz].resize(0);
       cntbin[z_rsz].resize(0);
       evtno[z_rsz].resize(0);
       runno[z_rsz].resize(0);
       evtcount[z_rsz].resize(0);
       
       for(Int_t k = 0; k < iiev; k++)
	 
	 {
	   vect_V0id_trg[z_rsz][k].resize(0);
	   V0_vect_trg[z_rsz][k].resize(0);
	   daup_vect_trg[z_rsz][k].resize(0);
	   daun_vect_trg[z_rsz][k].resize(0);
	   
	   vect_V0id_ass[z_rsz][k].resize(0);
	   V0_vect_ass[z_rsz][k].resize(0);
	   daup_vect_ass[z_rsz][k].resize(0);
	   daun_vect_ass[z_rsz][k].resize(0);
	 }
     }
   
   TString Str;
   ifstream fpr("2018_HM.txt", ios::in);
   
   for(Int_t ifile = 0; ifile < kFile; ifile++)
     {
       fpr >> Str;
       
       if(ifile % 500 ==0)
	 {
	   
	   cout << ifile <<"  files run"<< endl;
	   cout<<"file name: "<<Str<<std::endl;
	 }
       
       TFile *file = new TFile(Str, "READ");
       TDirectory *dir1 = (TDirectory*)file->Get("RAGHUV0");
       TDirectory *dir2 = (TDirectory*)dir1->Get("V0tree");
       
       TTree *evt_tree =(TTree*)dir2->Get("Events");
       Double_t        zvtx;
       Float_t         centbin;
       Int_t           run_no;
       Int_t           event_no;
       
       TBranch        *b_zBestVtx1_;
       TBranch        *b_cent;
       TBranch        *b_run;
       TBranch        *b_event;
       
       evt_tree->SetBranchAddress("zvtx", &zvtx, &b_zBestVtx1_);
       evt_tree->SetBranchAddress("centbin", &centbin, &b_cent);
       evt_tree->SetBranchAddress("run_no", &run_no, &b_run);
       evt_tree->SetBranchAddress("event_no", &event_no, &b_event);
       
       
       TTree *V0_tree =(TTree*)dir2->Get("Kshort");
       
       vector<double>  *V0mass_kshort=0;
       vector<double>  *V0pt_kshort=0;
       vector<double>  *V0eta_kshort=0;
       vector<double>  *V0rapidity_kshort=0;
       vector<double>  *V0phi_kshort=0;
       vector<double>  *V0id_kshort=0;
       vector<double>  *V0Chi2_kshort=0;
       vector<double>  *V0lxyz_kshort=0;
       vector<double>  *V0costhetaXYZ_kshort=0;
       vector<double>  *V0vtxDecaySigXYZ_kshort=0;
       vector<double>  *dca_kshort=0;
       vector<double>  *mva_kshort=0;
       
       
       vector<double>  *pdau_kshort_pt=0;
       vector<double>  *pdau_kshort_pterror=0;
       vector<double>  *pdau_kshort_eta=0;
       vector<double>  *pdau_kshort_phi=0;
       vector<double>  *pdau_kshort_nchi2=0;
       vector<double>  *pdau_kshort_nhits=0;
       vector<double>  *pdau_kshort_NPxlayer=0;
       vector<double>  *pdau_kshort_DCAsigXY=0;
       vector<double>  *pdau_kshort_DCAsigZ=0;
       
       
       vector<double>  *ndau_kshort_pt=0;
       vector<double>  *ndau_kshort_pterror=0;
       vector<double>  *ndau_kshort_eta=0;
       vector<double>  *ndau_kshort_phi=0;
       vector<double>  *ndau_kshort_nchi2=0;
       vector<double>  *ndau_kshort_nhits=0;
       vector<double>  *ndau_kshort_NPxlayer=0;
       vector<double>  *ndau_kshort_DCAsigXY=0;
       vector<double>  *ndau_kshort_DCAsigZ=0;
       
       TBranch        *b_V0mass;
       TBranch        *b_V0pt;
       TBranch        *b_V0eta;
       TBranch        *b_V0rapidity;
       TBranch        *b_V0phi;
       TBranch        *b_V0id;
       TBranch        *b_V0Chi2;
       TBranch        *b_V0lxyz;
       TBranch        *b_V0costhetaXYZ;
       TBranch        *b_V0vtxDecaySigXYZ;
       TBranch        *b_dca;
       TBranch        *b_mva;
       
       
       TBranch        *b_pdau_pt;
       TBranch        *b_pdau_pterror;
       TBranch        *b_pdau_eta;
       TBranch        *b_pdau_phi;
       TBranch        *b_pdau_nchi2;
       TBranch        *b_pdau_nhits;
       TBranch        *b_pdau_NPxlayer;
       TBranch        *b_pdau_DCAsigXY;
       TBranch        *b_pdau_DCAsigZ;
       
       TBranch        *b_ndau_pt;
       TBranch        *b_ndau_pterror;
       TBranch        *b_ndau_eta;
       TBranch        *b_ndau_phi;
       TBranch        *b_ndau_nchi2;
       TBranch        *b_ndau_nhits;
       TBranch        *b_ndau_NPxlayer;
       TBranch        *b_ndau_DCAsigXY;
       TBranch        *b_ndau_DCAsigZ;
       
      
      V0_tree->SetBranchAddress("V0mass_kshort", &V0mass_kshort, &b_V0mass);
      V0_tree->SetBranchAddress("V0pt_kshort", &V0pt_kshort, &b_V0pt);
      V0_tree->SetBranchAddress("V0eta_kshort", &V0eta_kshort, &b_V0eta);
      V0_tree->SetBranchAddress("V0rapidity_kshort", &V0rapidity_kshort, &b_V0rapidity);
      V0_tree->SetBranchAddress("V0phi_kshort", &V0phi_kshort, &b_V0phi);
      V0_tree->SetBranchAddress("V0id_kshort", &V0id_kshort, &b_V0id);
      V0_tree->SetBranchAddress("V0Chi2_kshort", &V0Chi2_kshort, &b_V0Chi2);
      V0_tree->SetBranchAddress("V0lxyz_kshort", &V0lxyz_kshort, &b_V0lxyz);
      V0_tree->SetBranchAddress("V0costhetaXYZ_kshort", &V0costhetaXYZ_kshort, &b_V0costhetaXYZ);
      V0_tree->SetBranchAddress("V0vtxDecaySigXYZ_kshort", &V0vtxDecaySigXYZ_kshort, &b_V0vtxDecaySigXYZ);
      V0_tree->SetBranchAddress("dca_kshort", &dca_kshort, &b_dca);
      V0_tree->SetBranchAddress("mva_kshort", &mva_kshort, &b_mva);


      V0_tree->SetBranchAddress("pdau_kshort_pt", &pdau_kshort_pt, &b_pdau_pt);
      V0_tree->SetBranchAddress("pdau_kshort_pterror", &pdau_kshort_pterror, &b_pdau_pterror);
      V0_tree->SetBranchAddress("pdau_kshort_eta", &pdau_kshort_eta, &b_pdau_eta);
      V0_tree->SetBranchAddress("pdau_kshort_phi", &pdau_kshort_phi, &b_pdau_phi);
      V0_tree->SetBranchAddress("pdau_kshort_nchi2", &pdau_kshort_nchi2, &b_pdau_nchi2);
      V0_tree->SetBranchAddress("pdau_kshort_nhits", &pdau_kshort_nhits, &b_pdau_nhits);
      V0_tree->SetBranchAddress("pdau_kshort_NPxlayer", &pdau_kshort_NPxlayer, &b_pdau_NPxlayer);
      V0_tree->SetBranchAddress("pdau_kshort_DCAsigXY", &pdau_kshort_DCAsigXY, &b_pdau_DCAsigXY);
      V0_tree->SetBranchAddress("pdau_kshort_DCAsigZ", &pdau_kshort_DCAsigZ, &b_pdau_DCAsigZ);

      
      V0_tree->SetBranchAddress("ndau_kshort_pt", &ndau_kshort_pt, &b_ndau_pt);
      V0_tree->SetBranchAddress("ndau_kshort_pterror", &ndau_kshort_pterror, &b_ndau_pterror);
      V0_tree->SetBranchAddress("ndau_kshort_eta", &ndau_kshort_eta, &b_ndau_eta);
      V0_tree->SetBranchAddress("ndau_kshort_phi", &ndau_kshort_phi, &b_ndau_phi);
      V0_tree->SetBranchAddress("ndau_kshort_nchi2", &ndau_kshort_nchi2, &b_ndau_nchi2);
      V0_tree->SetBranchAddress("ndau_kshort_nhits", &ndau_kshort_nhits, &b_ndau_nhits);
      V0_tree->SetBranchAddress("ndau_kshort_NPxlayer", &ndau_kshort_NPxlayer, &b_ndau_NPxlayer);
      V0_tree->SetBranchAddress("ndau_kshort_DCAsigXY", &ndau_kshort_DCAsigXY, &b_ndau_DCAsigXY);
      V0_tree->SetBranchAddress("ndau_kshort_DCAsigZ", &ndau_kshort_DCAsigZ, &b_ndau_DCAsigZ);


      //--------------------------------------------------------------------

      Int_t nevent = evt_tree->GetEntries();
      
      Int_t idx  = 0;
      for(Int_t ievt = 0; ievt < nevent; ievt++)
	{
	  evt_tree->GetEntry(ievt);
	  cbin = centbin;

	  if (cbin >= 120.) continue;
	  
	  hcentbin_dist->Fill(cbin);
	  
	  Long64_t tentry = V0_tree->LoadTree(ievt);
	  if(tentry < 0) break;
	  
	  b_V0mass->GetEntry(tentry);
	  b_V0pt->GetEntry(tentry);
	  b_V0eta->GetEntry(tentry);
	  b_V0rapidity->GetEntry(tentry);
	  b_V0phi->GetEntry(tentry);
	  b_V0id->GetEntry(tentry);
	  b_V0Chi2->GetEntry(tentry);
	  b_V0lxyz->GetEntry(tentry);
	  b_V0costhetaXYZ->GetEntry(tentry);
	  b_V0vtxDecaySigXYZ->GetEntry(tentry);
	  b_dca->GetEntry(tentry);
	  b_mva->GetEntry(tentry);
	  
	  b_pdau_pt->GetEntry(tentry);
	  b_pdau_pterror->GetEntry(tentry);
	  b_pdau_eta->GetEntry(tentry);
	  b_pdau_phi->GetEntry(tentry);
	  b_pdau_nchi2->GetEntry(tentry);
	  b_pdau_nhits->GetEntry(tentry);
	  b_pdau_NPxlayer->GetEntry(tentry);
	  b_pdau_DCAsigXY->GetEntry(tentry);
	  b_pdau_DCAsigZ->GetEntry(tentry);
	  
	  
	  b_ndau_pt->GetEntry(tentry);
	  b_ndau_pterror->GetEntry(tentry);
	  b_ndau_eta->GetEntry(tentry);
	  b_ndau_phi->GetEntry(tentry);
	  b_ndau_nchi2->GetEntry(tentry);
	  b_ndau_nhits->GetEntry(tentry);
	  b_ndau_NPxlayer->GetEntry(tentry);
	  b_ndau_DCAsigXY->GetEntry(tentry);
	  b_ndau_DCAsigZ->GetEntry(tentry);
	  
	  std::vector<double> V0pt_vect;
	  std::vector<double> V0eta_vect;
	  std::vector<double> V0phi_vect;
	  std::vector<int> V0id_vect;
	  std::vector<double> V0mass_vect;
	  
	  std::vector<double> pdau_pt_vect;
	  std::vector<double> pdau_eta_vect;
	  std::vector<double> pdau_phi_vect;
	  std::vector<double> pdau_chi2_vect;
	  
	  std::vector<double> ndau_pt_vect;
	  std::vector<double> ndau_eta_vect;
	  std::vector<double> ndau_phi_vect;
	  std::vector<double> ndau_chi2_vect;

	  Int_t npar = 0; 
	  for (unsigned int i = 0; i< V0pt_kshort->size(); i++)
	    {
	      V0_mass = (*V0mass_kshort)[i];
              V0_pt = (*V0pt_kshort)[i];
              V0_eta = (*V0eta_kshort)[i];
              V0_rapidity = (*V0rapidity_kshort)[i];
              V0_phi = (*V0phi_kshort)[i];
              V0_id = (*V0id_kshort)[i];
              V0_dca = (*dca_kshort)[i];
              V0_mva = (*mva_kshort)[i];
              V0_Chi2 = (*V0Chi2_kshort)[i];
              V0_lxyz = (*V0lxyz_kshort)[i];
              V0_costhetaXYZ = (*V0costhetaXYZ_kshort)[i];
              V0_vtxDecaySigXYZ = (*V0vtxDecaySigXYZ_kshort)[i];

	      daup_pt = (*pdau_kshort_pt)[i];
              daup_eta = (*pdau_kshort_eta)[i];
              daup_phi = (*pdau_kshort_phi)[i];
              daup_chi2 = (*pdau_kshort_nchi2)[i];
              daup_pterror = (*pdau_kshort_pterror)[i];
              daup_nhits = (*pdau_kshort_nhits)[i];
              daup_NPxlayer = (*pdau_kshort_NPxlayer)[i];
              daup_DCAsigXY = (*pdau_kshort_DCAsigXY)[i];
              daup_DCAsigZ = (*pdau_kshort_DCAsigZ)[i];


              daun_pt = (*ndau_kshort_pt)[i];
              daun_eta = (*ndau_kshort_eta)[i];
              daun_phi = (*ndau_kshort_phi)[i];
              daun_chi2 = (*ndau_kshort_nchi2)[i];
              daun_pterror = (*ndau_kshort_pterror)[i];
              daun_nhits = (*ndau_kshort_nhits)[i];
              daun_NPxlayer = (*ndau_kshort_NPxlayer)[i];
              daun_DCAsigXY = (*ndau_kshort_DCAsigXY)[i];
              daun_DCAsigZ = (*ndau_kshort_DCAsigZ)[i];

	      
	      if ( V0_mva <= bdt_cut) continue;
	      if ( TMath::Abs(V0_id) != 310) continue;
	      if ( TMath::Abs(V0_rapidity) > Rapmax) continue; //|| V0_rapidity <= Rapmin ) continue;
	      if ( V0_mass > Massmax || V0_mass < Massmin ) continue;
	      if ( V0_vtxDecaySigXYZ <= decaylength ) continue;
	      if ( V0_costhetaXYZ <= costheta) continue;
	      
	      TVector3 pV0vector(V0_pt * TMath::Cos(V0_phi), V0_pt * TMath::Sin(V0_phi), V0_pt * TMath::SinH(V0_eta));
	      TVector3 pd1vector(daup_pt * TMath::Cos(daup_phi), daup_pt * TMath::Sin(daup_phi), daup_pt * TMath::SinH(daup_eta));
	      TVector3 pd2vector(daun_pt * TMath::Cos(daun_phi), daun_pt * TMath::Sin(daun_phi), daun_pt * TMath::SinH(daun_eta));

	      TVector3 dauvec1(daup_pt * TMath::Cos(daup_phi), daup_pt * TMath::Sin(daup_phi), daup_pt * TMath::SinH(daup_eta));
	      TVector3 dauvec2(daun_pt * TMath::Cos(daun_phi), daun_pt * TMath::Sin(daun_phi), daun_pt * TMath::SinH(daun_eta));
	      TVector3 dauvecsum(dauvec1+dauvec2);
	      
	      double pd1 = TMath::Sqrt(pow(daup_pt * TMath::Cos(daup_phi),2)+ pow(daup_pt * TMath::Sin(daup_phi),2)+ pow(daup_pt * TMath::SinH(daup_eta),2));
	      double pd2 = TMath::Sqrt(pow(daun_pt * TMath::Cos(daun_phi),2)+ pow(daun_pt * TMath::Sin(daun_phi),2)+ pow(daun_pt * TMath::SinH(daun_eta),2));
	      	      
	      double v0masspiproton1 = sqrt((sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))*(sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))-dauvecsum.Mag2());
	      double v0masspiproton2 = sqrt((sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))*(sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))-dauvecsum.Mag2());

	     hpip1->Fill(v0masspiproton1);
	     hpip2->Fill(v0masspiproton2);

	     if((v0masspiproton1>=(1.115683-mis_la_range) && v0masspiproton1<=(1.115683+mis_la_range)) || (v0masspiproton2>=(1.115683-mis_la_range) && v0masspiproton2<=(1.115683+mis_la_range)) ) continue;
	     

	     int ptbin = GetPtbin(V0_pt);
             int ctbin_pt = GetCentbin(cbin);
	     
             hV0mass_KS_[ctbin_pt][ptbin]->Fill(V0_mass);
	     
             if ( V0_pt > pTmax || V0_pt < pTmin ) continue;
	     
	     
	      double P1l = (pV0vector.Dot(pd1vector))/pV0vector.Mag();
	      double P2l = (pV0vector.Dot(pd2vector))/pV0vector.Mag();
	      
	      double qt = ((pd1vector.Cross(pd2vector)).Mag())/pV0vector.Mag();
	      double alpha = (P1l - P2l)/(P1l+P2l);
	      
	      V0pt_vect.push_back(V0_pt);
	      V0eta_vect.push_back(V0_eta);
	      V0phi_vect.push_back(V0_phi);
	      V0id_vect.push_back(V0_id);
	      V0mass_vect.push_back(V0_mass);
	      
	      pdau_pt_vect.push_back(daup_pt);
	      pdau_eta_vect.push_back(daup_eta);
	      pdau_phi_vect.push_back(daup_phi);
	      pdau_chi2_vect.push_back(daup_chi2);
	      
	      ndau_pt_vect.push_back(daun_pt);
	      ndau_eta_vect.push_back(daun_eta);
	      ndau_phi_vect.push_back(daun_phi);
	      ndau_chi2_vect.push_back(daun_chi2);
	      
	      
	      h_V0mass->Fill(V0_mass);
	      h_V0pt->Fill(V0_pt);
	      h_V0eta->Fill(V0_eta);
	      h_V0phi->Fill(V0_phi);
	      h_V0rapidity->Fill(V0_rapidity);
	      if (V0_mass >= mass_sig_min && V0_mass <= mass_sig_max)
		{
		  harmen->Fill(alpha, qt);
		}
	      npar++;
	    }//particle loop end-------------------------------------
	  
	  if (npar < 2)
            {
	      V0pt_vect.clear();
              V0eta_vect.clear();
              V0phi_vect.clear();
              V0id_vect.clear();
              V0mass_vect.clear();

              pdau_pt_vect.clear();
              pdau_eta_vect.clear();
              pdau_phi_vect.clear();
              pdau_chi2_vect.clear();

              ndau_pt_vect.clear();
              ndau_eta_vect.clear();
              ndau_phi_vect.clear();
              ndau_chi2_vect.clear();

              continue;
	    }
	  

	  int ctbin = GetCentbin(cbin);
	  
	  vzvtx[ctbin].push_back(zvtx);
	  evtno[ctbin].push_back(event_no);
	  runno[ctbin].push_back(run_no);
	  cntbin[ctbin].push_back(cbin);

	  for (unsigned int ntrk_1 = 0; ntrk_1 < V0pt_vect.size(); ntrk_1 ++)
	    {
	      bool is_dcut = kFALSE;
	      
	      double V0pt_pure = V0pt_vect[ntrk_1];
	      double V0eta_pure = V0eta_vect[ntrk_1];
	      double V0phi_pure = V0phi_vect[ntrk_1];
	      double V0id_pure = V0id_vect[ntrk_1];
	      double V0mass_pure = V0mass_vect[ntrk_1];
	      
	      double ppt_1 = pdau_pt_vect[ntrk_1];	  
	      double peta_1 = pdau_eta_vect[ntrk_1];
	      double pphi_1 = pdau_phi_vect[ntrk_1];
	      double pchi2_1 = pdau_chi2_vect[ntrk_1];
	      
	      double npt_1 = ndau_pt_vect[ntrk_1];
	      double neta_1 = ndau_eta_vect[ntrk_1];
	      double nphi_1 = ndau_phi_vect[ntrk_1];
	      double nchi2_1 = ndau_chi2_vect[ntrk_1];

	      TVector3 vec_daup_1;
              vec_daup_1.SetPtEtaPhi(ppt_1, peta_1, pphi_1);

              TVector3 vec_daun_1;
              vec_daun_1.SetPtEtaPhi(npt_1, neta_1, nphi_1);

	      
	      for (unsigned int ntrk_2 = 0; ntrk_2 < V0pt_vect.size(); ntrk_2 ++)
		{
		  if(ntrk_1 == ntrk_2) continue;

		  double ppt_2 = pdau_pt_vect[ntrk_2];
		  double peta_2 = pdau_eta_vect[ntrk_2];
		  double pphi_2 = pdau_phi_vect[ntrk_2];
		  double pchi2_2 = pdau_chi2_vect[ntrk_2];

		  double npt_2 = ndau_pt_vect[ntrk_2];
		  double neta_2 = ndau_eta_vect[ntrk_2];
		  double nphi_2 = ndau_phi_vect[ntrk_2];
		  double nchi2_2 = ndau_chi2_vect[ntrk_2];
		  
		  double p_dphi = fabs(pphi_1 - pphi_2);
		  double n_dphi = fabs(nphi_1 - nphi_2);
		  
		  double p_deta = fabs(peta_1 - peta_2);
		  double n_deta = fabs(neta_1 - neta_2);
		  
		  double p_chi2 = fabs(pchi2_1 - pchi2_2);
		  double n_chi2 = fabs(nchi2_1 - nchi2_2);

		  double p_dpt = fabs(ppt_1 - ppt_2);
                  double n_dpt = fabs(npt_1 - npt_2);

		  TVector3 vec_daup_2;
                  vec_daup_2.SetPtEtaPhi(ppt_2, peta_2, pphi_2);

                  TVector3 vec_daun_2;
                  vec_daun_2.SetPtEtaPhi(npt_2, neta_2, nphi_2);


                  double dcosthetadaun = fabs((vec_daun_1.Dot(vec_daun_2))/(vec_daun_1.Mag()*vec_daun_2.Mag()));
                  double dcosthetadaup = fabs((vec_daup_1.Dot(vec_daup_2))/(vec_daup_1.Mag()*vec_daup_2.Mag()));


                  hdpt_dcostheta_daun_wo->Fill(n_dpt, dcosthetadaun);
                  hdpt_dcostheta_daup_wo->Fill(p_dpt, dcosthetadaup);


                  hdeta_dphi_daun_wo->Fill(n_deta, n_dphi);
                  hdeta_dphi_daup_wo->Fill(p_deta, p_dphi);

		  hdchi2_daup->Fill(p_chi2);
                  hdchi2_daun->Fill(n_chi2);
		 
		  //if ((p_deta == dau_etaphi_ && p_dphi == dau_etaphi_) || (n_deta == dau_etaphi_ && n_dphi == dau_etaphi_))
		  if ((p_chi2 <= dau_cut) || (n_chi2 <= dau_cut))
		    {
		      is_dcut = kTRUE;
		      break;
		    }
		  else is_dcut = kFALSE;
		}
	      if (!is_dcut)
		{
		  AssignpTbins(V0pt_pure, V0eta_pure, V0phi_pure, V0id_pure, V0mass_pure, ppt_1, peta_1, pphi_1, npt_1, neta_1, nphi_1, idx, ctbin);
		}
	      
	    }// vector particle loop end----------------

	  
	  V0pt_vect.clear();
	  V0eta_vect.clear();
	  V0phi_vect.clear();
	  V0id_vect.clear();
	  V0mass_vect.clear();
	  
	  pdau_pt_vect.clear();
	  pdau_eta_vect.clear();
	  pdau_phi_vect.clear();
	  pdau_chi2_vect.clear();
	  
	  ndau_pt_vect.clear();
	  ndau_eta_vect.clear();
	  ndau_phi_vect.clear();
	  ndau_chi2_vect.clear();

	  evtcount[ctbin].push_back(idx);
	  
	  idx ++;
	  
	} // event loop end-------------
      
      AnalyseEvents(vzvtx, evtcount, cntbin, evtno, runno, vect_V0id_trg, vect_V0id_ass, V0_vect_trg, daup_vect_trg, daun_vect_trg, V0_vect_ass, daup_vect_ass, daun_vect_ass);
      

      for(Int_t z = 0; z < NBINS_cent; z++)
	{
	  vzvtx[z].clear();
	  cntbin[z].clear();
	  evtno[z].clear();
	  runno[z].clear();
	  evtcount[z].clear();
	  
	  for(Int_t k = 0; k < iiev; k++)
	    
	    {
	      vect_V0id_trg[z][k].clear();
              V0_vect_trg[z][k].clear();
              daup_vect_trg[z][k].clear();
              daun_vect_trg[z][k].clear();

              vect_V0id_ass[z][k].clear();
              V0_vect_ass[z][k].clear();
              daup_vect_ass[z][k].clear();
              daun_vect_ass[z][k].clear();


	    }
	}
      
      delete evt_tree;
      delete V0_tree;
      file->Close();
      delete file;
      
      
    }//text file loop------------------------
  
   TFile *fout = new TFile("/home/matapita/cernbox/V0_analysis/plot_Anarootfile/new_plot/Final_dec8_2020/new_systematics/kshort/BDT/new_rootfiles_july6/f2mass_check_HM0_1_2_3_Oct2.root","recreate");
  
  TDirectory *Devent_hist = fout->mkdir("Event_hists");
  Devent_hist->cd();

  hzvtx->Write();
  hcentbin->Write();
  hcentbin_dist->Write();
  for (int iic = 0; iic< NBINS_cent; iic++)
    {
      hcent_dist[iic]->Write();
    }
      
  Devent_hist->cd();

  TDirectory *DV0_hist = fout->mkdir("V0_hists");
  DV0_hist->cd();
  
  h_V0mass->Write();
  h_V0mass_pure->Write();
  h_V0pt->Write();
  h_V0kt->Write();
  h_V0mt->Write();
  h_V0eta->Write();
  h_V0phi->Write();
  h_V0rapidity->Write();
  harmen->Write();
  hpip1->Write();
  hpip2->Write();
  hdptcostheta_daup->Write();
  hdptcostheta_daun->Write();
  hdpt_dcostheta_daun_wo->Write();
  hdpt_dcostheta_daup_wo->Write();
  hdeta_dphi_daun_wo->Write();
  hdeta_dphi_daup_wo->Write();
  hdchi2_daun->Write();
  hdchi2_daup->Write();
  
  for (int iic = 0; iic< NBINS_cent; iic++)
    {
      h_V0kt_[iic]->Write();
      h_V0mt_[iic]->Write();
      hf2_sig_[iic]->Write();
      hf2_mix_[iic]->Write();
    }
  hf2_sig->Write();
  hf2_mix->Write();
  
  DV0_hist->cd();
  
  TDirectory *DV0mass_hist = fout->mkdir("V0mass_hists");
  DV0mass_hist->cd();
  
  for (int i = 0; i< NBINS_cent; i++)
    {
      hV0mass_KS_cent_[i]->Write();

      for (int j = 0; j< NBINS_pt; j++)
        {
          hV0mass_KS_[i][j]->Write();
        }
    }
  
  DV0mass_hist->cd();

  TDirectory *Dqinvmass_hist = fout->mkdir("qinvV0mass_hists");
  Dqinvmass_hist->cd();
  
  hV0mass_KS_qinv0p8->Write();
  hV0mass_KS_qinv1p0->Write();
  hV0mass_KS_qinv1p5->Write();
  hV0mass_KS_qinv2p0->Write();
  for (int ii = 0; ii< NBINS_cent; ii++)
    {
      hV0mass_KS_cent_qinv0p8_[ii]->Write();
      hV0mass_KS_cent_qinv1p0_[ii]->Write();
      hV0mass_KS_cent_qinv1p5_[ii]->Write();
      hV0mass_KS_cent_qinv2p0_[ii]->Write();

      /*
      for(int iqbin = 0; iqbin <NBINS_qbin; iqbin++)
	{
	  hV0mass_KS_cbin_qbin_[ii][iqbin]->Write();
	}
      */
    }
  for(int iqbin = 0; iqbin <NBINS_qbin; iqbin++)
    {
      hV0mass_KS_cbin_qbin_[iqbin]->Write();
    }
  
  Dqinvmass_hist->cd();

  
  TDirectory *DV0corr_hist = fout->mkdir("V0corr_hists");
  DV0corr_hist->cd();

  hsig_KS_cent_qinv->Write();
  hmix_KS_cent_qinv->Write();
  
  hsig_daup->Write();
  hsig_daun->Write();
  
  hsig_KS->Write();
  hsig_KS_rotate->Write();
  hsig_KS_invert->Write();
  
  hmix_KS->Write();
  hmix_KS_rotate->Write();
  hmix_KS_invert->Write();
  
  hsig_KS_AB->Write();
  hsig_KS_left_AB->Write();
  hsig_KS_right_AB->Write();
  hsig_KS_rotate_AB->Write();
  hsig_KS_invert_AB->Write();

  hmix_KS_AB->Write();
  hmix_KS_left_AB->Write();
  hmix_KS_right_AB->Write();
  hmix_KS_rotate_AB->Write();
  hmix_KS_invert_AB->Write();

  hsig_KS_BB->Write();
  hsig_KS_BB1->Write();
  hsig_KS_BB2->Write();
  hsig_KS_left_BB->Write();
  hsig_KS_right_BB->Write();
  hsig_KS_rotate_BB->Write();
  hsig_KS_invert_BB->Write();

  hmix_KS_BB->Write();
  hmix_KS_BB1->Write();
  hmix_KS_BB2->Write();
  hmix_KS_left_BB->Write();
  hmix_KS_right_BB->Write();
  hmix_KS_rotate_BB->Write();
  hmix_KS_invert_BB->Write();

  for(int y = 0; y < NBINS_cent; y++)
    {
      hsig_KS_[y]->Write();
      hsig_KS_rotate_[y]->Write();
      hsig_KS_invert_[y]->Write();
      
      hmix_KS_[y]->Write();
      hmix_KS_rotate_[y]->Write();
      hmix_KS_invert_[y]->Write();

      hsig_KS_AB_[y]->Write();
      hsig_KS_left_AB_[y]->Write();
      hsig_KS_right_AB_[y]->Write();
      hsig_KS_rotate_AB_[y]->Write();
      hsig_KS_invert_AB_[y]->Write();

      hmix_KS_AB_[y]->Write();
      hmix_KS_left_AB_[y]->Write();
      hmix_KS_right_AB_[y]->Write();
      hmix_KS_rotate_AB_[y]->Write();
      hmix_KS_invert_AB_[y]->Write();

      hsig_KS_BB_[y]->Write();
      hsig_KS_BB1_[y]->Write();
      hsig_KS_BB2_[y]->Write();
      hsig_KS_left_BB_[y]->Write();
      hsig_KS_right_BB_[y]->Write();
      hsig_KS_rotate_BB_[y]->Write();
      hsig_KS_invert_BB_[y]->Write();

      hmix_KS_BB_[y]->Write();
      hmix_KS_BB1_[y]->Write();
      hmix_KS_BB2_[y]->Write();
      hmix_KS_left_BB_[y]->Write();
      hmix_KS_right_BB_[y]->Write();
      hmix_KS_rotate_BB_[y]->Write();
      hmix_KS_invert_BB_[y]->Write();
    }

  DV0corr_hist->cd();

  TDirectory *DV0corr_ktbin_hist = fout->mkdir("V0corr_ktbin_hists");
  DV0corr_ktbin_hist->cd();

  for(int z = 0; z < NBINS_kt; z++)
    {
      hsig_KS_kt_[z]->Write();
      hmix_KS_kt_[z]->Write();  
    }
  
  for(int y = 0; y < NBINS_cent; y++)
    {
      for(int z = 0; z < NBINS_kt; z++)
	{
	  hsig_KS_ckt_[y][z]->Write();
	  hmix_KS_ckt_[y][z]->Write();
	}
    }

  DV0corr_ktbin_hist->cd();
  
  fout->Write();
  fout->Close();
  delete fout;
  
  auto stop = high_resolution_clock::now();
  
  auto duration = duration_cast<minutes>(stop - start); 
  
  cout << "Total time taken: "
       << duration.count() << "minutes" << endl; 
}//main function-----------------------------

//----------------------------------------------------------------------------
void AssignpTbins(double V0pt, double V0eta, double V0phi, double V0id, double V0mass, double pdau_pt, double pdau_eta, double pdau_phi, double ndau_pt, double ndau_eta, double ndau_phi, int idx, int ctbin)
  
{

  h_V0mass_pure->Fill(V0mass);
  hV0mass_KS_cent_[ctbin]->Fill(V0mass);
   
  TLorentzVector V0vector;
  V0vector.SetPtEtaPhiM(V0pt, V0eta, V0phi, V0mass);
  //std::cout<<"eta is = "<<V0vector.Eta()<<std::endl;
  TLorentzVector daupvector;
  daupvector.SetPtEtaPhiM(pdau_pt, pdau_eta, pdau_phi, 0.139570 );
  
  TLorentzVector daunvector;
  daunvector.SetPtEtaPhiM(ndau_pt, ndau_eta, ndau_phi, 0.139570 );
  
  vect_V0id_trg[ctbin][idx].push_back(V0id);
  V0_vect_trg[ctbin][idx].push_back(V0vector);
  daup_vect_trg[ctbin][idx].push_back(daupvector);
  daun_vect_trg[ctbin][idx].push_back(daunvector);

  vect_V0id_ass[ctbin][idx].push_back(V0id);
  V0_vect_ass[ctbin][idx].push_back(V0vector);
  daup_vect_ass[ctbin][idx].push_back(daupvector);
  daun_vect_ass[ctbin][idx].push_back(daunvector);
}
//----------------------------------------------------------------------------
void AnalyseEvents(std::vector<double> vzvtx[NBINS_cent],std::vector<int> evtcount[NBINS_cent], std::vector<float> cntbin[NBINS_cent], std::vector<double> evtno[NBINS_cent], std::vector<double> runno[NBINS_cent], std::vector<double> vect_V0id_trg[NBINS_cent][iiev], std::vector<double> vect_V0id_ass[NBINS_cent][iiev], std:: vector<TLorentzVector> V0_vect_trg[NBINS_cent][iiev], std:: vector<TLorentzVector> daup_vect_trg[NBINS_cent][iiev], std:: vector<TLorentzVector> daun_vect_trg[NBINS_cent][iiev], std:: vector<TLorentzVector> V0_vect_ass[NBINS_cent][iiev], std:: vector<TLorentzVector> daup_vect_ass[NBINS_cent][iiev], std:: vector<TLorentzVector> daun_vect_ass[NBINS_cent][iiev])
  
{
  
  for(unsigned int ic = 0; ic < NBINS_cent; ic++)
    {
      //std::cout<< "Total of " << vzvtx[ic].size() << " events are selected in"<<"  " << ic <<" "<<"centbin"<< std::endl;
      //std::cout<< "Total of " << evtcount[ic].size() << " events are selected in"<<"  " << ic <<" "<<"centbin"<< std::endl;

      //std::cout<< "Started running correlation analysis in"<<"  "<< ic << "th centbin" <<std::endl;
      
      for(unsigned int iev = 0; iev < evtcount[ic].size(); iev++)
        {
          unsigned int ntrgsize_i =  V0_vect_trg[ic][evtcount[ic][iev]].size();
          unsigned int nasssize_i =  V0_vect_ass[ic][evtcount[ic][iev]].size();

	  if (ntrgsize_i != nasssize_i)
            {
              std::cout<< "ntrgsize_i and nasssize_i are not same" << std::endl;
              continue;
            }
	  
          if (ntrgsize_i < 2 || nasssize_i < 2) continue;
	  
	  hzvtx->Fill(vzvtx[ic][iev]);
	  hcentbin->Fill(cntbin[ic][iev]);
	  hcent_dist[ic]->Fill(cntbin[ic][iev]);
	  
	  SignalTwoPartCorr(evtcount[ic][iev], ic, cntbin[ic][iev]);
	  //TwoPartCorr(evtcount[ic][iev], ic);
	  

	  unsigned int mixstart = iev +1;
          unsigned int mixend = evtcount[ic].size();

	  
          /*if (iev < (unsigned int)(evtcount[ic].size()*0.6))
            {
              mixstart = iev+1;
              mixend = evtcount[ic].size();
            }
	  
          else if (iev >= (unsigned int)(evtcount[ic].size()*0.6))
            {
              mixstart = 0;
              mixend = (unsigned int)(evtcount[ic].size());
            }
	  */
	  
          Int_t nmix=0;
	  
          for( unsigned int jev = mixstart; jev < mixend; jev++ )
	    {
	      if(iev == jev) continue;
	      if(evtcount[ic][iev] == evtcount[ic][jev]) continue;
	      
              double deltazvtx = vzvtx[ic][iev]-vzvtx[ic][jev];
	      
              if(fabs(deltazvtx) >= 2.0) continue;
	      
              if ((evtno[ic][iev] == evtno[ic][jev]) && (runno[ic][iev] == runno[ic][jev])) continue;
	      
	      int cntbin_i = GetCent(cntbin[ic][iev]);
              int cntbin_j = GetCent(cntbin[ic][jev]);

	      if(cntbin_i != cntbin_j) continue;
	      
	      unsigned int ntrgsize_j =  V0_vect_trg[ic][evtcount[ic][iev]].size();
              unsigned int nasssize_j =  V0_vect_ass[ic][evtcount[ic][jev]].size();

	      
	      if (ntrgsize_j ==0 || nasssize_j ==0 ) continue;
	      //if( nasssize_j <= 1 && ntrgsize_j <= 1) continue;
	      
	      nmix++;
              if (nmix >= bkgFactor) break;
	      
	      MixedTwoPartCorr( evtcount[ic][iev], evtcount[ic][jev], ic, cntbin[ic][iev]);
	      
            }
	  //std::cout<<"nmix for the "<<iev<<" th event is: "<<nmix<<std::endl;
	}
      //std::cout<< "Finished running correlation analysis in"<<"  "<< ic << "th centbin" <<std::endl;
    }
}
//----------------------------------------------------------------------------
void SignalTwoPartCorr(unsigned int iev, unsigned int ic, float centbin)

{
  
  unsigned int ntrgsize =  V0_vect_trg[ic][iev].size();
  unsigned int nasssize =  V0_vect_ass[ic][iev].size();

  for( unsigned int ntrg = 0; ntrg < ntrgsize; ntrg++ )
    {
           
      TLorentzVector pvector_daup_trg = daup_vect_trg[ic][iev][ntrg];
      TLorentzVector pvector_daun_trg = daun_vect_trg[ic][iev][ntrg];
      int V0id_trg = vect_V0id_trg[ic][iev][ntrg];
      
      TLorentzVector pvector_V0_trg = V0_vect_trg[ic][iev][ntrg];
      
      double px_daup_trg = pvector_daup_trg.Px();
      double py_daup_trg = pvector_daup_trg.Py();
      double pz_daup_trg = pvector_daup_trg.Pz();
      double E_daup_trg = pvector_daup_trg.E();
      double Pt_daup_trg = pvector_daup_trg.Pt();
      double Eta_daup_trg = pvector_daup_trg.Eta();
      double Phi_daup_trg = pvector_daup_trg.Phi();
      
      double px_daun_trg = pvector_daun_trg.Px();
      double py_daun_trg = pvector_daun_trg.Py();
      double pz_daun_trg = pvector_daun_trg.Pz();
      double E_daun_trg = pvector_daun_trg.E();
      double Pt_daun_trg = pvector_daun_trg.Pt();
      double Eta_daun_trg = pvector_daun_trg.Eta();
      double Phi_daun_trg = pvector_daun_trg.Phi();

      TVector3 daup_trg(px_daup_trg, py_daup_trg, pz_daup_trg);
      TVector3 daun_trg(px_daun_trg, py_daun_trg, pz_daun_trg);
      
      double px_V0_trg = pvector_V0_trg.Px();
      double py_V0_trg = pvector_V0_trg.Py();
      double pz_V0_trg = pvector_V0_trg.Pz();
      double E_V0_trg = pvector_V0_trg.E();
      double M_V0_trg = pvector_V0_trg.M();
      double pt_V0_trg = pvector_V0_trg.Pt();
      double eta_V0_trg = pvector_V0_trg.Eta();
      double phi_V0_trg = pvector_V0_trg.Phi();

      //std::cout<<"px_V0_trg is: "<<px_V0_trg<<"  "<<"py_V0_trg is: "<<py_V0_trg<<"  "<<"pz_V0_trg is: "<<pz_V0_trg<<"  "<<"pt_V0_trg is: "<<pt_V0_trg<<std::endl;
      
      for( unsigned int nass = ntrg+1 ; nass < nasssize; nass++ )
	{
	  TLorentzVector pvector_daup_ass = daup_vect_ass[ic][iev][nass];
	  TLorentzVector pvector_daun_ass = daun_vect_ass[ic][iev][nass];
	  int V0id_ass = vect_V0id_ass[ic][iev][nass];
	  
	  TLorentzVector pvector_V0_ass = V0_vect_ass[ic][iev][nass];
      
	  
	  double px_daup_ass = pvector_daup_ass.Px();
	  double py_daup_ass = pvector_daup_ass.Py();
	  double pz_daup_ass = pvector_daup_ass.Pz();
	  double E_daup_ass = pvector_daup_ass.E();
	  double Pt_daup_ass = pvector_daup_ass.Pt();
	  double Eta_daup_ass = pvector_daup_ass.Eta();
	  double Phi_daup_ass = pvector_daup_ass.Phi();
      
	  double px_daun_ass = pvector_daun_ass.Px();
	  double py_daun_ass = pvector_daun_ass.Py();
	  double pz_daun_ass = pvector_daun_ass.Pz();
	  double E_daun_ass = pvector_daun_ass.E();
	  double Pt_daun_ass = pvector_daun_ass.Pt();
	  double Eta_daun_ass = pvector_daun_ass.Eta();
	  double Phi_daun_ass = pvector_daun_ass.Phi();
	  
	  TVector3 daup_ass(px_daup_ass, py_daup_ass, pz_daup_ass);
	  TVector3 daun_ass(px_daun_ass, py_daun_ass, pz_daun_ass);

	  
	  double px_V0_ass = (pvector_V0_ass.Px());
	  double py_V0_ass = (pvector_V0_ass.Py());
	  double pz_V0_ass = pvector_V0_ass.Pz();
	  double E_V0_ass = pvector_V0_ass.E();
	  double M_V0_ass = pvector_V0_ass.M();
	  double pt_V0_ass = pvector_V0_ass.Pt();
	  double eta_V0_ass = pvector_V0_ass.Eta();
	  double phi_V0_ass = pvector_V0_ass.Phi();

	  if (px_V0_trg == px_V0_ass && py_V0_trg == py_V0_ass && pz_V0_trg == pz_V0_ass &&  pt_V0_trg == pt_V0_ass && eta_V0_trg == eta_V0_ass && phi_V0_trg == phi_V0_ass && M_V0_trg == M_V0_ass) continue;

	  //double kt = 0.5*(fabs(pt_V0_trg + pt_V0_ass));

	  TLorentzVector total_V0 = pvector_V0_trg + pvector_V0_ass ;
          double kt_V0 = 0.5*(fabs(total_V0.Pt()));
          double mu = (2*(M_V0_trg * M_V0_ass))/((M_V0_trg + M_V0_ass));
          double mt = TMath::Sqrt(pow(mu,2) + pow(kt_V0,2));
          h_V0mt->Fill(mt);
          h_V0kt->Fill(kt_V0);


	  if(kt_V0 >= maxkt || kt_V0 < 0) continue;
	  
	  int iktbin = GetKtbin(kt_V0);
	  //if (iktbin == 4) continue;

	  double Dpt_p = fabs(Pt_daup_trg - Pt_daup_ass);
	  double Dpt_n = fabs(Pt_daun_trg - Pt_daun_ass);
	  double costheta_daup = fabs(daup_trg.Dot(daup_ass))/((daup_trg.Mag())*(daup_ass.Mag()));
	  double costheta_daun = fabs(daun_trg.Dot(daun_ass))/((daun_trg.Mag())*(daun_ass.Mag()));
	  
	  hdptcostheta_daup->Fill(Dpt_p, costheta_daup);
	  hdptcostheta_daun->Fill(Dpt_n, costheta_daun);
	  
	  double Qinv_daup = GetQInv(px_daup_ass,py_daup_ass,pz_daup_ass,E_daup_ass,px_daup_trg,py_daup_trg,pz_daup_trg,E_daup_trg);
	  double Qinv_daun = GetQInv(px_daun_ass,py_daun_ass,pz_daun_ass,E_daun_ass,px_daun_trg,py_daun_trg,pz_daun_trg,E_daun_trg);
	  
	  double Qinv_V0 = GetQInv(px_V0_ass,py_V0_ass,pz_V0_ass,E_V0_ass,px_V0_trg,py_V0_trg,pz_V0_trg,E_V0_trg);

	  double Qinv_V0_rotate = GetQInv(-px_V0_ass,-py_V0_ass,pz_V0_ass,E_V0_ass,px_V0_trg,py_V0_trg,pz_V0_trg,E_V0_trg);

	  double Qinv_V0_invert = GetQInv(-px_V0_ass,-py_V0_ass,-pz_V0_ass,E_V0_ass,px_V0_trg,py_V0_trg,pz_V0_trg,E_V0_trg);

	  double f2mass = GetInvMass(E_V0_ass, px_V0_ass, py_V0_ass, pz_V0_ass, E_V0_trg, px_V0_trg, py_V0_trg, pz_V0_trg);
	  
	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && (M_V0_ass >= mass_sig_min && M_V0_ass <= mass_sig_max))
	    {
	      hf2_sig_[ic]->Fill(f2mass);
	      hf2_sig->Fill(f2mass);
	      
	      h_V0kt_[ic]->Fill(kt_V0);
	      h_V0mt_[ic]->Fill(mt);
	      
	      hsig_KS_cent_qinv->Fill(centbin, Qinv_V0);
	      
	      hsig_daup->Fill(Qinv_daup);
	      hsig_daun->Fill(Qinv_daun);
	      hsig_KS->Fill(Qinv_V0);
	      hsig_KS_rotate->Fill(Qinv_V0_rotate);
              hsig_KS_invert->Fill(Qinv_V0_invert);
	      
	      hsig_KS_[ic]->Fill(Qinv_V0);
	      hsig_KS_rotate_[ic]->Fill(Qinv_V0_rotate);
	      hsig_KS_invert_[ic]->Fill(Qinv_V0_invert);

	      hsig_KS_kt_[iktbin]->Fill(Qinv_V0);
	      hsig_KS_ckt_[ic][iktbin]->Fill(Qinv_V0);
	    }
	  
	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && (( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high ) || ( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high ))) 
	    {
	      hsig_KS_AB->Fill(Qinv_V0);
	      hsig_KS_rotate_AB->Fill(Qinv_V0_rotate);
              hsig_KS_invert_AB->Fill(Qinv_V0_invert);
	      
	      hsig_KS_AB_[ic]->Fill(Qinv_V0);
	      hsig_KS_rotate_AB_[ic]->Fill(Qinv_V0_rotate);
	      hsig_KS_invert_AB_[ic]->Fill(Qinv_V0_invert);
	    }
	  
	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && ( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high))	    
            {
	      hsig_KS_left_AB->Fill(Qinv_V0);
	      hsig_KS_left_AB_[ic]->Fill(Qinv_V0);
	    }

	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && ( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high))
	    {
	      hsig_KS_right_AB->Fill(Qinv_V0);
	      hsig_KS_right_AB_[ic]->Fill(Qinv_V0);
	    }

	  if((( M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high ) && ( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high)) || ((M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high )&&( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high )))
	    {
	      hsig_KS_BB1->Fill(Qinv_V0);
              hsig_KS_BB1_[ic]->Fill(Qinv_V0);
	    }
	  
	  if ((( M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high ) || ( M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high )) && (( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high )||( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high )))
	    {
	      hsig_KS_BB2->Fill(Qinv_V0);
	      hsig_KS_BB2_[ic]->Fill(Qinv_V0);
	    }
	  
	  if (((M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high) && (( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high )||( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high))) || ((M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high) && ( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high)))
            {
	      hsig_KS_BB->Fill(Qinv_V0);
	      hsig_KS_rotate_BB->Fill(Qinv_V0_rotate);
              hsig_KS_invert_BB->Fill(Qinv_V0_invert);
	      
	      hsig_KS_BB_[ic]->Fill(Qinv_V0);
	      hsig_KS_rotate_BB_[ic]->Fill(Qinv_V0_rotate);
	      hsig_KS_invert_BB_[ic]->Fill(Qinv_V0_invert);
	    }

	  if ((M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high )&&( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high ))
	    {
	      hsig_KS_left_BB->Fill(Qinv_V0);
	      hsig_KS_left_BB_[ic]->Fill(Qinv_V0);
            }
	  
	  if ((M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high )&&( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high ))
	    {
	      hsig_KS_right_BB->Fill(Qinv_V0);
	      hsig_KS_right_BB_[ic]->Fill(Qinv_V0);
	    }

	}
    }
}
//----------------------------------------------------------------------------
void MixedTwoPartCorr(unsigned int iev, unsigned int jev, int ic, float centbin)

{
  
  unsigned int ntrgsize =  V0_vect_trg[ic][iev].size();
  unsigned int nasssize =  V0_vect_ass[ic][jev].size();

  for( unsigned int ntrg = 0; ntrg < ntrgsize; ntrg++ )
    {
           
      TLorentzVector pvector_V0_trg = V0_vect_trg[ic][iev][ntrg];
      int V0id_trg = vect_V0id_trg[ic][iev][ntrg];
      
      double px_V0_trg = pvector_V0_trg.Px();
      double py_V0_trg = pvector_V0_trg.Py();
      double pz_V0_trg = pvector_V0_trg.Pz();
      double E_V0_trg = pvector_V0_trg.E();
      double M_V0_trg = pvector_V0_trg.M();
      double pt_V0_trg = pvector_V0_trg.Pt();
      double eta_V0_trg = pvector_V0_trg.Eta();
      double phi_V0_trg = pvector_V0_trg.Phi();
      
      for( unsigned int nass = 0 ; nass < nasssize; nass++ )
	{
	  
	  TLorentzVector pvector_V0_ass = V0_vect_ass[ic][jev][nass];
	  int V0id_ass = vect_V0id_ass[ic][jev][nass];

	  double px_V0_ass = (pvector_V0_ass.Px());
	  double py_V0_ass = (pvector_V0_ass.Py());
	  double pz_V0_ass = pvector_V0_ass.Pz();
	  double E_V0_ass = pvector_V0_ass.E();
	  double M_V0_ass = pvector_V0_ass.M();
	  double pt_V0_ass = pvector_V0_ass.Pt();
	  double eta_V0_ass = pvector_V0_ass.Eta();
	  double phi_V0_ass = pvector_V0_ass.Phi();

	  if (px_V0_trg == px_V0_ass && py_V0_trg == py_V0_ass && pz_V0_trg == pz_V0_ass &&  pt_V0_trg == pt_V0_ass && eta_V0_trg == eta_V0_ass && phi_V0_trg == phi_V0_ass && M_V0_trg == M_V0_ass) continue;
	  
	  //double kt = 0.5*(fabs(pt_V0_trg + pt_V0_ass));
	  TLorentzVector total_V0 = pvector_V0_trg + pvector_V0_ass ;
          double kt_V0 = 0.5*(fabs(total_V0.Pt()));

	  if(kt_V0 >= maxkt || kt_V0 < 0) continue;
	  	  	  
	  int iktbin = GetKtbin(kt_V0);
	  //if (iktbin == 4) continue;
	  
	  double Qinv_V0 = GetQInv(px_V0_ass,py_V0_ass,pz_V0_ass,E_V0_ass,px_V0_trg,py_V0_trg,pz_V0_trg,E_V0_trg);

	  double Qinv_V0_rotate = GetQInv(-px_V0_ass,-py_V0_ass,pz_V0_ass,E_V0_ass,px_V0_trg,py_V0_trg,pz_V0_trg,E_V0_trg);
	  
	  double Qinv_V0_invert = GetQInv(-px_V0_ass,-py_V0_ass,-pz_V0_ass,E_V0_ass,px_V0_trg,py_V0_trg,pz_V0_trg,E_V0_trg);

	  double hf2mass =  GetInvMass(E_V0_ass, px_V0_ass, py_V0_ass, pz_V0_ass, E_V0_trg, px_V0_trg, py_V0_trg, pz_V0_trg);
	  
	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && (M_V0_ass >= mass_sig_min && M_V0_ass <= mass_sig_max))	  
	    {

	      hf2_mix_[ic]->Fill(hf2mass);
	      hf2_mix->Fill(hf2mass);
	      
	      hmix_KS_cent_qinv->Fill(centbin, Qinv_V0);
	      
	      hmix_KS->Fill(Qinv_V0);
              hmix_KS_rotate->Fill(Qinv_V0_rotate);
              hmix_KS_invert->Fill(Qinv_V0_invert);
	      
	      hmix_KS_[ic]->Fill(Qinv_V0);
	      hmix_KS_rotate_[ic]->Fill(Qinv_V0_rotate);
	      hmix_KS_invert_[ic]->Fill(Qinv_V0_invert);

	      hmix_KS_kt_[iktbin]->Fill(Qinv_V0);
	      hmix_KS_ckt_[ic][iktbin]->Fill(Qinv_V0);
	      
	    }

	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && (( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high ) || ( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high ))) 
	    {
	      hmix_KS_AB->Fill(Qinv_V0);
	      hmix_KS_rotate_AB->Fill(Qinv_V0_rotate);
              hmix_KS_invert_AB->Fill(Qinv_V0_invert);
	      
	      hmix_KS_AB_[ic]->Fill(Qinv_V0);
	      hmix_KS_rotate_AB_[ic]->Fill(Qinv_V0_rotate);
	      hmix_KS_invert_AB_[ic]->Fill(Qinv_V0_invert);
	    }
	  
	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && ( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high))
            {
		hmix_KS_left_AB->Fill(Qinv_V0);
                hmix_KS_left_AB_[ic]->Fill(Qinv_V0);
	    }

	  if ((M_V0_trg >= mass_sig_min && M_V0_trg <= mass_sig_max) && ( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high))
	    {
	      hmix_KS_right_AB->Fill(Qinv_V0);
	      hmix_KS_right_AB_[ic]->Fill(Qinv_V0);
	    }

	  if((( M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high ) && ( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high)) || ((M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high )&&( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high )))
	    {
	      hmix_KS_BB1->Fill(Qinv_V0);
              hmix_KS_BB1_[ic]->Fill(Qinv_V0);
	    }
	  
	  if ((( M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high ) || ( M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high )) && (( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high )||( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high )))
	    {
	      hmix_KS_BB2->Fill(Qinv_V0);
	      hmix_KS_BB2_[ic]->Fill(Qinv_V0);
	    }
	  
	  if (((M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high) && (( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high )||( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high))) || ((M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high) && ( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high)))
	    {
	      hmix_KS_BB->Fill(Qinv_V0);
	      hmix_KS_rotate_BB->Fill(Qinv_V0_rotate);
              hmix_KS_invert_BB->Fill(Qinv_V0_invert);
	      
	      hmix_KS_BB_[ic]->Fill(Qinv_V0);
	      hmix_KS_rotate_BB_[ic]->Fill(Qinv_V0_rotate);
	      hmix_KS_invert_BB_[ic]->Fill(Qinv_V0_invert);
	      
	    }
	  
	  if ((M_V0_trg >= mass_leftside_low && M_V0_trg <= mass_leftside_high )&&( M_V0_ass >= mass_leftside_low && M_V0_ass <= mass_leftside_high ))
	    
	    {
	      hmix_KS_left_BB->Fill(Qinv_V0);
	      hmix_KS_left_BB_[ic]->Fill(Qinv_V0);
	    }
	  
	  if ((M_V0_trg >= mass_rightside_low && M_V0_trg <= mass_rightside_high )&&( M_V0_ass >= mass_rightside_low && M_V0_ass <= mass_rightside_high ))
	    {
	      hmix_KS_right_BB->Fill(Qinv_V0);
	      hmix_KS_right_BB_[ic]->Fill(Qinv_V0);
	    }
	}
    }
}
//----------------------------------------------------------------------------
void TwoPartCorr(unsigned int iev, unsigned int ic)
{

  unsigned int ntrgsize =  V0_vect_trg[ic][iev].size();
  unsigned int nasssize =  V0_vect_ass[ic][iev].size();

  for( unsigned int ntrg = 0; ntrg < ntrgsize; ntrg++ )
    {

      int V0id_trg = vect_V0id_trg[ic][iev][ntrg];

      TLorentzVector pvector_V0_trg = V0_vect_trg[ic][iev][ntrg];

      double px_V0_trg = pvector_V0_trg.Px();
      double py_V0_trg = pvector_V0_trg.Py();
      double pz_V0_trg = pvector_V0_trg.Pz();
      double E_V0_trg = pvector_V0_trg.E();
      double M_V0_trg = pvector_V0_trg.M();
      double pt_V0_trg = pvector_V0_trg.Pt();
      double eta_V0_trg = pvector_V0_trg.Eta();
      double phi_V0_trg = pvector_V0_trg.Phi();


      int ii_1p0 = 0, ii_0p8 = 0, ii_1p5 = 0; int ii_2p0 = 0;
      
      int Count_qinv[NBINS_qbin] ={};
      
      for( unsigned int nass = 0 ; nass < nasssize; nass++ )
        {
          if (ntrg == nass) continue;

          int V0id_ass = vect_V0id_ass[ic][iev][nass];

          TLorentzVector pvector_V0_ass = V0_vect_ass[ic][iev][nass];
	  double px_V0_ass = pvector_V0_ass.Px();
          double py_V0_ass = pvector_V0_ass.Py();
          double pz_V0_ass = pvector_V0_ass.Pz();
          double E_V0_ass = pvector_V0_ass.E();
          double M_V0_ass = pvector_V0_ass.M();
          double pt_V0_ass = pvector_V0_ass.Pt();
	  double eta_V0_ass = pvector_V0_ass.Eta();
	  double phi_V0_ass = pvector_V0_ass.Phi();


	  if (px_V0_trg == px_V0_ass && py_V0_trg == py_V0_ass && pz_V0_trg == pz_V0_ass &&  pt_V0_trg == pt_V0_ass && eta_V0_trg == eta_V0_ass && phi_V0_trg == phi_V0_ass && M_V0_trg == M_V0_ass) continue;
	  
          double Qinv_V0 = GetQInv(px_V0_ass,py_V0_ass,pz_V0_ass,E_V0_ass,px_V0_trg,py_V0_trg,pz_V0_trg,E_V0_trg);
	  //double kt = 0.5*(fabs(pt_V0_trg + pt_V0_ass));

	  TLorentzVector total_V0 = pvector_V0_trg + pvector_V0_ass ;
          double kt_V0 = 0.5*(fabs(total_V0.Pt()));

	  if(kt_V0 >= maxkt || kt_V0 < 0) continue;

	  if(Qinv_V0 > 4. ) continue;

	  if(V0id_trg == 310 && V0id_ass == 310)
	    {
	      int qbin = Getqinvbin(Qinv_V0);
	      Count_qinv[qbin] += 1; 	  
	    }
	  
	  if(Qinv_V0 <= 0.8)
	    {
	      ii_0p8++;
	    }
	  
	  if(Qinv_V0 <= 1.0)
	    {
	      ii_1p0++;
	    }
	  
	  if (Qinv_V0 <= 1.5)
	    {
	      ii_1p5++;
	    }
	  
	  if (Qinv_V0 <= 2.0)
	    {
	      ii_2p0++;
	    }
	}

      if(ii_0p8 > 0 )
	{
	  hV0mass_KS_qinv0p8->Fill(M_V0_trg);
	  hV0mass_KS_cent_qinv0p8_[ic]->Fill(M_V0_trg);
        }

      if(ii_1p0 > 0 )
	{
	  hV0mass_KS_qinv1p0->Fill(M_V0_trg);
	  hV0mass_KS_cent_qinv1p0_[ic]->Fill(M_V0_trg);
        }
      
      if(ii_1p5 > 0 )
	{
	  hV0mass_KS_qinv1p5->Fill(M_V0_trg);
	  hV0mass_KS_cent_qinv1p5_[ic]->Fill(M_V0_trg);
        }
      
      if(ii_2p0 > 0 )
	{
	  hV0mass_KS_qinv2p0->Fill(M_V0_trg);
	  hV0mass_KS_cent_qinv2p0_[ic]->Fill(M_V0_trg);
        }

      
      for(int iqbin = 0; iqbin < NBINS_qbin; iqbin++)
	{
	  if(Count_qinv[iqbin] > 0)
	    {
	      //hV0mass_KS_cbin_qbin_[ic][iqbin]->Fill(M_V0_trg);
	      hV0mass_KS_cbin_qbin_[iqbin]->Fill(M_V0_trg);
	    }
	}
    }
}
//-----------------------------------------------------------------------------
double GetQInv(double px1,double py1,double pz1,double E1,double px2,double py2,double pz2,double E2)
{
  double Qinv;

  Qinv = sqrt(fabs(pow(px1-px2,2)+pow(py1-py2,2)+pow(pz1-pz2,2) - pow(E1-E2,2)));
  return Qinv;
}
//-----------------------------------------------------------------------------
int GetCentbin(Float_t Cent)
{
  //int cbin = hcbin->FindBin(Cent/2.)-1;
  int cbin = (hcbin->FindBin(Cent))-1;
  return cbin;
}

//-----------------------------------------------------------------------------
int GetCent(Float_t Cent)
{
  //int ct = hcent->FindBin(Cent/2.)-1;
  int ct = (hcent->FindBin(Cent))-1;
  return ct;
}
//-----------------------------------------------------------------------------                                     
int GetPtbin(Float_t pt)
{
  int ptbin = (hptbin->FindBin(pt))-1;
  return ptbin;
}
//-----------------------------------------------------------------------------                                     
int GetKtbin(Float_t kt)
{
  int ktbin = (hktbin->FindBin(kt))-1;
  return ktbin;
}
//-----------------------------------------------------------------------------
int Getqinvbin(Float_t qinv)
{
  int qinvbin = (hqinvbin->FindBin(qinv))-1;
  return qinvbin;
}
//-----------------------------------------------------------------------------                                     
double GetInvMass(double e1, double px1, double py1, double pz1, double e2, double px2, double py2, double pz2)
{
  double p1p2 = sqrt(pow((px1+px2),2)+pow((py1+py2),2)+pow((pz1+pz2),2));
  double Invmass = sqrt(pow((e1+e2),2) - pow((p1p2),2));
  return Invmass;
}
//-----------------------------------------------------------------------------                                     
