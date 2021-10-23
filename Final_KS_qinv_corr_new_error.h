#include "TLorentzVector.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TString.h"

const int iiev = 1000000;

//------------eta, pt, kt, centrality binning----------------- 
const Int_t ibin = 4;
Double_t edge_eta[ibin+1] = {-1.2 , -0.4, 0, 0.4, 1.2};
TH1D *h_eta = new TH1D("V0eta","", ibin, edge_eta);

const int NBINS_cent = 5;
//double cbin_edge[NBINS_cent+1]= {0., 20., 40., 60., 80., 100., 120., 200};
//double cbin_edge[NBINS_cent+1]= {0., 10., 30., 50., 80., 100.};
double cbin_edge[NBINS_cent+1]= {0, 20, 60, 100, 160, 200};
TH1D *hcbin  = new TH1D("hcbin", "", NBINS_cent, cbin_edge);
TH1D *hcent  = new TH1D("hcent", "", 20, 0, 200);
TH1D *hcent_dist[NBINS_cent];

const Int_t NBINS_pt = 13;
Double_t edges_pt[NBINS_pt + 1] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 9.0};
TH1D *hptbin  = new TH1D("hptbin", "", NBINS_pt, edges_pt);

const Int_t NBINS_kt = 5;
Double_t edges_kt[NBINS_kt + 1] = {0.0, 1.2, 1.5, 1.8, 2.5, 10.0};
TH1D *hktbin  = new TH1D("hktbin", "", NBINS_kt, edges_kt);

const Int_t NBINS_qbin = 55;
TH1D *hqinvbin = new TH1D("hqinvbin","", NBINS_qbin, 0., 4.4);
//-----------Histograms----------------------------------
TH1D *hzvtx = new TH1D("hzvtx distribution", "", 400 , -20, 20);
TH1D *hzvtx_bin = new TH1D("hzvtxbin", "", 8 , 0, 80);
TH1D *hcentbin  = new TH1D("hcentbin distribution", "", 200 , 0, 200);
TH1D *hcentbin_dist  = new TH1D("hcentbin dist", "", 200 , 0, 200);

TH1D *h_V0pt = new TH1D("V0pt distribution","", 90, 0.0, 9.0);
TH1D *h_V0kt = new TH1D("V0kt distribution","", 95, -0.5, 9.0);
TH1D *h_V0mt = new TH1D("V0mt distribution","", 95, -0.5, 9.0);
TH1D *h_V0eta = new TH1D("V0eta distribution","",480 , -2.4, 2.4);
TH1D *h_V0rapidity = new TH1D("V0rapidity distribution","",480 , -2.4, 2.4);
TH1D *h_V0phi = new TH1D("V0phi distribution","",628 , -TMath::Pi(),  TMath::Pi());
TH1D *h_V0mass = new TH1D("V0mass distribution","", 260 , 0.43, 0.56);
TH1D *h_V0mass_pure = new TH1D("V0mass_pure ","", 260 , 0.43, 0.56);
TH2D *harmen = new TH2D("harmen", "", 200, -1., 1., 300, 0., 0.3 );
TH2D *hdptcostheta_daup = new TH2D("hdptcostheta_daup", "", 100, 0, 0.1, 100, 0.99, 1.);
TH2D *hdptcostheta_daun = new TH2D("hdptcostheta_daun", "", 100, 0, 0.1, 100, 0.99, 1.);
TH2D *hdpt_dcostheta_daun_wo = new TH2D("hdpt_dcostheta_daun_wo", "", 100, 0., 0.1, 100, 0.9, 1. );
TH2D *hdpt_dcostheta_daup_wo = new TH2D("hdpt_dcostheta_daup_wo", "", 100, 0., 0.1, 100, 0.9, 1. );
TH2D *hdeta_dphi_daun_wo = new TH2D("hdeta_dphi_daun_wo", "", 100, 0., 0.1, 100, 0., 0.1 );
TH2D *hdeta_dphi_daup_wo = new TH2D("hdeta_dphi_daup_wo", "", 100, 0., 0.1, 100, 0., 0.1 );
TH1D *hdchi2_daup  = new TH1D("hdchi2_daup", "", 100, 0., 0.1);
TH1D *hdchi2_daun  = new TH1D("hdchi2_daun", "", 100, 0., 0.1);
TH1D* hpip1 = new TH1D("hpip1","", 160, 1.08, 1.16);
TH1D* hpip2 = new TH1D("hpip2","", 160, 1.08, 1.16);
TH1D* hf2_sig = new TH1D("hf2_sig","", 75, 1.2, 1.8);
TH1D* hf2_mix = new TH1D("hf2_mix","", 75, 1.2, 1.8);
TH1D* hf2_sig_[NBINS_cent];
TH1D* hf2_mix_[NBINS_cent];

TH1D* hsig_daup = new TH1D("hsig_daup","",100 ,0, 2) ;
TH1D* hsig_daun = new TH1D("hsig_daun","",100 ,0, 2) ;
TH1D* hV0mass_bais = new TH1D("hV0mass_bais","",500 ,0, 10) ;

TH1D* h_V0kt_[NBINS_cent];
TH1D* h_V0mt_[NBINS_cent];

std::vector< std::vector<TH1D*> > hV0mass_KS_;
//std::vector< std::vector<TH1D*> > hV0mass_KS_cbin_qbin_;
TH1D* hV0mass_KS_cbin_qbin_[NBINS_qbin];
TH1D* hV0mass_KS_cent_[NBINS_cent];

TH1D* hV0mass_KS_cent_qinv0p8_[NBINS_cent];
TH1D* hV0mass_KS_cent_qinv1p0_[NBINS_cent];
TH1D* hV0mass_KS_cent_qinv1p5_[NBINS_cent];
TH1D* hV0mass_KS_cent_qinv2p0_[NBINS_cent];

TH1D *hV0mass_KS_qinv0p8 = new TH1D("hV0mass_KS_qinv0p8","", 260 , 0.43, 0.56);
TH1D *hV0mass_KS_qinv1p0 = new TH1D("hV0mass_KS_qinv1p0","", 260 , 0.43, 0.56);
TH1D *hV0mass_KS_qinv1p5 = new TH1D("hV0mass_KS_qinv1p5","", 260 , 0.43, 0.56);
TH1D *hV0mass_KS_qinv2p0 = new TH1D("hV0mass_KS_qinv2p0","", 260 , 0.43, 0.56);


TH2D* hsig_KS_cent_qinv = new TH2D("hsig_KS_cent_qinv","",200, 0, 200, 500, 0, 10);
TH2D* hmix_KS_cent_qinv = new TH2D("hmix_KS_cent_qinv","",200, 0, 200, 500, 0, 10);

TH1D* hsig_KS = new TH1D("hsig_KS","",500, 0, 10);
TH1D* hsig_KS_rotate = new TH1D("hsig_KS_rotate","",500, 0, 10);
TH1D* hsig_KS_invert = new TH1D("hsig_KS_invert","",500, 0, 10);

TH1D* hmix_KS = new TH1D("hmix_KS","",500, 0, 10);
TH1D* hmix_KS_rotate = new TH1D("hmix_KS_rotate","",500, 0, 10);
TH1D* hmix_KS_invert = new TH1D("hmix_KS_invert","",500, 0, 10);

TH1D* hsig_KS_AB = new TH1D("hsig_KS_AB","",500, 0, 10);
TH1D* hsig_KS_left_AB = new TH1D("hsig_KS_left_AB","",500, 0, 10);
TH1D* hsig_KS_right_AB = new TH1D("hsig_KS_right_AB","",500, 0, 10);
TH1D* hsig_KS_rotate_AB = new TH1D("hsig_KS_rotate_AB","",500, 0, 10);
TH1D* hsig_KS_invert_AB = new TH1D("hsig_KS_invert_AB","",500, 0, 10);

TH1D* hmix_KS_AB = new TH1D("hmix_KS_AB","",500, 0, 10);
TH1D* hmix_KS_left_AB = new TH1D("hmix_KS_left_AB","",500, 0, 10);
TH1D* hmix_KS_right_AB = new TH1D("hmix_KS_right_AB","",500, 0, 10);
TH1D* hmix_KS_rotate_AB = new TH1D("hmix_KS_rotate_AB","",500, 0, 10);
TH1D* hmix_KS_invert_AB = new TH1D("hmix_KS_invert_AB","",500, 0, 10);

TH1D* hsig_KS_BB = new TH1D("hsig_KS_BB","",500, 0, 10);
TH1D* hsig_KS_BB1 = new TH1D("hsig_KS_BB1","",500, 0, 10);
TH1D* hsig_KS_BB2 = new TH1D("hsig_KS_BB2","",500, 0, 10);
TH1D* hsig_KS_left_BB = new TH1D("hsig_KS_left_BB","",500, 0, 10);
TH1D* hsig_KS_right_BB = new TH1D("hsig_KS_right_BB","",500, 0, 10);
TH1D* hsig_KS_rotate_BB = new TH1D("hsig_KS_rotate_BB","",500, 0, 10);
TH1D* hsig_KS_invert_BB = new TH1D("hsig_KS_invert_BB","",500, 0, 10);

TH1D* hmix_KS_BB = new TH1D("hmix_KS_BB","",500, 0, 10);
TH1D* hmix_KS_BB1 = new TH1D("hmix_KS_BB1","",500, 0, 10);
TH1D* hmix_KS_BB2 = new TH1D("hmix_KS_BB2","",500, 0, 10);
TH1D* hmix_KS_left_BB = new TH1D("hmix_KS_left_BB","",500, 0, 10);
TH1D* hmix_KS_right_BB = new TH1D("hmix_KS_right_BB","",500, 0, 10);
TH1D* hmix_KS_rotate_BB = new TH1D("hmix_KS_rotate_BB","",500, 0, 10);
TH1D* hmix_KS_invert_BB = new TH1D("hmix_KS_invert_BB","",500, 0, 10);

TH1D* hsig_KS_[NBINS_cent];
TH1D* hsig_KS_rotate_[NBINS_cent];
TH1D* hsig_KS_invert_[NBINS_cent];

TH1D* hmix_KS_[NBINS_cent];
TH1D* hmix_KS_rotate_[NBINS_cent];
TH1D* hmix_KS_invert_[NBINS_cent];


TH1D* hsig_KS_AB_[NBINS_cent];
TH1D* hsig_KS_left_AB_[NBINS_cent];
TH1D* hsig_KS_right_AB_[NBINS_cent];
TH1D* hsig_KS_rotate_AB_[NBINS_cent];
TH1D* hsig_KS_invert_AB_[NBINS_cent];

TH1D* hmix_KS_AB_[NBINS_cent];
TH1D* hmix_KS_left_AB_[NBINS_cent];
TH1D* hmix_KS_right_AB_[NBINS_cent];
TH1D* hmix_KS_rotate_AB_[NBINS_cent];
TH1D* hmix_KS_invert_AB_[NBINS_cent];

TH1D* hsig_KS_BB_[NBINS_cent];
TH1D* hsig_KS_BB1_[NBINS_cent];
TH1D* hsig_KS_BB2_[NBINS_cent];
TH1D* hsig_KS_left_BB_[NBINS_cent];
TH1D* hsig_KS_right_BB_[NBINS_cent];
TH1D* hsig_KS_rotate_BB_[NBINS_cent];
TH1D* hsig_KS_invert_BB_[NBINS_cent];

TH1D* hmix_KS_BB_[NBINS_cent];
TH1D* hmix_KS_BB1_[NBINS_cent];
TH1D* hmix_KS_BB2_[NBINS_cent];
TH1D* hmix_KS_left_BB_[NBINS_cent];
TH1D* hmix_KS_right_BB_[NBINS_cent];
TH1D* hmix_KS_rotate_BB_[NBINS_cent];
TH1D* hmix_KS_invert_BB_[NBINS_cent];

TH1D* hsig_KS_kt_[NBINS_kt];
TH1D* hmix_KS_kt_[NBINS_kt];

//--------Global variables------------------
std::vector< std::vector<TH1D*> > hsig_KS_ckt_;
std::vector< std::vector<TH1D*> > hmix_KS_ckt_;;

std::vector<double> vzvtx[NBINS_cent];
std::vector<float> cntbin[NBINS_cent];
std::vector<double> evtno[NBINS_cent];
std::vector<double> runno[NBINS_cent];
std::vector<int> evtcount[NBINS_cent];

std::vector<double> vect_V0id_trg[NBINS_cent][iiev];
std:: vector<TLorentzVector> V0_vect_trg[NBINS_cent][iiev];
std:: vector<TLorentzVector> daup_vect_trg[NBINS_cent][iiev];
std:: vector<TLorentzVector> daun_vect_trg[NBINS_cent][iiev];

std::vector<double> vect_V0id_ass[NBINS_cent][iiev];
std:: vector<TLorentzVector> V0_vect_ass[NBINS_cent][iiev];
std:: vector<TLorentzVector> daup_vect_ass[NBINS_cent][iiev];
std:: vector<TLorentzVector> daun_vect_ass[NBINS_cent][iiev];


Float_t V0_mass;
Float_t V0_pt;
Float_t V0_eta;
Float_t V0_rapidity;
Float_t V0_phi;
Float_t V0_id;
Float_t V0_dca;
Float_t V0_mva;
Float_t V0_Chi2;
Float_t V0_lxyz;
Float_t V0_costhetaXYZ;
Float_t V0_vtxDecaySigXYZ;


Float_t daup_pt;
Float_t daup_eta;
Float_t daup_phi;
Float_t daup_chi2;
Float_t daup_pterror;
Float_t daup_nhits;
Float_t daup_NPxlayer;
Float_t daup_DCAsigXY;
Float_t daup_DCAsigZ;

Float_t daun_pt;
Float_t daun_eta;
Float_t daun_phi;
Float_t daun_chi2;
Float_t daun_pterror;
Float_t daun_nhits;
Float_t daun_NPxlayer;
Float_t daun_DCAsigXY;
Float_t daun_DCAsigZ;
Float_t cbin;
//--------------kinematics cuts-------------------
double dau_cut = 0.000001;
double bdt_cut = 0.15;
double maxkt = 2.5;
double pTmax = 8.5;
double pTmin = 1.0;
//double pTmin = 0.4;
double Rapmax = 1.0;
double Rapmin = -1.0;
double Massmax = 0.56;
double Massmin = 0.43;
const Int_t bkgFactor = 20;
double mass_sig_min = 0.486;
double mass_sig_max = 0.509;
double mass_leftside_low = 0.435;
double mass_leftside_high = 0.480;
double mass_rightside_low = 0.515;
double mass_rightside_high = 0.56;
double mis_la_range = 0.010;
double decaylength = 2.5;
double costheta = 0.999;
/*
double mass_sig_min = 0.483;
double mass_sig_max = 0.512;
double mass_leftside_low = 0.435;
double mass_leftside_high = 0.480;
double mass_rightside_low = 0.515;
double mass_rightside_high = 0.56;
*/
//--------------------------------------
