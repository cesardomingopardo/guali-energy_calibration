//gROOT->Reset();

// c/c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <time.h>

// root specific headers
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

const int npeaks = 10;

Double_t fpeaks(Double_t *x, Double_t *par);
const std::string currentDateTime();
int write_peakfits_to_file(string cal_file, string file_peakfit_pars);
void make_calibration(const string inputfilename, const string outputfilename, const int npoints);

void e_calibration(void)
{
   
  int n_cal_files = 2;
  const int max_cal_files = 5;
  string cal_files[max_cal_files];
  
  cal_files[0] = "./Ecal_files/137Cs_00_00_1200s_45cm.root";
  cal_files[1] = "./Ecal_files/60Co_00_00_300s_3cm.root";

  // cal_files[0] = "./root_files_ecal_old/137Cs_Ecal.root";
  // cal_files[1] = "./root_files_ecal_old/60Co_Ecal.root";
  
  string file_cal_pars = "calibrations";
  string file_peakfit_pars = "peakfits_";
  const string date_time_string = currentDateTime();

  file_peakfit_pars += date_time_string + ".txt";
  string cal_dir = "calibrations/" + date_time_string;
  gSystem->Exec(Form("mkdir ./%s",cal_dir.c_str()));
  file_cal_pars = cal_dir + "/calibration.txt";
  file_peakfit_pars = cal_dir + "/peakfits.txt";

 

  cout << " peakfits file = " << file_peakfit_pars << endl;
  cout << " calibration file = " << file_cal_pars << endl;

  int ndata = 0;
  for(int i=0;i<n_cal_files;i++)
    {
      cout << " Analysing Calibration File: " << cal_files[i] << endl;
      ndata += write_peakfits_to_file(cal_files[i], file_peakfit_pars);
    }

  make_calibration(file_peakfit_pars,file_cal_pars,ndata); // linear energy fit --> Energy calibration parameters

  //  gSystem->Exec(Form("gv %s.137Cs.ps &",file_peakfit_pars.c_str()));
  gSystem->Exec(Form("inkscape %s.137Cs.ps &",file_peakfit_pars.c_str()));

  //  gSystem->Exec(Form("gv %s.60Co.ps &",file_peakfit_pars.c_str()));
  gSystem->Exec(Form("inkscape %s.60Co.ps &",file_peakfit_pars.c_str()));
  
  return;
}

int write_peakfits_to_file(string cal_file, string file_peakfit_pars) // returns 1 for Cs, 2 for Co
{
  bool cs_flag = 0;
  bool co_flag = 0;
  int ndata = 0;
  
  // Checking if root file exists
  TFile *f = new TFile(cal_file.c_str());
  if(!(f->IsOpen()))
    {
      cout << " PROBLEM OPENING FILE: " << cal_file << " (Check if file exists) " << endl;
      exit(0);
    }
 

  // Checking for the expected calibration source, only used to return one or two peaks
  if(strstr(cal_file.c_str(),"Cs"))
    cs_flag = 1;
  else if(strstr(cal_file.c_str(),"Co"))
    co_flag = 1;
  else
    {
      cout << " Check consistency of calibration file names" << endl;
      cout << " the filename must contain either Cs or Co " << endl;
      exit(0);
    }
    
  //  f->ls();
  TTree *t1 = (TTree*)f->Get("T");
  //  t1->Print();
  Float_t E;
  Long64_t T;
  
  t1->SetBranchAddress("E",&E);
  t1->SetBranchAddress("T",&T);
  
  Int_t nentries = (Int_t) t1->GetEntries();

  const float bin_units = 20;
  const int nbins = 1000;
  const float hrange = bin_units*nbins;

  // TCanvas *cT = new TCanvas("cT","",1);
  // cT->SetFillColor(10);
  // t1->Draw("T>>hT");
  t1->GetEntry(0);
  Long64_t t_0 = T;
  t1->GetEntry(nentries-1);
  Long64_t t_f = T;
  Long64_t Delta_t = (t_f - t_0)*4E-9; // s
  cout << "          Measurement time = " << Delta_t << " (s) " << endl;
    
  TCanvas *cE = new TCanvas("cE","",1);
  cE->SetFillColor(10);
  // t1->Draw(Form("E>>hE(%d,0,%d)",nbins,nbins*bin_units));
  TH1F *hE = new TH1F("hE","Energy Histogram",nbins,0,hrange);
  t1->Draw("E>>hE");

  //  hE->GetXaxis()->SetRange(0,nbins);
  //  hE->GetXaxis()->SetRangeUser(0,nbins);

  const int reb_factor = 10;
  
  hE->Rebin(reb_factor);

  // doing peak search:
  //-----------  parameters for peak search --------------------------------
  float default_sigma = 50/reb_factor; // approx. sigma peaks
  float threshold = 0.1;// amplitude_smallest_peak/amplitude_highest_peak
  float resolution = 1;// 1:= default value, higher values higher resolution (but not used)

  TSpectrum *my_spectrum = new TSpectrum(npeaks,resolution);
  my_spectrum->Search(hE,0.5*default_sigma,"new",threshold);
  Int_t n_found = my_spectrum->GetNPeaks();
  // Trying to get positions from polymarkers, because  GetPositionX() fails like mad
  //  //  Double_t *xpeaks = my_spectrum->GetPositionX();
  TList *functions = hE->GetListOfFunctions();
  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
  Double_t *xpeaks = pm->GetX();
  
  for(int i=0;i<n_found;i++)
    cout << "****************     DEBUG             peak " << i << " at " << xpeaks[i] << endl;
  // no ordered data:  first peak:
  
  // identifying the last two peaks in the spectrum
  const int max_cal_peaks = 2;

  Double_t lastpeak[max_cal_peaks] = xpeaks[0];
  
  for(int j=0;j<n_found;j++)
    if(xpeaks[j]>lastpeak[0])
      lastpeak[0] = xpeaks[j]; // lastpeak[0] = very last peak
  
  Float_t delta_x[npeaks]; // distance with respect to last one
  Float_t min_delta_x = 1E8;
  Int_t index_second_last = 0;
  
  for(int j=0;j<n_found;j++)
    {
      delta_x[j] = lastpeak[0] - xpeaks[j];
      if(min_delta_x > delta_x[j] && delta_x[j] > 0.0)
	{
	  min_delta_x = delta_x[j];
	  lastpeak[1] = xpeaks[j];
	}
    }

  if(co_flag)
    {
      n_found = 2; // taking last two
      cout << " 60Co peaks found at " << endl;
      cout << " peak 1 = " << lastpeak[1]<< endl;
      cout << " peak 2 = " << lastpeak[0]<< endl;
    }
  else if(cs_flag)
    {
      n_found = 1; // taking last one
      cout << " 137Cs peak found at " << endl;
      cout << " peak 1 = " << lastpeak[0]<< endl;
    }

  
  printf("Found %d candidate peaks to fit\n",n_found);
  for(int i=0;i<n_found;i++)
    cout << " peak# " << i << " at " << lastpeak[i] << endl;

    
  //-------------------------------------------------------------------------
  // REPLACE THIS BY BACKGROUND MEASUREMENT, NORMALIZATION AND SUBTRACTION
  //
  // doing background estimation (https://root.cern.ch/doc/v608/classTSpectrum.html)
  //-------------------------------------------------------------------------
  // //  hE->GetXaxis()->SetRange(firstpeak/(bin_units*reb_factor),nbins*bin_units/(bin_units*reb_factor));
  // TH1F* hbkg = my_spectrum->Background(hE,20,"");
  // hbkg->SetLineColor(2);
  // hbkg->Scale(hbkg->Integral()/hE->Integral());
  // hbkg->Draw("same");
  // TH1F* hE_bs = (TH1F*) hE->Clone();
  // hE_bs->Reset();
  // hE_bs->Add(hE,hbkg,1,-1);
  // hE_bs->SetLineColor(1);
  // hE_bs->Draw("same");
  //-----------------------------------------

  //-------------------------------------------
  // Fitting the peaks:
  //-------------------------------------------
  TF1 *func_peaks[npeaks];

  for(int i=0;i<n_found;i++) // 0:= last one; 1 := previous last one
    {
      
      float peak_pos = lastpeak[i];
      float peak_amp = hE->GetBinContent(hE->GetXaxis()->FindBin(peak_pos));
      float llimit = lastpeak[i]-default_sigma*reb_factor*bin_units;
      float ulimit = lastpeak[i]+default_sigma*reb_factor*bin_units;

      cout << " peak# = " << i << endl;
      cout << "  ll = " << llimit << endl;
      cout << "  ul = " << ulimit << endl;
      cout << "  peak_pos = "<< peak_pos << endl;
      cout << "  peak_amp = " << peak_amp << endl;
      
      func_peaks[i] = new TF1(Form("fit_func[%d]",i),fpeaks,llimit,ulimit,5);
      func_peaks[i]->SetParameter(2,peak_amp);
      func_peaks[i]->SetParameter(3,peak_pos);
      func_peaks[i]->SetParLimits(3,llimit,ulimit);
      
      func_peaks[i]->SetParameter(4,default_sigma);
      func_peaks[i]->SetLineColor(2);
      
      hE->Fit(Form("fit_func[%d]",i),"","same",llimit,ulimit);
    }
  
  for(int i=0;i<n_found;i++)
    func_peaks[i]->Draw("same");


  
  ofstream fout(file_peakfit_pars.c_str(),ios_base::app);

  // I store both numerical value of peak_search algorithm (which is quite reliable) and the fit-value (which is less reliable but may be more accurate)

  cout << "  Source         E(keV)            numerical_pos          fit_pos          Err_fit_pos   " << endl;
  if(cs_flag)
    {
      fout << "137Cs         662               " << lastpeak[0] << "  " << func_peaks[0]->GetParameter(3) << "   " << func_peaks[0]->GetParError(3) << endl;
       cout << "137Cs         662               " << lastpeak[0] << "  " << func_peaks[0]->GetParameter(3) << "   " << func_peaks[0]->GetParError(3) << endl;
      string psfilename = file_peakfit_pars + ".137Cs.png";
      cE->Print(psfilename.c_str());
      psfilename = file_peakfit_pars + ".137Cs.ps";
      cE->Print(psfilename.c_str());
      ndata = 1;
    }

  if(co_flag)
    {
      fout << "60Co         1173               " << lastpeak[1] << "  "<< func_peaks[1]->GetParameter(3) << "  " << func_peaks[1]->GetParError(3)<<  endl;
      cout << "60Co         1173               " << lastpeak[1] << "  "<< func_peaks[1]->GetParameter(3) << "  " << func_peaks[1]->GetParError(3)<<  endl;
      fout << "60Co         1332               " << lastpeak[0] << "  "<< func_peaks[0]->GetParameter(3) << "  " << func_peaks[0]->GetParError(3)<< endl;
      cout << "60Co         1332               " << lastpeak[0] << "  "<< func_peaks[0]->GetParameter(3) << "  " << func_peaks[0]->GetParError(3)<< endl;
      string psfilename = file_peakfit_pars + ".60Co.png";
      cE->Print(psfilename.c_str());
      psfilename = file_peakfit_pars + ".60Co.ps";
      cE->Print(psfilename.c_str());
      ndata = 2;
    }

  cE=0;
  delete cE;
  my_spectrum=0;
  delete my_spectrum;
  //  delete xpeaks; // this fucks up everything, do not understand why yet...
  func_peaks[] = 0;
  delete [] func_peaks;
  
  f->Close();
  
  return ndata;
}

Double_t fpeaks(Double_t *x, Double_t *par) { // linear background + gaussian peak
   Double_t result = par[0] + par[1]*x[0];

   Double_t norm  = par[2];
   Double_t mean  = par[3];
   Double_t sigma = par[4];
   
   result += norm*TMath::Gaus(x[0],mean,sigma);

   return result;
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "d%dm%my%Y_%OHh%Mm%Ss", &tstruct);

    return buf;
}

void make_calibration(const string inputfilename, const string outputfilename, const int npoints)
{
  bool bad_peak_fit_flag = 0;
  const int max_ndata = 5;
  Double_t Eg[max_ndata] = {0};
  Double_t err_Eg[max_ndata] = {0};
  Double_t peak_x[max_ndata] = {0};
  Double_t peak_x_num[max_ndata] = {0};
  Double_t err_peak_x_num[max_ndata] = {0};
  Double_t err_peak_x[max_ndata] = {0};
  string dummy;
  
  ifstream finp(inputfilename.c_str());
  if(!finp.good())
    {
      cout << " PROBLEMS OPENING FILE WITH PEAK-FIT RESULTS " << endl;
      cout << "   filename = " << inputfilename << endl;
      exit(0);
    }
  int i=0;
  while(!finp.eof())
    {
      finp >> dummy >> Eg[i] >> peak_x_num[i] >> peak_x[i] >> err_peak_x[i];
      err_Eg[i] = 0;
      err_peak_x_num[i] = 0.01*peak_x_num[i];
      i++;
    }

  
  int np = i-1;
  cout << " " << endl;
  cout << " " << endl;
  cout << np << " calibration points read from " << inputfilename << endl;
  cout << " " << endl;
  cout << " " << endl;


  for(int i=1;i<np;i++)
    {
      if(peak_x[i] < peak_x[i-1])
	bad_peak_fit_flag = 1;
    }   

  if(bad_peak_fit_flag)
    {
      cout << "                                                                   " << endl;
      cout << "                                                                   " << endl;
      cout << " Peak-fits are not in increasing values, use numerical values fit  " << endl;
      cout << "                                                                   " << endl;
      cout << "                                                                   " << endl;
    }
  
  TGraphErrors *gr = new TGraphErrors(np,peak_x,Eg,err_peak_x,err_Eg);
  gr->SetTitle("Calibration using Gaussian peak fit");
  gr->GetXaxis()->SetTitle("Channels");
  gr->GetYaxis()->SetTitle("Energy (keV)");

  TGraphErrors *grn = new TGraphErrors(np,peak_x_num,Eg,err_peak_x_num,err_Eg);
  grn->SetTitle("Calibration using numerical peak search algorithm");
  grn->GetXaxis()->SetTitle("Channels");
  grn->GetYaxis()->SetTitle("Energy (keV)");

        

  
  TCanvas *cCal = new TCanvas("cCal","Energy Calibration", 1);
  cCal->SetFillColor(10);
  cCal->Divide(1,2);
  cCal->cd(1);
  // calibration with fitted peaks
  gr->SetMarkerStyle(21);
  gr->Draw("ALP");
  TF1 *f1 = new TF1("f1","pol1",peak_x[0],peak_x[np-1]);
  // initial parameters estimate:
  float slope_i = (Eg[1] - Eg[0])/(peak_x[1] - peak_x[0]);
  float offset_i = Eg[0] - slope_i*peak_x[0];
  f1->SetParameters(offset_i,slope_i);
  TFitResultPtr fit_results = gr->Fit("f1","RS");

  // calibration with numerical peak values
  cCal->cd(2);
  grn->SetMarkerStyle(22);
  grn->SetMarkerColor(2);
  grn->Draw("ALP");
  TF1 *f1n = new TF1("f1n","pol1",peak_x_num[0],peak_x_num[np-1]);
  // initial parameters estimate:
  slope_i = (Eg[1] - Eg[0])/(peak_x_num[1] - peak_x_num[0]);
  offset_i = Eg[0] - slope_i*peak_x_num[0];
  f1n->SetParameters(offset_i,slope_i);
  TFitResultPtr fit_results_n = grn->Fit("f1n","RS");

  
  string fig_file_name = outputfilename + ".png";
  cCal->Print(fig_file_name.c_str());
  fig_file_name = outputfilename + ".ps";
  cCal->Print(fig_file_name.c_str());
  
  ofstream fout(outputfilename.c_str(),ios_base::app);
  fout<< f1n->GetParameter(0) << "        " << f1n->GetParameter(1) << endl;
  fout << " =========================================" << endl;
  fout << " above:   num_offset  num_slope  " << endl;
  fout << " below more details:  " << endl;
  fout << " =========================================" << endl;
  fout << " offset_keV  " << f1->GetParameter(0) << "  " << f1->GetParError(0)<< endl;
  fout << " slope_keV_per_chan  " << f1->GetParameter(1) << "  " << f1->GetParError(1)<< endl;
  fout << " Chi2 = " << fit_results->Chi2() << endl;
  fout << " Below Numerical Values" << endl;
  fout << " offset_keV  " << f1n->GetParameter(0) << "  " << f1n->GetParError(0)<< endl;
  fout << " slope_keV_per_chan  " << f1n->GetParameter(1) << "  " << f1n->GetParError(1)<< endl;
  fout << " Chi2 = " << fit_results_n->Chi2() << endl;
  
  fout << "  " << endl;

  fout.close();
}

