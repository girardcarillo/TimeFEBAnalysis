// Author Cloé Girard-Carillo girardcarillo@lal.in2p3.fr

#include <limits>
#include <string>
#include <iostream>
#include <fstream>

#include <TTree.h>
#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>

#include "/home/girardcarillo/Workdir/SNPlot/RootDisplay.h"
#include "/home/girardcarillo/Workdir/SNPlot/EventTree.h"

// root 'RootAnalyzer_syncFEB.cc("Cut_2e_DeltaT.root",6,0,1)' -b -q


using namespace std ;

// typedef numeric_limits<double> dbl ;

const int row_tot_number=16 ;
const int column_tot_number=20 ;
const double DEFAULT_CALO_SAMPLING_PERIOD = 400/1024 ;
const int CALO_MAX_NB_OF_SAMPLES = 1024 ;

void fit_delta_t(TH1F* histo, double *fit_parameters) ;
bool select_PM(vector<int> *vect_column,vector<int> *vect_row,int selected_column,int selected_row,int *hit) ;
double distance_OM(int col1,int row1,int col2,int row2) ;
// void Correct_CaloTimeINL(float inputDataArray[], int nbOFSamples,int slotIndex, int channel, int fcr, float *correctedDataArray) ;

void RootAnalyzer_syncFEB(string filename,int col_ref,int row_ref, bool enable_drawing = 0){

  TFile *theInFile = new TFile(filename.c_str(),"READ") ;
  TTree *theTree = nullptr ;
  theInFile->GetObject("T",theTree) ;

  if (theInFile->IsOpen()) {
    cout << "File " << filename << " opened sucessfully" << endl ;
  }

  theTree = (TTree*)theInFile->Get("DataCut") ;

  theTree->SetBranchAddress("trigger_id",&trigger_id) ;
  theTree->SetBranchAddress("calo_row",&calo_row) ;
  theTree->SetBranchAddress("calo_column",&calo_column) ;
  theTree->SetBranchAddress("calo_id",&calo_id) ;
  theTree->SetBranchAddress("calo_module",&calo_module) ;
  theTree->SetBranchAddress("channel_raw_tdc", &channel_raw_tdc) ;
  theTree->SetBranchAddress("calo_time",&calo_time) ;
  theTree->SetBranchAddress("calo_energy",&calo_energy) ;
  theTree->SetBranchAddress("cut_multiplicity", &cut_multiplicity) ;
  theTree->SetBranchAddress("calo_peak", &calo_peak) ;
  theTree->SetBranchAddress("channel_raw_tdc", &channel_raw_tdc) ;
  theTree->SetBranchAddress("calo_number", &calo_number) ;
  // theTree->SetBranchAddress("waveform_fcr", &waveform_fcr) ;

  cout << theTree->GetEntries() << " entries"<< endl ;

  TH2D *h2mean = new TH2D ("mean","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
  TH2D *h2meanError = new TH2D ("mean_error","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
  TH2D *h2sigma = new TH2D ("sigma","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
  TH2D *h2sigmaError = new TH2D ("sigma_error","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
  TH2D *h2counts = new TH2D ("counts","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
  TH2D *h2timeEnergy = new TH2D ("Delta t vs energy","", 100, 0.03, 0.17, 100, 3.9, 4.2) ;
  TH2D *h2timePeak = new TH2D ("Delta t vs amplitude","", 100, 0.03, 0.17, 100, -425, -380) ;
  TH2D *h2timet0 = new TH2D ("Delta t vs t0","", 500, 94, 101, 500, -0.22, -0.06) ;
  TH2D *h2t0Event = new TH2D ("t0 vs event number","", 5000, 0, 5000, 500, -0.22, -0.06) ;
  TH2D *h2t0t1 = new TH2D ("t0 vs t1","", 500, 94, 101, 500, 94, 101) ;

  TProfile *hmean_distance = new TProfile ("mean vs distance","", 100, 0, 10, -2, 2) ;

  TH1F *hmean = new TH1F("mean","",50, -2, 2) ;
  TH1F *hsigma = new TH1F("sigma","",50, 0 , 0.1) ;

  TH1F *htime[column_tot_number][row_tot_number] ;
  for (int i=0 ;i<column_tot_number ;i++) {
    for (int j=0 ;j<row_tot_number ;j++) {
      htime[i][j] = new TH1F(Form("[%d:%d]&[%d,%d]",col_ref,row_ref,i,j),Form("[%d:%d]&[%d,%d]",col_ref,row_ref,i,j),500,-2,2) ;
    }
  }


  for (Long64_t i=0 ;i<theTree->GetEntries() ;i++) {
    theTree->GetEntry(i) ;


    if (i%100000==0) cout << "event " << i << endl ;
    if (isnan(calo_energy->at(0)) || isnan(calo_energy->at(1))){
      continue ;
    }

    else {

      //tableau coïncidences en temps
      int hit0 = -1 ;
      bool selected_PMT = select_PM(calo_column,calo_row,col_ref,row_ref,&hit0) ;
      int hit1 = abs(hit0-1);
      double Delta_t = 0. ;

      if (selected_PMT) {
        Delta_t = calo_time->at(hit0)-calo_time->at(hit1) ;
        htime[calo_column->at(hit1)][calo_row->at(hit1)]->Fill(Delta_t) ;

        if (calo_column->at(hit0) == col_ref && calo_row->at(hit1) == 1) {
          h2t0Event->Fill(i,Delta_t) ;
          h2t0t1->Fill(calo_time->at(hit0)-channel_raw_tdc->at(hit0)*6.25,calo_time->at(hit1)-channel_raw_tdc->at(hit1)*6.25) ;
          h2timeEnergy->Fill(Delta_t,calo_energy->at(hit0)+calo_energy->at(hit1)) ;
          h2timePeak->Fill(Delta_t,calo_peak->at(hit0)) ;
          h2timet0->Fill(calo_time->at(hit0)-channel_raw_tdc->at(hit0)*6.25,Delta_t) ;
        }

      }
    }
  }


  int max_counts = -1 ;
  int counter_OM = 0 ;
  for (int i = 0 ; i < column_tot_number ; ++i) {
    for (int j = 0 ; j< row_tot_number ; ++j) {
  // for (int i = 0 ; i < 1 ; ++i) {
  //   for (int j = 0 ; j< row_tot_number ; ++j) {
      if (htime[i][j]->GetEntries() > max_counts) {
        max_counts = htime[i][j]->GetEntries() ;
      }
      h2counts->SetBinContent(i+2,j+2,htime[i][j]->GetEntries()) ;
      if (htime[i][j]->GetEntries()) {
        double fit_parameters[4] ;
        fit_delta_t(htime[i][j],fit_parameters) ;
        h2mean->SetBinContent(i+2,j+2,fit_parameters[0]) ;
        h2meanError->SetBinContent(i+2,j+2,fit_parameters[1]) ;
        h2sigma->SetBinContent(i+2,j+2,fit_parameters[2]) ;
        h2sigmaError->SetBinContent(i+2,j+2,fit_parameters[3]) ;
        hmean_distance->Fill(distance_OM(9,8,i,j),fit_parameters[0]) ;

        hmean->Fill(fit_parameters[0]) ;
        hsigma->Fill(fit_parameters[2]) ;
      }
      counter_OM++ ;
    }
  }


  // ///Drawing

  if (enable_drawing) {

    TCanvas *c2 = new TCanvas("c2","c2",10,10,2000,1000) ;

    int counter_color = 0 ;
    config_histo1D(htime[col_ref][1],"","Delta t (ns)","#counts",1,1,MultiPlotColors(counter_color)) ;
    htime[col_ref][1]->GetYaxis()->SetRangeUser(0,max_counts/50.) ;
    for (int i = 0 ; i < column_tot_number ; ++i) {
      for (int j = 0 ; j< row_tot_number ; ++j) {
        if (i != col_ref || j != 1) {
          if (htime[i][j]->GetEntries()) {

            counter_color++ ;
            config_histo1D(htime[i][j],"SAME","Delta t (ns)","#counts",1,1,MultiPlotColors(counter_color)) ;

          }
        }
      }
    }

    gStyle->SetLegendBorderSize(0) ;
    c2->BuildLegend(0.76,0.35,0.89,0.89) ;
    htime[col_ref][1]->SetTitle("Delta t for different channels") ;
    htime[col_ref][1]->GetYaxis()->SetRangeUser(0,max_counts/30.) ;
    gStyle->SetOptStat(0) ;
    c2->SaveAs("plots/delta_t.pdf") ;

    // TCanvas *c3 = new TCanvas("c3","c3",10,10,2000,1000) ;
    // c3->Divide(2,2) ;

    // c3->cd(1) ;
    // config_histo1D(htime[col_ref][2],"SAME","Delta t (ns)","#counts",1,1,4) ;
    // cout << htime[col_ref][2]->GetXaxis()->GetXmin() << endl ;
    // cout << htime[col_ref][2]->GetXaxis()->GetXmax() << endl ;
    // cout << htime[col_ref][2]->GetMinimum() << endl ;
    // cout << htime[col_ref][2]->GetMaximum() << endl ;

    // c3->cd(2) ;
    // config_histo2D(h2t0Event,"", "event number","Delta t","BOX") ;

    // c3->cd(3) ;
    // config_histo2D(h2timet0,"Delta t vs t0", "Delta t","t0","BOX") ;

    // c3->cd(4) ;
    // config_histo2D(h2t0t1,"t0 vs t1", "t0","t1","BOX") ;
    // c3->SaveAs("plots/4_plots.pdf") ;


    TCanvas *c1 = new TCanvas("c1","c1",10,10,2000,1000) ;
    gStyle->SetPaintTextFormat("1.f") ;
    config_histo2D(h2counts,"Total counts for each channel", "Slot","Channel","COLZTEXT") ; c1->SaveAs("plots/counts.pdf") ;

    gStyle->SetPaintTextFormat("1.4f") ;
    config_histo2D(h2mean,"Mean of fitted Delta t ditribution for each channel", "Slot","Channel","COLZTEXT") ; c1->SaveAs("plots/mean.pdf") ;
    config_histo2D(h2sigma,"Sigma of fitted Delta t ditribution for each channel", "Slot","Channel","COLZTEXT") ; c1->SaveAs("plots/sigma.pdf") ;
    config_histo2D(h2sigmaError,"Sigma error of fitted Delta t ditribution for each channel", "Slot","Channel","COLZTEXT") ; c1->SaveAs("plots/sigmaError.pdf") ;
    config_histo2D(h2meanError,"Mean error of fitted Delta t ditribution for each channel", "Slot","Channel","COLZTEXT") ; c1->SaveAs("plots/meanError.pdf") ;
    config_histo2D(h2timeEnergy,"Delta t vs charge", "Delta t","charge","COLZ") ; c1->SaveAs("plots/delta_t_charge.pdf") ;
    config_histo2D(h2timePeak,"Delta t vs amplitude", "Delta t","amplitude (mV)","COLZ") ; c1->SaveAs("plots/delta_t_ampl.pdf") ;

    TCanvas *c4 = new TCanvas("c4","c4",10,10,2000,1000) ;
    c4->Divide(1,2) ;
    c4->cd(1) ;
    config_histo1D(hmean,"","Mean of Delta t distribution (ns)","",2,1,1) ;
    c4->cd(2) ;
    config_histo1D(hsigma,"","Sigma of Delta t distribution (ns)","",2,1,1) ;
    c4->SaveAs("plots/sigma_and_mean.pdf") ;

  }

  gStyle->SetOptStat(0);
  hmean_distance->Draw() ; hmean_distance->GetXaxis()->SetTitle("Distance with CB") ; hmean_distance->GetYaxis()->SetTitle("Mean of delta t distributions") ;

  theTree->ResetBranchAddresses() ;
}












bool select_PM(vector<int> *vect_column,vector<int> *vect_row,int selected_column,int selected_row,int *hit){
  bool flag_test=0 ;
  for (int i = 0 ; i < vect_row->size() ; ++i) {
    if (vect_column->at(i)==selected_column) {
      if (vect_row->at(i)==selected_row) {
        *hit=i ;
        flag_test=1 ;
      }
    }
  }
  return flag_test ;
}

void fit_delta_t(TH1F* histo, double fit_parameters[4]){
  meanhisto = histo->GetMean() ;
  sigmahisto = histo->GetRMS() ;
  TF1 *f1 = new TF1("f1","gaus",(meanhisto-0.1),(meanhisto+0.1)) ;
  f1->SetParameter(1,meanhisto) ;
  f1->SetParameter(2,sigmahisto) ;
  histo->Fit("f1","RQ") ;
  fit_parameters[0] = f1->GetParameter(1) ;
  fit_parameters[1] = f1->GetParError(1) ;
  fit_parameters[2] = f1->GetParameter(2) ;
  fit_parameters[3] = f1->GetParError(2) ;
  delete f1 ;
}


double distance_OM(int col1,int row1,int col2,int row2){

  return sqrt(pow((col2-col1),2)+pow((row2-row1),2)) ;

}


// void Correct_CaloTimeINL(float inputDataArray[], int slotIndex, int channel, int fcr, float *correctedDataArray){
//   float prevCoef, currentCoef, nextCoef, x, x1, x2, x3, y1, y2, y3;
//   int i, n, prev_i, next_i, current_i, prev_n, next_n, current_n;
//   float samplingPeriod;
//   int nbOFSamples = 1024;

//   samplingPeriod = DEFAULT_CALO_SAMPLING_PERIOD*1000.0;   // en picosecond

//   for(n=0; n<nbOFSamples; n++){  // On corrige toutes les cellules
//     // Interpolation polynomiale de Lagrange

//     i = (n + fcr) % CALO_MAX_NB_OF_SAMPLES - 1; //  On se recale dans le tableau d'INL qui est décalé de 1

//     if(i < 0) i += CALO_MAX_NB_OF_SAMPLES;

//     if((n != 0 ) && (n!= nbOFSamples -1)){
//       prev_i = (i - 1 + CALO_MAX_NB_OF_SAMPLES)%CALO_MAX_NB_OF_SAMPLES ;
//       current_i = i;
//       next_i = (i+1) % CALO_MAX_NB_OF_SAMPLES;

//       prev_n = (n - 1 + CALO_MAX_NB_OF_SAMPLES)%CALO_MAX_NB_OF_SAMPLES ;
//       current_n = n;
//       next_n = (n+1) % CALO_MAX_NB_OF_SAMPLES;

//       x = current_n;
//     }
//     else if (n == 0){
//       prev_i = i;
//       current_i = (i+1)% CALO_MAX_NB_OF_SAMPLES;
//       next_i = (i+2)% CALO_MAX_NB_OF_SAMPLES;

//       prev_n = n;
//       current_n = (n+1)% CALO_MAX_NB_OF_SAMPLES;
//       next_n = (n+2)% CALO_MAX_NB_OF_SAMPLES;

//       x = prev_n;
//     }
//     else{ // n = Nb_Of_Samples -1
//       prev_i = (i-2 + CALO_MAX_NB_OF_SAMPLES)%CALO_MAX_NB_OF_SAMPLES ;
//       current_i = (i-1 + CALO_MAX_NB_OF_SAMPLES)% CALO_MAX_NB_OF_SAMPLES;
//       next_i = i;

//       prev_n = (n-2 + CALO_MAX_NB_OF_SAMPLES)%CALO_MAX_NB_OF_SAMPLES ;
//       current_n = (n-1 + CALO_MAX_NB_OF_SAMPLES)% CALO_MAX_NB_OF_SAMPLES;
//       next_n = n;

//       x = next_n;
//     }

//     prevCoef =     CaloCalibParams[slotIndex].TimeINLPedestals[channel][prev_i] / samplingPeriod;
//     currentCoef = CaloCalibParams[slotIndex].TimeINLPedestals[channel][current_i] / samplingPeriod;
//     nextCoef = CaloCalibParams[slotIndex].TimeINLPedestals[channel][next_i] /samplingPeriod;

//     x1 = prev_n + prevCoef;
//     x2 = current_n + currentCoef;
//     x3 = next_n + nextCoef;

//     y1 = inputDataArray[prev_n];
//     y2 = inputDataArray[current_n];
//     y3 = inputDataArray[next_n];

//     correctedDataArray[n] = (y1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3)))+
//       (y2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))) +
//       (y3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2)));
//   }
// }



// TGraph *gsigmaError = new TGraph(counter_OM,stat,sigmaError) ;
// gsigmaError->SetTitle("Error on sigma vs stat") ;
// gsigmaError->GetXaxis()->SetTitle("stat") ;
// gsigmaError->GetYaxis()->SetTitle("Error on sigma (ns)") ;
// gsigmaError->Draw("A*") ;
// gPad->Print("plots/statsigmaError.pdf") ;
