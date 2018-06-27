/*
The Goal of this plotting script is to consider the 1D case in each parameter (except for likelihood) and find the minimum likelihood in each iteration for that variable while varying the other parameters.
To do this, we make a series of histogram arrays in each of the four parameters where the histogram arrays are filled with the appropriate likelihood.
In particular, we will look at an ntuple value, and look at the ith parameter, then go to the appropriate hist array, find the corresponding histogram and fill that hist with the loglikelihood.
This will be significantly faster than iterating through the ntuple 4 times to distribute the appropriate loglikelihood values. Moreover, getting the minimum will be a simple ->GetMinimum() from each hist, which will then be plotted ina new histogram.
Requires that we know the precision for the spatial variables and the time variable and know where they start and end so that the histograms can be made apriori. Accordingly, will set constant variables in the beginning to fix that.
*/
{
#include <vector>     

  const int ntuplecount = 1;

  TObjArray xtemphists,ytemphists,ztemphists,ttemphists;    // Histogram array that will hold the hists that will be filled when each loglikelihood is distributed.
  
  /*TH1F* xhist;
  TH1F* yhist;
  TH1F* zhist;
  TH1F* thist;*/

  TH1F* xhist_final;
  TH1F* yhist_final;
  TH1F* zhist_final;
  TH1F* thist_final;

  Int_t check,count,mincount,spatial_steps,time_steps,count2,x_index,y_index,z_index,t_index,ntuple_entries;
  Float_t xpos,ypos,zpos,time,loglikelihood;

  // Read in the file and make an output file
  TFile TheFile("./fiTQun_015ns_delay.root");
  TFile NewFile("./Likelihood_Plots.root","RECREATE");

  // The Ntuple to be read into
  TNtuple* Ntuple; 

  // Read the First ntuple to get an idea of how large the ntuples are
  TheFile.cd();
  Ntuple = (TNtuple*)TheFile.Get("likelihoodprefit_Event000");

  // Values set apriori according to the fiTQun fitting parameters.
  const Int_t inf = 1000000000;
  const Int_t time_precision = 1;
  const Int_t time_min = 840;
  const Int_t time_max = 1000;
  const Int_t spatial_min = -500;
  const Int_t spatial_max = 500;
  const Int_t spatial_precision = 25;

  // Figure out which index has the most indices to run through.
  time_steps = (time_max - time_min)/time_precision + 1;
  spatial_steps = (spatial_max - spatial_min)/spatial_precision + 1;

  if (time_steps > spatial_steps) {
    mincount = time_steps;
    count2 = spatial_steps;
    cout << "More Time steps than spatial steps!" << endl;
  } else { 
    mincount = spatial_steps;
    count2 = time_steps;
    cout << "More spatial steps than time steps!" << endl;
  }

  check = 0;
  count = 0;
  // Allocate the matrices to store the data
  vector<Float_t> time_mat(time_steps);
  vector<Float_t> xpos_mat(spatial_steps);
  vector<Float_t> ypos_mat(spatial_steps);
  vector<Float_t> zpos_mat(spatial_steps);
  
  for (int m = 0; m < time_steps; m++) {
    switch(

  /*vector <Float_t> temptime;
  vector <Float_t> tempxpos;
  vector <Float_t> tempypos;
  vector <Float_t> tempzpos;*/
  
  // Loop through the ntuples
  for (int i = 0; i < ntuplecount; i++) {
    TheFile->cd();
    cout << "Working on Ntuple number " << i << endl;
    Ntuple = (TNtuple*)TheFile.Get(Form("likelihoodprefit_Event%03d",i));    // Read from the next ntuple
    Ntuple->SetBranchAddress("xpos",&xpos);
    Ntuple->SetBranchAddress("ypos",&ypos);
    Ntuple->SetBranchAddress("zpos",&zpos);
    Ntuple->SetBranchAddress("time",&time);
    Ntuple->SetBranchAddress("loglikelihood",&loglikelihood);
    xhist_final = new TH1F(Form("loglikelihood_vs_xpos_Event%03d",i), "Loglikelihood vs xpos", spatial_steps, spatial_min, spatial_max);
    yhist_final = new TH1F(Form("loglikelihood_vs_ypos_Event%03d",i), "Loglikelihood vs ypos", spatial_steps, spatial_min, spatial_max);
    zhist_final = new TH1F(Form("loglikelihood_vs_zpos_Event%03d",i), "Loglikelihood vs zpos", spatial_steps, spatial_min, spatial_max);
    thist_final = new TH1F(Form("loglikelihood_vs_time_Event%03d",i), "Loglikelihood vs time", time_steps, time_min, time_max);
    ntuple_entries = Ntuple->GetEntries();
    check = 0;
    for (int k = 0; k < ntuple_entries; k++) {
      Ntuple->GetEntry(k);
      if (loglikelihood != loglikelihood || loglikelihood > inf) {loglikelihood = inf;} // Set large/nan loglikelihoods to inf
      cout << "Loglikelihood at " << k << " is " << loglikelihood << endl;
      // Fill the appropriate vectors with the loglikelihood if it is smaller than the prior one
      x_index = ((int)xpos - spatial_min)/spatial_precision;
      y_index = ((int)ypos - spatial_min)/spatial_precision;
      z_index = ((int)zpos - spatial_min)/spatial_precision;
      t_index = ((int)time - time_min)/time_precision;
      if (xpos
	  xpos_mat[x_index] = loglikelihood;
	  ypos_mat[y_index] = loglikelihood;
	  zpos_mat[z_index] = loglikelihood;
	  time_mat[t_index] = loglikelihood;
    }
    // Write the minimum of each column vector to the appropriate 
    check = 0;
    for (int j = 0; j < time_steps; j++) {
      switch(check) {
      case 0:
	xhist_final->SetBinContent(j,xpos_mat[j].
    NewFile.cd();
    cout << "Saving plots ... " << endl;
    // Save the plots that we have made
    NewFile.cd();
    xhist_final->Write();
    yhist_final->Write();
    zhist_final->Write();
    thist_final->Write();
  }
}

