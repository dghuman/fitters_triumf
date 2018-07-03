/*
The Goal of this plotting script is to consider the 1D case in each parameter (except for likelihood) and find the minimum likelihood in each iteration for that variable while varying the other parameters.
To do this, we make a series of histogram arrays in each of the four parameters where the histogram arrays are filled with the appropriate likelihood.
In particular, we will look at an ntuple value, and look at the ith parameter, then go to the appropriate hist array, find the corresponding histogram and fill that hist with the loglikelihood.
This will be significantly faster than iterating through the ntuple 4 times to distribute the appropriate loglikelihood values. Moreover, getting the minimum will be a simple ->GetMinimum() from each hist, which will then be plotted ina new histogram.
Requires that we know the precision for the spatial variables and the time variable and know where they start and end so that the histograms can be made apriori. Accordingly, will set constant variables in the beginning to fix that.
*/
{
#include <vector>     

  const int ntuplecount = 2;

  TObjArray xtemphists,ytemphists,ztemphists,ttemphists;    // Histogram array that will hold the hists that will be filled when each loglikelihood is distributed.
  
  TH1F* xhist_final;
  TH1F* yhist_final;
  TH1F* zhist_final;
  TH1F* thist_final;

  TF1* time_norm;
  TF1* spatial_norm;

  Int_t check,count,mincount,spatial_steps,time_steps,count2,x_index,y_index,z_index,t_index,ntuple_entries,gmin_xpos,gmin_ypos,gmin_zpos,gmin_time;
  Float_t xpos,ypos,zpos,time,loglikelihood,gmin_loglikelihood;

  // Read in the file and make an output file
  TFile TheFile("./fiTQun_single_events.root");
  TFile NewFile("./TESTING.root","RECREATE");

  // The Ntuples to be read from and written to
  TNtuple* Ntuple; 
  TNtuple* EventGlobalMin = new TNtuple("Global_Min", "Global Minimum","event:xpos:ypos:zpos:time:loglikelihood");

  // Values set apriori according to the fiTQun fitting parameters.
  const Int_t inf = 1000000000;
  const Int_t time_precision = 1;
  const Int_t time_min = 950;
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
    cout << "More spatial steps than time steps! Stopping here." << endl;
    break;
  }

  check = 0;
  count = 0;
  // Allocate the matrices to store the data
  vector<Float_t> time_mat;
  vector<Float_t> xpos_mat;
  vector<Float_t> ypos_mat;
  vector<Float_t> zpos_mat;

  // Loop through the ntuples
  for (int i = 1; i < ntuplecount; i++) {
    // Clear the vectors of the previous ntuple values
    time_mat.clear();
    xpos_mat.clear();
    ypos_mat.clear();
    zpos_mat.clear();
    // refill the vectors with infs so that we can compare values
    for (int m = 0; m < time_steps; m++) {
      switch(check) {
      case 0: 
	xpos_mat.push_back(inf);
	ypos_mat.push_back(inf);
	zpos_mat.push_back(inf);
      case 1:
	time_mat.push_back(inf);
	break;
      }
      if (m == spatial_steps - 1) { check = 1; }
    }
    TheFile->cd();
    cout << "Working on Ntuple number " << i << endl;
    Ntuple = (TNtuple*)TheFile.Get(Form("likelihoodprefit_Event%03d",i));    // Read from the next ntuple
    Ntuple->SetBranchAddress("xpos",&xpos);
    Ntuple->SetBranchAddress("ypos",&ypos);
    Ntuple->SetBranchAddress("zpos",&zpos);
    Ntuple->SetBranchAddress("time",&time);
    Ntuple->SetBranchAddress("loglikelihood",&loglikelihood);
    xhist_final = new TH1F(Form("loglikelihood_vs_xpos_Event%03d",i), "Loglikelihood vs xpos", spatial_steps+1, (float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    yhist_final = new TH1F(Form("loglikelihood_vs_ypos_Event%03d",i), "Loglikelihood vs ypos", spatial_steps+1, (float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    zhist_final = new TH1F(Form("loglikelihood_vs_zpos_Event%03d",i), "Loglikelihood vs zpos", spatial_steps+1, (float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    thist_final = new TH1F(Form("loglikelihood_vs_time_Event%03d",i), "Loglikelihood vs time", time_steps+1, (float)time_min - ((float)time_precision)/2.0, (float)time_max + ((float)spatial_precision)/2.0);
    ntuple_entries = Ntuple->GetEntries();
    check = 0;
    gmin_loglikelihood = inf;
    for (int k = 0; k < ntuple_entries; k++) {
      Ntuple->GetEntry(k);
      if (loglikelihood != loglikelihood) {loglikelihood = inf;} // Set nan loglikelihoods to inf
      //cout << "Loglikelihood at " << k << " is " << loglikelihood << endl; // Debugging
      // Keep track of the minimum loglikelihood observed
      if (loglikelihood < gmin_loglikelihood) {
	cout << "New minimum Loglikelihood is " << loglikelihood << endl;
	gmin_loglikelihood = loglikelihood;
	gmin_xpos = xpos;
	gmin_ypos = ypos;
	gmin_zpos = zpos;
	gmin_time = time;
      }
      // Fill the appropriate vectors with the loglikelihood if it is smaller than the prior one
      x_index = ((int)xpos - spatial_min)/spatial_precision;
      y_index = ((int)ypos - spatial_min)/spatial_precision;
      z_index = ((int)zpos - spatial_min)/spatial_precision;
      t_index = ((int)time - time_min)/time_precision;
      if (xpos_mat[x_index] > loglikelihood) {
	xpos_mat[x_index] = loglikelihood;
      }
      if (ypos_mat[y_index] > loglikelihood) {
	ypos_mat[y_index] = loglikelihood;
      } 
      if (zpos_mat[z_index] > loglikelihood) {
	zpos_mat[z_index] = loglikelihood;
      }
      if (time_mat[t_index] > loglikelihood) {
	time_mat[t_index] = loglikelihood;
      }
    }
    // Write the minimum of each column vector to the appropriate histogram bin
    check = 0;
    for (int j = 0; j < time_steps; j++) {
      switch(check) {
      case 0:
	xhist_final->SetBinContent(j,xpos_mat[j]);
	yhist_final->SetBinContent(j,ypos_mat[j]);
	zhist_final->SetBinContent(j,zpos_mat[j]);
      case 1:
	thist_final->SetBinContent(j,time_mat[j]);
	break;
      }
      if (j == spatial_steps - 1) { check = 1; }
    }
    NewFile.cd();
    cout << "Filling ntuple with global minimum position ..." << endl;
    // Save to the ntuple
    EventGlobalMin->Fill((float)i,(float)gmin_xpos,(float)gmin_ypos,(float)gmin_zpos,(float)gmin_time,(float)gmin_loglikelihood);
    // Normalize the plots that were made to the global minimum experienced during the process
    time_norm = new TF1(Form("time_norm%03d",i),"-1",time_min,time_max);
    spatial_norm = new TF1(Form("spatial_norm%03d",i),"-1",spatial_min,spatial_max);
    xhist_final->Add(spatial_norm,gmin_loglikelihood);
    yhist_final->Add(spatial_norm,gmin_loglikelihood);
    zhist_final->Add(spatial_norm,gmin_loglikelihood);
    thist_final->Add(time_norm,gmin_loglikelihood);    
    // Save the plots that we have made
    cout << "Saving plots ... " << endl;
    NewFile.cd();
    xhist_final->Write();
    yhist_final->Write();
    zhist_final->Write();
    thist_final->Write();
  }
  EventGlobalMin->Write();
}

