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
  
  TH1F* xhist_final;
  TH1F* yhist_final;
  TH1F* zhist_final;
  TH1F* thist_final;

  // 2Dim plots of the same form as the 1D ones, except that the final plot will hold 2 parameters constant and will then spit out the projection among the other 2
  TH2F* XYhist;
  TH2F* XZhist;
  TH2F* YZhist;
  TH2F* TXhist;
  TH2F* TYhist;
  TH2F* TZhist;

  TF1* time_norm;
  TF1* spatial_norm;

  Int_t check,count,mincount,spatial_steps,time_steps,count2,x_index,y_index,z_index,t_index,ntuple_entries,gmin_xpos,gmin_ypos,gmin_zpos,gmin_time;
  Float_t xpos,ypos,zpos,time,loglikelihood,gmin_loglikelihood;

  // Read in the file and make an output file
  TFile TheFile("./fiTQun_3events_15ns.root");
  TFile NewFile("./testing.root","RECREATE");

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
  // Allocate the vectors to store the data
  vector<Float_t> time_mat;
  vector<Float_t> xpos_mat;
  vector<Float_t> ypos_mat;
  vector<Float_t> zpos_mat;

  // Initialize the matrices that store the 2 dim data
  Float_t** xy_mat = new Float_t*[spatial_steps];
  Float_t** xz_mat = new Float_t*[spatial_steps];
  Float_t** yz_mat = new Float_t*[spatial_steps];
  Float_t** tx_mat = new Float_t*[time_steps];
  Float_t** ty_mat = new Float_t*[time_steps];
  Float_t** tz_mat = new Float_t*[time_steps];

  xy_mat[0] = new Float_t[spatial_steps*spatial_steps];
  xz_mat[0] = new Float_t[spatial_steps*spatial_steps];
  yz_mat[0] = new Float_t[spatial_steps*spatial_steps];
  tx_mat[0] = new Float_t[spatial_steps*time_steps];
  ty_mat[0] = new Float_t[spatial_steps*time_steps];
  tz_mat[0] = new Float_t[spatial_steps*time_steps];

  for (int i = 1; i < time_steps; i++) {
    if (i < spatial_steps) {
      xy_mat[i] = xy_mat[0] + i*spatial_steps;
      xz_mat[i] = xz_mat[0] + i*spatial_steps;
      yz_mat[i] = yz_mat[0] + i*spatial_steps;
    }
    tx_mat[i] = tx_mat[0] + i*spatial_steps;
    ty_mat[i] = ty_mat[0] + i*spatial_steps;
    tz_mat[i] = tz_mat[0] + i*spatial_steps;
  }

  // Loop through the ntuples
  for (int i = 0; i < ntuplecount; i++) { 
    // Clear the vectors of the previous ntuple values
    time_mat.clear();
    xpos_mat.clear();
    ypos_mat.clear();
    zpos_mat.clear();
    // refill the vectors with infs so that we can compare values
    
    for (int m = 0; m < time_steps; m++) {
      if (m < spatial_steps) {
	xpos_mat.push_back(inf);
	ypos_mat.push_back(inf);
	zpos_mat.push_back(inf);
      }
      time_mat.push_back(inf);
    }

    // fill with infs
    for (int k = 0; k < time_steps; k++) {
      for (int l = 0; l < spatial_steps; l++) {
	if (k < spatial_steps) {
	  xy_mat[k][l] = inf;
	  xz_mat[k][l] = inf;
	  yz_mat[k][l] = inf;
	}
	tx_mat[k][l] = inf;
	ty_mat[k][l] = inf;
	tz_mat[k][l] = inf;
      }
    }
    TheFile.cd();
    cout << "Working on Ntuple number " << i << endl;
    // Make the 1D plots for the new ntuple
    xhist_final = new TH1F(Form("loglikelihood_vs_xpos_Event%03d",i), "Loglikelihood vs xpos", spatial_steps, (float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    yhist_final = new TH1F(Form("loglikelihood_vs_ypos_Event%03d",i), "Loglikelihood vs ypos", spatial_steps, (float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    zhist_final = new TH1F(Form("loglikelihood_vs_zpos_Event%03d",i), "Loglikelihood vs zpos", spatial_steps, (float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    thist_final = new TH1F(Form("loglikelihood_vs_time_Event%03d",i), "Loglikelihood vs time", time_steps, (float)time_min - ((float)time_precision)/2.0, (float)time_max + ((float)spatial_precision)/2.0);
    // Make the 2D plots for the new ntuple
    XYhist = new TH2F(Form("loglikelihood_XY_Event%03d",i), "XY Plot",spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0,spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    XZhist = new TH2F(Form("loglikelihood_XZ_Event%03d",i), "XZ Plot",spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0,spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    YZhist = new TH2F(Form("loglikelihood_YZ_Event%03d",i), "YZ Plot",spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0,spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    TXhist = new TH2F(Form("loglikelihood_TX_Event%03d",i), "TX Plot",time_steps,(float)time_min - ((float)time_precision)/2.0, (float)time_max + ((float)time_precision)/2.0,spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    TYhist = new TH2F(Form("loglikelihood_TY_Event%03d",i), "TY Plot",time_steps,(float)time_min - ((float)time_precision)/2.0, (float)time_max + ((float)time_precision)/2.0,spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    TZhist = new TH2F(Form("loglikelihood_TZ_Event%03d",i), "TZ Plot",time_steps,(float)time_min - ((float)time_precision)/2.0, (float)time_max + ((float)time_precision)/2.0,spatial_steps,(float)spatial_min - ((float)spatial_precision)/2.0, (float)spatial_max + ((float)spatial_precision)/2.0);
    // Set-up to read from the ntuple
    Ntuple = (TNtuple*)TheFile.Get(Form("likelihoodprefit_Event%03d",i));
    Ntuple->SetBranchAddress("xpos",&xpos);
    Ntuple->SetBranchAddress("ypos",&ypos);
    Ntuple->SetBranchAddress("zpos",&zpos);
    Ntuple->SetBranchAddress("time",&time);
    Ntuple->SetBranchAddress("loglikelihood",&loglikelihood);
    ntuple_entries = Ntuple->GetEntries();
    check = 0;
    gmin_loglikelihood = inf;
    for (int k = 0; k < ntuple_entries; k++) {
      Ntuple->GetEntry(k);
      if (loglikelihood != loglikelihood) {loglikelihood = inf;} // Set nan loglikelihoods to inf
      //cout << "Loglikelihood at " << k << " is " << loglikelihood << endl; // Debugging
      // Keep track of the minimum loglikelihood observed
      if (loglikelihood < gmin_loglikelihood) {
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
      // Now to fill the matrices with the minimum loglikelihood
      if (xy_mat[x_index][y_index] > loglikelihood) {
	xy_mat[x_index][y_index] = loglikelihood;
      }
      if (xz_mat[x_index][z_index] > loglikelihood) {
	xz_mat[x_index][z_index] = loglikelihood;
      }
      if (yz_mat[y_index][z_index] > loglikelihood) {
	yz_mat[y_index][z_index] = loglikelihood;
      }
      if (tx_mat[t_index][x_index] > loglikelihood) {
	tx_mat[t_index][x_index] = loglikelihood;
      }
      if (ty_mat[t_index][y_index] > loglikelihood) {
	ty_mat[t_index][y_index] = loglikelihood;
      }
      if (tz_mat[t_index][z_index] > loglikelihood) {
	tz_mat[t_index][z_index] = loglikelihood;
      }
    }
    // Write the minimum of each column vector to the appropriate histogram bin
    check = 0;
    // First the 1D hists
    for (int j = 0; j < time_steps; j++) {
      if (j < spatial_steps) {
	if (xpos_mat[j] == inf) { 
	  xhist_final->SetBinContent(j,0.0);
	} else {
	  xhist_final->SetBinContent(j,xpos_mat[j] - gmin_loglikelihood);
	}
	if (ypos_mat[j] == inf) {
	  yhist_final->SetBinContent(j,0.0);
	} else {
	  yhist_final->SetBinContent(j,ypos_mat[j] - gmin_loglikelihood);
	}
	if (zpos_mat[j] == inf) {	
	  zhist_final->SetBinContent(j,0.0);
	} else {
	zhist_final->SetBinContent(j,zpos_mat[j] - gmin_loglikelihood);
	}
      }
      if (time_mat[j] == inf) {
	      thist_final->SetBinContent(j,0.0);
      } else {
	thist_final->SetBinContent(j,time_mat[j] - gmin_loglikelihood);
      }
    }
    // Now to draw into the 2D hists
    for (int row = 0; row < time_steps; row++) {
      for (int col = 0; col < spatial_steps; col++) {
	if (row < spatial_steps) {
	  if (xy_mat[row][col] == inf) {
	    XYhist->Fill(spatial_min + row*spatial_precision,spatial_min + col*spatial_precision,0.0);
	  } else {
	    XYhist->Fill(spatial_min + row*spatial_precision,spatial_min + col*spatial_precision,xy_mat[row][col] - gmin_loglikelihood);
	  }
	  if (xz_mat[row][col] == inf) {	  
	    XZhist->Fill(spatial_min + row*spatial_precision,spatial_min + col*spatial_precision,0.0);
	  } else {
	    XZhist->Fill(spatial_min + row*spatial_precision,spatial_min + col*spatial_precision,xz_mat[row][col] - gmin_loglikelihood);
	  }
	  if (yz_mat[row][col] == inf) {	  	  
	    YZhist->Fill(spatial_min + row*spatial_precision,spatial_min + col*spatial_precision,0.0);
	  } else {
	    YZhist->Fill(spatial_min + row*spatial_precision,spatial_min + col*spatial_precision,yz_mat[row][col] - gmin_loglikelihood);
	  }
	}
	if (tx_mat[row][col] == inf) {
	  TXhist->Fill(time_min + row*time_precision,spatial_min + col*spatial_precision,0.0);
	} else {
	  TXhist->Fill(time_min + row*time_precision,spatial_min + col*spatial_precision,tx_mat[row][col] - gmin_loglikelihood);
	}
	if (ty_mat[row][col] == inf) {
	  TYhist->Fill(time_min + row*time_precision,spatial_min + col*spatial_precision,0.0);
	} else {
	  TYhist->Fill(time_min + row*time_precision,spatial_min + col*spatial_precision,ty_mat[row][col] - gmin_loglikelihood);
	}
	if (tz_mat[row][col] == inf) {
	  TZhist->Fill(time_min + row*time_precision,spatial_min + col*spatial_precision,0.0);
	} else {
	  TZhist->Fill(time_min + row*time_precision,spatial_min + col*spatial_precision,tz_mat[row][col] - gmin_loglikelihood);
	}
      }
    }
    NewFile.cd();
    cout << "Filling ntuple with global minimum position ..." << endl;
    // Save to the ntuple
    EventGlobalMin->Fill((float)i,(float)gmin_xpos,(float)gmin_ypos,(float)gmin_zpos,(float)gmin_time,gmin_loglikelihood);
    // Save the plots that we have made
    cout << "Saving plots ... " << endl;
    NewFile.cd();
    xhist_final->Write();
    yhist_final->Write();
    zhist_final->Write();
    thist_final->Write();
    XYhist->Write();
    XZhist->Write();
    YZhist->Write();
    TXhist->Write();
    TYhist->Write();
    TZhist->Write();
  }
  EventGlobalMin->Write();
  // Clean up the matrices
  delete [] xy_mat[0];
  delete [] xy_mat;
  delete [] xz_mat[0];
  delete [] xz_mat;
  delete [] yz_mat[0];
  delete [] yz_mat;
  delete [] tx_mat[0];
  delete [] tx_mat;
  delete [] ty_mat[0];
  delete [] ty_mat;
  delete [] tz_mat[0];
  delete [] tz_mat;
}

