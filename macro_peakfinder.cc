{
  // Peakfinding macro sfcript for n events
#include <vector>

  const Int_t n = 100; // The number of events in the rootfile
  const Int_t bins = 160;
  const Double_t ratio = 0.4;
  const Int_t inf = 100000000000000000;

  TFile TheFile("fiTQun_output_nopileup.root");
  TFile NewFile("NoPileUp_TrueLoglikelihood.root","RECREATE");
  NewFile->mkdir("Time_likelihood");
  NewFile->mkdir("xpos_likelihood");
  NewFile->mkdir("ypos_likelihood");
  NewFile->mkdir("zpos_likelihood");
  //TH1F* TempHist = new TH1F("temphist","temphist", 50000, 0.0, 10000000.0);
  TH1F* timeHist,xHist,yHist,zHist;                                 // Take the loglikelihood distribution in all parameters to get a better idea of how many events we are experiencing and hence cut this way


  TNtuple* Ntuple;
  TNtuple* peaks = new TNtuple("Peak_Ntuple", "Peak Ntuple","event:time:loglikelihood:left:right");
  
  Float_t xpos,ypos,zpos,time,loglikelihood,temptime,val,max,min,templog;
  Int_t cut = 50,bin,index,minbin,check,tempcheck,tempflag,flip,tempflip,tempahead,tempbehind,tempindex,temppush;

  vector <Int_t> partition,peakflag,ahead,behind;
  
  temptime = 840;
  for (int i = 0; i < n; i++) {
    partition.clear();
    peakflag.clear();
    ahead.clear();
    behind.clear();
    //TempHist->Reset();
    TheFile->cd();
    Ntuple = (TNtuple*)TheFile.Get(Form("likelihoodprefit_Event%03d",i));
    Ntuple->SetBranchAddress("xpos",&xpos);
    Ntuple->SetBranchAddress("ypos",&ypos);
    Ntuple->SetBranchAddress("zpos",&zpos);
    Ntuple->SetBranchAddress("time",&time);
    Ntuple->SetBranchAddress("loglikelihood",&loglikelihood);
    timeHist = new TH1F(Form("loglikelihood_vs_time_Event%03d",i), "Loglikelihood vs Time",bins,840.0,1000.0);
    xHist = new TH1F(Form("loglikelihood_vs_x_Event%03d",i), "Loglikelihood vs x-Pos",32,-400.0,400.0);
    yHist = new TH1F(Form("loglikelihood_vs_y_Event%03d",i), "Loglikelihood vs y-Pos",32,-400.0,400.0);
    zHist = new TH1F(Form("loglikelihood_vs_z_Event%03d",i), "Loglikelihood vs z-Pos",32,-400.0,400.0);
    cout << "File number: " << i << endl;
    for (int j = 0; j < Ntuple->GetEntries(); j++) {
      Ntuple->GetEntry(j);
      if (temptime != time) {
	bin = timeHist->Fill(temptime);
	/*timeHist->SetBinContent(bin,TempHist->GetMean());
	timeHist->SetBinError(bin,TempHist->GetMeanError());
	TempHist->Reset();*/
	temptime = time;
      }
      if (loglikelihood == loglikelihood && loglikelihood < inf) {
	//TempHist->Fill(loglikelihood);
	bin = xHist->Fill(xpos);                                                    // Recursively adds the loglikelihoods to the position bin
	xHist->SetBinContent(bin,xHist->GetBinContent(bin) + loglikelihood);
	bin = yHist->Fill(ypos);                                              
	yHist->SetBinContent(bin,yHist->GetBinContent(bin) + loglikelihood);
	bin = zHist->Fill(zpos);                                              
	zHist->SetBinContent(bin,zHist->GetBinContent(bin) + loglikelihood);
      }
    }
    cout << "TimeHist Entries:  " << timeHist->GetEntries() << endl;
    NewFile->cd("Time_likelihood");
    timeHist->Write();
    NewFile->cd("xpos_likelihood");
    xHist->Write();
    NewFile->cd("ypos_likelihood");
    yHist->Write();
    NewFile->cd("zpos_likelihood");
    zHist->Write();
  }
}
/*
    ////////////////////////////////////////// Stopping the code here while figuring out the plotting //////////////////////////////////////////////////////

    minbin = timeHist->GetMinimumBin();
    min = timeHist->GetBinContent(minbin);
    
    cut = min*ratio;
    
    check = 0;
    tempcheck = 0;
    
    for (int j = 0; j < bins; j++) {
      templog = timeHist->GetBinContent(j);
      if (templog <= cut) {
	if (j == minbin) { minbin = index; }
	check = 1;
      } else { check = 0; }
      if (tempcheck != check) { 
	tempcheck = check;
	partition.push_back(j);
      }
    }
    // With two vectors that are storing our values that passed the threshold, we can now start to actually define the places we will fit
    // First thing to check will be the minimum
    
    // Check if the point following and before the min are one time interval away
    
    
    if ( timeHist->GetBinContent(partition[0]) < timeHist->GetBinContent(partition[0] + 1)) { 
      tempflag = partition[0];
      tempflip = 0;
      behind.push_back(0);
      tempbehind = 0;
      tempahead = 1;
    } else if (timeHist->GetBinContent(partition[0]) > timeHist->GetBinContent(partition[0] + 1)) {
      tempflip = 1;
      tempflag = partition[0] + 1;
      tempbehind = 1;
      tempahead = 0;
    } else {
      tempflip = 2;
      tempflag = partition[0] + 1;
      tempbehind = 0;
      tempahead = 0;
    }
    index = 0;
    tempindex = index;
    check = partition[index] + 1;
    temppush = 0;

    //cout << "Partition size is " << partition.size() << endl;

    cout << "Peak(s) at :" << endl;   
    while (index < partition.size()) { 
      if (tempindex != index) {
	check = partition[index];
	tempindex = index;
      }
      if (check > 1000) {
	cout << "Peak loop is broken. Terminating ... " << endl;
	break;
      }
      if (timeHist->GetBinContent(check) < timeHist->GetBinContent(check+1)) {
	flip = 0;
      } else if (timeHist->GetBinContent(check) > timeHist->GetBinContent(check+1)) {
	flip = 1;
	tempflag = check+1;
      } else {
	flip = 2;
	tempflag = check+1;
	temppush = 0;
      }
      //cout << "First Bin: " << Hist->GetBinContent(check) << " Second Bin: " << Hist->GetBinContent(check+1) << " Flip: " << flip << " tempflip: " << tempflip << endl; // -Debugging
      //cout << check << " , partition: " << partition[index + 1] << endl;
      switch(flip) {
      case 0: if (tempflip != 0) {
	  behind.push_back(tempbehind);
	  tempbehind = 0;
	  tempflip = 0;
	  peakflag.push_back(tempflag);
	  temppush = 1;
	} 
	tempahead++;
	break;
      case 1: if (tempflip != 1) {
	  ahead.push_back(tempahead);
	  tempahead = 0;
	  tempflip = 1;
	}
	tempbehind++;
	break;
      case 2: break;
      }
      check++;
      if (check - 1 == partition[index + 1]) { 
	index += 2; 
	if (temppush == 0) { peakflag.push_back(tempflag); }
      }
      if (timeHist->GetBinContent(check+2) > 0) { index = 1000; } // If the returned value of the histogram is above 0, we know that it must be the last point in our hist 
    }

    for (int m = 0; m < peakflag.size(); m++) {
      cout << "Time: " << 800 + 15*peakflag[m] << " Loglikelihood: " << timeHist->GetBinContent(peakflag[m]) << endl;
      peaks.Fill(i,800 + 15*peakflag[m],timeHist->GetBinContent(peakflag[m]),0.,0.);
    }
  }
  NewFile->cd();
  peaks.Write();
}
