{
  // Used to plot the time vs likelihood using the means of each time bin
  // Also finds the maxima/minima using a simple binary comparison and likelihood cut
  // Modified to work with the true loglikelihood rather than the pseudo-likelihood made with a sum rather than a product

#include <vector>
  
  int bins, t;
  
  bins = 160;

  TFile TheFile("fiTQun_newlikelihood.root");
  TH1F* TempHist;
  TH1F* FinalHist = new TH1F("lglklihood","Loglikelihood vs Time", bins, 840.0, 1000.0);
  TNtuple* Ntuple;

  Float_t time,loglikelihood,temptime,templog,min,cut;
  Int_t bin,index,minbin,check,tempcheck;
  const Int_t inf = 100000000000000000;

  vector <Float_t> time_vec,likelihood_vec;
  vector <Int_t> partition;

  Ntuple = (TNtuple*)TheFile.Get("likelihoodprefit_Event000");
  Ntuple->SetBranchAddress("time",&time);
  Ntuple->SetBranchAddress("loglikelihood",&loglikelihood);
  

  temptime = 840;
  TempHist = new TH1F(Form("time_%d",temptime),Form("Time %d",temptime), 50000, 0.0, 10000000.0);
  check = 0;
  tempcheck = 0;
  for (int i = 0; i < Ntuple->GetEntries(); i++) {
    Ntuple->GetEntry(i);
    if (loglikelihood == loglikelihood && temptime == time && loglikelihood < inf) {
      TempHist->Fill(loglikelihood);
    } else if (temptime != time) {
      cout << "Time: " << time << " Mean: " << TempHist->GetMean() << " Entries: " << TempHist->GetEntries() << endl;
      check = 0;
      tempcheck = 0;
      bin = FinalHist->Fill(temptime);
      FinalHist->SetBinContent(bin,TempHist->GetMean());
      FinalHist->SetBinError(bin,TempHist->GetMeanError());
      TempHist = new TH1F(Form("time_%03d",(int) time),Form("Time %03d",(int) time), 50000, 0.0, 10000000.0);
      temptime = time;
    }
  }
  TFile NewFile("./Converted_Out.root", "RECREATE");
  FinalHist->Write();
  break;
}  

  // Now to see if we can fit the plot
  // Pick 0.5 the max value to be the threshold for if a pulse is truly an event. Could probably find a better choice for such a number but that is for later
  
  minbin = FinalHist->GetMinimumBin();
  min = FinalHist->GetBinContent(minbin);
  
  cut = min*0.5;
  
  check = 0;
  tempcheck = 0;

  for (int j = 0; j < bins; j++) {
    templog = FinalHist->GetBinContent(j);
    if (templog <= cut) {
      if (j == minbin) { minbin = index; }
      check = 1;
    } else { check = 0; }
    if (tempcheck != check) { 
      tempcheck = check;
      partition.push_back(j);
      cout << "Partition at time " << 800 + 15*j << " with likelihood " << templog << endl;
    }
  }
  
  // With two vectors that are storing our values that passed the threshold, we can now start to actually define the places we will fit
  // First thing to check will be the minimum
 
  // Check if the point following and before the min are one time interval away


  vector<Int_t> peakflag,ahead,behind;
  Int_t tempflag,flip,tempflip,tempahead,tempbehind,tempindex;                   // 'flip' used to check if we changing from increasing or decreasing

  if ( FinalHist->GetBinContent(partition[0]) < FinalHist->GetBinContent(partition[0] + 1)) { 
    tempflag = partition[0];
    tempflip = 0;
    behind.push_back(0);
    tempbehind = 0;
    tempahead = 1;
  } else { 
    tempflip = 1;
    tempflag = partition[0] + 1;
    tempbehind = 1;
    tempahead = 0;
  }
  cout << "First tempflip is " << tempflip << endl;
  index = 0;
  tempindex = index;
  check = partition[index] + 1;
  
  while (index < partition.size()) { 
    if (tempindex != index) {
      check = partition[index];
      tempindex = index;
    }
    if (FinalHist->GetBinContent(check) < FinalHist->GetBinContent(check+1)) {
      flip = 0;
    } else {
      flip = 1;
      tempflag = check+1;
    }
    cout << "Flip = " << flip << " at " << check << endl;
    switch(flip) {
    case 0: if (tempflip != 0) {
	behind.push_back(tempbehind);
	tempbehind = 0;
	tempflip = 0;
	peakflag.push_back(tempflag);
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
    }
      check++;
      if (check - 1 == partition[index + 1]) { index += 2; }
  }

  cout << "Peaks at: " << endl;
  for (int n = 0; n < peakflag.size(); n++) {
    cout << "Time: " << 800 + 15*(peakflag[n]) << " Likelihood: " << FinalHist->GetBinContent(peakflag[n]) << endl;
  }
}
