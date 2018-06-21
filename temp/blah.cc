{
  // A script meant to look at the variance in the variance of the loglikelihood at a given time and see if we can define a cut based on this value
  // Also used as a general script for looking at particular statistics


  TFile NewFile("BlahFile.root","RECREATE");
  TFile TheFile("./fiTQun_output_nopileup.root");
  TH1D* FinalHist = new TH1D("lglklihood","Loglkelihood vs Time", 600, 800.0, 1400.0, 150, -15.0, 0.0);
  TH1D* MeanHist = new TH1D("MHist", "Histogram of Loglikelihood Pulse Mean", 150, -15.0, 0.0);
  TH1D* RelHist = new TH1D("RelInt", "Relative Integral to cut-off", 50, 0.0, 1.0);  

  //TH1D* TempHist;
  TF1* func;

  TNtuple* Ntuple;
  TObjArray TempHists;
  Float_t xpos,ypos,zpos,time,loglikelihood,temptime,val,max;
  Int_t cut = 50;
  
  temptime = 800;
  for (int i = 0; i < 100; i++) {
    TheFile->cd();
    Ntuple = (TNtuple*)TheFile.Get(Form("likelihoodprefit_Event%03d",i));
    Ntuple->SetBranchAddress("xpos",&xpos);
    Ntuple->SetBranchAddress("ypos",&ypos);
    Ntuple->SetBranchAddress("zpos",&zpos);
    Ntuple->SetBranchAddress("time",&time);
    Ntuple->SetBranchAddress("loglikelihood",&loglikelihood);
    TempHist = new TH1D(Form("loglikelihooddist_Event%03d",i), "Loglikelihood Distribution",100,-1.0,0.0);
    cout << "File number: " << i << endl;
    max = -Ntuple->GetMinimum("loglikelihood");
    for (int j = 0; j < Ntuple->GetEntries(); j++) {
      Ntuple->GetEntry(j);
      if (loglikelihood == loglikelihood && loglikelihood < 0.01 && temptime == time) {
	TempHist->Fill(loglikelihood/max);
      /*if (temptime != time) {
	TempHists->Add(TempHist);
	TempHist = new TH1D(Form);
	}*/
      }
    }
    RelHist->Fill((TempHist->Integral(0,cut))/(TempHist->Integral(0,101)));
    TempHists.Add(TempHist);
    NewFile->cd();
    TempHist->Write();
    //TempHist->Fit("gaus","","",-15,cut);   // Fit a Gaussian based on some cut of the loglikelihood
    //func = TempHist->GetFunction("gaus");  // Use a function to call and read the fit values into a histogram
    //MeanHist->Fill(func->GetParameter(1));
  }


  //MeanHist->Write();
  RelHist->Write();
}
