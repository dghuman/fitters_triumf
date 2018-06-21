{
  /* A Short Script for plotting from the fiTQun NTuple written by me as well.
     Makes 3 plots per time interval; XY, XZ and YZ, which fill the z-coordinate in the plot with the likelihood value.
   */
  Float_t xpos,ypos,zpos,time,loglikelihood,tmptime;
  Int_t counter,index,entries;

  TH2D* hist0;
  TH2D* hist1;
  TH2D* hist2;
  TFile blah("./wcsim_output_fiTQun.root");
  TObjArray histsXY;
  TObjArray histsXZ;
  TObjArray histsYZ;
  TObjArray canvas0;
  TObjArray canvas1;
  TObjArray canvas2;
  TCanvas* c0;
  TCanvas* c1;
  TCanvas* c2;

  gStyle->SetOptStat(0);

  float win_scale=0.75;
  int n_wide=2;
  int n_high=3;
  
  // Setting the ntuple read in values
  likelihoodprefit_Event000->SetBranchAddress("xpos",&xpos);
  likelihoodprefit_Event000->SetBranchAddress("ypos",&ypos);
  likelihoodprefit_Event000->SetBranchAddress("zpos",&zpos);
  likelihoodprefit_Event000->SetBranchAddress("time",&time);
  likelihoodprefit_Event000->SetBranchAddress("loglikelihood",&loglikelihood);

  index = 0;
  entries = likelihoodprefit_Event000->GetEntries();
  // Now to fill the histogram with the appropriate values
  for (int i=0;  i<entries; i++) {
    likelihoodprefit_Event000->GetEntry(i);                               
    if (loglikelihood != loglikelihood) { index++;continue;}
    if (i == index) {
      hist0 = new TH2D(Form("XY_time_%d",(int) time),Form("XY_time_%d",(int) time),32,-400.0,400.0,32,-400.0,400.0);
      hist1 = new TH2D(Form("XZ_time_%d",(int) time),Form("XZ_time_%d",(int) time),32,-400.0,400.0,32,-400.0,400.0);
      hist2 = new TH2D(Form("YZ_time_%d",(int) time),Form("YZ_time_%d",(int) time),32,-400.0,400.0,32,-400.0,400.0);
      tmptime = time;
      cout << "Time: " << time << " Entries Removed: " << index << endl;
    }
    hist0->Fill(xpos, ypos, loglikelihood);
    hist1->Fill(xpos, zpos, loglikelihood);
    hist2->Fill(ypos, zpos, loglikelihood);
    if ((tmptime != time || i == entries - 1) && i != index) {
      cout << "Time: " << time << " Term: " << i << " Entries Removed: " << index << " Entry: " << i << endl;
      tmptime = time; 
      histsXY.Add(hist0);
      histsXZ.Add(hist1);
      histsYZ.Add(hist2);
      hist0 = new TH2D(Form("XY_time_%d",(int) time),Form("XY_time_%d",(int) time),32,-400.0,400.0,32,-400.0,400.0);
      hist1 = new TH2D(Form("XZ_time_%d",(int) time),Form("XZ_time_%d",(int) time),32,-400.0,400.0,32,-400.0,400.0);
      hist2 = new TH2D(Form("YZ_time_%d",(int) time),Form("YZ_time_%d",(int) time),32,-400.0,400.0,32,-400.0,400.0);
    }
  }

  // Make a new File to save the plots to
  TFile NewFile("PlotFile.root","RECREATE");
  gStyle->SetOptStat(0);

  // Pull from the histogram TObjectArray and draw into a canvas 6 times before makeing a new canvas to draw into, and writing everything to the RootFile.
  index = 0;
  cout << "HISTCOUNT: XY = " << histsXY.GetEntries() << " XZ = " << histsXZ.GetEntries() << " YZ = " << histsYZ.GetEntries() << endl;
  for (int j=0; j<histsXY.GetEntries(); j++) {
    if (j % 6 == 0) {
      canvas0.Add(c0);
      canvas1.Add(c1);
      canvas2.Add(c2);
      c0 = new TCanvas(Form("XY_canvas_%d",index),Form("XY_canvas_%d",index),700*n_wide*win_scale,500*n_high*win_scale);
      c1 = new TCanvas(Form("XZ_canvas_%d",index),Form("XZ_canvas_%d",index),700*n_wide*win_scale,500*n_high*win_scale);
      c2 = new TCanvas(Form("YZ_canvas_%d",index),Form("YZ_canvas_%d",index),700*n_wide*win_scale,500*n_high*win_scale);
      c0->Divide(n_wide,n_high);
      c1->Divide(n_wide,n_high);
      c2->Divide(n_wide,n_high);
      index++;
    }
    c0->cd((j % 6) + 1);
    histsXY[j]->Draw("colz");
    histsXY[j].Write();
    c1->cd((j % 6) + 1);
    histsXZ[j]->Draw("colz");
    histsXZ[j].Write();
    c2->cd((j % 6) + 1);
    histsYZ[j]->Draw("colz");
    histsYZ[j].Write();
    }


  for (int k=0; k<canvas0.GetEntries(); k++) {
    canvas0[k+1].Write();
    canvas1[k+1].Write();
    canvas2[k+1].Write();
  }
}
