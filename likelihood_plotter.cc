/*
The Goal of this plotting script is to consider the 1D case in each parameter (except for likelihood) and find the minimum likelihood in each iteration for that variable while varying the other parameters.
To do this, we make a series of histogram arrays in each of the four parameters where the histogram arrays are filled with the appropriate likelihood.
In particular, we will look at an ntuple value, and look at the ith parameter, then go to the appropriate hist array, find the corresponding histogram and fill that hist with the loglikelihood.
This will be significantly faster than iterating through the ntuple 4 times to distribute the appropriate loglikelihood values. Moreover, getting the minimum will be a simple ->GetMinimum() from each hist, which will then be plotted ina new histogram.
Requires that we know the precision for the spatial variables and the time variable and know where they start and end so that the histograms can be made apriori. Accordingly, will set constant variables in the beginning to fix that.
*/
{
#include <vector>     

  TObjArray xtemphists;    // Histogram array that will hold the hists that will be filled when each loglikelihood is distributed.
  TObjArray ytemphists;
  TObjArray ztemphists;
  TObjArray ttemphists;
  
  Int_t check,count,mincount;
  Float_t xpos,ypos,zpos,time,loglikelihood;

  // Read in the ntuple
  TNtuple Ntuple;

  // Values set apriori according to the fiTQun fitting parameters.
  const Int_t time_precision = 1;
  const Int_t time_min = 950;
  const Int_t time_max = 1000;
  const Int_t spatial_min = -400;
  const Int_t spatial_max = 400;
  const Float_t spatial_precision = 25.;

  // Figure out which index has the most indices to run through.
  check = ((float) time_max - time_min)/((float) time_precision);
  count = ((float) spatial_max - spatial_min)/((float) spatial_precision);

  if (check > count) {
    mincount = check;
  } else { mincount = count;}

  // Make the histograms and fill the hist array
  check = 0;
  count = 0;
  while (count < mincount) {
    
