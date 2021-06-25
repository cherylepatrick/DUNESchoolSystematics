// To run this, type: cafe Systematics3Solution.C

// These are standard header files from the CAFAna analysis tool
// They allow you to load and plot variables for each interaction event
// in your simulation file (and later, in data files)
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/ISyst.h" // *** this one is new for systematics!

// As we are working with simulation, we have access to information about the
// TRUE event - what GENIE, the simulation program, simulated for the event
// this is separate to the RECONSTRUCTED information - what the simulation 
// program thinks the DUNE detector would have seen
#include "CAFAna/Vars/Vars.h" // Variables
#include "CAFAna/Cuts/TruthCuts.h" // Cuts

#include "StandardRecord/SRProxy.h" // A wrapper for the CAF format

// These files come from the ROOT data analysis package
// This is used in many particle-physics experiments to make plots
// and do some basic statistics, cuts etc. The CAF simulation and data files
// are a special DUNE-specific format or a ROOT data file
// ROOT classes all start with a T, for some reason. It makes them easy to search
// for, except for things like the unfortunate TAxis...
#include "TCanvas.h" // Plots are drawn on a "canvas"
#include "TH1.h" // 1-dimensional histogram
#include "TPad.h" // Canvases are divided into pads. That could let you draw more than one plot on a canvas, if you wanted to, by using multiple pads. We will not bother with that today.
#include "TLegend.h" // Lets us draw a legend
#include "TMath.h" // I'll use some basic math functions
#include "TRandom3.h" // *** Random number generator

// Standard C++ library for input and output
#include <iostream>


/* *****************
 Define some GENIE interaction modes.
 See full list at https://wiki.dunescience.org/wiki/Scattering_mode
 Use these to make your TRUTH cuts - this the interaction type GENIE simulated
 */

const int MODE_QE = 1;
const int MODE_RES = 4;
const int MODE_DIS = 3;
const int MODE_MEC = 10;

/* ********
 Define some PDG codes (particle identifiers from https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
 You can also find the ID's for things like protons, neutrons, pions and even whole nuclei in the list!
 */
const int PDG_MU=13;
const int PDG_E=11;
const int PDG_NUMU=14;
const int PDG_NUE=12;

/* *********
 Define some physical constants
 */
const double M_P = .938; // Proton mass in GeV
const double M_N = .939; // Neutron mass in GeV
const double M_MU = .106; // Muon mass in GeV
const double E_B = .028; // Binding energy for nucleons in argon-40 in GeV

using namespace ana;
using util::sqr;

// **** function to make a fractional plot
TH1D *MakeFractionalPlot( TH1D* shifted, TH1D *cv)
{
  TH1D *frac= (TH1D*) shifted->Clone();
  frac->Add(cv,-1);
  frac->Divide(cv);
  return frac;
}


// Define the quasi-elastic formula for neutrino energy
double QEFormula(double Emu, double cosmu) // Muon energy and cosine of muon angle
{
  //Muon momentum
  const double pmu = sqrt(sqr(Emu) - sqr(M_MU)); // Use the relativity formula E^2 = p^2 + m^2
  // This is the neutrino-mode version of the formula. For antineutrino mode, swap neutron and proton masses.
  const double num = sqr(M_P) - sqr(M_N - E_B) - sqr(M_MU) + 2 * (M_N - E_B) * Emu;
  const double denom = 2 * (M_N - E_B - Emu + pmu * cosmu);
  if (denom==0) return 0;
  return num/denom;
}

// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Systematics3Solution()
{
  // CAFs
  const std::string CAFS = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_FHC_90*.root"; //This wildcard gives 10 files!
  
  // Source of events - load them from the one of the sets of files
  SpectrumLoader loader(CAFS);

  // We want to plot a histogram with 40 bins, covering the range 0 to 10 GeV
  const Binning binsEnergy = Binning::Simple(40, 0, 10);

  // This time we will plot neutrino energy from the QE reconstruction formula
  const Var kRecoQEFormulaEnergy([](const caf::SRProxy* sr)
                                 {
                                   const double Emu = sr->Elep_reco;

                                    // NB reco doesn't always work out!
                                    // Our shifts might have made it negative, that doesn't seem real...
                                   if(Emu < M_MU) return 0.;
                                   const double cosmu = cos(sr->theta_reco);
                                   if (isnan(Emu) || isnan(cosmu))return 0.;
                                   return QEFormula(Emu, cosmu);
                                 });
  const HistAxis axRecoQEFormula("Reconstructed QE energy (GeV)", binsEnergy, kRecoQEFormulaEnergy);

  
  // Now we have defined an axis and a variable we want to plot on it, let's decide which
  // events to plot. So here we are telling the loader to load a spectrum with our defined axis
  // Here, we are defining the selections or "cuts" using the variable names in CAFAna https://wiki.dunescience.org/wiki/CAFAna_Cuts
  
  // We'll use the CC0pi cut from last class
  const Cut kHasCC0PiFinalState([](const caf::SRProxy* sr)
                                {
                                  const int totPi = sr->nipip + sr->nipim + sr->nipi0;
                                  return abs(sr->LepPDG) == 13 && sr->nP >= 1 && totPi == 0;
                                });
  

  // Define the central-value Spectrum
  Spectrum sCV(loader, axRecoQEFormula, kHasCC0PiFinalState && kRecoQEFormulaEnergy>0);

  // Define a class to make the systematic energy shift
   // In this case it scales the muon energy by +/- 20 %
   class EMuScale: public ISyst
   {
   public:
     EMuScale(): ISyst("muScale", "Muon energy scale") {}

     virtual void Shift(double sigma, // This represents how much we shift by. We'll pass in +1 to shift up and -1 to shift down
                        Restorer& restore, // To enable us to change back
                        caf::SRProxy* sr, // Standard record
                        double& weight) const override // You can reweight
     {
       // Tell CAFAna we're about to shift this variable, so that it can put it
       // back again afterwards.
       restore.Add(sr->Elep_reco);
       // Now we can modify the muon energy
       sr->Elep_reco *= (1 + 0.2 * sigma); // Increase (sigma=+1) or decrease (sigma = -1) muon energy by 20 %
     }
   };
   
   // Make an object of the EMuScale shifter type we just defined
   EMuScale kEMuScale;
   SystShifts ssScaleUp(&kEMuScale, +1); // Scale UP by 20%
   SystShifts ssScaleDn(&kEMuScale, -1); // Scale DOWN by 20%
   
   // Make Spectrum objects for the shifted energies.
   // Note the extra parameter to indicate the shift (the rest stays the same as for the central value)
   Spectrum sScaleUp(loader, axRecoQEFormula, kHasCC0PiFinalState && kRecoQEFormulaEnergy>0, ssScaleUp);
   Spectrum sScaleDn(loader, axRecoQEFormula, kHasCC0PiFinalState && kRecoQEFormulaEnergy>0, ssScaleDn);
  
  // 20% smear
  class EMuSmear: public ISyst
  {
  public:
    EMuSmear(): ISyst("muSmear", "Muon energy smearing") {}

    virtual void Shift(double sigma,
                       Restorer& restore,
                       caf::SRProxy* sr,
                       double& weight) const override
    {
      restore.Add(sr->Elep_reco);

      // NB - the way this syst works there's no sense in doing -1 sigma
      sr->Elep_reco *= 1 + sigma*gRandom->Gaus(0, 0.2);
    }
  };
  EMuSmear kEMuSmear;

  // Let's require the reconstructed energy > 0 so we only include events we really managed to reconstruct
  Spectrum sSmear(loader, axRecoQEFormula, kHasCC0PiFinalState && kRecoQEFormulaEnergy>0,  SystShifts(&kEMuSmear, +1));
  
  
  
  // Smear muon angle with sigma of 30 degrees (pi/6 radians)
  class ThetaSmear: public ISyst
  {
  public:
    ThetaSmear(): ISyst("thetaSmear", "Muon angle smearing") {}

    virtual void Shift(double sigma,
                       Restorer& restore,
                       caf::SRProxy* sr,
                       double& weight) const override
    {
      restore.Add(sr->theta_reco);

      // Theta is in radians so we smear with a width of pi/6
      sr->theta_reco += sigma*gRandom->Gaus(0, TMath::Pi()/6.0);
    }
  };
  ThetaSmear kThetaSmear;
    
  Spectrum sThetaSmear(loader, axRecoQEFormula, kHasCC0PiFinalState && kRecoQEFormulaEnergy>0,  SystShifts(&kThetaSmear, +1));
 
  // Fill all the Spectrum objects from the loader
  loader.Go();

  //   Set to the same exposure as before
  const double pot = 1e20;

  // Convert and draw
  TCanvas *canvas = new TCanvas; // Make a canvas
  
  TH1D *hCV = sCV.ToTH1(pot, kAzure-7);
  // ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.html
  // Set the height a bit higher so we can see all the error bands
  hCV->GetYaxis()->SetRangeUser(0,hCV->GetMaximum()*1.3);
  hCV->Draw("E");
  
  // Scale up and down
  TH1D *hScaleUp = sScaleUp.ToTH1(pot, kOrange-2);
  TH1D *hScaleDn = sScaleDn.ToTH1(pot, kOrange-2,7);
  hScaleUp->Draw("HIST SAME");
  hScaleDn->Draw("HIST SAME");
  
  
  //Smear
  TH1D *hSmear = sSmear.ToTH1(pot, kOrange+7);
  hSmear->Draw("HIST SAME");
    
  // Angle Smear
  TH1D *hThetaSmear = sThetaSmear.ToTH1(pot, kAzure-9);
  hThetaSmear->Draw("HIST SAME");
    
  
  auto legend = new TLegend(0.65,0.65,0.9,0.9); // x and y coordinates of corners
  legend->AddEntry(hCV,"Central value","l");
  legend->AddEntry(hScaleUp,"Scale up","l");
  legend->AddEntry(hScaleDn,"Scale down","l");
  legend->AddEntry(hSmear,"Smear","l");
  legend->AddEntry(hThetaSmear,"#theta smear","l");
  legend->Draw();
  
  canvas->SaveAs("Systematics3.png"); // Save the result
  
  // Calculate ratios
  // Convert and draw
  TCanvas *canvasFrac = new TCanvas; // Make a canvas
  
  TH1D *fracScaleUp= MakeFractionalPlot(hScaleUp,hCV);
  TH1D *fracScaleDn= MakeFractionalPlot(hScaleDn,hCV);
  TH1D *fracSmear= MakeFractionalPlot(hSmear,hCV);
  TH1D *fracThetaSmear= MakeFractionalPlot(hThetaSmear,hCV);
  fracScaleUp->GetYaxis()->SetTitle("Fractional shift");
  fracScaleUp->GetYaxis()->SetRangeUser(-1.5,1.5);
  
  hCV->Sumw2(); //Make sure the errors behave
  TH1D *fracCV= MakeFractionalPlot(hCV,hCV); // Hope this is all zero! But it should have stat uncertainties on it
  
  fracScaleUp->Draw("HIST");
  fracScaleDn->Draw("HIST SAME");
  fracSmear->Draw("HIST SAME");
  fracThetaSmear->Draw("HIST SAME");
  fracCV->Draw("E SAME");
  
  auto legend2 = new TLegend(0.65,0.11,0.9,0.36); // x and y coordinates of corners
  legend2->AddEntry(fracCV,"Central value","l");
  legend2->AddEntry(fracScaleUp,"Scale up","l");
  legend2->AddEntry(fracScaleDn,"Scale down","l");
  legend2->AddEntry(fracSmear,"Smear","l");
  legend2->AddEntry(fracThetaSmear,"#theta smear","l");
  legend2->Draw();
  
  canvasFrac->SaveAs("SystematicsFractions.png");
  
}
