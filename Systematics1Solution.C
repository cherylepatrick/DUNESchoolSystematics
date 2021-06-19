// To run this, type: cafe Systematics1Solution.C

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


// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Systematics1Solution()
{
  // Various input CAF samples
  
  const std::string FIRST_CAF = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_FHC_900.root"; //ND-LAr FHC - one file
  const std::string SECOND_CAF = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_FHC_901.root"; //ND-LAr FHC - one file
  const std::string ELEVEN_CAFS = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_FHC_90*.root"; //This wildcard gives 11 files!
  
  // Source of events - load them from the one of the sets of files
  SpectrumLoader lFirstCaf(FIRST_CAF);
  // ***** Add more loaders here (use the shortcut names defined above)

  // We want to plot a histogram with 40 bins, covering the range 0 to 10 GeV
  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  // ****** This is where you can change the number of bins

  // Define the label, binning, and contents that we want for our first histogram
  // The axis label can be whatever you like.
  // The binning needs to be a Binning object like the one we just made.
  // The Variable that we are plotting can be a single variable or function of variables from
  // the CAFs. See https://wiki.dunescience.org/wiki/CAFAna_Variables
  // We want to plot true energy
  const HistAxis axTrue("True neutrino energy (GeV)", binsEnergy, kTrueEnergy);

  // Now we have defined an axis and a variable we want to plot on it, let's decide which
  // events to plot. So here we are telling the loader to load a spectrum with our defined axis
  // Here, we are defining the selections or "cuts" using the variable names in CAFAna https://wiki.dunescience.org/wiki/CAFAna_Cuts
  
  
  // This cut selects neutrino-mode CC interactions
  const Cut kNuMuCC = kIsNumuCC && !kIsAntiNu; //The modes are defined at the top

  // Define the Spectrum
  Spectrum sFirstCaf(lFirstCaf, axTrue, kNuMuCC); // **** You'll be adding more Spectrum objectss

  
  // Fill all the Spectrum objects from the loader
  lFirstCaf.Go(); // **** You'll be adding more Loaders

  /* 
     Set to the same exposure as before
  */  
  const double pot = 1e20;
  // ******* This is where you change the POT

  // Convert and draw
  TCanvas *canvas = new TCanvas; // Make a canvas
  
  // Spectrum for CAF file 1
  TH1D *hFirstCaf = sFirstCaf.ToTH1(pot, kAzure-7);
  // ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.

  //  hFirstCaf->Print() // ***** Uncomment to see a print of the histogram values
  
  hFirstCaf->Draw("E"); // ****** Change this to show error bars

  
  auto legend = new TLegend(0.65,0.65,0.9,0.9); // x and y coordinates of corners
  legend->SetHeader("CAFs used","C"); // option "C" to center the header
  legend->AddEntry(hFirstCaf,"First CAF","l");
  legend->Draw();
  
  canvas->SaveAs("Systematics1.png"); // Save the result
}
