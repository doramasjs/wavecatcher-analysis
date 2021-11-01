#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>

#include "../ReadRun.cc"
using namespace std;

void read_cosmics() // main
{
	int which = 7; //select meas

	// better create separate file just for DC measurements
	bool isDC = false;

	string path;

	// edit for your fs
	path = "C:/SHiP/data/";

	switch (which) { //specify folder to run below, ALL bin files in this folder will be used
	case(0): {
		path += "176_calib_vb41_tune8180_switch1_nopz/"; // one SiPM
		break;
	}//
	case(1): {
		path += "182_calib_vb41_tune8180_switch2_nopz/"; // one SiPM
		break;
	}//
	case(2): {
		path += "13_1904_night/"; // old before fw upgrade
		break;
	}//
	case(3): {
		path += "112_2106_day/"; // new w/o non-active PMTs
		break;
	}//
	case(4): {
		path += "229_cosmics_vb58_newsupply/"; // different PSU
		break;
	}//
	case(5): {
		path += "119_2806_day2_mod/"; // different PSU
		break;
	}//
	case(6): {
		path += "316_cosmics_pcbj_0109_day_nopz/"; // different PSU
		break;
	}//
	case(7): {
		path += "1_Test_12_10_2021/"; // different PSU
		break;
	}//
	default: {
		path += "35_0505_night/"; // new w/ two non-active PMTs
		break;
	}
	}

	// read data
	ReadRun mymeas(path);

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	int which_blc = 0;
	if (which_blc == 0) {
		//mymeas.SmoothAll(1., true);
		mymeas.CorrectBaseline(0., 50.);	//
	}
	else if (which_blc == 1) {
		mymeas.CorrectBaselineMinSlopeRMS(20, true, 5, 352, 300, false);
	}
	else {
		mymeas.CorrectBaselineMin(20, true, 1., 372, 0, false, 8);
	}

	//mymeas.FractionEventsAboveThreshold(4, 7, true, true);

	////plotting

	//investigate individual waveforms
	//TCanvas* tstc = new TCanvas("tstc", "", 1600, 1000);
	//TH1F* histo = mymeas.Getwf(0, 0, 0);
	//histo->Draw();
	//tstc->BuildLegend(0.85, 0.70, .99, .95);

	// sums of all events per channel
	mymeas.PlotChannelSums(true);

	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 10.;	// lower integration window in ns rel. to max
	float intwindowplus = 15.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume signal from laser arrives between here ...
	float findmaxto = 125.;		// ... and here (depends on trigger delay setting)

	if (isDC) {
		findmaxfrom = 10 + intwindowminus;
		findmaxto = 280. - intwindowplus;
	}

	// plot all channels
	if (isDC) {
		mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -5, 100, 100);
	}
	else {
		mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -50, 5000, 400);
	}

	int event1 = 214;
	int event2 = 1421;
	int event3 = 1500;

	// plot waveforms for certain events with integration window
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event1, -15, 150.);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event2, -15, 150);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event3, -15, 150);
	//mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 500, -15, 15);

	bool printfft = false;
	if (printfft) {
		mymeas.PrintFFTWF(event1, 0, .6, 64);
		mymeas.PrintFFTWF(event2, 0, .6, 64);
		mymeas.PrintFFTWF(event3, 0, .6, 64);
	}
}