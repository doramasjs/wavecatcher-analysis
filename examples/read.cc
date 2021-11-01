#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>

#include "../ReadRun.cc"
using namespace std;

void read() // main
{
	int which = 10; //select meas

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
	}
	case(2): {
		path += "181_calib_vb40_tune8180_switch2_nopz/"; // one SiPM
		break;
	}//
	case(3): {
		path += "193_calib_vb41_tune8180_switch12_nopz/"; // one SiPM
		break;
	}//
	case(4): {
		path += "189_calib_vb40_tune8180_switch12_nopz/"; // one SiPM
		break;
	}//
	case(5): {
		path += "197_calib_vb41_tune8220_switch12_nopz/"; // one SiPM
		break;
	}//
	case(6): {
		path += "203_calib_vb41_tune8310_switch1_nopz/"; // one SiPM
		break;
	}//
	case(7): {
		path += "214_calib_vb41_tune8310_switch123_nopz/"; // one SiPM
		break;
	}//
	case(8): {
		path += "216_calib_vb41_tune8350_switch1_nopz/"; // one SiPM
		break;
	}//
	case(9): {
		path += "220_calib_vb41_tune8350_switch5_nopz/"; // one SiPM
		break;
	}//
	case(10): {
		path += "226_calib_vb58_tune8260_pcbd/"; // one SiPM
		break;
	}//
	case(11): {
		path += "7_calib_vb58_tune8700_pcbd/"; // one SiPM
		break;
	}//
	default: {
		path += "3_calib_vb56_tune8350/"; // ???
		break;
	}
	}

	// read data
	ReadRun mymeas(path);

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	int which_blc = 2;
	if (which_blc == 0) {
		mymeas.SmoothAll(5);
		mymeas.CorrectBaseline(95., 105.);	//
	}
	else if (which_blc == 1) {
		mymeas.CorrectBaselineMinSlopeRMS(20, true, 5, 352, 300, false);
	}
	else {
		mymeas.CorrectBaselineMin(20, false, 1., 372, 320, false);
	}

	////plotting

	//investigate individual waveforms
	//TCanvas* tstc = new TCanvas("tstc", "", 1600, 1000);
	//TH1F* histo = mymeas.Getwf(0, 0, 0);
	//histo->Draw();
	//tstc->BuildLegend(0.85, 0.70, .99, .95);

	// sums of all events per channel
	mymeas.PlotChannelSums(true);

	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 1.;	// lower integration window in ns rel. to max
	float intwindowplus = 1.;	// upper integration window in ns rel. to max
	float findmaxfrom = 110.;	// assume signal from laser arrives between here ...
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
		mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -20, 300, 320);
	}

	// plot waveforms for certain events with integration window
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 171, -2., 10.);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 1100, -2, 10);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 4, -2, 10);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 5, -2, 10);

	mymeas.PrintFFTWF(171, 0, .6, 64);
	mymeas.PrintFFTWF(1100, 0, .6, 64);
	mymeas.PrintFFTWF(4, 0, .6, 64);
	mymeas.PrintFFTWF(5, 0, .6, 64);
}