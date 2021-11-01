#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <TCanvas.h>

#include "../ReadRun.cc"
using namespace std;

void read_ps_sr90_source() // main
{
	int which = 3; //select meas

	string path;

	// edit for your fs
	path = "C:/SHiP/data/";

	switch (which) { //specify folder to run below, ALL bin files in this folder will be used
	case(0): {
		path += "18_Sr_ps/"; //
		break;
	}//
	case(1): {
		path += "16_Sr_nops/"; //
		break;
	}//
	case(2): {
		path += "20_Sr_nops_drilled/"; // first version of casing with drilled hole covered with tape
		break;
	}//
	case(3): {
		path += "21_Sr_ps_drilled/"; //
		break;
	}//
	case(4): {
		path += "22_test_box2/"; // The box used for previous measurements
		break;
	}//
	case(5): {
		path += "23_test_box2_reverse/"; //
		break;
	}//
	case(6): {
		path += "24_test_box1/"; // The other trigger box
		break;
	}//
	case(7): {
		path += "25_test_box1_reverse/"; //
		break;
	}//
	default: {
		path += "18_Sr_ps/"; // default
		break;
	}
	}

	// read data
	ReadRun mymeas(path, true);

	// only plot channel 14
	//int channel_to_plot = 14;
	//mymeas.plot_active_channels.push_back(channel_to_plot);
	//mymeas.plot_active_channels.push_back(1);

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	//mymeas.SmoothAll(3); // smoothing of waveforms. Caution, will bias results!!
	mymeas.CorrectBaseline(0., 50.);	// use mean from 0 ns to 50 ns

	// print events above a threshold to identify interesting events
	mymeas.FractionEventsAboveThreshold(4, true, true, 100, 150);

	mymeas.FractionEventsAboveThreshold(-4, false, false, 100, 150);
	mymeas.FractionEventsAboveThreshold(-6, false, false, 100, 150);
	mymeas.FractionEventsAboveThreshold(-8, false, false, 100, 150);
	mymeas.FractionEventsAboveThreshold(-12, false, false, 100, 150);

	////plotting

	//investigate individual waveforms
	//TCanvas* tstc = new TCanvas("tstc", "", 1600, 1000);
	//TH1F* histo = mymeas.Getwf(0, 0, 0);
	//histo->Draw();
	//tstc->BuildLegend(0.85, 0.70, .99, .95);

	// plot sums of all events per channel
	mymeas.PlotChannelSums(true);

	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 3.;	// lower integration window in ns rel. to max
	float intwindowplus = 5.;	// upper integration window in ns rel. to max
	float findmaxfrom = 100.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 150.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)

	// plot all charge spectrum of channels
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -1e2, 2.e3, 500, 0, 0, 0);

	mymeas.PrintChargeSpectrumPMT(0, 0, findmaxfrom, findmaxto, -2e1, 1.8e2, 202, 4.);

	// timing of maximum
	mymeas.PrintTimeDist(findmaxfrom, findmaxto, findmaxfrom - 5, findmaxto + 5, 60);

	// plot waveforms of individual events
	int event1 = 68;
	int event2 = 79;
	int event3 = 269;
	int event4 = 270;
	//plot range
	double ymin = -5;
	double ymax = 25;

	// plot waveforms for certain events with integration window
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event1, ymin, ymax);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event2, ymin, ymax);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event3, ymin, ymax);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, event4, ymin, ymax);
}