#include <iostream>
#include <cmath> 
#include <string.h>
#include <stdio.h>  
#include <math.h> 
#include <TCanvas.h>

#include "ReadRun.cc"
using namespace std;


void read_ps_Sr90_source() // main
{
	string path;
	
	// edit for your data folder
	path = "/home/djimeno/SHiP-Chris-djs/data/";

	switch (13) { //specify folder to run below, ALL bin files in this folder will be used
	case(0): {
		path += "10_beta_far_ps/";	//  PMT on channel 15, triggering on itself. Same for all beta_xx_xx (Cs source)
		break;
	}//
	case(1): {
		path += "11_beta_no_ps/"; 
		break;
	}//
	case(2): {
		path += "12_beta_near_ps/"; 
		break;
	}//
	case(3): {
		path += "13_beta_near_nops/"; 
		break;
	}//
	case(4): {
		path += "14_beta_no_nops/"; 
		break;
	}//
	case(5): {
		path += "15_beta_far_nops/"; 
		break;
	}//
	case(6): {
		path += "20_Sr_nops_drilled/";	//  PMT on channel 9, triggering on ch14,15 which corresponds to the coincidence triggering box. Same for Sr_xx_xx
		break;
	}//
	case(7): {
		path += "21_Sr_ps_drilled/"; 
		break;
	}//
	case(8): {
		path += "22_test_box2/"; 
		break;
	}//
	case(9): {
		path += "23_test_box2_reverse/"; 
		break;
	}//
	case(10): {
		path += "24_test_box1/"; 
		break;
	}//
	case(11): {
		path += "25_test_box1_reverse/"; 
		break;
	}//
	case(12): {
		path += "26_test_box1_repeat/"; 
		break;
	}//
	case(13): {
		path += "27_test_box1_reverse_repeat/"; 
		break;
	}//
	case(66): {
		path += "3_Test_Laser_13_10_2021/";	//  Test of DC triggering the laser(ch8), PMT on ch15
		break;
	}
	//DORAMAS
	default: {
		cout <<  "\n\n\n------------ You must specify a folder to run!!!!!!!" << endl;
		exit(0);
		break;
	}
	//
	}

	// read data
	ReadRun mymeas(path, true);
	
	//  DORAMAS: Add the channels to be shown in the plots
	vector<int> channels_plot = {9, 14, 15};	
	sort (channels_plot.begin(), channels_plot.end());	// It sorts the vector elements in ascending order

	for (int i = 0; i < channels_plot.size(); i++) {
		mymeas.plot_active_channels.push_back(channels_plot[i]);
	}
	//

	//apply baseline correction to ALL waveforms <- NEEDED but slow when not compiled
	//mymeas.SmoothAll(3); // smoothing of waveforms. Caution, will bias results!!
	mymeas.CorrectBaseline(0., 50.);	// use mean from 0 ns to 50 ns


	// print events above a threshold to identify interesting events
	mymeas.FractionEventsAboveThreshold(-4, false, false);

	////plotting
	
	//investigate individual waveforms
	//TCanvas* tstc = new TCanvas("tstc", "", 1600, 1000);
	//TH1F* histo = mymeas.Getwf(0, 0, 0);
	//histo->Draw();
	//tstc->BuildLegend(0.85, 0.70, .99, .95);


	// plot sums of all events per channel
	mymeas.PlotChannelSums(true);


	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 3.;	// lower integrationwidow in ns rel. to max
	float intwindowplus = 5.;	// upper integrationwidow in ns rel. to max
	float findmaxfrom = 50.;	// assume pulse after trigger arrives between here ...
	float findmaxto = 200.;		// ... and here (depends on trigger delay setting etc., for dark counts the signal is random so we look at the whole recorded time range)	
	

	// plot all charge spectrum of channels
	mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -5e1, 7.5e2, 250, 0, 0, 0);

	//int channel15_index = find(mymeas.active_channels.begin(), mymeas.active_channels.end(), 15) - mymeas.active_channels.begin(); 
	//cout << "\n index of channel 15: " << channel15_index << endl;
	//auto charge_spectrum = mymeas.ChargeSpectrum(static_cast<int>(channel15_index), intwindowminus, intwindowplus, findmaxfrom, findmaxto, -5e1, 7.5e2, 100);
	//auto cscanv = new TCanvas("charge spectrum");
	//charge_spectrum->GetYaxis()->SetTitle("#Entries");
	//charge_spectrum->GetXaxis()->SetTitle("integral in mV#timesns");
	//charge_spectrum->Draw();

	// plot waveforms of individual events
	vector<int> events_plot = {1, 15, 102, 123, 146, 155, 208, 222, 284, 299};	//  DORAMAS
	//plot range
	double ymin = -5;
	double ymax = 25;

	// plot waveforms for certain events with integration window
	//DORAMAS: for loop to create waveform_event pdf's
	for (int i = 0; i < events_plot.size(); i++) {
		mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, events_plot[i], ymin, ymax);
	}
	//
}
