//This is the code. You just need to adjust the path to the data in "read_ps.cc" and then load it
//In the shell -> go to directory where the code is, open root, and type ".x read_ps.cc"
//You have to close root (".q") and restart it if you want to re-run the code
//You can select the channels you want with mymeas.plot_active_channels.push_back(channel_to_plot); to add them to the list
//   if you uncomment it it plots all channels as before

//DORAMAS -> means I have made some changes
//Doramas -> means I have only added a comment

#include "ReadRun.h"

ClassImp(ReadRun)

ReadRun::ReadRun(int bla) {
	cout << "\ninit" << endl;
}

ReadRun::ReadRun(string path, bool changesignofPMTs) {
	bool out = false; //experimental, not working

	// Wavecatcher hardware/software properties
	SP = 0.3125;					// ns per bin
	pe = 47.46;					//mV*ns ????
	coef = 2.5 / (4096 * 10);	//?????
	binNumber = 1024;				//default: 1024, hard coded later on so it can be removed
	const int nChannelsWC = 64;			//max number of channels default: 32

	rundata = new TClonesArray("TH1F", 1e7); //raw data will be stored here as TH1F
	rundata->BypassStreamer(kFALSE);  //Doramas: I don't know why is it used, but it's better to use when working with TClonesArray
	TClonesArray& testrundata = *rundata;

	// verbosity
	bool debug_header = 0;
	bool debug_data = 0;

	unsigned short output_channel;
	unsigned int output_event;
	unsigned long long int output_tdc;
	unsigned short output_nbchannels;

	amplValuessum = new double* [nChannelsWC]; //sum of all wf for each channel
	for (int i = 0; i < nChannelsWC; i++) {//init
		amplValuessum[i] = new double[binNumber];
		for (int k = 0; k < binNumber; k++) amplValuessum[i][k] = 0.;
	}

	maxSumBin = new int[nChannelsWC];

	//Start reading the raw data from .bin files.
	stringstream inFileList = list_files(path.c_str(), ".bin"); //all *.bin* files in folder path
	int nitem = 1;
	string fileName;
	string outfileName;
	int file_counter = 0;
	int currentPrint = -1;
	int wfcounter = 0;
	int event_counter = 0;

	while (inFileList >> fileName) {
		// file loop

		fileName = path + fileName;
		outfileName = path + "mod/" + fileName;
		ifstream input_file(fileName.c_str(), std::ios::binary | std::ios::in);
		ofstream output_file(outfileName.c_str(), std::ios::binary | std::ios::out);

		bool has_measurement = false;

		if (!input_file.is_open()) {
			printf("*** failed to open '%s'\n", fileName.c_str());
			continue;
		}
		printf("+++ reading '%s' ...\n", fileName.c_str());

		// Header
		string header_line;
		// HEADER 1 //
		//
		// "=== DATA FILE SAVED WITH SOFTWARE VERSION: V?.??.? ==="
		//
		getline(input_file, header_line, '\n');
		if (out) output_file.write(header_line.c_str(), sizeof(char) * header_line.size());

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t header_version_first = header_line.find_last_of('V');
		size_t header_version_last = header_line.find_first_of(' ', header_version_first);
		string software_version = header_line.substr(header_version_first, header_version_last - header_version_first);
		if (debug_header) printf("    |- data version = '%s'\n", software_version.data());

		//if (software_version == "V2.9.13")
		//	;
		//else if (software_version == "V2.9.15")
		//	;
		//else if (debug_header) printf("*** unsupported data version\n");

		// HEADER 2 //
		// "=== WAVECATCHER SYSTEM OF TYPE ?? WITH ?? CHANNELS AND GAIN: ??? ==="
		getline(input_file, header_line, '\n');
		if (out) output_file.write(header_line.c_str(), sizeof(char) * header_line.size());

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 3 //
		// === Rate coincidence masks ... === Posttrig in ns for SamBlock ... ===
		getline(input_file, header_line, '\n');
		if (out) output_file.write(header_line.c_str(), sizeof(char) * header_line.size());

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 4 //
		// V2.9.13: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1
		// V2.9.15: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1 == MEASUREMENTS: 0 ===
		getline(input_file, header_line, '\n');
		if (out) output_file.write(header_line.c_str(), sizeof(char) * header_line.size());

		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t nsamples_first = 1 + header_line.find_last_of('[');
		size_t nsamples_last = header_line.find_first_of(']', nsamples_first);
		string nsamples_str = header_line.substr(nsamples_first, nsamples_last - nsamples_first);

		int nsamples = atoi(nsamples_str.data());
		if (debug_header) printf("    |- data sample  = %d\n", nsamples);

		size_t nchannels_first = 10 + header_line.find("ACQUIRED: ", nsamples_first);
		size_t nchannels_last = header_line.find_first_of(' ', nchannels_first);
		string nchannels_str = header_line.substr(nchannels_first, nchannels_last - nchannels_first);

		nchannels = atoi(nchannels_str.data());
		if (debug_header) printf("    |- nchannels    = %d\n", nchannels);

		if (software_version == "V2.9.15" || software_version == "V2.9.16" || software_version == "V2.10.1") {
			size_t has_measurement_first = 14 + header_line.find("MEASUREMENTS: ", nsamples_first);
			size_t has_measurement_last = header_line.find_first_of(' ', has_measurement_first);
			string has_measurement_str = header_line.substr(has_measurement_first, has_measurement_last - has_measurement_first);
			has_measurement = atoi(has_measurement_str.data());
		}
		else {
			//if (software_version == "V2.9.13") {
				// V2.9.13 has always measurement stored
				// (everything is set to 0 when disabled!)
			has_measurement = true;
		}

		if (debug_header) printf("    `- measurement  = %d\n", has_measurement);

		// end of header reader

		event_data an_event;

		while (input_file.read((char*)(&an_event), sizeof(an_event))) {
			//file loop

			if (out) output_file.write((char*)(&an_event), sizeof(an_event));

			if (debug_data) printf("%03d has %d channels\n", an_event.EventNumber, an_event.nchannelstored);

			output_event = an_event.EventNumber;
			output_tdc = an_event.TDCsamIndex;
			output_nbchannels = an_event.nchannelstored;

			if (debug_data && output_event % 200 == 0) printf("EventNr: %d, nCh: %d\n", output_event, output_nbchannels);

			for (int ch = 0; ch < output_nbchannels; ++ch) {
				//
				channel_data_with_measurement a_channel_data;
				channel_data_without_measurement a_channel_data_without_measurement;

				if (has_measurement) {
					// read with 'channel_data_with_measurement' struct
					input_file.read((char*)(&a_channel_data), sizeof(channel_data_with_measurement));
				}
				else {
					// read with 'channel_data_without_measurement' struct
					input_file.read((char*)(&a_channel_data_without_measurement), sizeof(channel_data_without_measurement));

					// copy the content into 'channel_data_with_measurement' struct
					a_channel_data.channel = a_channel_data_without_measurement.channel;
					a_channel_data.EventIDsamIndex = a_channel_data_without_measurement.EventIDsamIndex;
					a_channel_data.FirstCellToPlotsamIndex = a_channel_data_without_measurement.FirstCellToPlotsamIndex;
					memcpy(a_channel_data.waveform, a_channel_data_without_measurement.waveform, 1024 * sizeof(short));
				}

				output_channel = a_channel_data.channel;
				if (debug_data) printf("- reading channel %d\n", output_channel);

				if (event_counter == 0) active_channels.push_back(static_cast<int>(output_channel));

				TString name(Form("channel_%02d, event %05d ", output_channel, an_event.EventNumber));
				TString title(Form("Channel %d, event %d data", output_channel, an_event.EventNumber));
				TH1F* hCh = (TH1F*)testrundata.ConstructedAt(wfcounter);
				hCh->SetName(name.Data());
				hCh->SetTitle(title.Data());
				hCh->SetBins(binNumber, -0.5 * SP, 1023.5 * SP);

				float val = 0.;
				for (int s = 0; s < binNumber; ++s) {
					val = a_channel_data.waveform[s] * coef * 1000.;
					if (changesignofPMTs && output_channel > 8) val *= -1.;
					hCh->SetBinContent(s + 1, val);
					//hCh->SetBinError(s, 0.5); //The error of each value in each bin is set to 0.5 mV -> Why??

					if (out) {
						a_channel_data.waveform[s] = 0.;
						if (s == 300) a_channel_data.waveform[s] = 1.;
					}

					// channel sums
					amplValuessum[ch][s] += static_cast<double>(val);
				}

				//hCh->SetLineColor(ch + 1); // gets a bit too colorful
				//hCh->SetMarkerColor(ch + 1);
				if (out) {
					if (has_measurement) {
						// read with 'channel_data_with_measurement' struct
						output_file.write((char*)(&a_channel_data), sizeof(channel_data_with_measurement));
					}
					else {
						output_file.write((char*)(&a_channel_data_without_measurement), sizeof(channel_data_without_measurement));
					}
				}

				wfcounter++;
			} // for ch

			eventnr_storage.push_back(output_event);	//  Adds the current event number(the one from the WaveCatcher) to the storage vector
			event_counter++;
		} // while an_event

		input_file.close();
		if (out) output_file.close();
		file_counter++;
	} // for file_id

	// in case there are empty channels, nchannels is the number of channels which contain
	nchannels = output_nbchannels;

	// get bins where the sum spectrum has its maximum for runs with fixed trigger delay and fixed integration window relative to the max of the sum spectrum (not working for DC measurement)
	for (int ch = 0; ch < nchannels; ch++) {
		double max = 0.;
		for (int i = 0; i < binNumber; i++) {
			if (amplValuessum[ch][i] > max) {
				max = amplValuessum[ch][i];
				maxSumBin[ch] = i;
			}
		}
	}

	nevents = event_counter;
	nwf = wfcounter;
}

ReadRun::~ReadRun() {
	// Destructor
	//rundata->Clear();
	//delete[] maxSumBin;
	//delete baseline_correction_result;
	plot_active_channels.clear();
	cout << "deleting nothing currently..." << endl;
}

// plot sums of all waveforms for each channel

void ReadRun::PlotChannelSums(bool doaverage) {
	// doaverage: if true it will plot the running average +/- 4 bins

	double* xv = getx();
	TMultiGraph* mgsums = new TMultiGraph();
	mgsums->SetTitle("channel sums; t [ns]; amplitude [arb.]");

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			double* yv = amplValuessum[i];
			if (doaverage) SmoothArray(yv, binNumber, 4);
			TGraph* gr = new TGraph(binNumber, xv, yv);
			delete[] yv;

			TString name(Form("channel_%02d", active_channels[i]));
			TString title(Form("Channel %d", active_channels[i]));
			gr->SetName(name.Data());
			gr->SetTitle(title.Data());
			gr->SetLineColor(i + 1);
			gr->SetMarkerColor(i + 1);
			mgsums->Add(gr);
		}
	}
	delete[] xv;

	TCanvas* sumc = new TCanvas("Sums", "", 1600, 1000);
	mgsums->Draw("APL");
	mgsums->GetYaxis()->SetRangeUser(-1e4, 10e4);
	sumc->BuildLegend(0.85, 0.70, .99, .95);
	sumc->SaveAs("channelsums.pdf");
}

// averaging all waveforms (for testing)

void ReadRun::SmoothAll(double sigma, bool doconv) { //deprecated since it can be done with baseline correction??
	// just for testing, not very efficient
	cout << "\nsmoothing wfs";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = ((TH1F*)rundata->At(j));
		double* yvals = gety(his);
		SmoothArray(yvals, binNumber, sigma, doconv);
		for (int i = 1; i < his->GetNbinsX(); i++) his->SetBinContent(i, yvals[i]);
		delete[] yvals;
		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

// baseline correction

void ReadRun::CorrectBaseline(float tCut, float tCutEnd) {
	// corrects the baseline (DC offset) of all waveformsf
	// tCut: time denoting the end or the beginning (if tCutEnd is set) of integration window
	// for no bl correction set tCut < 0

	int iCut, iCutEnd;
	float corr = 0;
	printf("\nBaseline correction (%d waveforms) :: ", nwf);

	for (int j = 0; j < nwf; j++) {
		TH1F* his = ((TH1F*)rundata->At(j));
		iCut = his->GetXaxis()->FindBin(tCut);

		// start of temp for strange PMT signals in cosmics setup
		//int currchannel = j - nchannels * floor(j / nchannels);
		//if (currchannel > 7 && his->GetMaximum() < 5) cout << "\nevent:\t" << 1 + floor(j / nchannels) << "\tchannel:\t" << currchannel;
		// end of temp for strange PMT signals in cosmics setup

		if (tCutEnd <= 0) { //
			corr = his->Integral(1, iCut) / (iCut - 1);
		}
		else {
			iCutEnd = his->GetXaxis()->FindBin(tCutEnd);
			corr = his->Integral(iCut, iCutEnd) / (iCutEnd - iCut);
		}

		// write corrected values to histograms
		if (tCut >= 0) {
			for (int i = 1; i < his->GetNbinsX(); i++) his->SetBinContent(i, his->GetBinContent(i) - corr);
		}

		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[j].push_back(corr);
		baseline_correction_result[j].push_back(0);
		baseline_correction_result[j].push_back(tCut);
		baseline_correction_result[j].push_back(tCutEnd);

		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

void ReadRun::CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool search_min, bool convolution, int skip_channel) {
	// corrects the baseline (DC offset) of all waveforms
	// determines the region of nIntegrationWindow bins where the root mean square of the slope of the smoothed waveform reaches its minimum
	// experimental, use with care and check output,
	// nIntegrationWindow: number of bins of integration window

	const int binNumberSlope = binNumber - 1;
	double* slope = new double[binNumberSlope];
	skip_channel += 1;

	if (start_at > max_bin_for_baseline - nIntegrationWindow) start_at = 0;

	int min_distance_from_max = 25 + nIntegrationWindow;
	float corr = 0;
	float minchange = 1.e9;
	float minsum = 0;
	float minsum0 = 0;
	float minsqsum = 0;
	int iintwindowstart = 0;

	float sum = 0.;
	float sum0 = 0.;
	float sqsum = 0.;
	float change = 0.;
	float sign = 1.;

	int imax = 0;
	int search_before = 0;

	printf("Baseline correction (%d waveforms) :: ", nwf);

	for (int j = 0; j < nwf; j++) {
		corr = 0;
		minchange = 1.e9;
		iintwindowstart = 0;

		sum = 0.;
		sum0 = 0.;
		sqsum = 0.;
		change = 0.;
		sign = 1.;

		imax = 0;
		search_before = 0;

		if (j == 0 || j != skip_channel - 1 || j % skip_channel != 0) { //eventnr * nchannels + i
			TH1F* his = ((TH1F*)rundata->At(j));
			double* yvals = gety(his); //find faster way
			SmoothArray(yvals, binNumber, sigma, convolution); // smoothing important to suppress variations in slope due to noise so the method is more sensitve to excluding peaks

			//calculate slope
			for (int i = 0; i < binNumberSlope; i++) slope[i] = yvals[i + 1] - yvals[i];

			if (max_bin_for_baseline != 0 && max_bin_for_baseline > nIntegrationWindow) {
				search_before = max_bin_for_baseline - nIntegrationWindow - 1;
			}
			else {
				imax = his->GetMaximumBin();
				search_before = imax - min_distance_from_max;
			}

			for (int i = start_at; i < search_before; i += 1) { // currently in steps of 3 bins (~1 ns) to make it faster
				sum = 0.;
				sum0 = 0.;
				sqsum = 0.;
				change = 0.;
				sign = 1.;
				for (int k = i; k < nIntegrationWindow + i; k++) {
					sum += slope[k];
					sum0 += yvals[k] / 100; // completely random choice
					sqsum += (slope[k] * slope[k]);
				}
				if (sum0 < 0) sign = -1.;

				if (search_min) change = sqsum + sum * sum + sum0 * sign;
				else change = sqsum + sum * sum;

				if (change < minchange) {
					minchange = change;
					iintwindowstart = i;
					minsum = sum * sum;
					minsqsum = sqsum;
					minsum0 = sum0 * sign;
				}
			}

			corr = 0.;
			if (!doaverage) {
				corr = his->Integral(iintwindowstart, iintwindowstart + nIntegrationWindow) / static_cast<float>(nIntegrationWindow);
			}
			else {
				for (int i = iintwindowstart; i < iintwindowstart + nIntegrationWindow; i++) corr += yvals[i];
				corr /= static_cast<float>(nIntegrationWindow);
			}

			for (int i = 0; i < binNumber; i++) {
				if (!doaverage) his->SetBinContent(i, his->GetBinContent(i) - corr);
				else his->SetBinContent(i, yvals[i] - corr);
			}
			delete[] yvals; //delete slow
		}

		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[j].push_back(corr);
		baseline_correction_result[j].push_back(minchange);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart) * SP);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart + nIntegrationWindow) * SP);
		baseline_correction_result[j].push_back(minsum);
		baseline_correction_result[j].push_back(minsum0);
		baseline_correction_result[j].push_back(minsqsum);

		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
	delete[] slope;

	TCanvas* sumc = new TCanvas("sumc", "sumc", 1600, 1000);
	TH1F* hiss = new TH1F("sum", "sum", 1e4, -2, 2);
	TH1F* hiss0 = new TH1F("sum0", "sum0", 1e4, -2, 2);
	TH1F* hisssq = new TH1F("sqsum", "sqsum", 1e4, -2, 2);
	for (int i = 0; i < nwf; i++) {
		hiss->Fill(baseline_correction_result[i][4]);
		hiss0->Fill(baseline_correction_result[i][5]);
		hisssq->Fill(baseline_correction_result[i][6]);
	}
	hiss0->SetLineColor(2); hisssq->SetLineColor(1);
	hiss->Draw(); hiss0->Draw("same"); hisssq->Draw("same");
}

void ReadRun::CorrectBaselineMin(int nIntegrationWindow, bool doaverage, double sigma, int max_bin_for_baseline, int start_at, bool convolution, int skip_channel) {
	// corrects the baseline (DC offset) of all waveforms
	// experimental - uses min(mean(nIntegrationWindow)) in range (start_at, max_bin_for_baseline)

	int binNumberSlope = binNumber - 1;
	skip_channel += 1;

	if (start_at > max_bin_for_baseline - nIntegrationWindow) start_at = 0;

	int min_distance_from_max = 25 + nIntegrationWindow;

	float corr = 0;
	float minchange = 1.e9;
	int iintwindowstart = 0;
	float sum0 = 0.;
	int imax = 0;
	int search_before = 0;

	printf("Baseline correction (%d waveforms) :: ", nwf);

	for (int j = 0; j < nwf; j++) {
		minchange = 1.e9;
		iintwindowstart = 0;

		if (j == 0 || j != skip_channel - 1 || j % skip_channel != 0) { //eventnr * nchannels + i
			TH1F* his = ((TH1F*)rundata->At(j));
			double* yvals = gety(his); //find faster way
			SmoothArray(yvals, binNumber, sigma, convolution); // smoothing

			if (max_bin_for_baseline != 0 && max_bin_for_baseline > nIntegrationWindow) {
				search_before = max_bin_for_baseline - nIntegrationWindow - 1;
			}
			else {
				imax = his->GetMaximumBin();
				search_before = imax - min_distance_from_max;
			}

			for (int i = start_at; i < search_before; i++) { // can also be done in coarser steps
				sum0 = 0.;
				for (int k = i; k < nIntegrationWindow + i; k += 2) { // can also be done in coarser steps
					sum0 += yvals[k];
				}

				if (sum0 < minchange) {
					minchange = sum0;
					iintwindowstart = i;
				}
			}

			corr = 0.;
			if (!doaverage) {
				corr = his->Integral(iintwindowstart, iintwindowstart + nIntegrationWindow) / static_cast<float>(nIntegrationWindow);
			}
			else {
				for (int i = iintwindowstart; i < iintwindowstart + nIntegrationWindow; i++) corr += yvals[i];
				corr /= static_cast<float>(nIntegrationWindow);
			}

			for (int i = 0; i < binNumber; i++) {
				if (!doaverage) his->SetBinContent(i, his->GetBinContent(i) - corr);
				else his->SetBinContent(i, yvals[i] - corr);
			}
			delete[] yvals; //delete slow
		}

		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[j].push_back(corr);
		baseline_correction_result[j].push_back(minchange);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart) * SP);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart + nIntegrationWindow) * SP);

		if ((j + 1) % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j + 1) / static_cast<float>(nwf) << "% -" << flush;
	}
}

void ReadRun::FractionEventsAboveThreshold(float threshold, bool max, bool greater, double from, double to) {
	// find events with max/min above/below a certain threshold
	// threshold -> in mV
	// max -> true uses max, false uses min
	// greater -> true looks for events with max/min>threshold, false looks for events with max/min<threshold
	int occurences = 0;
	int occurences2ch = 0;
	int o2ch = 0;
	int currchannel = 0;
	int currevent = 0;
	int lastevent = 0;
	if (plot_active_channels.empty()) plot_active_channels = active_channels;
	vector<int> counter_abovethr(plot_active_channels.size());	// DORAMAS: It stores a counter of events above threshold for each channel that will be plotted

	cout << "\n\n ------> ";
	if (max) cout << "max";
	else cout << "min";

	if (greater) cout << " > ";
	else cout << " < ";
	cout << threshold << " mV:\n";

	for (int j = 0; j < nwf; j++) {
		auto his = (TH1F*)((TH1F*)rundata->At(j))->Clone();

		// set range (changes all histograms..)
		if (from >= 0 && to > 0) {
			his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));
		}

		currchannel = j - nchannels * floor(j / nchannels);
		if (find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[currchannel]) != plot_active_channels.end()) {
			if ((max && greater && his->GetMaximum() > threshold) || (max && !greater && his->GetMaximum() < threshold) || (!max && greater && his->GetMinimum() > threshold) || (!max && !greater && his->GetMinimum() < threshold)) {
				currevent = eventnr_storage[floor(j / nchannels)];
				/*cout << "\nevent:\t" << currevent << "\tchannel:\t" << active_channels[currchannel];*/

				// We must use 'distance' to make sure the position in 'counter_above' matches with the corresponding channel's position at 'plot_active_channels'
				counter_abovethr[distance(plot_active_channels.begin(), find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[currchannel]))] += 1;
				// This is to detect events w/ at least two channels above threshold
				if (lastevent != currevent) occurences += 1;
				if (lastevent == currevent && o2ch != occurences) {
					occurences2ch += 1;
					o2ch = occurences;
				}
				lastevent = currevent;
				//
			}
		}
	}

	//  Loop to show the fraction of events above threshold for each channel that will be plotted
	for (int i = 0; i < plot_active_channels.size(); i++) {
		cout << "\nfraction of events in channel " << plot_active_channels[i] << " above threshold: " << 100. * static_cast<float>(counter_abovethr[i]) / static_cast<float>(nevents) << "%\n";
	}
	//
	cout << "\nfraction of events w/ at least 2 channels above threshold: " << 100. * static_cast<float>(occurences2ch) / static_cast<float>(nevents) << "%\n";
	cout << "\tfor a total of " << nevents << " events\n" << endl;
}

// functions for charge spectrum

int* ReadRun::GetIntWindow(TH1F* his, float windowlow, float windowhi, float start, float end, int channel) {
	// find maximum in range (start, end) and return bin numbers for [0] the max, [1] t_max - windowlow, and [2] t_max + windowhi
	// if (start < 0 || end < 0) doesn't return max and integration window is fixed t(max(sum_spectrum[channel])) +/- hi/lo
	// if (windowlow == start && windowwhi == end) doesn't return max and sets fixed integration window from start until end for all channels

	int istart, iend;
	int* foundindices = new int[3];//
	foundindices[0] = 0;

	if (start < 0 || end < 0) {									// fixed integration window relative to maximum of sum spectrum for each channel
		foundindices[1] = his->GetXaxis()->FindBin(static_cast<float>(maxSumBin[channel]) * SP - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(static_cast<float>(maxSumBin[channel]) * SP + windowhi);
	}
	else if (windowlow == start && windowhi == end) {				// fixed integration window for all channels
		foundindices[1] = his->GetXaxis()->FindBin(windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(windowhi);
	}
	else {															// fixed integration window relative to maximum of each individual waveform
		istart = his->GetXaxis()->FindBin(start);
		iend = his->GetXaxis()->FindBin(end);

		if (istart<1 || iend>his->GetNbinsX()) {
			cout << "\nError: start or end out of range" << endl;
			return 0;
		}

		float max = -1e3;
		float val = 0;
		for (int i = istart; i < iend; i++) {
			val = his->GetBinContent(i);
			if (val > max) {
				max = val;
				foundindices[0] = i;
			}
		}

		foundindices[1] = his->GetXaxis()->FindBin(static_cast<float>(foundindices[0]) * SP - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(static_cast<float>(foundindices[0]) * SP + windowhi);
	}
	return foundindices;
}

void ReadRun::PrintChargeSpectrumWF(float windowlow, float windowhi, float start, float end, int eventnr, float ymin, float ymax) {
	// plot waveforms of all channels for a given event number eventnr and add the determined integration windwos to the plot
	gStyle->SetOptStat(0);

	TString name(Form("waveforms_event__%05d", eventnr));
	TCanvas* intwinc = new TCanvas(name.Data(), name.Data(), 1600, 1000);

	if (plot_active_channels.empty()) intwinc->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	else intwinc->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);

	int current_canvas = 0;
	//  DORAMAS: We need to find the position where the event has been stored. This way we will go to the correct WC event when using rundata->At()
	eventnr = distance(eventnr_storage.begin(), find(eventnr_storage.begin(), eventnr_storage.end(), eventnr));

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;

			TH1F* his;
			his = ((TH1F*)rundata->At(eventnr * nchannels + i));
			int* windowind = GetIntWindow(his, windowlow, windowhi, start, end, i);
			// create lines to indicate the integration window
			TLine* low = new TLine(his->GetXaxis()->GetBinCenter(windowind[1]), -5, his->GetXaxis()->GetBinCenter(windowind[1]), 10);
			low->SetLineColor(2);
			TLine* hi = new TLine(his->GetXaxis()->GetBinCenter(windowind[2]), -2, his->GetXaxis()->GetBinCenter(windowind[2]), 3);
			hi->SetLineColor(2);
			TLine* zero = new TLine(0, 0, 320, 0); // draw line at x=0 to check if baseline correction worked
			zero->SetLineColor(1);
			delete[] windowind;

			TLine* baselinel = new TLine(baseline_correction_result[eventnr * nchannels + i][2], -1, baseline_correction_result[eventnr * nchannels + i][2], 1);
			baselinel->SetLineColor(6);
			baselinel->SetLineWidth(2);
			TLine* baselineh = new TLine(baseline_correction_result[eventnr * nchannels + i][3], -1, baseline_correction_result[eventnr * nchannels + i][3], 1);
			baselineh->SetLineColor(6);
			baselineh->SetLineWidth(2);
			TLine* baseline = new TLine(baseline_correction_result[eventnr * nchannels + i][2], 0, baseline_correction_result[eventnr * nchannels + i][3], 0);
			baseline->SetLineColor(6);

			// draw to canvas
			intwinc->cd(current_canvas);
			his->Draw();
			if (ymin != 0. && ymax != 0.) his->GetYaxis()->SetRangeUser(ymin, ymax); //for better comparison fix y range
			low->Draw("same");
			hi->Draw("same");
			zero->Draw("same");
			baselinel->Draw("same");
			baselineh->Draw("same");
			baseline->Draw("same");
		}
	}
	intwinc->Update();

	stringstream namess;
	namess << name.Data() << ".pdf";
	intwinc->SaveAs(namess.str().c_str());
}

TH1F* ReadRun::ChargeSpectrum(int channel_index, float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins) {
	// integrate all pulses in range (start, end) from t_max - windowlow to t_max + windowhi for a given channel and return the charge histogram with x range (rangestart, rangeend) and the number of bins nbins

	TString name(Form("channel__%02d", active_channels[channel_index]));
	TH1F* h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);

	for (int j = 0; j < nevents; j++) {
		TH1F* his = ((TH1F*)rundata->At(j * nchannels + channel_index));
		int* windowind = GetIntWindow(his, windowlow, windowhi, start, end, channel_index);	// find integration window
		h1->Fill(his->Integral(windowind[1], windowind[2]));					// fill charge spectrum
		delete[] windowind;
	}

	return h1;
}

void ReadRun::PrintChargeSpectrum(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend, int max_channel_nr_to_fit) {
	// print ReadRun::ChargeSpectrum for all channels optimized for SiPM signals
	// TODO: INCLUDE BETTER WAY TO READ STARTING VALUES OF FIT PARAMETERS FROM FILE STORED IN DATA DIRECTORY FOR EACH MEASUREMENT

	gStyle->SetOptStat("ne");
	gStyle->SetOptFit(1111);

	if (fitrangestart == 0.) fitrangestart = rangestart;
	if (fitrangeend == 0.) fitrangeend = rangeend;

	TCanvas* chargec = new TCanvas("charge spectra", "charge spectra", 1600, 1000);

	if (plot_active_channels.empty()) chargec->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	else chargec->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
	cout << "\n\nThere is data recorded in " << active_channels.size() << " channels \n\n\n";
	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;

			TH1F* his;
			his = ChargeSpectrum(i, windowlow, windowhi, start, end, rangestart, rangeend, nbins);
			chargec->cd(current_canvas);

			//Fitf fitf;
			//TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 7);
			//f->SetLineColor(3);
			//f->SetParName(0, "N0");					f->SetParameter(0, his->Integral()/3.);
			//f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
			//f->SetParName(2, "#lambda");			f->SetParameter(2, .15); //0.2 or 3
			//f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 3.2);//3.6
			//f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, .12);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
			//f->SetParName(5, "Gain");				f->SetParameter(5, 10.);	//f->FixParameter(5, 10.);
			//f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);//7.);									//f->FixParameter(6, 0.);

			Fitf_biased fitf_biased;
			TF1* f = new TF1("fitf", fitf_biased, fitrangestart, fitrangeend, 9);
			f->SetLineColor(3);

			//+-1ns 41V nopz 1 SiPM
			//f->SetParName(0, "N0");					f->SetParameter(0, his->Integral());
			//f->SetParName(1, "#mu");				f->SetParameter(1, 3.3);// 1.6);
			//f->SetParName(2, "#lambda");			f->SetParameter(2, .01); //0.2 or 3
			//f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 4.3);//3.6
			//f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 1.5);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
			//f->SetParName(5, "Gain");				f->SetParameter(5, 21.);	//f->FixParameter(5, 10.);
			//f->SetParName(6, "Pedestal");			f->SetParameter(6, 3.9);
			//f->SetParName(7, "norm_0");				f->SetParameter(7, 0.9); //f->FixParameter(7, 1.);
			//f->SetParName(8, "x_0");				f->SetParameter(8, 7.);

			//+-1ns 41V nopz 1 SiPM switch 5 tune 8350
			f->SetParName(0, "N0");					f->SetParameter(0, his->Integral());
			f->SetParName(1, "#mu");				f->SetParameter(1, 0.7);// 1.6);
			f->SetParName(2, "#lambda");			f->SetParameter(2, .04); //0.2 or 3
			f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 2.1);//3.6
			f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 3.4);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
			f->SetParName(5, "Gain");				f->SetParameter(5, 18.);	//f->FixParameter(5, 10.);
			f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
			f->SetParName(7, "norm_0");				f->SetParameter(7, 0.7); //f->FixParameter(7, 1.);
			f->SetParName(8, "x_0");				f->SetParameter(8, 5.);

			//+-1ns 41V nopz 2 SiPMs
			//f->SetParName(0, "N0");					f->SetParameter(0, his->Integral());
			//f->SetParName(1, "#mu");				f->SetParameter(1, 4.);
			//f->SetParName(2, "#lambda");			f->SetParameter(2, .1);		//f->FixParameter(2, .05);
			//f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 4.);
			//f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, .7);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
			//f->SetParName(5, "Gain");				f->SetParameter(5, 15.2);	//f->FixParameter(5, 10.);
			//f->SetParName(6, "Pedestal");			f->SetParameter(6, 15.);
			//f->SetParName(7, "norm_0");				f->SetParameter(7, 1.5); //f->FixParameter(7, 1.);
			//f->SetParName(8, "x_0");				f->SetParameter(8, 16.);

			if (i < max_channel_nr_to_fit) {
				cout << "\n\n---------------------- Fit for channel " << active_channels[i] << " ----------------------\n";
				TFitResultPtr fresults = his->Fit(f, "RS");
			}

			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("integral in mV#timesns");
			his->Draw();
		}
	}
	chargec->Update();
	chargec->SaveAs("ChargeSpectra.pdf");
}

void ReadRun::PrintChargeSpectrumPMT(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, double threshold) {
	// print ReadRun::ChargeSpectrum for all channels optimized for PMT signals

	gStyle->SetOptStat(0); // 11 is title + entries

	TCanvas* chargec = new TCanvas("charge spectra PMT", "charge spectra PMT", 1600, 1000);

	if (plot_active_channels.empty()) chargec->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	else chargec->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
	int current_canvas = 0;
	float threshold_bin_center = 0.;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;

			TH1F* his;
			his = ChargeSpectrum(i, windowlow, windowhi, start, end, rangestart, rangeend, nbins);
			chargec->cd(current_canvas);

			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("integral in mV#timesns");
			his->Draw();
			stringstream allname; allname << his->GetEntries() << " entries";
			his->SetTitle(allname.str().c_str());

			auto his_lo = (TH1F*)his->Clone();
			his_lo->GetXaxis()->SetRange(his_lo->GetXaxis()->FindBin(rangestart), his_lo->GetXaxis()->FindBin(threshold));
			his_lo->SetLineColor(2);
			his_lo->SetFillColor(2);
			his_lo->Draw("LF2 same");
			stringstream loname;
			loname << 100. * his->Integral(his->GetXaxis()->FindBin(rangestart), his->GetXaxis()->FindBin(threshold)) / his->GetEntries() << "% <= " << threshold << " mV";
			his_lo->SetTitle(loname.str().c_str());

			auto his_hi = (TH1F*)his->Clone();
			his_hi->GetXaxis()->SetRange(his_hi->GetXaxis()->FindBin(threshold), his_lo->GetXaxis()->FindBin(rangeend));
			his_hi->SetLineColor(1);
			his_hi->SetFillColor(1);
			his_hi->Draw("LF2 same");
			stringstream hiname;
			hiname << 100. * his->Integral(his->GetXaxis()->FindBin(threshold) + 1, his->GetXaxis()->FindBin(rangeend)) / his->GetEntries() << "% > " << threshold << " mV";
			his_hi->SetTitle(hiname.str().c_str());

			threshold_bin_center = his->GetXaxis()->GetBinCenter(his->GetXaxis()->FindBin(threshold) + 1);

			gPad->BuildLegend();
		}
	}
	cout << "\n PMT charge spectrum is counting events above threshold from bin center >= " << threshold_bin_center << " mV " << "for a threshold setting of " << threshold << " mV\n\n";

	chargec->Update();
	chargec->SaveAs("ChargeSpectraPMT.pdf");
}

// time distribution of max in a certain time window
// Maybe add functions for ToT and max+min amplitude in range?

TH1F* ReadRun::TimeDist(int channel_index, float from, float to, float rangestart, float rangeend, int nbins) {
	// find peak time for a given channel in time window [from, to] and return the peak time histogram with x range [rangestart, rangeend] and the number of bins nbins

	TString name(Form("timedist_ch%02d", active_channels[channel_index]));
	auto h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);

	for (int j = 0; j < nevents; j++) {
		auto his = (TH1F*)((TH1F*)rundata->At(j * nchannels + channel_index))->Clone();
		if (from >= 0 && to > 0) his->GetXaxis()->SetRange(his->GetXaxis()->FindBin(from), his->GetXaxis()->FindBin(to));
		h1->Fill(his->GetXaxis()->GetBinCenter(his->GetMaximumBin()));	// fill peak time histogram
	}
	return h1;
}

void ReadRun::PrintTimeDist(float from, float to, float rangestart, float rangeend, int nbins) {
	// print ReadRun::TimeDist for all channels

	gStyle->SetOptStat(1111); // 11 is title + entries

	TCanvas* time_dist_c = new TCanvas("timing of maximum", "timing of maximum", 1600, 1000);

	if (plot_active_channels.empty()) time_dist_c->Divide(TMath::Min(static_cast<double>(active_channels.size()), 4.), TMath::Max(TMath::Ceil(static_cast<double>(active_channels.size()) / 4.), 1.), 0, 0);
	else time_dist_c->Divide(TMath::Min(static_cast<double>(plot_active_channels.size()), 4.), TMath::Max(ceil(static_cast<double>(plot_active_channels.size()) / 4.), 1.), 0, 0);
	int current_canvas = 0;

	for (int i = 0; i < nchannels; i++) {
		if (plot_active_channels.empty() || find(plot_active_channels.begin(), plot_active_channels.end(), active_channels[i]) != plot_active_channels.end()) {
			current_canvas++;

			TH1F* his;
			his = TimeDist(i, from, to, rangestart, rangeend, nbins);
			time_dist_c->cd(current_canvas);

			his->GetYaxis()->SetTitle("#Entries");
			his->GetXaxis()->SetTitle("integral in mV#timesns");
			his->Draw();
			stringstream name; name << "t#_{max} for " << from << "<t<" << to << " ns";
			his->SetTitle(name.str().c_str());
		}
	}

	time_dist_c->Update();
	time_dist_c->SaveAs("TimeDist.pdf");
}

// helper functions

stringstream ReadRun::list_files(const char* dirname, const char* ext) {
	// helper creating list of all bin files in directory
	stringstream ss;
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	if (files) {
		TSystemFile* file;
		TString fname;
		TIter next(files);
		while ((file = (TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext)) {
				ss << fname.Data() << "\n";
				cout << fname.Data() << "\n";
			}
		}
		TIter next2(files);
		while ((file = (TSystemFile*)next2())) {
			fname = file->GetName();
			if (!file->IsDirectory() && !fname.EndsWith(ext) && fname.Contains(ext)) {
				ss << fname.Data() << "\n";
				cout << fname.Data() << "\n";
			}
		}
	}
	return ss;
}

TH1F* ReadRun::Getwf(int channelnr, int eventnr, int color) {
	TH1F* his;
	if (eventnr > 0) eventnr -= 1; //rundata->At() counter starts at 0 which contains event 1 (wavecatcher event numbering starts with 1)
	his = (TH1F*)rundata->At(eventnr * nchannels + channelnr);
	his->SetLineColor(color);
	his->SetMarkerColor(color);
	return his;
}

double* ReadRun::getx() {
	double* xvals = new double[1024];
	for (int i = 0; i < 1024; i++) {
		xvals[i] = static_cast<double>(SP) * static_cast<double>(i);
	}
	return xvals;
}

double* ReadRun::gety(int channelnr, int eventnr) {
	TH1F* his = Getwf(channelnr, eventnr);
	double* yvals = new double[1024];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}

double* ReadRun::gety(TH1F* his) {
	double* yvals = new double[1024];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}

void ReadRun::Convolute(double*& result, double* first, double* second, int size1, int size2) {
	// Include FFT convolution
	// faster if size1<size2

	//for (int i = 0; i < size2; i++) {
	//	result[i] = 0.;
	//	for (int j = 0; j < TMath::Min(size1, i); j++) {
	//		result[i] += first[j] * second[i - j];
	//	}
	//}

	double* refirst = new double[size1];
	double* imfirst = new double[size1];
	double* resecond = new double[size1];
	double* imsecond = new double[size1];
	double* reres = new double[size1];
	double* imres = new double[size1];

	TVirtualFFT* fftfirst = TVirtualFFT::FFT(1, &size1, "R2C ES");
	fftfirst->SetPoints(first);
	fftfirst->Transform();
	fftfirst->GetPointsComplex(refirst, imfirst);
	delete fftfirst;

	TVirtualFFT* fftsecond = TVirtualFFT::FFT(1, &size1, "R2C ES");
	fftsecond->SetPoints(second);
	fftsecond->Transform();
	fftsecond->GetPointsComplex(resecond, imsecond);
	delete fftsecond;

	TComplex cofirst;
	TComplex cosecond;
	TComplex cores;

	for (int i = 0; i < size1; i++) {
		cofirst(refirst[i], imfirst[i]);
		cosecond(resecond[i], imsecond[i]);

		cores = cofirst * cosecond / static_cast<double>(size1);

		reres[i] = cores.Re();
		imres[i] = cores.Im();
	}

	//cout << "performing IFFT ... ";
	TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &size1, "C2R ES");
	fft_back->SetPointsComplex(reres, imres);
	fft_back->Transform();
	fft_back->GetPoints(result);
	delete fft_back;
	delete[] imres; delete[] reres; delete[] refirst; delete[] imfirst; delete[] resecond; delete[] imsecond;
}

void ReadRun::SmoothArray(double*& ar, int nbins, double sigma, bool doconv) {
	//apply smoothing array of double with length nbins
	// very inefficient

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	if (doconv) {
		// convolution with gauss with sigma (experimental, not yet tested)
		double* gauss = new double[nbins];

		double sum = 0.;

		for (int i = 0; i < nbins; i++) {
			gauss[i] = TMath::Exp(-1. * TMath::Power((static_cast<double>(i) * SP - 5 * sigma/*-> offset?*/), 2.) / (2. * sigma * sigma)) / (sigma * 2.506628);
			sum += gauss[i];
		}

		for (int i = 0; i < nbins; i++) {
			gauss[i] /= sum;
		}

		Convolute(ar, artmp, gauss, nbins, nbins);
		delete[] gauss;
	}
	else {
		// calculate running average from -sigma until +sigma
		for (int i = 0; i < nbins; i++) {
			double mean1 = 0.;
			int nmn = 0;
			for (int k = -1 * static_cast<int>(floor(sigma)); k <= static_cast<int>(ceil(sigma)); k++) {
				if (i + k >= 0 && i + k < nbins) {
					mean1 += artmp[i + k];
					nmn++;
				}
			}
			if (nmn != 0.) {
				ar[i] = mean1 / static_cast<double>(nmn);
			}
		}
	}
	delete[] artmp;
}

void ReadRun::PrintFFTWF(int eventnr, float xmin, float xmax, int multiplier) {
	// plot waveforms of all channels for a given event number eventnr and add the determined integration windwos to the plot
	TString name(Form("fft_waveforms_event__%04d", eventnr));
	TCanvas* fftc = new TCanvas(name.Data(), name.Data(), 1600, 1000);
	fftc->Divide(4, ceil(nchannels / 3), 0, 0);

	TString imname(Form("fft_im_waveforms_event__%04d", eventnr));
	TCanvas* imfftc = new TCanvas(imname.Data(), imname.Data(), 1600, 1000);
	imfftc->Divide(4, ceil(nchannels / 3), 0, 0);

	if (eventnr > 0) eventnr -= 1; //rundata->At() counter starts at 0 which contains event 1 (wavecatcher event numbering starts with 1)

	int size = 1024 * multiplier;

	double* xvals = new double[size];
	for (int i = 0; i < size; i++) {
		xvals[i] = static_cast<double>(i) / (SP * static_cast<double>(size));
	}

	double* refft = new double[size];
	double* imfft = new double[size];

	double* yvals = new double[size];

	for (int i = 0; i < nchannels; i++) {
		TH1F* his;
		his = ((TH1F*)rundata->At(eventnr * nchannels + i));

		for (int j = 0; j < size; j++) {
			if (j < 1024) yvals[j] = his->GetBinContent(j + 1);
			else yvals[j] = 0.;
		}

		TVirtualFFT* ffthis = TVirtualFFT::FFT(1, &size, "R2C ES");
		ffthis->SetPoints(yvals);
		ffthis->Transform();
		ffthis->GetPointsComplex(refft, imfft);

		TGraph* re = new TGraph(size, xvals, refft);
		TGraph* im = new TGraph(size, xvals, imfft);

		// draw to canvas
		fftc->cd(i + 1);
		stringstream renamess; renamess << "Channel " << i << ", event " << eventnr + 1 << ", Re(FFT(data))";
		re->SetTitle(renamess.str().c_str());
		re->Draw("AL");
		gPad->Modified();
		if (xmin != 0. || xmax != 0.) re->GetXaxis()->SetLimits(xmin, xmax);
		//re->GetYaxis()->SetRangeUser(-1*size, size);

		imfftc->cd(i + 1);
		stringstream imnamess; imnamess << "Channel " << i << ", event " << eventnr + 1 << ", Im(FFT(data))";
		im->SetTitle(imnamess.str().c_str());
		im->Draw("AL");
		gPad->Modified();
		if (xmin != 0. || xmax != 0.) im->GetXaxis()->SetLimits(xmin, xmax);
		//im->GetYaxis()->SetRangeUser(-1*size, size);

		delete ffthis;
	}
	fftc->Update();
	imfftc->Update();

	stringstream namess;
	namess << name.Data() << ".pdf";
	fftc->SaveAs(namess.str().c_str());
	stringstream imnamess;
	imnamess << imname.Data() << ".pdf";
	imfftc->SaveAs(imnamess.str().c_str());

	delete[] yvals;
	delete[] refft;
	delete[] imfft;
	delete[] xvals;
}

//TODO:
// 1: DC probability																	<-
// 2: implement method to discard individual waveforms (needed???)						<-
// 3: store file with fit parameters in data directory									<-
// 4: compile as library																<- Priority