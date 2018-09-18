// =============================================================================
// 		The Histogram Class
// 		ADDED by : JAIDEEP
// 		21 Dec 2013
//
//		to compile, link with gsl: -lgsl, -lgslcblas
//
// =============================================================================

#include <gsl/gsl_histogram.h>
#include <iomanip>
#include <vector>

#include "../include/histogram.h"

Histogram::Histogram(){};

Histogram::Histogram(vector <float> &data, int nbins, float range_min, float range_max){

	h = gsl_histogram_alloc (nbins);

	if (range_min > 9e19) range_min = arrayMin(&data[0], data.size()); //cout << "min = " << data_min << '\n';
	if (range_max > 9e19) range_max = arrayMax(&data[0], data.size()); //cout << "min = " << data_min << '\n';

	if (range_min == range_max) {
		range_min -= 1e-3;
		range_max += 1e-3;
	}

	gsl_histogram_set_ranges_uniform (h, range_min, range_max);

	// create histogram
	gsl_histogram_reset(h);
	for (int i=0; i<data.size(); ++i) gsl_histogram_increment(h, data[i]);
}


Histogram::Histogram(vector <float> &data, vector <double> &breaks){
	h = gsl_histogram_alloc (breaks.size()-1);

	gsl_histogram_set_ranges (h, &breaks[0], breaks.size());

	// create histogram
	gsl_histogram_reset(h);
	for (int i=0; i<data.size(); ++i) gsl_histogram_increment(h, data[i]);
}


Histogram::Histogram(vector <float> &data, vector <float> &w, int nbins, float range_min, float range_max){

	h = gsl_histogram_alloc (nbins);

	if (range_min > 9e19) range_min = arrayMin(&data[0], data.size()); //cout << "min = " << data_min << '\n';
	if (range_max > 9e19) range_max = arrayMax(&data[0], data.size()); //cout << "min = " << data_min << '\n';

	gsl_histogram_set_ranges_uniform (h, range_min, range_max);

	// create histogram
	gsl_histogram_reset(h);
	for (int i=0; i<data.size(); ++i) gsl_histogram_accumulate(h, data[i], w[i]);
}


Histogram::Histogram(vector <float> &data, vector <float> &w, vector <double> &breaks){
	h = gsl_histogram_alloc (breaks.size()-1);

	gsl_histogram_set_ranges (h, &breaks[0], breaks.size());

	// create histogram
	gsl_histogram_reset(h);
	for (int i=0; i<data.size(); ++i) gsl_histogram_accumulate(h, data[i], w[i]);
}


Histogram::~Histogram(){
	gsl_histogram_free(h);
}

int Histogram::plot_console(){
	vector <float> b = getMids();
	cout << setprecision(2) << fixed;
	for (int i=0; i<h->n; ++i){
		cout << setw(7) << b[i] << " <- " << setw(7) << h->bin[i] << '\t';
		int max_count = gsl_histogram_max_val(h);
		int step = max_count/40 +1;
		for (int k=0; k<h->bin[i]; k+=step) cout << ".";
		cout << '\n';
	}
	cout.copyfmt(ios(NULL));	// reset cout state
	
	return 0;
}

vector <float> Histogram::getCounts(){
	vector <float> freq(h->n, 0);
	for (int i=0; i<h->n; ++i){
		freq[i] = h->bin[i];
	}
	return freq;
}

vector <float> Histogram::getMids(){
	vector <float> freq(h->n, 0);
	for (int i=0; i<h->n; ++i){
		freq[i] = (h->range[i] + h->range[i+1])/2;
	}
	return freq;
}

vector <float> Histogram::getMids_log(){
	vector <float> freq(h->n, 0);
	for (int i=0; i<h->n; ++i){
		freq[i] = sqrt(h->range[i] * h->range[i+1]);
	}
	return freq;
}

vector <float> Histogram::getBreaks(){
	vector <float> freq(h->n+1, 0);
	for (int i=0; i<h->n+1; ++i){
		freq[i] = h->range[i];
	}
	return freq;
}

int Histogram::convertToPdf(){
	float hsum = gsl_histogram_sum(h);
	for (int i=0; i<h->n; ++i){
		h->bin[i] /= (hsum+1e-6);
	}
	return 0;
	
}

// test program

//int main(){
//	int N = 10000;
//	vector <float> data(N), w(N,1);
//	for (int i=0; i<N; ++i) {
//		data[i] = norm_randf()*10;
//		if (data[i]<0) w[i] = 0;
//	}
//	double breaks_array[] = {-20,-15,-10,-5,0,5,10,15,20};
//	vector <double> breaks(breaks_array, breaks_array+9);
//	
//	Histogram h(data,w,32,-40,40);
//	
//	h.plot_console();

//	vector <float> counts = h.getMids();
//	for (int i=0; i<counts.size(); ++i) cout << counts[i] << " ";
//	cout << '\n';
//	
//	return 0;
//}


