#ifndef JAI_HISTOGRAM_H
#define JAI_HISTOGRAM_H

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

#include "simple_math.h"

class Histogram{
	public:
	gsl_histogram * h;
	
	Histogram();
	Histogram(vector <float> &data, int nbins, float range_min = 1e20, float range_max = 1e20);
	Histogram(vector <float> &data, vector <double> &breaks);

	Histogram(vector <float> &data, vector <float> &w, int nbins, float range_min = 1e20, float range_max = 1e20);
	Histogram(vector <float> &data, vector <float> &w, vector <double> &breaks);
	
	~Histogram();
	
	int plot_console();
	vector <float> getCounts();
	vector <float> getMids();
	vector <float> getMids_log();
	vector <float> getBreaks();
	int convertToPdf();
};



// =============================================================================
// 		Data summaries
// 		ADDED by : JAIDEEP
// 		11 May 2015
// =============================================================================

template <class T> 
void printSummary(T * data, int n, string s = ""){
	T dat_sum = data[0], dat_max = data[0], dat_min = data[0];
	for (int i=1; i<n; ++i){
		dat_sum += data[i];
		dat_min = min(dat_min, data[i]);
		dat_max = max(dat_max, data[i]);
	}
	float dat_mean = float(dat_sum)/n;

	cout << s << " summary: " << dat_min << " - " << dat_mean << " - " << dat_max << endl;
	vector <float> v(n);
	for (int i=0; i<n; ++i) v[i] = data[i];
	Histogram h(v, 20);

	h.plot_console();
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


#endif
