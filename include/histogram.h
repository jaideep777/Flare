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

/**
	@defgroup utils Utilities
	@brief Various utility functions and classes, such as vector math, colour palettes, histograms, and date-time arithmatic.
*/

/** @ingroup utils */

/** 
	\brief A histogram class based on gsl_histogram
*/
class Histogram{
	public:
	gsl_histogram * h;	//!< base GSL histogram 
	
	Histogram();		//!< Default constructor
	
	/** \brief Create histogram in one step by specifying number of breaks. 
	*/ 
	Histogram(vector <float> &data, 	//!< Data from which to create histogram
			  int nbins, 				//!< Number of bins (Equally spaced bins will be created)
			  float range_min = 1e20, 	//!< Min value (if not specified, will be calculated from the data)
			  float range_max = 1e20	//!< Max value (if not specified, will be calculated from the data)
			 );
			  
	/** \brief Create histogram in one step by specifying the breaks. 
	*/ 
	Histogram(vector <float> &data, 	//!< Data from which to create histogram
			  vector <double> &breaks	//!< Breaks
			 );

	/** \brief Create weighted histogram in one step by specifying number of breaks. 
	*/ 	
	Histogram(vector <float> &data, 	//!< Data from which to create histogram
			  vector <float> &w, 		//!< Weights (multiplied to the data before adding to counts)
			  int nbins, 				//!< Number of bins (Equally spaced bins will be created)
			  float range_min = 1e20, 	//!< Min value (if not specified, will be calculated from the data)
			  float range_max = 1e20	//!< Max value (if not specified, will be calculated from the data)
			 );
			 
	/** \brief Create weighted histogram in one step by specifying breaks. 
	*/ 
	Histogram(vector <float> &data, 	//!< Data from which to create histogram
			  vector <float> &w, 		//!< Weights (multiplied to the data before adding to counts)
			  vector <double> &breaks	//!< Breaks
			 );
	
	~Histogram();
	
	int plot_console();				//!< Plot the histogram to console 
	vector <float> getCounts();		//!< Get the counts vector from the histogram 
	vector <float> getMids();		//!< Get bin midpoints (Arithmatic mean of the bin ends)
	vector <float> getMids_log();	//!< Get bin midpoints (Geometric mean of the bin ends)
	vector <float> getBreaks();		//!< Get the breaks vector
	int convertToPdf();				//!< Normalize the counts to a probability distribution \f$\sum c = 1 \f$
};



/** @addtogroup utils 
	@{ 
*/

// =============================================================================
/** 	@brief Data summaries
 * 		@author Jaideep Joshi
 * 		@date 11 May 2015
 
 		Print the summary of given data (min, max, mean, and histogram) 
 */
// =============================================================================
template <class T> 
void printSummary(T * data, 		//!< Data array
				  int n, 			//!< Numeber of elements (array size)
				  string s = ""		//!< Name of the array to prefix the printed output
				 ){
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

/** @} */



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
