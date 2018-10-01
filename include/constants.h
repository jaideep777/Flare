#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <iostream>
using namespace std;

/**
	@brief Constants defined in Flare

	Constants predefined in the library.
	@author Jaideep Joshi
	@date Sept 2018	
	
*/

/** @defgroup constants Definitions
	@{
*/

extern ostream * gsm_log; 

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  DEFS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define CDEBUG if (gsm_debug_on) (*gsm_log) << "<GSM debug> "
#define CDEBUGC if (gsm_debug_on) (*gsm_log)
#define CINFO if (gsm_info_on) (*gsm_log) << "<GSM info> "
#define CINFOC if (gsm_info_on) (*gsm_log)
#define CWARN if (gsm_warnings_on) cout << "<GSM WARNING> "
#define CERR if (gsm_errors_on) cout << "<GSM ERROR> "

const double t_tol = 1e-3;	// tolerance when comparing gt values
							// corresponds to ~86 sec
const float std_missing_value = 9.9e20;

extern bool gsm_info_on;
extern bool gsm_debug_on;
extern bool gsm_warnings_on;
extern bool gsm_errors_on;

//const float glimits_globe[4] = {0, 360, -90, 90};
//const float glimits_india[4] = {66.5, 100.5, 6.5, 38.5};
//extern float glimits_custom[4];

/** 
	@}
*/

#endif


