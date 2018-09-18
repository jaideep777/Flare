#ifndef DATETIME_H
#define DATETIME_H

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  TIME    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//	Definitions:

//	gday = global day = integer days since 1-3-0000 AD
//	xhrs = hours elapsed over and above integer days as fraction of 24 hrs. eg: 18 hrs = 0.75 xhrs
//	gt = global time = decimal days since 1-3-0000 AD (gday + xhrs)


// notes: 
// fractional day is on base 10 not base 24
// global time is in "days since 1-3-0000 AD"

/*~~~~~~~~~~~~~~~~~~~~~~~~~  Date and Time functions  ~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string>
using namespace std;

int daysInMonth(int yr, int mon); // get days in month (31, 28/29, 31, 30 ...)
string xhrs2hms(double dayf); // convert fractional day (dayf) to hh-mm-ss.s string
double hms2xhrs(string time); // convert hh-mm-ss(.s) string to fractional day 
int    ymd2gday(string date); // convert date string yyyy-mm-dd to global day
int    ymd2gday(int year, int month, int day); // convert y, m, d to global day
string gday2ymd(int g); // convert global time to readable date string yyyy-mm-dd
string gt2string(double gt); // convert gday (including day fraction) to full date-time
string gt2string_date(double gt);	// " date only
string gt2string_time(double gt);   // " time only
string gtstr6d(double gt);		// print gt upto 6 decimals

int gt2year(double gt);	// calculate year only (non-decimal) from gday
int gt2month(double gt);	// calculate current month from gday (NOTE: month ranges from 1-12 and NOT from 0-11)
int gt2day(double g);		// calculate day in month
int gt2daynum(double gt);	// calculate day of year
int gt2dayOfYear(double gt);	// calculate day of year
void gt2array(double gt, int* tarr);	// get yyyy, MM, dd, hh, mm, ss in array from gt

// ------  lat lon functions ------- 
float sex2decLL(string s);	// convert "ll mm ss D" to decimal lat/lon
string dec2sexLL(float lon); // convert decimal lat/lon to "ll mm ss D"


#endif
