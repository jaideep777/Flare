/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    libGSM - library for Gridded Spatial Modelling
    Copyright (C) 2016 Jaideep Joshi

	This file is part of libGSM.

    libGSM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

	The author can be contacted at jaideep777@gmail.com 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "../include/gsm.h"

//// general GSM options functions
//void setDebugFlag(bool f) {debug_on = f;};
//void setInfoFlag(bool f) {info_on = f;};
//void setWarningFlag(bool f) {warnings_on = f;};
//void setErrorFlag(bool f) {errors_on = f;};


// gives time string in hh-mm-ss.s format
string xhrs2hms(double dayf){
	stringstream sss;
	int h, m; float s1;
	dayf = dayf*24;
	h = int(dayf);
	double r = dayf - int(dayf);
	m = int(r*60);
	r = r*60 - int(r*60);
	s1 = int(r*600)/10.0;
	sss.clear();
	sss << h << ":" << m << ":" << s1;
	return sss.str();
}

// gives fraction of day (on 24 hrs) given time string
// time format: hh:mm:ss.s OR hh:mm:ss
double hms2xhrs(string time){
	stringstream ss;
	int sep1 = time.find(':');
	int sep2 = time.find(':',sep1+1);
	time[sep1] = ' '; time[sep2] = ' ';
	ss.clear(); ss.str(time);
	float hours, mins, secs;
	ss >> hours >> mins >> secs;
	double dayfrac = (double(hours) + mins/60.0 + secs/3600.0)/24.0;
	return dayfrac;
}

// calculate days since 01-Mar-0000 AD (arbit reference date)
// see: http://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
// inputs: y, m, d
int ymd2gday(int y, int m, int d){
	m = (m+9)%12;
	y = y - m/10;
	return 365*y + y/4 - y/100 + y/400 + (m*306 + 5)/10 + (d-1);
}

// calculate days since 01-Mar-0000 AD (arbit reference date)
// see: http://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
// input format: yyyy-mm-dd
int ymd2gday(string date){
	stringstream ss;
	int y, m, d;
	int sep1 = date.find('-');
	int sep2 = date.find('-',sep1+1);
	date[sep1] = ' '; date[sep2] = ' ';
	//cout << "date given = " << date << '\n';
	ss.clear(); ss.str(date);
	ss >> y >> m >> d;
	//cout << y << "," << m << "," << d << '\n';
	m = (m+9)%12;
	y = y - m/10;
	return 365*y + y/4 - y/100 + y/400 + (m*306 + 5)/10 + (d-1);
}

// calculate date from days since 01-Mar-0000 AD
// see: http://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
// * note the typecasting to long long int
// output format: yyyy-mm-dd
string gday2ymd(int g){
	//long g = long(n);
	int ystr, mstr, dstr;
	int y, ddd, mm, mi;
	y = (10000*(long long int)g + 14780)/3652425;
	ddd = g - (365*y + y/4 - y/100 + y/400);
	if (ddd < 0){
		--y;
		ddd = g - (365*y + y/4 - y/100 + y/400);
	}
	mi = (52 + 100*ddd)/3060;
	ystr = y + (mi + 2)/12;
	mstr = (mi + 2)%12 + 1;
	dstr = ddd - (mi*306 + 5)/10 + 1;
	
	stringstream ss;
	ss << ystr << '-' << mstr << '-' << dstr;
	return ss.str();
}

// calculate integer year given gt
int gt2year(double gt){
	int ystr;
	int y, ddd, mm, mi;
	y = (10000*(long long int)gt + 14780)/3652425;
	ddd = gt - (365*y + y/4 - y/100 + y/400);
	if (ddd < 0){
		--y;
		ddd = gt - (365*y + y/4 - y/100 + y/400);
	}
	mi = (52 + 100*ddd)/3060;
	ystr = y + (mi + 2)/12;
	return ystr;
}

// calculate month given gt
int gt2month(double gt){
	int mstr;
	int y, ddd, mm, mi;
	y = (10000*(long long int)gt + 14780)/3652425;
	ddd = gt - (365*y + y/4 - y/100 + y/400);
	if (ddd < 0){
		--y;
		ddd = gt - (365*y + y/4 - y/100 + y/400);
	}
	mi = (52 + 100*ddd)/3060;
	mstr = (mi + 2)%12 + 1;
	return mstr;
}

// calculate day given GT g
int gt2day(double g){
	//long g = long(n);
	int ystr, mstr, dstr;
	int y, ddd, mm, mi;
	y = (10000*(long long int)g + 14780)/3652425;
	ddd = g - (365*y + y/4 - y/100 + y/400);
	if (ddd < 0){
		--y;
		ddd = g - (365*y + y/4 - y/100 + y/400);
	}
	mi = (52 + 100*ddd)/3060;
	ystr = y + (mi + 2)/12;
	mstr = (mi + 2)%12 + 1;
	dstr = ddd - (mi*306 + 5)/10 + 1;
	return dstr;
}

int gt2daynum(double g){
	//long g = long(n);
	int ystr, mstr, dstr;
	int y, ddd, mm, mi;
	y = (10000*(long long int)g + 14780)/3652425;
	ddd = g - (365*y + y/4 - y/100 + y/400);
	if (ddd < 0){
		--y;
		ddd = g - (365*y + y/4 - y/100 + y/400);
	}
	mi = (52 + 100*ddd)/3060;
	ystr = y + (mi + 2)/12;
	mstr = (mi + 2)%12 + 1;
	dstr = ddd - (mi*306 + 5)/10 + 1;

	return g - ymd2gday(ystr, 1, 1);
}

int gt2dayOfYear(double gt){
	return gt2daynum(gt);
}


string gt2string(double gt){
	int gday_day = int(gt);
	float gday_hours = gt - gday_day;	// use of double here gives 5:29:59.9 for 5:30:0!!
	// for the format hh:mm:ss.x float seems to suffice.. error occurs in the last decimal (x)
	// in conclusion, both float and double give error of 0.1s in different situations
	return (gday2ymd(gday_day) + " " + xhrs2hms(gday_hours));
}


string gt2string_date(double gt){
	int gday_day = int(gt);
	return (gday2ymd(gday_day));
}


string gt2string_time(double gt){
	int gday_day = int(gt);
	float gday_hours = gt - gday_day;	// use of double here gives 5:29:59.9 for 5:30:0!!
	// for the format hh:mm:ss.x float seems to suffice.. error occurs in the last decimal (x)
	// in conclusion, both float and double give error of 0.1s in different situations
	return (xhrs2hms(gday_hours));
}


void gt2array(double gt, int* tarr){
	int g = int(gt);
	double dayf = gt - g;	// use of double here gives 5:29:59.9 for 5:30:0!!

	// get day in yyyy, mm, dd
	int ystr, mstr, dstr;
	int y, ddd, mm, mi;
	y = (10000*(long long int)g + 14780)/3652425;
	ddd = g - (365*y + y/4 - y/100 + y/400);
	if (ddd < 0){
		--y;
		ddd = g - (365*y + y/4 - y/100 + y/400);
	}
	mi = (52 + 100*ddd)/3060;
	tarr[0] = ystr = y + (mi + 2)/12;
	tarr[1] = mstr = (mi + 2)%12 + 1;
	tarr[2] = dstr = ddd - (mi*306 + 5)/10 + 1;
	
	// get time in hh, mm, ss
	int h, m; float s1;
	dayf = dayf*24;
	h = int(dayf);
	double r = dayf - int(dayf);
	tarr[3] = m = int(r*60);
	tarr[4] = r = r*60 - int(r*60);
	tarr[5] = s1 = int(r*600)/10.0;
	 	
}


string gtstr6d(double gt){
	return (int2str(gt) + "." + int2str(int((gt-int(gt))*1e6)) );
}


/****************** FEW lat/lon functions.. included here for now. ***********/

// convert lat/lon from sexagesimal string to decimal float
float sex2decLL(string s){
	stringstream ss;
	ss.clear(); ss.str(s);
	float deg, min, sec;
	char dir;
	ss >> deg >> min >> sec >> dir;
	float sval = deg + min/60 + sec/3600;
	if (dir == 'W' || dir == 'S') sval = -sval;
	return sval;
}

// convert lat/lon from decimal float to sexagesimal string 
string dec2sexLL(float lon){
	stringstream ss;
	int d, m, s;
	d = int(lon);
	float r = lon - int(lon);
	m = int(r*60);
	r = r*60 - int(r*60);
	s = int(r*60);
	ss << d << " " << m << " " << s;
	return ss.str();
}

