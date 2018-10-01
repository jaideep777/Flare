#ifndef SIMPLE_INITIALIZER
#define SIMPLE_INITIALIZER

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;

/** \ingroup utils */

/**
	\brief A simple initializer that reads a parameter file and stores the values in a named map.
	
	The parameter file must have three sections - STRINGS, SCALARS, ARRAYS. Sections start with '>'. 
	Each section has name-value pairs separated by whitespace. Arrays have a name followed by values, ending in '-1'.
	Comments are allowed. Comments start with "# " (note the space) and can come either on a new line or on the same line after
	the name-value(s) pair. 
	
	Here is an example parameter file:
	
	~~~{.cpp}
	> STRINGS
	sim_name      mySimution
	output_file   ~/output/test.txt
	
	> SCALARS
	graphics      1           # Do we want graphics to be on? 
	timesteps     1000        # For how many timesteps do we run the simulation?
	dt            0.1
	# This is also a comment
	
	> ARRAYS
	parameter1    1 2 3 4 5 6 -1
	~~~
*/

class Initializer{
	private:
	string init_fname;
	ifstream fin;
	
	map <string, string> strings;
	map <string, float>  scalars;
	map <string, vector <float> > arrays;
	
	public:
	
	/** Default constructor */
	Initializer();
	
	/** The initializer can be created by specifying the parameter file name. */
	Initializer(string fname   //!< filename 
			   ); 
			   
	/** This function can be used to specify the parameter file. For example, if the initializer was created with the default constructor. */
	void setInitFile(string fname //!< filename 
					);

	/** Read values from the parameters file and store them in maps. */
	void readFile();

	/** Get string variable defined in the parameter file's STRINGS section. */
	string getString(string s //!< Name of the variable
					);

	/** Get a float variable defined in the parameter file's SCALARS section. */
	float getScalar(string s //!< Name of the variable 
				   );
	
	/** Read an array defined in the parameter file's ARRAYS section. 
		\returns vector of floating point numbers containing the array 
	*/
	vector <float> getArray(string s,	//!< name of the array 
							int size	//!< length of the array	
						   );

	/** Print out the values that have been read from the maps. */
	void printVars();
	
};



#endif

