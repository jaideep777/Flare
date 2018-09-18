#ifndef SIMPLE_INITIALIZER
#define SIMPLE_INITIALIZER

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;

class Initializer{
	private:
	string init_fname;
	ifstream fin;
	
	map <string, string> strings;
	map <string, float>  scalars;
	map <string, vector <float> > arrays;
	
	public:
	Initializer();
	Initializer(string fname);

	void setInitFile(string fname);
	void readFile();
	
	string getString(string s);
	float getScalar(string s);
	vector <float> getArray(string s, int size);

	void printVars();
	
};



#endif

