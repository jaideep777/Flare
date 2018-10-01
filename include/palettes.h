/*
 * A very simple library to create colours and palettes 
 *
 * Jaideep Joshi (24 Dec 2013)
 *
 */


#ifndef _PALETTES_H_
#define _PALETTES_H_


/** @addtogroup utils 
	@{
*/

#include <iostream>
#include <vector>
#include <cstdlib>		// for random numbers
using namespace std;

class Colour_rgb{
	public:
	float r;
	float g;
	float b;
	
	Colour_rgb();
	Colour_rgb(float rr, float gg, float bb);
};


//TODO: Make this into a Palette class

Colour_rgb HSVtoRGB(float h, float s, float v );

vector <Colour_rgb> createPalette_rainbow(int N, float start, float end);

vector <Colour_rgb> createPalette_random(int N, float start, float end);

vector <Colour_rgb> createPalette_grayscale(int N, float start, float end);

vector <Colour_rgb> createPalette_ramp(int N, Colour_rgb start, Colour_rgb end);

void printPalette(vector <Colour_rgb> &p);

#endif

/** @} */

// sample program
// ------------------

//int main(){
//	int N = 30;
//	
//	vector <Colour_rgb> p = createPalette_rainbow(30);

//	for (int i=0; i<N; ++i){
//		Colour_rgb c = p[i];
//		cout << int(c.r*255) << "\t" << int(c.g*255) << "\t" << int(c.b*255) << "\n";
//	}
//	return 0;	
//}


