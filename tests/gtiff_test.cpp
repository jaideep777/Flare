#include "gdal_priv.h"
//#include "cpl_conv.h" // for CPLMalloc()
using namespace std;

int main(){
	GDALDataset  *poDataset;
	GDALAllRegister();
	poDataset = (GDALDataset *) GDALOpen( "/media/jaideep/WorkData/Fire_G/GPP_modis/gpp_tiff/MOD17A2_GPP.2000.M01.tif", GA_ReadOnly );
	if( poDataset == NULL )
	{
		cout << "success\n";
	}

	return 0;
}

