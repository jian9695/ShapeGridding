#include "Rasterization.h"

#pragma warning( disable: 4049 ) 
#pragma warning( disable: 4201 ) 


void makeYilongGrid(std::string gridPath)
{
  //  dimensions  : 4260, 9060, 38595600  (nrow, ncol, ncell)
  //  resolution : 0.008333333, 0.008333333  (x, y)
  //  extent : -137.5, -62, 22, 57.5  (xmin, xmax, ymin, ymax)

  if (!QDir(gridPath.data()).exists())
    QDir(gridPath.data()).mkpath(".");

  OGRSpatialReference ogrSRS;
  ogrSRS.SetWellKnownGeogCS("WGS84");
  char projBuf[1024];
  char* pProjWKT = projBuf;
  ogrSRS.exportToWkt(&pProjWKT);

  OGREnvelope bound;
  bound.MinX = -137.5;
  bound.MaxX = -62;
  bound.MinY = 22;
  bound.MaxY = 57.5;
  int nrows = 4260;
  int ncols = 9060;
  GridConfig gridCFG;
  gridCFG.cellsize = (bound.MaxX- bound.MinX)/ ncols;
  gridCFG.bound = bound;
  gridCFG.nrows = nrows;
  gridCFG.ncols = ncols;
  gridCFG.proj = pProjWKT;

  gridCFG.exportConfig(gridPath + "grid.cfg");
  gridCFG.toRaster(gridPath + "grid.tif");

}

// Run BIN2RDS.R to convert binary spatial fractions files (.bin) to RDS format (.rds)

int main(int argc, char** argv)
{
	OGRRegisterAll();
	GDALAllRegister();

  makeYilongGrid("./Samples/grids/Yilong/");
  Rasterization::intersectFolder("./Samples/shapefiles/Vulcan/", "./Samples/grids/Yilong/grid.tif", "./Samples/shapefile_grid_intersection_fractions/Yilong/");

	return 0;

}
