#pragma once
#include "ShapeFile.h"
#include <vector>
#include "gdal_priv.h"
#include "Utils.h"
#include <iomanip>      // std::setprecision
#include <fstream>
#include "netcdf.h"
#define _NODATA 3.40282306074e+038

struct GridConfig
{
	OGREnvelope bound;
	double cellsize;
	std::string proj;
	int nrows;
	int ncols;
	GridConfig()
	{
		proj = "";
	}
	GridConfig(std::string configfile)
	{
		proj = "";
		fromConfig(configfile);
	}
	void getProjFromShapefile(std::string shapefile)
	{
		ShapeFile shp(shapefile);
		char buf[1024];
		char* pwkt = buf;
		shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		proj = pwkt;
	}
	void getProjFromRaster(std::string rasterfile)
	{
		GDALDataset* pDataset = (GDALDataset*)GDALOpen(rasterfile.data(), GA_ReadOnly);
		proj = "";
		if (pDataset->GetProjectionRef())
			proj = pDataset->GetProjectionRef();
		GDALClose((GDALDatasetH)pDataset);
	}

	void fromRaster(std::string rasterfile)
	{
		double adfGeoTransform[6];
		GDALDataset* pDataset = (GDALDataset*)GDALOpen(rasterfile.data(), GA_ReadOnly);
		pDataset->GetGeoTransform(adfGeoTransform);
		bound.MinX = adfGeoTransform[0];
		bound.MinY = adfGeoTransform[3] + adfGeoTransform[5] * pDataset->GetRasterYSize();
		bound.MaxX = adfGeoTransform[0] + adfGeoTransform[1] * pDataset->GetRasterXSize();
		bound.MaxY = adfGeoTransform[3];
		ncols = pDataset->GetRasterXSize();
		nrows = pDataset->GetRasterYSize();
		cellsize = adfGeoTransform[1];
		proj = "";
		if (pDataset->GetProjectionRef())
			proj = pDataset->GetProjectionRef();
		GDALClose((GDALDatasetH)pDataset);
	}
	void toRaster(std::string rasterfile)
	{
		GDALAllRegister();
		const char *pszFormat = "GTiff";
		char **papszOptions = NULL;
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		GDALDataset* pDataset = poDriver->Create(rasterfile.data(), ncols, nrows, 1, GDT_Byte, papszOptions);
		double adfGeoTransform[6];
		adfGeoTransform[0] = bound.MinX;
		adfGeoTransform[3] = bound.MaxY;
		adfGeoTransform[1] = cellsize;
		adfGeoTransform[5] = -cellsize;
		adfGeoTransform[2] = 0;
		adfGeoTransform[4] = 0;
		pDataset->SetGeoTransform(adfGeoTransform);
		GDALRasterBand *pBand = pDataset->GetRasterBand(1);
		pDataset->SetProjection(proj.data());
		GDALClose((GDALDatasetH)pDataset);
	}
	void exportConfig(std::string configfile)
	{
		std::ofstream ofs;
		ofs.open(configfile);
		ofs << std::setprecision(10);
		ofs << std::fixed << "minx" << " " << bound.MinX << std::endl;
		ofs << std::fixed << "miny" << " " << bound.MinY << std::endl;
		ofs << std::fixed << "maxx" << " " << bound.MaxX << std::endl;
		ofs << std::fixed << "maxy" << " " << bound.MaxY << std::endl;
		ofs << std::fixed << "ncols" << " " << ncols << std::endl;
		ofs << std::fixed << "nrows" << " " << nrows << std::endl;
		ofs << std::fixed << "cellsize" << " " << cellsize << std::endl;
		ofs << "projection" << " " << proj;
		ofs.close();
	}
	void fromConfig(std::string configfile)
	{
		std::ifstream ifs;
		ifs.open(configfile);
		std::string fieldname;
		ifs >> fieldname >> bound.MinX;
		ifs >> fieldname >> bound.MinY;
		ifs >> fieldname >> bound.MaxX;
		ifs >> fieldname >> bound.MaxY;
		ifs >> fieldname >> ncols;
		ifs >> fieldname >> nrows;
		ifs >> fieldname >> cellsize;
		std::string line;
		while (line.size() < 10 || line.substr(0, 10) != "projection")
		{
			std::getline(ifs, line);
		}
		int space = 10;
		for (size_t i = 10; i < line.size(); i++)
		{
			if (line[i] == ' ')
			{
				space = i;
				break;
			}
		}
		space = space + 1;
		fieldname = line.substr(0, 10);
		proj = line.substr(space, line.size() - space);
		//adfGeoTransform[5] = -adfGeoTransform[1];
		//adfGeoTransform[2] = 0;
		//adfGeoTransform[4] = 0;
		//adfGeoTransform[0] = bound.MaxY;
		//adfGeoTransform[3] = bound.MaxY;
		ifs.close();
	}
	GridConfig conformGrid(OGREnvelope newbound)
	{
		double minx = bound.MinX;
		if (newbound.MinX < bound.MinX)
		{
			while (minx > newbound.MinX)
			{
				minx -= cellsize;
			}
		}
		else
		{
			while (minx < newbound.MinX)
			{
				minx += cellsize;
			}
		}

		double maxy = bound.MaxY;
		if (newbound.MaxY < bound.MaxY)
		{
			while (maxy > newbound.MaxY)
			{
				maxy -= cellsize;
			}
		}
		else
		{
			while (maxy < newbound.MaxY)
			{
				maxy += cellsize;
			}
		}

		if (minx > newbound.MinX)
			minx -= cellsize;
		if (maxy < newbound.MaxY)
			maxy -= cellsize;

		double maxx = minx;
		double miny = maxy;
		GridConfig newgrid;
		newgrid.ncols = 0;
		while (newbound.MaxX > maxx)
		{
			maxx += cellsize;
			newgrid.ncols++;
		}
		newgrid.nrows = 0;
		while (newbound.MinY < miny)
		{
			miny -= cellsize;
			newgrid.nrows++;
		}


		newgrid.bound.MinX = minx;
		newgrid.bound.MinY = miny;
		newgrid.bound.MaxX = maxx;
		newgrid.bound.MaxY = maxy;
		newgrid.cellsize = cellsize;
		newgrid.proj = proj;

		printf("%f,%f\n", (newgrid.bound.MinY - bound.MinY) / cellsize, (bound.MinX - newgrid.bound.MinX) / cellsize);
		return newgrid;
	}

};

struct NetCDFAttribute
{
	//char name[255];
	std::string name;
	//char text[255];
	std::string text;
	size_t len;
	double value;
	nc_type vtype;
};

struct NetCDFVar
{
	//char name[255];
	std::string name;
	int id;
	int ndims;
	int dimids[4];
	int varnatt;
	nc_type vtype;
	std::vector<NetCDFAttribute> attributes;
};

struct NetCDFDim
{
	//char name[255];
	std::string name;
	size_t len;
	int id;
};

class NCFile
{
public:
	int _file;
	std::vector<NetCDFDim> _dims;
	std::vector<NetCDFVar> _vars;
	std::string _filename;
	void close()
	{
		nc_close(_file);
	}
	void create(const char* outputFile, int ncols, int nrows, int timedimsize, std::string timeunit, std::string spatialunit = "feet")
	{
		std::string xname = "X";
		std::string yname = "Y";
		if (spatialunit == "degree" || spatialunit == "degrees")
		{
			xname = "Lon";
			yname = "Lat";
		}
		std::vector<NetCDFDim> dims;
		NetCDFDim londim;
		londim.len = ncols;
		londim.name = xname;
		dims.push_back(londim);

		NetCDFDim latdim;
		latdim.len = nrows;
		latdim.name = yname;
		dims.push_back(latdim);

		if (timedimsize > 1)
		{
			NetCDFDim timedim;
			timedim.len = timedimsize;
			timedim.name = "Time";
			dims.push_back(timedim);
		}

		NetCDFAttribute attr;
		attr.len = 10;
		std::vector<NetCDFVar> vars;
		NetCDFVar lonvar;
		lonvar.dimids[0] = 0;
		lonvar.ndims = 1;
		lonvar.name = xname;
		lonvar.vtype = NC_DOUBLE;
		lonvar.varnatt = 2;
		attr.name = "units";
		attr.text = spatialunit.data();
		attr.vtype = 2;
		lonvar.attributes.push_back(attr);
		attr.name = xname;
		attr.text = xname;
		attr.vtype = 2;
		lonvar.attributes.push_back(attr);
		vars.push_back(lonvar);

		NetCDFVar latvar;
		latvar.dimids[0] = 1;
		latvar.ndims = 1;
		latvar.name = yname;
		latvar.vtype = NC_DOUBLE;
		latvar.varnatt = 2;
		attr.name = "units";
		attr.text = spatialunit.data();
		attr.vtype = 2;
		latvar.attributes.push_back(attr);
		attr.name = yname;
		attr.text = yname;
		attr.vtype = 2;
		latvar.attributes.push_back(attr);
		vars.push_back(latvar);

		NetCDFVar cavar;
		if (timedimsize > 1)
		{

			NetCDFVar timevar;
			timevar.dimids[0] = 2;
			timevar.ndims = 1;
			timevar.name = "Time";
			timevar.vtype = NC_INT;
			timevar.varnatt = 2;
			attr.name = "units";
			attr.text = timeunit.data();
			attr.vtype = NC_CHAR;
			timevar.attributes.push_back(attr);
			attr.name = "Time";
			attr.text = "Time";
			attr.vtype = NC_CHAR;
			timevar.attributes.push_back(attr);
			vars.push_back(timevar);


			cavar.dimids[0] = 2;
			cavar.dimids[1] = 1;
			cavar.dimids[2] = 0;
			cavar.ndims = 3;
		}
		else
		{
			cavar.dimids[0] = 1;
			cavar.dimids[1] = 0;
			cavar.ndims = 2;
		}
		cavar.name = "Carbon Emission";
		cavar.vtype = NC_DOUBLE;
		cavar.varnatt = 2;
		attr.name = "units";
		attr.text = "Metric Tonnes of Carbon (tC)";
		attr.vtype = NC_CHAR;
		cavar.attributes.push_back(attr);
		attr.name = "_FillValue";
		attr.value = _NODATA;
		attr.vtype = NC_DOUBLE;
		cavar.attributes.push_back(attr);
		vars.push_back(cavar);

		_dims = dims;
		_vars = vars;
		_filename = outputFile;

		for (size_t i = 0; i < _vars.size(); i++)
		{
			_vars[i].id = i;
		}
		//if (QFileInfo(outputFile).exists())
		//	return;
		int status = nc_create(outputFile, NC_CLOBBER, &(_file));

		//size_t start[] = { 0, 0 }; // start at first value 
		std::vector<size_t> latlondim;
		for (int i = 0; i < dims.size(); i++)
		{
			latlondim.push_back(dims[i].len);
			nc_def_dim(_file, dims[i].name.data(), dims[i].len, &i);
		}
		size_t* count = &latlondim[0];

		for (int i = 0; i < vars.size(); i++)
		{
			NetCDFVar var = vars[i];


			std::vector<size_t> dimlens;
			for (size_t j = 0; j < var.ndims; j++)
			{
				dimlens.push_back(dims[var.dimids[j]].len);
			}
			nc_def_var(_file, var.name.data(), var.vtype, var.ndims, var.dimids, &var.id);
			for (int j = 0; j < var.attributes.size(); j++)
			{
				NetCDFAttribute attr = var.attributes[j];
				if (attr.vtype == NC_DOUBLE)
					nc_put_att_double(_file, var.id, attr.name.data(), NC_DOUBLE, 1, &(attr.value));
				else
					nc_put_att_text(_file, var.id, attr.name.data(), attr.text.size(), attr.text.data());
			}
		}
		nc_close(_file);
		nc_open(_filename.data(), NC_WRITE | NC_SHARE, &_file);

	}
	void writeSlice(int timeSlice, double* data)
	{
		int err;
		//nc_open(_filename.data(), NC_WRITE | NC_SHARE, &_file);
		NetCDFVar var = _vars[_vars.size() - 1];
		err = nc_inq_varid(_file, var.name.data(), &var.id);
		std::vector<size_t> startArr;
		std::vector<size_t> latlondim;
		for (int i = 0; i < var.ndims; i++)
		{
			startArr.push_back(0);
			latlondim.push_back(_dims[var.dimids[i]].len);
		}
		startArr[0] = timeSlice;// *_dims[0].len * _dims[1].len;
		//latlondim[0] = 1;// *_dims[0].len * _dims[1].len;
		size_t* start = &startArr[0];
		size_t* count = &latlondim[0];
		err = nc_put_vara_double(_file, var.id, start, count, data);
		if (err != 0)
		{
			if (err == NC_NOERR)
			{
				printf("No error.\n");
			}
			else if (err == NC_ENOTVAR)
			{
				printf("Variable not found.\n");
			}
			else if (err == NC_EINVALCOORDS)
			{
				printf("Index exceeds dimension bound.\n");
			}
			else if (err == NC_EEDGE)
			{
				printf("Start + count exceeds dimension bound.\n");
			}
			else if (err == NC_ERANGE)
			{
				printf("One Orange more of the values are out of range.\n");
			}
			else if (err == NC_EINDEFINE)
			{
				printf("Operation not allowed in define mode.\n");
			}
			else if (err == NC_EBADID)
			{
				printf("Bad ncid.\n");
			}
		}
		//nc_close(_file);
	}
	void write(int varIndex, double* data)
	{
		int err;
		//nc_open(_filename.data(), NC_WRITE | NC_SHARE, &_file);
		NetCDFVar var = _vars[varIndex];
		err = nc_inq_varid(_file, var.name.data(), &var.id);
		std::vector<size_t> startArr;
		std::vector<size_t> latlondim;
		for (int i = 0; i < var.ndims; i++)
		{
			startArr.push_back(0);
			latlondim.push_back(_dims[var.dimids[i]].len);
		}
		size_t* start = &startArr[0];
		size_t* count = &latlondim[0];
		err = nc_put_vara_double(_file, var.id, start, count, data);
	}
};

class Grid
{
public:
	Grid();
	Grid(const OGREnvelope& bb, const double& gridcellsize, int expandByCells = 1);
	Grid(double * adfGeoTransform, int xsize, int yxize);
	void fromFishnetShape(std::string fishnetfile);
	void fromFishnetRaster(std::string fishnetfile, bool load = false);
	static void toGridConfig(std::string fishnetfile, std::string configfile);
	void fromGridConfig(std::string configfile);
	~Grid();
	int getNumberOfValidCells();
	void reset(std::string attributeName = "ca11");
	void resetValue(std::string attributeName, double initValue=1);
	void reset(std::vector<std::string> dimensions);
	static Grid* createFishnet(const std::string& shapefile, double gridcellsize);
	OGRPolygon* toOGRPolygon(OGRLayer* layer, const OGREnvelope& bb);
	void normalizedByArea();
	void toShape(const std::string& outputfile, bool writeAttribute = false);
	void toShape(const std::string& outputfile, int cellID);
	void toShape(const std::string& outputfile, OGREnvelope destBound, bool writeAttribute = false);
	void toBoundaryShape(const std::string & outputfile, const std::string & georeffile);
	void toRaster(std::string filename, std::string wkt, double nodata = -9999);
	void toNetCDF(std::string filename, std::string spatialUnits);
	Grid* upscale(int scale);
	void gatherCells(ShapeFile* input, const char* attributeID = "ca11", double scale = 1.0);
	void intersectShape(std::string inshapefile, std::string outshapefile);
	//void gatherCells(const std::string& inputfile);
private:
	double calFraction(OGRFeature* fea, const int& footprintIdx);
public:
	double* cells;
	OGREnvelope bound;
	int nrows;
	int ncols;
	int slice;
	double nodata;
  std::string proj;
	double _adfGeoTransform[6];
	std::vector<std::string> dims;
};

struct LayerInfo
{
	std::string sectorName;
	std::string spatialFileName;
	std::string ffco2Field;
	bool isVulcan;
	bool useMask;
	std::string attributeFileName;
	LayerInfo()
	{

	}
	LayerInfo(std::string line)
	{
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		sectorName = splits[0];
		attributeFileName = splits[1];
		ffco2Field = splits[2];
		spatialFileName = splits[3];
		isVulcan = (atoi(splits[4].data()) == 1);
		useMask = (atoi(splits[5].data()) == 1);
	}
};

class  GridMaker
{
public:
	 GridMaker();
	~ GridMaker();
	std::vector<double> loadAttribute(std::string dbffile, std::string attribute);
	void makeGrid(std::string sectorCFGFile, std::string gridCFGFile, std::string spatialunits, std::string sourceDir, std::string griddedDir, std::string outDir);
};

