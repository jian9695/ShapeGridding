#pragma once
#include "ShapeFile.h"
#include "netcdf.h"
#include "qdir.h"
#include "qfileinfo.h"
#include <fstream>
#include <sstream>
#include "TimestructTool.h"
#include "GDAL_DS.h"
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
	void create(const char* outputFile, int ncols, int nrows, int timedimsize, std::string timeunit,std::string spatialunit = "feet")
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

		if (timedimsize > 0)
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
		if (timedimsize > 0)
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
		attr.text = "Kilogram";
		attr.vtype = NC_CHAR;
		cavar.attributes.push_back(attr);
		attr.name = "_FillValue";
		attr.value = 0;
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
		NetCDFVar var = _vars[3];
		err = nc_inq_varid(_file, var.name.data(), &var.id);
		std::vector<size_t> startArr;
		std::vector<size_t> latlondim;
		for (int i = 0; i < var.ndims; i++)
		{
			startArr.push_back(0);
			latlondim.push_back(_dims[var.dimids[i]].len);
		}
		startArr[0] = timeSlice;// *_dims[0].len * _dims[1].len;
		latlondim[0] = 1;// *_dims[0].len * _dims[1].len;
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
		//nc_close(_file);
		//if (err == NC_NOERR)
		//{
		//	printf("No error.\n");
		//}
		//else if (err == NC_ENOTVAR)
		//{
		//	printf("Variable not found.\n");
		//}
		//else if (err == NC_EINVALCOORDS)
		//{
		//	printf("Index exceeds dimension bound.\n");
		//}
		//else if (err == NC_EEDGE)
		//{
		//	printf("Start + count exceeds dimension bound.\n");
		//}
		//else if (err == NC_ERANGE)
		//{
		//	printf("One Orange more of the values are out of range.\n");
		//}
		//else if (err = NC_EINDEFINE)
		//{
		//	printf("Operation not allowed in define mode.\n");
		//}
		//else if (err ==  NC_EBADID)
		//{
		//	printf("Bad ncid.\n");
		//}

	}
};
//struct HourlyFractions
//{
//public:
//	//int num;
//	char id[20];
//	double fractions[8760];
//	std::string idStr()
//	{
//		return id;
//	}
//	HourlyFractions()
//	{
//
//	}
//	void set(double* _fractions, std::string _id)
//	{
//		//num = _num;
//		memset(id, 0, 20);
//		memcpy(id, _id.data(), _id.size());
//		memcpy(fractions, _fractions, 8760*sizeof(double));
//	}
//	void normalizeFractions()
//	{
//		double sum = 0;
//		for (size_t i = 0; i < 8760; i++)
//		{
//			sum += fractions[i];
//		}
//		for (size_t i = 0; i < 8760; i++)
//		{
//			fractions[i] = fractions[i] / sum;
//		}
//	}
//	~HourlyFractions()
//	{
//
//	}
//};
struct FeatureCell
{
	int gridid;
	int fid;
	double* timestructure;
	double total;
};



class HestiaGridBase
{
public:
	HestiaGridBase()
	{
		numdays = 365;
		cells = NULL;
	}

	~HestiaGridBase()
	{
		if (cells != NULL)
			delete[] cells;
		//hourlyNCFile.close();
		//totalNCFile.close();
	}

public:
	NCFile hourlyNCFile;
	NCFile totalNCFile;
	NCFile dailyNCFile;
	std::string sectorname;
	double* cells;
	std::map<std::string, double*> timestructMap;
	OGREnvelope bound;
	int nrows;
	int ncols;
	int slice;
	int numhours;
	int numdays;
	std::string spatialunit;
	double _adfGeoTransform[6];
	std::string fieldname;
	void initializeHourlyNetCDF(std::string hourlyFilename)
	{
		hourlyNCFile.create(hourlyFilename.data(), ncols, nrows, numhours, "hour", spatialunit);
	}
	void initializeDailyNetCDF(std::string dailyFilename)
	{
		dailyNCFile.create(dailyFilename.data(), ncols, nrows, numdays, "daily", spatialunit);
	}
	void initializeAnnualNetCDF(std::string totalFilename)
	{
		totalNCFile.create(totalFilename.data(), ncols, nrows, 1, "year", spatialunit);
	}
	void fromGrid(HestiaGridBase* grid)
	{
		if (cells != NULL)
			delete[] cells;
		double gridcellsize = grid->_adfGeoTransform[1];
		ncols = grid->ncols;
		nrows = grid->nrows;
		slice = nrows * ncols;
		for (size_t i = 0; i < 6; i++)
		{
			_adfGeoTransform[i] = grid->_adfGeoTransform[i];
		}
		bound = grid->bound;
		cells = new double[slice];
		memset(cells, 0, slice*sizeof(double));
	}

	void fromFishnetShapeFile(std::string shapefile)
	{
		if (cells != NULL)
			delete[] cells;
		ShapeFile input(shapefile.data());
		input.poLayer->GetExtent(&bound);
		OGRFeature *poFeature;
		input.poLayer->ResetReading();
		OGREnvelope cellBB;
		while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
		{

			((OGRPolygon*)poFeature->GetGeometryRef())->getEnvelope(&cellBB);
			OGRFeature::DestroyFeature(poFeature);
			break;
		}
		double gridcellsize = (cellBB.MaxX - cellBB.MinX);
		ncols = (int)((bound.MaxX - bound.MinX + gridcellsize * 0.5) / gridcellsize);
		nrows = (int)((bound.MaxY - bound.MinY + gridcellsize * 0.5) / gridcellsize);
		slice = nrows * ncols;
		_adfGeoTransform[0] = bound.MinX;
		_adfGeoTransform[1] = gridcellsize;
		_adfGeoTransform[2] = 0;
		_adfGeoTransform[3] = bound.MaxY;
		_adfGeoTransform[4] = 0;
		_adfGeoTransform[5] = -gridcellsize;
		cells = new double[slice];
		memset(cells, 0, slice*sizeof(double));
	}
	void fromFishnetRaster(std::string rasterfile)
	{
		if (cells != NULL)
			delete[] cells;
		const char *pszFormat = "GTiff";
		char **papszOptions = NULL;
		GDALDataset* pDataset = (GDALDataset*)GDALOpen(rasterfile.data(), GA_ReadOnly);
		pDataset->GetGeoTransform(_adfGeoTransform);
		ncols = pDataset->GetRasterXSize();
		nrows = pDataset->GetRasterYSize();
		slice = nrows * ncols;
		GDALClose((GDALDatasetH)pDataset);
		bound.MinX = _adfGeoTransform[0];
		bound.MaxY = _adfGeoTransform[3];
		bound.MaxX = _adfGeoTransform[0] + _adfGeoTransform[1] * ncols;
		bound.MinY = _adfGeoTransform[3] + _adfGeoTransform[5] * nrows;
		cells = new double[slice];
		memset(cells, 0, slice*sizeof(double));
	}
	OGRPolygon* toOGRPolygon(OGRLayer* layer, const OGREnvelope& bb)
	{

		OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
		OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);
		linearRing->addPoint(bb.MinX, bb.MinY);
		linearRing->addPoint(bb.MinX, bb.MaxY);
		linearRing->addPoint(bb.MaxX, bb.MaxY);
		linearRing->addPoint(bb.MaxX, bb.MinY);
		linearRing->addPoint(bb.MinX, bb.MinY);

		poPolygon->addRing(linearRing);//also crashed
		return poPolygon;
	}
	void toShapefile(std::string spatialRefShape, std::string outputfile)
	{
		ShapeFile fishnet;
		if (spatialRefShape != "")
		{
			ShapeFile copyfrom(spatialRefShape.data());
			fishnet.create(outputfile.data(), copyfrom.poLayer->GetSpatialRef());
			copyfrom.close();
		}
		else
		{
			fishnet.create(outputfile.data());
		}

		OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();
		double gridcellsize = _adfGeoTransform[1];
		int idIdx = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("Id");
		if (idIdx < 0)
		{
			OGRFieldDefn def("Id", OGRFieldType::OFTInteger);
			fishnet.poLayer->CreateField(&def);
			idIdx = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("Id");
		}

		OGRFieldDefn def(fieldname.data(), OGRFieldType::OFTReal);
		fishnet.poLayer->CreateField(&def);
		int attributeIdx = fishnet.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data());
		int id = 0;
		int slice = nrows * ncols;
		for (size_t i = 0; i < nrows; i++)
		{
			//std::vector<GridCell>& newrow = gridcells[i];
			OGREnvelope tmpBound;
			tmpBound.MaxY = bound.MaxY - gridcellsize * i;
			tmpBound.MinY = bound.MaxY - gridcellsize * (i + 1);
			for (size_t j = 0; j < ncols; j++)
			{
				OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());
				OGREnvelope bb = tmpBound;
				bb.MinX = bound.MinX + gridcellsize * j;
				bb.MaxX = bound.MinX + gridcellsize * (j + 1);
				poFeaPolygon->SetField(idIdx, id);
				double val = cells[id];
				if (val > 0)
					poFeaPolygon->SetField(attributeIdx, cells[id]);
				OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bb);
				poFeaPolygon->SetGeometry(poPolygon);
				if (val > 0)
					fishnet.poLayer->CreateFeature(poFeaPolygon);
				OGRFeature::DestroyFeature(poFeaPolygon);
				id++;
			}

		}
		fishnet.close();
	}

};

class HestiaGrid : public HestiaGridBase
{
public:
	HestiaGrid()
		:HestiaGridBase()
	{

	}

public:

	std::vector<FeatureCell> features;
	std::vector<std::string> shapefiles;
	//void loadSectors(std::vector<std::string> sectorShapefiles, std::string emissionFieldName = "ca11", std::string timestructFieldName = "timestruct")
	//{
	//	fieldname = emissionFieldName;
	//	for (size_t i = 0; i < sectorShapefiles.size(); i++)
	//	{
	//		loadFeatures(sectorShapefiles[i],emissionFieldName,timestructFieldName);
	//	}
	//}
	virtual void loadAttribute(std::string _fieldname)
	{
		fieldname = _fieldname;
		features.clear();
		for (int ishapefile = 0; ishapefile < shapefiles.size(); ishapefile++)
		{
			std::string shapefile = shapefiles[ishapefile];
			ShapeFile input(shapefile.data());
			OGRFeature *poFeature;
			input.poLayer->ResetReading();

			int caIndex = input.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data());
			int idIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("Id");
	
			int timestructIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("timestruct");
			//if (!(caIndex == -1 || idIndex == -1 || timestructIndex == -1))
			//{
			//	printf("");
			//	continue
			//}
			if(timestructIndex < 0)
			   timestructIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("bt");

			bool isfloat = false;
			if (timestructIndex > -1 && input.poLayer->GetLayerDefn()->GetFieldDefn(timestructIndex)->GetType() == OGRFieldType::OFTReal)
				isfloat = true;
			OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
			while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
			{
				double ca = poFeature->GetFieldAsDouble(caIndex);
				if (ca > 0)
				{
					int id = poFeature->GetFieldAsInteger(idIndex);
					std::string tsid;
					if (timestructIndex > -1)
					{
						if (!isfloat)
						{
							tsid = poFeature->GetFieldAsString(timestructIndex);
						}
						else
						{
							std::stringstream ss;
							ss << (long)poFeature->GetFieldAsDouble(timestructIndex);
							tsid = ss.str();

						}
					}

		

					FeatureCell fc;
					fc.total = ca;
					fc.gridid = id;
					if (timestructIndex > -1)
					{
						if (timestructMap.find(tsid) != timestructMap.end())
						{
							fc.timestructure = timestructMap[tsid];

						}
						else
						{
							fc.timestructure = timestructMap["0"];
							if (timestructMap.size() > 2)
							{
								printf("%s\n", tsid.data());
							}

						}
					}
					features.push_back(fc);
				}
				OGRFeature::DestroyFeature(poFeature);
			}
		}
	}
	//void loadFeatures(std::string shapefile,std::string emissionFieldName = "ca11",std::string timestructFieldName = "timestruct")
	//{
	//	shapefiles.push_back(shapefile);
	//	ShapeFile input(shapefile.data());
	//	input.poLayer->GetExtent(&bound);
	//	OGRFeature *poFeature;
	//	input.poLayer->ResetReading();

	//	int caIndex = input.poLayer->GetLayerDefn()->GetFieldIndex(emissionFieldName.data());
	//	int idIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("Id");
	//	int timestructIndex = input.poLayer->GetLayerDefn()->GetFieldIndex(timestructFieldName.data());
	//	OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
	//	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	//	{
	//		double ca = poFeature->GetFieldAsDouble(caIndex);
	//		if (ca > 0)
	//		{
	//			int id = poFeature->GetFieldAsInteger(idIndex);
	//			std::string tsid = poFeature->GetFieldAsString(timestructIndex);
	//			FeatureCell fc;
	//			fc.total = ca;
	//			fc.gridid = id;
	//			if (timestructMap.find(tsid) == timestructMap.end())
	//			{
	//				printf("%s\n", tsid.data());
	//			}
	//			fc.timestructure = timestructMap[tsid];
	//			features.push_back(fc);
	//		}
	//		OGRFeature::DestroyFeature(poFeature);
	//	}
	//}
	virtual void getTimeSlice(double* gridcells, int houridx)
	{
		//memset(gridcells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			if (fc.total <= 0)
				continue;
			//printf("%d\n", fc.gridid);
			gridcells[fc.gridid] += (fc.total * fc.timestructure[houridx]);
		}
	}
	virtual void getTotal(double* gridcells)
	{
		memset(gridcells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			gridcells[fc.gridid] += fc.total;
		}
	}
	virtual void getTimeSlice(int houridx)
	{
		memset(cells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			if (fc.total <= 0)
				continue;
			//printf("%d\n", fc.gridid);
			cells[fc.gridid] += (fc.total * fc.timestructure[houridx]);
		}
	}
	virtual void getTotal()
	{
		memset(cells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			if (fc.gridid < 0 || fc.gridid > slice)
			{
				printf("");
			}
			cells[fc.gridid] += fc.total;
		}
	}

};

class FFDASGrid : public HestiaGrid
{
public:
	FFDASGrid()
		:HestiaGrid()
	{

	}

public:
	double* ffdasAnnual;
	std::vector<float*> ffdasHourly;
	virtual void loadAttribute(std::string _fieldname)
	{
		_fieldname = "area";
		fieldname = _fieldname;
		features.clear();
		for (int ishapefile = 0; ishapefile < shapefiles.size(); ishapefile++)
		{
			std::string shapefile = shapefiles[ishapefile];
			ShapeFile input(shapefile.data());
			OGRFeature *poFeature;
			input.poLayer->ResetReading();

			int caIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("area");
			int fidIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("Id");
			int gridIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("FID_fishne");
			OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
			while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
			{
				double ca = poFeature->GetFieldAsDouble(caIndex);
				if (ca > 0)
				{
					FeatureCell fc;
					fc.gridid = poFeature->GetFieldAsInteger(gridIndex);
					fc.fid = poFeature->GetFieldAsInteger(fidIndex);
					fc.total = ca;
					features.push_back(fc);
				}
				OGRFeature::DestroyFeature(poFeature);
			}
		}
	}

	virtual void getTimeSlice(double* gridcells, int houridx)
	{
		//memset(gridcells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			if (fc.total <= 0)
				continue;
			//printf("%d\n", fc.gridid);
			gridcells[fc.gridid] += (fc.total * ffdasHourly[houridx][fc.fid]);
		}
	}
	virtual void getTotal(double* gridcells)
	{
		memset(gridcells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			gridcells[fc.gridid] += (fc.total * ffdasAnnual[fc.fid]);
		}
	}
	virtual void getTimeSlice(int houridx)
	{
		memset(cells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			if (fc.total <= 0)
				continue;
			//printf("%d\n", fc.gridid);
			cells[fc.gridid] += (fc.total * ffdasHourly[houridx][fc.fid]);
		}
	}
	virtual void getTotal()
	{
		memset(cells, 0, sizeof(double)*slice);
		for (size_t i = 0; i < features.size(); i++)
		{
			FeatureCell& fc = features[i];
			cells[fc.gridid] += (fc.total * ffdasAnnual[fc.fid]);
			//printf("%d,%f,%f\n", fc.fid, fc.total, ffdasAnnual[fc.fid]);
		}
	}

};

class TemporalGridder : public HestiaGridBase
{
public:
	char* idbuf;
	double* fractionbuf;
	std::vector<HestiaGrid*> sectorGrids;
	std::string  outdir;
	TemporalGridder(std::string _outdir,std::string _spatialunitname = "feet")
		:HestiaGridBase()
	{
		outdir = _outdir;
		QDir qoutdir(outdir.data());
		//if (!qoutdir.exists())
		//	qoutdir.mkpath(".");
		outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
		spatialunit = _spatialunitname;
		idbuf = NULL;
		fractionbuf = NULL;
		numhours = 8760;
		numdays = 365;
	}
	~TemporalGridder()
	{
		//if (cells != NULL)
		//	delete[] cells;
		timestructMap.clear();
		if (idbuf)
		{
			delete[] idbuf;
			delete[] fractionbuf;
		}
		clearSectorGrids();
	}
	void normalizeFractions(double*& fractions)
	{
		double sum = 0;
		for (size_t i = 0; i < 8760; i++)
		{
			sum += fractions[i];
		}
		for (size_t i = 0; i < 8760; i++)
		{
			fractions[i] = fractions[i] / sum;
		}
	}
	void clearSectorGrids()
	{
		for (size_t i = 0; i < sectorGrids.size(); i++)
		{
			delete sectorGrids[i];
		}
		sectorGrids.clear();
	}
	void addSectorGrid(std::vector<std::string> shapefiles,std::string sectorname)
	{
		HestiaGrid* sectorGrid = new HestiaGrid;
		sectorGrid->fromGrid(this);
		sectorGrid->timestructMap = timestructMap;
		sectorGrid->shapefiles = shapefiles;
		sectorGrid->sectorname = sectorname;
		//sectorGrid->initializeNetCDF(numhours, outdir + totalFilename, outdir + hourlyFilename, spatialunit);
		sectorGrid->spatialunit = spatialunit;
		sectorGrid->numhours = numhours;
		sectorGrid->numdays = numdays;
		sectorGrids.push_back(sectorGrid);
	}

	void addFFDASGrid(std::vector<std::string> shapefiles, std::string sectorname, double* ffdasAnnual, std::vector<float*>& ffdasHourly)
	{
		FFDASGrid* sectorGrid = new FFDASGrid;
		sectorGrid->ffdasAnnual = ffdasAnnual;
		sectorGrid->ffdasHourly = ffdasHourly;
		sectorGrid->fromGrid(this);
		sectorGrid->timestructMap = timestructMap;
		sectorGrid->shapefiles = shapefiles;
		sectorGrid->sectorname = sectorname;
		//sectorGrid->initializeNetCDF(numhours, outdir + totalFilename, outdir + hourlyFilename, spatialunit);
		sectorGrid->spatialunit = spatialunit;
		sectorGrid->numhours = numhours;
		sectorGrid->numdays = numdays;
		sectorGrids.push_back(sectorGrid);
	}

	void loadAttribute(std::string fieldname)
	{
		this->fieldname = fieldname;
		for (int isector = 0; isector < sectorGrids.size(); isector++)
		{
			HestiaGrid* sectorGrid = sectorGrids[isector];
			sectorGrid->loadAttribute(fieldname);
		}
	}
	void loadtimestruct(std::string txtdir, std::string binaryfile)
	{
		if (!QFileInfo(binaryfile.data()).exists())
		{
			TimestructTool::txt2binary(txtdir, binaryfile,true);
		}
		
		if (idbuf)
		{
			delete[] idbuf;
			delete[] fractionbuf;
			timestructMap.clear();
		}
		std::ifstream filein(binaryfile.data(), std::ios::binary);
		int numstructs = 0;
		filein.read((char*)&numstructs, sizeof(int) * 1);
		idbuf = new char[numstructs * 20];
		memset(idbuf,0, (size_t)numstructs * 20);
		fractionbuf = new double[numstructs * 8760];
		filein.read(idbuf, numstructs * 20);
		filein.read((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * 8760);
		filein.close();
		//std::ofstream ofs("B:/LA_Version2/Time/2011.csv");
		char* pidbuf = idbuf;
		double* pfractionbuf = fractionbuf;
		//std::string outputdir = "C:/HestiaGridding/Baltimore/";
		for (size_t i = 0; i < numstructs; i++)
		{
			//ofs << pidbuf << std::endl;
			std::string tsid = pidbuf;
			/*double sum = 0;
			std::ofstream ofsts((outputdir + tsid + ".csv").data());
			for (size_t ihour = 0; ihour < 8760; ihour++)
			{
			sum += pfractionbuf[ihour];
			ofsts << pfractionbuf[ihour] << std::endl;
			}
			ofsts.close();
			printf("%f\n",sum);*/

			double sum = 0;
			
			for (size_t ihour = 0; ihour < 8760; ihour++)
			{
			    sum += pfractionbuf[ihour];
			}
			if (sum == 0 || (sum - 1) > 0.000001)
			{
				if (sum == 0)
				{
					double frac = 1.0 / 8760.0;
					for (size_t ihour = 0; ihour < 8760; ihour++)
					{
						pfractionbuf[ihour] = frac;
					}
				}
				else
				{
					for (size_t ihour = 0; ihour < 8760; ihour++)
					{
						pfractionbuf[ihour] = pfractionbuf[ihour] / sum;
					}
				}
				printf("%f\n", sum);
			}

			if (tsid[0] == '0' && tsid.length() > 1)
				tsid = tsid.substr(1, tsid.length() - 1);
			timestructMap[tsid] = pfractionbuf;
			pidbuf += 20;
			pfractionbuf += 8760;

		}
		//ofs.close();

	}
	void addGridCells(double* grid1, double* grid2, int numcells)
	{
		for (size_t i = 0; i < numcells; i++)
		{
			grid2[i] = grid1[i] + grid2[i];
		}
	}
	void makeHourlyTotal(std::string hourlyFilename)
	{
		if (timestructMap.size() == 0)
		{
			printf("Time structures not loaded!");
			return;
		}

		initializeHourlyNetCDF(hourlyFilename);

		double* xcoords = new double[ncols];
		double* ycoords = new double[nrows];
		double* tcoords = new double[numhours];
		for (int i = 0; i < ncols; i++)
		{
			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < nrows; i++)
		{
			int idx = i;// nrows - i - 1;
			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < numhours; i++)
		{
			tcoords[i] = i;
		}

		double* sectorCells = new double[slice];
		double sum = 0;
		for (int ihour = 0; ihour < numhours; ihour++)
		{
			memset(cells, 0, slice*sizeof(double));
			for (int isector = 0; isector < sectorGrids.size(); isector++)
			{
				memset(sectorCells, 0, slice*sizeof(double));
				HestiaGrid* sectorGrid = sectorGrids[isector];
				sectorGrid->getTimeSlice(sectorCells, ihour);
				addGridCells(sectorCells, cells,slice);
				//sectorGrid->hourlyNCFile.writeSlice(ihour, sectorCells);
			}
			hourlyNCFile.writeSlice(ihour, cells);
			for (size_t icell = 0; icell < slice; icell++)
			{
				sum += cells[icell];
			}
			if (ihour % 100 == 0)
			{
				printf("%d/%d\n", ihour, 8760);
			}
		}
		printf("total=%f\n", sum);
		hourlyNCFile.write(0, xcoords);
		hourlyNCFile.write(1, ycoords);
		hourlyNCFile.write(2, tcoords);
		hourlyNCFile.close();
		//for (int isector = 0; isector < sectorGrids.size(); isector++)
		//{
		//	HestiaGrid* sectorGrid = sectorGrids[isector];
		//	sectorGrid->hourlyNCFile.write(0, xcoords);
		//	sectorGrid->hourlyNCFile.write(1, ycoords);
		//	sectorGrid->hourlyNCFile.write(2, tcoords);
		//	sectorGrid->hourlyNCFile.close();
		//}
		for (int isector = 0; isector < sectorGrids.size(); isector++)
		{
			HestiaGrid* sectorGrid = sectorGrids[isector];
			memset(sectorGrid->cells, 0, slice*sizeof(double));
		}
		memset(cells, 0, slice*sizeof(double));
		delete[] sectorCells;
		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;

	}
	void makeDailyTotal(std::string dailyFilename)
	{
		if (timestructMap.size() == 0)
		{
			printf("Time structures not loaded!");
			return;
		}

		initializeDailyNetCDF(dailyFilename);

		double* xcoords = new double[ncols];
		double* ycoords = new double[nrows];
		double* tcoords = new double[numdays];
		for (int i = 0; i < ncols; i++)
		{
			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < nrows; i++)
		{
			int idx = i;// nrows - i - 1;
			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < numdays; i++)
		{
			tcoords[i] = i;
		}

		double* sectorCells = new double[slice];
		double sum = 0;
		int ihour = 0;
		for (int iday = 0; iday < numdays; iday++)
		{
			memset(cells, 0, slice*sizeof(double));
			for (size_t hourofday = 0; hourofday < 24; hourofday++)
			{
				for (int isector = 0; isector < sectorGrids.size(); isector++)
				{
					memset(sectorCells, 0, slice*sizeof(double));
					HestiaGrid* sectorGrid = sectorGrids[isector];
					sectorGrid->getTimeSlice(sectorCells, ihour);
					addGridCells(sectorCells, cells, slice);
				}
				ihour++;
			}
			dailyNCFile.writeSlice(iday, cells);
			for (size_t icell = 0; icell < slice; icell++)
			{
				sum += cells[icell];
			}
			printf("%d/%d\n", ihour, 8760);
		}
		printf("total=%f\n", sum);
		dailyNCFile.write(0, xcoords);
		dailyNCFile.write(1, ycoords);
		dailyNCFile.write(2, tcoords);
		dailyNCFile.close();
		//for (int isector = 0; isector < sectorGrids.size(); isector++)
		//{
		//	HestiaGrid* sectorGrid = sectorGrids[isector];
		//	sectorGrid->hourlyNCFile.write(0, xcoords);
		//	sectorGrid->hourlyNCFile.write(1, ycoords);
		//	sectorGrid->hourlyNCFile.write(2, tcoords);
		//	sectorGrid->hourlyNCFile.close();
		//}
		for (int isector = 0; isector < sectorGrids.size(); isector++)
		{
			HestiaGrid* sectorGrid = sectorGrids[isector];
			memset(sectorGrid->cells, 0, slice*sizeof(double));
		}
		memset(cells, 0, slice*sizeof(double));
		delete[] sectorCells;
		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;

	}
	void makeAnnualTotal(std::string totalFilename)
	{
		//if (timestructMap.size() == 0)
		//{
		//	printf("Time structures not loaded!");
		//	return;
		//}

		initializeAnnualNetCDF(totalFilename);

		double* xcoords = new double[ncols];
		double* ycoords = new double[nrows];
		double* tcoords = new double[numhours];
		for (int i = 0; i < ncols; i++)
		{
			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < nrows; i++)
		{
			int idx = i;// nrows - i - 1;
			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < numhours; i++)
		{
			tcoords[i] = i;
		}

		double* sectorCells = new double[slice];

		memset(cells, 0, slice*sizeof(double));
		for (int isector = 0; isector < sectorGrids.size(); isector++)
		{
			memset(sectorCells, 0, slice*sizeof(double));
			HestiaGrid* sectorGrid = sectorGrids[isector];
			sectorGrid->getTotal();
			addGridCells(sectorGrid->cells, cells, slice);
			//sectorGrid->initializeAnnualNetCDF()
			//sectorGrid->totalNCFile.writeSlice(0, sectorCells);
		}

		totalNCFile.writeSlice(0, cells);
		totalNCFile.write(0, xcoords);
		totalNCFile.write(1, ycoords);
		totalNCFile.write(2, tcoords);
		totalNCFile.close();
		//for (int isector = 0; isector < sectorGrids.size(); isector++)
		//{
		//	HestiaGrid* sectorGrid = sectorGrids[isector];

		//	sectorGrid->totalNCFile.write(0, xcoords);
		//	sectorGrid->totalNCFile.write(1, ycoords);
		//	sectorGrid->totalNCFile.write(2, tcoords);


		//	//sectorGrid->hourlyNCFile.write(0, xcoords);
		//	//sectorGrid->hourlyNCFile.write(1, ycoords);
		//	//sectorGrid->hourlyNCFile.write(2, tcoords);

		//	sectorGrid->totalNCFile.close();
		//	//sectorGrid->hourlyNCFile.close();
		//}

		delete[] sectorCells;
		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;

	}
	void makeAnnualTotal(std::string outdir,std::string cityname,std::string year,std::string version)
	{
		//if (timestructMap.size() == 0)
		//{
		//	printf("Time structures not loaded!");
		//	return;
		//}
		//gridder.makeAnnualTotal(outroot + cityname + ".total.annual." + year + "." + version + ".nc");
		initializeAnnualNetCDF(outdir + cityname + ".total.annual." + year + "." + version + ".nc");

		double* xcoords = new double[ncols];
		double* ycoords = new double[nrows];
		double* tcoords = new double[numhours];
		for (int i = 0; i < ncols; i++)
		{
			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < nrows; i++)
		{
			int idx = i;// nrows - i - 1;
			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < numhours; i++)
		{
			tcoords[i] = i;
		}

		//double* sectorCells = new double[slice];

		memset(cells, 0, slice*sizeof(double));
		for (int isector = 0; isector < sectorGrids.size(); isector++)
		{
			//memset(sectorCells, 0, slice*sizeof(double));
			HestiaGrid* sectorGrid = sectorGrids[isector];
			sectorGrid->getTotal();
			addGridCells(sectorGrid->cells, cells, slice);
			sectorGrid->initializeAnnualNetCDF(outdir + cityname + "." + sectorGrid->sectorname + +".annual." + year + "." + version + ".nc");
			sectorGrid->totalNCFile.writeSlice(0, sectorGrid->cells);
			sectorGrid->totalNCFile.write(0, xcoords);
			sectorGrid->totalNCFile.write(1, ycoords);
			sectorGrid->totalNCFile.write(2, tcoords);

			sectorGrid->totalNCFile.close();
		}

		totalNCFile.writeSlice(0, cells);
		totalNCFile.write(0, xcoords);
		totalNCFile.write(1, ycoords);
		totalNCFile.write(2, tcoords);
		totalNCFile.close();
		//for (int isector = 0; isector < sectorGrids.size(); isector++)
		//{
		//	HestiaGrid* sectorGrid = sectorGrids[isector];

		//	sectorGrid->totalNCFile.write(0, xcoords);
		//	sectorGrid->totalNCFile.write(1, ycoords);
		//	sectorGrid->totalNCFile.write(2, tcoords);


		//	//sectorGrid->hourlyNCFile.write(0, xcoords);
		//	//sectorGrid->hourlyNCFile.write(1, ycoords);
		//	//sectorGrid->hourlyNCFile.write(2, tcoords);

		//	sectorGrid->totalNCFile.close();
		//	//sectorGrid->hourlyNCFile.close();
		//}

		//delete[] sectorCells;
		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;

	}

	void tif2netcdf(std::vector<std::string> tifFiles, std::string netcdfFilename)
	{

		initializeAnnualNetCDF(netcdfFilename);

		double* xcoords = new double[ncols];
		double* ycoords = new double[nrows];
		double* tcoords = new double[numhours];
		for (int i = 0; i < ncols; i++)
		{
			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < nrows; i++)
		{
			int idx = i;// nrows - i - 1;
			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < numhours; i++)
		{
			tcoords[i] = i;
		}

		memset(cells, 0, slice*sizeof(double));
		for (int itifFile = 0; itifFile < tifFiles.size(); itifFile++)
		{
			GDAL_DS<double>* ds = new GDAL_DS<double>();
			ds->open(tifFiles[itifFile]);
			double* data = ds->readData(1);
			addGridCells(data, cells,  slice);
			delete[] data;
			delete ds;
		}

		totalNCFile.writeSlice(0, cells);
		totalNCFile.write(0, xcoords);
		totalNCFile.write(1, ycoords);
		totalNCFile.write(2, tcoords);
		totalNCFile.close();

		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;

	}
	void tif2netcdf(std::string tifFile, std::string netcdfFilename)
	{
		std::vector<std::string> tifFiles;
		tifFiles.push_back(tifFile);
		tif2netcdf(tifFiles, netcdfFilename);
	}
	NCFile* createNC(const char* outputFile, int ncols,int nrows,int timedimsize)
	{

		std::vector<NetCDFDim> dims;
		NetCDFDim londim;
		londim.len = ncols;
		londim.name = "X";
		dims.push_back(londim);

		NetCDFDim latdim;
		latdim.len = nrows;
		latdim.name = "Y";
		dims.push_back(latdim);

		if (timedimsize > 0)
		{
			NetCDFDim timedim;
			timedim.len = timedimsize;
			timedim.name = "Time";
			dims.push_back(timedim);
		}
	/*	NetCDFAttribute attr;
		attr.len = 10;
		std::vector<NetCDFVar> vars;
		NetCDFVar lonvar;
		lonvar.dimids[0] = 0;
		lonvar.ndims = 1;
		lonvar.name = "Lon";
		lonvar.vtype = NC_DOUBLE;
		lonvar.varnatt = 2;
		attr.name = "units";
		attr.text = "degree";
		attr.vtype = 2;
		lonvar.attributes.push_back(attr);
		attr.name = "long_name";
		attr.text = "Lon";
		attr.vtype = 2;
		lonvar.attributes.push_back(attr);
		vars.push_back(lonvar);

		NetCDFVar latvar;
		latvar.dimids[0] = 1;
		latvar.ndims = 1;
		latvar.name = "Lat";
		latvar.vtype = NC_DOUBLE;
		latvar.varnatt = 2;
		attr.name = "units";
		attr.text = "degree";
		attr.vtype = 2;
		latvar.attributes.push_back(attr);
		attr.name = "long_name";
		attr.text = "Lat";
		attr.vtype = 2;
		latvar.attributes.push_back(attr);
		vars.push_back(latvar);*/
		NetCDFAttribute attr;
		attr.len = 10;
		std::vector<NetCDFVar> vars;
		NetCDFVar lonvar;
		lonvar.dimids[0] = 0;
		lonvar.ndims = 1;
		lonvar.name = "X";
		lonvar.vtype = NC_DOUBLE;
		lonvar.varnatt = 2;
		attr.name = "units";
		attr.text = "feet";
		attr.vtype = 2;
		lonvar.attributes.push_back(attr);
		attr.name = "X";
		attr.text = "X";
		attr.vtype = 2;
		lonvar.attributes.push_back(attr);
		vars.push_back(lonvar);

		NetCDFVar latvar;
		latvar.dimids[0] = 1;
		latvar.ndims = 1;
		latvar.name = "Y";
		latvar.vtype = NC_DOUBLE;
		latvar.varnatt = 2;
		attr.name = "units";
		attr.text = "feet";
		attr.vtype = 2;
		latvar.attributes.push_back(attr);
		attr.name = "Y";
		attr.text = "Y";
		attr.vtype = 2;
		latvar.attributes.push_back(attr);
		vars.push_back(latvar);

		NetCDFVar cavar;
		if (timedimsize > 0)
		{

			NetCDFVar timevar;
			timevar.dimids[0] = 2;
			timevar.ndims = 1;
			timevar.name = "Time";
			timevar.vtype = NC_INT;
			timevar.varnatt = 2;
			attr.name = "units";
			attr.text = "hour";
			attr.vtype = NC_CHAR;
			timevar.attributes.push_back(attr);
			attr.name = "long_name";
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
		attr.text = "Kilogram";
		attr.vtype = NC_CHAR;
		cavar.attributes.push_back(attr);
		attr.name = "_FillValue";
		attr.value = NAN;
		attr.vtype = NC_DOUBLE;
		cavar.attributes.push_back(attr);
		vars.push_back(cavar);

		NCFile* ncfile = new NCFile;
		ncfile->_dims = dims;
		ncfile->_vars = vars;
		ncfile->_filename = outputFile;
		int status = nc_create(outputFile, NC_CLOBBER, &(ncfile->_file));

		//size_t start[] = { 0, 0 }; // start at first value 
		std::vector<size_t> latlondim;
		for (int i = 0; i < dims.size(); i++)
		{
			latlondim.push_back(dims[i].len);
			nc_def_dim(ncfile->_file, dims[i].name.data(), dims[i].len, &i);
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
			nc_def_var(ncfile->_file, var.name.data(), var.vtype, var.ndims, var.dimids, &var.id);
			for (int j = 0; j < var.attributes.size(); j++)
			{
				NetCDFAttribute attr = var.attributes[j];
				if (attr.vtype == NC_DOUBLE)
					nc_put_att_double(ncfile->_file, var.id, attr.name.data(), NC_DOUBLE, 1, &(attr.value));
				else
					nc_put_att_text(ncfile->_file, var.id, attr.name.data(), attr.text.size(), attr.text.data());
			}
		}
		nc_close(ncfile->_file);
		//NetCDFVar refvar = vars[vars.size() - 1];

		//std::vector<size_t> dimlens;
		//for (size_t j = 0; j < refvar.ndims; j++)
		//{
		//	dimlens.push_back(dims[refvar.dimids[j]].len);
		//}

		//size_t ncells = nrows*ncols;
		//for (size_t i = 0; i < 24; i++)
		//{
		//	std::stringstream ss;
		//	if (i > 8)
		//		ss << refvar.name << "_h" << i + 1;
		//	else
		//		ss << refvar.name << "_h0" << i + 1;
		//	nc_def_var(ncidDest, ss.str().data(), refvar.vtype, refvar.ndims, refvar.dimids, &refvar.id);
		//	NetCDFAttribute attr = refvar.attributes[0];
		//	nc_put_att_text(ncidDest, refvar.id, attr.name, attr.len, attr.text);
		//	//size_t offset = ncells * i;
		//	//nc_put_vara_float(ncidSrc, refvar.id, start, count, data);
		//}
		/*if (nc_enddef(ncidDest) != NC_NOERR) return;*/
		//delete[] data;
		//nc_close(ncidDest);
		return ncfile;


	}
	void getTotal()
	{
		double* xcoords = new double[ncols];
		double* ycoords = new double[nrows];
		double* tcoords = new double[numhours];
		for (int i = 0; i < ncols; i++)
		{
			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < nrows; i++)
		{
			int idx = i;// nrows - i - 1;
			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < numhours; i++)
		{
			tcoords[i] = i;
		}

		double* sectorCells = new double[slice];
		memset(cells, 0, slice*sizeof(double));
		for (int isector = 0; isector < sectorGrids.size(); isector++)
		{
			memset(sectorCells, 0, slice*sizeof(double));
			HestiaGrid* sectorGrid = sectorGrids[isector];
			sectorGrid->getTotal();
			addGridCells(sectorGrid->cells, cells, slice);
		}

		delete[] sectorCells;
		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;
	}
};


//class TemporalGridderWithFFDAS : public HestiaGridBase
//{
//public:
//	char* idbuf;
//	double* fractionbuf;
//	std::vector<HestiaGrid*> sectorGrids;
//	std::string outdir;
//	double* ffdasAnnual;
//	std::vector<float*> ffdasHourly;
//	TemporalGridderWithFFDAS(std::string _outdir, std::string _spatialunitname = "feet")
//		:HestiaGridBase()
//	{
//		outdir = _outdir;
//		QDir qoutdir(outdir.data());
//		if (!qoutdir.exists())
//			qoutdir.mkpath(".");
//		outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
//		spatialunit = _spatialunitname;
//		idbuf = NULL;
//		fractionbuf = NULL;
//		numhours = 8760;
//	}
//	~TemporalGridderWithFFDAS()
//	{
//		//if (cells != NULL)
//		//	delete[] cells;
//		timestructMap.clear();
//		if (idbuf)
//		{
//			delete[] idbuf;
//			delete[] fractionbuf;
//		}
//		clearSectorGrids();
//	}
//	void normalizeFractions(double*& fractions)
//	{
//		double sum = 0;
//		for (size_t i = 0; i < 8760; i++)
//		{
//			sum += fractions[i];
//		}
//		for (size_t i = 0; i < 8760; i++)
//		{
//			fractions[i] = fractions[i] / sum;
//		}
//	}
//	void clearSectorGrids()
//	{
//		for (size_t i = 0; i < sectorGrids.size(); i++)
//		{
//			delete sectorGrids[i];
//		}
//		sectorGrids.clear();
//	}
//	void addSectorGrid(std::vector<std::string> shapefiles, std::string sectorname)
//	{
//		//if (timestructMap.size() == 0)
//		//{
//		//	printf("Time structures not loaded!");
//		//	return;
//		//}
//		HestiaGrid* sectorGrid = new HestiaGrid;
//		sectorGrid->fromGrid(this);
//		sectorGrid->timestructMap = timestructMap;
//		sectorGrid->shapefiles = shapefiles;
//		sectorGrid->sectorname = sectorname;
//		//sectorGrid->initializeNetCDF(numhours, outdir + totalFilename, outdir + hourlyFilename, spatialunit);
//		sectorGrid->spatialunit = spatialunit;
//		sectorGrid->numhours = numhours;
//		sectorGrids.push_back(sectorGrid);
//	}
//	void loadAttribute(std::string fieldname)
//	{
//		this->fieldname = fieldname;
//		for (int isector = 0; isector < sectorGrids.size(); isector++)
//		{
//			HestiaGrid* sectorGrid = sectorGrids[isector];
//			sectorGrid->loadAttribute(fieldname);
//		}
//	}
//	void loadtimestruct(std::string txtdir, std::string binaryfile)
//	{
//		if (!QFileInfo(binaryfile.data()).exists())
//		{
//			TimestructTool::txt2binary(txtdir, binaryfile);
//		}
//
//		if (idbuf)
//		{
//			delete[] idbuf;
//			delete[] fractionbuf;
//			timestructMap.clear();
//		}
//		std::ifstream filein(binaryfile.data(), std::ios::binary);
//		int numstructs = 0;
//		filein.read((char*)&numstructs, sizeof(int) * 1);
//		idbuf = new char[numstructs * 20];
//		memset(idbuf, 0, (size_t)numstructs * 20);
//		fractionbuf = new double[numstructs * 8760];
//		filein.read(idbuf, numstructs * 20);
//		filein.read((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * 8760);
//		filein.close();
//		//std::ofstream ofs("B:/LA_Version2/Time/2011.csv");
//		char* pidbuf = idbuf;
//		double* pfractionbuf = fractionbuf;
//		//std::string outputdir = "C:/HestiaGridding/Baltimore/";
//		for (size_t i = 0; i < numstructs; i++)
//		{
//			//ofs << pidbuf << std::endl;
//			std::string tsid = pidbuf;
//			/*double sum = 0;
//			std::ofstream ofsts((outputdir + tsid + ".csv").data());
//			for (size_t ihour = 0; ihour < 8760; ihour++)
//			{
//			sum += pfractionbuf[ihour];
//			ofsts << pfractionbuf[ihour] << std::endl;
//			}
//			ofsts.close();
//			printf("%f\n",sum);*/
//
//			double sum = 0;
//
//			for (size_t ihour = 0; ihour < 8760; ihour++)
//			{
//				sum += pfractionbuf[ihour];
//			}
//			if ((sum - 1) > 0.0001)
//				printf("%f\n", sum);
//			timestructMap[tsid] = pfractionbuf;
//			pidbuf += 20;
//			pfractionbuf += 8760;
//
//		}
//		//ofs.close();
//
//	}
//	void addGridCells(double* grid1, double* grid2, int numcells)
//	{
//		for (size_t i = 0; i < numcells; i++)
//		{
//			grid2[i] = grid1[i] + grid2[i];
//		}
//	}
//	void makeHourlyTotal(std::string hourlyFilename)
//	{
//		if (timestructMap.size() == 0)
//		{
//			printf("Time structures not loaded!");
//			return;
//		}
//
//		initializeHourlyNetCDF(outdir + hourlyFilename);
//
//		double* xcoords = new double[ncols];
//		double* ycoords = new double[nrows];
//		double* tcoords = new double[numhours];
//		for (int i = 0; i < ncols; i++)
//		{
//			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
//		}
//		for (int i = 0; i < nrows; i++)
//		{
//			int idx = i;// nrows - i - 1;
//			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
//		}
//		for (int i = 0; i < numhours; i++)
//		{
//			tcoords[i] = i;
//		}
//
//		double* sectorCells = new double[slice];
//		double sum = 0;
//		for (int ihour = 0; ihour < numhours; ihour++)
//		{
//			memset(cells, 0, slice*sizeof(double));
//			for (int isector = 0; isector < sectorGrids.size(); isector++)
//			{
//				memset(sectorCells, 0, slice*sizeof(double));
//				HestiaGrid* sectorGrid = sectorGrids[isector];
//				sectorGrid->getTimeSlice(sectorCells, ihour);
//				addGridCells(sectorCells, cells, slice);
//				//sectorGrid->hourlyNCFile.writeSlice(ihour, sectorCells);
//			}
//			hourlyNCFile.writeSlice(ihour, cells);
//			for (size_t icell = 0; icell < slice; icell++)
//			{
//				sum += cells[icell];
//			}
//			if (ihour % 100 == 0)
//			{
//				printf("%d/%d\n", ihour, 8760);
//			}
//		}
//		printf("total=%f\n", sum);
//		hourlyNCFile.write(0, xcoords);
//		hourlyNCFile.write(1, ycoords);
//		hourlyNCFile.write(2, tcoords);
//		hourlyNCFile.close();
//		//for (int isector = 0; isector < sectorGrids.size(); isector++)
//		//{
//		//	HestiaGrid* sectorGrid = sectorGrids[isector];
//		//	sectorGrid->hourlyNCFile.write(0, xcoords);
//		//	sectorGrid->hourlyNCFile.write(1, ycoords);
//		//	sectorGrid->hourlyNCFile.write(2, tcoords);
//		//	sectorGrid->hourlyNCFile.close();
//		//}
//		for (int isector = 0; isector < sectorGrids.size(); isector++)
//		{
//			HestiaGrid* sectorGrid = sectorGrids[isector];
//			memset(sectorGrid->cells, 0, slice*sizeof(double));
//		}
//		memset(cells, 0, slice*sizeof(double));
//		delete[] sectorCells;
//		delete[] xcoords;
//		delete[] ycoords;
//		delete[] tcoords;
//
//	}
//	void makeAnnualTotal(std::string totalFilename)
//	{
//		//if (timestructMap.size() == 0)
//		//{
//		//	printf("Time structures not loaded!");
//		//	return;
//		//}
//
//		initializeAnnualNetCDF(outdir + totalFilename);
//
//		double* xcoords = new double[ncols];
//		double* ycoords = new double[nrows];
//		double* tcoords = new double[numhours];
//		for (int i = 0; i < ncols; i++)
//		{
//			xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
//		}
//		for (int i = 0; i < nrows; i++)
//		{
//			int idx = i;// nrows - i - 1;
//			ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
//		}
//		for (int i = 0; i < numhours; i++)
//		{
//			tcoords[i] = i;
//		}
//
//		double* sectorCells = new double[slice];
//
//		memset(cells, 0, slice*sizeof(double));
//		for (int isector = 0; isector < sectorGrids.size(); isector++)
//		{
//			memset(sectorCells, 0, slice*sizeof(double));
//			HestiaGrid* sectorGrid = sectorGrids[isector];
//			sectorGrid->getTotal();
//			addGridCells(sectorGrid->cells, cells, slice);
//			//sectorGrid->totalNCFile.writeSlice(0, sectorCells);
//		}
//
//		totalNCFile.writeSlice(0, cells);
//		totalNCFile.write(0, xcoords);
//		totalNCFile.write(1, ycoords);
//		totalNCFile.write(2, tcoords);
//		totalNCFile.close();
//		//for (int isector = 0; isector < sectorGrids.size(); isector++)
//		//{
//		//	HestiaGrid* sectorGrid = sectorGrids[isector];
//
//		//	sectorGrid->totalNCFile.write(0, xcoords);
//		//	sectorGrid->totalNCFile.write(1, ycoords);
//		//	sectorGrid->totalNCFile.write(2, tcoords);
//
//
//		//	sectorGrid->hourlyNCFile.write(0, xcoords);
//		//	sectorGrid->hourlyNCFile.write(1, ycoords);
//		//	sectorGrid->hourlyNCFile.write(2, tcoords);
//
//		//	sectorGrid->totalNCFile.close();
//		//	sectorGrid->hourlyNCFile.close();
//		//}
//
//		delete[] sectorCells;
//		delete[] xcoords;
//		delete[] ycoords;
//		delete[] tcoords;
//
//	}
//	NCFile* createNC(const char* outputFile, int ncols, int nrows, int timedimsize)
//	{
//
//		std::vector<NetCDFDim> dims;
//		NetCDFDim londim;
//		londim.len = ncols;
//		londim.name = "X";
//		dims.push_back(londim);
//
//		NetCDFDim latdim;
//		latdim.len = nrows;
//		latdim.name = "Y";
//		dims.push_back(latdim);
//
//		if (timedimsize > 0)
//		{
//			NetCDFDim timedim;
//			timedim.len = timedimsize;
//			timedim.name = "Time";
//			dims.push_back(timedim);
//		}
//		/*	NetCDFAttribute attr;
//		attr.len = 10;
//		std::vector<NetCDFVar> vars;
//		NetCDFVar lonvar;
//		lonvar.dimids[0] = 0;
//		lonvar.ndims = 1;
//		lonvar.name = "Lon";
//		lonvar.vtype = NC_DOUBLE;
//		lonvar.varnatt = 2;
//		attr.name = "units";
//		attr.text = "degree";
//		attr.vtype = 2;
//		lonvar.attributes.push_back(attr);
//		attr.name = "long_name";
//		attr.text = "Lon";
//		attr.vtype = 2;
//		lonvar.attributes.push_back(attr);
//		vars.push_back(lonvar);
//
//		NetCDFVar latvar;
//		latvar.dimids[0] = 1;
//		latvar.ndims = 1;
//		latvar.name = "Lat";
//		latvar.vtype = NC_DOUBLE;
//		latvar.varnatt = 2;
//		attr.name = "units";
//		attr.text = "degree";
//		attr.vtype = 2;
//		latvar.attributes.push_back(attr);
//		attr.name = "long_name";
//		attr.text = "Lat";
//		attr.vtype = 2;
//		latvar.attributes.push_back(attr);
//		vars.push_back(latvar);*/
//		NetCDFAttribute attr;
//		attr.len = 10;
//		std::vector<NetCDFVar> vars;
//		NetCDFVar lonvar;
//		lonvar.dimids[0] = 0;
//		lonvar.ndims = 1;
//		lonvar.name = "X";
//		lonvar.vtype = NC_DOUBLE;
//		lonvar.varnatt = 2;
//		attr.name = "units";
//		attr.text = "feet";
//		attr.vtype = 2;
//		lonvar.attributes.push_back(attr);
//		attr.name = "X";
//		attr.text = "X";
//		attr.vtype = 2;
//		lonvar.attributes.push_back(attr);
//		vars.push_back(lonvar);
//
//		NetCDFVar latvar;
//		latvar.dimids[0] = 1;
//		latvar.ndims = 1;
//		latvar.name = "Y";
//		latvar.vtype = NC_DOUBLE;
//		latvar.varnatt = 2;
//		attr.name = "units";
//		attr.text = "feet";
//		attr.vtype = 2;
//		latvar.attributes.push_back(attr);
//		attr.name = "Y";
//		attr.text = "Y";
//		attr.vtype = 2;
//		latvar.attributes.push_back(attr);
//		vars.push_back(latvar);
//
//		NetCDFVar cavar;
//		if (timedimsize > 0)
//		{
//
//			NetCDFVar timevar;
//			timevar.dimids[0] = 2;
//			timevar.ndims = 1;
//			timevar.name = "Time";
//			timevar.vtype = NC_INT;
//			timevar.varnatt = 2;
//			attr.name = "units";
//			attr.text = "hour";
//			attr.vtype = NC_CHAR;
//			timevar.attributes.push_back(attr);
//			attr.name = "long_name";
//			attr.text = "Time";
//			attr.vtype = NC_CHAR;
//			timevar.attributes.push_back(attr);
//			vars.push_back(timevar);
//
//
//			cavar.dimids[0] = 2;
//			cavar.dimids[1] = 1;
//			cavar.dimids[2] = 0;
//			cavar.ndims = 3;
//		}
//		else
//		{
//			cavar.dimids[0] = 1;
//			cavar.dimids[1] = 0;
//			cavar.ndims = 2;
//		}
//		cavar.name = "Carbon Emission";
//		cavar.vtype = NC_DOUBLE;
//		cavar.varnatt = 2;
//		attr.name = "units";
//		attr.text = "Kilogram";
//		attr.vtype = NC_CHAR;
//		cavar.attributes.push_back(attr);
//		attr.name = "_FillValue";
//		attr.value = NAN;
//		attr.vtype = NC_DOUBLE;
//		cavar.attributes.push_back(attr);
//		vars.push_back(cavar);
//
//		NCFile* ncfile = new NCFile;
//		ncfile->_dims = dims;
//		ncfile->_vars = vars;
//		ncfile->_filename = outputFile;
//		int status = nc_create(outputFile, NC_CLOBBER, &(ncfile->_file));
//
//		//size_t start[] = { 0, 0 }; // start at first value 
//		std::vector<size_t> latlondim;
//		for (int i = 0; i < dims.size(); i++)
//		{
//			latlondim.push_back(dims[i].len);
//			nc_def_dim(ncfile->_file, dims[i].name.data(), dims[i].len, &i);
//		}
//		size_t* count = &latlondim[0];
//
//		for (int i = 0; i < vars.size(); i++)
//		{
//			NetCDFVar var = vars[i];
//
//
//			std::vector<size_t> dimlens;
//			for (size_t j = 0; j < var.ndims; j++)
//			{
//				dimlens.push_back(dims[var.dimids[j]].len);
//			}
//			nc_def_var(ncfile->_file, var.name.data(), var.vtype, var.ndims, var.dimids, &var.id);
//			for (int j = 0; j < var.attributes.size(); j++)
//			{
//				NetCDFAttribute attr = var.attributes[j];
//				if (attr.vtype == NC_DOUBLE)
//					nc_put_att_double(ncfile->_file, var.id, attr.name.data(), NC_DOUBLE, 1, &(attr.value));
//				else
//					nc_put_att_text(ncfile->_file, var.id, attr.name.data(), attr.text.size(), attr.text.data());
//			}
//		}
//		nc_close(ncfile->_file);
//		//NetCDFVar refvar = vars[vars.size() - 1];
//
//		//std::vector<size_t> dimlens;
//		//for (size_t j = 0; j < refvar.ndims; j++)
//		//{
//		//	dimlens.push_back(dims[refvar.dimids[j]].len);
//		//}
//
//		//size_t ncells = nrows*ncols;
//		//for (size_t i = 0; i < 24; i++)
//		//{
//		//	std::stringstream ss;
//		//	if (i > 8)
//		//		ss << refvar.name << "_h" << i + 1;
//		//	else
//		//		ss << refvar.name << "_h0" << i + 1;
//		//	nc_def_var(ncidDest, ss.str().data(), refvar.vtype, refvar.ndims, refvar.dimids, &refvar.id);
//		//	NetCDFAttribute attr = refvar.attributes[0];
//		//	nc_put_att_text(ncidDest, refvar.id, attr.name, attr.len, attr.text);
//		//	//size_t offset = ncells * i;
//		//	//nc_put_vara_float(ncidSrc, refvar.id, start, count, data);
//		//}
//		/*if (nc_enddef(ncidDest) != NC_NOERR) return;*/
//		//delete[] data;
//		//nc_close(ncidDest);
//		return ncfile;
//
//
//	}
//
//};

