#include "Grid.h"
#include "SparseFractionGrid.h"
Grid::Grid(const OGREnvelope& bb, const double& gridcellsize, int expandByCells)
{
	nodata = FLT_MAX;
	//ShapeFile input(inputfile.data());
	ncols = (bb.MaxX - bb.MinX) / gridcellsize + expandByCells;
	if (gridcellsize * ncols < (bb.MaxX - bb.MinX))
	{
		ncols = ncols + 1;
	}

	nrows = (bb.MaxY - bb.MinY) / gridcellsize + expandByCells;

	//printf("%f,%f\n", bb.MaxX - bb.MinX, bb.MaxY - bb.MinY)
	if (gridcellsize * nrows < (bb.MaxY - bb.MinY))
	{
		nrows = nrows + 1;
	}

	OGREnvelope newbound;
	newbound.MinX = bb.MinX - ((ncols * gridcellsize) - (bb.MaxX - bb.MinX)) * 0.5;
	newbound.MaxX = bb.MaxX + ((ncols * gridcellsize) - (bb.MaxX - bb.MinX)) * 0.5;
	newbound.MinY = bb.MinY - ((nrows * gridcellsize) - (bb.MaxY - bb.MinY)) * 0.5;
	newbound.MaxY = bb.MaxY + ((nrows * gridcellsize) - (bb.MaxY - bb.MinY)) * 0.5;
	bound = newbound;
	_adfGeoTransform[0] = newbound.MinX;
	_adfGeoTransform[1] = gridcellsize;
	_adfGeoTransform[2] = 0;
	_adfGeoTransform[3] = newbound.MaxY;
	_adfGeoTransform[4] = 0;
	_adfGeoTransform[5] = -gridcellsize;
	cells = NULL;
	//reset();
}
Grid::Grid(double* adfGeoTransform,int xsize,int yxize)
{
	nodata = FLT_MAX;
	//ShapeFile input(inputfile.data());
	ncols = xsize;
	nrows = yxize;
	for (size_t i = 0; i < 6; i++)
	{
		_adfGeoTransform[i] = adfGeoTransform[i];
	}
	bound.MinX = _adfGeoTransform[0];
	bound.MaxX = bound.MinX + xsize*_adfGeoTransform[1];
	bound.MaxY = _adfGeoTransform[3];
	bound.MinY = bound.MaxY + yxize*_adfGeoTransform[5];
	cells = NULL;
	//reset();
}

void Grid::fromFishnetShape(std::string fishnetfile)
{
	nodata = FLT_MAX;
	ShapeFile input(fishnetfile.data());
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

	char wkt[512];
	char* pwkt = wkt;
	if (input.poLayer->GetSpatialRef())
	{
		input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		proj = pwkt;
	}


	double gridcellsize = (cellBB.MaxX - cellBB.MinX);
	ncols = (int)((bound.MaxX - bound.MinX + gridcellsize * 0.5) / gridcellsize);
	nrows = (int)((bound.MaxY - bound.MinY + gridcellsize * 0.5) / gridcellsize);

	_adfGeoTransform[0] = bound.MinX;
	_adfGeoTransform[1] = gridcellsize;
	_adfGeoTransform[2] = 0;
	_adfGeoTransform[3] = bound.MaxY;
	_adfGeoTransform[4] = 0;
	_adfGeoTransform[5] = -gridcellsize;
	cells = NULL;
}

void Grid::fromFishnetRaster(std::string fishnetfile,bool load)
{
	nodata = FLT_MAX;
	cells = NULL;
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDataset* pDataset = (GDALDataset*)GDALOpen(fishnetfile.data(), GA_ReadOnly);
	pDataset->GetGeoTransform(_adfGeoTransform);
	ncols = pDataset->GetRasterXSize();
	nrows = pDataset->GetRasterYSize();
	proj = pDataset->GetProjectionRef();
	if (load)
	{
		reset();
		pDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, ncols, nrows, cells, ncols, nrows, GDT_Float64, 0, 0);
		nodata = pDataset->GetRasterBand(1)->GetNoDataValue();
	}
	GDALClose((GDALDatasetH)pDataset);
	bound.MinX = _adfGeoTransform[0];
	bound.MaxY = _adfGeoTransform[3];
	bound.MaxX = _adfGeoTransform[0] + _adfGeoTransform[1] * ncols;
	bound.MinY = _adfGeoTransform[3] + _adfGeoTransform[5] * nrows;
	slice = ncols * nrows;
	
	//reset();
}

void Grid::toGridConfig(std::string fishnetfile, std::string configfile)
{
	double adfGeoTransform[6];
	GDALDataset* pDataset = (GDALDataset*)GDALOpen(fishnetfile.data(), GA_ReadOnly);
	pDataset->GetGeoTransform(adfGeoTransform);
	std::ofstream ofs;
	ofs.open(configfile);

	ofs << std::setprecision(10);
	ofs << std::fixed << "minx" << " " << adfGeoTransform[0] << std::endl;
	ofs << std::fixed << "miny" << " " << adfGeoTransform[3] + adfGeoTransform[5] * pDataset->GetRasterYSize() << std::endl;
	ofs << std::fixed << "maxx" << " " << adfGeoTransform[0] + adfGeoTransform[1] * pDataset->GetRasterXSize() << std::endl;
	ofs << std::fixed << "maxy" << " " << adfGeoTransform[3] << std::endl;
	ofs << std::fixed << "ncols" << " " << pDataset->GetRasterXSize() << std::endl;
	ofs << std::fixed << "nrows" << " " << pDataset->GetRasterYSize() << std::endl;
	ofs << std::fixed << "cellsize" << " " << adfGeoTransform[1] << std::endl;
	std::string proj = "";
	if (pDataset->GetProjectionRef())
		proj = pDataset->GetProjectionRef();
	ofs << "projection" << " " << proj;
	GDALClose((GDALDatasetH)pDataset);
	ofs.close();
}

void Grid::fromGridConfig(std::string configfile)
{
	nodata = FLT_MAX;
	cells = NULL;
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	std::ifstream ifs;
	ifs.open(configfile);
	std::string fieldname;
	ifs >> fieldname >> bound.MinX;
	ifs >> fieldname >> bound.MinY;
	ifs >> fieldname >> bound.MaxX;
	ifs >> fieldname >> bound.MaxY;
	ifs >> fieldname >> ncols;
	ifs >> fieldname >> nrows;
	ifs >> fieldname >> _adfGeoTransform[1];
	std::string line;
	while (line.size() < 10 || line.substr(0,10) != "projection")
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
	_adfGeoTransform[5] = -_adfGeoTransform[1];
	_adfGeoTransform[2] = 0;
	_adfGeoTransform[4] = 0;
	_adfGeoTransform[0] = bound.MaxY;
	_adfGeoTransform[3] = bound.MaxY;
	ifs.close();
	slice = ncols * nrows;
}

Grid::Grid()
{
	cells = NULL;
}


Grid::~Grid()
{
	if (cells)
		delete[] cells;
}
double Grid::calFraction(OGRFeature* fea, const int& footprintIdx)
{

	OGRGeometry* geo = fea->GetGeometryRef();
	OGRwkbGeometryType gtype = geo->getGeometryType();
	if (gtype == wkbLineString || gtype == wkbLineString25D || gtype == wkbMultiLineString)
	{
		double len1 = Utils::calPolylineLength(fea->GetGeometryRef());// *FOOTPRINT_SCALE_FACTOR;
		double len2 = fea->GetFieldAsDouble(footprintIdx);
		return len1 / len2;
	}
	else if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		double area1 = Utils::calPolygonArea(fea->GetGeometryRef());
		double area2 = fea->GetFieldAsDouble(footprintIdx);
		return area1 / area2;
	}

	return 1;
}

int Grid::getNumberOfValidCells()
{
	int count = 0;
	for (size_t i = 0; i < slice; i++)
	{
		double val = cells[i];
		if (val <= 0 || isinf(val) || val != val)
			continue;
		count++;
	}
	return count;
}
void Grid::reset(std::string attributeName)
{
	if (cells)
		delete[] cells;
	cells = new double[ncols*nrows];
	dims.clear();
	dims.push_back(attributeName);
	memset(cells, 0, sizeof(double)*ncols*nrows);
	slice = ncols*nrows;
}
void Grid::resetValue(std::string attributeName,double initValue)
{
	if (cells)
		delete[] cells;
	cells = new double[ncols*nrows];
	dims.clear();
	dims.push_back(attributeName);
	memset(cells, 0, sizeof(double)*ncols*nrows);
	slice = ncols*nrows;

	for (size_t i = 0; i < slice; i++)
	{
		cells[i] = initValue;
	}
}
void Grid::reset(std::vector<std::string> dimensions)
{
	dims = dimensions;
	if (cells)
		delete[] cells;
	cells = new double[ncols*nrows*dims.size()];
	memset(cells, 0, sizeof(double)*ncols*nrows*dims.size());
	slice = ncols*nrows;
}

Grid* Grid::createFishnet(const std::string& shapefile, double gridcellsize)
{
	ShapeFile input(shapefile.data());
	OGREnvelope bound;
	input.poLayer->GetExtent(&bound, true);
	Grid* grid = new Grid(bound, gridcellsize);

	return grid;
}

OGRPolygon* Grid::toOGRPolygon(OGRLayer* layer, const OGREnvelope& bb)
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
void Grid::normalizedByArea()
{
	double frac = _adfGeoTransform[1] * _adfGeoTransform[1];
	for (size_t i = 0; i < slice; i++)
	{
		cells[i] = cells[i] / frac;
	}
}

void Grid::toShape(const std::string& outputfile, bool writeAttribute)
{
	
	ShapeFile fishnet;
	if (this->proj != "")
	{
		OGRSpatialReference spatialref;
		const char* pProj = proj.data();
		spatialref.importFromWkt(&pProj);
		fishnet.create(outputfile.data(), &spatialref);
	}
  else
  {
    fishnet.create(outputfile.data());
  }

	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();
	double gridcellsize =  _adfGeoTransform[1];
	int idIdx =  fishnet.getOrCreateField("Id", OGRFieldType::OFTInteger);
	//int colIdx = fishnet.getOrCreateField("col", OGRFieldType::OFTInteger);
	//int rowIdx = fishnet.getOrCreateField("row", OGRFieldType::OFTInteger);

	//int caIndexOutput = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
	std::vector<int> fields;
	if (writeAttribute)
	{
		//if (caIndexOutput < 0)
		//{

		//OGRFieldDefn def("ca11", OGRFieldType::OFTReal);
		//fishnet.poLayer->CreateField(&def);
		//caIndexOutput = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
		//}

		for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
		{
			OGRFieldDefn def(dims[iattribute].data(), OGRFieldType::OFTReal);
			fishnet.poLayer->CreateField(&def);
			fields.push_back(fishnet.poLayer->GetLayerDefn()->GetFieldIndex(dims[iattribute].data()));
		}
	}

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

			OGREnvelope bb = tmpBound;
			bb.MinX = bound.MinX + gridcellsize * j;
			bb.MaxX = bound.MinX + gridcellsize * (j + 1);

		/*	if (bb.MinX >(-133.177071456802 - 0.1) && bb.MaxX < (-104.855640115368 + 0.1) && bb.MinY >(25.1208667507254 - 0.1) && bb.MaxY < (46.2050395140432 + 0.1))
			{*/
				OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());

				poFeaPolygon->SetField(idIdx, id);
	/*			poFeaPolygon->SetField(rowIdx, (int)i);
				poFeaPolygon->SetField(colIdx, (int)j);*/
				bool hasnonzero = false;
				if (writeAttribute)
				{
					//poFeaPolygon->SetField(caIndexOutput, cells[id]);
					for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
					{
						double val = cells[slice*iattribute + id];
						if (val > 0)
							hasnonzero = true;
						poFeaPolygon->SetField(fields[iattribute], cells[slice*iattribute + id]);
					}
				}
				OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bb);
				poFeaPolygon->SetGeometry(poPolygon);

				if (!writeAttribute || (writeAttribute && hasnonzero))
					fishnet.poLayer->CreateFeature(poFeaPolygon);

				OGRFeature::DestroyFeature(poFeaPolygon);
			//}
			
			id++;
		}

	}
	fishnet.close();
}
void Grid::toShape(const std::string& outputfile, OGREnvelope destBound, bool writeAttribute)
{

	ShapeFile fishnet;
	if (proj != "")
	{
		OGRSpatialReference spatialref;
		char* pProj = (char*)proj.data();
		spatialref.importFromWkt(&pProj);
		fishnet.create(outputfile.data(), &spatialref);
	}
	else
	{
		fishnet.create(outputfile.data());
	}

	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();
	double gridcellsize = _adfGeoTransform[1];
	//int idIdx = fishnet.getOrCreateField("Id", OGRFieldType::OFTInteger);
	std::vector<int> fields;
	if (writeAttribute)
	{
		for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
		{
			OGRFieldDefn def(dims[iattribute].data(), OGRFieldType::OFTReal);
			fishnet.poLayer->CreateField(&def);
			fields.push_back(fishnet.poLayer->GetLayerDefn()->GetFieldIndex(dims[iattribute].data()));
		}
	}

	int id = 0;
	int slice = nrows * ncols;
	for (size_t i = 0; i < nrows; i++)
	{
		OGREnvelope tmpBound;
		tmpBound.MaxY = bound.MaxY - gridcellsize * i;
		tmpBound.MinY = bound.MaxY - gridcellsize * (i + 1);
		for (size_t j = 0; j < ncols; j++)
		{

			OGREnvelope bb = tmpBound;
			bb.MinX = bound.MinX + gridcellsize * j;
			bb.MaxX = bound.MinX + gridcellsize * (j + 1);
			if (destBound.Intersects(bb) || destBound.Contains(bb))
			{
				OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());
				//poFeaPolygon->SetField(idIdx, id);
				bool hasnonzero = false;
				if (writeAttribute)
				{
					for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
					{
						double val = cells[slice*iattribute + id];
						poFeaPolygon->SetField(fields[iattribute], cells[slice*iattribute + id]);
					}
				}
				OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bb);
				poFeaPolygon->SetGeometry(poPolygon);

				fishnet.poLayer->CreateFeature(poFeaPolygon);
				OGRFeature::DestroyFeature(poFeaPolygon);
			}
			id++;
		}

	}
	fishnet.close();
}
void Grid::toShape(const std::string& outputfile, int cellID)
{

  ShapeFile fishnet;
  if (proj != "")
  {
    OGRSpatialReference spatialref;
    char* pProj = (char*)proj.data();
    spatialref.importFromWkt(&pProj);
    fishnet.create(outputfile.data(), &spatialref);
  }
  else
  {
    fishnet.create(outputfile.data());
  }


	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();
	double gridcellsize = _adfGeoTransform[1];
	int idIdx = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("Id");

	OGRFieldDefn def("Id", OGRFieldType::OFTInteger);
	fishnet.poLayer->CreateField(&def);
	idIdx = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("Id");

	int nrow = cellID / ncols;
	int ncol = cellID % ncols;
	int id = 0;

	OGREnvelope bb;
	bb.MaxY = bound.MaxY - gridcellsize * nrow;
	bb.MinY = bound.MaxY - gridcellsize * (nrow + 1);
	bb.MinX = bound.MinX + gridcellsize * ncol;
	bb.MaxX = bound.MinX + gridcellsize * (ncol + 1);

	OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());

	poFeaPolygon->SetField(idIdx, id);
	OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bb);
	poFeaPolygon->SetGeometry(poPolygon);

	fishnet.poLayer->CreateFeature(poFeaPolygon);

	OGRFeature::DestroyFeature(poFeaPolygon);
	fishnet.close();
}
void Grid::toBoundaryShape(const std::string& outputfile,const std::string& georeffile)
{

	ShapeFile fishnet;
	if (georeffile != "")
	{
		ShapeFile copyfrom(georeffile.data());
		fishnet.create(outputfile.data(), copyfrom.poLayer->GetSpatialRef());
		copyfrom.close();
	}
	else
	{
		fishnet.create(outputfile.data());
	}
	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();

 
	OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());

	OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bound);
	poFeaPolygon->SetGeometry(poPolygon);
	fishnet.poLayer->CreateFeature(poFeaPolygon);
	OGRFeature::DestroyFeature(poFeaPolygon);
	fishnet.close();


}

void Grid::toRaster(std::string filename, std::string wkt,double nodata)
{
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	//papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

	GDALDataset* pDataset = poDriver->Create(filename.data(), ncols, nrows, 1, GDT_Float64, papszOptions);

	pDataset->SetGeoTransform(_adfGeoTransform);
	double* pValue = cells;
	GDALRasterBand *pBand = pDataset->GetRasterBand(1);
	if (cells != NULL)
	{
		//for (size_t i = 0; i < ncols * nrows; i++)
		//{
		//	double val = cells[i];
		//	//if (val <= 0)
		//	//	cells[i] = nodata;
		//}
		pBand->RasterIO(GF_Write, 0, 0, ncols, nrows, cells, ncols, nrows, GDT_Float64, 0, 0);
	}
	pBand->SetNoDataValue(nodata);
	//if (wkt != NULL)
	//{
		pDataset->SetProjection(wkt.data());
	//}
	GDALClose((GDALDatasetH)pDataset);
}

void Grid::toNetCDF(std::string filename, std::string spatialUnits)
{
	//NCFile totalNCFile;
	//totalNCFile.create(filename.data(), ncols, nrows, 1, "year", spatialUnits);
	//double* xcoords = new double[ncols];
	//double* ycoords = new double[nrows];
	//double* tcoords = new double[numhours];
	//for (int i = 0; i < ncols; i++)
	//{
	//	xcoords[i] = _adfGeoTransform[0] + _adfGeoTransform[1] * i + +_adfGeoTransform[1] * 0.5;
	//}
	//for (int i = 0; i < nrows; i++)
	//{
	//	int idx = i;// nrows - i - 1;
	//	ycoords[i] = _adfGeoTransform[3] + _adfGeoTransform[5] * idx + +_adfGeoTransform[5] * 0.5;
	//}
	//for (int i = 0; i < numhours; i++)
	//{
	//	tcoords[i] = i;
	//}

	//double* sectorCells = new double[slice];

	//memset(cells, 0, slice*sizeof(double));
	//for (int isector = 0; isector < sectorGrids.size(); isector++)
	//{
	//	memset(sectorCells, 0, slice*sizeof(double));
	//	HestiaGrid* sectorGrid = sectorGrids[isector];
	//	sectorGrid->getTotal();
	//	addGridCells(sectorGrid->cells, cells, slice);
	//	//sectorGrid->initializeAnnualNetCDF()
	//	//sectorGrid->totalNCFile.writeSlice(0, sectorCells);
	//}

	//totalNCFile.writeSlice(0, cells);
	//totalNCFile.write(0, xcoords);
	//totalNCFile.write(1, ycoords);
	//totalNCFile.write(2, tcoords);
	//totalNCFile.close();

	//delete[] sectorCells;
	//delete[] xcoords;
	//delete[] ycoords;
	//delete[] tcoords;

}

//Grid* Grid::upscale(int scale)
//{
//	Grid* upscaledGrid = new Grid(bound, _adfGeoTransform[1] * scale, 0);
//	upscaledGrid->reset();
//	int idx = 0;
//	for (int nrow = 0; nrow < nrows; nrow++)
//	{
//		int nrow2 = nrow / scale;
//		for (int ncol = 0; ncol < ncols; ncol++)
//		{
//			int ncol2 = ncol / scale;
//			upscaledGrid->cells[ncol2 + nrow2*upscaledGrid->ncols] += cells[idx++];
//		}
//	}
//	return upscaledGrid;
//}
void Grid::gatherCells(ShapeFile* input,const char* attributeID, double scale)
{

	OGRFeature *poFeature;
	input->poLayer->ResetReading();
	int caIndexOld = input->poLayer->GetLayerDefn()->GetFieldIndex(attributeID);
	int idIndex = input->poLayer->GetLayerDefn()->GetFieldIndex("Id");
	if (caIndexOld == -1)
		return;
	OGRwkbGeometryType gtype = input->poLayer->GetGeomType();
	int footprintIndexOld = -1;
	int footprintIndexNew = -1;
	if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
	{
		footprintIndexOld = input->poLayer->GetLayerDefn()->GetFieldIndex("length");
	}
	if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		footprintIndexOld = input->poLayer->GetLayerDefn()->GetFieldIndex("area");
	}
	while ((poFeature = input->poLayer->GetNextFeature()) != NULL)
	{
		

		double fraction = calFraction(poFeature, footprintIndexOld) * scale;
		double ca = poFeature->GetFieldAsDouble(caIndexOld);

		if (ca > 0)
		{
			int id = poFeature->GetFieldAsInteger(idIndex);
			cells[id] += ca;// *fraction;

		}
		OGRFeature::DestroyFeature(poFeature);
	}

}

//void Grid::gatherCells(const std::string& inputfile)
//{
//	ShapeFile input(inputfile.data());
//	gatherCells(&input);
//
//}
#include "qdir.h"
GridMaker::GridMaker()
{

}
GridMaker::~GridMaker()
{

}
std::vector<double> GridMaker::loadAttribute(std::string dbffile,std::string attribute)
{
	std::vector<double> result;
	ShapeFile shp(dbffile);
	OGRFeature* poFeature;
	int fid = shp.poLayer->GetLayerDefn()->GetFieldIndex(attribute.data());
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		double val = poFeature->GetFieldAsDouble(fid);
		if (val != val || val < 0 || isinf((float)val) || !isfinite(val))
		{
			val = 0;
		}
		result.push_back(val);
		OGRFeature::DestroyFeature(poFeature);
	}
	return result;
}

void GridMaker::makeGrid(std::string sectorCFGFile, std::string fishnetRasterFile, std::string spatialunits, std::string sourceDir, std::string griddedDir, std::string outDir)
{
	QDir qoutdir(outDir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	outDir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	sourceDir = (QDir(sourceDir.data()).absolutePath() + "/").toLocal8Bit().data();
	griddedDir = (QDir(griddedDir.data()).absolutePath() + "/").toLocal8Bit().data();
	Grid totalGrid;
	totalGrid.fromFishnetRaster(fishnetRasterFile);
	totalGrid.reset();
	for (int icell = 0; icell < totalGrid.slice; icell++)
	{
		totalGrid.cells[icell] = _NODATA;
	}
	//std::string spatialUnits = "degrees";
	//OGRSpatialReference oSourceSRS;

	/*

	char* buf = new char[totalGrid.proj.size()];
	memcpy(buf, totalGrid.proj.data(), totalGrid.proj.size());
	oSourceSRS.importFromWkt(&buf);
	if (oSourceSRS.IsProjected())
	{
		char buf[255];
		char* pbuf = buf;
		oSourceSRS.GetLinearUnits(&pbuf);
		spatialUnits = pbuf;
	}*/
	std::ifstream ifs;
	ifs.open(sectorCFGFile);
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int idx = 0;
	std::map<std::string, std::vector<LayerInfo>> sectorMap;
	std::map<std::string, std::string> spatialMap;
	std::map<std::string, std::vector<LayerInfo>>::iterator iter;

	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		LayerInfo layerinfo(line);

		/*std::vector<std::string> splits = Utils::splitCSV(',', line);
		std::string sectorName = splits[0];
		std::string sectorFileName = splits[1];
		std::string ffco2Field = splits[2];
		bool isVulcan = (atoi(splits[4].data()) == 1);*/
		iter = sectorMap.find(layerinfo.sectorName);
		//bool isVulcan = fe_2007_us_zcta500

			//fe_2007_us_zcta500
		if (iter == sectorMap.end()){
			std::vector<LayerInfo> fileset;
			fileset.push_back(layerinfo);
			sectorMap[layerinfo.sectorName] = fileset;
		}
		else{
			iter->second.push_back(layerinfo);
		}
		spatialMap[layerinfo.attributeFileName] = layerinfo.spatialFileName;
	}

	double* xcoords = new double[totalGrid.ncols];
	double* ycoords = new double[totalGrid.nrows];

	for (int i = 0; i < totalGrid.ncols; i++)
	{
		xcoords[i] = totalGrid._adfGeoTransform[0] + totalGrid._adfGeoTransform[1] * i + totalGrid._adfGeoTransform[1] * 0.5;
	}
	for (int i = 0; i < totalGrid.nrows; i++)
	{
		int idx = i;// nrows - i - 1;
		ycoords[i] = totalGrid._adfGeoTransform[3] + totalGrid._adfGeoTransform[5] * idx + totalGrid._adfGeoTransform[5] * 0.5;
	}

	iter = sectorMap.begin();
	while (iter != sectorMap.end())
	{
		std::vector<LayerInfo> fileset = iter->second;
		Grid sectorGrid;
		sectorGrid.fromFishnetRaster(fishnetRasterFile);
		sectorGrid.reset();
		for (int icell = 0; icell < sectorGrid.slice; icell++)
		{
			sectorGrid.cells[icell] = _NODATA;
		}
		std::string outtif = outDir + iter->first + ".tif";
		std::string outnc = outDir + iter->first + ".nc";
		bool overwrite = false;
		if (!QFile(outtif.data()).exists() || overwrite)
		{
			for (size_t i = 0; i < fileset.size(); i++)
			{
				LayerInfo layer = fileset[i];
				std::string srcFile = sourceDir + layer.attributeFileName + ".dbf";
				std::string griddedFile = griddedDir + spatialMap[layer.attributeFileName] + ".bin";
				std::vector<double> srcMap = loadAttribute(srcFile, layer.ffco2Field);

				SparseFractionGrid sparseGrid;
				sparseGrid.open(griddedFile);
				for (size_t n = 0; n < sparseGrid.cells.size(); n++)
				{
					SparseGridCell& cell = sparseGrid.cells[n];
					int gridId = cell.gridId;
					int feaId = cell.feaId;
					double frac = cell.fraction;
					double emissions = srcMap[feaId] * frac;
					if (emissions > 0 && emissions < 1000000000)
					{
						if (sectorGrid.cells[gridId] == _NODATA)
							sectorGrid.cells[gridId] = 0;
						sectorGrid.cells[gridId] += emissions;
					}
				}
				//ShapeFile shp(griddedFile);
				//int gridField = shp.poLayer->GetLayerDefn()->GetFieldIndex("gridId");
				//int feaField = shp.poLayer->GetLayerDefn()->GetFieldIndex("feaId");
				//int fracField = shp.poLayer->GetLayerDefn()->GetFieldIndex("fraction");
				//OGRFeature* poFeature;
				//while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
				//{
				//	int gridId = poFeature->GetFieldAsInteger(gridField);
				//	int feaId = poFeature->GetFieldAsInteger(feaField);
				//	double frac = poFeature->GetFieldAsDouble(fracField);
				//	double emissions = srcMap[feaId] * frac;
				//	if (emissions > 0 && emissions < 1000000000)
				//	{
				//		if (sectorGrid.cells[gridId] == _NODATA)
				//			sectorGrid.cells[gridId] = 0;
				//		sectorGrid.cells[gridId] += emissions;
				//	}

				//	OGRFeature::DestroyFeature(poFeature);
				//}
				printf("%s,%s\n", layer.sectorName.data(), spatialMap[layer.attributeFileName].data());
			}
			////sectorGrid.toShape(sectorGrid.proj, outDir + iter->first + ".shp", true);
			//NCFile sectorNCFile;
			//sectorNCFile.create(outnc.data(), sectorGrid.ncols, sectorGrid.nrows, 1, "year", spatialunits);
			//sectorNCFile.writeSlice(0, sectorGrid.cells);
			//sectorNCFile.write(0, xcoords);
			//sectorNCFile.write(1, ycoords);
			////totalNCFile.write(2, tcoords);
			//sectorNCFile.close();
			sectorGrid.toRaster(outtif, sectorGrid.proj, _NODATA);
			if (fileset[0].isVulcan)
			{
				for (int icell = 0; icell < totalGrid.slice; icell++)
				{
					double cellval = sectorGrid.cells[icell];
					if (totalGrid.cells[icell] == _NODATA)
						totalGrid.cells[icell] = 0;
					if (cellval == _NODATA || cellval < 0 || cellval > 1000000000)
						continue;
					totalGrid.cells[icell] += cellval;
				}
			}


		}
		else
		{

			GDALDataset* pDataset = (GDALDataset*)GDALOpen(outtif.data(), GA_ReadOnly);
			pDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, totalGrid.ncols, totalGrid.nrows, sectorGrid.cells, totalGrid.ncols, totalGrid.nrows, GDT_Float64, 0, 0);
			GDALClose((GDALDatasetH)pDataset);
			if (fileset[0].isVulcan)
			{
				for (int icell = 0; icell < totalGrid.slice; icell++)
				{
					double cellval = sectorGrid.cells[icell];
					if (cellval == _NODATA || cellval < 0 || cellval > 1000000000)
						continue;
					if (totalGrid.cells[icell] == _NODATA)
						totalGrid.cells[icell] = 0;
					totalGrid.cells[icell] += cellval;
				}
			}

		}
		iter++;
	}

	//NCFile totalNCFile;
	//totalNCFile.create((outDir + "total.nc").data(), totalGrid.ncols, totalGrid.nrows, 1, "year", spatialunits);
	//totalNCFile.writeSlice(0, totalGrid.cells);
	//totalNCFile.write(0, xcoords);
	//totalNCFile.write(1, ycoords);
	////totalNCFile.write(2, tcoords);
	//totalNCFile.close();

	delete[] xcoords;
	delete[] ycoords;
	//delete[] tcoords;
	//totalGrid.toShape(totalGrid.proj, outDir + iter->first + ".shp", true);
	totalGrid.toRaster(outDir +  "Total.tif", totalGrid.proj, _NODATA);
}
