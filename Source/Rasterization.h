#pragma once

#include <QtGlobal>
#include <QApplication>
#include <QDesktopWidget>
#include <QSplashScreen>
#include <QTextCodec>
#include <QFileinfo>
#include "qdir.h"
#include "qpainter.h"
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <gdal_priv.h>
#include "ShapeFile.h"
#include "Grid.h"
#include "GDAL_DS.h"
#include "SparseFractionGrid.h"
#include "ogrsf_frmts.h"

struct FeatureFootprintCount
{
  int fid;
  int cellid;
  double frac;
  double totalarea;
  FeatureFootprintCount()
  {

  }
  FeatureFootprintCount(int _fid, int _cellid, double _area)
  {
    fid = _fid; cellid = _cellid; frac = _area; totalarea = _area;
  }

};

class RasterMask
{

public:
  QImage * _image;
  double _adfGeoTransform[6];
  int _xsize;
  int _ysize;
  uchar _nodata;
  int _stride;
  RasterMask(QImage* img, double xmin, double ymin, double xmax, double ymax)
  {
    _image = img;
    _adfGeoTransform[0] = xmin;
    _adfGeoTransform[1] = (xmax - xmin) / (img->width());
    _adfGeoTransform[2] = 0;
    _adfGeoTransform[3] = ymax;
    _adfGeoTransform[4] = 0;
    _adfGeoTransform[5] = -(ymax - ymin) / (img->height());
    _nodata = 0;
    _xsize = img->width();
    _ysize = img->height();
    _stride = 4;
  }
  ~RasterMask()
  {
    delete _image;
  }
  void save(std::string filename)
  {
    GDALAllRegister();
    _image->save((filename + ".png").data(), "PNG");
    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    char **papszMetadata;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    GDALDataset *poSrcDS = (GDALDataset *)GDALOpen((filename + ".png").data(), GA_ReadOnly);
    GDALDataset *poDstDS = poDriver->CreateCopy(filename.data(), poSrcDS, FALSE,
      NULL, NULL, NULL);
    poDstDS->SetGeoTransform(_adfGeoTransform);

    GDALClose((GDALDatasetH)poSrcDS);
    GDALClose((GDALDatasetH)poDstDS);
  }
  //adfGeoTransform[0] /* top left x */
  //adfGeoTransform[1] /* w-e pixel resolution */
  //adfGeoTransform[2] /* 0 */
  //adfGeoTransform[3] /* top left y */
  //adfGeoTransform[4] /* 0 */
  //adfGeoTransform[5] /* n-s pixel resolution (negative value) */
};

class TileRenderer
{
public:
  OGREnvelope envelope;
  double topleftx;
  double toplefty;
  int tilesizex;
  int tilesizey;
  int tilecols;
  int tilerows;
  double cellsize;
  int  ncols;
  int  nrows;
  int numtiles;
  void init(const double _cellsize, OGREnvelope _envelope)
  {
    envelope = _envelope;
    cellsize = _cellsize;
    ncols = (int)(ceil((envelope.MaxX - envelope.MinX) / cellsize));
    nrows = (int)(ceil((envelope.MaxY - envelope.MinY) / cellsize));
    while (ncols > 8192 * 5 || nrows >  8192 * 5)
    {
      cellsize = cellsize * 2;
      ncols = (int)((envelope.MaxX - envelope.MinX) / cellsize);
      nrows = (int)((envelope.MaxY - envelope.MinY) / cellsize);
    }
    int minsize = 10;
    int maxsize = 8192;

    while (ncols < 100 || nrows < 100)
    {
      cellsize = cellsize / 2;
      ncols = (int)((envelope.MaxX - envelope.MinX) / cellsize);
      nrows = (int)((envelope.MaxY - envelope.MinY) / cellsize);
    }
    tilesizex = maxsize;
    tilesizey = maxsize;
    if (ncols <= maxsize && nrows <= maxsize)
    {
      tilesizex = ncols;
      tilesizey = nrows;
    }
    tilecols = (int)(ceil((double)ncols / (double)tilesizex));
    tilerows = (int)(ceil((double)nrows / (double)tilesizey));
    numtiles = tilecols * tilerows;
    topleftx = envelope.MinX;
    toplefty = envelope.MaxY;

  }

  void drawPolygon(QPainter* painter, const QBrush& brush, const int& xsize, const int& ysize, const OGRPolygon *polygon)
  {
    QPainterPath outerpath;
    QPainterPath innerpath;
    const OGRLinearRing* exteriorRing = polygon->getExteriorRing();
    QPolygonF outerPolygon;
    for (size_t i = 0; i < exteriorRing->getNumPoints(); i++)
    {
      OGRPoint pt;
      exteriorRing->getPoint(i, &pt);
      outerPolygon.push_back(QPointF(pt.getX(), pt.getY()));
      //printf("%f,%f\n", pt.getX(), pt.getY());
    }
    outerpath.addPolygon(outerPolygon);
    for (size_t n = 0; n < polygon->getNumInteriorRings(); n++)
    {
      QPolygonF innerPolygon;
      const OGRLinearRing* interiorRing = polygon->getInteriorRing(n);
      for (size_t i = 0; i < interiorRing->getNumPoints(); i++)
      {
        OGRPoint pt;
        interiorRing->getPoint(i, &pt);
        innerPolygon.push_back(QPointF(pt.getX(), pt.getY()));
      }
      innerpath.addPolygon(innerPolygon);
      //path.addPolygon(outerPolygon);
    }
    if (polygon->getNumInteriorRings())
      outerpath = outerpath.subtracted(innerpath);
    painter->fillPath(outerpath, brush);

  }

  void drawPolyline(QPainter* painter, const QBrush& brush, const int& xsize, const int& ysize, const OGRLineString *polyline)
  {
    QVector<QPointF> qpoints;
    painter->setPen(Qt::black);
    for (size_t i = 0; i < polyline->getNumPoints() - 1; i++)
    {
      OGRPoint pt;
      polyline->getPoint(i, &pt);
      qpoints.push_back(QPointF(pt.getX(), pt.getY()));
      polyline->getPoint(i + 1, &pt);
      qpoints.push_back(QPointF(pt.getX(), pt.getY()));
    }

    painter->drawLines(qpoints);

  }

  RasterMask* drawTile(int tileid, OGRGeometry* geom, bool ispolygon = false)
  {

    int tiley = tileid / tilecols;
    int tilex = tileid % tilecols;

    if (tileid > numtiles - 1)
      return NULL;
    OGREnvelope tilebound;
    tilebound.MinX = topleftx + tilex * tilesizex * cellsize;
    tilebound.MaxY = toplefty - tiley * tilesizey * cellsize;
    tilebound.MaxX = tilebound.MinX + tilesizex * cellsize;
    tilebound.MinY = tilebound.MaxY - tilesizey * cellsize;
    QImage* qImage = new QImage(tilesizex, tilesizey, QImage::Format_ARGB32);
    qImage->fill(qRgba(0, 0, 0, 0));
    QBrush br(Qt::white);
    QPainter painter(qImage);
    painter.begin(qImage);
    //painter.setRenderHint(QPainter::Antialiasing); //make it look nicer
    painter.setRenderHint(QPainter::Antialiasing); //make it look nicer
                                                   //Antialiasing = 0x01,
                                                   //TextAntialiasing = 0x02,
                                                   //SmoothPixmapTransform = 0x04,
                                                   //HighQualityAntialiasing = 0x08,
                                                   //NonCosmeticDefaultPen = 0x10

    QTransform transform = QTransform().translate(-tilebound.MinX, -tilebound.MaxY) * QTransform().scale((float)tilesizex / (tilebound.MaxX - tilebound.MinX), -(float)tilesizey / (tilebound.MaxY - tilebound.MinY));
    painter.setTransform(transform);
    if (ispolygon)
    {
      painter.setBrush(br);
      drawPolygon(&painter, br, tilesizex, tilesizey, (OGRPolygon*)geom);
      painter.end();
    }
    else
    {
      drawPolyline(&painter, br, tilesizex, tilesizey, (OGRLineString*)geom);
      painter.end();
    }
    RasterMask* mask = new RasterMask(qImage, tilebound.MinX, tilebound.MinY, tilebound.MaxX, tilebound.MaxY);
    return mask;
  }
};

static double RASTERIZATION_RESOLUTION; // = 10 * 3.28084

class Rasterization
{
public:
  Rasterization();
  ~Rasterization();

  static void intersectRaster(std::string gridfile, std::string inputfile, std::string outputfile, bool usedMask = true)
  {
    Grid mask;
    mask.fromFishnetRaster(gridfile, true);
    Grid grid;
    grid.fromFishnetRaster(gridfile);
    double cellsize = grid._adfGeoTransform[1];
    GDAL_DS<double>* ds = new GDAL_DS<double>();
    ds->open(inputfile);
    double gridcellsize = ds->adfGeoTransform[1];
    //RASTERIZATION_RESOLUTION = gridcellsize / 50;
    OGRSpatialReference oTargetSRS;
    char* buf = new char[grid.proj.size()];
    memcpy(buf, grid.proj.data(), grid.proj.size());
    oTargetSRS.importFromWkt(&buf);
    OGRSpatialReference oSourceSRS;
    buf = new char[ds->projection.size()];
    memcpy(buf, ds->projection.data(), ds->projection.size());
    oSourceSRS.importFromWkt(&buf);
    double* cellvalues = ds->readData(1);
    double* pcellvalues = cellvalues;
    double nodataval = ds->getNoData(1);
    OGRCoordinateTransformation *poCT = NULL;//
    OGRCoordinateTransformation *poCT_Inverse = NULL;//
    if (!oTargetSRS.IsSame(&oSourceSRS))
    {
      poCT = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
      poCT_Inverse = OGRCreateCoordinateTransformation(&oTargetSRS, &oSourceSRS);
    }

    int fid = -1;
    grid.reset();
    OGREnvelope bound = ds->bound;
    OGREnvelope gridbound = grid.bound;
    for (size_t i = 0; i < grid.slice; i++)
    {
      grid.cells[i] = _NODATA;
    }
    for (size_t row = 0; row < ds->nrows; row++)
    {
      OGREnvelope tmpBound;
      tmpBound.MaxY = bound.MaxY - gridcellsize * row;
      tmpBound.MinY = bound.MaxY - gridcellsize * (row + 1);
      for (size_t col = 0; col < ds->ncols; col++)
      {
        OGREnvelope cellbound = tmpBound;
        cellbound.MinX = bound.MinX + gridcellsize * col;
        cellbound.MaxX = bound.MinX + gridcellsize * (col + 1);
        OGREnvelope transformedbb = cellbound;
        fid++;
        double val = *pcellvalues++;
        if (val == nodataval)
          continue;
        if (poCT)
        {
          poCT->Transform(1, &transformedbb.MinX, &transformedbb.MaxY);
          poCT->Transform(1, &transformedbb.MaxX, &transformedbb.MinY);
        }
        //printf("%f,%f,%f,%f\n", transformedbb.MinX, &transformedbb.MaxY, transformedbb.MaxX, transformedbb.MinY);
        if (!gridbound.Intersects(transformedbb))
          continue;
        double srcResol = transformedbb.MaxX - transformedbb.MinX;
        double cellscale = srcResol / cellsize * 25;
        //srcResol = gridcellsize / cellscale;
        srcResol = srcResol / 20;
        //std::map<int, FeatureFootprintCount> feamap;
        gatherFromGridCell(transformedbb, srcResol, &mask, &grid, val, NULL);
        //std::map<int, FeatureFootprintCount>::iterator feamapIter = feamap.begin();
        //while (feamapIter != feamap.end())
        //{
        //	FeatureFootprintCount& feacount = feamapIter->second;
        //	int cellval = (int)grid.cells[feacount.cellid];
        //	if (!usedMask || cellval != mask.nodata)
        //	{
        //		if (grid.cells[feacount.cellid] == _NODATA)
        //			grid.cells[feacount.cellid] = 0;
        //		grid.cells[feacount.cellid] += cellval * feacount.frac;
        //	}
        //	feamapIter++;
        //}
      }
    }

    //for (size_t i = 0; i < grid.slice; i++)
    //{
    //	if (grid.cells[i] != _NODATA)
    //		grid.cells[i] = log(grid.cells[i]);
    //}
    grid.toRaster(outputfile, grid.proj, _NODATA);
  }

  static void intersectPoints(std::string fishnetfile, std::string inputfile, std::string outputfile, bool usedMask = true)
  {
    Grid grd;
    grd.fromFishnetRaster(fishnetfile, true);
    Grid* grid = &grd;
    ShapeFile output;
    output.create(outputfile, NULL, NULL, OGRwkbGeometryType::wkbPoint);
    int feaField = output.getOrCreateField("feaId", OGRFieldType::OFTInteger);
    int gridField = output.getOrCreateField("gridId", OGRFieldType::OFTInteger);
    int fracField = output.getOrCreateField("fraction", OGRFieldType::OFTReal);

    double cellsize = grid->_adfGeoTransform[1];
    OGREnvelope bound = grid->bound;


    ShapeFile input(inputfile.data());
    OGRFeature* poFeature;
    int fid = -1;
    TileRenderer renderer;
    int caField = input.getOrCreateField("ca11", OGRFieldType::OFTReal);
    while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
    {
      fid++;
      OGRGeometry* poGeometry = poFeature->GetGeometryRef();
      if (!poGeometry)
      {
        OGRFeature::DestroyFeature(poFeature);
        continue;
      }
      OGRwkbGeometryType gtype = poGeometry->getGeometryType();
      if (gtype == OGRwkbGeometryType::wkbPoint || gtype == OGRwkbGeometryType::wkbPoint25D)
      {
        OGRPoint* pt = (OGRPoint*)poGeometry;
        if (pt->getX() >= bound.MinX && pt->getX() <= bound.MaxX && pt->getY() >= bound.MinY && pt->getY() <= bound.MaxY)
        {
          unsigned int row = (unsigned int)((bound.MaxY - pt->getY()) / cellsize);
          unsigned int col = (unsigned int)((pt->getX() - bound.MinX) / cellsize);
          int cellidx = col + row * grid->ncols;
          int cellval = (int)grid->cells[cellidx];
          if (cellidx < grid->slice)
          {
            if (!usedMask || cellval != grid->nodata)
            {
              //mask->cells[cellidx] += value;
              OGRFeature* poFeaFrac = OGRFeature::CreateFeature(output.poLayer->GetLayerDefn());
              poFeaFrac->SetGeometry(pt);
              poFeaFrac->SetField(feaField, fid);
              poFeaFrac->SetField(gridField, cellidx);
              poFeaFrac->SetField(fracField, poFeature->GetFieldAsDouble(caField));
              output.poLayer->CreateFeature(poFeaFrac);
              OGRFeature::DestroyFeature(poFeaFrac);
            }
          }
        }

      }
      OGRFeature::DestroyFeature(poFeature);
      if (fid % 100 == 0)
        printf("%d\n", fid);
    }

  }

  static void intersectShape(Grid* grid, std::string inputfile, std::string outputfile, bool usedMask = true)
  {
    /*ShapeFile output;
    output.create(outputfile, NULL, NULL, OGRwkbGeometryType::wkbNone);
    int feaField = output.getOrCreateField("feaId", OGRFieldType::OFTInteger);
    int gridField = output.getOrCreateField("gridId", OGRFieldType::OFTInteger);
    int fracField = output.getOrCreateField("fraction", OGRFieldType::OFTReal);*/
    SparseFractionGrid sparseGrid;

    /*char wkt[512];
    char* pwkt = wkt;
    if (input.poLayer->GetSpatialRef())
    input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);*/
    double cellsize = grid->_adfGeoTransform[1];
    OGREnvelope bound = grid->bound;

    //}
    //else
    //{
    ShapeFile input(inputfile.data());
    OGRFeature* poFeature;
    int fid = -1;
    TileRenderer renderer;
    RASTERIZATION_RESOLUTION = cellsize / 40;
    double LINE_SAMPLE_DIST = cellsize / 40;
    while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
    {
      fid++;

      OGRGeometry* poGeometry = poFeature->GetGeometryRef();
      if (!poGeometry)
      {
        OGRFeature::DestroyFeature(poFeature);
        continue;
      }
      OGRwkbGeometryType gtype = poGeometry->getGeometryType();
      if (gtype == OGRwkbGeometryType::wkbPoint || gtype == OGRwkbGeometryType::wkbPoint25D)
      {
        OGRPoint* pt = (OGRPoint*)poGeometry;
        if (pt->getX() >= bound.MinX && pt->getX() <= bound.MaxX && pt->getY() >= bound.MinY && pt->getY() <= bound.MaxY)
        {
          unsigned int row = (unsigned int)((bound.MaxY - pt->getY()) / cellsize);
          unsigned int col = (unsigned int)((pt->getX() - bound.MinX) / cellsize);
          int cellidx = col + row * grid->ncols;
          int cellval = (int)grid->cells[cellidx];
          if (cellidx < grid->slice)
          {
            if (!usedMask || cellval != grid->nodata)
            {
              //mask->cells[cellidx] += value;
              //OGRFeature* poFeaFrac = OGRFeature::CreateFeature(output.poLayer->GetLayerDefn());
              //poFeaFrac->SetGeometry(pt);
              //poFeaFrac->SetField(feaField, fid);
              //poFeaFrac->SetField(gridField, cellidx);
              //poFeaFrac->SetField(fracField, 1);
              //output.poLayer->CreateFeature(poFeaFrac);
              SparseGridCell spcell(fid, cellidx, 1.0);
              sparseGrid.cells.push_back(spcell);
              //OGRFeature::DestroyFeature(poFeaFrac);
            }
          }
        }

      }
      else
      {
        bool ispolygon = false;
        if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon || gtype == wkbMultiPolygon25D)
          ispolygon = true;
        std::vector<OGRGeometry*> geomlist;
        if (dynamic_cast<OGRGeometryCollection*>(poGeometry))
        {
          OGRGeometryCollection* geomcollect = dynamic_cast<OGRGeometryCollection*>(poGeometry);
          for (int igeom = 0; igeom < geomcollect->getNumGeometries(); igeom++)
          {
            geomlist.push_back(geomcollect->getGeometryRef(igeom));
          }
        }
        else
        {
          geomlist.push_back(poGeometry);
        }
        std::map<int, FeatureFootprintCount> feamap;
        double totalarea = 0;
        for (int igeom = 0; igeom < geomlist.size(); igeom++)
        {
          OGRGeometry* poGeometry = geomlist[igeom];
          OGREnvelope geobb;
          poGeometry->getEnvelope(&geobb);
          if (!bound.Intersects(geobb))
            continue;
          if (ispolygon)
          {
            //double adjustedResol = RASTERIZATION_RESOLUTION;
            //double bbarea = (geobb.MaxX - geobb.MaxX) * (geobb.MaxY- geobb.MaxY);
            //if (bbarea / adjustedResol > 100000 * 100000)
            //{
            //	adjustedResol = adjustedResol / (100000 * 100000);
            //}
            renderer.init(RASTERIZATION_RESOLUTION, geobb);
            for (int itile = 0; itile < renderer.numtiles; itile++)
            {
              RasterMask* mask = renderer.drawTile(itile, poGeometry, ispolygon);
              totalarea += gatherFromPolygonMap(mask, grid, fid, feamap);

              //if (gtype == wkbMultiPolygon || gtype == wkbMultiPolygon25D)
              //{
              //	uchar* data = mask->_image->bits();
              //	for (int ipixel = 0; ipixel < mask->_xsize * mask->_ysize; ipixel++)
              //	{
              //		if (data[1] == 255)
              //		{
              //			data[1] = rand() % 255;
              //			data[2] = rand() % 255;
              //			data[3] = rand() % 255;
              //		}
              //		data += 4;
              //	}
              //	std::stringstream ssSrc;
              //	ssSrc << fid << "_" << igeom << "_" << itile << "_.tif";
              //	std::stringstream destSrc;
              //	destSrc << fid << "_" << igeom << "_" << itile << ".tif";
              //	mask->save(ssSrc.str());
              //	const char *pszFormat = "GTiff";
              //	GDALDriver *poDriver;
              //	char **papszMetadata;
              //	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
              //	GDALDataset *poSrcDS =
              //		(GDALDataset *)GDALOpen(ssSrc.str().data(), GA_ReadOnly);
              //	GDALDataset *poDstDS;
              //	poDstDS = poDriver->CreateCopy(destSrc.str().data(), poSrcDS, FALSE,
              //		NULL, NULL, NULL);
              //	/* Once we're done, close properly the dataset */
              //	poDstDS->SetGeoTransform(mask->_adfGeoTransform);
              //	poDstDS->SetProjection(mask->proj.data());

              //	if (poDstDS != NULL)
              //		GDALClose((GDALDatasetH)poDstDS);
              //	GDALClose((GDALDatasetH)poSrcDS);
              //	poDriver->Delete(ssSrc.str().data());
              //}

              delete mask;
            }
          }
          else
          {
            if (dynamic_cast<OGRLineString*>(poGeometry))
            {
              OGRLineString* linestr = (OGRLineString *)poGeometry;
              totalarea += gatherFromPolyline(linestr, LINE_SAMPLE_DIST, grid, fid, feamap);
            }
            else if (dynamic_cast<OGRCircularString*>(poGeometry))
            {
              OGRCircularString* circularstr = (OGRCircularString *)poGeometry;
              OGRLineString* linestr = dynamic_cast<OGRLineString*>(circularstr->getLinearGeometry());
              totalarea += gatherFromPolyline(linestr, LINE_SAMPLE_DIST, grid, fid, feamap);
            }

          }


        }
        std::map<int, FeatureFootprintCount>::iterator feamapIter = feamap.begin();
        /*

        while (feamapIter != feamap.end())
        {
        totalarea += feamapIter->second.frac;
        feamapIter++;
        }
        */
        //feamapIter = feamap.begin();
        double mytotal = 0;
        feamapIter = feamap.begin();
        while (feamapIter != feamap.end())
        {
          FeatureFootprintCount& feacount = feamapIter->second;
          mytotal = mytotal + feacount.frac;
          feamapIter++;
        }


        feamapIter = feamap.begin();
        while (feamapIter != feamap.end())
        {
          FeatureFootprintCount& feacount = feamapIter->second;
          feacount.totalarea = totalarea;
          double frac = feacount.frac / totalarea;
          //mask->cells[feacount.cellid] += (frac * value);
          //intersectedFeaturesCount.push_back(feacount);

          int cellval = (int)grid->cells[feacount.cellid];
          if (!usedMask || cellval != grid->nodata)
          {
            //OGRFeature* poFeaFrac = OGRFeature::CreateFeature(output.poLayer->GetLayerDefn());
            //poFeaFrac->SetField(feaField, fid);
            //poFeaFrac->SetField(gridField, feacount.cellid);
            //poFeaFrac->SetField(fracField, frac);
            //output.poLayer->CreateFeature(poFeaFrac);
            //OGRFeature::DestroyFeature(poFeaFrac);

            SparseGridCell spcell(fid, feacount.cellid, frac);
            sparseGrid.cells.push_back(spcell);
          }

          feamapIter++;
        }
      }

      OGRFeature::DestroyFeature(poFeature);
      if (fid % 100 == 0)
        printf("%d\n", fid);
    }
    printf("number of features = %d\n", fid);
    //}

    //sparseGrid.save(outputfile);
    sparseGrid.saveByDatatype(outputfile);
    //sparseGrid.toDBF(outputfile.substr(0, outputfile.size() - 3) + "dbf");
  }

  static void intersectFolder(std::string indir, std::string fishnetfile, std::string outdir, std::map<std::string, LayerInfo> spatialFileMap)
  {
    printf("%s\n", indir.data());
    QDir qoutdir(outdir.data());
    if (!qoutdir.exists())
      qoutdir.mkpath(".");
    outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();

    std::vector<std::string> files;
    QDir input_dir(indir.data());
    input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
    input_dir.setSorting(QDir::Name);
    indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
    QFileInfoList list = input_dir.entryInfoList();
    Grid grid;
    grid.fromFishnetRaster(fishnetfile, true);
    for (int i = 0; i < list.size(); ++i) {
      QFileInfo fileInfo = list.at(i);
      std::string inputfile = fileInfo.absoluteFilePath().toLocal8Bit().data();
      std::string outputfile = outdir + fileInfo.completeBaseName().toLocal8Bit().data() + ".bin";

      std::string name = fileInfo.completeBaseName().toLocal8Bit().data();
      if (spatialFileMap.size() > 0 && spatialFileMap.find(name) == spatialFileMap.end())
        continue;
      LayerInfo info = spatialFileMap[name];
      if (!fileInfo.fileName().endsWith(".shp", Qt::CaseInsensitive) && !fileInfo.fileName().endsWith(".tif", Qt::CaseInsensitive))
      {
        continue;
      }

      printf("%s\n", inputfile.data());
      if (fileInfo.fileName().endsWith(".tif", Qt::CaseInsensitive))
      {
        outputfile = outdir + fileInfo.completeBaseName().toLocal8Bit().data() + ".tif";
        if (QFileInfo(outputfile.data()).exists())
        {
          continue;
        }
        intersectRaster(fishnetfile,
          inputfile, outputfile, info.useMask);
        continue;
      }
      ShapeFile shps(fileInfo.absoluteFilePath().toLocal8Bit().data());
      bool isPoint = false;
      OGRwkbGeometryType geotp = shps.poLayer->GetLayerDefn()->GetGeomType();
      if (geotp == wkbPoint || geotp == wkbMultiPoint || geotp == wkbPointM || geotp == wkbMultiPointM || geotp == wkbPointZM || geotp == wkbMultiPointZM || geotp == wkbPoint25D || geotp == wkbMultiPoint25D)
      {
        isPoint = true;
      }
      if (!isPoint) {
        if (!fileInfo.fileName().endsWith(".shp", Qt::CaseInsensitive) || QFileInfo(outputfile.data()).exists()) {
          continue;
        }
      }

      intersectShape(&grid, inputfile, outputfile, info.useMask);
      printf("%s\n", outputfile.data());

    }

  }

  static void intersectFolder(std::string indir, std::string fishnetfile, std::string outdir)
  {
    QDir qoutdir(outdir.data());
    if (!qoutdir.exists())
      qoutdir.mkpath(".");
    outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
    std::vector<std::string> files;
    QDir input_dir(indir.data());
    input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
    input_dir.setSorting(QDir::Name);
    indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
    QFileInfoList list = input_dir.entryInfoList();
    Grid grid;
    grid.fromFishnetRaster(fishnetfile, true);
    for (int i = 0; i < list.size(); ++i) {
      QFileInfo fileInfo = list.at(i);
      std::string inputfile = fileInfo.absoluteFilePath().toLocal8Bit().data();
      std::string outputfile = outdir + fileInfo.completeBaseName().toLocal8Bit().data() + ".bin";

      std::string name = fileInfo.completeBaseName().toLocal8Bit().data();
      if (!fileInfo.fileName().endsWith(".shp", Qt::CaseInsensitive) || QFileInfo(outputfile.data()).exists())
      {
        continue;
      }
      intersectShape(&grid, inputfile, outputfile, false);
      printf("%s\n", outputfile.data());
    }

  }

  static void intersectShape(OGREnvelope bound, double gridcellsize, std::string inputfile, std::string outputfile, bool usedMask = true)
  {
    Grid grid(bound, gridcellsize, 0);
    grid.reset();
    intersectShape(&grid, inputfile, outputfile, usedMask);
  }

  static void intersectShape(std::string fishnet_rasterfile, std::string inputfile, std::string outputfile, bool usedMask = true)
  {

    Grid grid;
    grid.fromFishnetRaster(fishnet_rasterfile);
    grid.reset();
    intersectShape(&grid, inputfile, outputfile, usedMask);
  }

private:
  static double gatherFromPolygonMap(RasterMask* mask, Grid* grid, int fid, std::map<int, FeatureFootprintCount>& featureCountMap)
  {


    size_t stride = mask->_stride;
    double* geoTransform = mask->_adfGeoTransform;
    double maskCellSize = abs(geoTransform[5]);

    size_t num = 0;
    size_t numpixels = mask->_ysize*mask->_xsize;
    //uchar *pData = mask->_image->bits();
    unsigned int* pData = (unsigned int*)mask->_image->bits();
    //int whiteidx = 0;
    /*for (size_t i = 0; i < numpixels; i++)
    {
    if (pData[0] == 255 && pData[3] == 255)
    {
    num++;
    }
    pData += stride;
    }*/
    //4294967295
    //for (size_t i = 0; i < numpixels; i++)
    //{
    //	if (*pData == 4294967295)
    //		num++;
    //	pData++;
    //}

    pData = (unsigned int*)mask->_image->bits();
    //double fraction = (double)1 / (double)num;
    double total = 0;
    double cellsize = grid->_adfGeoTransform[1];
    int maxcellidx = grid->nrows * grid->ncols;
    double total2 = 0;
    for (size_t i = 0; i < mask->_ysize; i++)
    {
      float y = geoTransform[3] + geoTransform[5] * 0.5 + geoTransform[5] * i;
      unsigned int yindex = (unsigned int)((grid->bound.MaxY - y) / cellsize);

      for (size_t j = 0; j < mask->_xsize; j++)
      {
        float x = geoTransform[0] + geoTransform[1] * 0.5 + geoTransform[1] * j;
        unsigned int xindex = (unsigned int)((x - grid->bound.MinX) / cellsize);
        if (*pData == 4294967295)
        {
          total += maskCellSize;
          if (x > grid->bound.MaxX || x < grid->bound.MinX || y > grid->bound.MaxY || y < grid->bound.MinY)
          {
            pData++;
            continue;
          }
          int cellidx = xindex + yindex * grid->ncols;
          total2 += maskCellSize;
          if (cellidx > 0 && cellidx < maxcellidx)
          {

            std::map<int, FeatureFootprintCount>::iterator findIter = featureCountMap.find(cellidx);

            if (findIter != featureCountMap.end())
            {
              findIter->second.frac += maskCellSize;
            }
            else
            {
              FeatureFootprintCount fea;
              fea.cellid = cellidx;
              fea.frac = maskCellSize;
              fea.fid = fid;
              featureCountMap[cellidx] = fea;
            }
            //mask->cells[xindex + yindex * mask->ncols] += (fraction * total);
            //mask->cells[xindex + yindex * mask->ncols] ++;
            //mask->cells[xindex + yindex * mask->ncols] += fraction;

          }

        }
        pData++;
      }
    }
    printf("%f\n", total2);
    return total;
  }

  static double gatherFromPolyline(const OGRLineString  *polyline, double& sampleDist, Grid* grid, int fid, std::map<int, FeatureFootprintCount>& featureCountMap)
  {
    unsigned int maxcellidx = grid->nrows * grid->ncols;
    double cellsize = grid->_adfGeoTransform[1];
    double total = 0;
    //std::vector<OGRPoint> points;
    //for (int i = 0; i < polyline->getNumPoints() - 1; i++)
    //{
    //	OGRPoint pt;
    //	PT->getNextPoint(&pt);
    //}
    for (int i = 0; i < polyline->getNumPoints() - 1; i++)
    {
      OGRPoint pt1;
      OGRPoint pt2;
      polyline->getPoint(i, &pt1);
      polyline->getPoint(i + 1, &pt2);
      //OGRPoint pt1 = points[i];
      //OGRPoint pt2 = points[i + 1];
      double xdif = pt2.getX() - pt1.getX();
      double ydif = pt2.getY() - pt1.getY();
      double dist = sqrt(xdif*xdif + ydif * ydif);
      int ndivisions = ceil(dist / sampleDist);
      double adjustedSampleDist = dist / ndivisions;
      double xstep = xdif / ndivisions;
      double ystep = ydif / ndivisions;

      double p1x = pt1.getX();
      double p1y = pt1.getY();
      for (int n = 0; n < ndivisions + 1; n++)
      {
        double px = p1x + xstep * n;
        double py = p1y + ystep * n;
        total += adjustedSampleDist;

        if (px > grid->bound.MaxX || px < grid->bound.MinX || py > grid->bound.MaxY || py < grid->bound.MinY)
          continue;
        unsigned int xindex = (unsigned int)((px - grid->bound.MinX) / cellsize);
        unsigned int yindex = (unsigned int)((grid->bound.MaxY - py) / cellsize);
        int cellidx = xindex + yindex * grid->ncols;


        if (cellidx > 0 && cellidx < maxcellidx)
        {
          std::map<int, FeatureFootprintCount>::iterator findIter = featureCountMap.find(cellidx);

          if (findIter != featureCountMap.end())
          {
            findIter->second.frac += adjustedSampleDist;
          }
          else
          {
            FeatureFootprintCount fea;
            fea.cellid = cellidx;
            fea.frac = adjustedSampleDist;
            fea.fid = fid;
            featureCountMap[cellidx] = fea;
          }

        }
      }

    }
    return total;
  }

  static double gatherFromGridCell(OGREnvelope bound, double cellsize, Grid* mask, Grid* grid, double val, OGRCoordinateTransformation *poCT = NULL)
  {

    int rows = (int)(ceil((bound.MaxY - bound.MinY) / cellsize));
    int cols = (int)(ceil((bound.MaxX - bound.MinX) / cellsize));
    int num = rows * cols;
    double fraction = (double)1 / (double)num;
    int maxcellidx = mask->nrows * mask->ncols;
    double targetcellsize = mask->_adfGeoTransform[1];
    for (size_t i = 0; i < rows; i++)
    {
      double py = bound.MaxY - cellsize * 0.5 - cellsize * i;


      for (size_t j = 0; j < cols; j++)
      {
        double px = bound.MinX + cellsize * 0.5 + cellsize * j;
        double x = px;
        double y = py;
        if (poCT)
          poCT->Transform(1, &x, &y);
        if (x > mask->bound.MaxX || x < mask->bound.MinX || y > mask->bound.MaxY || y < mask->bound.MinY)
          continue;
        unsigned int yindex = (unsigned int)((mask->bound.MaxY - y) / targetcellsize);
        unsigned int xindex = (unsigned int)((x - mask->bound.MinX) / targetcellsize);
        int cellidx = xindex + yindex * mask->ncols;
        if (cellidx > 0 && cellidx < maxcellidx)
        {
          if (mask->cells[cellidx] == mask->nodata)
            continue;
          //std::map<int, FeatureFootprintCount>::iterator findIter = featureCountMap.find(cellidx);
          if (grid->cells[cellidx] == _NODATA)
            grid->cells[cellidx] = 0;
          grid->cells[cellidx] += (fraction * val);
          //if (findIter != featureCountMap.end())
          //{
          //	findIter->second.frac += fraction;
          //}
          //else
          //{
          //	FeatureFootprintCount fea;
          //	fea.cellid = cellidx;
          //	fea.frac = fraction;
          //	fea.fid = fid;
          //	featureCountMap[cellidx] = fea;
          //}

        }
      }
    }
    return 1;
  }

};

