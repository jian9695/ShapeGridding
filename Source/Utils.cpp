#include "Utils.h"
#include <qfileinfo.h>
#include <qdir.h>
#include <fstream>
#include <sstream>
Utils::Utils()
{
}


Utils::~Utils()
{
}

double Utils::calPolygonArea(OGRGeometry* geom)
{

	double frac = 0;
	if (geom->getGeometryType() == wkbPolygon || geom->getGeometryType() == wkbPolygon25D)
	{
		OGRPolygon* poly = (OGRPolygon*)geom;
		frac = poly->get_Area();
	}
	else if (geom->getGeometryType() == wkbMultiPolygon)
	{
		OGRMultiPolygon* multipoly = (OGRMultiPolygon*)geom;
		for (size_t i = 0; i < multipoly->getNumGeometries(); i++)
		{
			OGRPolygon* poly = (OGRPolygon*)multipoly->getGeometryRef(i);
			frac += poly->get_Area();
		}
	}
	return frac;
}


double Utils::calPolylineLength(OGRGeometry* geom)
{
	double len = 0;
	OGRwkbGeometryType gtype = geom->getGeometryType();
	if (geom->getGeometryType() == wkbLineString || geom->getGeometryType() == wkbLineString25D)
	{
		OGRLineString* poly = (OGRLineString*)geom;
		len = poly->get_Length();
	}
	else if (geom->getGeometryType() == wkbMultiLineString)
	{
		OGRMultiLineString* multipoly = (OGRMultiLineString*)geom;
		for (size_t i = 0; i < multipoly->getNumGeometries(); i++)
		{
			OGRLineString* poly = (OGRLineString*)multipoly->getGeometryRef(i);
			len += poly->get_Length();
		}
	}
	return len;
}

void Utils::updateArea(ShapeFile* file, double scale)
{

	OGRFeature *poFeature;


	OGRwkbGeometryType gtype = file->poLayer->GetGeomType();
	if (gtype != wkbPolygon && gtype != wkbPolygon25D && gtype != wkbMultiPolygon)
		return;

	int idx = file->poLayer->GetLayerDefn()->GetFieldIndex("area");
	if (idx < 0)
	{
		OGRFieldDefn def("area", OGRFieldType::OFTReal);
		file->poLayer->CreateField(&def);
		idx = file->poLayer->GetLayerDefn()->GetFieldIndex("area");
	}

	file->poLayer->ResetReading();
	while ((poFeature = file->poLayer->GetNextFeature()) != NULL)
	{
		if (poFeature->GetGeometryRef() == NULL)
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		double frac = calPolygonArea(poFeature->GetGeometryRef()) * scale;
		poFeature->SetField(idx, frac);
		file->poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}


}

void Utils::dbf2csv(std::string shpfile, std::string csvfile)
{
	std::ofstream ofs(csvfile.data());
	ShapeFile shp(shpfile);

	OGRFeature *poFeature;
	shp.poLayer->ResetReading();
	int id = 0;
	OGRFeatureDefn* def = shp.poLayer->GetLayerDefn();
	int numfields = def->GetFieldCount();
	for (size_t i = 0; i < numfields; i++)
	{
		ofs << def->GetFieldDefn(i)->GetNameRef();
		if (i == numfields - 1)
		{
			ofs << std::endl;
		}
		else
		{
			ofs << ",";
		}
	}

	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		for (size_t i = 0; i < numfields; i++)
		{
			ofs << poFeature->GetFieldAsString(i);
			if (i == numfields - 1)
			{
				ofs << std::endl;
			}
			else
			{
				ofs << ",";
			}
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	ofs.close();

}
void Utils::updateLengh(ShapeFile* file, double scale)
{

	OGRwkbGeometryType gtype = file->poLayer->GetGeomType();
	if (gtype != wkbLineString && gtype != wkbMultiLineString && gtype != wkbLineString25D)
		return;

	int idx = file->poLayer->GetLayerDefn()->GetFieldIndex("length");


	if (idx < 0)
	{
		OGRFieldDefn def("length", OGRFieldType::OFTReal);
		OGRErr er = file->poLayer->CreateField(&def);
		idx = file->poLayer->GetLayerDefn()->GetFieldIndex("length");
	}

	OGRFeature *poFeature;
	file->poLayer->ResetReading();
	int id = 0;
	while ((poFeature = file->poLayer->GetNextFeature()) != NULL)
	{

		double length = calPolylineLength(poFeature->GetGeometryRef()) * scale;
		poFeature->SetField(idx, length);
		file->poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}
void Utils::updateFootprint(std::string filename, double scale)
{
	ShapeFile* file = new ShapeFile(filename, 1);
	OGRwkbGeometryType gtype = file->poLayer->GetGeomType();
	int footprintIndex = -1;
	if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D || gtype == wkbMultiLineString25D)
	{
		footprintIndex = file->poLayer->GetLayerDefn()->GetFieldIndex("length");
		updateLengh(file, scale);
	}
	if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype != wkbPolygon25D || gtype == wkbMultiPolygon25D)
	{
		footprintIndex = file->poLayer->GetLayerDefn()->GetFieldIndex("area");
		updateArea(file, scale);
	}
	delete file;


}



/// <summary>
/// </summary>
/// <param name="lat1"></param>
/// <param name="lon1"></param>
/// <param name="lat2"></param>
/// <param name="lon2"></param>
/// <returns></returns>
double Utils::getDistanceFromLatLonInMeter(double lat1, double lon1, double lat2, double lon2)
{
	double R = 6371000; // Radius of the earth in km
	double dLat = (lat2 - lat1) * 0.0174533;  // deg2rad below
	double dLon = (lon2 - lon1) * 0.0174533;
	double a =
		sin(dLat / 2) * sin(dLat / 2) +
		cos(lat1 * 0.0174533) * cos(lat2 * 0.0174533) *
		sin(dLon / 2) * sin(dLon / 2);
	double c = 2 * atan2(sqrt(a), sqrt(1 - a));
	double d = R * c; // Distance in km
	return d;
}
double Utils::getDegreeToMeter(OGRLayer* layer)
{
	OGREnvelope env;
	layer->GetExtent(&env);
	double centerLat = abs((env.MinY + env.MaxY) * 0.5);
	double centerLon = abs((env.MinX + env.MaxX) * 0.5);
	double shift = 0.05;
	double dist1 = getDistanceFromLatLonInMeter(centerLat, centerLon - shift, centerLat, centerLon + shift);
	double dist2 = getDistanceFromLatLonInMeter(centerLat - shift, centerLon, centerLat + shift, centerLon);
	dist1 = (dist1 / (shift * 2));
	dist2 = (dist2 / (shift * 2));

	return (dist1 + dist2) * 0.5;
}


void Utils::updateFootprintForDir(std::string indir, bool force, double scale)
{

	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (isdigit(name[0]))
			continue;
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		//if (!fileInfo.fileName().contains("NonRoad", Qt::CaseInsensitive))
		//	continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
		std::string infile = indir + files[files.size() - 1];
		ShapeFile* file = new ShapeFile(infile, 1);
		OGRwkbGeometryType gtype = file->poLayer->GetGeomType();
		int footprintIndex = -1;
		int footprintIndexNew = -1;
		if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
		{
			footprintIndex = file->poLayer->GetLayerDefn()->GetFieldIndex("length");
			if (footprintIndex < 0)
				Utils::updateLengh(file,scale);
			else if (force)
				Utils::updateLengh(file, scale);
		}
		if (gtype == wkbPolygon || gtype != wkbPolygon25D || gtype == wkbMultiPolygon)
		{
			footprintIndex = file->poLayer->GetLayerDefn()->GetFieldIndex("area");
			if (footprintIndex < 0)
				Utils::updateArea(file, scale);
			else if (force)
				Utils::updateArea(file, scale);
		}
		delete file;
	}
}

std::vector<std::string> Utils::split(const char& delimiter, const std::string& line)
{
	std::string str = line;
	std::vector<std::string> splits;
	if (str.length() < 1)
		return splits;
	std::string::iterator iter = str.begin();
	int index = 0;
	int lastIndex = 0;
	while (iter != str.end() && index < str.length())
	{
		char& c = *iter;
		if (c == delimiter || index == str.length() - 1)
		{
			if (index == str.length() - 1)
				index = str.length();
			if (index - lastIndex>0)
				splits.push_back(str.substr(lastIndex, index - lastIndex));
			lastIndex = index + 1;
		}
		index++;
		iter++;
	}
	return splits;
}

std::vector<std::string> Utils::splitCSV(const char& delimiter, const std::string& line)
{
	std::string str = line;
	std::vector<std::string> splits;
	if (str.length() < 1)
		return splits;
	std::string::iterator iter = str.begin();
	int index = 0;
	int lastIndex = 0;
	bool inQuote = false;
	while (iter != str.end() && index < str.length())
	{
		char& c = *iter;
		if (c == '"' && !inQuote)
		{
			inQuote = true;
		}
		else if (c == '"' && inQuote)
		{
			inQuote = false;
		}


		if (!inQuote && (c == delimiter || index >= str.length() - 1))
		{

			if (index == str.length() - 1)
				index = str.length();
			if (index - lastIndex > 0)
			{

				std::string substr = str.substr(lastIndex, index - lastIndex);
				//if (substr.size() > 1 && substr[0] == ',')
				//	substr = "";
				if (index - lastIndex > 1)
				{
					if (substr[0] == '"' && substr[substr.size() - 1] == '"')
					{
						substr = substr.substr(1, substr.size() - 2);
					}
				}
	/*			printf("%d,%d\n", index, str.length());*/
				if (index >= str.length() - 1 && substr[substr.size() - 1] == ',')
				{
					substr = substr.substr(0, substr.size() - 1);
					splits.push_back(substr);
					splits.push_back("");
				}
				else
				{
					splits.push_back(substr);
				}

			}
			else
			{
				splits.push_back("");
			}

			lastIndex = index + 1;
			inQuote = false;
		}
		index++;
		iter++;
	}
	return splits;
}

std::vector<std::string> Utils::findSubdirectories(std::string indir)
{
	std::vector<std::string> subdirs;
	QDir rootdir(indir.data());
	rootdir.setFilter(QDir::Dirs | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	rootdir.setSorting(QDir::Name);
	std::vector<std::string> files;
	QFileInfoList list = rootdir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string dir = (fileInfo.absoluteFilePath() + "/").toLocal8Bit().data();
		subdirs.push_back(dir);
	}
	return subdirs;
}

std::vector<std::string> Utils::findFiles(std::string indir,std::string match)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (/*!fileInfo.fileName().contains(match.data(),Qt::CaseInsensitive) || */match != "" && !fileInfo.fileName().endsWith(match.data(), Qt::CaseInsensitive))
			continue;
		files.push_back(fileInfo.absoluteFilePath().toLocal8Bit().data());
	}
	return files;

}
std::vector<std::string> Utils::buildVector(std::string dir, std::string* filearr, int numoffiles)
{
	std::vector<std::string> files;
	for (size_t i = 0; i < numoffiles; i++)
	{
		files.push_back(dir + filearr[i]);
	}
	return files;
}

SubdirManager::SubdirManager(std::string dir)
{
	m_dirpath = dir.data();
}

std::vector<std::string> SubdirManager::findFilesMatch(std::string match)
{
	std::vector<std::string> results;
	std::vector<std::string> files = Utils::findFiles(m_dirpath,".shp");
	for (size_t ifile = 0; ifile < files.size(); ifile++)
	{
		QString qfilename = QFileInfo(files[ifile].data()).fileName();
		if (!qfilename.contains(match.data(), Qt::CaseInsensitive) && !qfilename.endsWith(match.data(), Qt::CaseInsensitive))
			continue;
		results.push_back(files[ifile]);
	}
	std::vector<std::string> subdirs = Utils::findSubdirectories(m_dirpath);
	for (size_t idir = 0; idir < subdirs.size(); idir++)
	{
		std::vector<std::string> files = Utils::findFiles(subdirs[idir], ".shp");
		for (size_t ifile = 0; ifile < files.size(); ifile++)
		{
			QString qfilename = QFileInfo(files[ifile].data()).fileName();
			if (!qfilename.contains(match.data(),Qt::CaseInsensitive) && !qfilename.endsWith(match.data(), Qt::CaseInsensitive))
				continue;
			results.push_back(files[ifile]);
		}
	}
	return results;
}

std::vector<std::string> SubdirManager::findFilesMatch(std::vector<std::string> matches)
{
	std::vector<std::string> results;
	for (size_t imatch = 0; imatch < matches.size(); imatch++)
	{
		std::vector<std::string> files = findFilesMatch(matches[imatch]);
		for (size_t ifile = 0; ifile < files.size(); ifile++)
		{
			results.push_back(files[ifile]);
		}
	}
	return results;
}

