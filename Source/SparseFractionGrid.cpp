#include "SparseFractionGrid.h"
#include "ShapeFile.h"
#include "qdir.h"
#include "qfileinfo.h"
SparseFractionGrid::SparseFractionGrid()
{

}


SparseFractionGrid::~SparseFractionGrid()
{

}

void SparseFractionGrid::save(std::string filename)
{
	std::ofstream ofs;
	ofs.open(filename, std::ofstream::binary);
	int num = cells.size();
	ofs.write((char*)(&num), sizeof(int));
	size_t bufsize = sizeof(SparseGridCell) * cells.size();
	ofs.write((char*)(&cells[0]), bufsize);
	ofs.close();
}

void SparseFractionGrid::saveByDatatype(std::string filename)
{
	std::ofstream ofs;
	ofs.open(filename, std::ofstream::binary);
	std::vector<int> feaIds;
	std::vector<int> gridIds;
	std::vector<double> fractions;
	SparseGridCell* pcell = &cells[0];
	for (size_t i = 0; i < cells.size(); i++)
	{
		feaIds.push_back(pcell->feaId);
		gridIds.push_back(pcell->gridId);
		fractions.push_back(pcell->fraction);
		pcell++;
	}
	int num = cells.size();
	ofs.write((char*)(&num), sizeof(int));
	ofs.write((char*)(&feaIds[0]), sizeof(int)*num);
	ofs.write((char*)(&gridIds[0]), sizeof(int)*num);
	ofs.write((char*)(&fractions[0]), sizeof(double)*num);
	ofs.close();
}

void SparseFractionGrid::toDBF(std::string filename)
{
	ShapeFile output;
	output.create(filename, NULL, NULL, OGRwkbGeometryType::wkbNone);
	int feaField = output.getOrCreateField("feaId", OGRFieldType::OFTInteger);
	int gridField = output.getOrCreateField("gridId", OGRFieldType::OFTInteger);
	int fracField = output.getOrCreateField("fraction", OGRFieldType::OFTReal);
	for (size_t i = 0; i < cells.size(); i++)
	{
		OGRFeature* poFeaFrac = OGRFeature::CreateFeature(output.poLayer->GetLayerDefn());
		poFeaFrac->SetField(feaField, cells[i].feaId);
		poFeaFrac->SetField(gridField, cells[i].gridId);
		poFeaFrac->SetField(fracField, cells[i].fraction);
		output.poLayer->CreateFeature(poFeaFrac);
		OGRFeature::DestroyFeature(poFeaFrac);
	}
	output.close();
}
void SparseFractionGrid::open(std::string filename)
{
	std::ifstream ifs;
	ifs.open(filename, std::ofstream::binary);
	int num = 0;
	ifs.read((char*)(&num), sizeof(int));
	cells.resize(num);
	size_t bufsize = sizeof(SparseGridCell) * cells.size();
	ifs.read((char*)(&cells[0]), bufsize);
	ifs.close();
}

void SparseFractionGrid::fromDBF(std::string filename)
{
	cells.clear();
	ShapeFile shp(filename);
	int gridField = shp.poLayer->GetLayerDefn()->GetFieldIndex("gridId");
	int feaField = shp.poLayer->GetLayerDefn()->GetFieldIndex("feaId");
	int fracField = shp.poLayer->GetLayerDefn()->GetFieldIndex("fraction");
	OGRFeature* poFeature;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		int gridId = poFeature->GetFieldAsInteger(gridField);
		int feaId = poFeature->GetFieldAsInteger(feaField);
		double frac = poFeature->GetFieldAsDouble(fracField);
		SparseGridCell spcell(feaId, gridId, frac);
		cells.push_back(spcell);
		OGRFeature::DestroyFeature(poFeature);
	}
}

void SparseFractionGrid::convertDBF2BIN(std::string dir)
{


	std::vector<std::string> files;
	QDir input_dir(dir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	dir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();

	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string outputfile = dir + fileInfo.completeBaseName().toLocal8Bit().data() + ".bin";
		std::string inputfile = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (inputfile.substr(inputfile.length() - 4,4) == ".dbf" && !QFileInfo(outputfile.data()).exists())
		{
			SparseFractionGrid grid;
			grid.fromDBF(inputfile);
			grid.save(outputfile);
		}
	}
}
