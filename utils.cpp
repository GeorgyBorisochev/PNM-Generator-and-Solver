//utils.cpp
//various functions for io etc

#include "utils.h"
#include "framework.h"
#include "CSVRead.h"

//
//UTILS
//

bool classUtils::is_number(const std::string& s)
{
	char* end = nullptr;
	double val = strtod(s.c_str(), &end);
	return end != s.c_str() && *end == '\0' && val != HUGE_VAL;
}

void classUtils::writeLine(string string)
{
	cout << string << "\n";
}

void classUtils::emptyLine()
{
	cout << "\n";
}

std::vector<std::vector<string> > classUtils::transpose(const std::vector<std::vector<string> > data) {
	// this assumes that all inner vectors have the same size and
	// allocates space for the complete result in advance
	std::vector<std::vector<string> > result(data[0].size(),
		std::vector<string>(data.size()));
	for (std::vector<string>::size_type i = 0; i < data[0].size(); i++)
		for (std::vector<string>::size_type j = 0; j < data.size(); j++) {
			result[i][j] = data[j][i];
		}
	return result;
}

double classUtils::Cardano_min(double a, double b, double c, double d)
{
	double Pi = 3.1415926;
	double p = (3 * a * c - pow(b, 2)) / 3 / pow(a, 2);
	double q = 2 * pow(b, 3) / 27 / pow(a, 3) - b * c / 3 / pow(a, 2) + d / a;
	double discrim = pow((q / 2), 2) + pow((p / 3), 3);
	double cardano;

	if (discrim < 0)
	{
		double phi = acos(-q / 2 / sqrt(pow((-p / 3), 3)));
		double y1 = 2 * sqrt(-p / 3) * cos(phi / 3);
		double y2 = 2 * sqrt(-p / 3) * cos(phi / 3 + 2 * Pi / 3);
		double y3 = 2 * sqrt(-p / 3) * cos(phi / 3 + 4 * Pi / 3);

		double x1 = y1 - b / 3 / a;
		cardano = x1;
		double x2 = y2 - b / 3 / a;
		if (x2 < x1)
		{
			cardano = x2;
		}
		double x3 = y3 - b / 3 / a;
		if (x3 < x2)
		{
			cardano = x3;
		}
	}
	else if (discrim > 0)
	{
		double u = pow((-q / 2 + sqrt(discrim)), (1.0 / 3));
		double v = pow((-q / 2 - sqrt(discrim)), (1.0 / 3));

		double y1 = u + v;
		double x1 = y1 - b / 3 / a;
		cardano = x1;
	}
	else
	{
		double y1 = pow(-4 * q, 1.0 / 3);
		double y2 = pow(q / 2, 1.0 / 3);

		double x1 = y1 - b / 3 / a;
		cardano = x1;
		double x2 = y2 - b / 3 / a;
		if (x2 < x1)
		{
			cardano = x2;
		}
	}

	return cardano;
}

std::vector<std::vector<std::string> > classUtils::csvRead(string filePath)
{
	CSVReader reader(filePath);
	// Get the data from CSV File
	std::vector<std::vector<std::string> > dataList = reader.getData();

	return dataList;
}