#pragma once

#include "framework.h"

class classUtils
{
public:
	bool is_number(const std::string& s);
	void writeLine(string);
	void emptyLine();
	std::vector<std::vector<string> > transpose(const std::vector<std::vector<string> > data);
	std::vector<std::vector<std::string> > csvRead(string filePath);
	double Cardano_min(double a, double b, double c, double d);
};
