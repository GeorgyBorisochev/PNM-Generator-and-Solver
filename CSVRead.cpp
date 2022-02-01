#include "framework.h"
#include "CSVRead.h"

//
//CSVRead
//

std::vector<std::vector<std::string> > CSVReader::getData()
{
	std::ifstream file(fileName);
	if (!file.is_open()) {
		cout << " Failed to open" << endl;
	}
	else {
		//cout << "Opened OK" << endl;
	}
	std::vector<std::vector<std::string> > dataList;
	std::string line = "";
	// Iterate through each line and split the content using delimeter
	while (getline(file, line))
	{
		std::vector<std::string> vec;
		boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
		dataList.push_back(vec);
	}
	// Close the File
	file.close();
	return dataList;
}