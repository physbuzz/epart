
#include <cstdarg> 
#include <string>
#include <vector>


CSVData(std::string filename) { }
void newColumn(std::string name) { }
void newRow(std::vector<float> arg) { }
void newRow(std::vector<std::string> arg) { }
void newRow(float arg, ...) { }
void save() { }
