#ifndef CSVDATA_H
#define CSVDATA_H

#include <cstdarg> 
#include <string>
#include <vector>
#include <sstream>
#include <fstream>


class CSVData { 
    std::stringstream data;
    std::string fname;
    bool firstCol;
    bool firstRow;
public: 
    CSVData(std::string filename) : data(), fname(filename), firstCol(true),firstRow(true) { }
    template <typename T> void newColumn(T arg) { 
        if(firstCol) {
            data<<arg;
            firstCol=false;
        } else {
            data<<", "<<arg;
        }
    }
    void newRow() {
        if(!firstCol || !firstRow){
            data<<"\n";
            firstRow=false;
        }
        firstCol=true;
    }
    template <typename T> void newRow(std::vector<T> arg) { 
        if(arg.size()==0)
            return;
        if(!firstRow || !firstCol)
            data<<"\n";
        firstRow=false;
        firstCol=false;
        for(size_t i=0;i<arg.size();i++){
            if(i>0)
                data<<", ";
            data<<arg[i];
        }
    }
    /*
    void newRow(std::vector<std::string> arg) { 
        newRow<std::string>(arg);
    }
    void newRow(std::vector<float> arg) { 
        newRow<float>(arg);
    }*/
    /*
    void newRow(float arg, ...) { 
    }*/
    void save() { 
        std::ofstream output(fname.c_str());
        output<<data.str().c_str();
        output.close();
    }
};
#endif
