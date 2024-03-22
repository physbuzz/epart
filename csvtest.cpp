
#include <iostream>
#include <vector>
#include "CSVData.h"

using namespace std;

int main(){
    CSVData csv("csvtest.csv");
    csv.newColumn("n");
    csv.newColumn("f(n)");
    for(int n=0;n<10;n++){
        vector<float> arg({float(n),float(n*n)});
        csv.newRow(arg);
    }
    csv.save();
    return 0;

}
