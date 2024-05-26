#include <iostream>
using namespace std;
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <tuple>
#include "columnExtractionFunctions.cpp"




int main(int argc, char *argv[2]){
float thetaOffset = 0;

string file = "CuFoil_Cu_exafs.dat";
if (argc == 2){
string file = argv[1];
}


string fluoCounter = "xmap_roi00";
string monPattern = "mon_";
string ion1Pattern = "ion_1";

vector<string> counterNames = {"ZapEnergy","TwoTheta","mon_3","mon_4","mon_1","ion_1_2","ion_1_3","ion_1_1",fluoCounter};
vector<string> counterNames_NF; //NF - no fluorescence

for (int i = 0; i < counterNames.size(); i++){
    string item = counterNames[i];
    if (item != fluoCounter){
        counterNames_NF.push_back(item);
    }
}
//tuple<vector<vector<float>> , vector<string>> arrayandHeaders = datToVector(file, counterNames);
vector<vector<vector<float>>> allArrays;
int noscans = getLastScan(file);
vector<string> scanLines;
vector<string> dtLines;
string headerstring;
for (int i=1; i < noscans+1; i++){
auto  [dataArray, scanLine, dtLine, hstring] = datToVector(file, counterNames, i);
headerstring = hstring;
allArrays.push_back(dataArray);
scanLines.push_back(scanLine);
dtLines.push_back(dtLine);
}
string basefile = file.replace(file.find(".dat"),4,"");
for (int i=0;i<allArrays.size();i++){
    vector<vector<float>> tArray = transposeVector(allArrays[i]);
    ofstream outfile;
    string number = numberFormat(i,4);
    
    string filename = basefile + "_" + number + ".dat";
    cout << filename << "\n";
    outfile.open(filename);
    outfile << scanLines[i] << "\n" << dtLines[i]<< "\n" << headerstring << "\n";
    
    //print2dfloatVector(tArray);
    for (int j=0; j< tArray.size();j++ ){
        
        for (int k=0 ; k< tArray[0].size(); k++){
            outfile << tArray[j][k] << " ";
        }
        outfile << "\n";
    }
    outfile.close();
}
cout << allArrays[0].size() << "\n";
cout << transposeVector(allArrays[0]).size() << "\n";
}