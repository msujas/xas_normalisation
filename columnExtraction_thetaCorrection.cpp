#include <iostream>
using namespace std;
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <tuple>
#include "columnExtractionFunctions.cpp"
#include <filesystem>



int main(int argc, char *argv[2]){
float thetaOffset = 0;

string file = "10_ramp_O2_to480C_Cu_xanes.dat";
if (argc == 2){
string file = argv[1];
}
string file2 = file;




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

int noscans = getLastScan(file);
cout << "total scans " << noscans <<"\n";
datToVector(file, counterNames);

}