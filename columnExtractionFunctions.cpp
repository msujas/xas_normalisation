#include <iostream>
using namespace std;
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <tuple>
#include <filesystem>

float angle_to_kev(float angle){
    float pi = 3.14159265;
    float dspacing = 3.133789;
    float planck = 6.62607015e-34;
    float charge = 1.60217663e-19;
    int speedOfLight = 299792458;
    float wavelength = 2*dspacing*sin(angle*pi/(180));
    float wavelength_m = wavelength*pow(10,-10);
    float energy_kev = planck*speedOfLight/(wavelength_m*charge*1000);
    return energy_kev;
}

string numberFormat(int n, int width){
    string outstring(width,'0');
    string strDigit = to_string(n);
    for (int i = strDigit.length(); i--; i==0){
        outstring[width+i-strDigit.length()] = strDigit[i];
    }
    return outstring;
}

bool isIn(string s, string pattern){
    if (s.find(pattern)!= -1){
        return true;
    }
    else {
        return false;
    }
}



vector<string> splitString(string inputString, string delimiter){
    int pos;
    string substring;
    vector<string> outVector;
    string  tempString = inputString;
    if (delimiter ==" "){
        while (isIn(tempString,"  ")){ // check if 2 consecutive spaces in string, replace with 1
            tempString = tempString.replace(tempString.find("  "),2," ");

        }
        if (tempString.find(" ") == 0){
            tempString = tempString.replace(0,1,"");
        }
    }
    while (true){
        pos = tempString.find(delimiter);
        substring = tempString.substr(0,pos);
        outVector.push_back(substring);
        if (tempString.find(delimiter) == -1){
            break;
        }
        tempString.erase(0,pos+delimiter.length());    
    }
    return outVector;
}

vector<float> splitStringFloat(string inputString, string delimiter){
    int pos;
    string substring;
    vector<float> outVector;
    string tempString = inputString;
    if (delimiter ==" "){
        while (isIn(tempString,"  ")){
            tempString = tempString.replace(tempString.find("  "),2," ");
        }
        if (tempString.find(" ") == 0){
            tempString = tempString.replace(0,1,"");
        }
    }
    while (true){
        pos = tempString.find(delimiter);
        substring = tempString.substr(0,pos);
        outVector.push_back(stof(substring));
        if (tempString.find(delimiter) == -1){
            break;
        }
        tempString.erase(0,pos+delimiter.length());    
    }
    return outVector;
}

vector<vector<float>> transposeVector(vector<vector<float>> inputVector){

    vector<vector<float>> transposedVector;
    for (int i=0;i < inputVector[0].size();i++){
        transposedVector.push_back({});
        for (int j = 0; j < inputVector.size();j++){
            transposedVector[i].push_back(inputVector[j][i]);
        }
    }
    return transposedVector;
}

void print1dStringVector(vector<string> inputVector){
    cout << "{";
    for (int i = 0; i < inputVector.size(); i++){
        cout << inputVector[i];
        if (i != inputVector.size()-1){
            cout << ',';}
    }
    cout << "}" << endl;
}

void print2dfloatVector(vector<vector<float>> vec){

    cout << "{";
    for (int i = 0; i < vec.size(); i++){
        cout << "{";
        for (int j=0; j < vec[0].size(); j++){
        cout << vec[i][j];
        if (j != vec[0].size()-1){
            cout << ',';}
            }
        cout << "}";
        if (i != vec.size()-1){
            cout << ','<< endl;
        }
    }
    cout << "}" << endl;
}

int getLastScan (string file){
    ifstream input;
    vector<string> scanlines;
    input.open(file);
    if (input.is_open()){
        string s;
        while (getline(input,s)){
        if (isIn(s,"#S ")){
            scanlines.push_back(s);
        }
        }
    }

    int lastscan = stoi(splitString(scanlines.back(), " ")[1]);
    return lastscan;
}

void datToVector (string filename, vector<string> columnNames){

    ifstream input;
    string file2 = filename;
    string basefile = file2.replace(file2.find(".dat"),4,"");
    filesystem::create_directory(basefile);
    input.open(filename);

    string s;
   
    int startLine;
    int endLine;
    bool onscan = false;
    int spectrumCount = -1;

    vector<string> fileLineSplit;
    vector<string> arrayHeader;
    while (getline(input,s)){
        fileLineSplit.push_back(s);
    }
    input.close();
    int line = 0;
    string scanLine;
    string dtLine;
    string headString;
    vector<vector<float>> dataArray;
    bool scanStart = false;
    string columnString = "";
    for (int i = 0; i< columnNames.size(); i++){
        columnString += columnNames[i] + " ";
    }
    for (int i = 0; i < fileLineSplit.size(); i++){

        s = fileLineSplit[i];
        if (isIn(s,"#S") && isIn(s,"zapline")){
            dataArray = {};
            scanLine = s;
            onscan = true;
            spectrumCount++;
            }
        else if (isIn(s,"#D") && onscan){
            dtLine = s;
        }
        else if (isIn(s,"#L") && onscan){
            startLine = line+1;
            headString = s.replace(s.find("#L "),3,"");
            arrayHeader = splitString(headString," ");
            scanStart = true;
        }
        else if (scanStart && !isIn(s,"#")){
            dataArray.push_back(splitStringFloat(s," "));
        }
        else if (isIn(s,"#C") && onscan){
            onscan = false;
            scanStart = false;
            endLine = line -1;
            vector<vector<float>> transposedArray = transposeVector(dataArray);
            vector<vector<float>> columnSelectArray;



            for (int j = 0; j < columnNames.size(); j++){
                for (int k=0; k < arrayHeader.size(); k++) {
                    if (isIn(columnNames[j], arrayHeader[k])){
                        columnSelectArray.push_back(transposedArray[k]);
                    }
                }
            }
                vector<vector<float>> tArray = transposeVector(columnSelectArray);
                ofstream outfile;
                string number = numberFormat(spectrumCount,4);
                
                string newfilename = basefile + "/"+ basefile + "_" + number + ".dat";
                cout << newfilename << "\n";
                outfile.open(newfilename);

                outfile << scanLine << "\n" << dtLine<< "\n" << columnString << "\n";
                
                for (int j=0; j< tArray.size();j++ ){
                    
                    for (int k=0 ; k< tArray[0].size(); k++){
                        outfile << tArray[j][k] << " ";
                    }
                    outfile << "\n";
                }
                outfile.close();

        }
        /*
        if (onscan == true){
        cout << s << endl;
        }
        */
        line++;
    }
    
    
    //tuple<vector<vector<vector<float>>>,vector<string>, vector<string>, string> outTuple = {allArrays, scanLines, dtLines, columnString};
    
    
}



