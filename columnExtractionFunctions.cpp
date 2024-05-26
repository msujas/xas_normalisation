#include <iostream>
using namespace std;
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <tuple>

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

tuple<vector<vector<float>>,string, string, string> datToVector (string filename, vector<string> columnNames, int scanno){

    ifstream input;
    input.open(filename);

    string s;
   
    int startLine;
    int endLine;
    bool onscan = false;
    int spectrumCount = 0;
    vector<string> fileLineSplit;
    vector<string> arrayHeader;
    while (getline(input,s)){
        fileLineSplit.push_back(s);
    }
    int line = 0;
    vector<string> lineVector;
    string scanLine;
    string dtLine;
    string headString;
    for (int i = 0; i < fileLineSplit.size(); i++){
        s = fileLineSplit[i];
        string scanPattern = "#S " + to_string(scanno);
        if (isIn(s,scanPattern) && isIn(s,"zapline")){
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
        }
        else if (isIn(s,"#C") && onscan){
            onscan = false;
            endLine = line -1;
            //cout << startLine << " " << endLine << endl;
            for (int j = startLine; j < endLine; j++){
                lineVector.push_back(fileLineSplit[j]);
                //cout << lineVector[j - startLine] << endl;
            }
        break;
        }
        /*
        if (onscan == true){
        cout << s << endl;
        }
        */
        line++;
    }
    input.close();
    
    vector<vector<float>> dataArray;
    for (int i = 0; i < lineVector.size();i++){

        dataArray.push_back(splitStringFloat(lineVector[i]," "));
    }

    vector<vector<float>> transposedArray = transposeVector(dataArray);
    vector<vector<float>> columnSelectArray;
    string columnString = "";
    for (int i = 0; i< columnNames.size(); i++){
        columnString += columnNames[i] + " ";
    }

        //cout << arrayHeader[i] << endl;
    for (int j = 0; j < columnNames.size(); j++){
        for (int i=0; i < arrayHeader.size(); i++) {
            if (isIn(columnNames[j], arrayHeader[i])){
                columnSelectArray.push_back(transposedArray[i]);
            }
        }
    }

    tuple<vector<vector<float>>,string, string, string> outTuple = {columnSelectArray, scanLine, dtLine, columnString};
    return  outTuple;
    
}



