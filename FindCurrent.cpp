#include <iostream>
#include <string>
#include <ausa/util/FileUtil.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Default.h>
#include <Math/Vector3D.h>
#include <TROOT.h>
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include <fstream>
#include "TFitResult.h"
#include "include/runner.h"
#include <regex>
#include <array>
#include <iterator>
#include <TGraph.h>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace ROOT;

double linstep(double *x, double *par){
    return x[0]*par[0] + par[1];
}

vector<double> findCurrent(string in){
    string analyzed = "Al"+regex_replace(in, regex(R"([\D])"), "") + ".root";
    TFile *myFile = TFile::Open(analyzed.c_str());
    TTree *t = (TTree*)myFile->Get("h101");

    UInt_t vcharge;
    UInt_t clock;
    t ->SetBranchAddress("VCHARGE",&vcharge);
    t ->SetBranchAddress("CLOCK",&clock);

    t->GetEntry(0);
    double firstClock = clock;
    double firstCharge = vcharge;

    int entries = t -> GetEntries();

    t->GetEntry(entries - 1);

    UInt_t stepCharge;
    UInt_t lastClock = 0;
    UInt_t lastCharge;
    UInt_t stepClock;
    stepClock += 38227290460;

    auto gr = new TGraph();

    for(int i = 0; i < entries; i++){
        t->GetEntry(i);
        if(vcharge > firstCharge + 10){
            if(lastClock < clock){
                lastClock = clock;
                lastCharge = vcharge;
            }
            if(stepClock > clock){
                stepClock = clock;
                stepCharge = vcharge;
            }
        }
    }

    double deltaCharge = lastCharge - stepCharge;
    double deltaClock =  lastClock - stepClock;
    double current = deltaCharge/(deltaClock/1000);

    return {deltaClock, deltaCharge};
}