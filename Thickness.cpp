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
#include <regex>
#include <Math/Vector3D.h>
#include <TROOT.h>
#include <TVector3.h>
#include <include/runner.h>
#include <vector>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace ROOT;

double onAxisSolidAng(double a, double b, double d){
    return 4*TMath::ATan(a*b/TMath::Sqrt(1+a*a+b*b));
}

double offAxisSolidAng(double d, double A, double B, double a, double b){
    double s = 0.25*(onAxisSolidAng(2*(A+a),2*(B+b),d) - onAxisSolidAng(2*A,2*(B+b),d)
             - onAxisSolidAng(2*(A+a),2*B,d) + onAxisSolidAng(2*A,2*B,d));
    return s;
}

//løber igennem listen af pixels og returnerer true hvis pixelen er i listen
tuple<bool,int> findPixel(UInt_t toSearch[1000][3], UInt_t FI, UInt_t BI, UInt_t id, int loopuntil){
    for (int i = 0; i<loopuntil; i++){
        if(toSearch[i][0] == FI && toSearch[i][1] == BI && toSearch[i][2] == id){
            return make_tuple(true,i);
        }
    }
    return make_tuple(false,-1);
}

std::vector<double> thickness(string in){
    //energien denne fil blev optaget ved er givet i dens titel
    int energy = stoi(regex_replace(in, regex(R"([\D])"), ""));
    //skab en pointer til root-filen der er blevet lavet af analyse.
    string analyzed = "analyzed/"+regex_replace(in, regex(R"([\D])"), "") + "a.root";
    TFile *myFile = TFile::Open(analyzed.c_str());
    //Hent træet
    TTree *t = (TTree*)myFile->Get("a");

    //Skab variable som værdier kan loades ned i. Associer dem med branches i træet, så den i. entry bliver lagt
    //deri når vi kalder getEntry(i). Print antallet af entries i træet.
    double_t E[1000];
    double_t scatterAngle[1000];
    Short_t BI[1000];
    Short_t FI[1000];
    Short_t id[1000];
    double_t posx[1000];
    double_t posy[1000];
    double_t posz[1000];

    double_t dirx[1000];
    double_t diry[1000];
    double_t dirz[1000];
    double_t solid[1000];
    UInt_t mul;
    t->SetBranchAddress("pos.fX",posx);
    t->SetBranchAddress("pos.fY",posy);
    t->SetBranchAddress("pos.fZ",posz);
    t->SetBranchAddress("solidAngle",solid);

    t->SetBranchAddress("dir.fX",dirx);
    t->SetBranchAddress("dir.fY",diry);
    t->SetBranchAddress("dir.fZ",dirz);
    t->SetBranchAddress("id",id);
    t->SetBranchAddress("BI",BI);
    t->SetBranchAddress("FI",FI);
    t->SetBranchAddress("E",E);
    t->SetBranchAddress("scatterAngle",scatterAngle);
    t->SetBranchAddress("mul",&mul);
    auto entries = t->GetEntries();

    //skriv peak positioner til en .txt fil

    //Lav et array af histogrammer og vinkler på tilsvarende indekser
    TH1I *histograms[1000] = {};
    double_t solidangles[1000] = {};
    UInt_t pixelInfo[1000][3] = {};
    vector<double_t> positions[1000] = {};
    vector<double_t> directions[1000] = {};

    int lastPrinted = 0;
    for (Int_t i = 0; i < entries; i++) {
        if(i%300000 == 0){
            cout << "Finding 110 degrees in " << energy << "keV data: " << float(i)/float(entries)*100<< "%" << " \r";
            cout.flush();
        }

        //hent entry
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            //hvis vi ikke har set denne vinkel før skal vi lave et nyt histogram for denne vinkel.
            double currentAngle = 0;
            currentAngle += scatterAngle[j];
            //vi gider kun at bruge dem der er scatteret ved ca. 110 grader.
            if(fabs(currentAngle - 110) < 0.5) {
                short currentFI = 0;
                short currentBI = 0;
                short currentid = 0;
                double currentE = 0;
                double currentSolid = 0;
                vector<double_t> currentPos = {};
                currentPos = {posx[j],posy[j],posz[j]};
                vector<double_t> currentDir = {};
                currentDir = {dirx[j],diry[j],dirz[j]};
                currentFI += FI[j];
                currentBI += BI[j];
                currentid += id[j];
                currentE += E[j];
                currentSolid += solid[j];
                auto boolAndIndex = findPixel(pixelInfo, currentFI, currentBI, currentid, lastPrinted);
                //case for hvis der endnu ikke findes et histogram for denne pixel.
                if (!get<0>(boolAndIndex)) {
                    //skab nyt histogram til at indeholde events ved denne vinkel
                    string name = "ID: " + to_string(currentid) + " FI: " + to_string(currentFI) + " BI: " +
                                  to_string(currentBI) + " angle: " + to_string(currentAngle);
                    histograms[lastPrinted] = new TH1I(name.c_str(), name.c_str(), energy, 0.0, energy - 1);
                    pixelInfo[lastPrinted][0] = currentFI;
                    pixelInfo[lastPrinted][1] = currentBI;
                    pixelInfo[lastPrinted][2] = currentid;
                    solidangles[lastPrinted] = currentSolid;
                    positions[lastPrinted] = currentPos;
                    directions[lastPrinted] = currentDir;
                    //fyld energien ind i det skabte histogram
                    histograms[lastPrinted]->Fill(currentE);
                    //læg en til lastPrinted, så vi er klar til næste gang der er en ny vinkel
                    lastPrinted++;
                } else {
                    //der fandtes allerede et histogram for dette pixel
                    histograms[get<1>(boolAndIndex)]->Fill(E[j]);
                }
            }
        }
    }
    string histRoot = "summedhistograms.root";
    TH1I *summedHist = new TH1I("summed", "summed", energy, 0.0, energy - 1);
    TFile output(histRoot.c_str(), "RECREATE");
    double totalSolid = 0;
    double totalSolid2 = 0;
    for(int i = 0; i < lastPrinted; i++){
        TH1I *currentHist = histograms[i];
        //jeg tager kun fra detektor 1
        if(pixelInfo[i][2] == 1 && currentHist -> GetEntries() > 300){
            summedHist -> Add(currentHist);
            totalSolid += solidangles[i];
            //hvis jeg i stedet forsøger at beregne solidangle for denne pixel ved egen kraft:
            double a = 50/16;
            double b = 50/16;
            double d = abs(positions[i][2]);
            double A = abs(positions[i][0]);
            double B = abs(positions[i][1]);
            totalSolid2 += offAxisSolidAng(d, A, B, a, b);
        }
    }
    auto *s = new TSpectrum(100);
    Int_t nfound = s->Search(summedHist, 10, "", 0.05);

    auto xpeaks = s->GetPositionX();

    //For hvert peak skal jeg bruge 3 parametre til at bestemme en gauss. (og så skal jeg give en mere til
    //funktionen for at fortælle den, hvor mange peaks der er, men den fixes til npeaks).
    double par[nfound * 3 + 1];
    //gider ikke dem der har mange peaks
    for (int p = 0; p < nfound; p++) {
        Double_t xp = xpeaks[p];
        Int_t bin = summedHist->GetXaxis()->FindBin(xp);
        Double_t yp = summedHist->GetBinContent(bin);
        par[1 + 3 * p] = yp;
        par[2 + 3 * p] = xp;
        par[3 + 3 * p] = 10;
    }
    TF1 *fit = new TF1("fit", gaussSum, 0, energy, 1 + 3 * nfound);

    fit->SetParameters(par);
    fit->FixParameter(0, nfound);
    fit->SetNpx(2000);
    TFitResultPtr fp = summedHist->Fit("fit","s && Q");
    output.cd();
    summedHist -> Write();
    double area = 0;
    for (int p = 0; p < nfound; p++) {
        if(par[1+3*p] > area){
            area = fp -> Parameter(1+3*p);
        }
    }
    double en = energy;
    return {area,totalSolid,totalSolid2,en};
}