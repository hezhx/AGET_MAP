// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <fstream>
// #include <vector>

// #include "vector_define.h"
// #include "TrackFit.h"

using namespace std;

struct Reading{
    int pid;
    int bid;
    int cid;
    int chid;

    double yc;
    double zc;

    int rid;

    double yl;
    double zl;
};

struct ReadingRoot{
    int eid;
    int bid;
    int cid;
    int chid;

    double lev;
    double riv;
    double pep;
    double amp;
};

struct FitData{
    double yc;
    double zc;
    double amp;
    int rid;
};

TFile *f1 = new TFile("N850atm.root");
ifstream mapfile{"cert.map"};

