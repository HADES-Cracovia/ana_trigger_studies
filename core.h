#ifndef FWDET_TESTS_H
#define FWDET_TESTS_H

#include "hades.h"
#include "hloop.h"
#include "htaskset.h"
#include "TString.h"
#include <iostream>

#include "TStopwatch.h"

#include "TH1I.h"
#include "TH2F.h"

#include "hruntimedb.h"
#include "hrecevent.h"
#include "hreconstructor.h"
#include "hcategorymanager.h"
#include "hcategory.h"


//--------category definitions---------
#include "hparticledef.h"
#include "hstartdef.h"
#include "hgeantdef.h"
#include "hpiontrackerdef.h"
//-------------------------------------

//-------objects-----------------------
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticleevtinfo.h"
#include "hstart2hit.h"
#include "hgeantkine.h"
#include "heventheader.h"
#include "hpiontrackertrack.h"
#include "hparticletracksorter.h"

#include "hstart2cal.h"
//-------------------------------------
#include "hphysicsconstants.h"
#include "hparticletool.h"

#include "henergylosscorrpar.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

struct AnaParameters
{
    int verbose_flag;
    TString outfile, config;
    int events;
    float beam_momentum;
    bool sim;
};

typedef std::vector<int> SPVector;

struct GoodTrack
{
    int idx = -1;
    int pid;
    int track_id;
    bool req;
    bool found;
    SPVector search_path;
    string search_str;
    TH1 * hist_theta_all = NULL;
    TH1 * hist_theta_trig = NULL;
    TCanvas * can = NULL;
    void print() const
    {
        printf("%2d  pid=%2d  req=%d  found=%d  ss=%s  id_track=%d\n",
               idx, pid, req, found, search_str.c_str(), track_id);
    }
    void reset()
    {
        track_id = -1;
        found = false;
    }

    void create_gt_histogram()
    {
        char buff[200];
        stringstream ss, ss_can;
        ss << "h_gt_" << idx << "_pid_" << pid;
        ss_can << "c_gt_" << idx << "_pid_" << pid;
        for (auto & x : search_path)
        {
            ss << "_" << (x == -1 ? 0 : x);
            ss_can << "_" << (x == -1 ? 0 : x);
        }
        hist_theta_all = new TH1I(ss.str().c_str(), ";Theta / deg;counts", 360, 0, 180);
        ss << "_trig";
        hist_theta_trig = new TH1I(ss.str().c_str(), ";Theta / deg;counts", 360, 0, 180);

        can = new TCanvas(ss_can.str().c_str(), ss_can.str().c_str(), 800, 600);
    }
};

struct TrackInfo
{
    int pid;
    int parent_pid;
    int track_id;
    int parent_track_id;
};

Int_t core(HLoop * loop, const AnaParameters & anapars);

typedef std::vector<GoodTrack> GTVector;

Bool_t is_good_event(const GTVector & gtv);

#endif /* FWDET_TESTS_H */
