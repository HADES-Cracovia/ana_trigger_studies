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
    int idx = -1;           // #index       | fixed

    int pid;                // pid          | fixed
    int track_id;           // track id     | by-event
    bool req;               // required     | fixed
    bool found;             // is found     | by-event

    SPVector search_path;   // decay chain  | fixed
    string search_str;      // decay chain  | fixed
    TH1 * hist_theta_all = NULL;    // hist | accumulated
    TH1 * hist_theta_had = NULL;    // hist | accumulated
    TH1 * hist_theta_fwd = NULL;    // hist | accumulated
    TH1 * hist_theta_acc = NULL;    // host | accumulated
    TCanvas * can = NULL;   // canvas       | accumulated

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
        hist_theta_all = new TH1I(ss.str().c_str(), ";Theta / deg;counts", 180, 0, 180);
        hist_theta_had = new TH1I((ss.str() + "_had").c_str(), ";Theta / deg;counts", 180, 0, 180);
        hist_theta_fwd = new TH1I((ss.str() + "_fwd").c_str(), ";Theta / deg;counts", 180, 0, 180);
        hist_theta_acc = new TH1I((ss.str() + "_acc").c_str(), ";Theta / deg;counts", 180, 0, 180);

        can = new TCanvas(ss_can.str().c_str(), ss_can.str().c_str(), 800, 600);
    }
};

struct TrackInfo
{
    int pid = -1;
    int parent_pid = -1;
    int track_id = -1;
    int parent_track_id = -1;

    bool is_hades_hit = false;
    bool is_fwdet_hit = false;

    bool is_in_acc = false;
};

Int_t core(HLoop * loop, const AnaParameters & anapars);

typedef std::vector<GoodTrack> GTVector;
typedef std::vector<TrackInfo> TIVector;

Bool_t is_good_event(const GTVector & gtv);
Bool_t is_good_event_in_acc(const GTVector & gtv, const TIVector & tiv);

Bool_t is_full_straw_track(HGeantKine * pKine);

#endif /* FWDET_TESTS_H */
