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
    int nomdc_flag;
    int nofdrpc_flag;
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
    bool required;          // required     | fixed
    bool found;             // is found     | by-event

    SPVector search_path;   // decay chain  | fixed
    string search_str;      // decay chain  | fixed
    TH1 * hist_theta_all = NULL;    // hist | accumulated
    TH1 * hist_theta_had = NULL;    // hist | accumulated
    TH1 * hist_theta_fwd = NULL;    // hist | accumulated
    TH1 * hist_theta_acc = NULL;    // host | accumulated

    TH1 * hist_p_all = NULL;    // hist | accumulated
    TH1 * hist_p_had = NULL;    // hist | accumulated
    TH1 * hist_p_fwd = NULL;    // hist | accumulated
    TH1 * hist_p_acc = NULL;    // host | accumulated

    TH2 * hist_p_theta_all = NULL;    // hist | accumulated
    TH2 * hist_p_theta_had = NULL;    // hist | accumulated
    TH2 * hist_p_theta_fwd = NULL;    // hist | accumulated
    TH2 * hist_p_theta_acc = NULL;    // host | accumulated

    TCanvas * can_theta = NULL;   // canvas       | accumulated
    TCanvas * can_p = NULL;   // canvas       | accumulated
    TCanvas * can_p_theta = NULL;   // canvas       | accumulated

    void print() const
    {
        printf("%2d  pid=%2d  req=%d  found=%d  ss=%s  id_track=%d\n",
               idx, pid, required, found, search_str.c_str(), track_id);
    }
    void reset()
    {
        track_id = -1;
        found = false;
    }

    void create_gt_histogram()
    {
        stringstream ss_theta, ss_theta_can;
        stringstream ss_p, ss_p_can;
        stringstream ss_p_theta, ss_p_theta_can;

        ss_theta << "h_gt_theta_" << idx << "_pid_" << pid;
        ss_theta_can << "c_gt_theta_" << idx << "_pid_" << pid;

        ss_p << "h_gt_p_" << idx << "_pid_" << pid;
        ss_p_can << "c_gt_p_" << idx << "_pid_" << pid;

        ss_p_theta << "h_gt_p_theta_" << idx << "_pid_" << pid;
        ss_p_theta_can << "c_gt_p_theta_" << idx << "_pid_" << pid;

        for (auto & x : search_path)
        {
            ss_theta << "_" << (x == -1 ? 0 : x);
            ss_theta_can << "_" << (x == -1 ? 0 : x);

            ss_p << "_" << (x == -1 ? 0 : x);
            ss_p_can << "_" << (x == -1 ? 0 : x);

            ss_p_theta << "_" << (x == -1 ? 0 : x);
            ss_p_theta_can << "_" << (x == -1 ? 0 : x);
        }
        hist_theta_all = new TH1I(ss_theta.str().c_str(), ";Theta / deg;counts", 90, 0, 90);
        hist_theta_had = new TH1I((ss_theta.str() + "_had").c_str(), ";Theta / deg;counts", 90, 0, 90);
        hist_theta_fwd = new TH1I((ss_theta.str() + "_fwd").c_str(), ";Theta / deg;counts", 90, 0, 90);
        hist_theta_acc = new TH1I((ss_theta.str() + "_acc").c_str(), ";Theta / deg;counts", 90, 0, 90);
        can_theta = new TCanvas(ss_theta_can.str().c_str(), ss_theta_can.str().c_str(), 800, 600);

        hist_p_all = new TH1I(ss_p.str().c_str(), ";p / MeV;counts", 300, 0, 3000);
        hist_p_had = new TH1I((ss_p.str() + "_had").c_str(), ";p / MeV;counts", 300, 0, 3000);
        hist_p_fwd = new TH1I((ss_p.str() + "_fwd").c_str(), ";p / MeV;counts", 300, 0, 3000);
        hist_p_acc = new TH1I((ss_p.str() + "_acc").c_str(), ";p / MeV;counts", 300, 0, 3000);
        can_p = new TCanvas(ss_p_can.str().c_str(), ss_p_can.str().c_str(), 800, 600);

        hist_p_theta_all = new TH2I(ss_p_theta.str().c_str(), ";p / MeV;Theta / deg", 300, 0, 3000, 90, 0, 90);
        hist_p_theta_had = new TH2I((ss_p_theta.str() + "_had").c_str(), ";p / MeV;Theta / deg", 300, 0, 3000, 90, 0, 90);
        hist_p_theta_fwd = new TH2I((ss_p_theta.str() + "_fwd").c_str(), ";p / MeV;Theta / deg", 300, 0, 3000, 90, 0, 90);
        hist_p_theta_acc = new TH2I((ss_p_theta.str() + "_acc").c_str(), ";p / MeV;Theta / deg", 300, 0, 3000, 90, 0, 90);
        can_p_theta = new TCanvas(ss_p_theta_can.str().c_str(), ss_p_theta_can.str().c_str(), 800, 600);
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

#endif /* FWDET_TESTS_H */
