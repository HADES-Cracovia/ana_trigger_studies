#include "core.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetrpchit.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"

#include "TH2I.h"
#include "TCanvas.h"
#include "fstream"
#include "iostream"

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

Int_t core(HLoop * loop, const AnaParameters & anapars)
{
    if (!loop->setInput(""))
    {   // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    ifstream ifs;
    std::vector<GoodTrack> goodTracks;  //pattertn

    if (anapars.config.Length()>0)
    {
        ifs.open(anapars.config.Data(), std::fstream::in);
        if (!ifs.is_open())
        {
            std::cerr << "Cannot open config file." << endl;
            std::exit(1);
        }
        string path;
        int pid;
        TString path2;

        while (true)
        {
            ifs >> pid >> path;
            if (ifs.eof()) break;
            path2 = path;
            TObjArray *arr = path2.Tokenize(":");

            GoodTrack gd;
            gd.idx = goodTracks.size();
            gd.pid = abs(pid);
            gd.req = pid < 0 ? false : true;
            gd.reset();
            gd.search_str = path;

            int num = arr -> GetEntries();
            for (int i = 0; i < num; ++i)
            {
                TObjString * s = (TObjString *)arr -> At(i);
                int a = s -> GetString().Atoi();
                gd.search_path.push_back(a);
            }
            gd.create_gt_histogram();
            goodTracks.push_back(gd);
        }
    }

    int n_found, n_found_hf, n_req, n_nreq, n_all;

    //    float n_nra, n_fa, h_pim, f_pim, h_p, f_p;
    int n_h_p, n_h_p_acc, n_f_p, n_f_p_acc, n_h_pim, n_h_pim_acc, n_f_pim, n_f_pim_acc;
    n_h_p = 0; n_h_p_acc = 0; n_f_p = 0; n_f_p_acc = 0; n_h_pim = 0; n_h_pim_acc = 0; n_f_pim = 0; n_f_pim_acc = 0;

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    //////////////////////////////////////////////////////////////////////////////
    //      Fast tree builder for creating of ntuples                            //
    //////////////////////////////////////////////////////////////////////////////

    loop->printCategories();    // print all categories found in input + status

    HCategory * fCatGeantKine = nullptr;
    fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");

    if (!fCatGeantKine)
    {
        cout << "No catGeantKine!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    //
    Int_t entries = loop->getEntries();
    //     //setting numbers of events regarding the input number of events by the user
    if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;

    //     // specify output file
    TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
    output_file->cd();
    //

    TH1I * h_tr_mult = new TH1I("h_tr_mult", "Tracks mult;multiplicity", 1000, 0, 1000);
    TH1I * h_hit_s0 = new TH1I("h_hit_s0", "Tracks mult;multiplicity", 2, 0, 2);
    TH1I * h_hit_s1 = new TH1I("h_hit_s1", "Tracks mult;multiplicity", 2, 0, 2);
    TH1I * h_hit_s01 = new TH1I("h_hit_s01", "Tracks mult;multiplicity", 3, 0, 3);
    TH1I * h_hit_str = new TH1I("h_hit_str", "Tracks mult;multiplicity", 20, 0, 20);
    TH1I * h_hit_rpc = new TH1I("h_hit_rpc", "Tracks mult;multiplicity", 10, 0, 10);
    TH2I * h_hit_s01_str = new TH2I("h_hit_s01_str", "Tracks mult;multiplicity", 3, 0, 3, 20, 0, 20);

    TH1 * h_hit_mult_all_p = new TH1I("h_hit_mult_all_p", "Protons all mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_all_pip = new TH1I("h_hit_mult_all_pip", "pi-plus all mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_all_pim = new TH1I("h_hit_mult_all_pim", "pi-minus all mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_all_Kp = new TH1I("h_hit_mult_all_Kp", "K-plus all mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_all_Km = new TH1I("h_hit_mult_all_Km", "K-minus all mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_all = new TH1I("h_hit_mult_all", "Charged all mult;multiplicity", 10, 0, 10);

    // TH1 * h_hit_mult_fwdet_p = new TH1I("h_hit_mult_fwdet_p", "Protons in FT mult;multiplicity", 10, 0, 10);
    // TH1 * h_hit_mult_fwdet_pip = new TH1I("h_hit_mult_fwdet_pip", "pi-plus in FT mult;multiplicity", 10, 0, 10);
    // TH1 * h_hit_mult_fwdet_pim = new TH1I("h_hit_mult_fwdet_pim", "pi-minus in FT mult;multiplicity", 10, 0, 10);
    // TH1 * h_hit_mult_fwdet_Kp = new TH1I("h_hit_mult_fwdet_Kp", "K-plus in FT mult;multiplicity", 10, 0, 10);
    // TH1 * h_hit_mult_fwdet_Km = new TH1I("h_hit_mult_fwdet_Km", "K-minus in FT mult;multiplicity", 10, 0, 10);
    // TH1 * h_hit_mult_fwdet = new TH1I("h_hit_mult_fwdet", "Charged in FT mult;multiplicity", 10, 0, 10);

    // TH2 * h_hit_mult_hades_fwdet = new TH2I("h_hit_mult_hades_fwdet", "Charged in Hades/FwDet mult;multiplicity;multiplicity", 10, 0, 10, 10, 0, 10);

    //only good/accepted events
    TH1I * h_tr_mult_acc = new TH1I("h_tr_mult_acc", "Tracks mult_acc;multiplicity", 1000, 0, 1000);
    TH1I * h_hit_s0_acc = new TH1I("h_hit_s0_acc", "Tracks mult_acc;multiplicity", 2, 0, 2);
    TH1I * h_hit_s1_acc = new TH1I("h_hit_s1_acc", "Tracks mult_acc;multiplicity", 2, 0, 2);
    TH1I * h_hit_s01_acc = new TH1I("h_hit_s01_acc", "Tracks mult_acc;multiplicity", 3, 0, 3);
    TH1I * h_hit_str_acc = new TH1I("h_hit_str_acc", "Tracks mult_acc;multiplicity", 20, 0, 20);
    TH1I * h_hit_rpc_acc = new TH1I("h_hit_rpc_acc", "Tracks mult_acc;multiplicity", 10, 0, 10);
    TH2I * h_hit_s01_str_acc = new TH2I("h_hit_s01_str_acc", "Tracks mult_acc;multiplicity", 3, 0, 3, 20, 0, 20);

    TH1 * h_hit_mult_hades_p_acc = new TH1I("h_hit_mult_hades_p_acc", "Protons in Hades mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_pip_acc = new TH1I("h_hit_mult_hades_pip_acc", "pi-plus in Hades mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_pim_acc = new TH1I("h_hit_mult_hades_pim_acc", "pi-minus in Hades mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_Kp_acc = new TH1I("h_hit_mult_hades_Kp_acc", "K-plus in Hades mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_Km_acc = new TH1I("h_hit_mult_hades_Km_acc", "K-minus in Hades mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_acc = new TH1I("h_hit_mult_hades_acc", "Charged in Hades mult_acc;multiplicity", 10, 0, 10);

    TH1 * h_hit_mult_fwdet_p_acc = new TH1I("h_hit_mult_fwdet_p_acc", "Protons in FT mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_pip_acc = new TH1I("h_hit_mult_fwdet_pip_acc", "pi-plus in FT mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_pim_acc = new TH1I("h_hit_mult_fwdet_pim_acc", "pi-minus in FT mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_Kp_acc = new TH1I("h_hit_mult_fwdet_Kp_acc", "K-plus in FT mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_Km_acc = new TH1I("h_hit_mult_fwdet_Km_acc", "K-minus in FT mult_acc;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_acc = new TH1I("h_hit_mult_fwdet_acc", "Charged in FT mult_acc;multiplicity", 10, 0, 10);

    TH2 * h_hit_mult_hades_fwdet_acc = new TH2I("h_hit_mult_hades_fwdet_acc", "Charged in Hades/FwDet mult_acc;multiplicity;multiplicity", 10, 0, 10, 10, 0, 10);

    // theta of good particles
    TH1I * h_gt_p_theta = new TH1I("h_gt_p_theta", ";proton Theta / rad;counts", 360, 0, 180.0);
    TH1I * h_gt_Kp_theta = new TH1I("h_gt_Kp_theta", ";K+ Theta / rad;counts", 360, 0, 180.0);
    TH1I * h_gt_Km_theta = new TH1I("h_gt_Km_theta", ";K- Theta / rad;counts", 360, 0, 180.0);
    TH1I * h_gt_pip_theta = new TH1I("h_gt_pip_theta", ";#{pi}+ Theta / rad;counts", 360, 0, 180.0);
    TH1I * h_gt_pim_theta = new TH1I("h_gt_pim_theta", ";#{pi}- Theta / rad;counts", 360, 0, 180.0);

    // theta of good particles in acceptance
    TH1I * h_gt_p_theta_acc = new TH1I("h_gt_p_theta_acc", ";proton Theta / rad;counts", 180, 0, 90.0);
    TH1I * h_gt_Kp_theta_acc = new TH1I("h_gt_Kp_theta_acc", ";K+ Theta / rad;counts", 180, 0, 90.0);
    TH1I * h_gt_Km_theta_acc = new TH1I("h_gt_Km_theta_acc", ";K- Theta / rad;counts", 180, 0, 90.0);
    TH1I * h_gt_pip_theta_acc = new TH1I("h_gt_pip_theta_acc", ";#{pi}+ Theta / rad;counts", 180, 0, 90.0);
    TH1I * h_gt_pim_theta_acc = new TH1I("h_gt_pim_theta_acc", ";#{pi}- Theta / rad;counts", 180, 0, 90.0);

    // canvases to comapre both above
    TCanvas * c_gt_p_theta = new TCanvas("c_gt_p_theta", "proton Theta: acc vs all", 800, 600);
    TCanvas * c_gt_Kp_theta = new TCanvas("c_gt_Kp_theta", "K+ Theta: acc vs all", 800, 600);
    TCanvas * c_gt_Km_theta = new TCanvas("c_gt_Km_theta", "K- Theta: acc vs all", 800, 600);
    TCanvas * c_gt_pip_theta = new TCanvas("c_gt_pip_theta", "pi+ Theta: acc vs all", 800, 600);
    TCanvas * c_gt_pim_theta = new TCanvas("c_gt_pim_theta", "pi- Theta: acc vs all", 800, 600);

    const int gtnum = goodTracks.size();

    for (Int_t i = 0; i < entries; ++i)                    // event loop
    {
        /*Int_t nbytes =*/  loop -> nextEvent(i);         // get next event. categories will be cleared before
        //cout << fCatGeantFwDet->getEntries() << endl;

        int tracks_num = fCatGeantKine -> getEntries();
        n_all = tracks_num;
        h_tr_mult -> Fill(tracks_num);

        Int_t cnt_h_p_acc = 0, cnt_f_p_acc = 0;
        Int_t cnt_h_pip_acc = 0, cnt_f_pip_acc = 0;
        Int_t cnt_h_pim_acc = 0, cnt_f_pim_acc = 0;
        Int_t cnt_h_Kp_acc = 0, cnt_f_Kp_acc = 0;
        Int_t cnt_h_Km_acc = 0, cnt_f_Km_acc = 0;

        Int_t cnt_p = 0;
        Int_t cnt_pip = 0;
        Int_t cnt_pim = 0;
        Int_t cnt_Kp = 0;
        Int_t cnt_Km = 0;

        Int_t cnt_h = 0, cnt_h_acc = 0;
        Int_t cnt_f = 0, cnt_f_acc = 0;
    Int_t cnt_all = 0;

        n_found = 0;
        n_found_hf = 0;
        n_req = 0;
        n_nreq = 0;

        std::vector<TrackInfo> trackInf;

        for (int j = 0; j < gtnum; ++j)
            goodTracks[j].reset();


        if (anapars.verbose_flag) printf("Event %d (tracks=%d) ( good tracks: '*' - found, '#' - required )\n", i, tracks_num);

        for (int j = 0; j < tracks_num; ++j)
        {
            HGeantKine * pKine = (HGeantKine *)fCatGeantKine -> getObject(j);
            if (!pKine)
                continue;

            int check_nr = 0;

            TrackInfo ti;
            ti.track_id = j;
            ti.parent_track_id = pKine -> getParentTrack()-1;
            ti.pid = pKine -> getID();
            ti.parent_pid = ti.parent_track_id != -1 ? trackInf[ti.parent_track_id].pid : -1;
            trackInf.push_back(ti);

            if (pKine -> getID() == 1)   //no gamma
                continue;

            GoodTrack gt;
            gt.reset();

            // go over all good tracks (gt)
            for (auto & _gt : goodTracks)
            {
                if (_gt.found) continue;                // if gt was already found, skip it
                if (_gt.pid != ti.pid) continue;        // if gt has PID different that current track (ti), skip
                TrackInfo t1 = ti;                      // local copy of track for reverse search
                for (auto sp_pid : _gt.search_path)     // loop over all pid search path for given gt
                {
                    if (t1.parent_pid == sp_pid)         // if ti has parent ID the same searched PID
                    {
                        if (t1.parent_pid == -1)
                        {
                            _gt.found = true;
                            _gt.track_id = j;
                            break;
                        }
                        else
                        {
                            t1 = trackInf[t1.parent_track_id];
                        }
                    }
                    else if (t1.parent_pid == -1)
                    {
                        break;
                    }
                    else
                    {
                        if (ti.parent_pid == sp_pid)
                            t1 = trackInf[ti.parent_track_id];
                        else
                            break;
                    }
                }
                if (_gt.found)
                {
                    gt = _gt;
                    break;
                }
            }

//             if (anapars.verbose_flag)
//                 printf("  [%03d/%0d]  %c  pid=%2d parent=%d\n",
//                        j, tracks_num,
//                        gt.found ? (gt.req ? '#' : '*' ) : ' ',
//                        pKine->getID(), pKine->getParentTrack()-1);

            if (!gt.found)
                continue;

            if (anapars.verbose_flag)
                printf("  [%03d/%0d]  %c  pid=%2d parent=%d\n",
                       j, tracks_num,
                       gt.found ? (gt.req ? '#' : '*' ) : ' ',
                       pKine->getID(), pKine->getParentTrack()-1);

            gt.hist_theta_all->Fill(pKine->getThetaDeg());

            Int_t pid = pKine -> getID();

            Int_t m0 = 0, m1 = 0, m2 = 0, m3 = 0, s0 = 0, s1 = 0, str = 0, rpc = 0;

            pKine->getNHitsDecayBit(m0, m1, m2, m3, s0, s1);
            pKine->getNHitsFWDecayBit(str, rpc);

            //fill histos for good events in acceptance
            // if (gt.found)
            // {
        // h_hit_s0_acc -> Fill(s0);
        // h_hit_s1_acc -> Fill(s1);
        // h_hit_s01_acc -> Fill(s0+s1);
        // h_hit_str_acc -> Fill(str);
        // h_hit_rpc_acc -> Fill(rpc);
        // h_hit_s01_str_acc -> Fill(s0+s1, str);
        // h_tr_mult_acc -> Fill(tracks_num);

        if (s0 or s1)
        {
        switch (pid)
        {
        case 8: ++cnt_h_pip_acc; ++cnt_h_acc; break;       // pip
        case 9: ++cnt_h_pim_acc; ++cnt_h_acc; break;       // pim
        case 11: ++cnt_h_Kp_acc; ++cnt_h_acc; break;       // Kp
        case 12: ++cnt_h_Km_acc; ++cnt_h_acc; break;       // Km
        case 14: ++cnt_h_p_acc; ++cnt_h_acc; break;        // p
        }
        }

        if (str and rpc)
        {
        switch (pid)
        {
        case 8: ++cnt_f_pip_acc; ++cnt_f_acc; break;       // pip
        case 9: ++cnt_f_pim_acc; ++cnt_f_acc; break;       // pim
        case 11: ++cnt_f_Kp_acc; ++cnt_f_acc; break;       // Kp
        case 12: ++cnt_f_Km_acc; ++cnt_f_acc; break;       // Km
        case 14: ++cnt_f_p_acc; ++cnt_f_acc; break;        // p
        }
        }
        if (gt.req) n_req++;
        if (!gt.req) n_nreq++;
        n_found++;
        // }

            //fill histos for all good events
            h_hit_s0 -> Fill(s0);
            h_hit_s1 -> Fill(s1);
            h_hit_s01 -> Fill(s0+s1);
            h_hit_str -> Fill(str);
            h_hit_rpc -> Fill(rpc);
            h_hit_s01_str -> Fill(s0+s1, str);

        switch (pid)
        {
        case 8: ++cnt_pip; ++cnt_all; break;       // pip
        case 9: ++cnt_pim; ++cnt_all; break;       // pim
        case 11: ++cnt_Kp; ++cnt_all; break;       // Kp
        case 12: ++cnt_Km; ++cnt_all; break;       // Km
        case 14: ++cnt_p; ++cnt_all; break;        // p
        }
       
            if (anapars.verbose_flag)
        printf("  [%03d/%0d] pid=%2d parent=%d  s0=%d s1=%d  f=%d  rpc=%d found=%d\n", j, tracks_num, pKine->getID(), pKine->getParentTrack()-1, s0, s1, str, rpc, gt.found);

            Float_t theta = pKine->getThetaDeg();

        //theta of good particles in acceptance
        if ((str and rpc) or (s0 or s1))
        {
        switch (pid)
        {
        case 8: h_gt_pip_theta_acc->Fill(theta); break;       // pip
                case 9: h_gt_pim_theta_acc->Fill(theta); break;       // pim
                case 11: h_gt_Kp_theta_acc->Fill(theta); break;       // Kp
                case 12: h_gt_Km_theta_acc->Fill(theta); break;       // Km
                case 14: h_gt_p_theta_acc->Fill(theta); break;        // p
        }
        }

        // theta of good particles
        switch (pid)
            {
        case 8: h_gt_pip_theta->Fill(theta); break;       // pip
        case 9: h_gt_pim_theta->Fill(theta); break;       // pim
        case 11: h_gt_Kp_theta->Fill(theta); break;       // Kp
        case 12: h_gt_Km_theta->Fill(theta); break;       // Km
        case 14: h_gt_p_theta->Fill(theta); break;        // p
            }
        }

        Bool_t good_event = is_good_event(goodTracks);

        if (good_event)
        {
            for (auto & gt : goodTracks)
            {
                HGeantKine * pKine = (HGeantKine *)fCatGeantKine -> getObject(gt.track_id);
                if (!pKine)
                    continue;

                gt.hist_theta_trig->Fill(pKine->getThetaDeg());
            }
        }
        
        //fill histos for good events
        h_hit_mult_hades_p_acc -> Fill(cnt_h_p_acc);
        h_hit_mult_hades_pip_acc -> Fill(cnt_h_pip_acc);
        h_hit_mult_hades_pim_acc -> Fill(cnt_h_pim_acc);
        h_hit_mult_hades_Kp_acc -> Fill(cnt_h_Kp_acc);
        h_hit_mult_hades_Km_acc -> Fill(cnt_h_Km_acc);
        h_hit_mult_hades_acc -> Fill(cnt_h_acc);

        h_hit_mult_fwdet_p_acc -> Fill(cnt_f_p_acc);
        h_hit_mult_fwdet_pip_acc -> Fill(cnt_f_pip_acc);
        h_hit_mult_fwdet_pim_acc -> Fill(cnt_f_pim_acc);
        h_hit_mult_fwdet_Kp_acc -> Fill(cnt_f_Kp_acc);
        h_hit_mult_fwdet_Km_acc -> Fill(cnt_f_Km_acc);
        h_hit_mult_fwdet_acc -> Fill(cnt_f_acc);

        h_hit_mult_hades_fwdet_acc -> Fill(cnt_h_acc, cnt_f_acc);

        //foll histots for all events
        h_hit_mult_all_p -> Fill(cnt_p);
        h_hit_mult_all_pip -> Fill(cnt_pip);
        h_hit_mult_all_pim -> Fill(cnt_pim);
        h_hit_mult_all_Kp -> Fill(cnt_Kp);
        h_hit_mult_all_Km -> Fill(cnt_Km);
        h_hit_mult_all -> Fill(cnt_all);

        // h_hit_mult_fwdet_p -> Fill(cnt_f_p);
        // h_hit_mult_fwdet_pip -> Fill(cnt_f_pip);
        // h_hit_mult_fwdet_pim -> Fill(cnt_f_pim);
        // h_hit_mult_fwdet_Kp -> Fill(cnt_f_Kp);
        // h_hit_mult_fwdet_Km -> Fill(cnt_f_Km);
        // h_hit_mult_fwdet -> Fill(cnt_f);

    // h_hit_mult_hades_fwdet -> Fill(cnt_h, cnt_f);

        if (anapars.verbose_flag)
            printf("\n");

        //all counts
        // n_h_p += cnt_h_p;
        // n_h_p_acc += cnt_h_p_acc;
        // n_f_p += cnt_f_p;
        // n_f_p_acc += cnt_f_p_acc;
        // n_h_pim += cnt_h_pim;
        // n_h_pim_acc += cnt_h_pim_acc;

        //counts in event
        if (anapars.verbose_flag)
        {
            printf("h_pi_good=%d, f_pi_good=%d, pi_all=%d\n", cnt_h_pim_acc, cnt_f_pim_acc, cnt_pim);
             printf("h_p_good=%d, f_p_good=%d, p_all=%d\n", cnt_h_p_acc, cnt_f_p_acc, cnt_p);
            
            for (int j = 0; j < gtnum; ++j)
            {
                GoodTrack gd = goodTracks[j];
                gd.print();
            }
            printf(" Good event ? %s\n", good_event ? "yes" : "no");

            printf("\n********************************************************\n");
        }
    } // end eventloop

    // if (anapars.verbose_flag)
    //     printf("n_h_p=%d, n_h_p_good=%d\n n_f_p=%d, n_f_p_good=%d\n n_h_pim=%d, n_h_pim_good=%d\n", n_h_p, n_h_p_acc, n_f_p, n_f_p_acc, n_h_pim, n_h_pim_acc);

    h_tr_mult->Write();
    h_hit_s0->Write();
    h_hit_s1->Write();
    h_hit_s01->Write();
    h_hit_str->Write();
    h_hit_rpc->Write();
    h_hit_s01_str->Write();

    h_hit_mult_all_p->Write();
    h_hit_mult_all_pip->Write();
    h_hit_mult_all_pim->Write();
    h_hit_mult_all_Kp->Write();
    h_hit_mult_all_Km->Write();
    h_hit_mult_all->Write();

    // h_hit_mult_fwdet_p->Write();
    // h_hit_mult_fwdet_pip->Write();
    // h_hit_mult_fwdet_pim->Write();
    // h_hit_mult_fwdet_Kp->Write();
    // h_hit_mult_fwdet_Km->Write();
    // h_hit_mult_fwdet->Write();

    // h_hit_mult_hades_fwdet->Write();

    h_tr_mult_acc -> Write();
    h_hit_s0_acc -> Write();
    h_hit_s1_acc -> Write();
    h_hit_s01_acc -> Write();
    h_hit_str_acc -> Write();
    h_hit_rpc_acc -> Write();
    h_hit_s01_str_acc -> Write();

    h_hit_mult_hades_p_acc -> Write();
    h_hit_mult_hades_pip_acc -> Write();
    h_hit_mult_hades_pim_acc -> Write();
    h_hit_mult_hades_Kp_acc -> Write();
    h_hit_mult_hades_Km_acc -> Write();
    h_hit_mult_hades_acc -> Write();

    h_hit_mult_fwdet_p_acc -> Write();
    h_hit_mult_fwdet_pip_acc -> Write();
    h_hit_mult_fwdet_pim_acc -> Write();
    h_hit_mult_fwdet_Kp_acc -> Write();
    h_hit_mult_fwdet_Km_acc -> Write();
    h_hit_mult_fwdet_acc -> Write();

    h_hit_mult_hades_fwdet_acc -> Write();

    h_gt_p_theta->Write();
    h_gt_Kp_theta->Write();
    h_gt_Km_theta->Write();
    h_gt_pip_theta->Write();
    h_gt_pim_theta->Write();

    h_gt_p_theta_acc->SetLineColor(kRed);
    h_gt_Kp_theta_acc->SetLineColor(kRed);
    h_gt_Km_theta_acc->SetLineColor(kRed);
    h_gt_pip_theta_acc->SetLineColor(kRed);
    h_gt_pim_theta_acc->SetLineColor(kRed);

    h_gt_p_theta_acc->Write();
    h_gt_Kp_theta_acc->Write();
    h_gt_Km_theta_acc->Write();
    h_gt_pip_theta_acc->Write();
    h_gt_pim_theta_acc->Write();

    c_gt_p_theta->cd();
    h_gt_p_theta->Draw();
    h_gt_p_theta_acc->Draw("same");
    c_gt_p_theta->Write();

    c_gt_Kp_theta->cd();
    h_gt_Kp_theta->Draw();
    h_gt_Kp_theta_acc->Draw("same");
    c_gt_Kp_theta->Write();

    c_gt_Km_theta->cd();
    h_gt_Km_theta->Draw();
    h_gt_Km_theta_acc->Draw("same");
    c_gt_Km_theta->Write();

    c_gt_pip_theta->cd();
    h_gt_pip_theta->Draw();
    h_gt_pip_theta_acc->Draw("same");
    c_gt_pip_theta->Write();

    c_gt_pim_theta->cd();
    h_gt_pim_theta->Draw();
    h_gt_pim_theta_acc->Draw("same");
    c_gt_pim_theta->Write();

    for (auto & x : goodTracks)
    {
        x.hist_theta_all->SetLineColor(kBlack);
        x.hist_theta_all->Write();
        x.hist_theta_trig->SetLineColor(kRed);
        x.hist_theta_trig->Write();
        TCanvas * c = x.can;
        c->cd();
        x.hist_theta_all->Draw();
        x.hist_theta_trig->Draw("same");
        c->Write();
    }

    output_file -> Close();
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();

    return 0;
}

Bool_t is_good_event(const GTVector& gtv)
{
    for (auto & x : gtv)
    {
        if (x.req)
            if (!x.found)
                return kFALSE;
    }
    return kTRUE;
}
