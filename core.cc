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

    TH1 * h_hit_mult_hades_p = new TH1I("h_hit_mult_hades_p", "Protons in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_pip = new TH1I("h_hit_mult_hades_pip", "pi-plus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_pim = new TH1I("h_hit_mult_hades_pim", "pi-minus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_Kp = new TH1I("h_hit_mult_hades_Kp", "K-plus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades_Km = new TH1I("h_hit_mult_hades_Km", "K-minus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_hades = new TH1I("h_hit_mult_hades", "Charged in Hades mult;multiplicity", 10, 0, 10);

    TH1 * h_hit_mult_fwdet_p = new TH1I("h_hit_mult_fwdet_p", "Protons in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_pip = new TH1I("h_hit_mult_fwdet_pip", "pi-plus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_pim = new TH1I("h_hit_mult_fwdet_pim", "pi-minus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_Kp = new TH1I("h_hit_mult_fwdet_Kp", "K-plus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet_Km = new TH1I("h_hit_mult_fwdet_Km", "K-minus in Hades mult;multiplicity", 10, 0, 10);
    TH1 * h_hit_mult_fwdet = new TH1I("h_hit_mult_fwdet", "Charged in Hades mult;multiplicity", 10, 0, 10);

    TH2 * h_hit_mult_hades_fwdet = new TH2I("h_hit_mult_hades_fwdet", "Charged in Hades/FwDet mult;multiplicity;multiplicity", 10, 0, 10, 10, 0, 10);

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

        Int_t cnt_h_p = 0, cnt_h_p_acc = 0;
        Int_t cnt_h_pip = 0, cnt_h_pip_acc = 0;
        Int_t cnt_h_pim = 0, cnt_h_pim_acc = 0;
        Int_t cnt_h_Kp = 0, cnt_h_Kp_acc = 0;
        Int_t cnt_h_Km = 0, cnt_h_Km_acc = 0;

        Int_t cnt_f_p = 0, cnt_f_p_acc = 0;
        Int_t cnt_f_pip = 0, cnt_f_pip_acc = 0;
        Int_t cnt_f_pim = 0, cnt_f_pim_acc = 0;
        Int_t cnt_f_Kp = 0, cnt_f_Kp_acc = 0;
        Int_t cnt_f_Km = 0, cnt_f_Km_acc = 0;

        Int_t cnt_h = 0, cnt_h_acc = 0;
        Int_t cnt_f = 0, cnt_f_acc = 0;

        n_found = 0;
        n_found_hf = 0;
        n_req = 0;
        n_nreq = 0;

        std::vector<TrackInfo> trackInf;

        for (int j = 0; j < gtnum; ++j)
            goodTracks[j].reset();

//        if (anapars.verbose_flag)
	    printf("Event %d (tracks=%d)\n", i, tracks_num);

        for (int j = 0; j < tracks_num; ++j)
        {
            HGeantKine * pKine = (HGeantKine *)fCatGeantKine -> getObject(j);
            if (!pKine)
                continue;

            int check_nr = 0;

            if (pKine -> getID() == 1)   //no gamma
                continue;

            TrackInfo ti;
            ti.pid = pKine -> getID();
            ti.track_id = j;
            ti.parent_id = pKine -> getParentTrack()-1;
            trackInf.push_back(ti);

            GoodTrack gt;
            gt.found = false;

            for (int k = 0; k < gtnum; ++k)
            {
		if (goodTracks[k].found) continue;
		gt = goodTracks[k];
                gt.track_id = j;
                if (gt.pid != ti.pid) continue;   //search all p,pi-,K+
		TrackInfo t1 = ti;
                for (int l = 0; l < gt.search_path.size(); ++l)
                {
                    if (t1.parent_id == gt.search_path[l])
                    {	
                        if (t1.parent_id == -1)
                        {
                            gt.found = true;
                            goodTracks[k] = gt;
                        }
                        else
                        {
                            t1 = trackInf[t1.parent_id];
                        }
                    }
                    else if (t1.parent_id == -1)
                        break;
		    else
                    {	
                        if (trackInf[ti.parent_id].pid == gt.search_path[l])
                            t1 = trackInf[ti.parent_id];
                        else
                            break;
                    }
                }
                if (gt.found)
                    break;
            }

            //printf("j=%d: allTracks: pid=%d, track_id=%d, req=%d, found=%d\n", j, allTracks.back().pid, allTracks.back().track_id, allTracks.back().req, allTracks.back().found);
            // printf("at: pid=%d, track_id=%d, req=%d, found=%d\n", at.pid, at.track_id, at.req, at.found);

            if (!gt.found)
                continue;

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

            // if (s0 or s1)
            // {
	    switch (pid)
	    {
	    case 8: ++cnt_h_pip; ++cnt_h; break;       // pip
	    case 9: ++cnt_h_pim; ++cnt_h; break;       // pim
	    case 11: ++cnt_h_Kp; ++cnt_h; break;       // Kp
	    case 12: ++cnt_h_Km; ++cnt_h; break;       // Km
	    case 14: ++cnt_h_p; ++cnt_h; break;        // p
	    }
	    // }

            // if (str and rpc)
            // {
	    switch (pid)
	    {
	    case 8: ++cnt_f_pip; ++cnt_f; break;       // pip
	    case 9: ++cnt_f_pim; ++cnt_f; break;       // pim
	    case 11: ++cnt_f_Kp; ++cnt_f; break;       // Kp
	    case 12: ++cnt_f_Km; ++cnt_f; break;       // Km
	    case 14: ++cnt_f_p; ++cnt_f; break;        // p
	    }
	    //   }

//             if (anapars.verbose_flag)
	    printf("  [%03d/%0d] pid=%2d parent=%d  s0=%d s1=%d  f=%d  rpc=%d found=%d\n", j, tracks_num, pKine->getID(), pKine->getParentTrack()-1, s0, s1, str, rpc, gt.found);

            Float_t theta = pKine->getThetaDeg();

	    //theta of good particles in acceptance
	    if ((str and rpc) and (s0 or s1))
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
        h_hit_mult_hades_p -> Fill(cnt_h_p);
        h_hit_mult_hades_pip -> Fill(cnt_h_pip);
        h_hit_mult_hades_pim -> Fill(cnt_h_pim);
        h_hit_mult_hades_Kp -> Fill(cnt_h_Kp);
        h_hit_mult_hades_Km -> Fill(cnt_h_Km);
        h_hit_mult_hades -> Fill(cnt_h);

        h_hit_mult_fwdet_p -> Fill(cnt_f_p);
        h_hit_mult_fwdet_pip -> Fill(cnt_f_pip);
        h_hit_mult_fwdet_pim -> Fill(cnt_f_pim);
        h_hit_mult_fwdet_Kp -> Fill(cnt_f_Kp);
        h_hit_mult_fwdet_Km -> Fill(cnt_f_Km);
        h_hit_mult_fwdet -> Fill(cnt_f);

        h_hit_mult_hades_fwdet -> Fill(cnt_h, cnt_f);

        if (anapars.verbose_flag)
            printf("\n");

        //all counts
        n_h_p += cnt_h_p;
        n_h_p_acc += cnt_h_p_acc;
        n_f_p += cnt_f_p;
        n_f_p_acc += cnt_f_p_acc;
        n_h_pim += cnt_h_pim;
        n_h_pim_acc += cnt_h_pim_acc;

        //counts in event
	// if (anapars.verbose_flag)
        {
            printf("h_pi_good=%d, h_pi_all=%d\n", cnt_h_pim_acc, cnt_h_pim);
            printf("f_pi_good=%d, f_pi_all=%d\n", cnt_f_pim_acc, cnt_f_pim);
            printf("h_p_good=%d,  h_p_all=%d\n", cnt_h_p_acc, cnt_h_p);
            printf("f_p_good=%d,  f_p_all=%d\n", cnt_f_p_acc, cnt_f_p);

            for (int j = 0; j < gtnum; ++j)
            {
                GoodTrack gd = goodTracks[j];
                gd.print(j);
            }
            printf(" Good event ? %s\n", good_event ? "yes" : "no");

            printf("\n********************************************************\n");
        }
    } // end eventloop

    //if (anapars.verbose_flag)
        printf("n_h_p=%d, n_h_p_good=%d\n n_f_p=%d, n_f_p_good=%d\n n_h_pim=%d, n_h_pim_good=%d\n", n_h_p, n_h_p_acc, n_f_p, n_f_p_acc, n_h_pim, n_h_pim_acc);

    h_tr_mult->Write();
    h_hit_s0->Write();
    h_hit_s1->Write();
    h_hit_s01->Write();
    h_hit_str->Write();
    h_hit_rpc->Write();
    h_hit_s01_str->Write();

    h_hit_mult_hades_p->Write();
    h_hit_mult_hades_pip->Write();
    h_hit_mult_hades_pim->Write();
    h_hit_mult_hades_Kp->Write();
    h_hit_mult_hades_Km->Write();
    h_hit_mult_hades->Write();

    h_hit_mult_fwdet_p->Write();
    h_hit_mult_fwdet_pip->Write();
    h_hit_mult_fwdet_pim->Write();
    h_hit_mult_fwdet_Kp->Write();
    h_hit_mult_fwdet_Km->Write();
    h_hit_mult_fwdet->Write();

    h_hit_mult_hades_fwdet->Write();

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

    output_file -> Close();
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();

    return 0;
}

Bool_t is_good_event(const GTVector& gtv)
{
    for (auto x : gtv)
    {
        if (x.req)
            if (!x.found)
                return kFALSE;
    }
    return kTRUE;
}
