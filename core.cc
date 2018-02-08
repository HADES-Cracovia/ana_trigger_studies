#include "core.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetrpchit.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"

#include "TH2I.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

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
            gd.required = pid < 0 ? false : true;
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

    int n_found, n_req, n_nreq;

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

    TH2 * h_hit_mult_hades_fwdet_req_acc = new TH2I("h_hit_mult_hades_fwdet_req_acc", "Charged in Hades/FwDet mult_acc;Hades;FwDet", 8, 0, 8, 8, 0, 8);
    TH2 * h_hit_mult_hades_fwdet_acc = new TH2I("h_hit_mult_hades_fwdet_acc", "Charged in Hades/FwDet mult_acc;Hades;FwDet", 8, 0, 8, 8, 0, 8);

    TCanvas * c_hit_mult_hades_fwdet_req_acc = new TCanvas("c_hit_mult_hades_fwdet_req_acc", "Required tracks multiplicity in acceptance Hades/FwDet", 800, 600);
    TCanvas * c_hit_mult_hades_fwdet_acc = new TCanvas("c_hit_mult_hades_fwdet_acc", "Tracks multiplicity in acceptance Hades/FwDet", 800, 600);

    TH1 * h_gt_mult_acc = new TH1I("h_gt_mult_acc", ";Tracks multiplicity in acceptance;counts", 10, 0, 10);
    TCanvas * c_gt_mult_acc = new TCanvas("c_gt_mult_acc", "Tracks multiplicity in acceptance", 800, 600);

    const int gtnum = goodTracks.size();
    int n_all_ge = 0, n_acc_ge = 0;
    // int n_all_0 = 0, n_all_1 = 0, n_all_2 = 0;
    // int n_acc_0 = 0, n_acc_1 = 0, n_acc_2 = 0;
    // int n_facc_0 = 0, n_facc_1 = 0, n_facc_2 = 0;
    // int n_hacc_0 = 0, n_hacc_1 = 0, n_hacc_2 = 0;
    
    for (Int_t i = 0; i < entries; ++i)                    // event loop
    {
        /*Int_t nbytes =*/  loop -> nextEvent(i);         // get next event. categories will be cleared before
        //cout << fCatGeantFwDet->getEntries() << endl;

        int tracks_num = fCatGeantKine -> getEntries();
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

        Int_t cnt_h_acc = 0;
        Int_t cnt_f_acc = 0;

        Int_t cnt_h_req_acc = 0;
        Int_t cnt_f_req_acc = 0;

        Int_t cnt_all = 0;

        n_found = 0;
        n_req = 0;
        n_nreq = 0;

        std::vector<TrackInfo> trackInf;
        trackInf.reserve(10000);

        for (int j = 0; j < gtnum; ++j)
            goodTracks[j].reset();


        if (anapars.verbose_flag) printf("Event %d (tracks=%d) ( good tracks: '*' - found, '#' - required )\n", i, tracks_num);

        for (int j = 0; j < tracks_num; ++j)
        {
            HGeantKine * pKine = (HGeantKine *)fCatGeantKine -> getObject(j);
            if (!pKine)
                continue;

            TrackInfo ti;
            ti.track_id = pKine->getTrack()-1;
            ti.parent_track_id = pKine -> getParentTrack()-1;
            ti.pid = pKine -> getID();
            ti.parent_pid = ti.parent_track_id != -1 ? trackInf[ti.parent_track_id].pid : -1;

            // if gamma, push it to the vector and continue
            if (pKine -> getID() == 1)   //no gamma
            {
                trackInf.push_back(ti);
                continue;
            }

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
                            _gt.track_id = pKine->getTrack()-1;
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

//             if (anapars.verbose_flag)
//                 printf("  [%03d/%0d] pid=%2d parent=%d  found=%d\n", ti.track_id, tracks_num, pKine->getID(), pKine->getParentTrack()-1, gt.found);

            // if not good track, push it to the vector and continue
            if (!gt.found)
            {
                trackInf.push_back(ti);
                continue;
            }

            // check if track in aceptance
            Int_t m0 = 0, m1 = 0, m2 = 0, m3 = 0, s0 = 0, s1 = 0, str = 0, rpc = 0;
            pKine->getNHitsDecayBit(m0, m1, m2, m3, s0, s1);
            pKine->getNHitsFWDecayBit(str, rpc);

            pKine->fillAcceptanceBit();
            Int_t nrpc = 0;
            Bool_t is_good_fwdet_acc = pKine->isInTrackAcceptanceFWDecayBit(nrpc);

            // is hit in hades and fwdet
            //                |------------------- HADES --------------------|
            //                |- META -|     |------------- MDC -------------|
            if (anapars.nomdc_flag)
                ti.is_hades_hit = (s0 or s1);
            else
                ti.is_hades_hit = (s0 or s1) and (m0>0 and m1>0 and m2>0 and m3>0);

            if (anapars.nofdrpc_flag)
                ti.is_fwdet_hit = is_good_fwdet_acc;
            else
                ti.is_fwdet_hit = is_good_fwdet_acc and (nrpc > 0);

            ti.is_in_acc = ti.is_hades_hit or (ti.is_fwdet_hit);
            trackInf.push_back(ti);

            if (anapars.verbose_flag)
                printf("  [%03d/%0d]  %c  pid=%2d parent=%d\n",
                       j, tracks_num,
                       gt.found ? (gt.required ? '#' : '*' ) : ' ',
                       pKine->getID(), pKine->getParentTrack()-1);

            Float_t theta = pKine->getThetaDeg();
            Float_t p = pKine->getTotalMomentum();

	    /*  // fill all tracks that are good tracks
            gt.hist_theta_all->Fill(theta);
            gt.hist_p_all->Fill(p);
            gt.hist_p_theta_all->Fill(p, theta);
	    */
//             if (ti.hades_hit)
//                 gt.hist_theta_had->Fill(pKine->getThetaDeg());
//             if (is_good_fwdet_acc)
//                 gt.hist_theta_fwd->Fill(pKine->getThetaDeg());

            Int_t pid = pKine -> getID();

            pKine->getNHitsDecayBit(m0, m1, m2, m3, s0, s1);
            pKine->getNHitsFWDecayBit(str, rpc);

            if (gt.required) n_req++;
            if (!gt.required) n_nreq++;
            n_found++;

            //fill histos for all good events
            h_hit_s0 -> Fill(s0);
            h_hit_s1 -> Fill(s1);
            h_hit_s01 -> Fill(s0+s1);
            h_hit_str -> Fill(str);
            h_hit_rpc -> Fill(rpc);
            h_hit_s01_str -> Fill(s0+s1, str);

	    /* switch (pid)
            {
                case 8: ++cnt_pip; ++cnt_all; break;       // pip
                case 9: ++cnt_pim; ++cnt_all; break;       // pim
                case 11: ++cnt_Kp; ++cnt_all; break;       // Kp
                case 12: ++cnt_Km; ++cnt_all; break;       // Km
                case 14: ++cnt_p; ++cnt_all; break;        // p
		}
	    */
            // if (anapars.verbose_flag)
            //     printf("  [%03d/%0d] pid=%2d parent=%d  s0=%d s1=%d  f=%d  rpc=%d found=%d\n", j, tracks_num, pKine->getID(), pKine->getParentTrack()-1, s0, s1, str, rpc, gt.found);
        }

        // Good event is one, where all required tracks are found
        Bool_t good_event = is_good_event(goodTracks);

	if (good_event)
        {
	    Int_t cnt_all_ge = 0;
	    for (auto & gt : goodTracks)
	    {
		HGeantKine * pKine = (HGeantKine *)fCatGeantKine -> getObject(gt.track_id);
		if (!pKine)
		    continue;

		Float_t theta = pKine->getThetaDeg();
		Float_t p = pKine->getTotalMomentum();
		Int_t pid = pKine -> getID();
		
		// fill all tracks that are good tracks
		if(gt.found)
		{
		    gt.hist_theta_all->Fill(theta);
		    gt.hist_p_all->Fill(p);
		    gt.hist_p_theta_all->Fill(p, theta);

		    switch (pid)
		    {
		    case 8: ++cnt_pip; ++cnt_all; break;       // pip
		    case 9: ++cnt_pim; ++cnt_all; break;       // pim
		    case 11: ++cnt_Kp; ++cnt_all; break;       // Kp
		    case 12: ++cnt_Km; ++cnt_all; break;       // Km
		    case 14: ++cnt_p; ++cnt_all; break;        // p
		    }
		    ++n_all_ge;
		    ++cnt_all_ge;
		}
		
		}

	   // If we have good event, check if all required tracks are in the acceptance
            Bool_t all_req_in_acc = is_good_event_in_acc(goodTracks, trackInf);

            if (all_req_in_acc)
            {
                Int_t cnt = 0;
                for (auto & gt : goodTracks)
                {
                    HGeantKine * pKine = (HGeantKine *)fCatGeantKine -> getObject(gt.track_id);
                    if (!pKine)
                        continue;

                    Float_t theta = pKine->getThetaDeg();
                    Float_t p = pKine->getTotalMomentum();

                    TrackInfo & ti = trackInf[gt.track_id];
		    if (ti.is_fwdet_hit)
                    {
                        gt.hist_theta_fwd->Fill(theta);
                        gt.hist_p_fwd->Fill(p);
                        gt.hist_p_theta_fwd->Fill(p, theta);

                        ++cnt_f_acc;
                        if (gt.required) ++cnt_f_req_acc;
                    }
		    else if (ti.is_hades_hit)
                    {
                        gt.hist_theta_had->Fill(theta);
                        gt.hist_p_had->Fill(p);
                        gt.hist_p_theta_had->Fill(p, theta);

                        ++cnt_h_acc;
                        if (gt.required) ++cnt_h_req_acc;
                    }
                   
                    gt.hist_theta_acc->Fill(theta);
                    gt.hist_p_acc->Fill(p);
                    gt.hist_p_theta_acc->Fill(p, theta);

                    if (ti.is_in_acc)
		    {
                        ++cnt;
			++n_acc_ge;
		    }
                }
                h_gt_mult_acc->Fill(cnt);

                h_hit_mult_hades_fwdet_req_acc -> Fill(cnt_h_req_acc, cnt_f_req_acc);
                h_hit_mult_hades_fwdet_acc -> Fill(cnt_h_acc, cnt_f_acc);
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

        //fill histots for all events
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

    /* h_tr_mult->Write();
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
    h_hit_mult_fwdet_acc -> Write();*/

    h_hit_mult_hades_fwdet_req_acc->SetMarkerColor(kWhite);
    h_hit_mult_hades_fwdet_req_acc->SetMarkerSize(2);
    h_hit_mult_hades_fwdet_req_acc->Write();

    h_hit_mult_hades_fwdet_acc->SetMarkerColor(kWhite);
    h_hit_mult_hades_fwdet_acc->SetMarkerSize(2);
    h_hit_mult_hades_fwdet_acc->Write();

    c_hit_mult_hades_fwdet_req_acc->cd();
    h_hit_mult_hades_fwdet_req_acc->Draw("colz,text30");
    c_hit_mult_hades_fwdet_req_acc->Write();

    c_hit_mult_hades_fwdet_acc->cd();
    h_hit_mult_hades_fwdet_acc->Draw("colz,text30");
    c_hit_mult_hades_fwdet_acc->Write();

    TLatex * tex = new TLatex;
    tex->SetNDC(kTRUE);

    for (auto & x : goodTracks)
    {
        // theta
        x.hist_theta_all->SetLineColor(kBlack);
        x.hist_theta_all->Write();
        x.hist_theta_acc->SetLineWidth(0);
        x.hist_theta_acc->SetLineColor(0);
        if (x.required)
            x.hist_theta_acc->SetFillColor(46);
        else
            x.hist_theta_acc->SetFillColor(41);
        x.hist_theta_acc->Write();
        x.hist_theta_had->SetLineColor(30);
        x.hist_theta_had->SetLineWidth(2);
        x.hist_theta_had->Write();
        x.hist_theta_fwd->SetLineColor(38);
        x.hist_theta_fwd->SetLineWidth(2);
        x.hist_theta_fwd->Write();
        TCanvas * c = x.can_theta;
        c->cd();
        x.hist_theta_all->Draw();
        x.hist_theta_acc->Draw("same");
        x.hist_theta_fwd->Draw("same");
        x.hist_theta_had->Draw("same");

        TLegend * leg = new TLegend(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(x.hist_theta_all, "Good tracks", "lp");
        leg->AddEntry(x.hist_theta_had, "required in acceptance in hades", "lp");
        leg->AddEntry(x.hist_theta_fwd, "required in acceptance in fwd", "lp");
        leg->AddEntry(x.hist_theta_acc, "required in acceptance", "lp");
        leg->Draw();
        tex->DrawLatex(0.6, 0.70, TString::Format("# all = %.0f", x.hist_theta_all->Integral()));
        tex->DrawLatex(0.6, 0.65, TString::Format("# acc = %.0f", x.hist_theta_acc->Integral()));
        
        tex->DrawLatex(0.6, 0.55, TString::Format("# fwd = %.0f", x.hist_theta_fwd->Integral()));
        tex->DrawLatex(0.6, 0.50, TString::Format("# had = %.0f", x.hist_theta_had->Integral()));
	// tex->DrawLatex(0.6, 0.70, TString::Format("# all = %.0f", cnt_all_ge));
        // tex->DrawLatex(0.6, 0.65, TString::Format("# acc = %.0f", cnt));

        // tex->DrawLatex(0.6, 0.55, TString::Format("# fwd = %.0f", cnt_f_acc));
	// tex->DrawLatex(0.6, 0.50, TString::Format("# had = %.0f", cnt_h_acc));
        c->SetLogy();
        c->Write();

        // p
        x.hist_p_all->SetLineColor(kBlack);
        x.hist_p_all->Write();
        x.hist_p_acc->SetLineWidth(0);
        x.hist_p_acc->SetLineColor(0);
        if (x.required)
            x.hist_p_acc->SetFillColor(46);
        else
            x.hist_p_acc->SetFillColor(41);
        x.hist_p_acc->Write();
        x.hist_p_had->SetLineColor(30);
        x.hist_p_had->SetLineWidth(2);
        x.hist_p_had->Write();
        x.hist_p_fwd->SetLineColor(38);
        x.hist_p_fwd->SetLineWidth(2);
        x.hist_p_fwd->Write();
        c = x.can_p;
        c->cd();
        x.hist_p_all->Draw();
        x.hist_p_acc->Draw("same");
        x.hist_p_fwd->Draw("same");
        x.hist_p_had->Draw("same");

        leg->Draw();
        tex->DrawLatex(0.6, 0.70, TString::Format("# all = %.0f", x.hist_p_all->Integral()));
        tex->DrawLatex(0.6, 0.65, TString::Format("# acc = %.0f", x.hist_p_acc->Integral()));
        
        tex->DrawLatex(0.6, 0.55, TString::Format("# fwd = %.0f", x.hist_p_fwd->Integral()));
        tex->DrawLatex(0.6, 0.50, TString::Format("# had = %.0f", x.hist_p_had->Integral()));
        c->SetLogy();
        c->Write();

        // p-theta
        x.hist_p_theta_all->Write();
        x.hist_p_theta_acc->Write();
        x.hist_p_theta_had->Write();
        x.hist_p_theta_fwd->Write();
        c = x.can_p_theta;
        c->cd();
        x.hist_p_theta_acc->Draw("colz");
        c->Write();
    }
    h_gt_mult_acc->SetMarkerSize(2);
    h_gt_mult_acc->Write();
    c_gt_mult_acc->cd();
    h_gt_mult_acc->Draw("h,text30");
    c_gt_mult_acc->Write();

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
        if (x.required)
            if (!x.found)
                return kFALSE;
    }
    return kTRUE;
}

// This function assumes that event is good, no further check. Call is_good_event() before.
Bool_t is_good_event_in_acc(const GTVector& gtv, const TIVector& ti)
{
    for (auto & x : gtv)
    {
        if (x.required)
            if (!ti[x.track_id].is_in_acc)
                return kFALSE;
    }
    return kTRUE;
}
