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
    {                                                    // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    ifstream ifs;
    std::vector<GoodTrack> goodTracks;
    if(anapars.config.Length()>0){
      ifs.open(anapars.config.Data(), std::fstream::in);
      if(!ifs.is_open()){
	std::cerr << "Cannot open config file." << endl;
	std::exit(1);
      }
      string path;
      int pid;
      TString path2;

      while(true){
	ifs >> pid >> path;
	if(ifs.eof()) break;
	path2 = path;
	TObjArray *arr = path2.Tokenize(":");

	GoodTrack gd;
	gd.pid = abs(pid);
	gd.req = pid < 0 ? false : true;
	gd.reset();
	gd.search_str = path;

	int num = arr -> GetEntries();
	for(int i = 0; i < num; i++){
	  TObjString * s = (TObjString *)arr->At(i);
	  int a = s->GetString().Atoi();
	  gd.search_path.push_back(a);
	}
	goodTracks.push_back(gd);
      }
    }

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

    for (Int_t i = 0; i < entries; i++)                    // event loop
    {
        /*Int_t nbytes =*/  loop->nextEvent(i);         // get next event. categories will be cleared before
        //cout << fCatGeantFwDet->getEntries() << endl;

        int tracks_num = fCatGeantKine->getEntries();
        h_tr_mult->Fill(tracks_num);

        Int_t cnt_h_p = 0;
        Int_t cnt_h_pip = 0;
        Int_t cnt_h_pim = 0;
        Int_t cnt_h_Kp = 0;
        Int_t cnt_h_Km = 0;

        Int_t cnt_f_p = 0;
        Int_t cnt_f_pip = 0;
        Int_t cnt_f_pim = 0;
        Int_t cnt_f_Kp = 0;
        Int_t cnt_f_Km = 0;

        Int_t cnt_h = 0;
        Int_t cnt_f = 0;

	std::vector<TrackInfo> trackInf;
	const int gdnum = goodTracks.size();

	for (int j = 0; j < gdnum; ++j)
	  goodTracks[j].reset();

        printf("Event %d\n", i);
        for (int j = 0; j < tracks_num; ++j)
        {
            HGeantKine * pKine = (HGeantKine *)fCatGeantKine->getObject(j);
            if (!pKine)
                continue;
	    
	    TrackInfo ti;
	    ti.pid = pKine->getID();
	    ti.track_id = j;
	    ti.parent_id = pKine->getParentTrack()-1;
	    trackInf.push_back(ti);
	    
	    if(pKine->getID() == 1)
	      continue;

	    for(int k = 0; k < gdnum; k++){
	      GoodTrack gd = goodTracks[k];
	      //gd.print(k);
	      if(gd.found) continue;
	      if(gd.pid != ti.pid) continue;
	      TrackInfo t1 = ti;
	      for(int l = 0; l < gd.search_path.size(); l++){
		printf("%d    %d\n", t1.parent_id, gd.search_path[l]);
		if(t1.parent_id == gd.search_path[l]){
		  if(t1.parent_id == -1){
		    gd.found = true;
		    gd.track_id = j;
		    goodTracks[k]= gd;
		  }
		  else{
		    t1 = trackInf[t1.parent_id];
		  }
		}
		else if (t1.parent_id == -1)
		{
		  break;
		}
		else
		{
		  if (trackInf[t1.parent_id].pid == gd.search_path[l])
		    t1 = trackInf[t1.parent_id];
		  else
		    break;
		}
	      }
	      if (gd.found) break;
	    }
	    
            Int_t m0 = 0, m1 = 0, m2 = 0, m3 = 0, s0 = 0, s1 = 0, str = 0, rpc = 0;

            pKine->getNHitsDecayBit(m0, m1, m2, m3, s0, s1);
            pKine->getNHitsFWDecayBit(str, rpc);
            h_hit_s0->Fill(s0);
            h_hit_s1->Fill(s1);
            h_hit_s01->Fill(s0+s1);
            h_hit_str->Fill(str);
            h_hit_rpc->Fill(rpc);
            h_hit_s01_str->Fill(s0+s1, str);

            if (s0 or s1)
            {
                Int_t pid = pKine->getID();
                switch (pid)
                {
                    case 8: ++ cnt_h_pip; ++cnt_h; break;       // pip
                    case 9: ++ cnt_h_pim; ++cnt_h; break;       // pim
                    case 11: ++ cnt_h_Kp; ++cnt_h; break;       // Kp
                    case 12: ++ cnt_h_Km; ++cnt_h; break;       // Km
                    case 14: ++ cnt_h_p; ++cnt_h; break;        // p
                }
            }

            if (str and rpc)
            {
                Int_t pid = pKine->getID();
                switch (pid)
                {
                    case 8: ++ cnt_f_pip; ++cnt_f; break;       // pip
                    case 9: ++ cnt_f_pim; ++cnt_f; break;       // pim
                    case 11: ++ cnt_f_Kp; ++cnt_f; break;       // Kp
                    case 12: ++ cnt_f_Km; ++cnt_f; break;       // Km
                    case 14: ++ cnt_f_p; ++cnt_f; break;        // p
                }
            }

            printf("  [%03d/%0d] pid=%2d parent=%d  s0=%d s1=%d  f=%d  rpc=%d\n", j, tracks_num, pKine->getID(), pKine->getParentTrack()-1, s0, s1, str, rpc);
        }

        h_hit_mult_hades_p->Fill(cnt_h_p);
        h_hit_mult_hades_pip->Fill(cnt_h_pip);
        h_hit_mult_hades_pim->Fill(cnt_h_pim);
        h_hit_mult_hades_Kp->Fill(cnt_h_Kp);
        h_hit_mult_hades_Km->Fill(cnt_h_Km);
        h_hit_mult_hades->Fill(cnt_h);

        h_hit_mult_fwdet_p->Fill(cnt_f_p);
        h_hit_mult_fwdet_pip->Fill(cnt_f_pip);
        h_hit_mult_fwdet_pim->Fill(cnt_f_pim);
        h_hit_mult_fwdet_Kp->Fill(cnt_f_Kp);
        h_hit_mult_fwdet_Km->Fill(cnt_f_Km);
        h_hit_mult_fwdet->Fill(cnt_f);

        h_hit_mult_hades_fwdet->Fill(cnt_h, cnt_f);

        printf("\n");

	for(int j = 0; j < gdnum; j++){
	  GoodTrack gd = goodTracks[j];
	  gd.print(j);
	}
	printf("\n\n\n\n");
    } // end eventloop

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

    output_file->Close();
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();

    return 0;
}
