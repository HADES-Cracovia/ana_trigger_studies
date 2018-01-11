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

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <string>

struct AnaParameters
{
  TString outfile, config;
  int events;
  float beam_momentum;
  bool sim;
};

struct GoodTrack{
  int pid;
  int track_id;
  bool req;
  bool found;
  std::vector<int> search_path;
  string search_str;
  void print(int idx = -1)
  {
    printf("%2d  pid=%2d  req=%d  found=%d  ss=%s  id_track=%d\n",
		 idx, pid, req, found, search_str.c_str(), track_id);
  }
};

struct TrackInfo{
  int pid;
  int parent_id;
  int track_id;
};

Int_t core(HLoop * loop, const AnaParameters & anapars);

#endif /* FWDET_TESTS_H */
