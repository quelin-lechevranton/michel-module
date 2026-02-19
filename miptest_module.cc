////////////////////////////////////////////////////////////////////////
// Class:       miptest
// Plugin Type: analyzer (Unknown Unknown)
// File:        miptest_module.cc
//
// Generated at Thu Feb 12 09:05:58 2026 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

namespace ana {
  class miptest;

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  template<typename Data> simb::MCParticle const* trk2mcp(
    art::Ptr<recob::Track> const& pt,
    detinfo::DetectorClocksData const& clockData,
    art::FindManyP<recob::Hit, Data> const& fmp_trk2hit
  ) {
    std::unordered_map<int, float> map_tid_ene;
    for (art::Ptr<recob::Hit> const& p_hit : fmp_trk2hit.at(pt.key()))
      for (sim::TrackIDE ide : bt_serv->HitToTrackIDEs(clockData, p_hit))
        map_tid_ene[ide.trackID] += ide.energy;

    float max_ene = util::kBogusF;
    int tid_max = 0;
    for (std::pair<int, float> p : map_tid_ene)
      if (p.second > max_ene)
        max_ene = p.second, tid_max = p.first;
    return max_ene == util::kBogusF ? nullptr : pi_serv->TrackIdToParticle_P(tid_max);
  }
}


class ana::miptest : public art::EDAnalyzer {
public:
  explicit miptest(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  miptest(miptest const&) = delete;
  miptest(miptest&&) = delete;
  miptest& operator=(miptest const&) = delete;
  miptest& operator=(miptest&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  const detinfo::DetectorPropertiesService* asDetProp;
  const detinfo::DetectorClocksService* asDetClocks;
};


ana::miptest::miptest(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
{
  asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>{};    
  asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>{};
}

void ana::miptest::analyze(art::Event const& e)
{
  auto const clockData = asDetClocks->DataFor(e);
  auto const detProp = asDetProp->DataFor(e,clockData);

  auto const& vh_trk = e.getValidHandle<std::vector<recob::Track>>(art::InputTag("pandoraTrack"));
  auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(art::InputTag("michelprod"));

  if (vh_hit->empty()) return;

  std::vector<art::Ptr<recob::Track>> vp_trk;
  art::fill_ptr_vector(vp_trk, vh_trk);

  std::cout << "\t" << vh_hit->size() << " michel hits in event" << std::endl;

  size_t nhit=0;
  art::FindManyP<recob::Hit> fmp_trk2mi(vh_trk, e, art::InputTag("michelprod"));
  art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, art::InputTag("pandoraTrack"));

  for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
    std::vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk2mi.at(p_trk.key());
    simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);

    if (mcp) std::cout << " track associated with a MCParticle" << std::endl;
    if (vp_hit.empty()) continue;
    std::cout << "\t\t" << vp_hit.size() << " michel hits associated with a track" << std::endl;
    nhit += vp_hit.size();
  }
  if (vh_hit->size() == nhit) {
    std::cout << "\033[92;1m" "all hits were associated" "\033[0m" << std::endl;
  } else {
    std::cout << "\033[91;1m" "number mismatch: " "\033[0m" << vh_hit->size() << " vs. " << nhit << std::endl;
  }
}

DEFINE_ART_MODULE(ana::miptest)
