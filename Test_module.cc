////////////////////////////////////////////////////////////////////////
// Class:       Test
// Plugin Type: analyzer (Unknown Unknown)
// File:        Test_module.cc
//
// Generated at Tue Apr 21 10:11:48 2026 by Jeremy Quelin Lechevranton using cetskelgen
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


#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"

#include <TTree.h>
#include <TBranch.h>

#include <cstdio>
#include <algorithm>




namespace ana {
  class Test;
}


class ana::Test : public art::EDAnalyzer {
public:
  explicit Test(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Test(Test const&) = delete;
  Test(Test&&) = delete;
  Test& operator=(Test const&) = delete;
  Test& operator=(Test&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

  art::InputTag
    tag_pfp{"pandora", ""},
    tag_hit{"hitpdune", ""},
    tag_trk{"pandoraTrack", ""},
    tag_shw{"pandoraShower", ""};
};


ana::Test::Test(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

}

void ana::Test::analyze(art::Event const& e)
{
  auto const& particleHandle = e.getHandle<std::vector<recob::PFParticle>>(tag_pfp);
  if (!particleHandle.isValid()) {
    std::cout << "MiAna: " "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
    return;
  }

  auto const & hitHandle = e.getHandle<std::vector<recob::Hit>>(tag_hit);
  if (!hitHandle.isValid()) {
    std::cout << "MiAna: " "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::Hit>> pHits;
  art::fill_ptr_vector(pHits, hitHandle);

  auto const & trackHandle = e.getHandle<std::vector<recob::Track>>(tag_trk);
  if (!trackHandle.isValid()) {
    std::cout << "MiAna: " "\033[1;91m" "No valid recob::Track handle" "\033[0m" << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::Track>> pTracks;
  art::fill_ptr_vector(pTracks, trackHandle);

  auto const & showerHandle = e.getHandle<std::vector<recob::Shower>>(tag_shw);
  if (!showerHandle.isValid()) {
    std::cout << "MiAna: " "\033[1;91m" "No valid recob::Shower handle" "\033[0m" << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::Shower>> pShowers;
  art::fill_ptr_vector(pShowers, showerHandle);

  art::FindManyP<recob::Hit, recob::TrackHitMeta> trackToManyHits(trackHandle, e, tag_trk);
  art::FindOneP<recob::Track> hitToOneTrack(hitHandle, e, tag_trk);

  art::FindOneP<recob::PFParticle> trackToOneParticle(trackHandle, e, tag_trk);

  art::FindOneP<anab::T0> particleToOneT0(particleHandle, e, tag_pfp);
  art::FindManyP<anab::T0> particleToManyT0(particleHandle, e, tag_pfp);

  art::FindManyP<recob::Hit> showerToManyHits(showerHandle, e, tag_shw);
  art::FindOneP<recob::Shower> hitToOneShower(hitHandle, e, tag_shw);

  art::FindManyP<recob::Wire> hitToManyWires(hitHandle, e, tag_hit);

  // std::vector<size_t> sizes;
  // for (auto const& pHit : pHits) {
  //   size_t size = hitToManyWires.at(pHit.key()).size();

  //   if (std::find(
  //     sizes.begin(), sizes.end(), size
  //   ) != sizes.end()) continue;
  //   sizes.push_back(size);
  // }

  // for (size_t size : sizes) {
  //   std::cout << " " << size;
  // }
  // std::cout << std::endl;



  size_t count=0;

  for (auto const& pTrack : pTracks) {
    std::cout << "track#" << pTrack.id();

    auto const& pParticle = trackToOneParticle.at(pTrack.key());
    if (pParticle.isNull()) {
      std::cout << "\t" "has no pfp" << std::endl;
      continue;
    }
    auto t0s = particleToManyT0.at(pParticle.key());
    if (t0s.empty()) {
      std::cout << "\t" "no t0" << std::endl;
      continue;
    }
    std::cout << "\t" "with " << t0s.size() << " t0s" << std::endl;
    count++;
  }

  std::cout << "tracks with T0: " << count << " among " << pTracks.size() << " tracks" << std::endl;
}

DEFINE_ART_MODULE(ana::Test)
