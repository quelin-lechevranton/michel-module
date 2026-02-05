////////////////////////////////////////////////////////////////////////
// Class:       MichelProd
// Plugin Type: producer (Unknown Unknown)
// File:        MichelProd_module.cc
//
// Generated at Tue Feb  3 04:31:11 2026 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "TTree.h"

#include <memory>

using PtrHit    = art::Ptr<recob::Hit>;
using VecPtrHit = std::vector<art::Ptr<recob::Hit>>;
using PtrTrk    = art::Ptr<recob::Track>;
using VecPtrTrk = std::vector<art::Ptr<recob::Track>>;

namespace ana {
  class MichelProd;


  // class Michel {
  // private:
  //   art::Ptr<recob::Hit> fMuonLastHit;

  // public:
  //   Michel(art::Ptr<recob::Hit> muonLastHit): fMuonLastHit(muonLastHit) {}

  //   inline art::Ptr<recob::Hit> MuonLastHit() const { return fMuonLastHit; } 
  // };
}

class ana::MichelProd : public art::EDProducer {
public:
  explicit MichelProd(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelProd(MichelProd const&) = delete;
  MichelProd(MichelProd&&) = delete;
  MichelProd& operator=(MichelProd const&) = delete;
  MichelProd& operator=(MichelProd&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  const geo::GeometryCore* asGeo;
  const geo::WireReadoutGeom* asWire;
  const detinfo::DetectorPropertiesService* asDetProp;
  const detinfo::DetectorClocksService* asDetClocks;
  ana::Det_t geoDet;
  art::InputTag tag_hit, tag_trk, tag_cal;
  float fTick2cm;

  TTree *tuple;
  std::vector<float> residual_range, dEds;

  float GetDistance(PtrHit const& ph1, PtrHit const& ph2);
};

void ana::MichelProd::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  tuple = tfs->make<TTree>("calo", "");
  tuple->Branch("ResidualRange", &residual_range);
  tuple->Branch("dEds", &dEds);
}


ana::MichelProd::MichelProd(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , tag_hit{"hitpdune"}
  , tag_trk{"pandoraTrack"}
  , tag_cal{"pandoraGnocalo"}
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  // produces<std::vector<ana::Michel>>();
  // produces<art::Assns<recob::Track, ana::Michel>>();
  produces<std::vector<recob::Hit>>();

  asGeo = &*art::ServiceHandle<geo::Geometry>{};
  asWire = &art::ServiceHandle<geo::WireReadout>{}->Get();
  asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>{};    
  asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>{};
  if (asGeo->DetectorName().find("vd") != std::string::npos)
      geoDet = kPDVD;
  else if (asGeo->DetectorName().find("hd") != std::string::npos)
      geoDet = kPDHD;
  else {
      std::cout << "\033[1;91m" "unknown geometry: "
          << asGeo->DetectorName() << "\033[0m" << std::endl;
      exit(1);
  }

  auto const clockData = asDetClocks->DataForJob();
  auto const detProp = asDetProp->DataForJob(clockData);
  fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
}

void ana::MichelProd::produce(art::Event& e)
{
  // Output
  // auto outMichels = std::make_unique<std::vector<ana::Michel>>();
  // auto outAssns = std::make_unique<art::Assns<recob::Track, ana::Michel>>();
  // art::PtrMaker<ana::Michel> michelPtrMaker(e);
  auto outMichelHits = std::make_unique<std::vector<recob::Hit>>();

  // Input
  auto const& vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
  if (!vh_hit.isValid()) return;
  VecPtrHit vph_ev; art::fill_ptr_vector(vph_ev, vh_hit);

  auto const& vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
  if (!vh_trk.isValid()) return;
  VecPtrTrk vpt_ev; art::fill_ptr_vector(vpt_ev, vh_trk);

  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
  art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
  art::FindManyP<anab::Calorimetry> fop_trk2cal(vh_trk, e, tag_cal);

  for (PtrTrk const& pt_ev: vpt_ev) {
    if (pt_ev->Length() > 30) continue;

    VecPtrHit vph_trk = fmp_trk2hit.at(pt_ev.key());
    if (vph_trk.empty()) continue;

    std::vector<recob::TrackHitMeta const*> const& vhm_trk = fmp_trk2hit.data(pt_ev.key());
    if (vph_trk.size() != vhm_trk.size()) continue;

    std::vector<art::Ptr<anab::Calorimetry>> vpc_trk = fop_trk2cal.at(pt_ev.key());
    std::cout << "number of cal in tracks: " << vpc_trk.size() << std::endl;
    if (vpc_trk.empty()) continue;
    art::Ptr<anab::Calorimetry> cal_trk = vpc_trk.front();
    residual_range = cal_trk->ResidualRange();

    // get Track Index of the calo record with minimum residual range
    size_t lastCal = std::distance(residual_range.begin(), std::min_element(
      residual_range.begin(), residual_range.end()
    ));
    if (lastCal == residual_range.size()) continue;
    size_t lastIndex = cal_trk->TpIndices().at(lastCal);
    if (!pt_ev->HasValidPoint(lastIndex)) continue;

    // get collection Hit Index corresponding to the Track Index
    size_t lastHitIndex = vph_trk.size();
    for (size_t i=0; i<vph_trk.size(); i++) {
      if (vph_trk.at(i)->View() != geo::kW) continue;
      if (vhm_trk.at(i)->Index() != lastIndex) continue;
      lastHitIndex = i;
      break;
    }
    if (lastHitIndex == vph_trk.size()) continue;
    PtrHit const& lastHit = vph_trk.at(lastHitIndex);

    // NTuple output
    // residual_range = cal_trk->ResidualRange();
    dEds = cal_trk->dEdx();
    tuple->Fill();
    
    // Producer outputs
    // ana::Michel outMichel(lastHit);
    // outMichels->emplace_back(std::move(outMichel));
    // art::Ptr<ana::Michel> ptr_outMichel = michelPtrMaker(outMichels->size()-1);
    // outAssns->addSingle(pt_ev, ptr_outMichel);

    for (PtrHit const& ph_ev: vph_ev) {
      if (GetDistance(ph_ev, lastHit) > 30) continue;

      PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
      if (pt_hit && pt_hit->Length() > 30) continue;

      PtrHit outMichelHit(ph_ev);
      outMichelHits->emplace_back(std::move(outMichelHit));
    }
  }

  // Put outputs in event
  // e.put(std::move(outMichels));
  // e.put(std::move(outAssns));
  e.put(std::move(outMichelHits));
}

float ana::MichelProd::GetDistance(PtrHit const& ph1, PtrHit const& ph2) {
  ana::Sec_t sec1 = ana::tpc2sec.at(geoDet).at(ph1->WireID().TPC);
  ana::Sec_t sec2 = ana::tpc2sec.at(geoDet).at(ph2->WireID().TPC);
  if (sec1 != sec2) return std::numeric_limits<float>::max();
  float z1 = asWire->Wire(ph1->WireID()).GetCenter().Z();
  float z2 = asWire->Wire(ph2->WireID()).GetCenter().Z();
  float t1 = ph1->PeakTime() * fTick2cm;
  float t2 = ph2->PeakTime() * fTick2cm;
  return sqrt(pow(z1-z2, 2) + pow(t1-t2, 2));
}

DEFINE_ART_MODULE(ana::MichelProd)
