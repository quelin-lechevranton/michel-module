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

class ana::MichelProd : public art::EDProducer, private ana::MichelModule {
public:
  explicit MichelProd(fhicl::ParameterSet const& p);
  MichelProd(MichelProd const&) = delete;
  MichelProd(MichelProd&&) = delete;
  MichelProd& operator=(MichelProd const&) = delete;
  MichelProd& operator=(MichelProd&&) = delete;

  void produce(art::Event& e) override;
  void beginJob() override;
private:
  ana::Bounds<float> wireWindow;
  ana::Bounds3D<float> geoTop, geoBot;
  float geoCathodeGap;

  bool inLog;
  TTree *tuple;
  std::vector<float> residual_range, dEds;
};

void ana::MichelProd::beginJob() {
  // tuple = asFile->make<TTree>("calo", "");
  // tuple->Branch("ResidualRange", &residual_range);
  // tuple->Branch("dEds", &dEds);
}

ana::MichelProd::MichelProd(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , MichelModule{p}
  , inLog(p.get<bool>("Log", true))
{
  produces<std::vector<recob::Hit>>();
  produces<art::Assns<recob::Track, recob::Hit>>();

  auto const clockData = asDetClocks->DataForJob();
  auto const detProp = asDetProp->DataForJob(clockData);
  fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
  wireWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
  switch (geoDet) {
    case kPDVD:
      geoBot = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 0}).Min(),
        asGeo->TPC(geo::TPCID{0, 7}).Max()
      };
      geoTop = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 8}).Min(),
        asGeo->TPC(geo::TPCID{0, 15}).Max()
      };
      break;
    case kPDHD:
      geoBot = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 1}).Min(),
        asGeo->TPC(geo::TPCID{0, 5}).Max()
      };
      geoTop = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 2}).Min(),
        asGeo->TPC(geo::TPCID{0, 6}).Max()
      };
      break;
    default: break;
  }
  geoCathodeGap = geoTop.x.min - geoBot.x.max;
}

void ana::MichelProd::produce(art::Event& e)
{
  // Output
  auto outMichelHits = std::make_unique<std::vector<recob::Hit>>();
  auto outAssns = std::make_unique<art::Assns<recob::Track, recob::Hit>>();
  art::PtrMaker<recob::Hit> hitPtrMaker(e);

  auto const clockData = asDetClocks->DataFor(e);
  auto const detProp = asDetProp->DataFor(e,clockData);
  fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

  // Input
  auto const& vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
  if (!vh_hit.isValid()) return;
  VecPtrHit vph_ev; art::fill_ptr_vector(vph_ev, vh_hit);

  auto const& vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
  if (!vh_trk.isValid()) return;
  VecPtrTrk vpt_ev; art::fill_ptr_vector(vpt_ev, vh_trk);

  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
  art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

  for (PtrTrk const& pt_ev: vpt_ev) {
    ASSERT(pt_ev->Length() < 30);

    VecPtrHit vph_trk = fmp_trk2hit.at(pt_ev.key());
    std::vector<recob::TrackHitMeta const*> const& vhm_trk = fmp_trk2hit.data(pt_ev.key());
    std::map<size_t, unsigned> map_hitkey2idx;
    ASSERT(!vph_trk.empty())
    ASSERT(vph_trk.size() == vhm_trk.size())

    // remove hits not associated to a track point (hit misassociated to the track?)
    std::vector<unsigned> bad_hit_indices;
    for (unsigned i=0; i<vph_trk.size(); i++) {
      if (vph_trk[i]->View() != geo::kW) continue;
      if (!pt_ev->HasValidPoint(vhm_trk[i]->Index())) {
        bad_hit_indices.push_back(i);
      } else {
        map_hitkey2idx[vph_trk[i].key()] = i;
      }
    }
    for (int i=bad_hit_indices.size()-1; i>=0; i--)
      vph_trk.erase(vph_trk.begin() + bad_hit_indices[i]);

    ana::SortedHits sh_mu = GetSortedHits(vph_trk);
    ASSERT(bool(sh_mu))
    ASSERT(sh_mu.is_cc())

    // bool track_is_up =  IsUpright(*pt_ev);
    // geo::Point_t Start = track_is_up ? pt_ev->Start() : pt_ev->End();
    // geo::Point_t End = track_is_up ? pt_ev->End() : pt_ev->Start();

    bool is_up = true;
    if (geoDet == kPDVD) {
      Side_t front_side = ana::tpc2side.at(geoDet).at(sh_mu.start()->WireID().TPC);
      is_up = sh_mu.is_cc()
        ? front_side == kTop
        : ( front_side == kTop
          ? sh_mu.start()->PeakTime() < sh_mu.end()->PeakTime()
          : sh_mu.start()->PeakTime() > sh_mu.end()->PeakTime()
        );
    } else if (geoDet == kPDHD) {
      size_t front_hit_track_idx = vhm_trk[map_hitkey2idx.at(sh_mu.start().key())]->Index();
      float front_hit_y = pt_ev->HasValidPoint(front_hit_track_idx)
        ? pt_ev->LocationAtPoint(front_hit_track_idx).Y()
        : util::kBogusF;
      size_t back_hit_track_idx = vhm_trk[map_hitkey2idx.at(sh_mu.end().key())]->Index();
      float back_hit_y = pt_ev->HasValidPoint(back_hit_track_idx)
        ? pt_ev->LocationAtPoint(back_hit_track_idx).Y()
        : util::kBogusF;
      is_up = front_hit_y > back_hit_y;
    }
    if (!is_up) sh_mu.reverse();

    ana::Hit start_hit = GetHit(sh_mu.start());
    ana::Hit end_hit = GetHit(sh_mu.end());
    size_t end_track_idx = vhm_trk[map_hitkey2idx.at(sh_mu.end().key())]->Index();
    float end_y = pt_ev->HasValidPoint(end_track_idx)
      ? pt_ev->LocationAtPoint(end_track_idx).Y()
      : util::kBogusF;
    bool end_in_y = geoBot.y.isInside(end_y, 20) || geoTop.y.isInside(end_y, 20);
    bool end_in_z = geoBot.z.isInside(end_hit.space, 20) || geoTop.z.isInside(end_hit.space, 20);
    bool end_in_t = wireWindow.isInside(end_hit.tick, 20 / fTick2cm);

    Side_t end_side = ana::tpc2side.at(geoDet).at(end_hit.tpc);
    float end_x = !sh_mu.is_cc() ? util::kBogusF : end_side == kBot
      ? -(geoCathodeGap/2) - (sh_mu.cc_second()->PeakTime() - end_hit.tick) * fTick2cm
      : +(geoCathodeGap/2) + (sh_mu.cc_first()->PeakTime() - end_hit.tick) * fTick2cm;
    bool end_in_x = !sh_mu.is_cc() ? false 
      : (end_side==kBot ? geoBot : geoTop).x.isInside(end_x, 20);
    // LOG(end_in_x && end_in_y && end_in_z && end_in_t);
    ASSERT(end_in_x && end_in_y && end_in_z && end_in_t)


    std::vector<float> dQds = GetdQds(sh_mu.after_cathode_it(), sh_mu.vph.end(), 6);
    float ADC2fC = 1.60e-4 * 200; // el2fC * ADC2el
    // LOG(dQds.back() > 13.5 / ADC2fC);
    ASSERT(dQds.back() > 13.5 / ADC2fC)


    VecPtrHit vph_ev_endsec;
    for (PtrHit const& ph_ev : vph_ev) {
      if (ph_ev->View() != geo::kW) continue;
      Sec_t sec = ana::tpc2sec.at(geoDet).at(ph_ev->WireID().TPC);
      if (sec != sh_mu.end_sec()) continue;
      vph_ev_endsec.push_back(ph_ev);
    }

    ana::Hits bary_hits;
    for (PtrHit const& ph_ev : vph_ev_endsec) {
      if (GetDistance(ph_ev, sh_mu.end()) > 10) continue;

      PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
      if (pt_hit && pt_hit->Length() > 30) continue;

      bary_hits.push_back(GetHit(ph_ev));
    }

    ASSERT(bary_hits.size() > 4)
    ana::Vec2 end_to_bary = bary_hits.barycenter(fTick2cm) - end_hit.vec(fTick2cm);
    float mu_angle = sh_mu.regs
      .at(ana::sec2side.at(geoDet).at(end_hit.section))
      .theta(end_hit.space > start_hit.space ? 1 : -1);
    float da = end_to_bary.angle() - mu_angle;
    da = abs(da) > M_PI ? da - (da>0 ? 1 : -1)*2*M_PI : da;
    da = 90 - abs(abs(da * TMath::RadToDeg()) - 90);
    ASSERT(da > 30)

    if (inLog) std::cout << "\033[93;1m" "michel" "\033[0m" << std::endl;

    for (PtrHit const& ph_ev: vph_ev_endsec) {
      if (GetDistance(ph_ev, sh_mu.end()) > 30) continue;

      PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
      if (pt_hit && pt_hit->Length() > 30) continue;

      outMichelHits->emplace_back(*ph_ev);
      art::Ptr<recob::Hit> ph_out = hitPtrMaker(outMichelHits->size()-1);
      outAssns->addSingle(pt_ev, ph_out);
    }
  }
  
  // Put outputs in event
  e.put(std::move(outMichelHits));
  e.put(std::move(outAssns));
}

DEFINE_ART_MODULE(ana::MichelProd)
