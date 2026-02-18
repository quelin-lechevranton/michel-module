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

  bool      inLog;
  float     inTrackLengthCut; // in cm
  float     inFiducialLength; // in cm
  unsigned  inSmoothingLength;
  float     inBraggCut; // in fC/cm
  unsigned  inMinBaryHits;
  float     inBarycenterRadius; // in cm
  float     inMichelRadius; // in cm
  float     inAngleCut; // in deg
};

void ana::MichelProd::beginJob() {
  // tuple = asFile->make<TTree>("calo", "");
  // tuple->Branch("ResidualRange", &residual_range);
  // tuple->Branch("dEds", &dEds);
}

ana::MichelProd::MichelProd(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , MichelModule{p}
  , inLog             (p.get<bool>    ("Log", true))
  , inTrackLengthCut  (p.get<float>   ("TrackLengthCut", 30.F)) // in cm
  , inFiducialLength  (p.get<float>   ("FiducialLength", 20.F)) // in cm
  , inSmoothingLength (p.get<unsigned>("SmoothingLength", 6))
  , inBraggCut        (p.get<float>   ("BraggCut", 13.5F)) // in fC/cm
  , inMinBaryHits     (p.get<unsigned>("MinBaryHits", 4))
  , inBarycenterRadius(p.get<float>   ("BarycenterRadius", 10.F)) // in cm
  , inMichelRadius    (p.get<float>   ("MichelRadius", 20.F)) //in cm
  , inAngleCut        (p.get<float>   ("AngleCut", 30.F)) //in cm
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

  std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
    << "  Detector Geometry: " << asGeo->DetectorName()
    << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
    << "  Tick Window: " << wireWindow << std::endl
    << "  Top Bounds: " << geoTop << std::endl
    << "  Bot Bounds: " << geoBot << std::endl;
  std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
    << "  Track Length Cut: "     << inTrackLengthCut << " cm" << std::endl
    << "  Fiducial Length: "      << inFiducialLength << " cm" << std::endl
    << "  Smoothing Length: "     << inSmoothingLength << std::endl
    << "  Bragg Cut: "            << inBraggCut << " fC/cm" << std::endl
    << "  Min Bary Hits: "        << inMinBaryHits << std::endl
    << "  Barycenter Radius: "    << inBarycenterRadius << " cm" << std::endl
    << "  Bary Center Radius: "   << inBarycenterRadius << " cm" << std::endl
    << "  Michel Radius: "        << inMichelRadius << " cm" << std::endl
    << "  Michel Space Radius: "  << inMichelRadius << " cm" << std::endl
    << "  Angle Cut: "            << inAngleCut << " deg" << std::endl;
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
  auto const& vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
  art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
  if (!vh_hit.isValid()
   || !vh_trk.isValid()
   || !fmp_trk2hit.isValid()
   || !fop_hit2trk.isValid()
  ) {
    if (inLog) std::cout << "\033[91;1m" "missing required product" "\033[0m" << std::endl;
    e.put(std::move(outMichelHits));
    e.put(std::move(outAssns));
    return;
  }

  VecPtrHit vph_ev; art::fill_ptr_vector(vph_ev, vh_hit);
  VecPtrTrk vpt_ev; art::fill_ptr_vector(vpt_ev, vh_trk);

  size_t trk_idx=0;
  for (PtrTrk const& pt_ev: vpt_ev) {
    if (inLog) std::cout << "t" << ++trk_idx << "/" << vpt_ev.size() << "\r" << std::flush;
    ASSERT(pt_ev->Length() > inTrackLengthCut)

    VecPtrHit vph_trk = fmp_trk2hit.at(pt_ev.key());
    ASSERT(!vph_trk.empty())
    std::vector<recob::TrackHitMeta const*> const& vhm_trk = fmp_trk2hit.data(pt_ev.key());
    ASSERT(vph_trk.size() == vhm_trk.size())
    std::map<size_t, unsigned> map_hitkey2idx;

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
    bool end_in_y = geoBot.y.isInside(end_y, inFiducialLength) || geoTop.y.isInside(end_y, inFiducialLength);
    bool end_in_z = geoBot.z.isInside(end_hit.space, inFiducialLength) || geoTop.z.isInside(end_hit.space, inFiducialLength);
    bool end_in_t = wireWindow.isInside(end_hit.tick, inFiducialLength / fTick2cm);

    Side_t end_side = ana::tpc2side.at(geoDet).at(end_hit.tpc);
    float end_x = util::kBogusF;
    bool end_in_x = false;
    if (sh_mu.is_cc()) {
      end_x = end_side == kBot
        ? -(geoCathodeGap/2) - (sh_mu.cc_second()->PeakTime() - end_hit.tick) * fTick2cm
        : +(geoCathodeGap/2) + (sh_mu.cc_first()->PeakTime() - end_hit.tick) * fTick2cm;
      end_in_x = (end_side==kBot ? geoBot : geoTop).x.isInside(end_x, inFiducialLength);
    }
    ASSERT(end_in_x && end_in_y && end_in_z && end_in_t)

    bool track_is_up =  IsUpright(*pt_ev);
    geo::Point_t Start = track_is_up ? pt_ev->Start() : pt_ev->End();
    geo::Point_t End = track_is_up ? pt_ev->End() : pt_ev->Start();
    float cosY = 1 / sqrt(1 + pow(Start.Y()-End.Y(), 2) / (pow(Start.X()-End.X(), 2) + pow(Start.Z()-End.Z(), 2)));

    size_t queue = 0;
    float dist;
    do {
      PtrHit const& ph = sh_mu.vph.at(sh_mu.vph.size() - ++queue);
      dist = GetDistance(ph, sh_mu.end());
    } while (queue < sh_mu.vph.size() && dist < 2);

    std::vector<float> dQds = GetdQds(sh_mu.after_cathode_it(), sh_mu.vph.end(), inSmoothingLength);
    float ADC2fC = 1.60e-4 * 200; // el2fC * ADC2el
    float end_dQds = cosY * ADC2fC * std::accumulate(
      dQds.end()-queue, dQds.end(), 0.F
    ) / queue;
    if (inLog) std::cout << "\tqueue: " << queue << std::endl;
    ASSERT(end_dQds > inBraggCut)

    VecPtrHit vph_ev_endsec;
    for (PtrHit const& ph_ev : vph_ev) {
      Sec_t sec = ana::tpc2sec.at(geoDet).at(ph_ev->WireID().TPC);
      if (sec != sh_mu.end_sec()) continue;
      vph_ev_endsec.push_back(ph_ev);
    }

    ana::Hits bary_hits;
    for (PtrHit const& ph_ev : vph_ev_endsec) {
      if (ph_ev->View() != geo::kW) continue;
      if (GetDistance(ph_ev, sh_mu.end()) > inBarycenterRadius) continue;

      PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
      if (pt_hit && pt_hit->Length() > inTrackLengthCut) continue;

      bary_hits.push_back(GetHit(ph_ev));
    }

    ASSERT(bary_hits.size() >= inMinBaryHits)
    ana::Vec2 end_to_bary = bary_hits.barycenter(fTick2cm) - end_hit.vec(fTick2cm);
    float mu_angle = sh_mu.regs
      .at(ana::sec2side.at(geoDet).at(end_hit.section))
      .theta(end_hit.space > start_hit.space ? 1 : -1);
    float da = end_to_bary.angle() - mu_angle;
    da = abs(da) > M_PI ? da - (da>0 ? 1 : -1)*2*M_PI : da;
    da = 90 - abs(abs(da * TMath::RadToDeg()) - 90);
    ASSERT(da > inAngleCut)

    if (inLog) std::cout << "\t\033[93;1m" "michel" "\033[0m" << std::endl;

    for (PtrHit const& ph_ev: vph_ev_endsec) {
      if (GetDistance(ph_ev, end_side, end_y, end_hit.space, end_hit.tick) > inMichelRadius) continue;

      PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
      if (pt_hit && pt_hit->Length() > inTrackLengthCut) continue;

      outMichelHits->emplace_back(*ph_ev);
      art::Ptr<recob::Hit> ph_out = hitPtrMaker(outMichelHits->size()-1);
      outAssns->addSingle(pt_ev, ph_out);
    }
  }

  if (inLog) std::cout << "\033[93;1m" << outMichelHits->size() << " michel(s) found" "\033[0m" << std::endl;
  
  // Put outputs in event
  e.put(std::move(outMichelHits));
  e.put(std::move(outAssns));
}

DEFINE_ART_MODULE(ana::MichelProd)
