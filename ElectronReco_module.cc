#include "utils.h"

/* 
Naming conventions,
some variables are prefixed by their type followed by an underscore:
 - `ph_`  for `art::Ptr<recob::Hit>`
 - `vph_` for `std::vector<art::Ptr<recob::Hit>>`

 - `pt_`  for `art::Ptr<recob::Track>`
 - `vpt_` for `std::vector<art::Ptr<recob::Track>>`

 - `ps_`  for `art::Ptr<recob::Shower>`
 - `vps_` for `std::vector<art::Ptr<recob::Shower>>`

 - `tag_` for `art::InputTag`
 - `mcp_` for `simb::MCParticle*`

 - `sh_` for custom `ana::SortedHits`, a structure containing info on the sorted list of hits of a track
ElectronReco class properties are prefixed:
 - `as`  for art services
 - `geo` for geometric information
 - `in`  for input variable, taken from fhicl file

 - `ev`  for reco info about the full event (stored in `evTree`)
 - `mu`  for reco info about the track and michel electron (stored in `muTree`)
 - `tru` for truth info about the track (stored in `muTree`)
 - `mi`  for truth info about the michel electron (stored in `muTree`)

Top means X>0, Bot means X<0 (by reference to PDVD)
*/

namespace ana { class ElectronReco; }

class ana::ElectronReco: 
  public art::EDAnalyzer, 
  private ana::MichelModule 
{
public:
  explicit ElectronReco(fhicl::ParameterSet const& p);
  ElectronReco(ElectronReco const&) = delete;
  ElectronReco(ElectronReco&&) = delete;
  ElectronReco& operator=(ElectronReco const&) = delete;
  ElectronReco& operator=(ElectronReco&&) = delete;

  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
private:
  // Geometry Information
  ana::Bounds<float> wireWindow;
  ana::Bounds3D<float> geoTop, geoBot;
  float geoCathodeGap; // cm

  // Input Parameters
  bool        inLog;

  TTree       *_tree;

  // Output Variables
  float _muEnergy;
  ana::Point _muDir;
  float _muTheta;
  float _muPhi;
  float _muEndEnergy;
  ana::Point _muEndDir;
  float _muEndTheta;
  float _muEndPhi;
  float _miEnergy;
  ana::Point _miDir;
  ana::Hit _muEndHit;
  ana::Point _muEndPoint;

  bool _muHasTrack;
  enum EnumMichelHasTrack: int { 
    kIsMuonTrack = -1, kNoTrack = 0, kHasTrack = 1
  };
  EnumMichelHasTrack _miHasTrack;
  bool _miHasShower;
  float _miTrackLength;
  float _miShowerLength;

  ana::Hits _miHits;
  std::vector<float> _miHitEnergyFrac;
  std::vector<float> _miHitDistance;
  std::vector<float> _miHitAngle;
  std::vector<bool> _miHitIsPrime;
  std::vector<bool> _miHitFromMichelTrack;
  std::vector<bool> _miHitFromMichelShower;
  std::vector<bool> _miHitFromMuonTrack;
};

ana::ElectronReco::ElectronReco(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} 
  , MichelModule{p}
  , inLog(p.get<bool>("Log", true))
{
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
    case kPDSP:
      geoBot = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 1}).Min(),
        asGeo->TPC(geo::TPCID{0, 9}).Max()
      };
      geoTop = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 2}).Min(),
        asGeo->TPC(geo::TPCID{0, 10}).Max()
      };
      break;
    default:
      std::cout << "ElectronRecoModule: " << "\033[1;93m" << "geometry not handled: " << ana::det_name[geoDet] << "\033[0m" << std::endl;
  }
  geoCathodeGap = geoTop.x.min - geoBot.x.max;

  std::cout << "ElectronRecoModule: " "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
    << "  Detector Geometry: " << asGeo->DetectorName()
    << "  (" << ana::det_name[geoDet] << ")" << std::endl
    << "  Tick Window: " << wireWindow << std::endl
    << "  Top Bounds: " << geoTop << std::endl
    << "  Bot Bounds: " << geoBot << std::endl
  ;
  std::cout << "ElectronRecoModule: " "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
  ;

  _tree = asFile->make<TTree>("michel","");

  _tree->Branch("muEnergy",              &_muEnergy);
  SetBranches(_tree, "muDir",            &_muDir);
  _tree->Branch("muTheta",               &_muTheta);
  _tree->Branch("muPhi",                 &_muPhi);
  _tree->Branch("muEndEnergy",           &_muEndEnergy);
  SetBranches(_tree, "muEndDir",         &_muEndDir);
  _tree->Branch("muEndTheta",            &_muEndTheta);
  _tree->Branch("muEndPhi",              &_muEndPhi);
  _tree->Branch("miEnergy",              &_miEnergy);
  SetBranches(_tree, "miDir",            &_miDir);
  SetBranches(_tree, "muEnd",            &_muEndHit);
  SetBranches(_tree, "muEnd",            &_muEndPoint);

  _tree->Branch("muHasTrack",            &_muHasTrack);
  _tree->Branch("miHasTrack",      (int*)&_miHasTrack);
  _tree->Branch("miHasShower",           &_miHasShower);
  _tree->Branch("miTrackLength",         &_miTrackLength);
  _tree->Branch("miShowerLength",        &_miShowerLength);

  SetBranches(_tree, "mi",               &_miHits);
  _tree->Branch("miHitEnergyFrac",       &_miHitEnergyFrac);
  _tree->Branch("miHitDistance",         &_miHitDistance);
  _tree->Branch("miHitAngle",            &_miHitAngle);
  _tree->Branch("miHitIsPrime",          &_miHitIsPrime);
  _tree->Branch("miHitFromMichelTrack",  &_miHitFromMichelTrack);
  _tree->Branch("miHitFromMichelShower", &_miHitFromMichelShower);
  _tree->Branch("miHitFromMuonTrack",    &_miHitFromMuonTrack);
}

void ana::ElectronReco::analyze(art::Event const& e) {
  auto const clockData = asDetClocks->DataFor(e);
  auto const detProp = asDetProp->DataFor(e,clockData);
  fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

  auto const& vh_mcp = e.getHandle<std::vector<simb::MCParticle>>(tag_mcp);
  if (!vh_mcp.isValid()) { std::cout<<"ElectronRecoModule: No valid simb::MCParticle handle"<<std::endl; return; }

  auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
  if (!vh_hit.isValid()) { std::cout<<"ElectronRecoModule: No valid recob::Hit handle"<<std::endl; return; }
  VecPtrHit vph_ev;
  art::fill_ptr_vector(vph_ev, vh_hit);

  auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
  if (!vh_trk.isValid()) { std::cout<<"ElectronRecoModule: No valid recob::Track handle"<<std::endl; return; }
  VecPtrTrk vpt_ev;
  art::fill_ptr_vector(vpt_ev, vh_trk);

  auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
  if (!vh_shw.isValid()) { std::cout<<"ElectronRecoModule: No valid recob::Shower handle"<<std::endl; return; }
  VecPtrShw vps_ev;
  art::fill_ptr_vector(vps_ev, vh_shw);

  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
  art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
  art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
  art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

  // dump event information
  // evRun = e.run();
  // evSubRun = e.subRun();
  // evEvent = e.event();
  // evIsData = e.isRealData();
  // for (PtrHit p_hit : vph_ev)
  //     if (p_hit->View() == geo::kW)
  //         evHits.push_back(GetHit(p_hit));

  for (simb::MCParticle const& muon : *vh_mcp) {
    simb::MCParticle const* michel = GetMichelMCP(&muon);
    if (!michel) continue;

    VecPtrHit muonHits = ana::mcp2hits(&muon, vph_ev, clockData, false);
    if (muonHits.empty()) continue;

    _muEnergy = (muon.E() - muon.Mass()) * 1e3;
    _muEndEnergy = (muon.EndE() - muon.Mass()) * 1e3;
    _miEnergy = (michel->E() - michel->Mass()) * 1e3;

    TVector3 const& muVect = muon.Momentum().Vect();
    geo::Vector_t muonDir;
    switch (geoDet) {
      case kPDVD:
        muonDir.SetCoordinates(muVect.Y(), muVect.Z(), muVect.X());
        break;
      case kPDHD:
      case kPDSP:
        muonDir.SetCoordinates(muVect.Z(), muVect.X(), muVect.Y());
        break;
      default:
    }
    _muDir = ana::Point(muVect.Unit());
    _muTheta = muonDir.Theta();
    _muPhi = muonDir.Phi();

    _muEndPoint = ana::Point(muon.EndPosition().Vect());
    TVector3 const& muEndVect = muon.EndMomentum().Vect();
    geo::Vector_t muonEndDir;
    switch (geoDet) {
      case kPDVD:
        muonEndDir.SetCoordinates(muEndVect.Y(), muEndVect.Z(), muEndVect.X());
        break;
      case kPDHD:
      case kPDSP:
        muonEndDir.SetCoordinates(muEndVect.Z(), muEndVect.X(), muEndVect.Y());
        break;
      default:
    }
    _muEndDir = ana::Point(muEndVect.Unit());
    _muEndTheta = muonEndDir.Theta();
    _muEndPhi = muonEndDir.Phi();

    PtrTrk const& muonTrack = ana::mcp2trk(&muon, vpt_ev, clockData, fmp_trk2hit);
    _muHasTrack = muonTrack.isNonnull();
    // geo::Point_t endPoint = geo::Point_t(muon.EndPosition().Vect());

    PtrTrk const& michelTrack = ana::mcp2trk(michel, vpt_ev, clockData, fmp_trk2hit);

    _miHasTrack = michelTrack.isNull() ? kNoTrack : (_muHasTrack && muonTrack.key()==michelTrack.key() ? kIsMuonTrack : kHasTrack);
    _miTrackLength = _miHasTrack ? michelTrack->Length() : util::kBogusF;

    PtrShw const& michelShower = ana::mcp2shw(michel, vps_ev, clockData, fmp_shw2hit);
    _miHasShower = michelShower.isNonnull();
    _miShowerLength = _miHasShower ? michelShower->Length() : util::kBogusF;

    _miDir = ana::Point(michel->Momentum().Vect().Unit());
    ana::Vec2 michelVec2(_miDir.z, _miDir.x);
    float michelVec2Angle = michelVec2.angle();

    _miHitEnergyFrac.clear();
    VecPtrHit eveHits = ana::mcp2hits(michel, vph_ev, clockData, true, &_miHitEnergyFrac);
    VecPtrHit primeHits = ana::mcp2hits(michel, vph_ev, clockData, false);

    int increasingX = muVect.X() > 0 ? +1 : -1;
    PtrHit const& endHit = *std::max_element(
      muonHits.begin(), muonHits.end(), 
      [&](PtrHit const& hit1, PtrHit const& hit2) {
        if (GetSide(hit1) == kTop && GetSide(hit2) == kTop)
          return increasingX * (hit1->PeakTime() - hit2->PeakTime()) > 0;
        if (GetSide(hit1) == kBot && GetSide(hit2) == kBot)
          return increasingX * (hit1->PeakTime() - hit2->PeakTime()) < 0;
        return GetSide(hit1) == (increasingX>0 ? kBot : kTop);
      }
    );
    _muEndHit = GetHit(endHit);

    _miHits.clear();
    _miHits.reserve(eveHits.size());
    _miHitDistance.clear();
    _miHitDistance.reserve(eveHits.size());
    _miHitAngle.clear();
    _miHitAngle.reserve(eveHits.size());
    _miHitIsPrime.clear();
    _miHitIsPrime.reserve(eveHits.size());
    _miHitFromMichelTrack.clear();
    _miHitFromMichelTrack.reserve(eveHits.size());
    _miHitFromMichelShower.clear();
    _miHitFromMichelShower.reserve(eveHits.size());
    _miHitFromMuonTrack.clear();
    _miHitFromMuonTrack.reserve(eveHits.size());

    for (PtrHit const& hitPtr : eveHits) {
      ana::Hit hit = GetHit(hitPtr);
      _miHits.push_back(hit);

      _miHitIsPrime.push_back(std::find_if(
        primeHits.begin(), primeHits.end(), 
        [key=hitPtr.key()](PtrHit const& hit1){ return hit1.key()==key; }
      ) != primeHits.end());
      _miHitDistance.push_back(GetDistance(hit, _muEndHit));

      ana::Vec2 decay2hit(hit.space - _muEndHit.space, (hit.tick - _muEndHit.tick) * fTick2cm);
      float angle = decay2hit.angle() - michelVec2Angle;
      angle = abs(angle) > M_PI ? angle - (angle>0 ? 1 : -1) * 2*M_PI : angle;
      _miHitAngle.push_back(angle);

      PtrTrk const& hitTrack = fop_hit2trk.at(hitPtr.key());
      _miHitFromMichelTrack.push_back(michelTrack.isNonnull() && hitTrack.isNonnull() && hitTrack.key() == michelTrack.key());
      _miHitFromMuonTrack.push_back(muonTrack.isNonnull() && hitTrack.isNonnull() && hitTrack.key() == muonTrack.key());

      PtrShw const& hitShower = fop_hit2shw.at(hitPtr.key());
      _miHitFromMichelShower.push_back(michelShower.isNonnull() && hitShower.isNonnull() && hitShower.key() == michelShower.key());
    }

    _tree->Fill();
  }
}

void ana::ElectronReco::beginJob() {}
void ana::ElectronReco::endJob() {}

DEFINE_ART_MODULE(ana::ElectronReco)