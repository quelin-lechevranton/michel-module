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
SNElectron class properties are prefixed:
 - `as`  for art services
 - `geo` for geometric information
 - `in`  for input variable, taken from fhicl file

 - `ev`  for reco info about the full event (stored in `evTree`)
 - `mu`  for reco info about the track and michel electron (stored in `muTree`)
 - `tru` for truth info about the track (stored in `muTree`)
 - `mi`  for truth info about the michel electron (stored in `muTree`)

Top means X>0, Bot means X<0 (by reference to PDVD)
*/

namespace ana { class SNElectron; }

class ana::SNElectron: 
  public art::EDAnalyzer, 
  private ana::MichelModule 
{
public:
  explicit SNElectron(fhicl::ParameterSet const& p);
  SNElectron(SNElectron const&) = delete;
  SNElectron(SNElectron&&) = delete;
  SNElectron& operator=(SNElectron const&) = delete;
  SNElectron& operator=(SNElectron&&) = delete;

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

  TTree       *_evTree;
  TTree       *_partTree;

  unsigned                _evIndex=0;
  unsigned                _evPartNumber;
  std::vector<unsigned>   _evPartIndices;
  ana::Hits               _evHits;
  unsigned                _partIndex=0;

  int                     _pdg;
  float                   _energy;
  ana::Point              _startPoint;
  ana::Point              _endPoint;
};

ana::SNElectron::SNElectron(fhicl::ParameterSet const& p)
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
    case kFDHD:
      geoBot = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 0}).Min(),
        asGeo->TPC(geo::TPCID{0, 22}).Max()
      };
      geoTop = ana::Bounds3D<float>{
        asGeo->TPC(geo::TPCID{0, 1}).Min(),
        asGeo->TPC(geo::TPCID{0, 23}).Max()
      };
      break;
    default:
      std::cout << "SNElecModule: " << "\033[1;91m" << "geometry not handled (" << ana::det_name[geoDet] << ")" << "\033[0m" << std::endl;
  }
  geoCathodeGap = geoTop.x.min - geoBot.x.max;

  std::cout << "SNElecModule: " "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
    << "  Detector Geometry: " << asGeo->DetectorName()
    << "  (" << ana::det_name[geoDet] << ")" << std::endl
    << "  Tick Window: " << wireWindow << std::endl
    << "  Top Bounds: " << geoTop << std::endl
    << "  Bot Bounds: " << geoBot << std::endl
  ;

  _evTree = asFile->make<TTree>("event", "");

  _evTree->Branch("Index",       &_evIndex);
  _evTree->Branch("PartNumber",  &_evPartNumber);
  _evTree->Branch("PartIndices", &_evPartIndices);
  SetBranches(_evTree, "",       &_evHits);

  _partTree = asFile->make<TTree>("particle","");

  _partTree->Branch("EventIndex",  &_evIndex);
  _partTree->Branch("IndexInEvent",&_evPartNumber);
  _partTree->Branch("Index",       &_partIndex);
  _partTree->Branch("pdg",         &_pdg);
  _partTree->Branch("energy",      &_energy);
  SetBranches(_partTree, "start",  &_startPoint);
  SetBranches(_partTree, "end",    &_endPoint);
}

void ana::SNElectron::analyze(art::Event const& e) {
  auto const clockData = asDetClocks->DataFor(e);
  auto const detProp = asDetProp->DataFor(e,clockData);
  fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

  auto const& vh_mct = e.getHandle<std::vector<simb::MCTruth>>(tag_mct);
  if (!vh_mct.isValid()) { std::cout<<"SNElecModule: No valid simb::MCTruth handle"<<std::endl; return; }

  auto const& vh_mcp = e.getHandle<std::vector<simb::MCParticle>>(tag_mcp);
  if (!vh_mcp.isValid()) { std::cout<<"SNElecModule: No valid simb::MCParticle handle"<<std::endl; return; }

  auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
  if (!vh_hit.isValid()) { std::cout<<"SNElecModule: No valid recob::Hit handle"<<std::endl; return; }
  VecPtrHit vph_ev;
  art::fill_ptr_vector(vph_ev, vh_hit);

  auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
  if (!vh_trk.isValid()) { std::cout<<"SNElecModule: No valid recob::Track handle"<<std::endl; return; }
  VecPtrTrk vpt_ev;
  art::fill_ptr_vector(vpt_ev, vh_trk);

  auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
  if (!vh_shw.isValid()) { std::cout<<"SNElecModule: No valid recob::Shower handle"<<std::endl; return; }
  VecPtrShw vps_ev;
  art::fill_ptr_vector(vps_ev, vh_shw);

  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
  art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
  art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
  art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);


  std::cout << "mct (" << vh_mct->size() << "): {" << std::endl;
  for (auto const& t : *vh_mct) {
    std::cout << "\t";
    for (int i=0; i<t.NParticles(); i++) 
      std::cout << t.GetParticle(i) << ", ";
    std::cout << std::endl;
  }
  std::cout << "}" << std::endl;

  std::cout << "mcp (" << vh_mcp->size() << "): { ";
  for (auto const& p : *vh_mcp) std::cout << p.PdgCode() << ", ";
  std::cout << "}" << std::endl;


  // dump event information
  // evRun = e.run();
  // evSubRun = e.subRun();
  // evEvent = e.event();
  // evIsData = e.isRealData();

  _evPartNumber=0;
  _evPartIndices.clear();
  _evHits.clear();
  for (PtrHit p_hit : vph_ev)
    if (p_hit->View() == geo::kW)
      _evHits.push_back(GetHit(p_hit));

  // std::cout << vh_mcp->size() << " particles" << std::endl;
  if (vh_mcp->empty()) return;
  // int np=0;
  for (simb::MCParticle const& part : *vh_mcp) {

    // std::cout << "#" << ++np << ": " << part.PdgCode() << "   ";

    _pdg = part.PdgCode();
    _energy = (part.E() - part.Mass()) * 1e3;
    _startPoint = part.Position().Vect();
    _endPoint = part.EndPosition().Vect();

    _partTree->Fill();
    _evPartIndices.push_back(_partIndex);
    _partIndex++;
    _evPartNumber++;
  }
  // std::cout << std::endl;

  _evTree->Fill();
  _evIndex++;
}

void ana::SNElectron::beginJob() {}
void ana::SNElectron::endJob() {}

DEFINE_ART_MODULE(ana::SNElectron)