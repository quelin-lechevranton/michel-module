////////////////////////////////////////////////////////////////////////
// Class:       Muchecks
// Plugin Type: analyzer
// File:        Muchecks_module.cc
//
// Generated at Wed Dec 18 10:31:28 2024 by Jeremy Quelin Lechevranton
////////////////////////////////////////////////////////////////////////

// Art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "larcorealg/Geometry/Exceptions.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// #include "lardata/ArtDataHelper/TrackUtils.h"

// ProtoDUNE includes
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

// ROOT includes
#include "TLorentzVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEllipse.h"


// #include "ROOT/RVec.hxx"
// #include "RVec.hxx"

// std includes
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <map>
#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <utility>

// using std::std::cout;
// using std::std::endl;
// using std::std::flush;
// using std::vector;
// using std::string;
// using std::unordered_map;
// using std::pair;
// using std::min;
// using std::max;


namespace ana {
  class Muchecks;
  struct Binning {
        int n;
        double min, max;
  };
  struct Hit {
        int channel;
        float tick;
        float adc;
  };
}


class ana::Muchecks : public art::EDAnalyzer {
public:
    explicit Muchecks(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Muchecks(Muchecks const&) = delete;
    Muchecks(Muchecks&&) = delete;
    Muchecks& operator=(Muchecks const&) = delete;
    Muchecks& operator=(Muchecks&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

private:

    // Utilities
    art::ServiceHandle<art::TFileService> tfs;

    const geo::GeometryCore* asGeo;
    const geo::WireReadoutGeom* asWire;
    const detinfo::DetectorPropertiesService* asDetProp;
    const detinfo::DetectorClocksService* asDetClocks;

    protoana::ProtoDUNETruthUtils truthUtil;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Conversion factors
    float fADCtoE = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    float fChannelPitch;
    float fDriftVelocity; // cm/µs
    float fSamplingRate; // µs/tick

    // Verbosity
    int iLogLevel;
    enum EnumFlag { kImportant, kBasics, kInfos, kDetails };

    // Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag   tag_mcp,
                    tag_sed,
                    tag_wir,
                    tag_hit,
                    tag_clu,
                    tag_trk,
                    tag_spt;

    // Input Variables
    float fMichelTimeRadius, // in µs
          fMichelSpaceRadius, // in cm
          fTrackLengthCut; // in cm

    bool fKeepOutside,
         fKeepNonDecaying;
    
    // Output Variables

    std::vector<std::string> anomalies;

    enum EnumTrackInside { kOutside, kInUpper, kInLower };
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside };
    enum EnumTickInside { kNoTick = -1, kTickOutside = 0, kTickInside = 1 };

    TTree* tEvent;

    size_t iEvent=0;

    size_t EventNHit;
    std::vector<int> EventHitChannel; //[hit]
    std::vector<float> EventHitTick; //[hit]
    std::vector<float> EventHitADC; //[hit]

    size_t EventNMuon;
    std::vector<int> EventiMuon; //[muon]

    // enum Tag { kMCP, kDaughters, kSphere, kTrue, kNTag };
    TTree* tMuon;
    TTree* tMichel;

    size_t iMuon=0;

    int MuonIsAnti;
    int MuonDoesDecay;
    int MuonHasMichel;

    float MuonTrackLength;

    size_t MuonNTrackPoint;
    std::vector<float> MuonTrackPointX, MuonTrackPointY, MuonTrackPointZ;

    int MuonEndChannel;
    float MuonEndTick;
    int MuonEndIsInWindow;
    int MuonEndIsInTPC;
    int MuonEndIsInVolumeYZ;

    size_t MuonNHit;
    std::vector<ana::Hit> MuonHit; //[hit]

    size_t MicheliMuon;

    size_t MichelNHit;
    float MichelTrackLength;
    // int MichelChannelBary;
    // float MichelTickBary;
    std::vector<int> MichelHitChannel; //[hit]
    std::vector<float> MichelHitTick; //[hit]
    std::vector<float> MichelHitADC; //[hit]
    float MichelTrueEnergy;
    float MichelHitEnergy;

    size_t MichelSphereNHit;
    std::vector<int> MichelSphereChannel; //[hit]
    std::vector<float> MichelSphereTick; //[hit]
    std::vector<float> MichelSphereADC; //[hit]
    float MichelSphereEnergy;

    size_t MichelSphereTruePositive;
    size_t MichelSphereFalsePositive;

    // Diagnostic Variables
    unsigned int n_section = 4, n_plane = 3;
    // size_t iNEventDiagnostic;
    ana::Binning binTick;
    std::vector<ana::Binning> binChan; //[section]
    std::vector<ana::Binning> chanTPC; //[tpc]

    // Functions
    void resetEvent();
    void resetMuon();
    void resetMichel();
    bool Log(bool cond, int flag, int tab, std::string msg, std::string succ, std::string fail);

    size_t GetSection(int ch);
    geo::View_t GetPlane(raw::ChannelID_t ch, int sec);
    geo::WireID GetWireID(geo::Point_t const& P, geo::View_t plane);
    raw::ChannelID_t GetChannel(geo::Point_t const& P, geo::View_t plane);

    bool IsInVolume(double x, double y, double z, float eps);
    bool IsInUpperVolume(double x, double y, double z, float eps);
    bool IsInLowerVolume(double x, double y, double z, float eps);
    bool IsInVolume(float x, float y, float z, float eps);
    bool IsInUpperVolume(float x, float y, float z, float eps);
    bool IsInLowerVolume(float x, float y, float z, float eps);
    template<class Vec> bool IsInVolume(Vec const& V, float eps);
    template<class Vec> bool IsInUpperVolume(Vec const& V, float eps);
    template<class Vec> bool IsInLowerVolume(Vec const& V, float eps);

    double Dist(TLorentzVector u, TLorentzVector v);

    bool IsInsideWindow(float tick, float eps);
    bool IsInsideTPC(int channel, int eps);

    bool IsUpright(recob::Track const& T);

    // int iCorr(int n, int i, int j);
};


ana::Muchecks::Muchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),

    fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("TrackLengthCut")), // in cm

    fKeepOutside(p.get<bool>("KeepOutside")),
    fKeepNonDecaying(p.get<bool>("KeepNonDecaying"))
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::Muchecks: =============================================================" << "\033[0m" << std::endl;
    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    geo::WireGeo const wiregeo1 = asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 0});
    geo::WireGeo const wiregeo2 = asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 1});
    fChannelPitch = geo::WireGeo::WirePitch(wiregeo1, wiregeo2);

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    // Retrieving product tags
    for (std::vector<std::string> prod : vvsProducts) {

        const std::string   process     = prod[0],
                            label       = prod[1],
                            instance    = prod[2],
                            type        = prod[3];

        const art::InputTag tag = art::InputTag(label,instance);

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
    }



    tEvent = tfs->make<TTree>("Event","");

    tEvent->Branch("iEvent", &iEvent);

    tEvent->Branch("NHit", &EventNHit);
    tEvent->Branch("HitChannel", &EventHitChannel);
    tEvent->Branch("HitTick", &EventHitTick);
    tEvent->Branch("HitADC", &EventHitADC);

    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);

    tMuon = tfs->make<TTree>("Muon","");

    tMuon->Branch("iEvent", &iEvent);
    tMuon->Branch("iMuon", &iMuon);

    tMuon->Branch("IsAnti", &MuonIsAnti);
    tMuon->Branch("DoesDecay", &MuonDoesDecay);
    tMuon->Branch("HasMichel", &MuonHasMichel);

    tMuon->Branch("TrackLength", &MuonTrackLength);

    tMuon->Branch("NTrackPoint", &MuonNTrackPoint);
    tMuon->Branch("TrackPointX", &MuonTrackPointX);
    tMuon->Branch("TrackPointY", &MuonTrackPointY);
    tMuon->Branch("TrackPointZ", &MuonTrackPointZ);

    tMuon->Branch("EndChannel", &MuonEndChannel);
    tMuon->Branch("EndTick", &MuonEndTick);
    tMuon->Branch("EndIsInWindow", &MuonEndIsInWindow);
    tMuon->Branch("EndIsInTPC", &MuonEndIsInTPC);
    tMuon->Branch("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ);

    tMuon->Branch("NHit", &MuonNHit);
    tMuon->Branch("HitChannel", &MuonHit.channel);
    tMuon->Branch("HitTick", &MuonHit.tick);
    tMuon->Branch("HitADC", &MuonHit.adc);

    // tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy);
    // tMuon->Branch("MichelTrackLength", &MichelTrackLength);
    // // tMuon->Branch("MichelChannelBary", &MichelChannelBary);
    // // tMuon->Branch("MichelTickBary", &MichelTickBary);

    // tMuon->Branch("MichelNHit", &MichelNHit);
    // tMuon->Branch("MichelHitChannel", &MichelHitChannel);
    // tMuon->Branch("MichelHitTick", &MichelHitTick);
    // tMuon->Branch("MichelHitADC", &MichelHitADC);
    // tMuon->Branch("MichelHitEnergy", &MichelHitEnergy);

    // tMuon->Branch("MichelSphereNHit", &MichelSphereNHit);
    // tMuon->Branch("MichelSphereChannel", &MichelSphereChannel);
    // tMuon->Branch("MichelSphereTick", &MichelSphereTick);
    // tMuon->Branch("MichelSphereADC", &MichelSphereADC);
    // tMuon->Branch("MichelSphereEnergy", &MichelSphereEnergy);

    tMichel = tfs->make<TTree>("Michel","");

    tMichel->Branch("iEvent", &iEvent);
    tMichel->Branch("iMuon", &MicheliMuon);

    tMichel->Branch("MichelTrueEnergy", &MichelTrueEnergy);
    tMichel->Branch("MichelTrackLength", &MichelTrackLength);

    tMichel->Branch("MichelNHit", &MichelNHit);
    tMichel->Branch("MichelHitChannel", &MichelHitChannel);
    tMichel->Branch("MichelHitTick", &MichelHitTick);
    tMichel->Branch("MichelHitADC", &MichelHitADC);
    tMichel->Branch("MichelHitEnergy", &MichelHitEnergy);

    tMichel->Branch("MichelSphereNHit", &MichelSphereNHit);
    tMichel->Branch("MichelSphereChannel", &MichelSphereChannel);
    tMichel->Branch("MichelSphereTick", &MichelSphereTick);
    tMichel->Branch("MichelSphereADC", &MichelSphereADC);
    tMichel->Branch("MichelSphereEnergy", &MichelSphereEnergy);

    tMichel->Branch("SphereTruePositive", &MichelSphereTruePositive);
    tMichel->Branch("SphereFalsePositive", &MichelSphereFalsePositive);


    // Diagnostic Variables
    binTick.n = detProp.ReadOutWindowSize()/4;
    binTick.min = 0;
    binTick.max = detProp.ReadOutWindowSize();

    binChan = std::vector<ana::Binning>(n_section);
    for (unsigned int s=0; s<n_section; s++) {
        geo::PlaneID planeid1{geo::TPCID{0, n_section*s}, geo::kW};
        geo::PlaneID palenid2{geo::TPCID{0, n_section*(s+1)-1}, geo::kW};
        binChan[s].min = asWire->PlaneWireToChannel(geo::WireID{planeid1, 0});
        binChan[s].max = asWire->PlaneWireToChannel(geo::WireID{palenid2, asWire->Nwires(palenid2)-1});
        binChan[s].n = binChan[s].max - binChan[s].min + 1;
        if (iLogLevel >= kDetails) std::cout << "S#" << s << " channel range: " << binChan[s].min << " - " << binChan[s].max << std::endl;
    }

    chanTPC = std::vector<ana::Binning>(asGeo->NTPC());
    for (unsigned int tpc=0; tpc<asGeo->NTPC(); tpc++) {
        geo::TPCID tpcid{0, tpc};
        geo::PlaneID planeid{tpcid, geo::kW};
        chanTPC[tpc].min = asWire->PlaneWireToChannel(geo::WireID{planeid, 0});
        chanTPC[tpc].max = asWire->PlaneWireToChannel(geo::WireID{planeid, asWire->Nwires(planeid)-1});
        chanTPC[tpc].n = chanTPC[tpc].max - chanTPC[tpc].min + 1;
        if (iLogLevel >= kDetails) std::cout << "TPC#" << tpc << " channel range: " << chanTPC[tpc].min << " - " << chanTPC[tpc].max << std::endl;
    }

    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    float fMichelTickRadius = fMichelTimeRadius / fSamplingRate;
    int fMichelChannelRadius = int(fMichelSpaceRadius / fChannelPitch);

    if (iLogLevel >= kDetails) {
        std::cout << "fMichelTimeRadius: " << fMichelTimeRadius << " µs" << std::endl;
        std::cout << "fSamplingRate: " << fSamplingRate << " µs/tick" << std::endl;
        std::cout << "fMichelTickRadius: " << fMichelTickRadius << " ticks" << std::endl;
        std::cout << "fMichelSpaceRadius: " << fMichelSpaceRadius << " cm" << std::endl;
        std::cout << "fChannelPitch: " << fChannelPitch << " cm/channel" << std::endl;
        std::cout << "fMichelChannelRadius: " << fMichelChannelRadius << " channels" << std::endl;
    }


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Muchecks::Muchecks =======================================================" << "\033[0m" << std::endl;
}

void ana::Muchecks::analyze(art::Event const& e)
{

    if (iLogLevel >= kBasics) {
        std::cout << "\033[93m" << "Muchecks::analyze: Initialization evt#" << iEvent << "\r" << std::flush;
        std::cout << std::string(5,'\t') << " ======================================" << "\033[0m" << std::endl;
    }

    resetEvent();

    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    float fMichelTickRadius = fMichelTimeRadius / fSamplingRate;
    int fMichelChannelRadius = int(fMichelSpaceRadius / fChannelPitch);

    std::vector<int> visited;
        
    auto const & vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    for (recob::Hit const& hit : *vh_hit) {
        if (hit.View() != geo::kW) continue;

        EventNHit++;
        EventHitChannel.push_back(hit.Channel());
        EventHitTick.push_back(hit.PeakTime());
        EventHitADC.push_back(hit.Integral());
    } // end loop over hits

    auto const vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindManyP<recob::SpacePoint> fmp_hit2spt(vh_hit, e, tag_spt);
    art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Muon Ends ===================================================" << "\033[0m" << std::endl;

    std::vector<int> MuonEndChannelAll;
    std::vector<float> MuonEndTickAll;
    std::vector<simb::MCParticle const*> MCPMichelAll;


    if (iLogLevel >= kInfos) std::cout << "looping over " << vp_trk.size() << " tracks..." << std::endl;
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        if (iLogLevel >= kInfos) std::cout << "trk#" << p_trk->ID() << "\r" << std::flush;

        if (!Log(p_trk->Length() > fTrackLengthCut, kDetails, 1, "is long enough...", "yes", "no")
        ) continue;

        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());

        if (!Log(mcp, kDetails, 1, "trk to mcp...", "done", "failed")
        ) continue;    

        if (!Log(abs(mcp->PdgCode()) == 13, kDetails, 1, "is muon...", "yes", "no")
        ) continue;

        // if (Log(mcp->EndProcess() == "Decay", kDetails, 1, "is decaying...", "yes", "no")
        // ) continue;

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());

        if (!Log(vp_hit_muon.size() > 0, kDetails, 1, "has hits...", "yes", "no")
        ) continue;

        resetMuon();

        float TickUpMax = 0;
        float TickLowMin = detProp.ReadOutWindowSize();
        int ChanUpMax = -1;
        int ChanLowMin = -1;

        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            if (p_hit_muon->View() != geo::kW) continue;

            if (p_hit_muon->Channel() > binChan[1].max) {
                if (p_hit_muon->PeakTime() > TickUpMax) {
                    TickUpMax = p_hit_muon->PeakTime();
                    ChanUpMax = p_hit_muon->Channel();
                }
            }
            if (p_hit_muon->Channel() < binChan[2].min) {
                if (p_hit_muon->PeakTime() < TickLowMin) {
                    TickLowMin = p_hit_muon->PeakTime();
                    ChanLowMin = p_hit_muon->Channel();
                }
            }

            MuonNHit++;
            MuonHit.push_back({p_hit_muon->Channel(), p_hit_muon->PeakTime(), p_hit_muon->Integral()});
        }

        if (iLogLevel >= kDetails) {
            if (ChanUpMax == -1) std::cout << "\t" << "\033[91m" << "no hit in upper volume" << "\033[0m" << std::endl;
            else {
                std::cout << "\t" << "up max channel: " << "\033[93m" << ChanUpMax << "\033[0m" << std::endl;
                std::cout << "\t" << "up max tick: " << "\033[93m" << TickUpMax << "\033[0m" << std::endl;
            }
            std::cout << "\t" << "------------------------------------------------" << std::endl;
            if (ChanLowMin == -1) std::cout << "\t" << "\033[91m" << "no hit in lower volume" << "\033[0m" << std::endl;
            else {
                std::cout << "\t" << "low min channel: " << "\033[93m" << ChanLowMin << "\033[0m" << std::endl;
                std::cout << "\t" << "low min tick: " << "\033[93m" << TickLowMin << "\033[0m" << std::endl;
            }
            std::cout << "\t" << "------------------------------------------------" << std::endl;
        }

        if (ChanLowMin == -1) {
            if (ChanUpMax == -1) continue;
            else {
                MuonEndChannel = ChanUpMax;
                MuonEndTick = TickUpMax;
            }
        } else {
            MuonEndChannel = ChanLowMin;
            MuonEndTick = TickLowMin;
        }

        MuonEndChannelAll.push_back(MuonEndChannel);
        MuonEndTickAll.push_back(MuonEndTick);
    
        // Analyzing muon

        MuonEndIsInWindow = IsInsideWindow(MuonEndTick, fMichelTickRadius);
        MuonEndIsInTPC = IsInsideTPC(MuonEndChannel, fMichelChannelRadius);

        geo::Point_t end;
        if (IsUpright(*p_trk)) end = p_trk->End();
        else end = p_trk->Start();
        MuonEndIsInVolumeYZ = IsInUpperVolume(50., end.Y(), end.Z(), fMichelSpaceRadius);

        Log(MuonEndIsInWindow, kInfos, 1, Form("is in window (±%.1f ticks)...", fMichelTickRadius), "yes", "no");
        Log(MuonEndIsInTPC, kInfos, 1, Form("is in TPC (±%d chans)...", fMichelChannelRadius), "yes", "no");
        Log(MuonEndIsInVolumeYZ, kInfos, 1, Form("is in volumeYZ (±%.1f cm)...", fMichelSpaceRadius), "yes", "no");

        if (!fKeepOutside && !MuonEndIsInWindow && !MuonEndIsInTPC && !MuonEndIsInVolumeYZ) continue;

        if (iLogLevel >= kInfos) std::cout << "\t" << "mu" << (MuonIsAnti ? "+" : "-") << "\033[93m" << " #" << EventNMuon << " (" << iMuon << ")" << "\033[0m" << " found" << std::endl;

        if (iLogLevel >= kInfos) std::cout << "\t" << "muon end @ channel:" << "\033[93m" << MuonEndChannel << "\033[0m" 
            << " tick:" << "\033[93m" << MuonEndTick << "\033[0m"
            << " YZ:" << "\033[93m(" << end.Y() << "," << end.Z() << ")\033[0m" << std::endl;

        EventNMuon++;
        EventiMuon.push_back(iMuon);

        MuonIsAnti = int(mcp->PdgCode() < 0);
        MuonTrackLength = p_trk->Length();
        MuonDoesDecay = mcp->EndProcess() == "Decay";

        if (!fKeepNonDecaying && !MuonDoesDecay) continue;

        for (size_t i_tpt=0; i_tpt<p_trk->NumberTrajectoryPoints(); i_tpt++) {
            if (!p_trk->HasValidPoint(i_tpt)) continue;
            MuonTrackPointX.push_back(p_trk->LocationAtPoint(i_tpt).X());
            MuonTrackPointY.push_back(p_trk->LocationAtPoint(i_tpt).Y());
            MuonTrackPointZ.push_back(p_trk->LocationAtPoint(i_tpt).Z());
            MuonNTrackPoint++;
        }

        int i_dau = mcp->NumberDaughters() - 1;
        if (iLogLevel >= kInfos) std::cout << "\t" << "looping over muon's " << mcp->NumberDaughters() << " daughters..." << std::endl;
        for (; i_dau >= 0; i_dau--) {
            if (iLogLevel >= kDetails) std::cout << "\t" << "dau#" << i_dau+1 << "\r" << std::flush;
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            if (mcp_dau->Mother() != mcp->TrackId()) anomalies.push_back(Form("e%ld t%d µ%ld: TrackID: mcp_dau->Mother() (%d) != mcp->TrackId() (%d)", iEvent, p_trk->ID(), iMuon, mcp_dau->Mother(), mcp->TrackId()));

            if (!Log(mcp_dau, kDetails, 2, "id to mcp...", "done", "failed")
            ) continue;

            if (!Log(abs(mcp_dau->PdgCode()) == 11 && mcp_dau->Process() == "Decay", kDetails, 2, "is michel...", "yes", "no")
            ) continue;

            break;
        } // end loop over muon daughters
        if (i_dau == -1) {
            if (iLogLevel >= kInfos) std::cout << "\t" << "no michel found" << std::endl;

            MuonHasMichel = kNoMichel;

            MCPMichelAll.push_back(nullptr);
        }
        else {
            if (iLogLevel >= kInfos) std::cout << "\t" << "michel found" << i_dau << std::endl;

            simb::MCParticle const * mcp_michel = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));

            MCPMichelAll.push_back(mcp_michel);

            // recob::Track const * trk_michel = truthUtil.GetRecoTrackFromMCParticle(clockData, *mcp_michel, e, tag_trk.label());

            // if (trk_michel) {
            //     MichelTrackLength = trk_michel->Length();
            //     if (iLogLevel >= kInfos) std::cout << "\t" << "michel track length: " << "\033[93m" << MichelTrackLength << "\033[0m" << std::endl;
            // }

            // if (std::find(visited.begin(), visited.end(), mcp_michel->TrackId()) != visited.end()) anomalies.push_back(Form("e%ld t%d µ%ld: TrackID: michel already visited (%d)", iEvent, p_trk->ID(), iMuon, mcp_michel->TrackId()));
            // visited.push_back(mcp_michel->TrackId());

            // MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;

            // double dist = Dist(mcp->EndPosition(), mcp_michel->Position(0));
            // if (dist > 0) anomalies.push_back(Form("e%ld t%d µ%ld: distance of %f", iEvent, p_trk->ID(), iMuon, dist)); 
            // if (MichelTrueEnergy > 100) anomalies.push_back(Form("e%ld t%d µ%ld: Michel True Energy %.1f MeV", iEvent, p_trk->ID(), iMuon, MichelTrueEnergy)); 

            // if (iLogLevel >= kInfos) std::cout << "\t" << "michel true energy: " << "\033[93m" << MichelTrueEnergy << " MeV" << "\033[0m" << std::endl;

            if (Log(IsInVolume(mcp_michel->Position(0), fMichelSpaceRadius), kInfos, 1, "michel is inside...", "yes", "no")
            ) {
                MuonHasMichel = kHasMichelInside; 
            } else { 
                MuonHasMichel = kHasMichelOutside;
            }

            // std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());

	
            // float prev = 0;
            // for (const recob::Hit* hit_michel : v_hit_michel) {
            //     if (hit_michel->View() != geo::kW) continue;

            //     MichelNHit++;
            //     // MichelChannelBary += hit_michel->Channel();
            //     // MichelTickBary += hit_michel->PeakTime();
            //     MichelHitChannel.push_back(hit_michel->Channel());
            //     MichelHitTick.push_back(hit_michel->PeakTime());
            //     MichelHitADC.push_back(hit_michel->Integral());
            //     if (hit_michel->ROISummedADC() != prev) {
            //         MichelHitEnergy += hit_michel->ROISummedADC();
            //         prev = hit_michel->ROISummedADC();
            //     }
            // }
            // MichelHitEnergy *= fADCtoE;

            // if (MichelNHit > 0) {
            //     MichelChannelBary /= MichelNHit;
            //     MichelTickBary /= MichelNHit;
            // }

            // if (iLogLevel >= kInfos) {
            //     if (MichelNHit == 0) std::cout << "\t" << "\033[91m" << "no hit in michel" << "\033[0m" << std::endl;
            //     else {
            //         std::cout << "\t" << "michel hits: " << "\033[93m" << MichelNHit << "\033[0m" << std::endl;
            //         // std::cout << "\t" << "michel barycenter channel: " << "\033[93m" << MichelChannelBary << "\033[0m" << std::endl;
            //         // std::cout << "\t" << "michel barycenter tick: " << "\033[93m" << MichelTickBary << "\033[0m" << std::endl;
            //         std::cout << "\t" << "michel energy: " << "\033[93m" << MichelHitEnergy << "\033[0m" << std::endl;
            //     }
            //     std::cout << "\t" << "michel end process: " << "\033[93m" << mcp_michel->EndProcess() << "\033[0m" << std::endl;
            // }
        
            // prev = 0;
            // for (const art::Ptr<recob::Hit> p_hit : vp_hit) {
            //     if (p_hit->View() != geo::kW) continue;


            //     std::vector<art::Ptr<recob::Track>> vp_trk = fmp_hit2trk.at(p_hit.key());
            //     auto it = vp_trk.begin();
            //     while (it != vp_trk.end() && (*it)->Length() > fTrackLengthCut) it++;
            //     if (it != vp_trk.end()) continue;


            //     float X = (float(p_hit->Channel()) - MuonEndChannel) / fMichelChannelRadius;
            //     float Y = (p_hit->PeakTime() - MuonEndTick) / fMichelTickRadius;
            //     if (X*X + Y*Y > 1) continue;

            //     MichelSphereNHit++;
            //     MichelSphereChannel.push_back(p_hit->Channel());
            //     MichelSphereTick.push_back(p_hit->PeakTime());
            //     MichelSphereADC.push_back(p_hit->Integral());
            //     if (p_hit->ROISummedADC() != prev) {
            //         MichelSphereEnergy += p_hit->ROISummedADC();
            //         prev = p_hit->ROISummedADC();
            //     }
            // }
            // MichelSphereEnergy *= fADCtoE;

            // if (iLogLevel >= kInfos) {
            //     std::cout << "\t" << "michel sphere hits: " << "\033[93m" << MichelSphereNHit << "\033[0m" << std::endl;
            //     std::cout << "\t" << "michel sphere energy: " << "\033[93m" << MichelSphereEnergy << "\033[0m" << std::endl;
            // }
        }

        tMuon->Fill();
        iMuon++;
    } // end loop over tracks

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Michel Hits =================================================" << "\033[0m" << std::endl;

    std::vector<std::vector<const recob::Hit*>> tpMichelHitAll;
    for (const simb::MCParticle * mcp_michel : MCPMichelAll) {
        if (!mcp_michel) {
            tpMichelHitAll.push_back({});
            continue;
        }
        tpMichelHitAll.push_back(truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label()));
    }

    std::vector<int> tpMSNHit(EventNMuon, 0);
    std::vector<std::vector<int>> tpMSChannel(EventNMuon);
    std::vector<std::vector<float>> tpMSTick(EventNMuon);
    std::vector<std::vector<float>> tpMSADC(EventNMuon);
    std::vector<float> tpMSEnergy(EventNMuon, 0);
    
    std::vector<size_t> tpMSTruePositive(EventNMuon, 0);
    std::vector<size_t> tpMSFalsePositive(EventNMuon, 0);
    
    std::vector<float> tp_prev(EventNMuon, 0);
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;

        std::vector<art::Ptr<recob::Track>> vp_trk = fmp_hit2trk.at(p_hit.key());
        auto it = vp_trk.begin();
        while (it != vp_trk.end() && (*it)->Length() > fTrackLengthCut) it++;
        if (it != vp_trk.end()) continue;

        for (size_t m = 0; m<EventNMuon; m++) {
            if (!MCPMichelAll[m]) continue;

            float X = (float(p_hit->Channel()) - MuonEndChannelAll[m]) / fMichelChannelRadius;
            float Y = (p_hit->PeakTime() - MuonEndTickAll[m]) / fMichelTickRadius;

            if (X*X + Y*Y > 1) continue;

            if (std::find(tpMichelHitAll[m].begin(), tpMichelHitAll[m].end(), &*p_hit) != tpMichelHitAll[m].end()) {
                tpMSTruePositive[m]++;
            } else {
                tpMSFalsePositive[m]++;
            }

            tpMSNHit[m]++;
            tpMSChannel[m].push_back(p_hit->Channel());
            tpMSTick[m].push_back(p_hit->PeakTime());
            tpMSADC[m].push_back(p_hit->Integral());

            if (p_hit->ROISummedADC() != tp_prev[m]) {
                tpMSEnergy[m] += p_hit->ROISummedADC();
                tp_prev[m] = p_hit->ROISummedADC();
            }
        }
    }
    tp_prev.clear();

    for (size_t m=0; m<EventNMuon; m++) {


        if (iLogLevel >= kInfos) std::cout << "mu#" << m << "\r" << std::flush;

        resetMichel();

        MicheliMuon = EventiMuon[m];

        const simb::MCParticle * mcp_michel = MCPMichelAll[m];

        if (!mcp_michel) continue;

        recob::Track const * trk_michel = truthUtil.GetRecoTrackFromMCParticle(clockData, *mcp_michel, e, tag_trk.label());

        if (trk_michel) {
            MichelTrackLength = trk_michel->Length();
            if (iLogLevel >= kInfos) std::cout << "\t" << "michel track length: " << "\033[93m" << MichelTrackLength << "\033[0m" << std::endl;
        }

        if (std::find(visited.begin(), visited.end(), mcp_michel->TrackId()) != visited.end()) anomalies.push_back(Form("e%ld µ%ld: TrackID: michel already visited (%d)", iEvent, iMuon, mcp_michel->TrackId()));
        visited.push_back(mcp_michel->TrackId());

        MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;
        if (iLogLevel >= kInfos) std::cout << "\t" << "michel true energy: " << "\033[93m" << MichelTrueEnergy << " MeV" << "\033[0m" << std::endl;

        // double dist = Dist(mcp->EndPosition(), mcp_michel->Position(0));
        // if (dist > 0) anomalies.push_back(Form("e%ld t%d µ%ld: distance of %f", iEvent, p_trk->ID(), iMuon, dist)); 
        if (MichelTrueEnergy > 100) anomalies.push_back(Form("e%ld µ%ld: Michel True Energy %.1f MeV", iEvent, iMuon, MichelTrueEnergy)); 



        float prev = 0;
        for (const recob::Hit* hit_michel : tpMichelHitAll[m]) {
            if (hit_michel->View() != geo::kW) continue;

            MichelNHit++;
            MichelHitChannel.push_back(hit_michel->Channel());
            MichelHitTick.push_back(hit_michel->PeakTime());
            MichelHitADC.push_back(hit_michel->Integral());
            if (hit_michel->ROISummedADC() != prev) {
                MichelHitEnergy += hit_michel->ROISummedADC();
                prev = hit_michel->ROISummedADC();
            }
        }
        MichelHitEnergy *= fADCtoE;

        if (iLogLevel >= kInfos) {
            if (MichelNHit == 0) std::cout << "\t" << "\033[91m" << "no hit in michel" << "\033[0m" << std::endl;
            else {
                std::cout << "\t" << "michel hits: " << "\033[93m" << MichelNHit << "\033[0m" << std::endl;
                std::cout << "\t" << "michel hit energy: " << "\033[93m" << MichelHitEnergy << "\033[0m" << std::endl;
            }
            std::cout << "\t" << "michel end process: " << "\033[93m" << mcp_michel->EndProcess() << "\033[0m" << std::endl;
        }

        MichelSphereNHit = tpMSNHit[m];
        MichelSphereChannel = tpMSChannel[m];
        MichelSphereTick = tpMSTick[m];
        MichelSphereADC = tpMSADC[m];
        MichelSphereEnergy = tpMSEnergy[m] * fADCtoE;

        MichelSphereTruePositive = tpMSTruePositive[m];
        MichelSphereFalsePositive = tpMSFalsePositive[m];

        if (iLogLevel >= kInfos) {
            std::cout << "\t" << "michel sphere hits: " << "\033[93m" << MichelSphereNHit << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere energy: " << "\033[93m" << MichelSphereEnergy << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere true positive: " << "\033[93m" << MichelSphereTruePositive << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere false positive: " << "\033[93m" << MichelSphereFalsePositive << "\033[0m" << std::endl;
        }

        tMichel->Fill();
    } // end loop over michels

    tEvent->Fill();
    iEvent++;

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Muchecks::analyze =======================================================" << "\033[0m" << std::endl;
} // end analyze

void ana::Muchecks::beginJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::beginJob: ============================================================" << "\033[0m" << std::endl;
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Muchecks::beginJob ======================================================" << "\033[0m" << std::endl;
} // end beginJob


void ana::Muchecks::endJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::endJob: ==============================================================" << "\033[0m" << std::endl;

    if (iLogLevel >= kImportant) {
        std::cout << "\033[91m" << "Anomalies:" << "\033[0m" << std::endl;
        for (std::string anomaly : anomalies) std::cout << anomaly << std::endl;
    }

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Muchecks::endJob ========================================================" << "\033[0m" << std::endl;
} // end endJob


void ana::Muchecks::resetEvent() {
    if (iLogLevel >= kDetails) std::cout << "resetting event branches...";

    EventNHit = 0;
    EventHitChannel.clear();
    EventHitTick.clear();
    EventHitADC.clear();

    EventNMuon = 0;
    EventiMuon.clear();

    if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
}
void ana::Muchecks::resetMuon() {
    if (iLogLevel >= kDetails) std::cout << "\t" << "resetting muon branches...";

    MuonNTrackPoint = 0;
    MuonTrackPointX.clear();
    MuonTrackPointY.clear();
    MuonTrackPointZ.clear();

    MuonNHit = 0;
    MuonHit.clear();

    // MichelNHit = 0;
    // MichelTrackLength = 0;
    // // MichelChannelBary = 0;
    // // MichelTickBary = 0;
    // MichelHitChannel.clear();
    // MichelHitTick.clear();
    // MichelHitADC.clear();
    // MichelTrueEnergy = 0;
    // MichelHitEnergy = 0;

    // MichelSphereNHit = 0;
    // MichelSphereChannel.clear();
    // MichelSphereTick.clear();
    // MichelSphereADC.clear();
    // MichelSphereEnergy = 0;

    if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
}
void ana::Muchecks::resetMichel() {
    if (iLogLevel >= kDetails) std::cout << "\t" << "resetting michel branches...";
;
    MichelTrueEnergy = 0;
    MichelTrackLength = 0;

    MichelNHit = 0;
    MichelHitChannel.clear();
    MichelHitTick.clear();
    MichelHitADC.clear();
    MichelHitEnergy = 0;

    MichelSphereNHit = 0;
    MichelSphereChannel.clear();
    MichelSphereTick.clear();
    MichelSphereADC.clear();
    MichelSphereEnergy = 0;

    if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
}

bool ana::Muchecks::Log(bool cond, int flag, int tab, std::string msg, std::string succ, std::string fail) {
    if (iLogLevel >= flag) {
        std::cout << std::string(tab,'\t') << msg << " ";
        if (cond) std::cout << "\033[92m" << succ << "\033[0m" << std::endl;
        else std::cout << "\033[91m" << fail << "\033[0m" << std::endl;
    }
    return cond;
}




size_t ana::Muchecks::GetSection(int ch) {
    return int(4.*ch / asWire->Nchannels());
}
geo::View_t ana::Muchecks::GetPlane(raw::ChannelID_t ch, int sec) {
    return static_cast<geo::View_t>(12.*ch / asWire->Nchannels() - 3*sec);
}
geo::WireID ana::Muchecks::GetWireID(geo::Point_t const& P, geo::View_t plane) {
    geo::TPCID tpcid = asGeo->FindTPCAtPosition(P);
    if (!tpcid.isValid) return geo::WireID();
    geo::PlaneGeo const& planegeo = asWire->Plane(tpcid,plane);
    geo::WireID wireid;
    try {
        wireid = planegeo.NearestWireID(P);
    } catch (geo::InvalidWireError const& e) {
        return e.suggestedWireID();
    }
    return wireid;
}
raw::ChannelID_t ana::Muchecks::GetChannel(geo::Point_t const& P, geo::View_t plane) {
    geo::WireID wireid = GetWireID(P, plane);
    if (!wireid.isValid) return raw::InvalidChannelID;
    return asWire->PlaneWireToChannel(wireid);
}

bool ana::Muchecks::IsInLowerVolume(double x, double y, double z, float eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpcid0{cryoid, 0}, tpcidN{cryoid, 7};
    return (x-eps > asGeo->TPC(tpcid0).MinX() and x+eps < asGeo->TPC(tpcidN).MaxX() and
            y-eps > asGeo->TPC(tpcid0).MinY() and y+eps < asGeo->TPC(tpcidN).MaxY() and
            z-eps > asGeo->TPC(tpcid0).MinZ() and z+eps < asGeo->TPC(tpcidN).MaxZ());   
}
bool ana::Muchecks::IsInUpperVolume(double x, double y, double z, float eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpcid0{cryoid, 8}, tpcidN{cryoid, asGeo->NTPC(cryoid)-1};
    return (x-eps > asGeo->TPC(tpcid0).MinX() and x+eps < asGeo->TPC(tpcidN).MaxX() and
            y-eps > asGeo->TPC(tpcid0).MinY() and y+eps < asGeo->TPC(tpcidN).MaxY() and
            z-eps > asGeo->TPC(tpcid0).MinZ() and z+eps < asGeo->TPC(tpcidN).MaxZ());   
}
bool ana::Muchecks::IsInVolume(double x, double y, double z, float eps) {
    return IsInLowerVolume(x, y, z, eps) or IsInUpperVolume(x, y, z, eps);
}
bool ana::Muchecks::IsInLowerVolume(float x, float y, float z, float eps) {
    return IsInLowerVolume(double(x), double(y), double(z), eps);
    return IsInUpperVolume(double(x), double(y), double(z), eps);
}
bool ana::Muchecks::IsInVolume(float x, float y, float z, float eps) {
    return IsInVolume(double(x), double(y), double(z), double(eps));
}
template<class Vec>
bool ana::Muchecks::IsInUpperVolume(Vec const& V, float eps) {
    return IsInUpperVolume(V.X(), V.Y(), V.Z(), eps);
}
template<class Vec>
bool ana::Muchecks::IsInLowerVolume(Vec const& V, float eps) {
    return IsInLowerVolume(V.X(), V.Y(), V.Z(), eps);
}
template<class Vec>
bool ana::Muchecks::IsInVolume(Vec const& V, float eps) {
    return IsInVolume(V.X(), V.Y(), V.Z(), eps);
}


double ana::Muchecks::Dist(TLorentzVector u, TLorentzVector v) {
    double X = u.X() - v.X();
    double Y = u.Y() - v.Y();
    double Z = u.Z() - v.Z();
    return X*X + Y*Y + Z*Z;
}

bool ana::Muchecks::IsInsideWindow(float tick, float eps) {
    return binTick.min + eps < tick and tick < binTick.max - eps;
}
bool ana::Muchecks::IsInsideTPC(int ch, int eps) {
    bool isIn = false;
    for (unsigned int tpc=0; tpc<asGeo->NTPC(); tpc++) {
        isIn = chanTPC[tpc].min + eps < ch and ch < chanTPC[tpc].max - eps;
        if (isIn) break;
    }
    return isIn;
}


bool ana::Muchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}

DEFINE_ART_MODULE(ana::Muchecks)
