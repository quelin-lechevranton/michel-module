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
        float min, max;
    };
    struct Hit {
        unsigned int slice;
        float Z;
        int channel;
        float tick;
        float adc;

        Hit() : 
            slice(0),
            Z(0), 
            channel(0),
            tick(0), 
            adc(0) {}
        Hit(unsigned int s, float z, int c, float t, float a) : 
            slice(s), 
            Z(z), 
            channel(c),
            tick(t), 
            adc(a) {}
    };
    struct Hits {
        std::vector<unsigned int> slice;
        std::vector<float> Z;
        std::vector<int> channel;
        std::vector<float> tick;
        std::vector<float> adc;

        Hits() : slice(), Z(), channel(), tick(), adc() {}
        void push_back(Hit const& hit) {
            slice.push_back(hit.slice);
            Z.push_back(hit.Z);
            channel.push_back(hit.channel);
            tick.push_back(hit.tick);
            adc.push_back(hit.adc);
        }
        void clear() {
            slice.clear();
            Z.clear();
            channel.clear();
            tick.clear();
            adc.clear();
        }
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
    float fMichelSpaceRadius; // in cm
    float fTrackLengthCut; // in cm
    float fMichelTimeRadius; // in µs
    float fMichelTickRadius; // in ticks
    int fMichelChannelRadius; // in channels

    bool fKeepOutside,
         fKeepNonDecaying;
    
    // Output Variables

    std::vector<std::string> anomalies;

    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside };

    TTree* tEvent;

    size_t iEvent=0;

    size_t EventNHit;
    ana::Hits EventHits; //[hit]

    size_t EventNMuon;
    std::vector<int> EventiMuon; //[muon]

    // enum Tag { kMCP, kDaughters, kSphere, kTrue, kNTag };
    TTree* tMuon;

    size_t iMuon=0;

    int MuonIsAnti;
    int MuonDoesDecay;
    int MuonHasMichel;

    float MuonTrackLength;

    size_t MuonNTrackPoint;
    std::vector<float> MuonTrackPointX, MuonTrackPointY, MuonTrackPointZ;

    ana::Hit MuonEnd;
    int MuonEndIsInWindow;
    int MuonEndIsInVolumeZ;
    int MuonEndIsInVolumeYZ;

    size_t MuonNHit;
    ana::Hits MuonHits; //[hit]

    float MichelTrackLength;
    float MichelTrueEnergy;
    float MichelHitEnergy;

    size_t MichelNHit;
    ana::Hits MichelHits; //[hit]

    size_t SphereNHit;
    ana::Hits SphereHits; //[hit]
    float SphereEnergy;

    size_t SphereTruePositive;
    size_t SphereFalsePositive;

    float SphereEnergyTruePositive;
    float SphereEnergyFalsePositive;

    // Diagnostic Variables
    unsigned int n_section = 4;
    ana::Binning binTick;
    std::vector<ana::Binning> binChan; //[section]
    std::vector<ana::Binning> chanTPC; //[tpc]

    std::map<int, float> mapChanZ;


    // Functions
    void resetEvent();
    void resetMuon();
    bool Log(bool cond, int flag, int tab, std::string msg, std::string succ, std::string fail);

    size_t GetSection(int ch);
    unsigned int GetSlice(int ch);
    ana::Hit GetHit(int ch, float t, float a);
    ana::Hit GetHit(recob::Hit const& hit);
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
    bool IsInsideZ(float z, float eps);

    bool IsUpright(recob::Track const& T);

    // int iCorr(int n, int i, int j);
};



ana::Muchecks::Muchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),

    // fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
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
    tEvent->Branch("HitSlice", &EventHits.slice);
    tEvent->Branch("HitZ", &EventHits.Z);
    tEvent->Branch("HitChannel", &EventHits.channel);
    tEvent->Branch("HitTick", &EventHits.tick);
    tEvent->Branch("HitADC", &EventHits.adc);

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

    tMuon->Branch("EndSlice", &MuonEnd.slice);
    tMuon->Branch("EndZ", &MuonEnd.Z);
    tMuon->Branch("EndChannel", &MuonEnd.channel);
    tMuon->Branch("EndTick", &MuonEnd.tick);

    tMuon->Branch("EndIsInWindow", &MuonEndIsInWindow);
    tMuon->Branch("EndIsInVolumeZ", &MuonEndIsInVolumeZ);
    tMuon->Branch("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ);

    tMuon->Branch("NHit", &MuonNHit);
    tMuon->Branch("HitSlice", &MuonHits.slice);
    tMuon->Branch("HitZ", &MuonHits.Z);
    tMuon->Branch("HitChannel", &MuonHits.channel);
    tMuon->Branch("HitTick", &MuonHits.tick);
    tMuon->Branch("HitADC", &MuonHits.adc);

    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy);
    tMuon->Branch("MichelTrackLength", &MichelTrackLength);
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy);

    tMuon->Branch("MichelNHit", &MichelNHit);
    tMuon->Branch("MichelHitSlice", &MichelHits.slice);
    tMuon->Branch("MichelHitZ", &MichelHits.Z);
    tMuon->Branch("MichelHitChannel", &MichelHits.channel);
    tMuon->Branch("MichelHitTick", &MichelHits.tick);
    tMuon->Branch("MichelHitADC", &MichelHits.adc);

    tMuon->Branch("SphereNHit", &SphereNHit);
    tMuon->Branch("SphereHitSlice", &SphereHits.slice);
    tMuon->Branch("SphereHitZ", &SphereHits.Z);
    tMuon->Branch("SphereHitChannel", &SphereHits.channel);
    tMuon->Branch("SphereHitTick", &SphereHits.tick);
    tMuon->Branch("SphereHitADC", &SphereHits.adc);

    tMuon->Branch("SphereEnergy", &SphereEnergy);

    tMuon->Branch("SphereTruePositive", &SphereTruePositive);
    tMuon->Branch("SphereFalsePositive", &SphereFalsePositive);

    tMuon->Branch("SphereEnergyTruePositive", &SphereEnergyTruePositive);
    tMuon->Branch("SphereEnergyFalsePositive", &SphereEnergyFalsePositive);


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

        for (unsigned int w=0; w<asWire->Nwires(planeid); w++) {
            geo::WireID const wireid = geo::WireID{planeid, w};
            geo::WireGeo const wiregeo = asWire->Wire(wireid);
            mapChanZ[asWire->PlaneWireToChannel(wireid)] = wiregeo.GetStart().Z();

            if (iLogLevel >= kDetails) {
                std::cout << "ch:" << asWire->PlaneWireToChannel(wireid) << " Z:" << wiregeo.GetStart().Z() << " ";
                if (w % 8 == 7) std::cout << std::endl;
            }
        }
    }



    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fMichelTimeRadius = fMichelSpaceRadius / fDriftVelocity;
    fMichelTickRadius = fMichelTimeRadius / fSamplingRate;
    fMichelChannelRadius = int(fMichelSpaceRadius / fChannelPitch);

    if (iLogLevel >= kDetails) {
        std::cout << "fMichelSpaceRadius: " << fMichelSpaceRadius << " cm" << std::endl;
        std::cout << "fChannelPitch: " << fChannelPitch << " cm/channel" << std::endl;
        std::cout << "fMichelChannelRadius: " << fMichelChannelRadius << " channels" << std::endl; std::cout << "fMichelTimeRadius: " << fMichelTimeRadius << " µs" << std::endl;
        std::cout << "fSamplingRate: " << fSamplingRate << " µs/tick" << std::endl;
        std::cout << "fMichelTickRadius: " << fMichelTickRadius << " ticks" << std::endl;
    }


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Muchecks::Muchecks =======================================================" << "\033[0m" << std::endl;
}

void ana::Muchecks::analyze(art::Event const& e)
{

    if (iLogLevel >= kBasics) {
        std::cout << "\033[93m" << "Muchecks::analyze: Initialization evt#" << iEvent << "\r" << std::flush;
        std::cout << std::string(5,'\t') << " ======================================" << "\033[0m" << std::endl;
    }

    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fMichelTickRadius = fMichelTimeRadius / fSamplingRate;
    fMichelChannelRadius = int(fMichelSpaceRadius / fChannelPitch);

    auto const & vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    auto const & vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);



    resetEvent();

    for (recob::Hit const& hit : *vh_hit) {
        if (hit.View() != geo::kW) continue;

        EventNHit++;
        EventHits.push_back(GetHit(hit));
    } // end loop over hits


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Searching for muons' end ====================================" << "\033[0m" << std::endl;


    std::vector<ana::Hit> MuonsEnd;
    std::vector<size_t> MuonsTrkID;
    std::vector<simb::MCParticle const*> MichelsMCP;
    std::vector<simb::MCParticle const*> MuonsMCP;



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

        if (!fKeepNonDecaying && mcp->EndProcess() != "Decay") continue;

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());

        if (!Log(vp_hit_muon.size() > 0, kDetails, 1, "has hits...", "yes", "no")
        ) continue;


        // tacking min of ticks in lower volume and max of ticks in upper volume
        float TickUpMax = 0;
        float TickLowMin = detProp.ReadOutWindowSize();
        int ChanUpMax = -1;
        int ChanLowMin = -1;

        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            if (p_hit_muon->View() != geo::kW) continue;

            // upper volume (section 2 and 3) 
            if (p_hit_muon->Channel() > binChan[1].max) {
                if (p_hit_muon->PeakTime() > TickUpMax) {
                    TickUpMax = p_hit_muon->PeakTime();
                    ChanUpMax = p_hit_muon->Channel();
                }
            }

            // lower volume (section 0 and 1)
            if (p_hit_muon->Channel() < binChan[2].min) {
                if (p_hit_muon->PeakTime() < TickLowMin) {
                    TickLowMin = p_hit_muon->PeakTime();
                    ChanLowMin = p_hit_muon->Channel();
                }
            }
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


        // Muon End is 
        if (ChanLowMin == -1) {
            // if no hit in lower volume, muon end is in upper volume

            // if no hit at all, continue
            if (ChanUpMax == -1) continue;

            else {
                MuonEnd = GetHit(ChanUpMax, TickUpMax, 0);
            }
        } else {
            // if there is hit in lower volume, muon end is in lower volume
            MuonEnd = GetHit(ChanLowMin, TickLowMin, 0);
        }

        // getting track end point
        geo::Point_t end;
        if (IsUpright(*p_trk)) end = p_trk->End();
        else end = p_trk->Start();

        // three fiducial checks
        bool isInWindow = IsInsideWindow(MuonEnd.tick, fMichelTickRadius);
        bool isInVolumeZ = IsInsideZ(MuonEnd.Z, fMichelSpaceRadius);
        // dummy X to put it in upper volume, only check Y and Z
        bool isInVolumeYZ = IsInUpperVolume(50., end.Y(), end.Z(), fMichelSpaceRadius);

        Log(isInWindow, kDetails, 1, Form("is in window (±%.1f ticks)...", fMichelTickRadius), "yes", "no");
        Log(isInVolumeZ, kDetails, 1, Form("is in volumeZ (±%.1f cm)...", fMichelSpaceRadius), "yes", "no");
        Log(isInVolumeYZ, kDetails, 1, Form("is in volumeYZ (±%.1f cm)...", fMichelSpaceRadius), "yes", "no");

        if (!fKeepOutside && !isInWindow && !isInVolumeZ && !isInVolumeYZ) continue;

        // from now on, we keep the muon
    
        EventNMuon++;
        MuonsMCP.push_back(mcp);    
        MuonsEnd.push_back(MuonEnd);
        MuonsTrkID.push_back(p_trk->ID());

        // looking at daughters to find a michel
        int i_dau = mcp->NumberDaughters() - 1;
        if (iLogLevel >= kDetails) std::cout << "\t" << "looping over muon's " << mcp->NumberDaughters() << " daughters..." << std::endl;
        for (; i_dau >= 0; i_dau--) {
            if (iLogLevel >= kDetails) std::cout << "\t" << "dau#" << i_dau+1 << "\r" << std::flush;
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            // if (mcp_dau->Mother() != mcp->TrackId()) anomalies.push_back(Form("e%ld t%d µ%ld: TrackID: mcp_dau->Mother() (%d) != mcp->TrackId() (%d)", iEvent, p_trk->ID(), iMuon, mcp_dau->Mother(), mcp->TrackId()));

            if (!Log(mcp_dau, kDetails, 2, "id to mcp...", "done", "failed")
            ) continue;

            if (!Log(abs(mcp_dau->PdgCode()) == 11 && mcp_dau->Process() == "Decay", kDetails, 2, "is michel...", "yes", "no")
            ) continue;

            break;
        } // end loop over muon daughters
        if (i_dau == -1) {
            // loop over daughters ended without finding a michel
            if (iLogLevel >= kDetails) std::cout << "\t" << "no michel found" << std::endl;
            MichelsMCP.push_back(nullptr);
        }
        else {
            // loop over daughters ended with a michel at the index i_dau
            if (iLogLevel >= kDetails) std::cout << "\t" << "michel found at i_dau=" << i_dau << std::endl;
            MichelsMCP.push_back(pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau)));
        }
    } // end loop over tracks

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Michel Hits =================================================" << "\033[0m" << std::endl;

    std::vector<int> visited;
        
    
    // retrieving hits associated to michels MCParticles
    std::vector<std::vector<const recob::Hit*>> tpMichelHitAll;
    for (const simb::MCParticle * mcp_michel : MichelsMCP) {
        if (!mcp_michel) {
            tpMichelHitAll.push_back({});
            continue;
        }
        tpMichelHitAll.push_back(truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label()));
    }

    std::vector<int> tpSphereNHit(EventNMuon, 0);
    std::vector<ana::Hits> tpSphereHit(EventNMuon);
    std::vector<float> tpSphereEnergy(EventNMuon, 0);
    
    std::vector<size_t> tpSphereTruePositive(EventNMuon, 0);
    std::vector<size_t> tpSphereFalsePositive(EventNMuon, 0);

    std::vector<float> tpSphereEnergyTruePositive(EventNMuon, 0);
    std::vector<float> tpSphereEnergyFalsePositive(EventNMuon, 0);
    
    std::vector<float> tp_prev(EventNMuon, 0);

    // looping over hits and checking if they are in a sphere around the muons ends
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;

        // check if the hit is associated to a long track
        std::vector<art::Ptr<recob::Track>> vp_trk = fmp_hit2trk.at(p_hit.key());
        auto it = vp_trk.begin();
        while (it != vp_trk.end() && (*it)->Length() < fTrackLengthCut) it++;
        if (it != vp_trk.end()) continue;

        // looping over muons ends
        for (size_t m = 0; m<EventNMuon; m++) {
            if (!MichelsMCP[m]) continue;

            if (GetSlice(p_hit->Channel()) != MuonsEnd[m].slice) continue;

            float Z = (mapChanZ[p_hit->Channel()] - MuonsEnd[m].Z) / fMichelSpaceRadius;
            float T = (p_hit->PeakTime() - MuonsEnd[m].tick) / fMichelTickRadius;

            // ellipse around muon's end
            if (Z*Z + T*T > 1) continue;

            // checking if the hit is associated to the michel MCParticle
            if (std::find(tpMichelHitAll[m].begin(), tpMichelHitAll[m].end(), &*p_hit) != tpMichelHitAll[m].end()) {
                tpSphereTruePositive[m]++;
                tpSphereEnergyTruePositive[m] += p_hit->Integral();
            } else {
                tpSphereFalsePositive[m]++;
                tpSphereEnergyFalsePositive[m] += p_hit->Integral();
            }

            tpSphereNHit[m]++;
            tpSphereHit[m].push_back(GetHit(*p_hit));

            // avoiding double counting
            if (p_hit->ROISummedADC() != tp_prev[m]) {
                tpSphereEnergy[m] += p_hit->ROISummedADC();
                tp_prev[m] = p_hit->ROISummedADC();
            }
        }
    }
    tp_prev.clear();

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Filling Muon Tree ===========================================" << "\033[0m" << std::endl;

    // filling the muon tree
    for (size_t m=0; m<EventNMuon; m++) {

        resetMuon();

        EventiMuon.push_back(iMuon);

        MuonEnd = MuonsEnd[m];

        // getting track end point
        const art::Ptr<recob::Track> & p_trk = vp_trk[MuonsTrkID[m]];
        geo::Point_t end;
        if (IsUpright(*p_trk)) end = p_trk->End();
        else end = p_trk->Start();

        // three fiducial checks
        MuonEndIsInWindow = IsInsideWindow(MuonEnd.tick, fMichelTickRadius);
        MuonEndIsInVolumeZ = IsInsideZ(MuonEnd.Z, fMichelSpaceRadius);
        MuonEndIsInVolumeYZ = IsInUpperVolume(50., end.Y(), end.Z(), fMichelSpaceRadius);

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());

        // getting all muon hits
        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            if (p_hit_muon->View() != geo::kW) continue;

            MuonNHit++;
            MuonHits.push_back(GetHit(*p_hit_muon));
        }

        if (iLogLevel >= kInfos) std::cout << "\t" << "mu" << (MuonIsAnti ? "+" : "-") << "\033[93m" << " #" << EventNMuon << " (" << iMuon << ")" << "\033[0m" << " found" << std::endl;

        if (iLogLevel >= kInfos) std::cout << "\t" << "muon end @ Z:" << "\033[93m" << MuonEnd.Z << "\033[0m" 
            << " tick:" << "\033[93m" << MuonEnd.tick << "\033[0m"
            << " YZ:" << "\033[93m(" << end.Y() << "," << end.Z() << ")\033[0m" << std::endl;

        // some properties of the muon
        MuonIsAnti = int(MuonsMCP[m]->PdgCode() < 0);
        MuonTrackLength = p_trk->Length();
        MuonDoesDecay = MuonsMCP[m]->EndProcess() == "Decay";

        // getting track points
        for (size_t i_tpt=0; i_tpt<p_trk->NumberTrajectoryPoints(); i_tpt++) {
            if (!p_trk->HasValidPoint(i_tpt)) continue;
            MuonTrackPointX.push_back(p_trk->LocationAtPoint(i_tpt).X());
            MuonTrackPointY.push_back(p_trk->LocationAtPoint(i_tpt).Y());
            MuonTrackPointZ.push_back(p_trk->LocationAtPoint(i_tpt).Z());
            MuonNTrackPoint++;
        }

        simb::MCParticle const * mcp_michel = MichelsMCP[m];

        // checking if the muon has a michel
        if (!mcp_michel) {
            MuonHasMichel = kNoMichel;
            continue;
        }
        // checking if the michel MCParticle is inside the detector
        if (Log(IsInVolume(mcp_michel->Position(0), fMichelSpaceRadius), kInfos, 1, "michel is inside...", "yes", "no")
        ) {
            MuonHasMichel = kHasMichelInside; 
        } else { 
            MuonHasMichel = kHasMichelOutside;
        }

        recob::Track const * trk_michel = truthUtil.GetRecoTrackFromMCParticle(clockData, *mcp_michel, e, tag_trk.label());

        // getting the track length of the michel (if any)
        if (trk_michel) {
            MichelTrackLength = trk_michel->Length();
            if (iLogLevel >= kInfos) std::cout << "\t" << "michel track length: " << "\033[93m" << MichelTrackLength << "\033[0m" << std::endl;
        }

        // checking if the michel MCParticle has already been visited
        if (std::find(visited.begin(), visited.end(), mcp_michel->TrackId()) != visited.end()) anomalies.push_back(Form("e%ld µ%ld: TrackID: michel already visited (%d)", iEvent, iMuon, mcp_michel->TrackId()));
        visited.push_back(mcp_michel->TrackId());

        MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;
        if (iLogLevel >= kInfos) std::cout << "\t" << "michel true energy: " << "\033[93m" << MichelTrueEnergy << " MeV" << "\033[0m" << std::endl;

        // double dist = Dist(MuonsMCP[m]->EndPosition(), mcp_michel->Position(0));
        // if (dist > 0) anomalies.push_back(Form("e%ld t%d µ%ld: distance of %f", iEvent, p_trk->ID(), iMuon, dist)); 

        // saving anomalous michel energies
        if (MichelTrueEnergy > 100) anomalies.push_back(Form("e%ld µ%ld: Michel True Energy %.1f MeV", iEvent, iMuon, MichelTrueEnergy)); 

        // getting hits associated to the michel MCParticle
        float prev = 0;
        for (const recob::Hit* hit_michel : tpMichelHitAll[m]) {
            if (hit_michel->View() != geo::kW) continue;

            MichelNHit++;
            MichelHits.push_back(GetHit(*hit_michel));
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

        // getting all the hits selected within the michel sphere
        SphereNHit = tpSphereNHit[m];
        SphereHits = tpSphereHit[m];
        SphereEnergy = tpSphereEnergy[m] * fADCtoE;

        SphereTruePositive = tpSphereTruePositive[m];
        SphereFalsePositive = tpSphereFalsePositive[m];

        SphereEnergyTruePositive = tpSphereEnergyTruePositive[m] * fADCtoE;
        SphereEnergyFalsePositive = tpSphereEnergyFalsePositive[m] * fADCtoE;

        if (iLogLevel >= kInfos) {
            std::cout << "\t" << "michel sphere hits: " << "\033[93m" << SphereNHit << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere energy: " << "\033[93m" << SphereEnergy << " MeV" << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere true positive: " << "\033[93m" << SphereTruePositive << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere false positive: " << "\033[93m" << SphereFalsePositive << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere energy true positive: " << "\033[93m" << SphereEnergyTruePositive << " MeV" << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere energy false positive: " << "\033[93m" << SphereEnergyFalsePositive << " MeV" << "\033[0m" << std::endl;
        }

        tMuon->Fill();
        iMuon++;
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
    EventHits.clear();

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
    MuonHits.clear();

    MichelTrueEnergy = 0;
    MichelTrackLength = 0;

    MichelNHit = 0;
    MichelHits.clear();
    MichelHitEnergy = 0;

    SphereNHit = 0;
    SphereHits.clear();
    SphereEnergy = 0;

    SphereTruePositive = 0;
    SphereFalsePositive = 0;
    SphereEnergyTruePositive = 0;
    SphereEnergyFalsePositive = 0;

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


ana::Hit ana::Muchecks::GetHit(int c, float t, float a) {
    unsigned int s = GetSlice(c);
    float z = mapChanZ[c];
    return ana::Hit(s, z, c, t, a);
}
ana::Hit ana::Muchecks::GetHit(recob::Hit const& hit) {
    return ana::Hit(GetSlice(hit.Channel()), mapChanZ[hit.Channel()], hit.Channel(), hit.PeakTime(), hit.Integral());
}

unsigned int ana::Muchecks::GetSlice(int ch) {
    for (unsigned int tpc=0; tpc<asGeo->NTPC(); tpc++) {
        if (chanTPC[tpc].min <= ch and ch <= chanTPC[tpc].max) {
            unsigned int slice = tpc/4 * 2 + tpc % 2;
            return slice;
        }
    }
    std::cerr << "Error: channel " << ch << " not in any slice" << std::endl;
    return 8;
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
bool ana::Muchecks::IsInsideZ(float z, float eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpcid0{cryoid, 0}, tpcidN{cryoid, asGeo->NTPC(cryoid)-1};
    return z - eps > asGeo->TPC(tpcid0).MinZ() and z + eps < asGeo->TPC(tpcidN).MaxZ();
}


bool ana::Muchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}

DEFINE_ART_MODULE(ana::Muchecks)
