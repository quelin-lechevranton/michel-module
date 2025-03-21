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
        
        friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
            os << "sl:" << hit.slice << " Z:" << hit.Z << " ch:" << hit.channel << " tick:" << hit.tick << " ADC:" << hit.adc;
            return os;
        }
    };
    struct Hits {
        unsigned int N;
        std::vector<unsigned int> slice;
        std::vector<float> Z;
        std::vector<int> channel;
        std::vector<float> tick;
        std::vector<float> adc;

        Hits() : N(0), slice(), Z(), channel(), tick(), adc() {}
        void push_back(Hit const& hit) {
            N++;
            slice.push_back(hit.slice);
            Z.push_back(hit.Z);
            channel.push_back(hit.channel);
            tick.push_back(hit.tick);
            adc.push_back(hit.adc);
        }
        void clear() {
            N = 0;
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

    unsigned iEvent=0;

    ana::Hits EventHits;

    unsigned EventNMuon;
    std::vector<unsigned> EventiMuon; //[muon]

    // enum Tag { kMCP, kDaughters, kSphere, kTrue, kNTag };
    TTree* tMuon;
    std::vector<TBranch*> brMuon, brMichel;

    unsigned iMuon=0;

    int MuonIsAnti;
    std::string MuonEndProcess;
    int MuonHasMichel;

    float MuonTrackLength;

    unsigned MuonNTrackPoint;
    std::vector<float> MuonTrackPointX, MuonTrackPointY, MuonTrackPointZ;

    float MuonEndPointY, MuonEndPointZ;
    float MuonEndSpacePointY, MuonEndSpacePointZ;

    ana::Hit MuonEnd;
    int MuonEndIsInWindow;
    int MuonEndIsInVolumeZ;
    int MuonEndIsInVolumeYZ;

    ana::Hits MuonHits;

    float MichelTrackLength;
    float MichelTrueEnergy;
    float MichelHitEnergy;

    ana::Hits MichelHits;

    ana::Hits SphereHits;
    float SphereEnergy;

    unsigned SphereTruePositive;
    unsigned SphereFalsePositive;

    float SphereEnergyTruePositive;
    float SphereEnergyFalsePositive;

    ana::Hits NearbyHits;

    // Diagnostic Variables
    unsigned n_section = 4;
    ana::Binning binTick;
    std::vector<ana::Binning> binChan; //[section]
    std::vector<ana::Binning> chanTPC; //[tpc]

    std::map<int, float> mapChanZ;


    // Functions
    void resetEvent();
    void resetMuon();
    void resetMichel();
    bool Log(bool cond, int flag, int tab, std::string msg, std::string succ, std::string fail);

    // unsigned int GetSection(int ch);
    unsigned int GetSlice(raw::ChannelID_t ch);
    // ana::Hit GetHit(int ch, float t, float a);
    ana::Hit GetHit(recob::Hit const& hit);
    // geo::View_t GetPlane(raw::ChannelID_t ch, int sec);
    // geo::WireID GetWireID(geo::Point_t const& P, geo::View_t plane);
    // raw::ChannelID_t GetChannel(geo::Point_t const& P, geo::View_t plane);

    bool IsInUpperVolume(raw::ChannelID_t ch);
    bool IsInLowerVolume(raw::ChannelID_t ch);

    bool IsInUpperVolume(double x, double y, double z, float eps);
    bool IsInLowerVolume(double x, double y, double z, float eps);
    bool IsInVolume(double x, double y, double z, float eps);
    bool IsInVolume(float x, float y, float z, float eps);
    // bool IsInUpperVolume(float x, float y, float z, float eps);
    // bool IsInLowerVolume(float x, float y, float z, float eps);
    template<class Vec> bool IsInVolume(Vec const& V, float eps);
    // template<class Vec> bool IsInUpperVolume(Vec const& V, float eps);
    // template<class Vec> bool IsInLowerVolume(Vec const& V, float eps);

    bool IsInsideWindow(float tick, float eps);
    // bool IsInsideTPC(raw::ChannelID_t channel, int eps);
    bool IsInsideZ(float z, float eps);
    bool IsInsideYZ(geo::Point_t const& P, float eps);

    bool IsUpright(recob::Track const& T);
};



ana::Muchecks::Muchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),

    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("TrackLengthCut")), // in cm

    fKeepOutside(p.get<bool>("KeepOutside")),
    fKeepNonDecaying(p.get<bool>("KeepNonDecaying"))
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::Muchecks: ============================================================" << "\033[0m" << std::endl;
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
    tEvent->Branch("NHit", &EventHits.N);
    tEvent->Branch("HitSlice", &EventHits.slice);
    tEvent->Branch("HitZ", &EventHits.Z);
    tEvent->Branch("HitChannel", &EventHits.channel);
    tEvent->Branch("HitTick", &EventHits.tick);
    tEvent->Branch("HitADC", &EventHits.adc);

    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);

    tMuon = tfs->make<TTree>("Muon","");

    brMuon.push_back(tMuon->Branch("iEvent", &iEvent));
    brMuon.push_back(tMuon->Branch("iMuon", &iMuon));

    brMuon.push_back(tMuon->Branch("IsAnti", &MuonIsAnti));
    brMuon.push_back(tMuon->Branch("TrackLength", &MuonTrackLength));
    brMuon.push_back(tMuon->Branch("EndProcess", &MuonEndProcess));
    brMuon.push_back(tMuon->Branch("HasMichel", &MuonHasMichel));

    brMuon.push_back(tMuon->Branch("NHit", &MuonHits.N));
    brMuon.push_back(tMuon->Branch("HitSlice", &MuonHits.slice));
    brMuon.push_back(tMuon->Branch("HitZ", &MuonHits.Z));
    brMuon.push_back(tMuon->Branch("HitChannel", &MuonHits.channel));
    brMuon.push_back(tMuon->Branch("HitTick", &MuonHits.tick));
    brMuon.push_back(tMuon->Branch("HitADC", &MuonHits.adc));

    brMuon.push_back(tMuon->Branch("NTrackPoint", &MuonNTrackPoint));
    brMuon.push_back(tMuon->Branch("TrackPointX", &MuonTrackPointX));
    brMuon.push_back(tMuon->Branch("TrackPointY", &MuonTrackPointY));
    brMuon.push_back(tMuon->Branch("TrackPointZ", &MuonTrackPointZ));

    brMuon.push_back(tMuon->Branch("EndPointY", &MuonEndPointY));
    brMuon.push_back(tMuon->Branch("EndPointZ", &MuonEndPointZ));
    brMuon.push_back(tMuon->Branch("EndSpacePointY", &MuonEndSpacePointY));
    brMuon.push_back(tMuon->Branch("EndSpacePointZ", &MuonEndSpacePointZ));

    brMuon.push_back(tMuon->Branch("EndHitSlice", &MuonEnd.slice));
    brMuon.push_back(tMuon->Branch("EndHitZ", &MuonEnd.Z));
    brMuon.push_back(tMuon->Branch("EndHitChannel", &MuonEnd.channel));
    brMuon.push_back(tMuon->Branch("EndHitTick", &MuonEnd.tick));
    brMuon.push_back(tMuon->Branch("EndHitADC", &MuonEnd.adc));

    brMuon.push_back(tMuon->Branch("EndIsInWindow", &MuonEndIsInWindow));
    brMuon.push_back(tMuon->Branch("EndIsInVolumeZ", &MuonEndIsInVolumeZ));
    brMuon.push_back(tMuon->Branch("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ));


    brMichel.push_back(tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy));
    brMichel.push_back(tMuon->Branch("MichelTrackLength", &MichelTrackLength));
    brMichel.push_back(tMuon->Branch("MichelHitEnergy", &MichelHitEnergy));

    brMichel.push_back(tMuon->Branch("MichelNHit", &MichelHits.N));
    brMichel.push_back(tMuon->Branch("MichelHitSlice", &MichelHits.slice));
    brMichel.push_back(tMuon->Branch("MichelHitZ", &MichelHits.Z));
    brMichel.push_back(tMuon->Branch("MichelHitChannel", &MichelHits.channel));
    brMichel.push_back(tMuon->Branch("MichelHitTick", &MichelHits.tick));
    brMichel.push_back(tMuon->Branch("MichelHitADC", &MichelHits.adc));


    brMichel.push_back(tMuon->Branch("SphereNHit", &SphereHits.N));
    brMichel.push_back(tMuon->Branch("SphereHitSlice", &SphereHits.slice));
    brMichel.push_back(tMuon->Branch("SphereHitZ", &SphereHits.Z));
    brMichel.push_back(tMuon->Branch("SphereHitChannel", &SphereHits.channel));
    brMichel.push_back(tMuon->Branch("SphereHitTick", &SphereHits.tick));
    brMichel.push_back(tMuon->Branch("SphereHitADC", &SphereHits.adc));

    brMichel.push_back(tMuon->Branch("SphereEnergy", &SphereEnergy));

    brMichel.push_back(tMuon->Branch("SphereTruePositive", &SphereTruePositive));
    brMichel.push_back(tMuon->Branch("SphereFalsePositive", &SphereFalsePositive));

    brMichel.push_back(tMuon->Branch("SphereEnergyTruePositive", &SphereEnergyTruePositive));
    brMichel.push_back(tMuon->Branch("SphereEnergyFalsePositive", &SphereEnergyFalsePositive));

    brMichel.push_back(tMuon->Branch("NearbyNHit", &NearbyHits.N));
    brMichel.push_back(tMuon->Branch("NearbyHitSlice", &NearbyHits.slice));
    brMichel.push_back(tMuon->Branch("NearbyHitZ", &NearbyHits.Z));
    brMichel.push_back(tMuon->Branch("NearbyHitChannel", &NearbyHits.channel));
    brMichel.push_back(tMuon->Branch("NearbyHitTick", &NearbyHits.tick));
    brMichel.push_back(tMuon->Branch("NearbyHitADC", &NearbyHits.adc));

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

            // if (iLogLevel >= kDetails) {
            //     std::cout << "ch:" << asWire->PlaneWireToChannel(wireid) << " Z:" << wiregeo.GetStart().Z() << " ";
            //     if (w % 8 == 7) std::cout << std::endl;
            // }
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
        std::cout << "fDriftVelocity: " << fDriftVelocity << " cm/µs" << std::endl; 
        std::cout << "fSamplingRate: " << fSamplingRate << " µs/tick" << std::endl;
        std::cout << "fMichelTickRadius: " << fMichelTickRadius << " ticks" << std::endl;
    }


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Muchecks::Muchecks ======================================================" << "\033[0m" << std::endl;
}

void ana::Muchecks::analyze(art::Event const& e)
{

    if (iLogLevel >= kBasics) {
        std::cout << "\033[93m" << "Muchecks::analyze: Initialization evt#" << std::setw(5) << iEvent << " ====================================" << "\033[0m" << std::endl;
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
    art::FindManyP<recob::SpacePoint> fmp_hit2spt(vh_hit, e, tag_spt);


    resetEvent();

    for (recob::Hit const& hit : *vh_hit) {
        if (hit.View() != geo::kW) continue;

        EventHits.push_back(GetHit(hit));
    } // end loop over hits


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Analyzing Muons =============================================" << "\033[0m" << std::endl;


    std::vector<ana::Hit> MuonsEnd;
    std::vector<simb::MCParticle const*> MichelsMCP;
    std::vector<simb::MCParticle const*> MuonsMCP;
    std::vector<size_t> MuonTrackKey;


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
        float TickUpMax = 0, TickLowMin = detProp.ReadOutWindowSize();
        art::Ptr<recob::Hit> HitUpMax, HitLowMin;
        // recob::Hit const* HitLowMin = new recob::Hit(raw::TDCtick_t(0), raw::TDCtick_t(0), raw::InvalidChannelID, detProp.ReadOutWindowSize(), 0., 0., 0., 0., 0., 0., 0., 0., 0, 0, 0., 0, geo::View_t(0), geo::SigType_t(0), geo::WireID());

        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            if (p_hit_muon->View() != geo::kW) continue;

            if (IsInUpperVolume(p_hit_muon->Channel()) && p_hit_muon->PeakTime() > TickUpMax) {
                TickUpMax = p_hit_muon->PeakTime();
                HitUpMax = p_hit_muon;
            }

            if (IsInLowerVolume(p_hit_muon->Channel()) && p_hit_muon->PeakTime() < TickLowMin) {
                TickLowMin = p_hit_muon->PeakTime();
                HitLowMin = p_hit_muon;
            }
        }


        if (iLogLevel >= kDetails) {
            if (!HitUpMax) std::cout << "\t" << "\033[91m" << "no hit in upper volume" << "\033[0m" << std::endl;
            else {
                std::cout << "\t" << "up max channel: " << "\033[93m" << HitUpMax->Channel() << "\033[0m" << std::endl;
                std::cout << "\t" << "up max tick: " << "\033[93m" << HitUpMax->PeakTime() << "\033[0m" << std::endl;
            }
            std::cout << "\t" << "------------------------------------------------" << std::endl;
            if (!HitLowMin) std::cout << "\t" << "\033[91m" << "no hit in lower volume" << "\033[0m" << std::endl;
            else {
                std::cout << "\t" << "low min channel: " << "\033[93m" << HitLowMin->Channel() << "\033[0m" << std::endl;
                std::cout << "\t" << "low min tick: " << "\033[93m" << HitLowMin->PeakTime() << "\033[0m" << std::endl;
            }
            std::cout << "\t" << "------------------------------------------------" << std::endl;
        }

        resetMuon();

        // Muon End is 
        // if no hit in lower volume, muon end is in upper volume
        // if no hit at all, continue
        // if there is hit in lower volume, muon end is in lower volume
        std::vector<art::Ptr<recob::SpacePoint>> vp_spt;

        if (!HitLowMin) {   
            if (!HitUpMax) {
                continue;
            } else {
                MuonEnd = GetHit(*HitUpMax);
                vp_spt = fmp_hit2spt.at(HitUpMax.key());
            }
        } else {
            MuonEnd = GetHit(*HitLowMin); 
            vp_spt = fmp_hit2spt.at(HitLowMin.key());
        }

        // getting track end point
        geo::Point_t EndPoint;
        if (IsUpright(*p_trk)) EndPoint = p_trk->End();
        else EndPoint = p_trk->Start();

        MuonEndPointY = EndPoint.Y();
        MuonEndPointZ = EndPoint.Z();

        if (vp_spt.size()) {
            MuonEndSpacePointY = vp_spt.front()->position().Y();
            MuonEndSpacePointZ = vp_spt.front()->position().Z();
        }

        // three fiducial checks
        MuonEndIsInWindow = IsInsideWindow(MuonEnd.tick, fMichelTickRadius);
        MuonEndIsInVolumeZ = IsInsideZ(MuonEnd.Z, fMichelSpaceRadius);
        // dummy X to put it in upper volume, only check Y and Z
        MuonEndIsInVolumeYZ = IsInUpperVolume(150., EndPoint.Y(), EndPoint.Z(), fMichelSpaceRadius);

        Log(MuonEndIsInWindow, kDetails, 1, Form("is in window (±%.1f ticks)...", fMichelTickRadius), "yes", "no");
        Log(MuonEndIsInVolumeZ, kDetails, 1, Form("is in volumeZ (±%.1f cm)...", fMichelSpaceRadius), "yes", "no");
        Log(MuonEndIsInVolumeYZ, kDetails, 1, Form("is in volumeYZ (±%.1f cm)...", fMichelSpaceRadius), "yes", "no");

        if (!fKeepOutside && !(MuonEndIsInWindow && MuonEndIsInVolumeYZ)) continue;

        // from now on, we keep the muon
    
        EventNMuon++;
        EventiMuon.push_back(iMuon);
        MuonTrackKey.push_back(p_trk.key());

        // some properties of the muon
        MuonIsAnti = int(mcp->PdgCode() < 0);
        MuonTrackLength = p_trk->Length();
        MuonEndProcess = mcp->EndProcess();


        // getting all muon hits
        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            if (p_hit_muon->View() != geo::kW) continue;

            MuonHits.push_back(GetHit(*p_hit_muon));
        }
        // getting track points
        for (unsigned i_tpt=0; i_tpt<p_trk->NumberTrajectoryPoints(); i_tpt++) {
            if (!p_trk->HasValidPoint(i_tpt)) continue;

            MuonNTrackPoint++;
            MuonTrackPointX.push_back(p_trk->LocationAtPoint(i_tpt).X());
            MuonTrackPointY.push_back(p_trk->LocationAtPoint(i_tpt).Y());
            MuonTrackPointZ.push_back(p_trk->LocationAtPoint(i_tpt).Z());
        }

        if (iLogLevel >= kInfos) std::cout << "\t" << "mu" << (MuonIsAnti ? "+" : "-") << "\033[93m" << " #" << EventNMuon << " (" << iMuon << ")" << "\033[0m" << " found" << std::endl;

        if (iLogLevel >= kInfos) {
            std::cout << "\t" << "muon end hit: " << "\033[93m" << MuonEnd << "\033[0m" << std::endl;
            std::cout << "\t" << "muon end point: " << "YZ:" << "\033[93m(" << MuonEndPointY << "," << MuonEndPointZ << ")\033[0m" << std::endl;
        }


        MuonsMCP.push_back(mcp);    
        MuonsEnd.push_back(MuonEnd);


        // looking at daughters to find a michel
        bool has_numu=false, has_nue=false;
        simb::MCParticle const *mcp_michel = nullptr;
        if (iLogLevel >= kDetails) std::cout << "\t" << "looping over muon's " << mcp->NumberDaughters() << " daughters..." << std::endl;
        for (int i_dau=mcp->NumberDaughters()-3; i_dau<mcp->NumberDaughters(); i_dau++) {
            if (iLogLevel >= kDetails) std::cout << "\t" << "dau#" << i_dau+1 << "\r" << std::flush;
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            if (!Log(mcp_dau, kDetails, 2, "id to mcp...", "done", "failed")
            ) continue;

            if (Log(abs(mcp_dau->PdgCode()) == 14, kDetails, 2, "is muon neutrino...", "yes", "no")
            ) has_numu = true;
            if (Log(abs(mcp_dau->PdgCode()) == 12, kDetails, 2, "is electron neutrino...", "yes", "no")
            ) has_nue = true;
            if (Log(abs(mcp_dau->PdgCode()) == 11, kDetails, 2, "is electron...", "yes", "no")
            ) mcp_michel = mcp_dau;
        } // end loop over muon daughters
        if (mcp_michel && has_numu && has_nue) {
            if (iLogLevel >= kDetails) std::cout << "\t" << "michel found" << std::endl;

            if (Log(IsInVolume(mcp_michel->Position(0), fMichelSpaceRadius), kInfos, 1, "michel is inside...", "yes", "no"))
                MuonHasMichel = kHasMichelInside; 
            else
                MuonHasMichel = kHasMichelOutside;
        }
        else {
            if (iLogLevel >= kDetails) std::cout << "\t" << "no michel found" << std::endl;
            MuonHasMichel = kNoMichel;
            mcp_michel = nullptr;
        }
        MichelsMCP.push_back(mcp_michel);

        for (TBranch * b : brMuon) b->Fill();
        iMuon++;
    } // end loop over tracks

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Searching Michel Hits in Sphere =============================" << "\033[0m" << std::endl;

    std::vector<int> visited;

    std::vector<ana::Hits> tpNearbyHits(EventNMuon);
    
    std::vector<ana::Hits> tpSphereHits(EventNMuon);
    std::vector<float> tpSphereEnergy(EventNMuon, 0);
    
    std::vector<unsigned> tpSphereTruePositive(EventNMuon, 0);
    std::vector<unsigned> tpSphereFalsePositive(EventNMuon, 0);

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
        bool from_track = (it != vp_trk.end());

        // looping over muons ends
        for (unsigned m = 0; m<EventNMuon; m++) {
            if (GetSlice(p_hit->Channel()) != MuonsEnd[m].slice) continue;

            float Z = (mapChanZ[p_hit->Channel()] - MuonsEnd[m].Z) / fMichelSpaceRadius;
            float T = (p_hit->PeakTime() - MuonsEnd[m].tick) / fMichelTickRadius;
            float sphere = Z*Z + T*T;

            bool from_this_muon = from_track ? (MuonTrackKey[m] == it->key()) : false;

            // not from on other track
            if (from_track && !from_this_muon) continue;
            // nearby muon's end
            if (sphere > 4) continue;

            tpNearbyHits[m].push_back(GetHit(*p_hit));

            // only integrate around muon with a michel electron
            if (!MichelsMCP[m]) continue;

            // not from any track
            if (from_track) continue;
            // sphere around muon's end
            if (sphere > 1) continue;

            // checking if the hit is associated to the michel MCParticle
            std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *MichelsMCP[m], e, tag_hit.label());
            if (std::find(v_hit_michel.begin(), v_hit_michel.end(), &*p_hit) != v_hit_michel.end()) {
                tpSphereTruePositive[m]++;
                tpSphereEnergyTruePositive[m] += p_hit->Integral();
            } else {
                tpSphereFalsePositive[m]++;
                tpSphereEnergyFalsePositive[m] += p_hit->Integral();
            }

            tpSphereHits[m].push_back(GetHit(*p_hit));

            // avoiding double counting
            if (p_hit->ROISummedADC() != tp_prev[m]) {
                tpSphereEnergy[m] += p_hit->ROISummedADC();
                tp_prev[m] = p_hit->ROISummedADC();
            }
        }
    }
    tp_prev.clear();

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Filling Michel Branches =====================================" << "\033[0m" << std::endl;

    // filling the michel tree
    for (unsigned m=0; m<EventNMuon; m++) {

        resetMichel();

        NearbyHits = tpNearbyHits[m];

        simb::MCParticle const * mcp_michel = MichelsMCP[m];

        if (!mcp_michel) {
            for (TBranch * b : brMichel) b->Fill();
            continue;
        }

        if (iLogLevel >= kInfos) std::cout << "mu#" << m+1 << "\r" << std::flush;
        recob::Track const * trk_michel = truthUtil.GetRecoTrackFromMCParticle(clockData, *mcp_michel, e, tag_trk.label());

        // getting the track length of the michel (if any)
        if (trk_michel) {
            MichelTrackLength = trk_michel->Length();
            if (iLogLevel >= kInfos) std::cout << "\t" << "michel track length: " << "\033[93m" << MichelTrackLength << "\033[0m" << std::endl;
        }

        // checking if the michel MCParticle has already been visited
        if (std::find(visited.begin(), visited.end(), mcp_michel->TrackId()) != visited.end()) anomalies.push_back(Form("e%u µ%u: TrackID: michel already visited (%d)", iEvent, m, mcp_michel->TrackId()));
        visited.push_back(mcp_michel->TrackId());

        MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;
        if (iLogLevel >= kInfos) std::cout << "\t" << "michel true energy: " << "\033[93m" << MichelTrueEnergy << " MeV" << "\033[0m" << std::endl;

        // double dist = Dist(MuonsMCP[m]->EndPosition(), mcp_michel->Position(0));
        // if (dist > 0) anomalies.push_back(Form("e%u t%d µ%u: distance of %f", iEvent, p_trk->ID(), iMuon, dist)); 

        // saving anomalous michel energies
        if (MichelTrueEnergy > 100) anomalies.push_back(Form("e%u µ%u: Michel True Energy %.1f MeV", iEvent, m, MichelTrueEnergy)); 

        // getting hits associated to the michel MCParticle
        std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());
        float prev = 0;
        for (const recob::Hit* hit_michel : v_hit_michel) {
            if (hit_michel->View() != geo::kW) continue;

            MichelHits.push_back(GetHit(*hit_michel));

            if (hit_michel->ROISummedADC() != prev) {
                MichelHitEnergy += hit_michel->ROISummedADC();
                prev = hit_michel->ROISummedADC();
            }
        }
        MichelHitEnergy *= fADCtoE;

        if (iLogLevel >= kInfos) {
            if (MichelHits.N == 0) std::cout << "\t" << "\033[91m" << "no hit in michel" << "\033[0m" << std::endl;
            else {
                std::cout << "\t" << "michel hits: " << "\033[93m" << MichelHits.N << "\033[0m" << std::endl;
                std::cout << "\t" << "michel hit energy: " << "\033[93m" << MichelHitEnergy << "\033[0m" << std::endl;
            }
            std::cout << "\t" << "michel process: " << "\033[93m" << mcp_michel->Process() << "\033[0m" << std::endl;
            std::cout << "\t" << "michel end process: " << "\033[93m" << mcp_michel->EndProcess() << "\033[0m" << std::endl;
        }

        // getting all the hits selected within the michel sphere
        SphereHits = tpSphereHits[m];
        SphereEnergy = tpSphereEnergy[m] * fADCtoE;

        SphereTruePositive = tpSphereTruePositive[m];
        SphereFalsePositive = tpSphereFalsePositive[m];

        SphereEnergyTruePositive = tpSphereEnergyTruePositive[m] * fADCtoE;
        SphereEnergyFalsePositive = tpSphereEnergyFalsePositive[m] * fADCtoE;

        // prev = 0;
        // float ene = std::accumulate(SphereHits.adc.begin(), SphereHits.adc.end(), 0, [&prev](float sum, float adc) {
        //     if (adc != prev) {
        //         prev = adc;
        //         return sum + adc;
        //     }
        //     return sum;
        // }) * fADCtoE;

        if (iLogLevel >= kInfos) {
            std::cout << "\t" << "michel sphere hits: " << "\033[93m" << SphereHits.N << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere energy: " << "\033[93m" << SphereEnergy << " MeV" << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere true positive: " << "\033[93m" << SphereTruePositive << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere false positive: " << "\033[93m" << SphereFalsePositive << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere energy true positive: " << "\033[93m" << SphereEnergyTruePositive << " MeV" << "\033[0m" << std::endl;
            std::cout << "\t" << "michel sphere energy false positive: " << "\033[93m" << SphereEnergyFalsePositive << " MeV" << "\033[0m" << std::endl;
        }

        // tMuon->Fill();
        // iMuon++;
        for (TBranch * b : brMichel) b->Fill();
    } // end loop over michels

    tMuon->SetEntries(brMuon.front()->GetEntries());
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

    EventHits.clear();

    EventNMuon = 0;
    EventiMuon.clear();

    if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
}
void ana::Muchecks::resetMuon() {
    if (iLogLevel >= kDetails) std::cout << "\t" << "resetting muon branches...";

    MuonIsAnti = 0;
    MuonTrackLength = 0;
    MuonHasMichel = 0;

    MuonHits.clear();

    MuonEndPointY = 0;
    MuonEndPointZ = 0;

    MuonEndSpacePointY = 0;
    MuonEndSpacePointZ = 0;

    MuonNTrackPoint = 0;
    MuonTrackPointX.clear();
    MuonTrackPointY.clear();
    MuonTrackPointZ.clear();

    MuonEnd = ana::Hit();

    MuonEndIsInWindow = false;
    MuonEndIsInVolumeZ = false;
    MuonEndIsInVolumeYZ = false;


    if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
}
void ana::Muchecks::resetMichel() {
    if (iLogLevel >= kDetails) std::cout << "\t" << "resetting michel branches...";

    MichelTrueEnergy = 0;
    MichelTrackLength = 0;

    MichelHits.clear();
    MichelHitEnergy = 0;

    SphereHits.clear();
    SphereEnergy = 0;

    NearbyHits.clear();

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


// ana::Hit ana::Muchecks::GetHit(int c, float t, float a) {
//     unsigned int s = GetSlice(c);
//     float z = mapChanZ[c];
//     return ana::Hit(s, z, c, t, a);
// }
ana::Hit ana::Muchecks::GetHit(recob::Hit const& hit) {
    return ana::Hit(
        GetSlice(hit.Channel()), 
        mapChanZ[hit.Channel()], 
        hit.Channel(), 
        hit.PeakTime(), 
        hit.Integral());
}
unsigned int ana::Muchecks::GetSlice(raw::ChannelID_t ch) {
    std::vector<geo::WireID> wids = asWire->ChannelToWire(ch);
    if (wids.size() == 0) return 8;
    geo::TPCID::TPCID_t tpc = wids.at(0).TPC;
    if (tpc == geo::TPCID::InvalidID) return 8;
    return (unsigned int) (tpc/4 * 2 + tpc % 2);
}
// unsigned int ana::Muchecks::GetSection(int ch) {
//     return (unsigned int) (4.*ch / asWire->Nchannels());
// }
// geo::View_t ana::Muchecks::GetPlane(raw::ChannelID_t ch, int sec) {
//     return geo::View_t(12.*ch / asWire->Nchannels() - 3*sec);
// }
// geo::WireID ana::Muchecks::GetWireID(geo::Point_t const& P, geo::View_t plane) {
//     geo::TPCID tpcid = asGeo->FindTPCAtPosition(P);
//     if (!tpcid.isValid) return geo::WireID();
//     geo::PlaneGeo const& planegeo = asWire->Plane(tpcid,plane);
//     geo::WireID wireid;
//     try {
//         wireid = planegeo.NearestWireID(P);
//     } catch (geo::InvalidWireError const& e) {
//         return e.suggestedWireID();
//     }
//     return wireid;
// }
// raw::ChannelID_t ana::Muchecks::GetChannel(geo::Point_t const& P, geo::View_t plane) {
//     geo::WireID wireid = GetWireID(P, plane);
//     if (!wireid.isValid) return raw::InvalidChannelID;
//     return asWire->PlaneWireToChannel(wireid);
// }


bool ana::Muchecks::IsInUpperVolume(raw::ChannelID_t ch) {
    if (ch == raw::InvalidChannelID) return false;
    return ch > binChan[1].max;
}
bool ana::Muchecks::IsInLowerVolume(raw::ChannelID_t ch) {
    if (ch == raw::InvalidChannelID) return false;
    return ch < binChan[2].min;
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
    return IsInUpperVolume(x, y, z, eps) or IsInLowerVolume(x, y, z, eps);
}
// bool ana::Muchecks::IsInLowerVolume(float x, float y, float z, float eps) {
//     return IsInLowerVolume(double(x), double(y), double(z), eps);
//     return IsInUpperVolume(double(x), double(y), double(z), eps);
// }
bool ana::Muchecks::IsInVolume(float x, float y, float z, float eps) {
    return IsInVolume(double(x), double(y), double(z), double(eps));
}
// template<class Vec>
// bool ana::Muchecks::IsInUpperVolume(Vec const& V, float eps) {
//     return IsInUpperVolume(V.X(), V.Y(), V.Z(), eps);
// }
// template<class Vec>
// bool ana::Muchecks::IsInLowerVolume(Vec const& V, float eps) {
//     return IsInLowerVolume(V.X(), V.Y(), V.Z(), eps);
// }
template<class Vec>
bool ana::Muchecks::IsInVolume(Vec const& V, float eps) {
    return IsInVolume(V.X(), V.Y(), V.Z(), eps);
}

bool ana::Muchecks::IsInsideWindow(float tick, float eps) {
    return binTick.min + eps < tick and tick < binTick.max - eps;
}
// bool ana::Muchecks::IsInsideTPC(raw::ChannelID_t ch, int eps) {
//     if (ch == raw::InvalidChannelID) return false;
//     bool isIn = false;
//     for (unsigned int tpc=0; tpc<asGeo->NTPC(); tpc++) {
//         isIn = chanTPC[tpc].min + eps < ch and ch < chanTPC[tpc].max - eps;
//         if (isIn) break;
//     }
//     return isIn;
// }
bool ana::Muchecks::IsInsideZ(float z, float eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpcid0{cryoid, 0}, tpcidN{cryoid, asGeo->NTPC(cryoid)-1};
    return z - eps > asGeo->TPC(tpcid0).MinZ() and z + eps < asGeo->TPC(tpcidN).MaxZ();
}
bool ana::Muchecks::IsInsideYZ(geo::Point_t const& P, float eps) {
    return IsInUpperVolume(50., P.Y(), P.Z(), eps);
}

bool ana::Muchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}

DEFINE_ART_MODULE(ana::Muchecks)
