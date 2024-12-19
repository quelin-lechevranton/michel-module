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
#include "lardata/ArtDataHelper/TrackUtils.h"

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
    float fChannelPitch = 0.5; // cm/channel
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
    
    // Output Variables


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

    size_t iMuon=0;

    int MuonIsAnti;
    int MuonDoesDecay;
    int MuonHasMichel;

    float MuonTrackLength;
    // int MuonIsUpright;

    size_t MuonNTrackPoint;
    std::vector<float> MuonTrackPointX, MuonTrackPointY, MuonTrackPointZ;

    geo::Point_t MuonEndPoint;
    geo::Point_t MuonStartPoint;

    // int MuonEndTrackInside;
    // int MuonStartTrackInside;

    // float MuonEndTickInside;
    int MuonEndChannel;
    float MuonEndTick;

    size_t MuonNHit;
    std::vector<int> MuonHitChannel; //[hit]
    std::vector<float> MuonHitTick; //[hit]
    std::vector<float> MuonHitADC; //[hit]

    // float MuonTickLowMin, MuonTickLowMax;
    // float MuonTickUpMin, MuonTickUpMax;



    // Diagnostic Variables
    size_t n_section = 4, n_plane = 3;
    // size_t iNEventDiagnostic;
    ana::Binning binTick;
    std::vector<ana::Binning> binChan; //[section]
    // std::vector<std::vector<TH2F*>> th2Hit; //[event][section]

    // Functions
    void resetEvent();
    void resetMuon();
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

    bool IsUpright(recob::Track const& T);

    // int iCorr(int n, int i, int j);
};


ana::Muchecks::Muchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),

    fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("TrackLengthCut")) // in cm

    // iNEventDiagnostic(p.get<int>("NEventDiagnostic"))
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::Muchecks: =============================================================" << "\033[0m" << std::endl;
    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

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
    // tMuon->Branch("IsUpright", &MuonIsUpright);

    tMuon->Branch("NTrackPoint", &MuonNTrackPoint);
    tMuon->Branch("TrackPointX", &MuonTrackPointX);
    tMuon->Branch("TrackPointY", &MuonTrackPointY);
    tMuon->Branch("TrackPointZ", &MuonTrackPointZ);


    tMuon->Branch("EndPoint", &MuonEndPoint);
    tMuon->Branch("StartPoint", &MuonStartPoint);

    // tMuon->Branch("EndTrackInside", &MuonEndTrackInside);
    // tMuon->Branch("StartTrackInside", &MuonStartTrackInside);

    // tMuon->Branch("EndTickInside", &MuonEndTickInside);

    tMuon->Branch("EndChannel", &MuonEndChannel);
    tMuon->Branch("EndTick", &MuonEndTick);

    tMuon->Branch("NHit", &MuonNHit);
    tMuon->Branch("HitChannel", &MuonHitChannel);
    tMuon->Branch("HitTick", &MuonHitTick);
    tMuon->Branch("HitADC", &MuonHitADC);

    // tMuon->Branch("TickLowMin", &MuonTickLowMin);
    // tMuon->Branch("TickLowMax", &MuonTickLowMax);
    // tMuon->Branch("TickUpMin", &MuonTickUpMin);
    // tMuon->Branch("TickUpMax", &MuonTickUpMax);



    // Diagnostic Variables
    binTick.n = detProp.ReadOutWindowSize()/4;
    binTick.min = 0;
    binTick.max = detProp.ReadOutWindowSize();

    int pp = -1; // previous plane
    size_t s = 0; // section
    binChan = std::vector<ana::Binning>(n_section);
    for (unsigned int c=0; c<asWire->Nchannels(); c++) {
        int p = asWire->View(raw::ChannelID_t(c));
        if (p == pp) continue;
        if (p == geo::kW) {
            binChan[s].min = c;
        } else if (pp == geo::kW) {
            binChan[s].max = c-1;
            binChan[s].n = binChan[s].max - binChan[s].min + 1;
            s++;
            if (s == n_section) {
                std::cerr << "\033[91m" << "ERROR: More than " << n_section << " sections found in the wire plane" << "\033[0m" << std::endl;
                break;
            }
        }
        pp = p;
    }
    binChan[n_section-1].max = asWire->Nchannels()-1;
    binChan[n_section-1].n = binChan[n_section-1].max - binChan[n_section-1].min + 1;


    if (iLogLevel >= kInfos) {
        for (size_t s=0; s<n_section; s++) {
            std::cout << "section#" << s << " plane#" << geo::kW << std::endl;
            std::cout << "\t" << "Tick: " << binTick.min << " -> " << binTick.max << " (" << binTick.n << ")" << std::endl;
            std::cout << "\t" << "Chan: " << binChan[s].min << " -> " << binChan[s].max << " (" << binChan[s].n << ")" << std::endl;
        }
    }

    // th2Hit = std::vector<std::vector<TH2F*>>(iNEventDiagnostic, std::vector<TH2F*>(n_section));
    // for (size_t e=0; e<iNEventDiagnostic; e++) {
    //     for (size_t s=0; s<n_section; s++) {
    //         th2Hit[e][s] = new TH2F(
    //             Form("th2Hit_%ld_%ld",e,s), 
    //             "",
    //             binTick.n, binTick.min, binTick.max,
    //             binChan[s].n, binChan[s].min, binChan[s].max
    //         );
    //     }

    //     for (int t=0; t<kTrue; t++) {
    //         tgMichelHit[e][t] = new TGraph();
    //     }
    // }


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

    auto const vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);
        
    auto const & vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);

    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;

        EventNHit++;
        EventHitChannel.push_back(p_hit->Channel());
        EventHitTick.push_back(p_hit->PeakTime());
        EventHitADC.push_back(p_hit->Integral());
    } // end loop over hits
    

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Muchecks::analyze: Muon Ends ===================================================" << "\033[0m" << std::endl;

    if (iLogLevel >= kInfos) std::cout << "looping over " << vp_trk.size() << " tracks..." << std::endl;
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        if (iLogLevel >= kInfos) std::cout << "trk#" << p_trk->ID()+1 << "\r" << std::flush;

        if (!Log(p_trk->Length() > fTrackLengthCut, kDetails, 1, "is long enough...", "yes", "no")
        ) continue;

        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());

        // simb::MCParticle const * mcp2 = GetMCParticle(clockData, p_trk, fmp_trk2hit);
        // if (mcp->TrackId() != mcp2->TrackId()) std::cout << "\033[91m" << "ERROR: mcp->TrackId() != mcp2->TrackId()" << "\033[0m" << " delta: " << (int) mcp2->TrackId() - (int) mcp->TrackId() << std::endl;

        if (!Log(mcp, kDetails, 1, "trk to mcp...", "done", "failed")
        ) continue;    

        if (!Log(abs(mcp->PdgCode()) == 13, kDetails, 1, "is muon...", "yes", "no")
        ) continue;

        // if (Log(mcp->EndProcess() == "Decay", kDetails, 1, "is decaying...", "yes", "no")
        // ) continue;



        // if (iLogLevel >= kInfos) std::cout << "\t" << p_trk->Start() << " -> " << p_trk->End();
        // if (IsUpright(*p_trk)) {
        //     MuonIsUpright = 1;

        //     MuonEndPoint = p_trk->End();
        //     MuonStartPoint = p_trk->Start();

        //     if (iLogLevel >= kInfos) std::cout << " " << "\033[93m" << "upright" << "\033[0m";
        // } else {
        //     MuonIsUpright = 0;

        //     MuonEndPoint = p_trk->Start();
        //     MuonStartPoint = p_trk->End();

        //     if (iLogLevel >= kInfos) std::cout << " " << "\033[93m" << "upside down" << "\033[0m";
        // }
        // if (iLogLevel >= kInfos) std::cout << std::endl;

        // if (IsInUpperVolume(MuonStartPoint, fMichelSpaceRadius)
        // ) {
        //     MuonStartTrackInside = kInUpper;
        // } else if (IsInLowerVolume(MuonStartPoint, fMichelSpaceRadius)
        // ) {
        //     MuonStartTrackInside = kInLower;
        // } else {
        //     MuonStartTrackInside = kOutside;
        // }

        // if (IsInUpperVolume(MuonEndPoint, fMichelSpaceRadius)
        // ) {
        //     MuonEndTrackInside = kInUpper;
        // } else if (IsInLowerVolume(MuonEndPoint, fMichelSpaceRadius)
        // ) {
        //     MuonEndTrackInside = kInLower;
        // } else {
        //     MuonEndTrackInside = kOutside;
        // }
        // if (MuonEndTrackInside == kOutside and MuonStartTrackInside == kOutside) {
        //     continue;
        // }



        std::vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk2hit.at(p_trk.key());

        if (!Log(vp_hit.size() > 0, kDetails, 1, "has hits...", "yes", "no")
        ) continue;

        float TickUpMax = 0;
        float TickLowMin = detProp.ReadOutWindowSize();
        int ChanUpMax = -1;
        int ChanLowMin = -1;

        for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;

            if (p_hit->Channel() > binChan[1].max) {
                if (p_hit->PeakTime() > TickUpMax) {
                    TickUpMax = p_hit->PeakTime();
                    ChanUpMax = p_hit->Channel();
                }
            }
            if (p_hit->Channel() < binChan[2].min) {
                if (p_hit->PeakTime() < TickLowMin) {
                    TickLowMin = p_hit->PeakTime();
                    ChanLowMin = p_hit->Channel();
                }
            }

            MuonNHit++;
            MuonHitChannel.push_back(p_hit->Channel());
            MuonHitTick.push_back(p_hit->PeakTime());
            MuonHitADC.push_back(p_hit->Integral());
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
    
        if (!Log(MuonEndTick > binTick.min + fMichelTimeRadius / fSamplingRate and
                MuonEndTick < binTick.max - fMichelTimeRadius / fSamplingRate, kDetails, 1, 
                Form("is in time window (%f -> %f)...", binTick.min + fMichelTimeRadius / fSamplingRate, binTick.max - fMichelTimeRadius / fSamplingRate),
                "yes", "no")
        ) continue;

        // if (iLogLevel >= kInfos) {
        //     std::cout << "\t" << "track end channel: " << "\033[93m" << MuonEndChannel << "\033[0m" << std::endl;
        //     std::cout << "\t" << "------------------------------------------------" << std::endl;
        //     std::cout << "\t" << "chan dist min: " << "\033[93m" << MuonEndChannel + ChanDistMin << "\033[0m" << std::endl;
        //     std::cout << "\t" << "tick dist min: " << "\033[93m" << TickDistMin << "\033[0m" << std::endl;
        //     if (iUpMin * iUpMax < 0) std::cerr << "\033[91m" << "ERROR: iUpMin * iUpMax < 0" << "\033[0m" << std::endl;
        //     else if (iUpMin == -1) {
        //         std::cout << "\t" << "\033[91m" << "no hit in upper volume" << "\033[0m" << std::endl;
        //     } else {
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "up min channel: " << "\033[93m" << MuonHitChannel[iUpMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "up min tick: " << "\033[93m" << MuonHitTick[iUpMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "up max channel: " << "\033[93m" << MuonHitChannel[iUpMax] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "up max tick: " << "\033[93m" << MuonHitTick[iUpMax] << "\033[0m" << std::endl;
        //     }
        //     if (iLowMin * iLowMax < 0) std::cerr << "\033[91m" << "ERROR: iLowMin * iLowMax < 0" << "\033[0m" << std::endl;
        //     else if (iLowMin == -1) {
        //         std::cout << "\t" << "\033[91m" << "no hit in lower volume" << "\033[0m" << std::endl;
        //     } else {
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "low min channel: " << "\033[93m" << MuonHitChannel[iLowMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "low min tick: " << "\033[93m" << MuonHitTick[iLowMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "low max channel: " << "\033[93m" << MuonHitChannel[iLowMax] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "low max tick: " << "\033[93m" << MuonHitTick[iLowMax] << "\033[0m" << std::endl;
        //     }
        // }
        // MuonTickLowMin = TickLowMin;
        // MuonTickLowMax = TickLowMax;
        // MuonTickUpMin = TickUpMin;
        // MuonTickUpMax = TickUpMax;


        // if (MuonHasMichel != kNoMichel) {
        //     simb::MCParticle const * mcp_michel = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));
        //     if (Log(mcp_michel, kInfos, 1, "id to mcp...", "done", "failed")
        //     ) {
        //         std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());

        //         size_t n_hit_michel = 0;
        //         int bary_channel_michel = 0;
        //         float bary_tick_michel = 0;

        //         for (const recob::Hit* hit_michel : v_hit_michel) {
        //             if (hit_michel->View() != geo::kW) continue;


        //             std::cout << "bonjour" << std::endl;

        //             n_hit_michel++;
        //             bary_channel_michel += hit_michel->Channel();
        //             bary_tick_michel += hit_michel->PeakTime();
        //         }

        //         // bary_channel_michel /= n_hit_michel;
        //         // bary_tick_michel /= n_hit_michel;

        //         if (iLogLevel >= kInfos) {
        //             std::cout << "\t" << "------------------------------------------------" << std::endl;
        //             std::cout << "\t" << "------------------------------------------------" << std::endl;
        //             if (n_hit_michel == 0) std::cout << "\t" << "\033[91m" << "no hit in michel" << "\033[0m" << std::endl;
        //             else {
        //                 std::cout << "\t" << "michel hits: " << "\033[93m" << n_hit_michel << "\033[0m" << std::endl;
        //                 std::cout << "\t" << "michel barycenter channel: " << "\033[93m" << bary_channel_michel << "\033[0m" << std::endl;
        //                 std::cout << "\t" << "michel barycenter tick: " << "\033[93m" << bary_tick_michel << "\033[0m" << std::endl;
        //             }
        //         }
        //     }
        // }


        // Analyzing muon

        resetMuon();

        MuonIsAnti = int(mcp->PdgCode() < 0);
        MuonTrackLength = p_trk->Length();
        MuonDoesDecay = mcp->Process() == "Decay";

        if (iLogLevel >= kInfos) std::cout << "\t" << "mu" << (MuonIsAnti ? "+" : "-") << "\033[93m" << " #" << EventNMuon << " (" << iMuon << ")" << "\033[0m" << " found" << std::endl;

        if (iLogLevel >= kInfos) std::cout << "\t" << "muon end @ channel:" << "\033[93m" << MuonEndChannel << "\033[0m" << " tick:" << "\033[93m" << MuonEndTick << "\033[0m" << std::endl;

        EventNMuon++;
        EventiMuon.push_back(iMuon++);
        // if (iLogLevel >= kInfos) {
        //     std::cout << "\t" << "start in: ";
        //     if (MuonStartTrackInside == kOutside) std::cout << "\033[91m" << "outside" << "\033[0m" << std::endl;
        //     else std::cout << "\033[93m" << (MuonStartTrackInside == kInUpper ? "upper" : "lower") << "\033[0m" << std::endl;
        //     std::cout << "\t" << "end in: ";
        //     if (MuonEndTrackInside == kOutside) std::cout << "\033[91m" << "outside" << "\033[0m" << std::endl;
        //     else std::cout << "\033[93m" << (MuonEndTrackInside == kInUpper ? "upper" : "lower") << "\033[0m" << std::endl;
        // }



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

            if (!Log(mcp_dau, kDetails, 2, "id to mcp...", "done", "failed")
            ) continue;

            if (!Log(abs(mcp_dau->PdgCode()) == 11 || mcp_dau->Process() != "Decay", kDetails, 2, "is michel...", "yes", "no")
            ) continue;

            if (!Log(IsInVolume(mcp_dau->Position(0), fMichelSpaceRadius), kDetails, 2, "is inside...", "yes", "no")
            ) {
                MuonHasMichel = kHasMichelInside;
            } else {
                MuonHasMichel = kHasMichelOutside;
            }
            if (iLogLevel >= kInfos) std::cout << "\t" << "michel found @i_dau:" << i_dau << "\033[0m" << std::endl;
            break;
        } // end loop over muon daughters
        if (i_dau == -1) MuonHasMichel = kNoMichel;

        // raw::ChannelID_t c = GetChannel(MuonEndPoint, geo::kW);
        // if (Log(c != raw::InvalidChannelID, kDetails, 1, "valid muon end channel...", "yes", "no")
        // ) MuonEndChannel = int(c);
        // else MuonEndChannel = -1;


        // std::vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk2hit.at(p_trk.key());

        // // std::vector<float> tpTick;
        // // std::vector<int> tpChanDist;
        // // std::vector<std::vector<float>> tpMuonHitVars(3); //[hitVar][hit]

        // float TickUpMin = detProp.ReadOutWindowSize(), TickUpMax = 0;
        // float TickLowMin = detProp.ReadOutWindowSize(), TickLowMax = 0;
        // int iUpMin=-1, iUpMax=-1, iLowMin=-1, iLowMax=-1;

        // int ChanDistMin = asWire->Nchannels();
        // float TickDistMin = 0;

        // int i_hit = 0;
        // bool error = false;
        // for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        //     if (p_hit->View() != geo::kW) continue;

        //     if (!error && MuonEndTrackInside == kInUpper) {
        //         if (p_hit->Channel() < binChan[2].min) {
        //             error = true;
        //             std::cout << "\033[91m" << "ERROR: MuonEndTrackInside == kInUpper and hits in lower volume around chan: " << p_hit->Channel() << "\033[0m" << std::endl;
        //         }
        //     }

        //     if (p_hit->Channel() > binChan[1].max) {
        //         if (p_hit->PeakTime() < TickUpMin) {
        //             TickUpMin = p_hit->PeakTime();
        //             iUpMin = i_hit;
        //         }
        //         if (p_hit->PeakTime() > TickUpMax) {
        //             TickUpMax = p_hit->PeakTime();
        //             iUpMax = i_hit;
        //         }
        //     }
        //     if (p_hit->Channel() < binChan[2].min) {
        //         if (p_hit->PeakTime() < TickLowMin) {
        //             TickLowMin = p_hit->PeakTime();
        //             iLowMin = i_hit;
        //         }
        //         if (p_hit->PeakTime() > TickLowMax) {
        //             TickLowMax = p_hit->PeakTime();
        //             iLowMax = i_hit;
        //         }
        //     }

        //     int dist = p_hit->Channel() - MuonEndChannel;
        //     if (abs(dist) < abs(ChanDistMin)) {
        //         ChanDistMin = dist;
        //         TickDistMin = p_hit->PeakTime();
        //     }

        //     MuonNHit++;
        //     MuonHitChannel.push_back(p_hit->Channel());
        //     MuonHitTick.push_back(p_hit->PeakTime());
        //     MuonHitADC.push_back(p_hit->Integral());

        //     i_hit++;
        // }

        // if (iLogLevel >= kInfos) {
        //     std::cout << "\t" << "track end channel: " << "\033[93m" << MuonEndChannel << "\033[0m" << std::endl;
        //     std::cout << "\t" << "------------------------------------------------" << std::endl;
        //     std::cout << "\t" << "chan dist min: " << "\033[93m" << MuonEndChannel + ChanDistMin << "\033[0m" << std::endl;
        //     std::cout << "\t" << "tick dist min: " << "\033[93m" << TickDistMin << "\033[0m" << std::endl;
        //     if (iUpMin * iUpMax < 0) std::cerr << "\033[91m" << "ERROR: iUpMin * iUpMax < 0" << "\033[0m" << std::endl;
        //     else if (iUpMin == -1) {
        //         std::cout << "\t" << "\033[91m" << "no hit in upper volume" << "\033[0m" << std::endl;
        //     } else {
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "up min channel: " << "\033[93m" << MuonHitChannel[iUpMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "up min tick: " << "\033[93m" << MuonHitTick[iUpMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "up max channel: " << "\033[93m" << MuonHitChannel[iUpMax] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "up max tick: " << "\033[93m" << MuonHitTick[iUpMax] << "\033[0m" << std::endl;
        //     }
        //     if (iLowMin * iLowMax < 0) std::cerr << "\033[91m" << "ERROR: iLowMin * iLowMax < 0" << "\033[0m" << std::endl;
        //     else if (iLowMin == -1) {
        //         std::cout << "\t" << "\033[91m" << "no hit in lower volume" << "\033[0m" << std::endl;
        //     } else {
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "low min channel: " << "\033[93m" << MuonHitChannel[iLowMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "low min tick: " << "\033[93m" << MuonHitTick[iLowMin] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "------------------------------------------------" << std::endl;
        //         std::cout << "\t" << "low max channel: " << "\033[93m" << MuonHitChannel[iLowMax] << "\033[0m" << std::endl;
        //         std::cout << "\t" << "low max tick: " << "\033[93m" << MuonHitTick[iLowMax] << "\033[0m" << std::endl;
        //     }
        // }
        // MuonTickLowMin = TickLowMin;
        // MuonTickLowMax = TickLowMax;
        // MuonTickUpMin = TickUpMin;
        // MuonTickUpMax = TickUpMax;


        if (MuonHasMichel != kNoMichel) {
            simb::MCParticle const * mcp_michel = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));
            if (Log(mcp_michel, kInfos, 1, "id to mcp...", "done", "failed")
            ) {
                std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());



                size_t n_hit_michel = 0;
                int bary_channel_michel = 0;
                float bary_tick_michel = 0;

                for (const recob::Hit* hit_michel : v_hit_michel) {
                    if (hit_michel->View() != geo::kW) continue;


                    std::cout << "bonjour" << std::endl;

                    n_hit_michel++;
                    bary_channel_michel += hit_michel->Channel();
                    bary_tick_michel += hit_michel->PeakTime();
                }

                // bary_channel_michel /= n_hit_michel;
                // bary_tick_michel /= n_hit_michel;

                if (iLogLevel >= kInfos) {
                    std::cout << "\t" << "------------------------------------------------" << std::endl;
                    if (n_hit_michel == 0) std::cout << "\t" << "\033[91m" << "no hit in michel" << "\033[0m" << std::endl;
                    else {
                        std::cout << "\t" << "michel hits: " << "\033[93m" << n_hit_michel << "\033[0m" << std::endl;
                        std::cout << "\t" << "michel barycenter channel: " << "\033[93m" << bary_channel_michel << "\033[0m" << std::endl;
                        std::cout << "\t" << "michel barycenter tick: " << "\033[93m" << bary_tick_michel << "\033[0m" << std::endl;
                    }
                }
            }
        }




        // if (iLogLevel >= kInfos) std::cout << "\t" << "muon nhit: " << "\033[93m" << tpMuonHitVars[kTick].size() << "\033[0m" << std::endl;
        // if (!Log(trk_end_chan != raw::InvalidChannelID, kInfos, 1, "valid muon end channel...", "yes", "no") ||
        //     !Log(tpTick.size() != 0, kInfos, 1, "muon hits in michel radius...", "yes", "no")) {
        //     AllMuonEndChan.push_back(-1);
        //     AllMuonEndTick.push_back(-1);
        //     continue;
        // }

        // float trk_end_tick = tpTick.at(std::distance(tpChanDist.begin(), std::min_element(tpChanDist.begin(), tpChanDist.end())));

        // AllMuonEndChan.push_back(trk_end_chan);
        // AllMuonEndTick.push_back(trk_end_tick);

        // if (iLogLevel >= kInfos) std::cout << "\t" << "muon end tick: " << "\033[93m" << trk_end_tick << "\033[0m" << std::endl; 
        // if (iLogLevel >= kInfos) std::cout << "\t" << "muon end channel: " << "\033[93m" << trk_end_chan << "\033[0m" << std::endl; 

        // float trk_min_tick=0, trk_max_tick=0;

        // if (iLogLevel >= kInfos) std::cout << "\t" << "looping over muon track's " << vp_hit.size() << " hits...";
        // for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        //     trk_min_tick = std::min(trk_min_tick, p_hit->PeakTime());
        //     trk_max_tick = std::max(trk_max_tick, p_hit->PeakTime());
        // }
        // if (iLogLevel >= kInfos) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;

        // if (iLogLevel >= kInfos) std::cout << "\t" << "muon ticks: " << trk_min_tick << " -> " << trk_max_tick;
        // float trk_end_tick;
        // if (abs(trk_end_pt.X()) < abs(trk_start_pt.X())) {
        //     trk_end_tick = trk_max_tick;

        //     if (iLogLevel >= kInfos) std::cout << " " << "\033[93m" << "max" << "\033[0m";
        // } else {
        //     trk_end_tick = trk_min_tick;

        //     if (iLogLevel >= kInfos) std::cout << " " << "\033[93m" << "min" << "\033[0m";
        // }
        // if (iLogLevel >= kInfos) std::cout << std::endl;
        // AllMuonEndTick.push_back(trk_end_tick);

        tMuon->Fill();
    } // end loop over tracks

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
    // if (iLogLevel >= kDetails) cout << "\033[93m" << "Muchecks::endJob: Plotting section =============================================" << "\033[0m" << endl;

    // auto const clockData = asDetClocks->DataForJob();
    // auto const detProp = asDetProp->DataForJob(clockData);

    // cout << "\t";
    // for (int k=0; k<kN; k++) cout << "| " << counterKeyNames[k] << "\t\t";
    // cout << endl;
    // cout << "--------";
    // for (int k=0; k<kN; k++) cout << "|-------" << "--------";
    // cout << endl;

    // for (int anti=0; anti<2; anti++) {
    //     cout << (anti ? "Anti" : "Muon") << "\t";
    //     for (int k=0; k<kN; k++) cout << "|" << viCounters[k][anti] << "\t" << fixed << setprecision(2) << 100.*viCounters[k][anti]/viCounters[kMuon][anti] << "%\t";
    //     cout << endl;
    // }

    // if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Muchecks::endJob ========================================================" << "\033[0m" << std::endl;
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
    MuonHitChannel.clear();
    MuonHitTick.clear();
    MuonHitADC.clear();

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


bool ana::Muchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}

DEFINE_ART_MODULE(ana::Muchecks)