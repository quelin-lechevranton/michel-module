////////////////////////////////////////////////////////////////////////
// Class:       Michecks
// Plugin Type: analyzer
// File:        Michecks_module.cc
//
// Generated at Wed Nov  6 10:31:28 2024 by Jeremy Quelin Lechevranton
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
using std::vector;
using std::string;
// using std::unordered_map;
// using std::pair;
// using std::min;
// using std::max;


namespace ana {
  class Michecks;
  struct Binning {
        int n;
        double min, max;
  };
}


class ana::Michecks : public art::EDAnalyzer {
public:
    explicit Michecks(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Michecks(Michecks const&) = delete;
    Michecks(Michecks&&) = delete;
    Michecks& operator=(Michecks const&) = delete;
    Michecks& operator=(Michecks&&) = delete;

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
    vector<vector<string>> vvsProducts;
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

    TTree* tEvent;

    size_t iEvent=0;

    // enum EnumTag { kMCP, kDaughters, kSphere, kTrue, kNTag };
    enum EnumTag { kMCP, kSphere, kTrue, kNTag };
    enum EnumHitVar { kTick, kChannel, kADC, kNVar };

    size_t EventNMichel;
    vector<size_t> EventiMichel; //[michel]
    vector<int> EventMichelPositron; //[michel]

    vector<int> EventMuonEndChan; //[michel]
    vector<float> EventMuonEndTick; //[michel]
    vector<float> EventMuonTrackLength; //[michel]

    vector<vector<float>> EventMichelEnergy; //[tag][michel]

    size_t EventNHit;
    vector<int> EventHitChannel; //[hit]
    vector<float> EventHitTick; //[hit]
    vector<float> EventHitADC; //[hit]

    TTree* tMichel;

    size_t iMichel=0;

    int MichelIsPositron;
    int MuonEndChannel;
    float MuonEndTick;
    float MuonTrackLength;

    vector<float> MichelEnergy; //[tag]

    vector<size_t> MichelNHit; //[tag]

    vector<vector<int>> MichelHitChannel; //[tag][hitFromMichel]
    vector<vector<float>> MichelHitTick; //[tag][hitFromMichel]
    vector<vector<float>> MichelHitADC; //[tag][hitFromMichel]

    size_t MuonNHit;
    vector<int> MuonHitChannel; //[hitFromMuon]
    vector<float> MuonHitTick; //[hitFromMuon]
    vector<float> MuonHitADC; //[hitFromMuon]


    // Diagnostic Variables
    size_t n_section = 4, n_plane = 3;
    size_t iNEventDiagnostic;
    ana::Binning binTick;
    vector<ana::Binning> binChan; //[section]
    vector<vector<TH2F*>> th2Hit; //[event][section]

    vector<vector<TGraph*>> tgMichelHit; //[event][tag]
    vector<TGraph*> tgMuonEnd; //[event]
    vector<TGraph*> tgMuonHit; //[event]

    vector<TGraph2D*> tg2Muon; //[event]
    vector<TGraph2D*> tg2MuonEnd; //[event]
    vector<TGraph2D*> tg2Michel; //[event]

    // Functions
    void resetEvent();
    void resetMichel();
    bool Log(bool cond, int flag, int tab, string msg, string succ, string fail);

    size_t GetSection(int ch);
    geo::View_t GetPlane(raw::ChannelID_t ch, int sec);
    geo::WireID GetWireID(geo::Point_t const& P, geo::View_t plane);
    raw::ChannelID_t GetChannel(geo::Point_t const& P, geo::View_t plane);
    // vector<art::Ptr<recob::Hit>> GetHits(detinfo::DetectorClocksData const& clockData, simb::MCParticle const& mpc, vector<art::Ptr<recob::Hit>> const& vp_hit);
    // simb::MCParticle const* GetMCParticle(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Track> const& p_trk, art::FindManyP<recob::Hit> const& fmp_trk2hit);

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


ana::Michecks::Michecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<vector<vector<string>>>("Products")),

    fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("TrackLengthCut")), // in cm

    iNEventDiagnostic(p.get<int>("NEventDiagnostic"))
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::Michecks: =============================================================" << "\033[0m" << std::endl;
    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    // Retrieving product tags
    for (vector<string> prod : vvsProducts) {

        const string    process     = prod[0],
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
    EventMichelEnergy = vector<vector<float>>(kNTag);

    tEvent->Branch("Event", &iEvent);

    tEvent->Branch("NMichel", &EventNMichel);
    tEvent->Branch("iMichel", &EventiMichel);

    tEvent->Branch("NHit", &EventNHit);
    tEvent->Branch("HitChannel", &EventHitChannel);
    tEvent->Branch("HitTick", &EventHitTick);
    tEvent->Branch("HitADC", &EventHitADC);


    tMichel = tfs->make<TTree>("Michel","");

    tMichel->Branch("Event", &iEvent);

    tMichel->Branch("IsPositron", &MichelIsPositron);

    tMichel->Branch("MuonEndChannel", &MuonEndChannel);
    tMichel->Branch("MuonEndTick", &MuonEndTick);
    tMichel->Branch("MuonTrackLength", &MuonTrackLength);

    MichelEnergy = vector<float>(kNTag);
    tMichel->Branch("TrueEnergy", (float*) &MichelEnergy[kTrue]);

    MichelNHit = vector<size_t>(kTrue);
    MichelHitChannel = vector<vector<int>>(kTrue);
    MichelHitTick = vector<vector<float>>(kTrue);
    MichelHitADC = vector<vector<float>>(kTrue);

    tMichel->Branch("MCPEnergy", (float*) &MichelEnergy[kMCP]);
    tMichel->Branch("MCPNHit", (size_t*) &MichelNHit[kMCP]);
    tMichel->Branch("MCPHitChannel", (vector<int>*) &MichelHitChannel[kMCP]);
    tMichel->Branch("MCPHitTick", (vector<float>*) &MichelHitTick[kMCP]);
    tMichel->Branch("MCPHitADC", (vector<float>*) &MichelHitADC[kMCP]);

    // tMichel->Branch("DaughtersEnergy", (float*) &MichelEnergy[kDaughters]);
    // tMichel->Branch("DaughtersNHit", (size_t*) &MichelNHit[kDaughters]);
    // tMichel->Branch("DaughtersHitChannel", (vector<int>*) &MichelHitChannel[kDaughters]);
    // tMichel->Branch("DaughtersHitTick", (vector<float>*) &MichelHitTick[kDaughters]);
    // tMichel->Branch("DaughtersHitADC", (vector<float>*) &MichelHitADC[kDaughters]);

    tMichel->Branch("SphereEnergy", (float*) &MichelEnergy[kSphere]);
    tMichel->Branch("SphereNHit", (size_t*) &MichelNHit[kSphere]);
    tMichel->Branch("SphereHitChannel", (vector<int>*) &MichelHitChannel[kSphere]);
    tMichel->Branch("SphereHitTick", (vector<float>*) &MichelHitTick[kSphere]);
    tMichel->Branch("SphereHitADC", (vector<float>*) &MichelHitADC[kSphere]);

    tMichel->Branch("MuonNHit", &MuonNHit);
    tMichel->Branch("MuonHitChannel", &MuonHitChannel);
    tMichel->Branch("MuonHitTick", &MuonHitTick);
    tMichel->Branch("MuonHitADC", &MuonHitADC);



    // Diagnostic Variables
    binTick.n = detProp.ReadOutWindowSize()/4;
    binTick.min = 0;
    binTick.max = detProp.ReadOutWindowSize();

    binChan.push_back({1175, 1900, 3075});
    binChan.push_back({1175, 4975, 6150});
    binChan.push_back({1175, 8045, 9220});
    binChan.push_back({1175, 11115, 12290});

    if (iLogLevel >= kInfos) {
        for (size_t s=0; s<n_section; s++) {
            std::cout << "section#" << s << " plane#" << geo::kW << std::endl;
            std::cout << "\t" << "Tick: " << binTick.min << " -> " << binTick.max << " (" << binTick.n << ")" << std::endl;
            std::cout << "\t" << "Chan: " << binChan[s].min << " -> " << binChan[s].max << " (" << binChan[s].n << ")" << std::endl;
        }
    }

    th2Hit = vector<vector<TH2F*>>(iNEventDiagnostic, vector<TH2F*>(n_section));
    tgMichelHit = vector<vector<TGraph*>>(iNEventDiagnostic, vector<TGraph*>(kTrue));
    tgMuonEnd = vector<TGraph*>(iNEventDiagnostic);
    tgMuonHit = vector<TGraph*>(iNEventDiagnostic);
    tg2Muon = vector<TGraph2D*>(iNEventDiagnostic);
    tg2MuonEnd = vector<TGraph2D*>(iNEventDiagnostic);
    tg2Michel = vector<TGraph2D*>(iNEventDiagnostic);
    for (size_t e=0; e<iNEventDiagnostic; e++) {
        for (size_t s=0; s<n_section; s++) {
            th2Hit[e][s] = new TH2F(
                Form("th2Hit_%ld_%ld",e,s), 
                "",
                binTick.n, binTick.min, binTick.max,
                binChan[s].n, binChan[s].min, binChan[s].max
            );
        }

        tgMuonEnd[e] = new TGraph();
        tgMuonHit[e] = new TGraph();
        tg2Muon[e] = new TGraph2D();
        tg2MuonEnd[e] = new TGraph2D();
        tg2Michel[e] = new TGraph2D();

        for (int t=0; t<kTrue; t++) {
            tgMichelHit[e][t] = new TGraph();
        }
    }


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Michecks::Michecks =======================================================" << "\033[0m" << std::endl;
}

void ana::Michecks::analyze(art::Event const& e)
{

    if (iLogLevel >= kBasics) {
        std::cout << "\033[93m" << "Michecks::analyze: Initialization evt#" << iEvent << "\r" << std::flush;
        std::cout << string(5,'\t') << " ======================================" << "\033[0m" << std::endl;
    }

    resetEvent();

    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();

    auto const vh_trk = e.getValidHandle<vector<recob::Track>>(tag_trk);
    vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);
        
    auto const & vh_hit = e.getValidHandle<vector<recob::Hit>>(tag_hit);
    vector<art::Ptr<recob::Hit>> vp_hit;
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
    

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::analyze: Muon Ends ===================================================" << "\033[0m" << std::endl;

    vector<simb::MCParticle const*> v_mcp_michel;
    vector<vector<vector<vector<float>>>> EventMichelHitVars(kTrue); //[tag][michel][hitVar][hit]
    vector<vector<vector<float>>> EventMuonHitVars; //[michel][hitVar][hit]
    vector<geo::Point_t> EventMuonEndPoint; //[michel]

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

        if (!Log(mcp->EndProcess() == "Decay", kDetails, 1, "is decaying...", "yes", "no")
        ) continue;

        geo::Point_t trk_end_pt, trk_start_pt;
        if (iLogLevel >= kInfos) std::cout << "\t" << p_trk->Start() << " -> " << p_trk->End();
        if (IsUpright(*p_trk)) {
            trk_end_pt = p_trk->End();
            trk_start_pt = p_trk->Start();
            if (iLogLevel >= kInfos) std::cout << " " << "\033[93m" << "upright" << "\033[0m";
        } else {
            trk_end_pt = p_trk->Start();
            trk_start_pt = p_trk->End();
            if (iLogLevel >= kInfos) std::cout << " " << "\033[93m" << "upside down" << "\033[0m";
        }
        if (iLogLevel >= kInfos) std::cout << std::endl;

        EventMuonEndPoint.push_back(trk_end_pt);

        if (!Log(IsInVolume(trk_start_pt, fMichelSpaceRadius), kInfos, 1, "starts in volume...", "yes", "no")
        ) continue;
        if (!Log(IsInVolume(trk_end_pt, fMichelSpaceRadius), kInfos, 1, "ends in volume...", "yes", "no")
        ) continue;

        raw::ChannelID_t trk_end_chan = GetChannel(trk_end_pt, geo::kW);

        if (!Log(trk_end_chan != raw::InvalidChannelID, kDetails, 1, "valid muon end channel...", "yes", "no")
        ) continue;



























        // if (!Log(IsInChannel(trk_end_chan, fMichelSpaceRadius), kDetails, 1, "ends in michel radius of channels boundaries...", "yes", "no")
        // ) continue;


        if (iLogLevel >= kInfos) std::cout << "\t" << "track length: " << "\033[93m" << p_trk->Length() << " cm" << "\033[0m" << std::endl;

        int i_dau = mcp->NumberDaughters() - 1;
        if (iLogLevel >= kInfos) std::cout << "\t" << "looping over " << "\033[93m" << "mu" << (mcp->PdgCode() < 0 ? "+" : "-") << "\033[0m" << "'s " << mcp->NumberDaughters() << " daughters..." << std::endl;
        for (; i_dau >= 0; i_dau--) {
            if (iLogLevel >= kDetails) std::cout << "\t" << "dau#" << i_dau+1 << "\r" << std::flush;
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            if (!Log(mcp_dau, kDetails, 2, "id to mcp...", "done", "failed")
            ) continue;

            if (!Log(abs(mcp_dau->PdgCode()) == 11 && mcp_dau->Process() == "Decay", kDetails, 2, "is michel...", "yes", "no")
            ) continue;

            if (!Log(IsInVolume(mcp_dau->Position(0), fMichelSpaceRadius), kDetails, 2, "is inside...", "yes", "no")
            ) i_dau = -1;
            break;
        } // end loop over muon daughters
        if (i_dau == -1) continue;
        if (iLogLevel >= kInfos) std::cout << "\t" << "michel found " << "\033[93m" << iMichel << (iMichel == 1 ? "st" : (iMichel == 2 ? "nd" : (iMichel == 3 ? "rd" : "th"))) << "\033[0m" << std::endl;

        EventNMichel++;
        EventiMichel.push_back(iMichel++);
        v_mcp_michel.push_back(pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau)));

        EventMuonTrackLength.push_back(p_trk->Length());

        !Log(IsInUpperVolume(trk_start_pt, 0.), kInfos, 1, "starts in upper volume...", "yes", "no");
        !Log(IsInLowerVolume(trk_start_pt, 0.), kInfos, 1, "starts in lower volume...", "yes", "no");
        !Log(IsInUpperVolume(trk_end_pt, 0.), kInfos, 1, "ends in upper volume...", "yes", "no");
        !Log(IsInLowerVolume(trk_end_pt, 0.), kInfos, 1, "ends in lower volume...", "yes", "no");


        if (iEvent < iNEventDiagnostic) {
            // for (size_t i = p_trk->FirstValidPoint(); i<p_trk->LastValidPoint(); i=p_trk->NextValidPoint(i)) {
            for (size_t i=0; i<p_trk->NumberTrajectoryPoints(); i++) {
                if (!p_trk->HasValidPoint(i)) continue;
                tg2Muon[iEvent]->AddPoint(p_trk->LocationAtPoint(i).Y(), p_trk->LocationAtPoint(i).Z(), p_trk->LocationAtPoint(i).X());
            }
            tg2MuonEnd[iEvent]->AddPoint(trk_end_pt.Y(), trk_end_pt.Z(), trk_end_pt.X());
            tg2Michel[iEvent]->AddPoint(v_mcp_michel.back()->Position(0).Y(), v_mcp_michel.back()->Position(0).Z(), v_mcp_michel.back()->Position(0).X());
            tg2Michel[iEvent]->AddPoint(v_mcp_michel.back()->Position(1).Y(), v_mcp_michel.back()->Position(1).Z(), v_mcp_michel.back()->Position(1).X());
        }
        

        vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk2hit.at(p_trk.key());

        vector<float> tpTick;
        vector<int> tpChanDist;
        vector<vector<float>> tpMuonHitVars(3); //[hitVar][hit]
        for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;

            tpMuonHitVars[kTick].push_back(p_hit->PeakTime());
            tpMuonHitVars[kChannel].push_back(p_hit->Channel());
            tpMuonHitVars[kADC].push_back(p_hit->Integral());

            if (trk_end_chan == raw::InvalidChannelID) continue;
            if (abs(int(p_hit->Channel() - trk_end_chan)) > fMichelSpaceRadius / fChannelPitch) continue;

            tpTick.push_back(p_hit->PeakTime());
            tpChanDist.push_back(abs(int(p_hit->Channel() - trk_end_chan)));
        }
        EventMuonHitVars.push_back(tpMuonHitVars);

        if (iLogLevel >= kInfos) std::cout << "\t" << "muon nhit: " << "\033[93m" << tpMuonHitVars[kTick].size() << "\033[0m" << std::endl;
        if (!Log(trk_end_chan != raw::InvalidChannelID, kInfos, 1, "valid muon end channel...", "yes", "no") ||
            !Log(tpTick.size() != 0, kInfos, 1, "muon hits in michel radius...", "yes", "no")) {
            EventMuonEndChan.push_back(-1);
            EventMuonEndTick.push_back(-1);
            continue;
        }

        float trk_end_tick = tpTick.at(std::distance(tpChanDist.begin(), std::min_element(tpChanDist.begin(), tpChanDist.end())));

        EventMuonEndChan.push_back(trk_end_chan);
        EventMuonEndTick.push_back(trk_end_tick);

        if (iLogLevel >= kInfos) std::cout << "\t" << "muon end tick: " << "\033[93m" << trk_end_tick << "\033[0m" << std::endl; 
        if (iLogLevel >= kInfos) std::cout << "\t" << "muon end channel: " << "\033[93m" << trk_end_chan << "\033[0m" << std::endl; 

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
        // EventMuonEndTick.push_back(trk_end_tick);
    } // end loop over tracks

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::analyze: Michel Energy ===============================================" << "\033[0m" << std::endl;

    if (iLogLevel >= kInfos) std::cout << "looping over " << v_mcp_michel.size() << " michels' mcp..." << std::endl;
    int i_mich = 0;
    for (simb::MCParticle const* mcp_mich : v_mcp_michel) {
        if (iLogLevel >= kInfos) std::cout << "mich#" << i_mich << "\r" << std::flush;

        EventMichelPositron.push_back(mcp_mich->PdgCode() < 0);
        EventMichelEnergy[kTrue].push_back((mcp_mich->E() - mcp_mich->Mass())*1e3);

        if (EventMichelEnergy[kTrue].back() > 60) {
            if (iLogLevel >= kImportant) std::cout << "\033[91m" << "too many energy" << "\033[0m" << " @ e" << iEvent << " m" << i_mich << std::endl;
        }

        if (iLogLevel >= kInfos) std::cout << "\t" << "\033[93m" << (mcp_mich->PdgCode() < 0 ? "positron" : "electron") << "\033[0m" << std::endl;
        if (iLogLevel >= kInfos) std::cout << "\t" << "true energy: " << "\033[93m" << EventMichelEnergy[kTrue].back() << " MeV" << "\033[0m" << std::endl;

        // vector<art::Ptr<recob::Hit>> vp_hit_michel = GetHits(clockData, *mcp_mich, vp_hit);
        vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_mich, e, tag_hit.label());

        if (v_hit_michel.size() == 0) {
            if (iLogLevel >= kImportant) std::cout << "\033[91m" << "0 hits" << "\033[0m" << " @ e" << iEvent << " m" << i_mich << std::endl;
        }

        float McpEnergy = 0;
        float previousADC = -1;
        vector<vector<float>> tpMichelHitVars(3); //[hitVar][hit]
        if (iLogLevel >= kInfos) std::cout << "\t" << "looping over " << v_hit_michel.size() << " hits...";
        for (recob::Hit const* hit_michel : v_hit_michel) {

            if (hit_michel->View() != geo::kW) continue;

            tpMichelHitVars[kTick].push_back(hit_michel->PeakTime());
            tpMichelHitVars[kChannel].push_back(hit_michel->Channel());
            tpMichelHitVars[kADC].push_back(hit_michel->Integral());

            // McpEnergy += hit_michel->HitSummedADC();

            if (hit_michel->ROISummedADC() == previousADC) continue;
            McpEnergy += hit_michel->ROISummedADC();
            previousADC = hit_michel->ROISummedADC();
        } // end loop over michel hits
        if (iLogLevel >= kInfos) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;

        EventMichelHitVars[kMCP].push_back(tpMichelHitVars);
        EventMichelEnergy[kMCP].push_back(McpEnergy*fADCtoE);

        if (iLogLevel >= kInfos) std::cout << "\t" << "MCP energy: " << "\033[93m" << McpEnergy << " ADC.Tick ~ " << EventMichelEnergy[kMCP].back() << " MeV" << "\033[0m" << std::endl;

        // float DaughtersEnergy = McpEnergy;
        // previousADC = -1;
        // tpMichelHitVars = vector<vector<float>>(3);
        // if (iLogLevel >= kInfos) std::cout << "\t" << "looping over " << mcp_mich->NumberDaughters() << " daughters...";
        // if (iLogLevel >= kDetails) std::cout << std::endl;
        // for (int i_dau=0; i_dau< mcp_mich->NumberDaughters(); i_dau++) {
        //     if (iLogLevel >= kDetails) std::cout << "\t" << "dau#" << i_dau+1 << "\r" << std::flush;

        //     simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp_mich->Daughter(i_dau));

        //     if (!Log(mcp_dau, kDetails, 2, "id to mcp...", "done", "failed")
        //     ) continue;

        //     const int pdg = mcp_dau->PdgCode();
        //     if (!Log(abs(pdg) == 11 || abs(pdg) == 22, kDetails, 2, "is electron or photon...", Form("%s", abs(pdg)==22 ? "photon" : (pdg>0 ? "elec" : "posi")), "no")
        //     ) continue;
            
        //     // vector<art::Ptr<recob::Hit>> vp_hit_dau = GetHits(clockData, *mcp_dau, vp_hit);
        //     vector<const recob::Hit*> v_hit_dau = truthUtil.GetMCParticleHits(clockData, *mcp_dau, e, tag_hit.label());

        //     if (iLogLevel >= kImportant && abs(pdg)==22 && v_hit_dau.size() > 0) std::cout << "\033[91m" << "photon with " << v_hit_dau.size() << " hits" << "\033[0m" << std::endl;

        //     if (iLogLevel >= kDetails) std::cout << "\t\t" << "looping over " << v_hit_dau.size() << " hits...";
        //     for (recob::Hit const* hit_dau : v_hit_dau) {

        //         std::cout << hit_dau->View() << std::endl;

        //         if (hit_dau->View() != geo::kW) continue;
                
        //         tpMichelHitVars[kTick].push_back(hit_dau->PeakTime());
        //         tpMichelHitVars[kChannel].push_back(hit_dau->Channel());
        //         tpMichelHitVars[kADC].push_back(hit_dau->Integral());

        //         std::cout << tpMichelHitVars[kTick].size() << std::endl;

        //         // DaughtersEnergy += hit_dau->HitSummedADC();

        //         if (hit_dau->ROISummedADC() == previousADC) continue;
        //         DaughtersEnergy += hit_dau->ROISummedADC();
        //         previousADC = hit_dau->ROISummedADC();
        //     } // end loop over electron hits
        //     if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
        // } // end loop over michel daughters
        // if (iLogLevel >= kInfos && iLogLevel < kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
        // EventMichelHitVars[kDaughters].push_back(tpMichelHitVars);
        // EventMichelEnergy[kDaughters].push_back(DaughtersEnergy*fADCtoE);

        // if (iLogLevel >= kInfos) std::cout << "\t" << "daughters energy: " << "\033[93m" << DaughtersEnergy << " ADC.Tick" << "\033[0m"
        //                                 << " ~ " << "\033[93m" << EventMichelEnergy[kDaughters].back() << " MeV" << "\033[0m"
        //                                 << " over " << "\033[93m" << tpMichelHitVars[kTick].size() << "\033[0m" << " hits" << std::endl;

        i_mich++;
    } // end loop over michel


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::analyze: Michel Sphere ===============================================" << "\033[0m" << std::endl;

    EventMichelEnergy[kSphere] = vector<float>(EventNMichel, 0);
    vector<vector<vector<float>>> tpSphereHitVars(EventNMichel, vector<vector<float>>(3));

    vector<vector<float>> tpADC(EventNMichel); //[michel][hit]
    if (iLogLevel >= kInfos) std::cout << "looping over all hits...";
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {

        if (p_hit->View() != geo::kW) continue;

        vector<art::Ptr<recob::Track>> vp_trk = fmp_hit2trk.at(p_hit.key());

        bool BelongsToLongTrack = false;
        for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
            BelongsToLongTrack = p_trk->Length() > fTrackLengthCut;
            if (BelongsToLongTrack) break;
        }
        if (BelongsToLongTrack) continue;

        for (size_t m=0; m<EventNMichel; m++) {

            if (EventMuonEndChan[m] == -1) continue;

            float X = (int(p_hit->Channel()) - EventMuonEndChan[m]) * fChannelPitch / fMichelSpaceRadius;
            float Y = (p_hit->PeakTime() - EventMuonEndTick[m]) * fSamplingRate / fMichelTimeRadius;

            if (X*X + Y*Y > 1) continue;

            tpSphereHitVars[m][kTick].push_back(p_hit->PeakTime());
            tpSphereHitVars[m][kChannel].push_back(p_hit->Channel());
            tpSphereHitVars[m][kADC].push_back(p_hit->Integral());

            tpADC[m].push_back(p_hit->ROISummedADC());
        } // end loop over michels
    } // end loop over hits
    if (iLogLevel >= kInfos) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
    EventMichelHitVars[kSphere] = tpSphereHitVars;

    if (iLogLevel >= kInfos) std::cout << "looping over " << EventNMichel << " michels..." << std::endl;
    for (size_t m=0; m<EventNMichel; m++) {
        if (iLogLevel >= kInfos) std::cout << "mich#" << m << "\r" << std::flush;

        // delete duplicates
        std::sort(tpADC[m].begin(), tpADC[m].end());
        tpADC[m].erase(std::unique(tpADC[m].begin(), tpADC[m].end()), tpADC[m].end());
        
        float SphereEnergy = std::accumulate(tpADC[m].begin(), tpADC[m].end(), 0.0);
        EventMichelEnergy[kSphere][m] = SphereEnergy * fADCtoE;

        if (iLogLevel >= kInfos) std::cout << "\t" << "sphere energy: " << "\033[93m" << SphereEnergy << " ADC.Tick ~ " << EventMichelEnergy[kSphere][m] << " MeV" << "\033[0m" << std::endl;
    } // end loop over michels


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::analyze: Filling trees ===============================================" << "\033[0m" << std::endl;

    tEvent->Fill();

    if (iLogLevel >= kDetails) std::cout << "Event tree filled" << std::endl;

    for (size_t m=0; m<EventNMichel; m++) {

        if (iLogLevel >= kDetails) std::cout << "mich#" << m << "\r" << std::flush;

        resetMichel();

        if (iLogLevel >= kDetails) std::cout << "\t" << "step 1" << std::endl;

        MichelIsPositron = EventMichelPositron[m];
        MuonEndChannel = EventMuonEndChan[m];
        MuonEndTick = EventMuonEndTick[m];
        MuonTrackLength = EventMuonTrackLength[m];

        if (iLogLevel >= kDetails) std::cout << "\t" << "step 2" << std::endl;

        MichelEnergy[kTrue] = EventMichelEnergy[kTrue][m];
        MichelEnergy[kMCP] = EventMichelEnergy[kMCP][m];
        // MichelEnergy[kDaughters] = EventMichelEnergy[kDaughters][m];
        MichelEnergy[kSphere] = EventMichelEnergy[kSphere][m];

        if (iLogLevel >= kDetails) std::cout << "\t" << "step 3" << std::endl;

        MichelNHit[kMCP] = EventMichelHitVars[kMCP][m][kTick].size();
        MichelHitChannel[kMCP] = vector<int>(EventMichelHitVars[kMCP][m][kChannel].begin(), EventMichelHitVars[kMCP][m][kChannel].end());
        MichelHitTick[kMCP] = EventMichelHitVars[kMCP][m][kTick];
        MichelHitADC[kMCP] = EventMichelHitVars[kMCP][m][kADC];

        // MichelNHit[kDaughters] = EventMichelHitVars[kDaughters][m][kTick].size();
        // MichelHitChannel[kDaughters] = vector<int>(EventMichelHitVars[kDaughters][m][kChannel].begin(), EventMichelHitVars[kDaughters][m][kChannel].end());
        // MichelHitTick[kDaughters] = EventMichelHitVars[kDaughters][m][kTick];
        // MichelHitADC[kDaughters] = EventMichelHitVars[kDaughters][m][kADC];

        if (iLogLevel >= kDetails) std::cout << "\t" << "step 4" << std::endl;

        MichelNHit[kSphere] = EventMichelHitVars[kSphere][m][kTick].size();
        MichelHitChannel[kSphere] = vector<int>(EventMichelHitVars[kSphere][m][kChannel].begin(), EventMichelHitVars[kSphere][m][kChannel].end());
        MichelHitTick[kSphere] = EventMichelHitVars[kSphere][m][kTick];
        MichelHitADC[kSphere] = EventMichelHitVars[kSphere][m][kADC];

        if (iLogLevel >= kDetails) std::cout << "\t" << "step 5" << std::endl;

        MuonNHit = EventMuonHitVars[m][kTick].size();
        MuonHitChannel = vector<int>(EventMuonHitVars[m][kChannel].begin(), EventMuonHitVars[m][kChannel].end());
        MuonHitTick = EventMuonHitVars[m][kTick];
        MuonHitADC = EventMuonHitVars[m][kADC];

        if (iLogLevel >= kDetails) std::cout << "\t" << "step 6" << std::endl;

        tMichel->Fill();
    }

    // Diagnostic
    if (iEvent >= iNEventDiagnostic) { iEvent++; return; }
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::analyze: Diagnostic ==================================================" << "\033[0m" << std::endl;

    if (iLogLevel >= kInfos) std::cout << "evt#" << iEvent << "\r" << std::flush;

    for (size_t i=0; i<EventNHit; i++) {
        for (size_t s=0; s<n_section; s++) {
            if (GetSection(EventHitChannel[i]) != s) continue;
            th2Hit[iEvent][s]->Fill(EventHitTick[i], EventHitChannel[i], EventHitADC[i]);
        }
    }

    for (size_t m=0; m<EventNMichel; m++) {
        if (iLogLevel >= kInfos) std::cout << "\t" << "mich#" << m << "\r" << std::flush;

        tgMuonEnd[iEvent]->AddPoint(EventMuonEndTick[m], EventMuonEndChan[m]);

        if (iLogLevel >= kInfos) {
            std::cout << "\t\t" << (EventMichelPositron[m] ? "positron" : "electron") << "\r" << std::flush;
            std::cout << string(5,'\t') << "@ " << EventMuonEndPoint[m] << std::endl;
            std::cout << string(5,'\t') << "@ tick:" << EventMuonEndTick[m] << " chan:" << EventMuonEndChan[m] << std::endl;
        }

        for (int t=0; t<kTrue; t++) {
            int n = EventMichelHitVars[t][m][kTick].size();

            for (int i=0; i<n; i++) tgMichelHit[iEvent][t]->AddPoint(EventMichelHitVars[t][m][kTick][i], EventMichelHitVars[t][m][kChannel][i]);

            if (iLogLevel >= kInfos) {
                // string tag = (t==kMCP ? "MCP" : (t==kDaughters ? "Daughters" : "Sphere"));
                string tag = (t==kMCP ? "MCP" : (t==kSphere ? "Sphere" : "???"));
                std::cout << "\t\t" << tag << ": " << n << " hits" << "\r" << std::flush;
                if (n > 0) std::cout << string(5,'\t') << "@ tick:" << EventMichelHitVars[t][m][kTick][0] << " chan:" << EventMichelHitVars[t][m][kChannel][0];
                std::cout << std::endl;
            }
        }

        for (unsigned int i=0; i<EventMuonHitVars[m][kTick].size(); i++) tgMuonHit[iEvent]->AddPoint(EventMuonHitVars[m][kTick][i], EventMuonHitVars[m][kChannel][i]);

    }
    iEvent++;

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Michecks::analyze =======================================================" << "\033[0m" << std::endl;
} // end analyze

void ana::Michecks::beginJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::beginJob: ============================================================" << "\033[0m" << std::endl;
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Michecks::beginJob ======================================================" << "\033[0m" << std::endl;
} // end beginJob


void ana::Michecks::endJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::endJob: Plotting =====================================================" << "\033[0m" << std::endl;

    auto const clockData = asDetClocks->DataForJob();
    // auto const detProp = asDetProp->DataForJob(clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Michecks::endJob: Diagnostic ===================================================" << "\033[0m" << std::endl;

    //Tag { kMCP, kDaughters, kSphere, kTrue, kNTag };
    vector<short> markerStyle = {kFullDoubleDiamond, kFullFourTrianglesX, kFullCircle};
    vector<float> markerSize = {2.5, 2.5, .8};
    // vector<short> markerColor = {kAzure-5, kGray+1, kMagenta-3};
    vector<short> markerColor = {kAzure-5, kMagenta-3};

    for (size_t e=0; e<iNEventDiagnostic; e++) {
        TCanvas* c2d = tfs->make<TCanvas>(Form("c2d_e%ld",e),Form("2D Event %ld",e));
        c2d->Divide(2,2);

        for (size_t s=0; s<n_section; s++) {
            c2d->cd(s+1);
            gPad->DrawFrame(
                binTick.min, binChan[s].min,
                binTick.max, binChan[s].max,
                Form("Section %ld, W plane;Tick;Channel",s)
            );
            gPad->SetLogz();

            th2Hit[e][s]->SetMinimum(0.1);
            th2Hit[e][s]->Draw("colz same");
            

            for (int m=0; m<tgMuonEnd[e]->GetN(); m++) {
                auto te = new TEllipse(
                    tgMuonEnd[e]->GetPointX(m),
                    tgMuonEnd[e]->GetPointY(m),
                    fMichelTimeRadius / fSamplingRate,
                    fMichelSpaceRadius / fChannelPitch 
                );
                te->SetFillStyle(0);
                te->SetLineWidth(2);
                te->SetLineColor(markerColor[kSphere]);
                te->Draw();
            }

            for (int t=0; t<kTrue; t++) {
                if (tgMichelHit[e][t]->GetN() == 0) continue;
                tgMichelHit[e][t]->SetEditable(kFALSE);
                tgMichelHit[e][t]->SetMarkerStyle(markerStyle[t]);
                tgMichelHit[e][t]->SetMarkerSize(markerSize[t]);
                tgMichelHit[e][t]->SetMarkerColor(markerColor[t]);
                // tgMichelHit[e][t]->SetMarkerColorAlpha(markerColor[t], 0.5);
                tgMichelHit[e][t]->Draw("P same");
            }
            if (tgMuonHit[e]->GetN() == 0) continue;
            tgMuonHit[e]->SetEditable(kFALSE);
            tgMuonHit[e]->SetMarkerStyle(kOpenSquare);
            tgMuonHit[e]->SetMarkerSize(1);
            tgMuonHit[e]->SetMarkerColor(kPink-5);
            tgMuonHit[e]->Draw("P same");
        }
        c2d->Write();

        TCanvas* c3d = tfs->make<TCanvas>(Form("c3d_e%ld",e),Form("3D Event %ld",e));

        TH3F *axes = new TH3F(Form("axes_e%ld",e), ";Y (cm);Z (cm);X (cm)", 1, -400, 400, 1, -10, 400, 1, -400, 400);
        axes->SetDirectory(0); // prevent it from being written to the file
        axes->SetStats(0); // prevent the stats box from being drawn



        c3d->cd();
        axes->Draw();

        tg2Muon[e]->SetMarkerStyle(kFullSquare);
        tg2Muon[e]->SetMarkerSize(.5);
        tg2Muon[e]->SetMarkerColor(kPink-5);
        tg2Muon[e]->Draw("P same");

        tg2Michel[e]->SetMarkerStyle(kFullDoubleDiamond);
        tg2Michel[e]->SetMarkerSize(2.5);
        tg2Michel[e]->SetMarkerColor(kAzure-5);
        tg2Michel[e]->Draw("P same");

        tg2MuonEnd[e]->SetMarkerStyle(kFullCircle);
        tg2MuonEnd[e]->SetMarkerSize(1);
        tg2MuonEnd[e]->SetMarkerColor(kOrange+5);
        tg2MuonEnd[e]->Draw("P same");

        c3d->Write();
    }

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Michecks::endJob ========================================================" << "\033[0m" << std::endl;
} // end endJob


void ana::Michecks::resetEvent() {
    if (iLogLevel >= kDetails) std::cout << "resetting event branches...";

    EventNMichel = 0;
    EventiMichel.clear();

    EventMichelPositron.clear();

    EventMuonEndChan.clear();
    EventMuonEndTick.clear();
    EventMuonTrackLength.clear();

    for (int t=0; t<kNTag; t++) EventMichelEnergy[t].clear();

    EventNHit = 0;
    EventHitChannel.clear();
    EventHitTick.clear();
    EventHitADC.clear();

    if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
}
void ana::Michecks::resetMichel() {
    if (iLogLevel >= kDetails) std::cout << "resetting michel branches...";

    // for (int t=0; t<kTrue; t++) {
    //     MichelHitChannel[t].clear();
    //     MichelHitTick[t].clear();
    //     MichelHitADC[t].clear();
    // }

    if (iLogLevel >= kDetails) std::cout << "\033[92m" << " done" << "\033[0m" << std::endl;
}

bool ana::Michecks::Log(bool cond, int flag, int tab, string msg, string succ, string fail) {
    if (iLogLevel >= flag) {
        std::cout << string(tab,'\t') << msg << " ";
        if (cond) std::cout << "\033[92m" << succ << "\033[0m" << std::endl;
        else std::cout << "\033[91m" << fail << "\033[0m" << std::endl;
    }
    return cond;
}




size_t ana::Michecks::GetSection(int ch) {
    return int(4.*ch / asWire->Nchannels());
}
geo::View_t ana::Michecks::GetPlane(raw::ChannelID_t ch, int sec) {
    return static_cast<geo::View_t>(12.*ch / asWire->Nchannels() - 3*sec);
}
geo::WireID ana::Michecks::GetWireID(geo::Point_t const& P, geo::View_t plane) {
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
raw::ChannelID_t ana::Michecks::GetChannel(geo::Point_t const& P, geo::View_t plane) {
    geo::WireID wireid = GetWireID(P, plane);
    if (!wireid.isValid) return raw::InvalidChannelID;
    return asWire->PlaneWireToChannel(wireid);
}

// https://github.com/DUNE/protoduneana/blob/develop/protoduneana/Utilities/ProtoDUNETruthUtils.cxx#L105
// vector<art::Ptr<recob::Hit>> ana::Michecks::GetHits(detinfo::DetectorClocksData const& clockData, simb::MCParticle const& mpc, vector<art::Ptr<recob::Hit>> const& vp_hit) {
//     vector<art::Ptr<recob::Hit>> vp_mcp_hit;
//     for (art::Ptr<recob::Hit> const& p_hit : vp_hit ) {
//         for (sim::TrackIDE const& ide : bt_serv->HitToTrackIDEs(clockData, p_hit)) {
//             if (ide.trackID == mpc.TrackId()) {
//             // if (pi_serv->TrackIdToParticle_P(ide.trackID) == &mpc) {
//             // if (pi_serv->TrackIdToParticle_P(ide.trackID) == pi_serv->TrackIdToParticle_P(mpc.TrackId())) {
//                 vp_mcp_hit.push_back(p_hit);
//             }
//         }
//     }
//     return vp_mcp_hit;
// }


// https://github.com/DUNE/protoduneana/blob/develop/protoduneana/Utilities/ProtoDUNETruthUtils.cxx#L610
// https://github.com/DUNE/protoduneana/blob/develop/protoduneana/Utilities/ProtoDUNETruthUtils.cxx#L407
// https://github.com/DUNE/protoduneana/blob/develop/protoduneana/Utilities/ProtoDUNETrackUtils.cxx#L95
// https://github.com/DUNE/protoduneana/blob/develop/protoduneana/Utilities/ProtoDUNETruthUtils.cxx#L238
// simb::MCParticle const* ana::Michecks::GetMCParticle(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Track> const& p_trk, art::FindManyP<recob::Hit> const& fmp_trk2hit) {
//     vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk2hit.at(p_trk.key());
//     // int i_max_ide_energy = 0;
//     // float max_ide_energy = 0;
//     std::unordered_map<simb::MCParticle const*, float> um_mcp_energy;
//     for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
//         for (sim::TrackIDE const& ide : bt_serv->HitToTrackIDEs(clockData, p_hit)) {
//             um_mcp_energy[pi_serv->TrackIdToParticle_P(ide.trackID)] += ide.energy;
            
//             // if (ide.energy > max_ide_energy) {
//             //     max_ide_energy = ide.energy;
//             //     i_max_ide_energy = ide.trackID;
//             // }
//         }
//     }
//     float max_energy = 0;
//     simb::MCParticle const* max_mcp = nullptr;
//     // vector<std::pair<simb::MCParticle const*, float>> v_mcp_energy;
//     for (std::pair<simb::MCParticle const*, float> p : um_mcp_energy) {
//         if (p.second > max_energy) {
//             max_energy = p.second;
//             max_mcp = p.first;
//         }
//         // v_mcp_energy.push_back(p);
//     }
//     // std::sort(v_mcp_energy.begin(), v_mcp_energy.end(), [](std::pair<simb::MCParticle const*, float> a, std::pair<simb::MCParticle const*, float> b) {return a.second > b.second;});
//     return max_mcp;
//     // return v_mcp_energy[0].first;
//     // return pi_serv->TrackIdToParticle_P(i_max_ide_energy);
// }



bool ana::Michecks::IsInLowerVolume(double x, double y, double z, float eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpcid0{cryoid, 0}, tpcidN{cryoid, 7};
    return (x-eps > asGeo->TPC(tpcid0).MinX() and x+eps < asGeo->TPC(tpcidN).MaxX() and
            y-eps > asGeo->TPC(tpcid0).MinY() and y+eps < asGeo->TPC(tpcidN).MaxY() and
            z-eps > asGeo->TPC(tpcid0).MinZ() and z+eps < asGeo->TPC(tpcidN).MaxZ());   
}
bool ana::Michecks::IsInUpperVolume(double x, double y, double z, float eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpcid0{cryoid, 8}, tpcidN{cryoid, asGeo->NTPC(cryoid)-1};
    return (x-eps > asGeo->TPC(tpcid0).MinX() and x+eps < asGeo->TPC(tpcidN).MaxX() and
            y-eps > asGeo->TPC(tpcid0).MinY() and y+eps < asGeo->TPC(tpcidN).MaxY() and
            z-eps > asGeo->TPC(tpcid0).MinZ() and z+eps < asGeo->TPC(tpcidN).MaxZ());   
}
bool ana::Michecks::IsInVolume(double x, double y, double z, float eps) {
    return IsInLowerVolume(x, y, z, eps) or IsInUpperVolume(x, y, z, eps);
}
bool ana::Michecks::IsInLowerVolume(float x, float y, float z, float eps) {
    return IsInLowerVolume(double(x), double(y), double(z), eps);
}
bool ana::Michecks::IsInUpperVolume(float x, float y, float z, float eps) {
    return IsInUpperVolume(double(x), double(y), double(z), eps);
}
bool ana::Michecks::IsInVolume(float x, float y, float z, float eps) {
    return IsInVolume(double(x), double(y), double(z), double(eps));
}
template<class Vec>
bool ana::Michecks::IsInUpperVolume(Vec const& V, float eps) {
    return IsInUpperVolume(V.X(), V.Y(), V.Z(), eps);
}
template<class Vec>
bool ana::Michecks::IsInLowerVolume(Vec const& V, float eps) {
    return IsInLowerVolume(V.X(), V.Y(), V.Z(), eps);
}
template<class Vec>
bool ana::Michecks::IsInVolume(Vec const& V, float eps) {
    return IsInVolume(V.X(), V.Y(), V.Z(), eps);
}


bool ana::Michecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}

DEFINE_ART_MODULE(ana::Michecks)
