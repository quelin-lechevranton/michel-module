////////////////////////////////////////////////////////////////////////
// Class:       Muchecks
// Plugin Type: analyzer
// File:        Muchecks_module.cc
//
// Generated at Thu Nov 28 13:09:28 2024 by Jeremy Quelin Lechevranton
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
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEllipse.h"

// std includes
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>

using std::cerr;
using std::cout;
using std::fixed;
using std::setprecision;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::min;
using std::max;

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

    const geo::Geometry* asGeo;
    const detinfo::DetectorPropertiesService* asDetProp;
    const detinfo::DetectorClocksService* asDetClocks;

    protoana::ProtoDUNETruthUtils truthUtil;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Conversion factors
    float fADCtoE = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    float fChannelPitch = 0.5; // cm/channel
    
    // Verbosity
    int iLogLevel;
    int iFlagSpecial = 1;
    int iFlagDetails = 2;

    // Products
    vector<vector<string>> vvsProducts;
    art::InputTag   tag_mcp,
                    tag_sed,
                    tag_wir,
                    tag_hit,
                    tag_clu,
                    tag_trk,
                    tag_spt;

    // Variables
    float fTrackLengthCut; // in cm
    float fMichelSpaceRadius; // in cm


    vector<vector<int>> viCounters;
    enum counterKeys { kMuon, kDecay, kInside, kDecayIn, kDIDauless, kDIMuDau, kN };
    vector<string> counterKeyNames = {"Muon", "Decay", "Inside", "DecIn", "w/0Dau", "w/µDau"};


    // Functions
    geo::WireID GetWireID(geo::Point_t const& P, int plane);
    raw::ChannelID_t GetChannel(geo::Point_t const& P, int plane);

    bool IsInVolume(double x, double y, double z, double eps);
    bool IsInVolume(float x, float y, float z, float eps);
    bool IsInVolume(TLorentzVector const& V, float eps);
    bool IsInVolume(geo::Point_t const& P);

    bool IsUpright(recob::Track const& T);

    int GetSection(unsigned int ch);
    int GetPlane(unsigned int ch, int sec);

    bool Logging(bool cond, int tab, string msg, string succ, string fail);
    bool LogFail(bool cond, int tab, string msg, string fail);
    bool LogSucc(bool cond, int tab, string msg, string succ);

    int iCorr(int n, int i, int j);
};


ana::Muchecks::Muchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel",1)),
    vvsProducts(p.get<vector<vector<string>>>("Products")),

    fTrackLengthCut(p.get<float>("LengthCut")), // in cm
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")) // in cm
{

    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
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

        const auto  tag = art::InputTag(label,instance);

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
    }


    // Initialize counters
    viCounters = vector<vector<int>>(kN, vector<int>(2, 0));
}

void ana::Muchecks::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    // auto const detProp = asDetProp->DataFor(e,clockData);

    if (iLogLevel >= iFlagDetails) cout << "evt#" << e.id().event() << "\r" << flush;

    auto const vh_trk = e.getValidHandle<vector<recob::Track>>(tag_trk);
    vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    if (iLogLevel >= iFlagDetails) cout << "\tlooping over " << vp_trk.size() << " tracks..." << endl;
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        if (iLogLevel >= iFlagDetails) cout << "\ttrk#" << p_trk->ID()+1 << "\r" << flush;

        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());

        if (LogFail(
            !mcp, 
            2, "trk to mcp...", "failed")
        ) continue;

        if (Logging(
            abs(mcp->PdgCode()) != 13,
            2, "is muon...", "yes", "no")
        ) continue;

        int anti = mcp->PdgCode() > 0 ? 0 : 1;
        viCounters[kMuon][anti]++;

        if (Logging(
            !IsInVolume(mcp->EndPosition(), fMichelSpaceRadius),
            2, "ends inside...", "yes", "no")
        ) {
            if (Logging(
                mcp->EndProcess() != "Decay",
                2, "is decaying...", "yes", "no")
            ) {} else viCounters[kDecay][anti]++;
            continue;
        }
        
        viCounters[kInside][anti]++;

        if (Logging(
            mcp->EndProcess() != "Decay",
            2, "is decaying...", "yes", "no")
        ) continue;
        
        viCounters[kDecay][anti]++;
        viCounters[kDecayIn][anti]++;

        if (Logging(
            mcp->NumberDaughters() == 0,
            2, Form("looping over muon's %d daughters...",mcp->NumberDaughters()), "", "none")
        ) {viCounters[kDIDauless][anti]++; continue;}

        bool has_muon = false;
        for (int i_dau=0; i_dau < mcp->NumberDaughters(); i_dau++) {
            if (iLogLevel >= iFlagDetails) cout << "\t\tdau#" << i_dau+1 << "\r" << flush;
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            if (LogFail(
                !mcp_dau,
                3, "id to mcp...", "failed")
            ) continue;

            if (LogSucc(
                abs(mcp_dau->PdgCode()) != 13,
                3, "is muon...", "yes")
            ) continue;

            if (iLogLevel >= iFlagDetails) cout << "\t\t\tMuon daughter energy: " << mcp_dau->E() << endl;

            has_muon = true;
            break;
        } // end loop over muon daughters

        if (!has_muon) continue;

        viCounters[kDIMuDau][anti]++;


        //     if (Logging(
        //         !IsInVolume(mcp_dau->Position(0), fMichelSpaceRadius),
        //         3, "is inside...", "yes", "no")
        //     ) i_dau = -1;
        //     break;
        // } // end loop over muon daughters
        // if (i_dau == -1) continue;

        // simb::MCParticle const * mcp_mich = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));

        // if (IsUpright(*p_trk)) {
        //     trk_end = p_trk->End();
        // } else {
        //     trk_end = p_trk->Start();
        // }

        // raw::ChannelID_t ch_mu_mcp_end = GetChannel(geo::Point_t(mcp_mich->EndPosition().Vect()), 2);
        // raw::ChannelID_t ch_mu_trk_end = GetChannel(trk_end, 2);
    } // end loop over tracks
} // end analyze

void ana::Muchecks::beginJob() {} // end beginJob


void ana::Muchecks::endJob()
{

    if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "Muchecks::endJob: Plotting section =============================================" << "\033[0m" << endl;

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    cout << "\t";
    for (int k=0; k<kN; k++) cout << "| " << counterKeyNames[k] << "\t\t";
    cout << endl;
    cout << "--------";
    for (int k=0; k<kN; k++) cout << "|-------" << "--------";
    cout << endl;

    for (int anti=0; anti<2; anti++) {
        cout << (anti ? "Anti" : "Muon") << "\t";
        for (int k=0; k<kN; k++) cout << "|" << viCounters[k][anti] << "\t" << fixed << setprecision(2) << 100.*viCounters[k][anti]/viCounters[kMuon][anti] << "%\t";
        cout << endl;
    }

    if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "End of Muchecks::endJob ========================================================" << "\033[0m" << endl;
} // end endJob


bool ana::Muchecks::IsInVolume(double x, double y, double z, double eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpc0{cryoid, 0}, tpcN{cryoid, asGeo->NTPC(cryoid)-1};
    return (x-eps > asGeo->TPC(tpc0).MinX() and x+eps < asGeo->TPC(tpcN).MaxX() and
            y-eps > asGeo->TPC(tpc0).MinY() and y+eps < asGeo->TPC(tpcN).MaxY() and
            z-eps > asGeo->TPC(tpc0).MinZ() and z+eps < asGeo->TPC(tpcN).MaxZ());   
}
bool ana::Muchecks::IsInVolume(float x, float y, float z, float eps) {
    return IsInVolume(double(x), double(y), double(z), double(eps));
}
bool ana::Muchecks::IsInVolume(TLorentzVector const& V, float eps) {
    return IsInVolume(V.X(), V.Y(), V.Z(), double(eps));
}
bool ana::Muchecks::IsInVolume(geo::Point_t const& P) {
    return IsInVolume(P.X(), P.Y(), P.Z(), 0.);
}
bool ana::Muchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}


int ana::Muchecks::GetSection(unsigned int ch) {
    return 4.*ch / asGeo->Nchannels();
}
int ana::Muchecks::GetPlane(unsigned int ch, int sec) {
    return (12.*ch / asGeo->Nchannels() - 3*sec);
}
geo::WireID ana::Muchecks::GetWireID(geo::Point_t const& P, int plane) {
    geo::TPCID tpc = asGeo->FindTPCAtPosition(P);
    if (!tpc.isValid) return geo::WireID();
    geo::WireID wire;
    try {
        wire = asGeo->NearestWireID(P, geo::PlaneID(tpc,2));
    } catch (geo::InvalidWireError const& e) {
        return e.suggestedWireID();
    }
    return wire;
}
raw::ChannelID_t ana::Muchecks::GetChannel(geo::Point_t const& P, int plane) {
    geo::WireID wire = GetWireID(P, plane);
    if (!wire.isValid) return raw::InvalidChannelID;
    return asGeo->PlaneWireToChannel(wire);
}


bool ana::Muchecks::Logging(bool cond, int tab, string msg, string succ, string fail) {
    if (iLogLevel >= iFlagDetails) {
        for (int i=0; i<tab; i++) cout << "\t";
        cout << msg << " ";
        if (!cond) cout << "\033[92m" << succ << "\033[0m" << endl;
        else cout << "\033[91m" << fail << "\033[0m" << endl;
    }
    return cond;
}
bool ana::Muchecks::LogFail(bool cond, int tab, string msg, string fail) {
    if (iLogLevel >= iFlagDetails) {
        if (cond) {
            for (int i=0; i<tab; i++) cout << "\t";
            cout << msg << " ";
            cout << "\033[91m" << fail << "\033[0m" << endl;
        }
    }
    return cond;
}
bool ana::Muchecks::LogSucc(bool cond, int tab, string msg, string succ) {
    if (iLogLevel >= iFlagDetails) {
        if (!cond) {
            for (int i=0; i<tab; i++) cout << "\t";
            cout << msg << " ";
            cout << "\033[92m" << succ << "\033[0m" << endl;
        }
    }
    return cond;
}
    

// Upper triangle of a symmetric matrix
int ana::Muchecks::iCorr(int n, int i, int j) {
    return j - 1 + i*n - (i*(i+3))/2;
}

DEFINE_ART_MODULE(ana::Muchecks)
