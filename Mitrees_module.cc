////////////////////////////////////////////////////////////////////////
// Class:       Mitrees
// Plugin Type: analyzer
// File:        Mitrees_module.cc
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
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEllipse.h"

// std includes
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::min;
using std::max;


namespace ana {
  class Mitrees;
  struct Binning {
        int n;
        double min, max;
  };
}


class ana::Mitrees : public art::EDAnalyzer {
public:
    explicit Mitrees(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Mitrees(Mitrees const&) = delete;
    Mitrees(Mitrees&&) = delete;
    Mitrees& operator=(Mitrees const&) = delete;
    Mitrees& operator=(Mitrees&&) = delete;

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
    // float fDriftVelocity = 0.16; // cm/µs
    float fSamplingRate;

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

    // Input Variables
    float fMichelTimeRadius, // in µs
          fMichelSpaceRadius, // in cm
          fTrackLengthCut; // in cm
    
    // Output Variables
    TTree* tTree;
    size_t iNMichel;
    vector<bool> vbIsPositron;
    vector<float> vfTrueE, vfRecoADC, vfFullRecoADC, vfSphereADC;

    size_t iNHit, iNHitMCP, iNHitMCPdau, iNHitSphere;
    vector<unsigned int> viChan, viChanMCP, viChanMCPdau, viChanSphere;
    vector<float> vfTick, vfTickMCP, vfTickMCPdau, vfTickSphere;
    vector<float> vfADC;

    // Functions
    void reset();

    geo::WireID GetWireID(geo::Point_t const& P, geo::View_t plane);
    raw::ChannelID_t GetChannel(geo::Point_t const& P, geo::View_t plane);

    bool IsInVolume(double x, double y, double z, double eps);
    bool IsInVolume(float x, float y, float z, float eps);
    bool IsInVolume(TLorentzVector const& V, float eps);
    bool IsInVolume(geo::Point_t const& P);

    bool IsUpright(recob::Track const& T);

    int GetSection(raw::ChannelID_t ch);
    geo::View_t GetPlane(raw::ChannelID_t ch, int sec);

    bool Logging(bool cond, int tab, string msg, string succ, string fail);

    int iCorr(int n, int i, int j);
};


ana::Mitrees::Mitrees(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<vector<vector<string>>>("Products")),

    fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("LengthCut")) // in cm
{

    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    fSamplingRate = detinfo::sampling_rate(clockData);

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


    tTree = tfs->make<TTree>("Tree","");
    tTree->Branch("NMichel", &iNMichel);
    tTree->Branch("IsPositron", &vbIsPositron);

    tTree->Branch("TrueE", &vfTrueE);
    tTree->Branch("RecoADC", &vfRecoADC);
    tTree->Branch("FullRecoADC", &vfFullRecoADC);
    tTree->Branch("SphereADC", &vfSphereADC);

    tTree->Branch("NHit", &iNHit);
    tTree->Branch("NHitMCP", &iNHitMCP);
    tTree->Branch("NHitMCPdau", &iNHitMCPdau);
    tTree->Branch("NHitSphere", &iNHitSphere);

    tTree->Branch("Chan", &viChan);
    tTree->Branch("ChanMCP", &viChanMCP);
    tTree->Branch("ChanMCPdau", &viChanMCPdau);
    tTree->Branch("ChanSphere", &viChanSphere);

    tTree->Branch("Tick", &vfTick);
    tTree->Branch("TickMCP", &vfTickMCP);
    tTree->Branch("TickMCPdau", &vfTickMCPdau);
    tTree->Branch("TickSphere", &vfTickSphere);

    tTree->Branch("ADC", &vfADC);
}

void ana::Mitrees::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    // auto const detProp = asDetProp->DataFor(e,clockData);

    reset();

    if (iLogLevel >= iFlagDetails) cout << "evt#" << e.id().event() << "\r" << flush;

    auto const vh_trk = e.getValidHandle<vector<recob::Track>>(tag_trk);
    vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);
        
    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);

    vector<raw::ChannelID_t> vch_mu_end;
    vector<float> vf_mu_end_tick;

    if (iLogLevel >= iFlagDetails) cout << "\tlooping over " << vp_trk.size() << " tracks..." << endl;
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        if (iLogLevel >= iFlagDetails) cout << "\ttrk#" << p_trk->ID()+1 << "\r" << flush;

        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());
        if (Logging(
            !mcp, 
            2, "trk to mcp...", "done", "failed")
        ) continue;    

        if (Logging(
            abs(mcp->PdgCode()) != 13,
            2, "is muon...", "yes", "no")
        ) continue;

        if (Logging(
            mcp->EndProcess() != "Decay",
            2, "is decaying...", "yes", "no")
        ) continue;

        int anti = mcp->PdgCode() < 0;
        if (Logging(
            mcp->NumberDaughters() == 0,
            2, Form("looping over \003[93mmu%s\003[0m daughters...",anti ? "+" : "-"), "", "none")
        ) continue;

        int i_dau = mcp->NumberDaughters() - 1;
        for (; i_dau >= 0; i_dau--) {
            if (iLogLevel >= iFlagDetails) cout << "\t\tdau#" << i_dau+1 << "\r" << flush;
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            if (Logging(
                !mcp_dau,
                3, "id to mcp...", "done", "failed")
            ) continue;

            if (Logging(
                abs(mcp_dau->PdgCode()) != 11 || mcp_dau->Process() != "Decay",
                3, "is michel...", "yes", "no")
            ) continue;

            if (Logging(
                !IsInVolume(mcp_dau->Position(0), fMichelSpaceRadius),
                3, "is inside...", "yes", "no")
            ) i_dau = -1;
            break;
        } // end loop over muon daughters
        if (i_dau == -1) continue;

        iNMichel++;

        simb::MCParticle const * mcp_mich = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));

        vbIsPositron.push_back(mcp_mich->PdgCode() < 0);
        vfTrueE.push_back((mcp_mich->E() - mcp_mich->Mass())*1e3);
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tTrueE: " << "\033[93m" << vfTrueE.back() << " MeV" << "\033[0m" << endl;

        vector<const recob::Hit*> hits_michel = truthUtil.GetMCParticleHits(clockData, *mcp_mich, e, tag_hit.label(), false);

        float RecoADC = 0;
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tlooping over michel's " << hits_michel.size() << " hits...";
        for (recob::Hit const* hit : hits_michel) {

            if (hit->View() != geo::kW) continue;
            // RecoADC += hit->HitSummedADC(); 
            RecoADC += hit->Integral();
            iNHitMCP++;
            viChanMCP.push_back(hit->Channel());
            vfTickMCP.push_back(hit->PeakTime());
        } // end loop over michel hits
        if (iLogLevel >= iFlagDetails) cout << "\033[92m" << " done" << "\033[0m" << endl;

        vfRecoADC.push_back(RecoADC*fADCtoE);
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tRecoADC: " << "\033[93m" << RecoADC << " ~ " << vfRecoADC.back() << " MeV" << "\033[0m" << endl;

        float FullRecoADC = RecoADC;
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tlooping over michel's " << mcp_mich->NumberDaughters() << " daughters..." << endl;
        for (int i_dau=0; i_dau< mcp_mich->NumberDaughters(); i_dau++) {
            if (iLogLevel >= iFlagDetails) cout << "\t\t\tdau#" << i_dau+1 << "\r" << flush;

            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp_mich->Daughter(i_dau));

            if (Logging(
                !mcp_dau,
                4, "id to mcp...", "done", "failed")
            ) continue;

            const int pdg = mcp_dau->PdgCode();
            if (Logging(
                abs(pdg) != 11 && abs(pdg) != 22,
                4, "is electron or photon...",
                Form("%s", abs(pdg)==22 ? "photon" : (pdg>0 ? "elec" : "posi")), "no")
            ) continue;
            
            vector<const recob::Hit*> hits_dau = truthUtil.GetMCParticleHits(clockData, *mcp_dau, e, tag_hit.label(), false);

            if (iLogLevel >= iFlagSpecial && abs(pdg)==22 && hits_dau.size() > 0) cout << "\t\t\t\tphoton with " << hits_dau.size() << " hits" << endl;

            if (iLogLevel >= iFlagDetails) cout << "\t\t\t\tlooping over " << hits_dau.size() << " hits...";
            for (recob::Hit const* hit : hits_dau) {

                if (hit->View() != geo::kW) continue;
                // FullRecoADC += hit->HitSummedADC();
                FullRecoADC += hit->Integral();
                iNHitMCPdau++;
                viChanMCPdau.push_back(hit->Channel());
                vfTickMCPdau.push_back(hit->PeakTime());
            } // end loop over electron hits
            if (iLogLevel >= iFlagDetails) cout << "\033[92m" << " done" << "\033[0m" << endl;
        } // end loop over michel daughters

        vfFullRecoADC.push_back(FullRecoADC*fADCtoE);
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tFullRecoADC: " << "\033[93m" << FullRecoADC << " ~ " << vfFullRecoADC.back() << " MeV" << "\033[0m" << endl;

        geo::Point_t trk_end_pt, trk_start_pt;
        if (IsUpright(*p_trk)) {
            trk_end_pt = p_trk->End();
            trk_start_pt = p_trk->Start();
        } else {
            trk_end_pt = p_trk->Start();
            trk_start_pt = p_trk->End();
        }

        vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk2hit.at(p_trk.key());
        float trk_min_tick=0, trk_max_tick=0;
        for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
            trk_min_tick = min(trk_min_tick, p_hit->PeakTime());
            trk_max_tick = max(trk_max_tick, p_hit->PeakTime());
        }

        // float trk_end_tick, trk_start_tick;
        float trk_end_tick;
        if (abs(trk_end_pt.X()) < abs(trk_start_pt.X())) {
            // trk_start_tick = trk_min_tick;
            trk_end_tick = trk_max_tick;
        } else {
            // trk_start_tick = trk_max_tick;
            trk_end_tick = trk_min_tick;
        }

        // raw::ChannelID_t ch_mu_mcp_end = GetChannel(geo::Point_t(mcp_mich->EndPosition().Vect()), geo::kW);
        raw::ChannelID_t ch_mu_trk_end = GetChannel(trk_end_pt, geo::kW);


        // if(Logging(
        //     ch_mu_trk_end == raw::InvalidChannelID || ch_mu_mcp_end == raw::InvalidChannelID,
        //     2, "mu end channel is valid...", "yes", "no")
        // ) {vch_mu_end.push_back(raw::InvalidChannelID); continue;}

        // if (iLogLevel >= iFlagDetails) cout << "\033[93m" << "\t\tMuon Track End Channel: " << "\033[0m" << ch_mu_trk_end << endl
        //                     << "\033[93m" << "\t\tMichel End Channel: " << "\033[0m" << ch_mu_mcp_end << endl;

        vch_mu_end.push_back(ch_mu_trk_end);
        vf_mu_end_tick.push_back(trk_end_tick);
    } // end loop over tracks

    vfSphereADC = vector<float>(vch_mu_end.size(), 0);

    auto const & vh_hit = e.getValidHandle<vector<recob::Hit>>(tag_hit);
    vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);

    if (iLogLevel >= iFlagDetails) cout << "looping over all hits...";
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {

        if (p_hit->View() != geo::kW) continue;

        iNHit++;
        viChan.push_back(p_hit->Channel());
        vfTick.push_back(p_hit->PeakTime());

        vector<art::Ptr<recob::Track>> vp_trk = fmp_hit2trk.at(p_hit.key());

        bool BelongsToLongTrack = false;
        for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
            BelongsToLongTrack = p_trk->Length() > fTrackLengthCut;
            if (BelongsToLongTrack) break;
        }
        if (BelongsToLongTrack) continue;

        for (unsigned int i=0; i<vch_mu_end.size(); i++) {

            if (vch_mu_end[i] == raw::InvalidChannelID) continue;

            if (abs(int(p_hit->Channel() - vch_mu_end[i])) > fMichelSpaceRadius / fChannelPitch) continue;
            if (abs(int(p_hit->PeakTime() - vf_mu_end_tick[i])) > fMichelTimeRadius / fSamplingRate) continue;

            // vfSphereADC[i] += p_hit->HitSummedADC() * fADCtoE;
            vfSphereADC[i] += p_hit->Integral() * fADCtoE;
            iNHitSphere++;
            viChanSphere.push_back(p_hit->Channel());
            vfTickSphere.push_back(p_hit->PeakTime());
            vfADC.push_back(p_hit->HitSummedADC());
        } // end loop over muon ends
    } // end loop over hits
    if (iLogLevel >= iFlagDetails) cout << "\033[92m" << " done" << "\033[0m" << endl;

    tTree->Fill();

} // end analyze

void ana::Mitrees::beginJob()
{

} // end beginJob


void ana::Mitrees::endJob()
{

    // if (iLogLevel >= iFlagDetails) cout << "\033[93m" << "Mitrees::endJob: Plotting section =============================================" << "\033[0m" << endl;

    // if (iLogLevel >= iFlagDetails) cout << "\033[93m" << "End of Mitrees::endJob ========================================================" << "\033[0m" << endl;
} // end endJob


void ana::Mitrees::reset()
{
    if (iLogLevel >= iFlagDetails) cout << "resetting...";

    iNMichel = 0;

    vbIsPositron.clear();
    vfTrueE.clear();
    vfRecoADC.clear();
    vfFullRecoADC.clear();
    vfSphereADC.clear();

    iNHit = 0;
    iNHitMCP = 0;
    iNHitMCPdau = 0;
    iNHitSphere = 0;

    viChan.clear();
    viChanMCP.clear();
    viChanMCPdau.clear();
    viChanSphere.clear();

    vfTick.clear();
    vfTickMCP.clear();
    vfTickMCPdau.clear();
    vfTickSphere.clear();

    if (iLogLevel >= iFlagDetails) cout << "\033[92m" << " done" << "\033[0m" << endl;
}



bool ana::Mitrees::IsInVolume(double x, double y, double z, double eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpcid0{cryoid, 0}, tpcidN{cryoid, asGeo->NTPC(cryoid)-1};
    return (x-eps > asGeo->TPC(tpcid0).MinX() and x+eps < asGeo->TPC(tpcidN).MaxX() and
            y-eps > asGeo->TPC(tpcid0).MinY() and y+eps < asGeo->TPC(tpcidN).MaxY() and
            z-eps > asGeo->TPC(tpcid0).MinZ() and z+eps < asGeo->TPC(tpcidN).MaxZ());   
}
bool ana::Mitrees::IsInVolume(float x, float y, float z, float eps) {
    return IsInVolume(double(x), double(y), double(z), double(eps));
}
bool ana::Mitrees::IsInVolume(TLorentzVector const& V, float eps) {
    return IsInVolume(V.X(), V.Y(), V.Z(), double(eps));
}
bool ana::Mitrees::IsInVolume(geo::Point_t const& P) {
    return IsInVolume(P.X(), P.Y(), P.Z(), 0.);
}
bool ana::Mitrees::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}


int ana::Mitrees::GetSection(raw::ChannelID_t ch) {
    return int(4.*ch / asWire->Nchannels());
}
geo::View_t ana::Mitrees::GetPlane(raw::ChannelID_t ch, int sec) {
    return static_cast<geo::View_t>(12.*ch / asWire->Nchannels() - 3*sec);
}
geo::WireID ana::Mitrees::GetWireID(geo::Point_t const& P, geo::View_t plane) {
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
raw::ChannelID_t ana::Mitrees::GetChannel(geo::Point_t const& P, geo::View_t plane) {
    geo::WireID wireid = GetWireID(P, plane);
    if (!wireid.isValid) return raw::InvalidChannelID;
    return asWire->PlaneWireToChannel(wireid);
}


bool ana::Mitrees::Logging(bool cond, int tab, string msg, string succ, string fail) {
    if (iLogLevel >= iFlagDetails) {
        for (int i=0; i<tab; i++) cout << "\t";
        cout << msg << " ";
        if (!cond) cout << "\033[92m" << succ << "\033[0m" << endl;
        else cout << "\033[91m" << fail << "\033[0m" << endl;
    }
    return cond;
}

// Upper triangle of a symmetric matrix
int ana::Mitrees::iCorr(int n, int i, int j) {
    return j - 1 + i*n - (i*(i+3))/2;
}

DEFINE_ART_MODULE(ana::Mitrees)
