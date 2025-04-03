////////////////////////////////////////////////////////////////////////
// Class:       PandoraChecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        PandoraChecks_module.cc
//
// Generated at Mon Mar 24 08:43:32 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <TTree.h>

namespace ana {
    class PandoraChecks;
    struct Point {
        float x, y, z;
        Point(void) : x(0), y(0), z(0) {}
        Point(geo::Point_t const& p) : x(p.X()), y(p.Y()), z(p.Z()) {}
        Point(TVector3 const& p) : x(p.X()), y(p.Y()), z(p.Z()) {}
    };
    struct TrackOut {
        int id;
        Point start;
        Point end;

        TrackOut(void) : id(-999), start(), end() {}
        TrackOut(recob::Track const& trk) : id(trk.ID()), start(trk.Start()), end(trk.End()) {}
        TrackOut(simb::MCParticle const& mcp) : id(mcp.TrackId()), start(mcp.Position().Vect()), end(mcp.EndPosition().Vect()) {}
        void SetBranches(TTree* t, const char* name) {
            t->Branch(Form("%sID", name), &id);
            t->Branch(Form("%sStartX", name), &start.x);
            t->Branch(Form("%sStartY", name), &start.y);
            t->Branch(Form("%sStartZ", name), &start.z);
            t->Branch(Form("%sEndX", name), &end.x);
            t->Branch(Form("%sEndY", name), &end.y);
            t->Branch(Form("%sEndZ", name), &end.z);
        }
    };
    struct HitsOut {
        std::vector<unsigned> channel;
        std::vector<float> tick;

        HitsOut() : channel(), tick() {}
        void push_back(recob::Hit const& hit) {
            this->channel.push_back(hit.Channel());
            this->tick.push_back(hit.PeakTime());
        }
        HitsOut(std::vector<art::Ptr<recob::Hit>> hits) : channel(), tick() {
            for (art::Ptr<recob::Hit> hit : hits) this->push_back(*hit);
        }
        void SetBranches(TTree* t, const char* name) {
            t->Branch(Form("%sChannel", name), &this->channel);
            t->Branch(Form("%sTick", name), &this->tick);
        }
    };
}


class ana::PandoraChecks : public art::EDAnalyzer {
public:
    explicit PandoraChecks(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    PandoraChecks(PandoraChecks const&) = delete;
    PandoraChecks(PandoraChecks&&) = delete;
    PandoraChecks& operator=(PandoraChecks const&) = delete;
    PandoraChecks& operator=(PandoraChecks&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
private:

    const detinfo::DetectorPropertiesService* asDetProp;
    const detinfo::DetectorClocksService* asDetClocks;

    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    art::InputTag tag_trk;
    art::InputTag tag_pdr;
    art::InputTag tag_hit;

    TTree* tTrack;
    TrackOut evt_trk, evt_mcp, mcp_trk, mich_mcp;
    double evt_trk_t0;
    HitsOut mcp_hits, mcp_hits_eve;

    simb::MCParticle const* trk2mcp(art::Ptr<recob::Track>, detinfo::DetectorClocksData, art::FindManyP<recob::Hit>);
    art::Ptr<recob::Track> mcp2trk(simb::MCParticle const*, std::vector<art::Ptr<recob::Track>> const&, detinfo::DetectorClocksData, art::FindManyP<recob::Hit>);
    std::vector<art::Ptr<recob::Hit>> mcp2hits(simb::MCParticle const*, std::vector<art::Ptr<recob::Hit>> const&, detinfo::DetectorClocksData, bool = true);
};


ana::PandoraChecks::PandoraChecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // ,
{
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    tag_trk = art::InputTag("pandoraTrack");
    tag_pdr = art::InputTag("pandora");
    tag_hit = art::InputTag("hitpdune");

    tTrack = tfs->make<TTree>("tracks", "");
    evt_trk.SetBranches(tTrack, "EventTrack");
    tTrack->Branch("EventTrackT0", &evt_trk_t0);
    mich_mcp.SetBranches(tTrack, "MichelMCP");
    evt_mcp.SetBranches(tTrack, "EventMCP");
    mcp_trk.SetBranches(tTrack, "MCPTrack");
    mcp_hits.SetBranches(tTrack, "MCPHits");
    mcp_hits_eve.SetBranches(tTrack, "MCPEveHits");
}

void ana::PandoraChecks::beginJob() {
    sim::ParticleList::AdoptEveIdCalculator(new sim::EmEveIdCalculator);
}

void ana::PandoraChecks::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);

    std::cout << "==== EVENT " << e.event() << "====" << std::endl;

    auto const& vh_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(tag_pdr);
    auto const& vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    art::FindOneP<recob::PFParticle> fop_trk2pfp(vh_trk, e, tag_trk);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    art::FindOneP<anab::T0> fop_pfp2t0(vh_pfp, e, tag_pdr);


    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);


    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        std::cout << "trk.ID(): " << p_trk->ID() << std::endl;

        simb::MCParticle const* mcp = trk2mcp(p_trk, clockData, fmp_trk2hit);
        if (!mcp) {
            std::cout << "\t" "no mcp found" << std::endl;
            continue;
        }


        if (abs(mcp->PdgCode()) != 13) {
            std::cout << "\t" "not a muon" << std::endl;
            continue;
        }

        art::Ptr<recob::PFParticle> p_pfp = fop_trk2pfp.at(p_trk.key());
        art::Ptr<anab::T0> p_t0 = fop_pfp2t0.at(p_pfp.key());

        art::Ptr<recob::Track> p_trk_from_mcp = mcp2trk(mcp, vp_trk, clockData, fmp_trk2hit);
        if (p_trk_from_mcp != p_trk) std::cout << "\033[1;91m" "mcp2trk: ID(): " << p_trk_from_mcp->ID() << " NOT SAME TRACK" "\033[0m" << std::endl;

        std::vector<art::Ptr<recob::Hit>> vp_hit_from_mcp_eve = mcp2hits(mcp, vp_hit, clockData, true);
        std::vector<art::Ptr<recob::Hit>> vp_hit_from_mcp_no_eve = mcp2hits(mcp, vp_hit, clockData, false);

        // michel
        bool has_numu = false, has_nue = false;
        simb::MCParticle const* mcp_michel = nullptr;
        for (int i_dau=mcp->NumberDaughters()-3; i_dau<mcp->NumberDaughters(); i_dau++) {
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    
            if (!mcp_dau) continue;

            switch (abs(mcp_dau->PdgCode())) {
                case 14: has_numu = true; break;
                case 12: has_nue = true; break;
                case 11: mcp_michel = mcp_dau; break;
            }
        }

        evt_trk = TrackOut(*p_trk);
        evt_trk_t0 = p_t0 ? p_t0->Time() : -999;

        evt_mcp = mcp ? TrackOut(*mcp) : TrackOut();
        mcp_trk = p_trk_from_mcp ? TrackOut(*p_trk_from_mcp) : TrackOut();
        mich_mcp = mcp_michel and has_numu and has_nue ? TrackOut(*mcp_michel) : TrackOut();

        mcp_hits_eve = HitsOut(vp_hit_from_mcp_eve);
        mcp_hits = HitsOut(vp_hit_from_mcp_no_eve);

        tTrack->Fill();
    }
}


simb::MCParticle const* ana::PandoraChecks::trk2mcp(art::Ptr<recob::Track> p_trk, detinfo::DetectorClocksData clockData, art::FindManyP<recob::Hit> fmp) {
    std::unordered_map<int, float> map_tid_ene;
    for (art::Ptr<recob::Hit> const& p_hit : fmp.at(p_trk.key()))
        for (sim::TrackIDE ide : bt_serv->HitToTrackIDEs(clockData, p_hit))
            map_tid_ene[ide.trackID] += ide.energy;

    float max_ene = -1;
    int tid_max = 0;
    for (std::pair<int, float> p : map_tid_ene)
        if (p.second > max_ene) max_ene = p.second, tid_max = p.first;
    return max_ene == -1 ? nullptr : pi_serv->TrackIdToParticle_P(tid_max);
}

art::Ptr<recob::Track> ana::PandoraChecks::mcp2trk(simb::MCParticle const* mcp, std::vector<art::Ptr<recob::Track>> const& vp_trk ,detinfo::DetectorClocksData clockData, art::FindManyP<recob::Hit> fmp) {
    std::unordered_map<art::Ptr<recob::Track>, unsigned> map_trk_nhit;
    for (art::Ptr<recob::Track> p_trk : vp_trk)
        map_trk_nhit[p_trk] += bt_serv->TrackIdToHits_Ps(clockData, mcp->TrackId(), fmp.at(p_trk.key())).size();

    unsigned max = 0;
    art::Ptr<recob::Track> p_trk_from_mcp;
    for (std::pair<art::Ptr<recob::Track>, unsigned> p : map_trk_nhit)
        if (p.second > max) max = p.second, p_trk_from_mcp = p.first;
    return p_trk_from_mcp;
}

std::vector<art::Ptr<recob::Hit>> ana::PandoraChecks::mcp2hits(simb::MCParticle const* mcp, std::vector<art::Ptr<recob::Hit>> const& vp_hit, detinfo::DetectorClocksData clockData, bool use_eve) {
    std::vector<art::Ptr<recob::Hit>> vp_hit_from_mcp;
    for (art::Ptr<recob::Hit> p_hit : vp_hit)
        for (sim::TrackIDE ide : (use_eve ? bt_serv->HitToEveTrackIDEs(clockData, p_hit) : bt_serv->HitToTrackIDEs(clockData, p_hit)))
            if (ide.trackID == mcp->TrackId()) vp_hit_from_mcp.push_back(p_hit);
    return vp_hit_from_mcp;
}

DEFINE_ART_MODULE(ana::PandoraChecks)
