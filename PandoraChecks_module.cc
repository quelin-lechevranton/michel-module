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
        void SetBranches(TTree* t, const char* name) const {
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
        void SetBranches(TTree* t, const char* name) const {
            t->Branch(Form("%sChannel"), &this->channel);
            t->Branch(Form("%sTick"), &this->tick);
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
    HitsOut mcp_hits;
};


ana::PandoraChecks::PandoraChecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // ,
{
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    tag_trk = art::InputTag("pandoraTrack");
    tag_pdr = art::InputTag("pandora");
    tag_hit = art::InputTag("hitpdune");

    tTrack = tfs->make<TTree>("tracks");
    evt_trk.SetBranches(tTrack, "EventTrack");
    tTrack->Branch("EventTrackT0", &evt_trk_t0);
    mich_mcp.SetBranches(tTrack, "MichelMCP");
    evt_mcp.SetBranches(tTrack, "EventMCP");
    mcp_trk.SetBranches(tTrack, "MCPTrack");
    mcp_hits.SetBranches(tTrack, "MCPHits");
}

void ana::PandoraChecks::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);

    std::cout << "==== EVENT " << e.event() << "====" << std::endl;

    auto const& vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    auto const& vh_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(tag_pdr);

    auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    art::FindOneP<recob::PFParticle> fop_trk2pfp(vh_trk, e, tag_trk);
    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<anab::T0> fop_pfp2t0(vh_pfp, e, tag_pdr);

    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        std::cout << "trk.ID(): " << p_trk->ID() << std::endl;

        // trk2mcp
        std::unordered_map<int, float> map_tid_ene;
        std::vector<art::Ptr<recob::Hit>> vp_hit = fmp_trk2hit.at(p_trk.key());
        for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
            std::vector<sim::TrackIDE> v_ide = bt_serv->HitToTrackIDEs(clockData, p_hit);
            for (sim::TrackIDE ide : v_ide) {
                map_tid_ene[ide.trackID] += ide.energy;
            }
        }
        float max_ene = -1;
        int tid_max;
        for (std::pair<int, float> p : map_tid_ene) {
            if (p.second > max_ene) {
                max_ene = p.second;
                tid_max = p.first;
            }
        }
        if (max_ene == -1) {
            std::cout << "\t" "no mcp found" << std::endl;
            continue;
        }
        simb::MCParticle const* mcp = pi_serv->TrackIdToParticle_P(tid_max);



        if (abs(mcp->PdgCode()) != 13) {
            std::cout << "\t" "not a muon" << std::endl;
            continue;
        }

        art::Ptr<recob::PFParticle> p_pfp = fop_trk2pfp.at(p_trk.key());
        art::Ptr<anab::T0> p_t0 = fop_pfp2t0.at(p_pfp.key());

        // mcp2trk
        std::unordered_map<art::Ptr<recob::Track>*, float> map_trk_ene;
        for (art::Ptr<recob::Track> p_trk2 : vp_trk) {
            std::vector<art::Ptr<recob::Hit>> vp_hit_from_trk2_mcp = bt_serv->TrackIdToHits_Ps(clockData, mcp->TrackId(), fmp_trk2hit.at(p_trk2.key()));
            for (art::Ptr<recob::Hit> p_hit : vp_hit_from_trk2_mcp) {
                std::vector<sim::IDE const*> v_ide = bt_serv->HitToSimIDEs_Ps(clockData, p_hit);
                for (sim::IDE const* ide : v_ide) {
                    if (ide->trackID == mcp->TrackId()) map_trk_ene[&p_trk2] += ide->energy;
                }
            }
        }
        max_ene = -1;
        art::Ptr<recob::Track> p_trk_from_mcp;
        for (std::pair<art::Ptr<recob::Track>*, float> p : map_trk_ene) {
            if (p.second > max_ene) {
                max_ene = p.second;
                p_trk_from_mcp = *p.first;
            }
        }

        if (p_trk_from_mcp != p_trk) std::cout << "\033[1;91m" "NOT SAME TRACK" "\033[0m" << std::endl;

        printf("ax.plot([%f,%f],[%f,%f],[%f,%f], c=\"red\")\n", p_trk->Start().X(), p_trk->End().X(), p_trk->Start().Y(), p_trk->End().Y(), p_trk->Start().Z(), p_trk->End().Z());
        printf("ax.plot([%f,%f],[%f,%f],[%f,%f], c=\"red\")\n", mcp->Position().X(), mcp->EndPosition().X(), mcp->Position().Y(), mcp->EndPosition().Y(), mcp->Position().Z(), mcp->EndPosition().Z());

        // mcp2hits
        std::vector<art::Ptr<recob::Hit>> vp_hit_from_mcp;
        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            std::vector<sim::TrackIDE> v_ide = bt_serv->HitToTrackIDEs(clockData, p_hit);
            for (sim::TrackIDE ide : v_ide) {
                if (ide.trackID == mcp->TrackId()) vp_hit_from_mcp.push_back(p_hit);
            }
        }


        // michel
        // for (int i=1; i<=3; i++) {
        //     pi_serv->
        // }

        evt_trk = TrackOut(*p_trk);
        evt_trk_t0 = p_t0 ? p_t0->Time() : -999;

        evt_mcp = mcp ? TrackOut(*mcp) : TrackOut();
        mcp_trk = p_trk_from_mcp ? TrackOut(*p_trk_from_mcp) : TrackOut();

        tTrack->Fill();
    }
}

DEFINE_ART_MODULE(ana::PandoraChecks)
