////////////////////////////////////////////////////////////////////////
// Class:       BackTrackerTest
// Plugin Type: analyzer (Unknown Unknown)
// File:        BackTrackerTest_module.cc
//
// Generated at Tue Oct 21 03:01:57 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"

namespace ana {
    class BackTrackerTest;
}


class ana::BackTrackerTest : public art::EDAnalyzer, private ana::MichelAnalyzer {
public:
    explicit BackTrackerTest(fhicl::ParameterSet const& p);
    BackTrackerTest(BackTrackerTest const&) = delete;
    BackTrackerTest(BackTrackerTest&&) = delete;
    BackTrackerTest& operator=(BackTrackerTest const&) = delete;
    BackTrackerTest& operator=(BackTrackerTest&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    ana::Bounds3D<float> geoBounds;
};


ana::BackTrackerTest::BackTrackerTest(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, MichelAnalyzer{p}
{
    switch (geoDet) {
        case kPDVD:
            geoBounds = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(), 
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoBounds = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        default: break;
    }
}

void ana::BackTrackerTest::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const& vh_mcp = e.getValidHandle<std::vector<simb::MCParticle>>(tag_mcp);
    auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);

    if (!vh_mcp || !vh_hit) {
        std::cout 
            << "\033[1;91m" "invalid handles" "\033[0m"
            << (vh_mcp?"" :" mcp ") 
            << (vh_hit?"" :" hit ") 
            << std::endl;
        return;
    }

    std::vector<art::Ptr<recob::Hit>> vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    std::cout << "sizes:" << std::endl;
    std::cout << "\tmcp: " << vh_mcp->size() << std::endl;
    std::cout << "\thit: " << vh_hit->size() << std::endl;

    // std::map<int, unsigned> pdg_counts;
    // float mean_muon_nhit=0;
    // unsigned muon_n=0;
    // for (simb::MCParticle const& mcp : *vh_mcp) {
    //     pdg_counts[mcp.PdgCode()]++;
    //     if (abs(mcp.PdgCode()) != 13) continue;
    //     mean_muon_nhit += bt_serv->TrackIdToHits_Ps(clockData, mcp.TrackId(), vph_ev).size();
    //     muon_n++;
    // }
    // mean_muon_nhit /= muon_n;

    // std::cout << "PDG\t\tCOUNT" << std::endl;
    // std::map<unsigned, int> counts_pdg;
    // for (auto const& [pdg, count] : pdg_counts) {
    //     counts_pdg[count] = pdg;
    // }
    // for (auto const& [count, pdg] : counts_pdg) {
    //     std::cout << GetParticleName(pdg) << ": \t\t" << count << std::endl;
    // }  
    // std::cout << "MEAN MUON NHIT: " << mean_muon_nhit << std::endl;

    // float mean_hit_nides=0;
    // for (art::Ptr<recob::Hit> const& ph_ev : vph_ev) {
    //     std::vector<sim::TrackIDE> ides = bt_serv->HitToTrackIDEs(clockData, ph_ev);
    //     mean_hit_nides += ides.size();  
    // }
    // mean_hit_nides /= vph_ev.size();
    // std::cout << "MEAN HIT NIDES: " << mean_hit_nides << std::endl;

    std::cout
        << "evt#" << e.event()
        << std::endl;

    unsigned n_mi=0;
    for (simb::MCParticle const& mcp : *vh_mcp) {
        if (abs(mcp.PdgCode()) != 13) continue;
        if (!geoBounds.isInside(mcp.EndPosition().Vect(), 10.F)) continue;
        simb::MCParticle const* mcp_mi = GetMichelMCP(&mcp);
        if (!mcp_mi) continue;

        std::cout << "\t"
            << "mi#" << n_mi++
            << " (" << GetParticleName(mcp_mi->PdgCode()) << ")"
            << "\ttrackID: " << mcp_mi->TrackId()
            << std::endl;

        VecPtrHit vph_mi = bt_serv->TrackIdToHits_Ps(clockData, mcp_mi->TrackId(), vph_ev);

        VecPtrHit vph2_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, false);
        VecPtrHit vph3_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true);

        std::cout << "\t\tnhits (bt_serv): " << vph_mi.size() << std::endl;
        std::cout << "\t\tnhits (ana no eve): " << vph2_mi.size() << std::endl;
        std::cout << "\t\tnhits (ana with eve): " << vph3_mi.size() << std::endl;


        VecPtrHit unseen_ph2, unseen_ph3;
        for (PtrHit const& ph_mi2 : vph2_mi) { unseen_ph2.push_back(ph_mi2); }
        for (PtrHit const& ph_mi3 : vph3_mi) { unseen_ph3.push_back(ph_mi3); }
        for (PtrHit const& ph_mi : vph_mi) {

            auto it2 = std::find_if(
                unseen_ph2.begin(), unseen_ph2.end(),
                [&ph_mi](PtrHit ph2){ return ph2.key() == ph_mi.key(); }
            );
            if (it2 != unseen_ph2.end()) {
                unseen_ph2.erase(it2);
            } else {
                std::cout << "\t\t\tmissing (bt) in (noeve): hit key " << ph_mi.key() << std::endl;
            }

            auto it3 = std::find_if(
                unseen_ph3.begin(), unseen_ph3.end(),
                [&ph_mi](PtrHit ph3){ return ph3.key() == ph_mi.key(); }
            );
            if (it3 != unseen_ph3.end()) {
                unseen_ph3.erase(it3);
            } else {
                std::cout << "\t\t\tmissing (bt) in (eve): hit key " << ph_mi.key() << std::endl;
            }
        }
        for (PtrHit ph2 : unseen_ph2) {
            std::cout << "\t\t\textra (noeve) hit key " << ph2.key() << "\tADC: " << ph2->HitSummedADC() << " MeV: " << ph2->HitSummedADC()*(200 * 23.6 * 1e-6 / 0.7) << std::endl;
            for(int tid : bt_serv->HitToTrackIds(clockData, *ph2)) {
                std::cout << "\t\t\t\ttrackID: " << tid << std::endl;
            }
            for (sim::TrackIDE tide : bt_serv->HitToTrackIDEs(clockData, ph2)) {
                std::cout << "\t\t\t\tide trackID: " << tide.trackID << "\tenergy: " << tide.energy << std::endl;
            }
        }
        for (PtrHit ph3 : unseen_ph3) {
            std::cout << "\t\t\textra (eve) hit key " << ph3.key() << std::endl;
        }

        // for (PtrHit const& ph_ev : vph_ev) {
        //     std::vector<sim::TrackIDE> tides = bt_serv->HitToEveTrackIDEs(clockData, ph_ev);
        //     for (sim::TrackIDE const& tide : tides) {
        //         tide.
        //     }
        //     std::vector<sim::IDE> ides = bt_serv->HitToAvgSimIDEs(clockData, ph_ev);
        //     for (sim::IDE const& ide : ides) {
        //         ide.
        //     }
        // }

    }
}

void ana::BackTrackerTest::beginJob() {}
void ana::BackTrackerTest::endJob() {}

DEFINE_ART_MODULE(ana::BackTrackerTest)
