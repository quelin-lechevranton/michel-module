#include "utils.h"

namespace ana {
    class TrackAnalysis;
}

class ana::TrackAnalysis : 
    public art::EDAnalyzer, 
    private ana::MichelAnalyzer 
{
public:
    explicit TrackAnalysis(fhicl::ParameterSet const& p);
    TrackAnalysis(TrackAnalysis const&) = delete;
    TrackAnalysis(TrackAnalysis&&) = delete;
    TrackAnalysis& operator=(TrackAnalysis const&) = delete;
    TrackAnalysis& operator=(TrackAnalysis&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    ana::Bounds<float> geoTickWindow;
    ana::Bounds3D<float> geoHighX, geoLowX;
    float geoCathodeGap; // cm

    bool fLog;

    struct {
        TTree* tree;
        ana::Hits hits;
        size_t  run, 
                subrun, 
                event; 
        bool is_data;
    } ev;

    struct {
        TTree* tree;
        ana::Hits hits, sec_crossing_hits;
        float length;
        ana::Hit start_hit, end_hit, top_last_hit, bottom_first_hit;

        bool    sh_error,
                is_upright,
                is_cathode_crossing,
                is_anode_crossing,
                is_section_jumping,
                is_section_misaligned,
                is_cathode_misaligned;
    } mu;
};

ana::TrackAnalysis::TrackAnalysis(fhicl::ParameterSet const& p) : 
    EDAnalyzer{p}, 
    MichelAnalyzer{p},
    fLog(p.get<bool>("Log", true))
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
    geoTickWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(),
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        default: break;
    }
    geoCathodeGap = geoHighX.x.min - geoLowX.x.max;

    //
    ev.tree = asFile->make<TTree>("event", "");

    ev.tree->Branch("run", &ev.run);
    ev.tree->Branch("subrun", &ev.subrun);
    ev.tree->Branch("event", &ev.event);
    ev.tree->Branch("is_data", &ev.is_data);

    //
    mu.tree = asFile->make<TTree>("muon", "");

    mu.tree->Branch("length", &mu.length);
    mu.tree->Branch("sh_error", &mu.sh_error);
    mu.tree->Branch("upright", &mu.is_upright);
    mu.tree->Branch("cathode_crossing", &mu.is_cathode_crossing);
    mu.tree->Branch("anode_crossing", &mu.is_anode_crossing);
    mu.tree->Branch("section_jumping", &mu.is_section_jumping);
    mu.tree->Branch("section_misaligned", &mu.is_section_misaligned);
    mu.tree->Branch("cathode_misaligned", &mu.is_cathode_misaligned);
    mu.start_hit.SetBranches(mu.tree, "start_");
    mu.end_hit.SetBranches(mu.tree, "end_");
    mu.hits.SetBranches(mu.tree);
}

void ana::TrackAnalysis::beginJob() {}
void ana::TrackAnalysis::endJob() {}

void ana::TrackAnalysis::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Track handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrTrk vpt_ev;
    art::fill_ptr_vector(vpt_ev, vh_trk);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    ev.run = e.run();
    ev.subrun = e.subRun();
    ev.event = e.event();
    ev.is_data = e.isRealData();

    ev.hits.clear();
    for (PtrHit ph : vph_ev)
        if (ph->View() == geo::kW)
            ev.hits.push_back(GetHit(ph));

    for (PtrTrk pt_ev : vpt_ev) {
        VecPtrHit vph_trk = fmp_trk2hit.at(pt_ev.key());
        ASSERT(vph_trk.size())

        mu.is_upright = IsUpright(*pt_ev);
        geo::Point_t Start = mu.is_upright ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = mu.is_upright ? pt_ev->End() : pt_ev->Start();

        ana::SortedHits sh = GetSortedHits(vph_trk, End.Z() > Start.Z() ? 1 : -1);

        mu.hits.clear();
        mu.sec_crossing_hits.clear();
        mu.is_section_jumping = false;
        mu.is_section_misaligned = false;
        mu.sh_error = !sh;
        LOG(!mu.sh_error);
        if (!mu.sh_error) {
            mu.start_hit = GetHit(sh.start);
            mu.end_hit = GetHit(sh.end);
            if (sh.is_cc()) {
                mu.top_last_hit = GetHit(sh.cc.first);
                mu.bottom_first_hit = GetHit(sh.cc.second);
            } else {
                mu.top_last_hit = ana::Hit{};
                mu.bottom_first_hit = ana::Hit{};
            }

            for (PtrHit ph : sh.vph)
                mu.hits.push_back(GetHit(ph));

            for (PtrHit ph : sh.sc)
                mu.sec_crossing_hits.push_back(GetHit(ph));

            mu.is_cathode_crossing = sh.is_cc();
            mu.is_cathode_misaligned = sh.is_cc() 
                && !(abs(sh.cc.first->PeakTime()-sh.cc.second->PeakTime())*fTick2cm < 3 * geoCathodeGap);
            LOG(mu.is_cathode_crossing);
            LOG(mu.is_cathode_misaligned);

            if (geoDet == kPDVD)
                mu.is_anode_crossing = mu.start_hit.section < 4
                    && geoHighX.z.isInside(mu.start_hit.space, 10.F)
                    && geoTickWindow.isInside(mu.start_hit.tick, 10.F/fTick2cm);
            else if (geoDet == kPDHD)
                mu.is_anode_crossing =
                    geoHighX.z.isInside(mu.start_hit.space, 10.F)
                    && geoTickWindow.isInside(mu.start_hit.tick, 10.F/fTick2cm);
            LOG(mu.is_anode_crossing);

            LOG(sh.secs.size() > 1);
            if (sh.secs.size() > 1) {
                VecPtrHit::iterator sc_it = sh.sc.begin();
                for (unsigned i=0; i<sh.secs.size()-1; i++) {
                    int sec_curr = sh.secs[i];
                    int sec_next = sh.secs[i+1];

                    if (ana::sec2side[geoDet][sec_curr] == ana::sec2side[geoDet][sec_next]) {
                        if (abs(sec_curr - sec_next) != 1)
                            mu.is_section_jumping = true;

                        if ((*sc_it)->PeakTime()-(*++sc_it)->PeakTime())*fTick2cm > 2.F))
                            mu.is_section_misaligned = true;
                    }
                }
                LOG(mu.is_section_jumping);
                LOG(mu.is_section_misaligned);
            }    
        } else {
            mu.start_hit = ana::Hit{};
            mu.end_hit = ana::Hit{};
            mu.top_last_hit = ana::Hit{};
            mu.bottom_first_hit = ana::Hit{};
            mu.is_cathode_crossing = false;
            mu.is_anode_crossing = false;
            mu.is_section_jumping = false;
            mu.is_section_misaligned = false;
            mu.is_cathode_misaligned = false;
        }
    }
}

DEFINE_ART_MODULE(ana::TrackAnalysis)