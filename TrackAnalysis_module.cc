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
    float misalignmentTolerance = 3.F; // cm

    bool fLog;

    struct {
        TTree* tree;
        size_t index=0;
        size_t muon_number=0;
        ana::Hits hits;
        size_t  run, 
                subrun, 
                event; 
        bool is_data;
    } ev;

    struct {
        TTree* tree;
        size_t index=0;
        ana::Hits hits, sec_crossing_hits;
        float length;
        float max_consecutive_dist;
        ana::Point start_point, end_point;
        ana::Hit start_hit, end_hit, top_last_hit, bottom_first_hit;
        ana::LinearRegression top_reg, bot_reg;
        std::vector<float> top_dQdx, bot_dQdx;

        bool    sh_error,
                is_cathode_crossing,
                is_anode_crossing,
                is_section_jumping,
                is_section_misaligned,
                is_cathode_misaligned;

        struct {
            int pdg;
            std::string end_process;
            ana::Hits hits, sec_crossing_hits;
            float max_consecutive_dist;
            ana::Point start_point, end_point;
            ana::Hit start_hit, end_hit, top_last_hit, bottom_first_hit;
            float end_hit_x;
            int dir_z;
            ana::LinearRegression top_reg, bot_reg;
            std::vector<float> top_dQdx, bot_dQdx;
            float end_energy;

            bool    sh_error,
                    is_cathode_crossing,
                    is_anode_crossing,
                    is_section_jumping,
                    is_section_misaligned,
                    is_cathode_misaligned;
            bool has_michel;
            float michel_energy;
        } tru;
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

    ev.tree->Branch("Run",      &ev.run);
    ev.tree->Branch("Subrun",   &ev.subrun);
    ev.tree->Branch("Event",    &ev.event);
    ev.tree->Branch("iEvent",   &ev.index);
    ev.tree->Branch("NMuon",    &ev.muon_number);
    ev.tree->Branch("IsData",   &ev.is_data);
    ev.hits.SetBranches(ev.tree);

    //
    mu.tree = asFile->make<TTree>("muon", "");

    mu.tree->Branch("Run",                  &ev.run);
    mu.tree->Branch("Subrun",               &ev.subrun);
    mu.tree->Branch("Event",                &ev.event);
    mu.tree->Branch("iEvent",               &ev.index);
    mu.tree->Branch("iMuon",                &mu.index);
    mu.tree->Branch("Length",               &mu.length);
    mu.tree->Branch("SortedHitsError",      &mu.sh_error);
    mu.tree->Branch("MaxConsecutiveDist",   &mu.max_consecutive_dist);
    mu.tree->Branch("CathodeCrossing",      &mu.is_cathode_crossing);
    mu.tree->Branch("AnodeCrossing",        &mu.is_anode_crossing);
    mu.tree->Branch("SectionJumping",       &mu.is_section_jumping);
    mu.tree->Branch("SectionMisaligned",    &mu.is_section_misaligned);
    mu.tree->Branch("CathodeMisaligned",    &mu.is_cathode_misaligned);
    mu.start_point.             SetBranches(mu.tree, "Start");
    mu.end_point.               SetBranches(mu.tree, "End");
    mu.start_hit.               SetBranches(mu.tree, "Start");
    mu.end_hit.                 SetBranches(mu.tree, "End");
    mu.top_last_hit.            SetBranches(mu.tree, "TopLast");
    mu.bottom_first_hit.        SetBranches(mu.tree, "BottomFirst");
    mu.hits.                    SetBranches(mu.tree);
    mu.sec_crossing_hits.       SetBranches(mu.tree, "SecCross");
    mu.top_reg.                 SetBranches(mu.tree, "Top");
    mu.bot_reg.                 SetBranches(mu.tree, "Bot");
    mu.tree->Branch("TopdQdx",                &mu.top_dQdx);
    mu.tree->Branch("BotdQdx",                &mu.bot_dQdx);

    mu.tree->Branch("TruPdg",                 &mu.tru.pdg);
    mu.tree->Branch("TruEndProcess",          &mu.tru.end_process);
    mu.tree->Branch("TruEndEnergy",           &mu.tru.end_energy);
    mu.tree->Branch("truSortedHitsError",     &mu.tru.sh_error);
    mu.tree->Branch("TruMaxConsecutiveDist",  &mu.tru.max_consecutive_dist);
    mu.tree->Branch("TruCathodeCrossing",     &mu.tru.is_cathode_crossing);
    mu.tree->Branch("TruAnodeCrossing",       &mu.tru.is_anode_crossing);
    mu.tree->Branch("TruSectionJumping",      &mu.tru.is_section_jumping);
    mu.tree->Branch("TruSectionMisaligned",   &mu.tru.is_section_misaligned);
    mu.tree->Branch("TruCathodeMisaligned",   &mu.tru.is_cathode_misaligned);
    mu.tru.start_point.         SetBranches(mu.tree, "TruStart");
    mu.tru.end_point.           SetBranches(mu.tree, "TruEnd");
    mu.tru.start_hit.           SetBranches(mu.tree, "TruStart");
    mu.tru.end_hit.             SetBranches(mu.tree, "TruEnd");
    mu.tru.top_last_hit.        SetBranches(mu.tree, "TruTopLast");
    mu.tru.bottom_first_hit.    SetBranches(mu.tree, "TruBottomFirst");
    mu.tru.hits.                SetBranches(mu.tree, "Tru");
    mu.tru.sec_crossing_hits.   SetBranches(mu.tree, "TruSecCross");
    mu.tru.top_reg.             SetBranches(mu.tree, "TruTop");
    mu.tru.bot_reg.             SetBranches(mu.tree, "TruBot");
    mu.tree->Branch("TruTopdQdx",          &mu.tru.top_dQdx);
    mu.tree->Branch("TruBotdQdx",          &mu.tru.bot_dQdx);
    mu.tree->Branch("TruHasMichel",           &mu.tru.has_michel);
    mu.tree->Branch("TruMichelEnergy",        &mu.tru.michel_energy);
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

        mu.length = pt_ev->Length();
        bool up = IsUpright(*pt_ev);
        mu.start_point = ana::Point(up ? pt_ev->Start() : pt_ev->End());
        mu.end_point = ana::Point(up ? pt_ev->End() : pt_ev->Start());

        ana::SortedHits sh = GetSortedHits(vph_trk, mu.end_point.z > mu.end_point.z ? 1 : -1);

        mu.hits.clear();
        mu.sec_crossing_hits.clear();
        mu.is_section_jumping = false;
        mu.is_section_misaligned = false;
        mu.top_dQdx.clear();
        mu.bot_dQdx.clear();
        mu.sh_error = !sh;
        LOG(!mu.sh_error);
        if (!mu.sh_error) {

            mu.top_reg = sh.regs[0];
            mu.bot_reg = sh.regs[1];
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

            mu.top_dQdx = GetdQdx(sh.vph.begin(), sh.bot_it(), 6);
            mu.bot_dQdx = GetdQdx(sh.bot_it(), sh.vph.end(), 6);

            mu.max_consecutive_dist = 0.F;
            for (VecPtrHit::iterator it=sh.vph.begin(); it!=sh.vph.end()-1; ++it) {
                Hit h1 = GetHit(*it), h2 = GetHit(*(it+1));
                if (h1.section != h2.section) continue;
                float dist = GetDistance(h1, h2);
                if (dist > mu.max_consecutive_dist)
                    mu.max_consecutive_dist = dist; 
            }

            mu.is_cathode_crossing = sh.is_cc();
            mu.is_cathode_misaligned = sh.is_cc() 
                && abs(sh.cc.first->PeakTime()-sh.cc.second->PeakTime()) > (geoCathodeGap + misalignmentTolerance)/fTick2cm;
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

                        if (((*sc_it)->PeakTime()-(*++sc_it)->PeakTime()) > misalignmentTolerance/fTick2cm)
                            mu.is_section_misaligned = true;
                    }
                }
                LOG(mu.is_section_jumping);
                LOG(mu.is_section_misaligned);
            }    
        } else {
            mu.max_consecutive_dist = -1.F;
            mu.top_reg = ana::LinearRegression{};
            mu.bot_reg = ana::LinearRegression{};
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


        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        if (mcp) {

            mu.tru.pdg = mcp->PdgCode();
            mu.tru.end_process = mcp->EndProcess();
            mu.tru.end_energy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV

            mu.tru.start_point = ana::Point(mcp->Position().Vect());
            mu.tru.end_point = ana::Point(mcp->EndPosition().Vect());

            VecPtrHit vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            ana::SortedHits sh_mcp = GetSortedHits(vph_mcp_mu, mu.tru.end_point.z > mu.tru.start_point.z ? 1 : -1);

            mu.tru.hits.clear();   
            mu.tru.sec_crossing_hits.clear();
            mu.tru.is_section_jumping = false;
            mu.tru.is_section_misaligned = false;
            mu.tru.top_dQdx.clear();
            mu.tru.bot_dQdx.clear();
            mu.tru.sh_error = !sh_mcp;
            LOG(!mu.tru.sh_error);
            if (!mu.tru.sh_error) {

                mu.tru.top_reg = sh_mcp.regs[0];
                mu.tru.bot_reg = sh_mcp.regs[1];
                mu.tru.start_hit = GetHit(sh_mcp.start);
                mu.tru.end_hit = GetHit(sh_mcp.end);
                if (sh_mcp.is_cc()) {
                    mu.tru.top_last_hit = GetHit(sh_mcp.cc.first);
                    mu.tru.bottom_first_hit = GetHit(sh_mcp.cc.second);
                } else {
                    mu.tru.top_last_hit = ana::Hit{};
                    mu.tru.bottom_first_hit = ana::Hit{};
                }

                for (PtrHit ph : sh_mcp.vph)
                    mu.tru.hits.push_back(GetHit(ph));

                for (PtrHit ph : sh_mcp.sc)
                    mu.tru.sec_crossing_hits.push_back(GetHit(ph));

                mu.tru.top_dQdx = GetdQdx(sh_mcp.vph.begin(), sh_mcp.bot_it(), 6);
                mu.tru.bot_dQdx = GetdQdx(sh_mcp.bot_it(), sh_mcp.vph.end(), 6);

                mu.tru.max_consecutive_dist = 0.F;
                for (VecPtrHit::iterator it=sh_mcp.vph.begin(); it!=sh_mcp.vph.end()-1; ++it) {
                    Hit h1 = GetHit(*it), h2 = GetHit(*(it+1));
                    if (h1.section != h2.section) continue;
                    float dist = GetDistance(h1, h2);
                    if (dist > mu.tru.max_consecutive_dist)
                        mu.tru.max_consecutive_dist = dist; 
                }

                mu.tru.is_cathode_crossing = sh_mcp.is_cc();
                mu.tru.is_cathode_misaligned = sh_mcp.is_cc() 
                    && !(abs(sh_mcp.cc.first->PeakTime()-sh_mcp.cc.second->PeakTime())*fTick2cm < 3 * geoCathodeGap);
                LOG(mu.tru.is_cathode_crossing);
                LOG(mu.tru.is_cathode_misaligned);

                if (geoDet == kPDVD)
                    mu.tru.is_anode_crossing = mu.tru.start_hit.section < 4
                        && geoHighX.z.isInside(mu.tru.start_hit.space, 10.F)
                        && geoTickWindow.isInside(mu.tru.start_hit.tick, 10.F/fTick2cm);
                else if (geoDet == kPDHD)
                    mu.tru.is_anode_crossing =
                        geoHighX.z.isInside(mu.tru.start_hit.space, 10.F)
                        && geoTickWindow.isInside(mu.tru.start_hit.tick, 10.F/fTick2cm);
                LOG(mu.tru.is_anode_crossing);

                LOG(sh_mcp.secs.size() > 1);
                if (sh_mcp.secs.size() > 1) {
                    VecPtrHit::iterator sc_it = sh_mcp.sc.begin();
                    for (unsigned i=0; i<sh_mcp.secs.size()-1; i++) {
                        int sec_curr = sh_mcp.secs[i];
                        int sec_next = sh_mcp.secs[i+1];

                        if (ana::sec2side[geoDet][sec_curr] == ana::sec2side[geoDet][sec_next]) {
                            if (abs(sec_curr - sec_next) != 1)
                                mu.tru.is_section_jumping = true;

                            if (((*sc_it)->PeakTime()-(*++sc_it)->PeakTime())*fTick2cm > 2.F)
                                mu.tru.is_section_misaligned = true;
                        }
                    }
                    LOG(mu.tru.is_section_jumping);
                    LOG(mu.tru.is_section_misaligned);
                }    
            } else {
                mu.tru.max_consecutive_dist = -1.F;
                mu.tru.top_reg = ana::LinearRegression{};
                mu.tru.bot_reg = ana::LinearRegression{};
                mu.tru.start_point = ana::Point{};
                mu.tru.end_point = ana::Point{};
                mu.tru.start_hit = ana::Hit{};
                mu.tru.end_hit = ana::Hit{};
                mu.tru.top_last_hit = ana::Hit{};
                mu.tru.bottom_first_hit = ana::Hit{};
                mu.tru.is_cathode_crossing = false;
                mu.tru.is_anode_crossing = false;
                mu.tru.is_section_jumping = false;
                mu.tru.is_section_misaligned = false;
                mu.tru.is_cathode_misaligned = false;
            }

            simb::MCParticle const* mcp_mi = GetMichelMCP(mcp);
            mu.tru.has_michel = !!mcp_mi;
            LOG(mu.tru.has_michel);
            if (mu.tru.has_michel)
                mu.tru.michel_energy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3; // MeV
            else 
                mu.tru.michel_energy = -1.F;
        }


        mu.tree->Fill();
        mu.index++;
        ev.muon_number++;
    }
    ev.tree->Fill();
    ev.index++;
}

DEFINE_ART_MODULE(ana::TrackAnalysis)