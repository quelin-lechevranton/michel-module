#include "utils.h"

namespace ana {
    class TrackAnalysis;
}

class ana::TrackAnalysis : 
    public art::EDAnalyzer, 
    private ana::MichelModule 
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

    bool inLog, fAssert, fRemoveBadHits;

    struct {
        TTree* tree;
        size_t index=0;
        size_t muon_number;
        ana::Hits hits;
        size_t  run, 
                subrun, 
                event; 
        bool is_data;
    } ev;

    struct {
        TTree* tree;
        size_t index=0;
        ana::Hits hits, 
                  sec_crossing_hits;
        // struct {
        //     ana::Hits hits;
        //     std::vector<unsigned> hit_indices;
        //     std::vector<float> hit_dxs;
        //     ana::Hit end_hit;
        // } trk;
        float length;
        float max_consecutive_dist;
        ana::Point  start_point, 
                    end_point;
        ana::Hit    start_hit,
                    end_hit, 
                    top_last_hit, 
                    bot_first_hit;
        float end_hit_y;
        ana::LinearRegression   top_reg, 
                                bot_reg;
        std::vector<float>  top_dQds, 
                            bot_dQds;
        bool    sh_error,
                cathode_crossing,
                anode_crossing,
                section_jumping,
                section_misaligned,
                cathode_misaligned;

        struct {
            int pdg;
            std::string end_process;
            ana::Hits   hits, 
                        sec_crossing_hits;
            float max_consecutive_dist;
            ana::Point  start_point, 
                        end_point;
            ana::Hit    start_hit, 
                        end_hit, 
                        top_last_hit, 
                        bot_first_hit;
            // float end_hit_x;
            // int dir_z;
            ana::LinearRegression   top_reg, 
                                    bot_reg;
            std::vector<float>  top_dQds, 
                                bot_dQds;
            float end_energy;
            bool    downward,
                    sh_error,
                    cathode_crossing,
                    anode_crossing,
                    section_jumping,
                    section_misaligned,
                    cathode_misaligned;
            bool has_michel;
            float michel_energy;
        } tru;
    } mu;

    void reset_ev(void);
    void reset_mu(void);
};

void ana::TrackAnalysis::reset_ev(void) {
    ev.muon_number = 0;
    ev.hits.clear();
}
void ana::TrackAnalysis::reset_mu(void) {
    // mu.trk.hits.clear();
    // mu.trk.hit_indices.clear();
    // mu.trk.hit_dxs.clear();
    // mu.trk.end_hit = ana::Hit{};

    mu.hits.clear();
    mu.sec_crossing_hits.clear();
    mu.top_dQds.clear();
    mu.bot_dQds.clear();

    mu.sh_error = false;
    mu.cathode_crossing = false;
    mu.anode_crossing = false;
    mu.section_jumping = false;
    mu.section_misaligned = false;
    mu.cathode_misaligned = false;

    mu.max_consecutive_dist = -1.F;
    mu.top_reg = ana::LinearRegression{};
    mu.bot_reg = ana::LinearRegression{};
    mu.start_hit = ana::Hit{};
    mu.end_hit = ana::Hit{};
    mu.end_hit_y = -500.F;
    mu.top_last_hit = ana::Hit{};
    mu.bot_first_hit = ana::Hit{};


    mu.tru.pdg = 0;
    mu.tru.end_process = "";
    mu.tru.end_energy = -1.F;
    mu.tru.hits.clear();
    mu.tru.sec_crossing_hits.clear();
    mu.tru.max_consecutive_dist = -1.F;
    mu.tru.top_dQds.clear();
    mu.tru.bot_dQds.clear();

    mu.tru.top_reg = ana::LinearRegression{};
    mu.tru.bot_reg = ana::LinearRegression{};
    mu.tru.start_point = ana::Point{};
    mu.tru.end_point = ana::Point{};
    mu.tru.start_hit = ana::Hit{};
    mu.tru.end_hit = ana::Hit{};
    mu.tru.top_last_hit = ana::Hit{};
    mu.tru.bot_first_hit = ana::Hit{};
    // mu.tru.end_hit_x = 0.F;
    // mu.tru.dir_z = 0.F;

    mu.tru.downward = false;
    mu.tru.sh_error = false;
    mu.tru.cathode_crossing = false;
    mu.tru.anode_crossing = false;
    mu.tru.section_jumping = false;
    mu.tru.section_misaligned = false;
    mu.tru.cathode_misaligned = false;
    mu.tru.has_michel = false;
    mu.tru.michel_energy = -1.F;
}

ana::TrackAnalysis::TrackAnalysis(fhicl::ParameterSet const& p) : 
    EDAnalyzer{p}, 
    MichelModule{p},
    inLog(p.get<bool>("Log", true)),
    fAssert(p.get<bool>("Assert", false)),
    fRemoveBadHits(p.get<bool>("RemoveBadHits", false))
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
    mu.tree->Branch("iMuonInEvent",         &ev.muon_number);
    mu.tree->Branch("iMuon",                &mu.index);
    mu.tree->Branch("Length",               &mu.length);
    mu.tree->Branch("SortedHitsError",      &mu.sh_error);
    mu.tree->Branch("MaxConsecutiveDist",   &mu.max_consecutive_dist);
    mu.tree->Branch("CathodeCrossing",      &mu.cathode_crossing);
    mu.tree->Branch("AnodeCrossing",        &mu.anode_crossing);
    mu.tree->Branch("SectionJumping",       &mu.section_jumping);
    mu.tree->Branch("SectionMisaligned",    &mu.section_misaligned);
    mu.tree->Branch("CathodeMisaligned",    &mu.cathode_misaligned);
    mu.start_point.             SetBranches(mu.tree, "Start");
    mu.end_point.               SetBranches(mu.tree, "End");
    mu.start_hit.               SetBranches(mu.tree, "Start");
    mu.end_hit.                 SetBranches(mu.tree, "End");
    mu.tree->Branch("EndHitY",             &mu.end_hit_y);
    mu.top_last_hit.            SetBranches(mu.tree, "TopLast");
    mu.bot_first_hit.           SetBranches(mu.tree, "BotFirst");
    mu.hits.                    SetBranches(mu.tree, "");
    mu.sec_crossing_hits.       SetBranches(mu.tree, "SecCross");
    mu.top_reg.                 SetBranches(mu.tree, "Top");
    mu.bot_reg.                 SetBranches(mu.tree, "Bot");
    mu.tree->Branch("TopdQds",              &mu.top_dQds);
    mu.tree->Branch("BotdQds",              &mu.bot_dQds);

    // mu.trk.hits.                SetBranches(mu.tree, "Trk");
    // mu.trk.end_hit.             SetBranches(mu.tree, "TrkEnd");
    // mu.tree->Branch("TrkHitIndex", &mu.trk.hit_indices);
    // mu.tree->Branch("TrkHitDx",    &mu.trk.hit_dxs);

    mu.tree->Branch("TruPdg",               &mu.tru.pdg);
    mu.tree->Branch("TruEndProcess",        &mu.tru.end_process);
    mu.tree->Branch("TruEndEnergy",         &mu.tru.end_energy);
    mu.tree->Branch("TruDownward",          &mu.tru.downward);
    mu.tree->Branch("TruSortedHitsError",   &mu.tru.sh_error);
    mu.tree->Branch("TruMaxConsecutiveDist",&mu.tru.max_consecutive_dist);
    mu.tree->Branch("TruCathodeCrossing",   &mu.tru.cathode_crossing);
    mu.tree->Branch("TruAnodeCrossing",     &mu.tru.anode_crossing);
    mu.tree->Branch("TruSectionJumping",    &mu.tru.section_jumping);
    mu.tree->Branch("TruSectionMisaligned", &mu.tru.section_misaligned);
    mu.tree->Branch("TruCathodeMisaligned", &mu.tru.cathode_misaligned);
    mu.tru.start_point.         SetBranches(mu.tree, "TruStart");
    mu.tru.end_point.           SetBranches(mu.tree, "TruEnd");
    mu.tru.start_hit.           SetBranches(mu.tree, "TruStart");
    mu.tru.end_hit.             SetBranches(mu.tree, "TruEnd");
    mu.tru.top_last_hit.        SetBranches(mu.tree, "TruTopLast");
    mu.tru.bot_first_hit.       SetBranches(mu.tree, "TruBotFirst");
    mu.tru.hits.                SetBranches(mu.tree, "Tru");
    mu.tru.sec_crossing_hits.   SetBranches(mu.tree, "TruSecCross");
    mu.tru.top_reg.             SetBranches(mu.tree, "TruTop");
    mu.tru.bot_reg.             SetBranches(mu.tree, "TruBot");
    mu.tree->Branch("TruTopdQds",           &mu.tru.top_dQds);
    mu.tree->Branch("TruBotdQds",           &mu.tru.bot_dQds);
    mu.tree->Branch("TruHasMichel",         &mu.tru.has_michel);
    mu.tree->Branch("TruMichelEnergy",      &mu.tru.michel_energy);
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

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    reset_ev();

    ev.run = e.run();
    ev.subrun = e.subRun();
    ev.event = e.event();
    ev.is_data = e.isRealData();

    for (PtrHit ph : vph_ev)
        if (ph->View() == geo::kW)
            ev.hits.push_back(GetHit(ph));

    for (PtrTrk pt_ev : vpt_ev) {
        VecPtrHit vph_trk = fmp_trk2hit.at(pt_ev.key());
        std::vector<recob::TrackHitMeta const*> const& vhm_trk = fmp_trk2hit.data(pt_ev.key());
        ASSERT(vph_trk.size())

        reset_mu();

        mu.length = pt_ev->Length();
        bool up = IsUpright(*pt_ev);
        mu.start_point = ana::Point(up ? pt_ev->Start() : pt_ev->End());
        mu.end_point = ana::Point(up ? pt_ev->End() : pt_ev->Start());

        // Trying TrackHitMeta
        ASSERT(vph_trk.size() == vhm_trk.size())
        // float min_dist = std::numeric_limits<float>::max();
        std::vector<unsigned> bad_hit_indices;
        std::map<size_t, unsigned> map_hitkey_to_metaidx;
        for (unsigned i=0; i<vph_trk.size(); i++) {
            // PtrHit const& ph = vph_trk[i];
            // recob::TrackHitMeta const* ihm = vhm_trk[i];
            // if (ph->View() != geo::kW) continue;
            // mu.trk.hits.push_back(GetHit(ph));
            // unsigned idx = ihm->Index();
            // mu.trk.hit_indices.push_back(idx);
            // mu.trk.hit_dxs.push_back(ihm->Dx());

            // if (pt_ev->HasValidPoint(idx)) {
            //     geo::Point_t pt = pt_ev->LocationAtPoint(idx);
            //     float dist = sqrt(
            //         pow(pt.x() - mu.end_point.x, 2) +
            //         pow(pt.y() - mu.end_point.y, 2) +
            //         pow(pt.z() - mu.end_point.z, 2)
            //     );
            //     if (dist < min_dist) {
            //         min_dist = dist;
            //         mu.trk.end_hit = GetHit(ph);
            //     }
            // } else bad_hit_indices.push_back(i);

            if (vph_trk[i]->View() != geo::kW) continue;
            if (!pt_ev->HasValidPoint(vhm_trk[i]->Index())) {
                bad_hit_indices.push_back(i);
            } else {
                map_hitkey_to_metaidx[vph_trk[i].key()] = i;
            }
        }
        if (fRemoveBadHits)
            for (int i=bad_hit_indices.size()-1; i>=0; i--)
                vph_trk.erase(vph_trk.begin() + bad_hit_indices[i]);

        // ana::SortedHits sh = GetSortedHits(vph_trk, mu.end_point.z > mu.end_point.z ? 1 : -1);
        ana::SortedHits sh = GetSortedHits(vph_trk);

        mu.sh_error = !sh;
        if (fAssert) { ASSERT(!mu.sh_error) }
        else         { LOG(!mu.sh_error); }
        if (!mu.sh_error) {

            mu.top_reg = sh.regs[0];
            mu.bot_reg = sh.regs[1];
            mu.start_hit = GetHit(sh.start());
            mu.end_hit = GetHit(sh.end());
            size_t end_track_idx = vhm_trk[map_hitkey_to_metaidx.at(sh.end().key())]->Index();
            mu.end_hit_y = pt_ev->HasValidPoint(end_track_idx)
                ? pt_ev->LocationAtPoint(end_track_idx).Y()
                : -500.F;

            if (sh.is_cc()) {
                mu.top_last_hit = GetHit(sh.cc_first());
                mu.bot_first_hit = GetHit(sh.cc_second());
            } else {
                mu.top_last_hit = ana::Hit{};
                mu.bot_first_hit = ana::Hit{};
            }

            for (PtrHit ph : sh.vph)
                mu.hits.push_back(GetHit(ph));

            for (PtrHit ph : sh.scs())
                mu.sec_crossing_hits.push_back(GetHit(ph));

            mu.top_dQds = GetdQds(sh.vph.begin(), sh.after_cathode_it(), 6);
            mu.bot_dQds = GetdQds(sh.after_cathode_it(), sh.vph.end(), 6);

            mu.max_consecutive_dist = 0.F;
            for (VecPtrHit::iterator it=sh.vph.begin(); it!=sh.vph.end()-1; ++it) {
                Hit h1 = GetHit(*it), h2 = GetHit(*(it+1));
                if (h1.section != h2.section) continue;
                float dist = GetDistance(h1, h2);
                if (dist > mu.max_consecutive_dist)
                    mu.max_consecutive_dist = dist; 
            }

            mu.cathode_crossing = sh.is_cc();
            mu.cathode_misaligned = sh.is_cc() 
                && abs(sh.cc_first()->PeakTime()-sh.cc_second()->PeakTime()) > (geoCathodeGap + misalignmentTolerance)/fTick2cm;
            if (fAssert) { ASSERT(mu.cathode_crossing) }
            else         { LOG(mu.cathode_crossing); }
            LOG(mu.cathode_misaligned);

            if (geoDet == kPDVD)
                mu.anode_crossing = mu.start_hit.section < 4
                    && geoHighX.z.isInside(mu.start_hit.space, 10.F)
                    && geoTickWindow.isInside(mu.start_hit.tick, 10.F/fTick2cm);
            else if (geoDet == kPDHD)
                mu.anode_crossing =
                    geoHighX.z.isInside(mu.start_hit.space, 10.F)
                    && geoTickWindow.isInside(mu.start_hit.tick, 10.F/fTick2cm);
            LOG(mu.anode_crossing);

            LOG(sh.secs.size() > 1);
            if (sh.secs.size() > 1) {
                VecPtrHit scs = sh.scs();
                VecPtrHit::iterator sc_it = scs.begin();
                for (unsigned i=0; i<sh.secs.size()-1; i++) {
                    int sec_curr = sh.secs[i];
                    int sec_next = sh.secs[i+1];

                    if (ana::sec2side.at(geoDet).at(sec_curr) == ana::sec2side.at(geoDet).at(sec_next)) {
                        if (abs(sec_curr - sec_next) != 1)
                            mu.section_jumping = true;

                        if (((*sc_it)->PeakTime()-(*++sc_it)->PeakTime()) > misalignmentTolerance/fTick2cm)
                            mu.section_misaligned = true;
                    }
                }
                LOG(mu.section_jumping);
                LOG(mu.section_misaligned);
            }    
        }

        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        if (mcp) {
            mu.tru.pdg = mcp->PdgCode();
            mu.tru.end_process = mcp->EndProcess();
            mu.tru.end_energy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV
            mu.tru.downward = geoDet == kPDVD ? mcp->Vx(0) > mcp->EndX() : mcp->Vy(0) > mcp->EndY();

            mu.tru.start_point = ana::Point(mcp->Position().Vect());
            mu.tru.end_point = ana::Point(mcp->EndPosition().Vect());

            VecPtrHit vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            ana::SortedHits sh_mcp = GetSortedHits(vph_mcp_mu);

            mu.tru.sh_error = !sh_mcp;
            LOG(!mu.tru.sh_error);
            if (!mu.tru.sh_error) {
                bool is_up = mcp->Vx() > mcp->EndX();
                if (!is_up) sh_mcp.reverse();

                mu.tru.top_reg = sh_mcp.regs[0];
                mu.tru.bot_reg = sh_mcp.regs[1];
                mu.tru.start_hit = GetHit(sh_mcp.start());
                mu.tru.end_hit = GetHit(sh_mcp.end());
                if (sh_mcp.is_cc()) {
                    mu.tru.top_last_hit = GetHit(sh_mcp.cc_first());
                    mu.tru.bot_first_hit = GetHit(sh_mcp.cc_second());
                } else {
                    mu.tru.top_last_hit = ana::Hit{};
                    mu.tru.bot_first_hit = ana::Hit{};
                }

                for (PtrHit ph : sh_mcp.vph)
                    mu.tru.hits.push_back(GetHit(ph));

                for (PtrHit ph : sh_mcp.scs())
                    mu.tru.sec_crossing_hits.push_back(GetHit(ph));

                mu.tru.top_dQds = GetdQds(sh_mcp.vph.begin(), sh_mcp.after_cathode_it(), 6);
                mu.tru.bot_dQds = GetdQds(sh_mcp.after_cathode_it(), sh_mcp.vph.end(), 6);

                mu.tru.max_consecutive_dist = 0.F;
                for (VecPtrHit::iterator it=sh_mcp.vph.begin(); it!=sh_mcp.vph.end()-1; ++it) {
                    Hit h1 = GetHit(*it), h2 = GetHit(*(it+1));
                    if (h1.section != h2.section) continue;
                    float dist = GetDistance(h1, h2);
                    if (dist > mu.tru.max_consecutive_dist)
                        mu.tru.max_consecutive_dist = dist; 
                }

                mu.tru.cathode_crossing = sh_mcp.is_cc();
                mu.tru.cathode_misaligned = sh_mcp.is_cc() 
                    && !(abs(sh_mcp.cc_first()->PeakTime()-sh_mcp.cc_second()->PeakTime())*fTick2cm < 3 * geoCathodeGap);
                LOG(mu.tru.cathode_crossing);
                LOG(mu.tru.cathode_misaligned);

                // if (geoDet == kPDVD)
                //     mu.tru.anode_crossing = mu.tru.start_hit.section < 4
                //         && geoHighX.z.isInside(mu.tru.start_hit.space, 10.F)
                //         && geoTickWindow.isInside(mu.tru.start_hit.tick, 10.F/fTick2cm);
                // else if (geoDet == kPDHD)
                //     mu.tru.anode_crossing =
                //         geoHighX.z.isInside(mu.tru.start_hit.space, 10.F)
                //         && geoTickWindow.isInside(mu.tru.start_hit.tick, 10.F/fTick2cm);
                float t_anode = (geoHighX.x.max - mcp->Vx(0)) / mcp->Px(0);
                float z_anode = mcp->Vz(0) + mcp->Pz(0) * t_anode;
                float y_anode = mcp->Vy(0) + mcp->Py(0) * t_anode;
                mu.tru.anode_crossing = geoHighX.z.isInside(z_anode, 10.F) && geoHighX.y.isInside(y_anode, 10.F);
                // the anode crossing might be out of the tick window!
                LOG(mu.tru.anode_crossing);



                LOG(sh_mcp.secs.size() > 1);
                if (sh_mcp.secs.size() > 1) {
                    VecPtrHit scs = sh.scs();
                    VecPtrHit::iterator sc_it = scs.begin();
                    for (unsigned i=0; i<sh_mcp.secs.size()-1; i++) {
                        int sec_curr = sh_mcp.secs[i];
                        int sec_next = sh_mcp.secs[i+1];

                        if (ana::sec2side.at(geoDet).at(sec_curr) == ana::sec2side.at(geoDet).at(sec_next)) {
                            if (abs(sec_curr - sec_next) != 1)
                                mu.tru.section_jumping = true;

                            if (((*sc_it)->PeakTime()-(*++sc_it)->PeakTime())*fTick2cm > 2.F)
                                mu.tru.section_misaligned = true;
                        }
                    }
                    LOG(mu.tru.section_jumping);
                    LOG(mu.tru.section_misaligned);
                }    
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