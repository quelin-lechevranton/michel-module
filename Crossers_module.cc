#include "utils.h"

/* 
Naming conventions,
some variables are prefixed by their type followed by an underscore:
 - `ph_`  for `art::Ptr<recob::Hit>`
 - `vph_` for `std::vector<art::Ptr<recob::Hit>>`

 - `pt_`  for `art::Ptr<recob::Track>`
 - `vpt_` for `std::vector<art::Ptr<recob::Track>>`

 - `ps_`  for `art::Ptr<recob::Shower>`
 - `vps_` for `std::vector<art::Ptr<recob::Shower>>`

 - `tag_` for `art::InputTag`
 - `mcp_` for `simb::MCParticle*`

 - `sh_` for custom `ana::SortedHits`, a structure containing info on the sorted list of hits of a track
MichelAnalysis class properties are prefixed:
 - `as`  for art services
 - `geo` for geometric information
 - `in`  for input variable, taken from fhicl file

 - `ev`  for reco info about the full event (stored in `evTree`)
 - `mu`  for reco info about the track and michel electron (stored in `muTree`)
 - `tru` for truth info about the track (stored in `muTree`)
 - `mi`  for truth info about the michel electron (stored in `muTree`)

Top means X>0, Bot means X<0 (by reference to PDVD)
*/

namespace ana { class Crossers; }

class ana::Crossers: 
    public art::EDAnalyzer, 
    private ana::MichelModule 
{
public:
    explicit Crossers(fhicl::ParameterSet const& p);
    Crossers(Crossers const&) = delete;
    Crossers(Crossers&&) = delete;
    Crossers& operator=(Crossers const&) = delete;
    Crossers& operator=(Crossers&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    // Geometry Information
    ana::Bounds<float> wireWindow;
    ana::Bounds3D<float> geoTop, geoBot;
    float geoCathodeGap; // cm

    // Input Parameters
    bool        inLog;
    float       inTrackLengthCut; // in cm
    float       inFiducialLength; // in cm
    unsigned    inRegN;

    // Output Variables
    TTree       *evTree;
    TTree       *trTree;

    // Event information
    unsigned                evRun;
    unsigned                evSubRun;
    unsigned                evEvent;
    unsigned                evIndex=0;
    unsigned                evTrackNumber;
    std::vector<unsigned>   evTrackIndices;
    bool                    evIsData;
    ana::Hits               evHits;

    // Track information: recob::Track
    unsigned                trIndex=0;
    float                   trLength;
    ana::Point              trStartPoint;
    ana::Point              trEndPoint;
    ana::Points             trPoints;

    // Track information: recob::Hit
    ana::Hits               trHits;
    std::vector<float>      trHitAX;
    std::vector<float>      trHitCX;
    std::vector<float>      trHitY;
    // float                   trEndAngle;
    bool                    trCathodeCrossing;
    float                   trCathodeAlignment;
    bool                    trAnodeCrossing;

    ana::LinearRegression   trStartReg;
    ana::LinearRegression   trGhostReg;

    bool                    trGhostTrack;

    std::vector<float>      trHitdQds;
    
    // Truth information: Muon
    int                     truPdg;
    float                   truEnergy;
    std::string             truEndProcess;
    ana::Point              truStartPoint;
    ana::Point              truEndPoint;
    ana::Points             truPoints;
    float                   truEndEnergy;

    ana::Point              truCathodePoint;
    bool                    truCathodeCrossing;
    ana::Point              truAnodePoint;
    bool                    truAnodeCrossing;

    // ana::Hit                truStartHit;
    // ana::Hit                truEndHit;
    // ana::LinearRegression   truReg;
    // float                   truEndAngle;
    // enum EnumHasMichel: int { 
    //     kHasNoMichel        = 0, 
    //     kHasMichelOutside   = 1, 
    //     kHasMichelInside    = 2, 
    //     kHasMichelFiducial  = 3
    // };
    // EnumHasMichel           truHasMichel;

    // Truth information: Michel electron
    // float               miTrueEnergy;
    // float               miTrackLength;
    // float               miShowerLength;
    // ana::Hits           miHits;
    // std::vector<float>  miHitEnergyFrac;
    // std::vector<float>  miHitMuonAngle;
    // float               miHitEnergy;

    // Truth information: Hits nearby muon's end
    // unsigned            miBaryNHit;
    // ana::Vec2           miBary;
    // float               miBaryAngle;
    // float               miBaryMuonAngle;

    void resetEvent(void);
    void resetMuon(void);
};

ana::Crossers::Crossers(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} 
    , MichelModule{p}
    , inLog(p.get<bool>("Log", true))
    , inTrackLengthCut(p.get<float>("TrackLengthCut", 30.F)) // in cm
    , inFiducialLength(p.get<float>("FiducialLength", 20.F)) // in cm
    , inRegN(p.get<unsigned>("RegN", 6))
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
    wireWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoBot = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(),
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoTop = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoBot = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoTop = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        case kPDSP:
            geoBot = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 9}).Max()
            };
            geoTop = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 10}).Max()
            };
            break;
        default: break;
    }
    geoCathodeGap = geoTop.x.min - geoBot.x.max;

    std::cout << "Crossers: " "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << ana::det_name[geoDet] << ")" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  Top Bounds: " << geoTop << std::endl
        << "  Bot Bounds: " << geoBot << std::endl;
    std::cout << "Crossers: " "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << inTrackLengthCut << " cm" << std::endl
        << "  Fiducial Length: " << inFiducialLength << " cm" << std::endl
        << "  Smoothing Length: " << inRegN << " points" << std::endl;

    evTree = asFile->make<TTree>("event","");

    evTree->Branch("EventRun",      &evRun);
    evTree->Branch("EventSubRun",   &evSubRun);
    evTree->Branch("EventEvent",    &evEvent);
    evTree->Branch("Index",         &evIndex);
    evTree->Branch("MuonNumber",    &evTrackNumber);
    evTree->Branch("MuonIndices",   &evTrackIndices);
    evTree->Branch("IsData",        &evIsData);
    SetBranches(evTree, "",         &evHits);

    trTree = asFile->make<TTree>("track","");

    // Event
    trTree->Branch("EventRun",      &evRun);
    trTree->Branch("EventSubRun",   &evSubRun);
    trTree->Branch("EventEvent",    &evEvent);
    trTree->Branch("IsData",        &evIsData);
    trTree->Branch("EventIndex",    &evIndex);
    trTree->Branch("IndexInEvent",  &evTrackNumber);
    trTree->Branch("Index",         &trIndex);

    // Track
    trTree->Branch("Length",        &trLength);
    SetBranches(trTree, "Start",    &trStartPoint);
    SetBranches(trTree, "End",      &trEndPoint);
    SetBranches(trTree, "",         &trPoints);

    // Hit
    trTree->Branch("CathodeCrossing",       &trCathodeCrossing);
    trTree->Branch("CathodeAlignment",      &trCathodeAlignment);
    trTree->Branch("AnodeCrossing",         &trAnodeCrossing);
    SetBranches(trTree, "Start",            &trStartReg);
    SetBranches(trTree, "Ghost",            &trGhostReg);
    trTree->Branch("GhostTrack",            &trGhostTrack);
    // trTree->Branch("EndAngle",              &trEndAngle);
    SetBranches(trTree, "",                 &trHits);
    trTree->Branch("HitAX",                 &trHitAX);
    trTree->Branch("HitCX",                 &trHitCX);
    trTree->Branch("HitY",                  &trHitY);
    trTree->Branch("HitdQds",               &trHitdQds);

    // Truth
    trTree->Branch("TruePdg",               &truPdg);
    trTree->Branch("TrueEnergy",            &truEnergy);
    trTree->Branch("TrueEndProcess",        &truEndProcess);
    SetBranches(trTree, "TrueStart",        &truStartPoint);
    SetBranches(trTree, "TrueEnd",          &truEndPoint);
    SetBranches(trTree, "True",             &truPoints);
    trTree->Branch("TrueEndEnergy",         &truEndEnergy);

    SetBranches(trTree, "TrueCathode",      &truCathodePoint);
    trTree->Branch("TrueCathodeCrossing",   &truCathodeCrossing);
    SetBranches(trTree, "TrueAnode",        &truAnodePoint);
    trTree->Branch("TrueAnodeCrossing",     &truAnodeCrossing);

    // SetBranches(trTree, "TrueStart",        &truStartHit);
    // SetBranches(trTree, "TrueEnd",          &truEndHit);
    // SetBranches(trTree, "True",             &truReg);
    // trTree->Branch("TrueEndAngle",          &truEndAngle);
    // trTree->Branch("TrueHasMichel",   (int*)&truHasMichel);
    // trTree->Branch("MichelTrueEnergy",      &miTrueEnergy); // MeV
    // trTree->Branch("MichelTrackLength",     &miTrackLength); // cm
    // trTree->Branch("MichelShowerLength",    &miShowerLength); // cm
    // SetBranches(trTree, "Michel",           &miHits);
    // trTree->Branch("MichelHitEnergyFrac",   &miHitEnergyFrac);
    // trTree->Branch("MichelHitMuonAngle",    &miHitMuonAngle);
    // trTree->Branch("MichelHitEnergy",       &miHitEnergy); // ADC
}

void ana::Crossers::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) {
        std::cout << "Crossers: " "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) {
        std::cout << "Crossers: " "\033[1;91m" "No valid recob::Track handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrTrk vpt_ev;
    art::fill_ptr_vector(vpt_ev, vh_trk);

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    resetEvent();

    // dump event information
    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    evIsData = e.isRealData();
    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            evHits.push_back(GetHit(p_hit));

    // loop over tracks to find stopping muons
    for (PtrTrk const& pt_ev : vpt_ev) {
        if (inLog) std::cout << "e" << evIndex << "t" << pt_ev->ID() << "\r" << std::flush;
        resetMuon();

        // get hits and metadata associated to the track
        VecPtrHit vph_mu_all = fmp_trk2hit.at(pt_ev.key());
        std::vector<recob::TrackHitMeta const*> const& vhm_mu = fmp_trk2hit.data(pt_ev.key());
        // std::map<size_t, unsigned> map_hitkey2metaidx;
        std::map<size_t, size_t> map_hitkey2trkidx;
        ASSERT(!vph_mu_all.empty())
        ASSERT(vph_mu_all.size() == vhm_mu.size())

        // remove hits not associated to a track point (hit misassociated to the track?)
        VecPtrHit vph_mu;
        for (unsigned i=0; i<vph_mu_all.size(); i++) {
            if (vph_mu_all[i]->View() != geo::kW) continue;
            if (!pt_ev->HasValidPoint(vhm_mu[i]->Index())) continue;
            map_hitkey2trkidx[vph_mu_all[i].key()] = vhm_mu[i]->Index();
            vph_mu.push_back(vph_mu_all[i]);
        }
        vph_mu_all.clear();
        ASSERT(!vph_mu.empty())

        bool track_is_up =  IsUpright(*pt_ev);
        geo::Point_t Start = track_is_up ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = track_is_up ? pt_ev->End() : pt_ev->Start();
        trStartPoint = ana::Point(Start);
        trEndPoint = ana::Point(End);


        if (inLog) std::cout << "\t" "\033[1;93m" "e" << evIndex << "m" << evTrackNumber << " (" << trIndex << ")" "\033[0m" << std::endl;

        // dump basic track information
        trLength = pt_ev->Length();
        ASSERT(trLength > inTrackLengthCut)

        std::sort(vph_mu.begin(), vph_mu.end(), [&map_hitkey2trkidx](PtrHit const& ph1, PtrHit const& ph2) {
            return map_hitkey2trkidx.at(ph1.key()) < map_hitkey2trkidx.at(ph2.key());
        });

        bool is_in_bot=false, is_in_top=false;
        for (PtrHit const& ph : vph_mu) {
            if (is_in_bot && is_in_top) break;
            if (!is_in_bot && GetSide(ph) == kBot) is_in_bot = true;
            if (!is_in_top && GetSide(ph) == kTop) is_in_top = true;
        }
        bool is_cc = is_in_bot && is_in_top;

        bool is_upright = false;
        if (geoDet == kPDVD) {
            Side_t front_side = GetSide(vph_mu.front());
            if (is_cc)
                is_upright = front_side == kTop;
            else if (front_side == kTop)
                is_upright = vph_mu.front()->PeakTime() < vph_mu.back()->PeakTime();
            else
                is_upright = vph_mu.front()->PeakTime() > vph_mu.back()->PeakTime();
        } else if (geoDet == kPDHD || geoDet == kPDSP) {
            float front_y = pt_ev->LocationAtPoint(map_hitkey2trkidx.at(vph_mu.front().key())).Y();
            float back_y = pt_ev->LocationAtPoint(map_hitkey2trkidx.at(vph_mu.back().key())).Y();
            is_upright = front_y > back_y;
        }
        if (!is_upright)
            std::reverse(vph_mu.begin(), vph_mu.end());

        Side_t first_side = GetSide(vph_mu.front());
        VecPtrHit::iterator cathode_it = std::find_if(
            vph_mu.begin(), vph_mu.end(), [&](PtrHit const& ph) {
                return GetSide(ph) != first_side;
            }
        );
        // if (is_cc) { ASSERT(cathode_it != vph_mu.end()) }
        PtrHit const cc_first = is_cc ? *(cathode_it-1) : PtrHit{};
        PtrHit const cc_second = is_cc ? *cathode_it : PtrHit{};

        // dump dQ/ds
        trHitdQds.clear();
        std::vector<float> tmp;
        tmp = GetdQds(vph_mu.begin(), cathode_it, inRegN);
        trHitdQds.insert(trHitdQds.end(), tmp.begin(), tmp.end());
        tmp = GetdQds(cathode_it, vph_mu.end(), inRegN);
        trHitdQds.insert(trHitdQds.end(), tmp.begin(), tmp.end());
        tmp.clear();

        // dump cathode crossing boolean: defined as track with 4+ hits on both sides of the cathode
        trCathodeCrossing = is_cc;
        trCathodeAlignment = is_cc ? abs(cc_first->PeakTime()-cc_second->PeakTime())*fTick2cm : -1.F;
        LOG(trCathodeCrossing);

        float start_t = vph_mu.front()->PeakTime();
        float start_y = pt_ev->HasValidPoint(map_hitkey2trkidx.at(vph_mu.front().key()))
            ? pt_ev->LocationAtPoint(map_hitkey2trkidx.at(vph_mu.front().key())).Y()
            : util::kBogusF;
        float start_z = GetSpace(vph_mu.front());
        switch (geoDet) {
        case kPDVD: /* ASSUMS DOWNWARD MUON */
            trAnodeCrossing = GetSide(vph_mu.front()) == kTop
                && geoTop.y.isInside(start_y, inFiducialLength)
                && geoTop.z.isInside(start_z, inFiducialLength)
                && wireWindow.isInside(start_t, inFiducialLength/fTick2cm);
            break;
        case kPDHD:
            trAnodeCrossing =
                geoTop.y.isInside(start_y, inFiducialLength)
                && geoTop.z.isInside(start_z, inFiducialLength)
                && wireWindow.isInside(start_t, inFiducialLength/fTick2cm);
            break;
        case kPDSP:
            trAnodeCrossing = false;
            break;
        }
        LOG(trAnodeCrossing);


        // dump hits
        for (PtrHit const& ph_mu : vph_mu) {
            ana::Hit hit = GetHit(ph_mu);
            trHits.push_back(hit);

            size_t hit_track_idx = map_hitkey2trkidx.at(ph_mu.key());
            trHitY.push_back(pt_ev->HasValidPoint(hit_track_idx)
                ? pt_ev->LocationAtPoint(hit_track_idx).Y()
                : util::kBogusF
            );
        }

        for (size_t index=pt_ev->NextValidPoint(0); index!=recob::TrackTrajectory::InvalidIndex; index=pt_ev->NextValidPoint(index)) {
            trPoints.push_back(pt_ev->LocationAtPoint(index));
        }

        if (trCathodeCrossing) {
            Side_t start_side = GetSide(vph_mu.front());
            PtrHit const& cc_bot = start_side == kBot ? cc_first : cc_second;
            PtrHit const& cc_top = start_side == kTop ? cc_first : cc_second;

            for (PtrHit const& ph_mu : vph_mu) {
                trHitCX.push_back(GetSide(ph_mu) == kBot
                    ? -(geoCathodeGap/2) - (cc_bot->PeakTime() - ph_mu->PeakTime()) * fTick2cm
                    : +(geoCathodeGap/2) + (cc_top->PeakTime() - ph_mu->PeakTime()) * fTick2cm
                );
            }
        }
        if (trAnodeCrossing) {
            Side_t start_side = GetSide(vph_mu.front());

            for (PtrHit const& ph_mu : vph_mu) {
                if (GetSide(ph_mu) == start_side)
                    trHitAX.push_back(GetSide(ph_mu) == kTop
                        ? geoTop.x.max - (ph_mu->PeakTime() - start_t) * fTick2cm
                        : geoBot.x.min + (ph_mu->PeakTime() - start_t) * fTick2cm
                    );
                else
                    trHitAX.push_back(util::kBogusF);
            }
        }

        // Search for ghost track
        int n = 0;
        for (PtrHit const& ph_mu : vph_mu) {
            if (n++ == 20) break;
            trStartReg.add( GetSpace(ph_mu), ph_mu->PeakTime()*fTick2cm );
        }
        trStartReg.compute();

        int induction_hits = 0;
        for (PtrHit const& ph_ev : vph_ev) {
            PtrTrk const& pt = fop_hit2trk.at(ph_ev.key());
            if (pt.isNonnull() && pt->Length() > inTrackLengthCut) continue;

            if (ph_ev->View() != geo::kW) {
                if (GetDistance(ph_ev, (Side_t)-1, start_y, start_z, start_t, true) < 20)
                    induction_hits++;
                continue;
            }
            if (GetDistance(ph_ev, vph_mu.front()) > 20) continue;
            trGhostReg.add( GetSpace(ph_ev), ph_ev->PeakTime()*fTick2cm );
        }
        trGhostReg.compute();

        trGhostTrack = trGhostReg.r2 > 0.8 
            && abs( (trGhostReg.m - trStartReg.m) / trStartReg.m ) < 10 
            && (trGhostReg.n - induction_hits) > 0;


        // Truth Information
        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        // simb::MCParticle const* mcp_mi = nullptr;
        // VecPtrHit vph_mcp_mu, vph_mi;
        // std::vector<float> energyFracs_mi;
        // if (mcp) {
        //     mcp_mi = GetMichelMCP(mcp);
        //     vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
        //     vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true, &energyFracs_mi);
        // }

        LOG(mcp);
        if (mcp) {
            truPdg = mcp->PdgCode();
            truEnergy = mcp->E();
            truEndProcess = mcp->EndProcess();
            truStartPoint = ana::Point(mcp->Position().Vect());
            truEndPoint = ana::Point(mcp->EndPosition().Vect());
            truEndEnergy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV


            // truCathodeCrossing 
            int before_anode=-1, before_cathode=-1;
            TVector3 prev_pt = mcp->Position().Vect();
            for (size_t i=1; i<mcp->NumberTrajectoryPoints(); i++) {
                TVector3 const& pt = mcp->Position(i).Vect();
                truPoints.push_back(ana::Point(pt));

                if (before_cathode != -1 && before_anode != -1) break;

                if (before_cathode == -1 && prev_pt.X() * pt.X() < 0) {
                    before_cathode = i-1;
                }

                if (before_anode != -1) continue;

                switch (geoDet) {
                case kPDVD:
                    if (prev_pt.X() > geoTop.x.max && geoTop.x.max > pt.X())
                        before_anode = i-1;
                    break;
                case kPDHD: 
                    std::cout << "\033[1;91m" "truCathodeCrossing not implemented for PDHD" "\033[0m" << std::endl;
                    break;
                case kPDSP: 
                    std::cout << "\033[1;91m" "truCathodeCrossing not implemented for PDSP" "\033[0m" << std::endl;
                    break;
                }
            }
            if (before_cathode != -1) {
                truCathodePoint = ana::Point(mcp->Position(before_cathode).Vect());
                truCathodeCrossing = geoTop.y.isInside(truCathodePoint.y, 5.F)
                    && geoTop.z.isInside(truCathodePoint.z, 5.F);
            }
            if (before_anode != -1) {
                truAnodePoint = ana::Point(mcp->Position(before_anode).Vect());
                truAnodeCrossing = geoTop.y.isInside(truAnodePoint.y, 5.F)
                    && geoTop.z.isInside(truAnodePoint.z, 5.F);
            }

            // LOG(mcp_mi);
            // if (mcp_mi) {
            //     truHasMichel = (
            //         geoTop.isInside(mcp_mi->Position().Vect(), 20.F)
            //         || geoBot.isInside(mcp_mi->Position().Vect(), 20.F)
            //     ) ? kHasMichelFiducial : (
            //         geoTop.isInside(mcp_mi->Position().Vect())
            //         || geoBot.isInside(mcp_mi->EndPosition().Vect())
            //         ? kHasMichelInside
            //         : kHasMichelOutside
            //     );
            //     miTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3;
            //     PtrTrk pt_mi = ana::mcp2trk(mcp_mi, vpt_ev, clockData, fmp_trk2hit);
            //     miTrackLength = pt_mi ? pt_mi->Length() : util::kBogusF;
            //     // PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
            //     // MichelShowerLength = ps_mi ? ps_mi->Length() : util::kBogusF;

            //     for (size_t i=0; i<vph_mi.size(); i++) {
            //         PtrHit const& ph_mi = vph_mi[i];
            //         float energyFrac = energyFracs_mi[i];
            //         if (ph_mi->View() != geo::kW) continue;
            //         miHits.push_back(GetHit(ph_mi));
            //         miHitEnergyFrac.push_back(energyFrac);
            //     }
            //     miHitEnergy = std::accumulate(miHits.adc.begin(), miHits.adc.end(), 0.F);
            // }
        }

        trTree->Fill();
        evTrackIndices.push_back(trIndex);
        trIndex++;
        evTrackNumber++;
    } // end of loop over tracks
    evTree->Fill();
    evIndex++;
}

void ana::Crossers::beginJob() {}
void ana::Crossers::endJob() {}

void ana::Crossers::resetEvent() {
    evTrackNumber = 0;
    evTrackIndices.clear();
    evHits.clear();
}
void ana::Crossers::resetMuon() {
    trHits.clear();
    trHitAX.clear();
    trHitCX.clear();
    trHitY.clear();
    // trEndAngle = util::kBogusF;
    trCathodeCrossing = false;
    trCathodeAlignment = util::kBogusF;
    trAnodeCrossing = false;
    trStartReg.clear();
    trGhostReg.clear();
    trGhostTrack = false;
    trHitdQds.clear();

    truPdg = 0;
    truEnergy = util::kBogusF;
    truEndProcess = "";
    truStartPoint = ana::Point{};
    truEndPoint = ana::Point{};
    truEndEnergy = util::kBogusF;

    truCathodePoint = ana::Point{};
    truCathodeCrossing = false;
    truAnodePoint = ana::Point{};
    truAnodeCrossing = false;

    // truStartHit = ana::Hit{};
    // truEndHit = ana::Hit{};
    // truReg = ana::LinearRegression{};
    // truEndAngle = util::kBogusF;
    // truHasMichel = kHasNoMichel;

    // miTrueEnergy = util::kBogusF;
    // miTrackLength = util::kBogusF;
    // miHits.clear();
    // miHitEnergyFrac.clear();
    // miHitEnergy = util::kBogusF;
}

DEFINE_ART_MODULE(ana::Crossers)