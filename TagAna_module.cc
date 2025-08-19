////////////////////////////////////////////////////////////////////////
// Class:       TagAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        Agnochecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"

namespace ana {
    class TagAna;
}

class ana::TagAna : public art::EDAnalyzer, private ana::MichelAnalyzer {
public:
    explicit TagAna(fhicl::ParameterSet const& p);
    TagAna(TagAna const&) = delete;
    TagAna(TagAna&&) = delete;
    TagAna& operator=(TagAna const&) = delete;
    TagAna& operator=(TagAna&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    ana::Bounds<float> wireWindow;
    ana::Bounds3D<float> geoHighX, geoLowX;
    float fCathodeGap; // cm

    // Input Parameters
    bool fLog;
    float fTrackLengthCut; // in cm
    float fMichelRadius; // in cm
    float fNearbyRadius; // in cm
    float fBodyDistance; // in cm
    unsigned fRegN;
    float fBraggThreshold; // in MIP

    // Output Variables
    unsigned evRun, evSubRun, evEvent;

    TTree* tEvent;
    bool EventIsReal;
    unsigned iEvent=0;
    unsigned EventNMuon;
    std::vector<unsigned> EventiMuon;

    ana::Hits EventHits;


    TTree* tMuon;
    unsigned iMuon=0;

    bool TagIsUpright; // Supposition: Muon is downward
    bool TagCathodeCrossing;
    bool TagAnodeCrossing;
    float CutTrackLength;
    bool TagEndInWindow;
    bool TagEndInVolume;
    bool TagHitCathodeCrossing;
    bool TagHitGoodCathodeCrossing;

    float CutdQdxMax;
    enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
    int TagBraggError;
    ana::Hits BraggMuonHits;

    bool AgnoTagTrkEndLowN,
         AgnoTagTrkEndBadCC,
         AgnoTagEndHitError;
    unsigned AgnoCountAll=0,
             AgnoCountNoHit=0,
             AgnoCountNoSec=0;

    bool TrueTagDownward;
    int TrueTagPdg;
    std::string TrueTagEndProcess;

    int TrueTagHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside, kHasMichelFiducial };

    ana::Hits MuonHits;
    ana::Hit MuonEndHit;
    ana::Point MuonEndPoint;
    ana::Hit MuonTrueEndHit;
    ana::Point MuonTrueEndPoint;
    ana::Hit BraggEndHit;

    float MichelTrackLength;
    ana::Hits MichelHits;
    float MichelTrueEnergy, MichelHitEnergy;
    float BraggSphereEnergy;
    float Bragg2SphereEnergy;
    float PandoraSphereEnergy;

    void resetEvent();
    void resetMuon();

    bool IsUpright(recob::Track const& T);
};

ana::TagAna::TagAna(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}, MichelAnalyzer{p},
    fLog(p.get<bool>("Log", true)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 20.F)), // in cm
    fMichelRadius(p.get<float>("MichelRadius", 20.F)), //in cm
    fNearbyRadius(p.get<float>("NearbyRadius", 40.F)), //in cm
    fBodyDistance(p.get<float>("BodyDistance", 20.F)), //in cm
    fRegN(p.get<unsigned>("RegN", 6)),
    fBraggThreshold(p.get<float>("BraggThreshold", 1.7)) // in MIP dE/dx
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
    wireWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
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
    fCathodeGap = geoHighX.x.min - geoLowX.x.max;

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX << std::endl
        << "  LowX Bounds: " << geoLowX << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl
        << "  Michel Space Radius: " << fMichelRadius << " cm";

    tEvent = asFile->make<TTree>("event","");

    tEvent->Branch("eventRun", &evRun);
    tEvent->Branch("eventSubRun", &evSubRun);
    tEvent->Branch("eventEvent", &evEvent);
    tEvent->Branch("isReal", &EventIsReal);

    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    EventHits.SetBranches(tEvent);

    tMuon = asFile->make<TTree>("muon","");

    tMuon->Branch("eventRun", &evRun);
    tMuon->Branch("eventSubRun", &evSubRun);
    tMuon->Branch("eventEvent", &evEvent);

    tMuon->Branch("iEvent", &iEvent);
    tMuon->Branch("iMuon", &iMuon);
    tMuon->Branch("iMuonInEvent", &EventNMuon);

    tMuon->Branch("TagIsUpright", &TagIsUpright);
    tMuon->Branch("TagCathodeCrossing", &TagCathodeCrossing);
    tMuon->Branch("TagAnodeCrossing", &TagAnodeCrossing);
    tMuon->Branch("CutTrackLength", &CutTrackLength);
    tMuon->Branch("TagEndInWindow", &TagEndInWindow);
    tMuon->Branch("TagEndInVolume", &TagEndInVolume);
    tMuon->Branch("TagHitCathodeCrossing", &TagHitCathodeCrossing);
    tMuon->Branch("TagHitGoodCathodeCrossing", &TagHitGoodCathodeCrossing);
    tMuon->Branch("CutdQdxMax", &CutdQdxMax);
    tMuon->Branch("TagBraggError", &TagBraggError);

    tMuon->Branch("AgnoTagTrkEndLowN", &AgnoTagTrkEndLowN);
    tMuon->Branch("AgnoTagTrkEndBadCC", &AgnoTagTrkEndBadCC);
    tMuon->Branch("AgnoTagEndHitError", &AgnoTagEndHitError);

    tMuon->Branch("TrueTagPdg", &TrueTagPdg);
    tMuon->Branch("TrueTagDownward", &TrueTagDownward);
    tMuon->Branch("TrueTagEndProcess", &TrueTagEndProcess);
    tMuon->Branch("TrueTagHasMichel", &TrueTagHasMichel);

    MuonHits.SetBranches(tMuon);
    MuonEndHit.SetBranches(tMuon, "End");
    MuonEndPoint.SetBranches(tMuon, "End");
    MuonTrueEndHit.SetBranches(tMuon, "TrueEnd");
    MuonTrueEndPoint.SetBranches(tMuon, "TrueEnd");
    BraggEndHit.SetBranches(tMuon, "BraggEnd");
    BraggMuonHits.SetBranches(tMuon, "BraggMuon");

    tMuon->Branch("MichelTrackLength", &MichelTrackLength); // cm
    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy); // MeV
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy); // ADC
    tMuon->Branch("BraggSphereEnergy", &BraggSphereEnergy); // ADC
    tMuon->Branch("Bragg2SphereEnergy", &Bragg2SphereEnergy); // ADC
    tMuon->Branch("PandoraSphereEnergy", &PandoraSphereEnergy); // ADC

    MichelHits.SetBranches(tMuon, "Michel");
}

void ana::TagAna::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) return;
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) return;
    VecPtrTrk vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    resetEvent();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    EventIsReal = e.isRealData();

    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            EventHits.push_back(GetHit(p_hit));

    // loop over tracks to find muons
    for (PtrTrk const& p_trk : vp_trk) {
        if (fLog) std::cout << "e" << iEvent << "t" << p_trk->ID() << "\r" << std::flush;
        resetMuon();
        AgnoCountAll++;

        VecPtrHit vph_mu = fmp_trk2hit.at(p_trk.key());
        if (vph_mu.empty()) AgnoCountNoHit++;
        ASSERT(vph_mu.size())

        CutTrackLength = p_trk->Length();

        TagIsUpright =  IsUpright(*p_trk);
        geo::Point_t Start = TagIsUpright ? p_trk->Start() : p_trk->End();
        geo::Point_t End = TagIsUpright ? p_trk->End() : p_trk->Start();

        TagEndInVolume = geoHighX.isInsideYZ(End, fMichelRadius);

        TagCathodeCrossing = (
            geoLowX.isInside(Start)
            && geoHighX.isInside(End)
        ) || (
            geoHighX.isInside(Start)
            && geoLowX.isInside(End)
        );

        // Anode Crossing: SUPPOSITION: Muon is downward
        if (geoDet == kPDVD)
            TagAnodeCrossing = geoHighX.isInsideYZ(Start, fMichelRadius);
        else if (geoDet == kPDHD)
            TagAnodeCrossing = geoHighX.isInsideYZ(Start, fMichelRadius) || geoLowX.isInsideYZ(Start, fMichelRadius);


        



        // Compare to Agnochecks
        std::vector<ana::LinearRegression> reg_side(2);
        std::vector<VecPtrHit> vph_side(2);
        for (PtrHit const& ph_hit : vph_mu) {
            if (ph_hit->View() != geo::kW) continue;
            int side = ana::tpc2side[geoDet][ph_hit->WireID().TPC];
            if (side == -1) continue;
            double z = GetSpace(ph_hit->WireID());
            double t = ph_hit->PeakTime() * fTick2cm;
            reg_side[side].add(z, t);
            vph_side[side].push_back(ph_hit);
        }
        AgnoTagTrkEndLowN = reg_side[0].n < ana::LinearRegression::nmin
            && reg_side[1].n < ana::LinearRegression::nmin;

        if (reg_side[0].n >= ana::LinearRegression::nmin
            && reg_side[1].n >= ana::LinearRegression::nmin
        ) {
            reg_side[0].compute();
            reg_side[1].compute();
            std::vector<std::pair<VecPtrHit::iterator, VecPtrHit::iterator>> ends_side(2);
            ends_side[0] = std::minmax_element(
                vph_side[0].begin(), vph_side[0].end(),
                [&](PtrHit const& a, PtrHit const& b) -> bool {
                    double sa = reg_side[0].projection(GetSpace(a->WireID()), a->PeakTime() * fTick2cm);
                    double sb = reg_side[0].projection(GetSpace(b->WireID()), b->PeakTime() * fTick2cm);
                    return sa < sb;
                } 
            );
            ends_side[1] = std::minmax_element(
                vph_side[1].begin(), vph_side[1].end(),
                [&](PtrHit const& a, PtrHit const& b) -> bool {
                    double sa = reg_side[1].projection(GetSpace(a->WireID()), a->PeakTime() * fTick2cm);
                    double sb = reg_side[1].projection(GetSpace(b->WireID()), b->PeakTime() * fTick2cm);
                    return sa < sb;
                } 
            );
            std::vector<std::pair<PtrHit, PtrHit>> pairs = {
                { *ends_side[0].first, *ends_side[1].first },
                { *ends_side[0].first, *ends_side[1].second },
                { *ends_side[0].second, *ends_side[1].first },
                { *ends_side[0].second, *ends_side[1].second }
            };
            AgnoTagTrkEndBadCC = std::find_if(
                pairs.begin(), pairs.end(),
                [&](std::pair<PtrHit, PtrHit> const& p) -> bool {
                    return GetDistance(p.first, p.second) < 2 * fCathodeGap;
                }
            ) == pairs.end();
        } else {
            AgnoTagTrkEndBadCC = false;
        }
        // End of Agnochecks



        if (fLog) std::cout << "\t" "\033[1;93m" "e" << iEvent << "m" << EventNMuon << " (" << iMuon << ")" "\033[0m" << std::endl;
        EventiMuon.push_back(iMuon);


        int dirz = End.Z() > Start.Z() ? 1 : -1;
        ana::SortedHits sh_mu = GetSortedHits(vph_mu, dirz);

        TagHitCathodeCrossing = sh_mu.is_cc();
        TagHitGoodCathodeCrossing = sh_mu.is_cc() && GetDistance(sh_mu.cc.first, sh_mu.cc.second) < 2 * fCathodeGap;

        if (sh_mu) {
            AgnoTagEndHitError = false; 
            // ASSERT(sh_mu)
            // bool TagHitCathodeCrossing = sh_mu.is_cc();

            // Last Hit: SUPPOSITION: Muon is downward
            // MuonEndHit = GetHit(sh_mu.lastHit(dirz));
            MuonEndHit = GetHit(sh_mu.end);
            MuonEndPoint = ana::Point(End);

            TagEndInWindow = wireWindow.isInside(MuonEndHit.tick, fMichelRadius / fTick2cm);

            ana::Bragg bragg = GetBragg(
                sh_mu.vph,
                sh_mu.end,
                p_trk,
                vph_ev,
                fop_hit2trk,
                { fBodyDistance, fRegN, fTrackLengthCut, fNearbyRadius }
            );
            TagBraggError = bragg.error;

            BraggMuonHits.clear();
            if (bragg) {
                CutdQdxMax = bragg.max_dQdx / bragg.mip_dQdx;
                BraggEndHit = GetHit(bragg.end);
                for (PtrHit const& ph_bragg_mu : bragg.vph_clu)
                    BraggMuonHits.push_back(GetHit(ph_bragg_mu));
            } else {
                CutdQdxMax = 0.F;
                BraggEndHit = ana::Hit{};
            }

            // getting all muon hits
            for (PtrHit const& ph_mu : sh_mu.vph)
                if (ph_mu->View() == geo::kW)
                    MuonHits.push_back(GetHit(ph_mu));

            VecPtrHit::iterator iph_bragg = std::find_if(
                bragg.vph_clu.begin(), bragg.vph_clu.end(),
                [&](PtrHit const& h) -> bool { return h.key() == bragg.end.key(); }
            );
            if (iph_bragg != bragg.vph_clu.end()) iph_bragg++;
            // integrate charges around muon endpoint
            for (PtrHit const& ph_ev : vph_ev) {
                if (ph_ev->View() != geo::kW) continue;
                if (ana::tpc2sec[geoDet][ph_ev->WireID().TPC]
                    != MuonEndHit.section) continue;

                PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());

                // ??????????????????????
                // art::Ptr<recob::Shower> ps_hit = fop_hit2shw.at(ph_ev.key());
                // if (ps_hit) continue;

                if ((
                    GetDistance(ph_ev, sh_mu.end) <= fMichelRadius
                ) && (
                    !pt_hit || pt_hit->Length() < fTrackLengthCut
                )) {
                    PandoraSphereEnergy += ph_ev->Integral();
                } 

                if ((
                    bragg
                ) && (
                    GetDistance(ph_ev, bragg.end) <= fMichelRadius
                ) && (
                    !pt_hit || pt_hit.key() == p_trk.key() || pt_hit->Length() < fTrackLengthCut
                ) && (
                    std::find_if(
                        bragg.vph_clu.begin(), iph_bragg,
                        [&](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                    ) == iph_bragg
                )) {
                    BraggSphereEnergy += ph_ev->Integral();
                }
            }
            Bragg2SphereEnergy = std::accumulate(
                iph_bragg, bragg.vph_clu.end(), 0.F,
                [&](float sum, PtrHit const& h) -> float {
                    if (GetDistance(h, bragg.end) < fMichelRadius)
                        return sum + h->Integral(); 
                    else 
                        return sum;
                }
            );
        } else AgnoTagEndHitError = true;
        
        // Truth Information
        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        if (mcp) {
            TrueTagPdg = mcp->PdgCode();
            TrueTagEndProcess = mcp->EndProcess();

            if (geoDet == kPDVD)
                TrueTagDownward = mcp->Position(0).X() > mcp->EndPosition().X();
            else if (geoDet == kPDHD)
                TrueTagDownward = mcp->Position(0).Y() > mcp->EndPosition().Y();

            VecPtrHit vph_mcp_mu;
            vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            ana::SortedHits sh_mcp = GetSortedHits(vph_mcp_mu, mcp->EndZ() > mcp->Vz() ? 1 : -1);
            if (sh_mcp) 
                MuonTrueEndHit = GetHit(sh_mcp.end);
            
            MuonTrueEndPoint = ana::Point(mcp->EndPosition().Vect());

            // a decaying muon has nu_mu, nu_e and elec as last daughters
            simb::MCParticle const* mcp_michel = GetMichelMCP(mcp);
            if (mcp_michel) {
                TrueTagHasMichel = (
                    geoHighX.isInside(mcp_michel->Position().Vect(), 20.F)
                    || geoLowX.isInside(mcp_michel->Position().Vect(), 20.F)
                ) ? kHasMichelFiducial : (
                    geoHighX.isInside(mcp_michel->Position().Vect())
                    || geoLowX.isInside(mcp_michel->EndPosition().Vect())
                    ? kHasMichelInside
                    : kHasMichelOutside
                );

                PtrTrk trk_michel = ana::mcp2trk(mcp_michel, vp_trk, clockData, fmp_trk2hit);
                MichelTrackLength = trk_michel ? trk_michel->Length() : 0;
                MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;

                VecPtrHit vph_michel = ana::mcp2hits(mcp_michel, vph_ev, clockData, true);
                for (PtrHit const& ph_michel : vph_michel)
                    if (ph_michel->View() == geo::kW)
                        MichelHits.push_back(GetHit(ph_michel));
                MichelHitEnergy = MichelHits.energy();
            }
        }
        tMuon->Fill();
        iMuon++;
        EventNMuon++;
    } // end of loop over tracks
    tEvent->Fill();
    iEvent++;
}

void ana::TagAna::beginJob() {}
void ana::TagAna::endJob() {

    std::cout << "DEBUG: ===================" << std::endl;
    std::cout << "  Tracks found: " << AgnoCountAll << std::endl;
    std::cout << "  Tracks without hits: " << AgnoCountNoHit << std::endl;
    std::cout << "  Tracks with < " << ana::LinearRegression::nmin << " hits on each sec: " << AgnoCountNoSec << std::endl;
    std::cout << "==========================" << std::endl;

}

void ana::TagAna::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
}
void ana::TagAna::resetMuon() {
    MuonHits.clear();
    MichelTrackLength = 0;
    MichelHits.clear();
    BraggMuonHits.clear();
    MichelHitEnergy = 0;
    BraggSphereEnergy = 0;
    Bragg2SphereEnergy = 0;
    PandoraSphereEnergy = 0;

    // Truth Information
    TrueTagPdg = 0;
    MuonTrueEndPoint = ana::Point{};
    TrueTagEndProcess = "";
    TrueTagDownward = false;
    MuonTrueEndHit = ana::Hit{};
    MuonTrueEndPoint = ana::Point{};
    // Michel Information
    TrueTagHasMichel = kNoMichel;
    MichelTrackLength = -1;
    MichelTrueEnergy = -1;
    MichelHitEnergy = -1;
}

bool ana::TagAna::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}

DEFINE_ART_MODULE(ana::TagAna)