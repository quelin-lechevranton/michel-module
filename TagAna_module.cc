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

using PtrHit = art::Ptr<recob::Hit>;
using VecPtrHit = std::vector<art::Ptr<recob::Hit>>;

class ana::TagAna : public art::EDAnalyzer, private MichelAnalyzer {
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
    // Detector Properties
    float fADC2MeV;
    float fTick2cm; // cm/tick
    float fSamplingRate; // µs/tick
    float fDriftVelocity; // cm/µs
    float fCathodeGap; // cm

    bounds3D<float> geoHighX, geoLowX;
    bounds<float> wireWindow;

    // Data Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag tag_mcp, tag_sed, tag_wir,
        tag_hit, tag_clu, tag_trk,
        tag_shw, tag_spt, tag_pfp;

    // Input Parameters
    bool fLog;
    bool fKeepOutside;
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

    float CutdQdxMax;
    enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
    int TagBraggError;
    ana::Hits BraggMuonHits;

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
    float PandoraSphereEnergy;

    void resetEvent();
    void resetMuon();

    bool IsUpright(recob::Track const& T);
    PtrHit GetBraggEnd(
        VecPtrHit const& vph_trk,
        PtrHit const& ph_trk_end,
        PtrTrk const& p_trk,
        VecPtrHit const& vph_ev,
        art::FindOneP<recob::Track> const& fop_hit2trk,
        VecPtrHit *vph_sec_bragg = nullptr,
        float *max_dQdx = nullptr,
        int *error = nullptr
    );
};

ana::TagAna::TagAna(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fLog(p.get<bool>("Log", false)),
    fKeepOutside(p.get<bool>("KeepOutside", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 40.F)), // in cm
    fMichelRadius(p.get<float>("MichelRadius", 20.F)), //in cm
    fNearbyRadius(p.get<float>("NearbyRadius", 40.F)), //in cm
    fBodyDistance(p.get<float>("BodyDistance", 20.F)), //in cm
    fRegN(p.get<unsigned>("RegN", 6)),
    fBraggThreshold(p.get<float>("BraggThreshold", 1.7)) // in MIP dE/dx
{
    asGeo = &*art::ServiceHandle<geo::Geometry>{};
    asWire = &art::ServiceHandle<geo::WireReadout>{}->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>{};    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>{};
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    for (std::vector<std::string> prod : vvsProducts) {
        const std::string process  = prod[0],
                          label    = prod[1],
                          instance = prod[2],
                          type     = prod[3];

        const art::InputTag tag = art::InputTag{label,instance};

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::Shower")           tag_shw = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
        else if (type == "recob::PFParticle")       tag_pfp = tag;
    }

    if (asGeo->DetectorName().find("vd") != std::string::npos)
        geoDet = kPDVD;
    else if (asGeo->DetectorName().find("hd") != std::string::npos)
        geoDet = kPDHD;
    else {
        std::cout << "\033[1;91m" "unknown geometry: "
            << asGeo->DetectorName() << "\033[0m" << std::endl;
        exit(1);
    }

    // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    // fADC2MeV = (geoDet == kPDVD ? 200 : ) * 23.6 * 1e-6 / 0.7;
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    wireWindow = bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoLowX = bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(), 
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoHighX = bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoLowX = bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoHighX = bounds3D<float>{
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
        << "  Sampling Rate: " << fSamplingRate << " µs/tick" << std::endl
        << "  Drift Velocity: " << fDriftVelocity << " cm/µs" << std::endl
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
    tMuon->Branch("CutdQdxMax", &CutdQdxMax);
    tMuon->Branch("TagBraggError", &TagBraggError);

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
    std::vector<PtrTrk> vp_trk;
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

        VecPtrHit vph_muon = fmp_trk2hit.at(p_trk.key());
        ASSERT(vph_muon.size())

        CutTrackLength = p_trk->Length();

        TagIsUpright =  IsUpright(*p_trk);
        geo::Point_t Start = TagIsUpright ? p_trk->Start() : p_trk->End();
        geo::Point_t End = TagIsUpright ? p_trk->End() : p_trk->Start();

        TagEndInVolume = geoHighX.isInsideYZ(End, 20.F);

        TagCathodeCrossing = (
            geoLowX.isInside(Start)
            && geoHighX.isInside(End)
        ) || (
            geoHighX.isInside(Start)
            && geoLowX.isInside(End)
        );

        // Anode Crossing: SUPPOSITION: Muon is downward
        if (geoDet == kPDVD)
            TagAnodeCrossing = geoHighX.isInsideYZ(Start, 20.);
        else if (geoDet == kPDHD)
            TagAnodeCrossing = geoHighX.isInsideYZ(Start, 20.) || geoLowX.isInsideYZ(Start, 20.);

        ana::SortedHits sh_muon = GetSortedHits(vph_muon);

        ASSERT(sh_muon)
        // bool TagHitCathodeCrossing = sh_muon.isCathodeCrossing();

        // Last Hit: SUPPOSITION: Muon is downward
        int dirz = End.Z() > Start.Z() ? 1 : -1;
        MuonEndHit = GetHit(sh_muon.lastHit(dirz));
        MuonEndPoint = ana::Point(End);

        TagEndInWindow = wireWindow.isInside(MuonEndHit.tick, fMichelRadius / fTick2cm);

        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        ASSERT(mcp)

        // we found a muon candidate!
        if (fLog) std::cout << "\t" "\033[1;93m" "e" << iEvent << "m" << EventNMuon << " (" << iMuon << ")" "\033[0m" << std::endl;
        EventiMuon.push_back(iMuon);
        resetMuon();

        VecPtrHit vph_bragg_muon;
        PtrHit ph_bragg = GetBraggEnd(
            sh_muon.vph, 
            sh_muon.lastHit(dirz),
            p_trk,
            vph_ev,
            fop_hit2trk,
            &vph_bragg_muon,
            &CutdQdxMax,
            &TagBraggError
        );
        BraggEndHit = ph_bragg ? GetHit(ph_bragg) : ana::Hit{};
        for (PtrHit const& ph_bragg_muon : vph_bragg_muon)
            BraggMuonHits.push_back(GetHit(ph_bragg_muon));

        // getting all muon hits
        for (PtrHit const& vph_muon : sh_muon.vph)
            if (vph_muon->View() == geo::kW)
                MuonHits.push_back(GetHit(vph_muon));

        // integrate charges around muon endpoint
        for (PtrHit const& ph_ev : vph_ev) {
            if (ph_ev->View() != geo::kW) continue;
            if (ana::tpc2sec[geoDet][ph_ev->WireID().TPC]
                != ana::tpc2sec[geoDet][MuonEndHit.tpc]) continue;
            ana::Hit hit = GetHit(ph_ev);

            float dz = (hit.space - MuonEndHit.space);
            float dt = (hit.tick - MuonEndHit.tick) * fTick2cm;
            float dr2 = dz*dz + dt*dt;

            if (dr2 > fMichelRadius * fMichelRadius) continue;

            // check if the hit is in an other track
            PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
            if (pt_hit
                && pt_hit.key() != p_trk.key()
                && pt_hit->Length() > fTrackLengthCut
            ) continue;

            // ??????????????????????
            art::Ptr<recob::Shower> ps_hit = fop_hit2shw.at(ph_ev.key());
            if (ps_hit) continue;

            if (std::find_if(
                vph_bragg_muon.begin(), 
                vph_bragg_muon.end(),
                [&](PtrHit const& h) -> bool {
                    return h.key() == ph_ev.key();
                }) == vph_bragg_muon.end()
            ) {
                BraggSphereEnergy += hit.adc;
            }

            if (!pt_hit) {
                PandoraSphereEnergy += hit.adc;
            }
        }
        
        // Truth Information

        TrueTagPdg = mcp->PdgCode();
        TrueTagEndProcess = mcp->EndProcess();

        if (geoDet == kPDVD)
            TrueTagDownward = mcp->Position(0).X() > mcp->EndPosition().X();
        else if (geoDet == kPDHD)
            TrueTagDownward = mcp->Position(0).Y() > mcp->EndPosition().Y();

        VecPtrHit vph_ev_mcp_muon;
        if (mcp) {
            vph_ev_mcp_muon = ana::mcp2hits(mcp, vph_ev, clockData, false);
            // vph_ev_mcp_muon = GetSortedHits(
            //     vph_ev_mcp_muon, 
            //     mcp->EndZ() > mcp->Vz() ? 1 : -1
            // );
            ana::SortedHits sh_mcp = GetSortedHits(vph_ev_mcp_muon);
            if (sh_mcp) 
                MuonTrueEndHit = GetHit(sh_mcp.lastHit(mcp->EndZ() > mcp->Vz() ? 1 : -1));
            
            MuonTrueEndPoint = ana::Point(mcp->EndPosition().Vect());
        }

        // a decaying muon has nu_mu, nu_e and elec as last daughters
        simb::MCParticle const* mcp_michel = nullptr;
        VecPtrHit vph_ev_mcp_michel;
        if (mcp && mcp->NumberDaughters() >= 3) {
            bool has_numu = false, has_nue = false;
            for (
                int i_dau=mcp->NumberDaughters()-3;
                i_dau<mcp->NumberDaughters();
                i_dau++
            ) {
                simb::MCParticle const* mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    
                if (!mcp_dau) continue;
                switch (abs(mcp_dau->PdgCode())) {
                    case 14: has_numu = true; break;
                    case 12: has_nue = true; break;
                    case 11: mcp_michel = mcp_dau; break;
                    default: break;
                }
            }

            if (mcp_michel and has_numu and has_nue) {
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

                vph_ev_mcp_michel = ana::mcp2hits(mcp_michel, vph_ev, clockData, true);
                for (PtrHit p_hit_michel : vph_ev_mcp_michel)
                    if (p_hit_michel->View() == geo::kW)
                        MichelHits.push_back(GetHit(p_hit_michel));
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
void ana::TagAna::endJob() {}

void ana::TagAna::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
}
void ana::TagAna::resetMuon() {
    MuonHits.clear();
    MuonTrueEndHit = ana::Hit{};
    MuonTrueEndPoint = ana::Point{};
    MichelTrackLength = 0;
    MichelTrueEnergy = 0;
    MichelHits.clear();
    MichelHitEnergy = 0;
    BraggMuonHits.clear();

    TrueTagHasMichel = kNoMichel;
    MichelTrackLength = -1;
    MichelHitEnergy = -1;
}

bool ana::TagAna::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}

PtrHit ana::TagAna::GetBraggEnd(
    VecPtrHit const& vph_trk,
    PtrHit const& ph_trk_end,
    PtrTrk const& p_trk,
    VecPtrHit const& vph_ev,
    art::FindOneP<recob::Track> const& fop_hit2trk,
    VecPtrHit *vph_sec_bragg,
    float *max_dQdx,
    int *error
) {
    VecPtrHit vph_sec_trk;
    for (PtrHit const& p_hit : vph_trk) {
        if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != ana::tpc2sec[geoDet][ph_trk_end->WireID().TPC]) continue;
        vph_sec_trk.push_back(p_hit);
    }
    if (vph_sec_trk.empty()
        || std::find_if(
            vph_sec_trk.begin(), vph_sec_trk.end(),
            [key=ph_trk_end.key()](PtrHit const& ph) -> bool {
                return ph.key() == key;
            }
        ) == vph_sec_trk.end()
    ) {
        if (error) *error = kEndNotFound;
        return PtrHit{};
    }

    auto dist2 = [&](PtrHit const& ph1, PtrHit const& ph2) -> double {
        return pow((ph1->PeakTime() - ph2->PeakTime()) * fTick2cm, 2)
            + pow(GetSpace(ph1->WireID()) - GetSpace(ph2->WireID()), 2);
    };
    std::sort(
        vph_sec_trk.begin(),
        vph_sec_trk.end(),
        [&](PtrHit const& ph1, PtrHit const& ph2) -> bool {
            return dist2(ph2, ph_trk_end) > dist2(ph1, ph_trk_end);
        }
    );

    VecPtrHit::iterator iph_body = std::find_if(
        vph_sec_trk.begin(),
        vph_sec_trk.end(),
        [&](PtrHit const& h) -> bool {
            return dist2(h, ph_trk_end) > fBodyDistance * fBodyDistance;
        }
    );
    if (std::distance(iph_body, vph_sec_trk.end()) < fRegN) {
        if (error) *error = kSmallBody;
        return PtrHit{};
    }

    VecPtrHit vph_reg{iph_body, iph_body+fRegN};
    auto orientation = [&](VecPtrHit const& vph) -> std::pair<double, double> {
        LinearRegression reg;
        for (PtrHit const& ph : vph) {
            double z = GetSpace(ph->WireID());
            double t = ph->PeakTime() * fTick2cm;
            reg.add(z, t);
        }
        reg.compute();
        // DEBUG(reg.corr == 0)
        int dirz = GetSpace(vph.back()->WireID())
            > GetSpace(vph.front()->WireID())
            ? 1 : -1;
        double sigma = TMath::Pi() / 4 / reg.corr;
        double theta = reg.theta(dirz);
        return std::make_pair(theta, sigma);
    };
    std::pair<double, double> reg = orientation(vph_reg);

    VecPtrHit vph_near;
    for (PtrHit const& ph_ev : vph_ev) {
        if (ana::tpc2sec[geoDet][ph_ev->WireID().TPC]
            != ana::tpc2sec[geoDet][ph_trk_end->WireID().TPC]) continue;

        // check if the hit is in a track
        PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
        if (pt_hit) {
            // check if the hit is on the body of the track
            if (pt_hit.key() == p_trk.key()
                && std::find_if(
                    iph_body, vph_sec_trk.end(),
                    [key=ph_ev.key()](PtrHit const& ph) -> bool {
                        return ph.key() == key;
                    }
                ) != vph_sec_trk.end()
            ) continue;
            // check if the hit is on a long track
            if (
                pt_hit.key() != p_trk.key()
                && pt_hit->Length() > fTrackLengthCut
            ) continue;
        }

        // check if the hit is close enough
        if (dist2(ph_ev, ph_trk_end) > fNearbyRadius) continue;

        vph_near.push_back(ph_ev);
    }
    DEBUG(vph_near.empty())

    auto score = [&](PtrHit const& ph1, PtrHit const& ph2, double theta, double sigma) -> double {
        double dz = GetSpace(ph2->WireID()) - GetSpace(ph1->WireID());
        double dt = (ph2->PeakTime() - ph1->PeakTime()) * fTick2cm;
        double da = atan2(dt, dz) - theta;
        da = abs(da) < TMath::Pi() ? da : da - (da>0?1:-1)*2*TMath::Pi();
        double r = sqrt(dt*dt + dz*dz);

        return TMath::Gaus(da, 0, sigma) / r;
    };

    VecPtrHit vph_sec{vph_reg.begin(), vph_reg.end()};
    std::reverse(vph_sec.begin(), vph_sec.end());
    PtrHit ph_prev = ph_trk_end;
    while (vph_near.size()) {
        VecPtrHit::iterator iph_max = std::max_element(
            vph_near.begin(), vph_near.end(),
            [&](PtrHit const& ph1, PtrHit const& ph2) -> bool {
                return score(ph_prev, ph1, reg.first, reg.second)
                    < score(ph_prev, ph2, reg.first, reg.second);
            }
        );

        vph_near.erase(iph_max);
        vph_reg.insert(vph_reg.begin(), *iph_max);
        vph_reg.pop_back();
        vph_sec.push_back(*iph_max);
        ph_prev = *iph_max;

        reg = orientation(vph_reg);

        // reg = LinearRegression{};
        // for (PtrHit const& ph : vph_reg) {
        //     double z = GetSpace(ph->WireID());
        //     double t = ph->PeakTime() * fTick2cm;
        //     reg.add(z, t);
        // }
        // reg.compute();
        // DEBUG(reg.corr == 0)
        // dirz = GetSpace(vph_reg.back()->WireID())
        //     > GetSpace(vph_reg.front()->WireID())
        //     ? 1 : -1;
        // sigma = TMath::Pi() / 4 / reg.corr;
        // theta = reg.theta(dirz);
    }

    unsigned const trailing_radius = 6;
    double max = std::numeric_limits<double>::lowest();
    PtrHit ph_max;
    for (auto iph_sec=vph_sec.begin(); iph_sec!=vph_sec.end(); iph_sec++) {
        VecPtrHit::iterator jph_sec = std::distance(iph_sec, vph_sec.end()) > trailing_radius
            ? iph_sec+trailing_radius : vph_sec.end();
        unsigned l = std::distance(iph_sec, jph_sec);

        double dQ = std::accumulate(
            iph_sec, jph_sec, 0.,
            [](double sum, PtrHit const& ph) {
                return sum+ph->Integral();
            }
        );
        dQ /= l;

        double dx = iph_sec == vph_sec.begin()
            ? sqrt(dist2(*iph_sec, *(iph_sec+1)))
            : ( iph_sec == vph_sec.end()-1
                ? sqrt(dist2(*(iph_sec-1), *iph_sec))
                : .5*sqrt(dist2(*(iph_sec-1), *(iph_sec+1)))
            );
        for (auto iph=iph_sec+1; iph!=jph_sec; iph++)
            dx += iph == vph_sec.end()-1
                ? sqrt(dist2(*(iph-1), *iph))
                : .5*sqrt(dist2(*(iph-1), *(iph+1)));
        dx /= l;

        double dQdx = dQ / dx;
        if (dQdx > max) {
            max = dQdx;
            ph_max = *iph_sec;
        }
    }
    if (vph_sec_bragg) *vph_sec_bragg = vph_sec;
    if (max_dQdx) *max_dQdx = float(max);
    if (error) *error = kNoError;
    return ph_max;
}   



DEFINE_ART_MODULE(ana::TagAna)