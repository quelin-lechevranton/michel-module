////////////////////////////////////////////////////////////////////////
// Class:       Tagchecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Agnochecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"

namespace ana {
    class Tagchecks;
}

using HitPtr = art::Ptr<recob::Hit>;
using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

class ana::Tagchecks : public art::EDAnalyzer {
public:
    explicit Tagchecks(fhicl::ParameterSet const& p);
    Tagchecks(Tagchecks const&) = delete;
    Tagchecks(Tagchecks&&) = delete;
    Tagchecks& operator=(Tagchecks const&) = delete;
    Tagchecks& operator=(Tagchecks&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:

    // Utilities
    art::ServiceHandle<art::TFileService> asFile;
    art::ServiceHandle<cheat::ParticleInventoryService> asPartInv;
    art::ServiceHandle<cheat::BackTrackerService> asBackTrack;

    const geo::GeometryCore* asGeo;
    const geo::WireReadoutGeom* asWire;
    const detinfo::DetectorPropertiesService* asDetProp;
    const detinfo::DetectorClocksService* asDetClocks;

    int geoDet;
    enum EnumDet { kPDVD, kPDHD };

    // Detector Properties
    float fADC2MeV;
    float fTick2cm; // cm/tick
    float fSamplingRate; // µs/tick
    float fDriftVelocity; // cm/µs
    float fCathodeGap; // cm

    geo::BoxBoundedGeo geoHighX, geoLowX;
    bounds<float> wireWindow;
    std::map<geo::PlaneID, ana::axis> plane2axis;
    std::map<geo::PlaneID, double> plane2pitch;

    // Data Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag tag_mcp, tag_sed, tag_wir,
        tag_hit, tag_clu, tag_trk,
        tag_spt, tag_pfp, tag_r3d;

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
    enum BraggError { kNoError, kEndNotFound, kSmallBody };
    BraggError TagBraggError;
    ana::Hits BraggMuonHits;

    bool TrueTagDownward;
    int TrueTagPdg;
    std::string TrueTagEndProcess;

    int TrueTagHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside };

    ana::Hits MuonHits;
    ana::Hit MuonEndHit;
    ana::Hit MuonTrueEndHit;
    ana::Hit BraggEndHit;

    float MichelTrackLength;
    ana::Hits MichelHits;
    float MichelTrueEnergy, MichelHitEnergy;

    void resetEvent();
    void resetMuon();

    bool IsUpright(recob::Track const& T);
    double GetSpace(geo::WireID);
    ana::Hit GetHit(HitPtr const p_hit);
    HitPtrVec GetSortedHits(
        HitPtrVec const& vp_hit,
        int dirz,
        HitPtrPair *pp_cathode_crossing = nullptr,
        HitPtrVec *vp_tpc_crossing = nullptr,
        std::vector<ana::LinearRegression> *p_side_reg = nullptr,
        geo::View_t view = geo::kW
    );
    std::pair<double, double> LinRegPCA(HitPtrVec const& vp_hit);
    HitPtr GetBraggEnd(
        HitPtrVec const& vph_trk,
        HitPtr const& ph_trk_end,
        art::Ptr<recob::Track> const& p_trk,
        HitPtrVec const& vph_ev,
        art::FindOneP<recob::Track> const& fop_hit2trk,
        HitPtrVec *vph_sec_bragg = nullptr,
        float *max_dQdx = nullptr,
        BraggError *error = nullptr
    );
};


ana::Tagchecks::Tagchecks(fhicl::ParameterSet const& p)
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
    // fADC2MeV = (geoDet == kPDVD ? 200 : 1000) * 23.6 * 1e-6 / 0.7;
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fTick2cm = fDriftVelocity * fSamplingRate;

    wireWindow = bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoLowX = geo::BoxBoundedGeo{
                asGeo->TPC(geo::TPCID{0, 0}).Min(), 
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoHighX = geo::BoxBoundedGeo{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoLowX = geo::BoxBoundedGeo{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoHighX = geo::BoxBoundedGeo{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        default: break;
    }
    fCathodeGap = geoHighX.MinX() - geoLowX.MaxX();

    for (unsigned t=0; t<asGeo->NTPC(); t++) {
        for (unsigned p=0; p<asWire->Nplanes(); p++) {
            unsigned c = 0; // cryostat
            geo::PlaneID pid{c, t, p};
            // assume constant wire pitch in one plane
            geo::WireGeo w0 = asWire->Wire(geo::WireID{pid, 0});
            geo::WireGeo w1 = asWire->Wire(geo::WireID{pid, 1});

            int dy = w1.GetCenter().Y() > w0.GetCenter().Y() ? 1 : -1;
            int dz = w1.GetCenter().Z() > w0.GetCenter().Z() ? 1 : -1;

            plane2axis[pid] = { dy * w0.CosThetaZ(), dz * w0.SinThetaZ() };
            plane2pitch[pid] = geo::WireGeo::WirePitch(w0, w1);
        }
    }

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
        << "  Sampling Rate: " << fSamplingRate << " µs/tick" << std::endl
        << "  Drift Velocity: " << fDriftVelocity << " cm/µs" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX.Min() << " -> " << geoHighX.Max() << std::endl
        << "  LowX Bounds: " << geoLowX.Min() << " -> " << geoLowX.Max() << std::endl;
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
    MuonTrueEndHit.SetBranches(tMuon, "TrueEnd");
    BraggEndHit.SetBranches(tMuon, "BraggEnd");
    BraggMuonHits.SetBranches(tMuon, "BraggMuon");

    tMuon->Branch("MichelTrackLength", &MichelTrackLength); // cm
    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy); // MeV
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy); // MeV

    MichelHits.SetBranches(tMuon, "Michel");
}

void ana::Tagchecks::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fTick2cm = fDriftVelocity * fSamplingRate;

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) return;
    HitPtrVec vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) return;
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    auto const & vh_pfp = e.getHandle<std::vector<recob::PFParticle>>(tag_pfp);
    if (!vh_pfp.isValid()) return;

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    art::FindOneP<recob::PFParticle> fop_trk2pfp(vh_trk, e, tag_trk);
    art::FindManyP<recob::SpacePoint> fmp_pfp2spt(vh_pfp, e, tag_pfp);

    resetEvent();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    EventIsReal = e.isRealData();

    for (HitPtr p_hit : vp_hit)
        if (p_hit->View() == geo::kW)
            EventHits.push_back(GetHit(p_hit));

    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        if (fLog) std::cout << "e" << iEvent << "t" << p_trk->ID() << "\r" << std::flush;

        HitPtrVec vph_muon = fmp_trk2hit.at(p_trk.key());
        ASSERT(vph_muon.size())

        CutTrackLength = p_trk->Length();

        TagIsUpright =  IsUpright(*p_trk);
        geo::Point_t Start = TagIsUpright ? p_trk->Start() : p_trk->End();
        geo::Point_t End = TagIsUpright ? p_trk->End() : p_trk->Start();

        TagEndInVolume = geoHighX.InFiducialY(End.Y(), 20.)
            && geoHighX.InFiducialZ(End.Z(), 20.);

        TagCathodeCrossing = (
            geoLowX.ContainsPosition(Start) 
            && geoHighX.ContainsPosition(End)
        ) || (
            geoHighX.ContainsPosition(Start)
            && geoLowX.ContainsPosition(End)
        );

        // Anode Crossing: SUPPOSITION: Muon is downward
        if (geoDet == kPDVD)
            TagAnodeCrossing = geoHighX.InFiducialY(Start.Y(), 20.)
                && geoHighX.InFiducialZ(Start.Z(), 20.);
        else if (geoDet == kPDHD)
            TagAnodeCrossing = (
                geoHighX.InFiducialY(Start.Y(), 20.)
                && geoHighX.InFiducialZ(Start.Z(), 20.)
            ) || (
                geoLowX.InFiducialY(Start.Y(), 20.)
                && geoLowX.InFiducialZ(Start.Z(), 20.)
            );

        // Last Hit: SUPPOSITION: Muon is downward
        HitPtrVec vph_muon_sorted = GetSortedHits(
            vph_muon, 
            End.Z() > Start.Z() ? 1 : -1
        );
        ASSERT(vph_muon_sorted.size())
        MuonEndHit = GetHit(vph_muon_sorted.back());

        TagEndInWindow = wireWindow.isInside(MuonEndHit.tick, fMichelRadius / fTick2cm);

        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        ASSERT(mcp)

        // we found a muon candidate!
        if (fLog) std::cout << "\t" "\033[1;93m" "e" << iEvent << "m" << EventNMuon << " (" << iMuon << ")" "\033[0m" << std::endl;
        EventiMuon.push_back(iMuon);
        resetMuon();

        HitPtrVec vph_bragg;
        HitPtr ph_bragg = GetBraggEnd(
            vph_muon_sorted, 
            vph_muon_sorted.back(),
            p_trk,
            vp_hit,
            fop_hit2trk,
            &vph_bragg,
            &CutdQdxMax,
            &TagBraggError
        );
        BraggEndHit = ph_bragg ? GetHit(ph_bragg) : ana::Hit{};
        for (HitPtr const& ph_bragg : vph_bragg)
            BraggMuonHits.push_back(GetHit(ph_bragg));

        // getting all muon hits
        for (HitPtr const& p_hit_muon : vph_muon_sorted)
            if (p_hit_muon->View() == geo::kW)
                MuonHits.push_back(GetHit(p_hit_muon));

        TrueTagPdg = mcp->PdgCode();
        TrueTagEndProcess = mcp->EndProcess();

        if (geoDet == kPDVD)
            TrueTagDownward = mcp->Position(0).X() > mcp->EndPosition().X();
        else if (geoDet == kPDHD)
            TrueTagDownward = mcp->Position(0).Y() > mcp->EndPosition().Y();

        HitPtrVec vp_hit_mcp_muon;
        MuonTrueEndHit = ana::Hit{};
        if (mcp) {
            vp_hit_mcp_muon = ana::mcp2hits(mcp, vp_hit, clockData, false);
            vp_hit_mcp_muon = GetSortedHits(
                vp_hit_mcp_muon, 
                mcp->EndZ() > mcp->Vz() ? 1 : -1
            );
            if (vp_hit_mcp_muon.size())
                MuonTrueEndHit = GetHit(vp_hit_mcp_muon.back());
        }

        // a decaying muon has nu_mu, nu_e and elec as last daughters
        simb::MCParticle const* mcp_michel = nullptr;
        HitPtrVec vp_hit_mcp_michel;
        if (mcp && mcp->NumberDaughters() >= 3) {
            bool has_numu = false, has_nue = false;
            for (int i_dau=mcp->NumberDaughters()-3; i_dau<mcp->NumberDaughters(); i_dau++) {
                simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    
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
                    geoHighX.ContainsPosition(mcp_michel->Position().Vect())
                    || geoLowX.ContainsPosition(mcp_michel->Position().Vect())
                )   ? kHasMichelInside
                    : kHasMichelOutside;

                art::Ptr<recob::Track> trk_michel = ana::mcp2trk(mcp_michel, vp_trk, clockData, fmp_trk2hit);
                MichelTrackLength = trk_michel ? trk_michel->Length() : 0;
                MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;

                vp_hit_mcp_michel = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);
                for (HitPtr p_hit_michel : vp_hit_mcp_michel)
                    if (p_hit_michel->View() == geo::kW)
                        MichelHits.push_back(GetHit(p_hit_michel));
                MichelHitEnergy = MichelHits.energy();
            } else TrueTagHasMichel = kNoMichel;
        } else {
            TrueTagHasMichel = -1;
            MichelTrackLength = -1;
            MichelHitEnergy = -1;
        }
        tMuon->Fill();
        iMuon++;
        EventNMuon++;
    } // end of loop over tracks
    tEvent->Fill();
    iEvent++;
}

void ana::Tagchecks::beginJob() {}
void ana::Tagchecks::endJob() {}


void ana::Tagchecks::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
}
void ana::Tagchecks::resetMuon() {
    MuonHits.clear();
    MichelTrackLength = 0;
    MichelTrueEnergy = 0;
    MichelHits.clear();
    MichelHitEnergy = 0;
}


double ana::Tagchecks::GetSpace(geo::WireID wid) {
    return plane2axis[wid].space(asWire->Wire(wid));
}

ana::Hit ana::Tagchecks::GetHit(HitPtr const p_hit) {
    geo::WireID wid = p_hit->WireID();
    // if (geoDet == kPDHD)
    //     for (int t : (int[]){0, 4, 3, 7})
    //         if (wireid.TPC == t)
    //             return ana::Hit{};
    
    return ana::Hit{
        wid.TPC,
        ana::tpc2sec[geoDet][wid.TPC],
        float(GetSpace(wid)),
        p_hit->Channel(),
        p_hit->PeakTime(),
        p_hit->Integral()
    };
}

bool ana::Tagchecks::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}

HitPtrVec ana::Tagchecks::GetSortedHits(
    HitPtrVec const& vp_hit,
    int dirz,
    HitPtrPair *pp_cathode_crossing,
    HitPtrVec *vp_tpc_crossing,
    std::vector<ana::LinearRegression> *p_side_reg,
    geo::View_t view
) {
    std::vector<ana::LinearRegression> side_reg(2);
    std::vector<HitPtrVec> side_hit(2);
    for (HitPtr const& p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        int side = ana::tpc2side[geoDet][p_hit->WireID().TPC];
        if (side == -1) continue;
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        side_reg[side].add(z, t);
        side_hit[side].push_back(p_hit);
    }
    if (
        side_reg[0].n < ana::LinearRegression::nmin 
        && side_reg[1].n < ana::LinearRegression::nmin
    ) return HitPtrVec{};
    for (ana::LinearRegression& reg : side_reg)
        if (reg.n >= ana::LinearRegression::nmin)
            reg.compute();
    if (p_side_reg) *p_side_reg = side_reg;
    for (int side=0; side<2; side++)
        if (side_reg[side].n >= ana::LinearRegression::nmin)
            std::sort(
                side_hit[side].begin(),
                side_hit[side].end(),
                [&, &reg=side_reg[side]](
                    HitPtr const& ph1, HitPtr const& ph2
                ) -> bool {
                    double const s1 = reg.projection(
                        GetSpace(ph1->WireID()),
                        ph1->PeakTime() * fTick2cm
                    );
                    double const s2 = reg.projection(
                        GetSpace(ph2->WireID()),
                        ph2->PeakTime() * fTick2cm
                    );
                    return (s2 - s1) * dirz > 0;
                }
            );
    if (side_reg[0].n < ana::LinearRegression::nmin)
        return side_hit[1];
    if (side_reg[1].n < ana::LinearRegression::nmin)
        return side_hit[0];
   
    std::pair<unsigned, unsigned> side_pair = 
        (side_reg[1].mx - side_reg[0].mx) * dirz > 0
        ? std::make_pair(0, 1) : std::make_pair(1, 0);
    
    if (pp_cathode_crossing) {
        pp_cathode_crossing->first = side_hit[side_pair.first].back();
        pp_cathode_crossing->second = side_hit[side_pair.second].front();
    }
    HitPtrVec vp_sorted_hit;
    if (vp_tpc_crossing) {
        vp_tpc_crossing->clear();
        HitPtr& prev_hit = side_hit[side_pair.first].front();
        unsigned prev_tpc = prev_hit->WireID().TPC;
        for (HitPtr const& p_hit : side_hit[side_pair.first]) {
            vp_sorted_hit.push_back(p_hit);
            if (p_hit->WireID().TPC != prev_tpc) {
                vp_tpc_crossing->push_back(prev_hit);
                vp_tpc_crossing->push_back(p_hit);
                prev_tpc = p_hit->WireID().TPC;
            }
            prev_hit = p_hit;
        }
        prev_hit = side_hit[side_pair.second].front();
        prev_tpc = prev_hit->WireID().TPC;
        for (HitPtr const& p_hit : side_hit[side_pair.second]) {
            vp_sorted_hit.push_back(p_hit);
            if (p_hit->WireID().TPC != prev_tpc) {
                vp_tpc_crossing->push_back(prev_hit);
                vp_tpc_crossing->push_back(p_hit);
                prev_tpc = p_hit->WireID().TPC;
            }
            prev_hit = p_hit;
        }
        return vp_sorted_hit;
    }
    vp_sorted_hit.insert(
        vp_sorted_hit.end(),
        side_hit[side_pair.first].begin(),
        side_hit[side_pair.first].end()
    );
    vp_sorted_hit.insert(
        vp_sorted_hit.end(),
        side_hit[side_pair.second].begin(),
        side_hit[side_pair.second].end()
    );
    return vp_sorted_hit;
}

HitPtr ana::Tagchecks::GetBraggEnd(
    HitPtrVec const& vph_trk,
    HitPtr const& ph_trk_end,
    art::Ptr<recob::Track> const& p_trk,
    HitPtrVec const& vph_ev,
    art::FindOneP<recob::Track> const& fop_hit2trk,
    HitPtrVec *vph_sec_bragg,
    float *max_dQdx,
    BraggError *error
) {
    HitPtrVec vph_sec_trk;
    for (HitPtr const& p_hit : vph_trk) {
        if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != ana::tpc2sec[geoDet][ph_trk_end->WireID().TPC]) continue;
        vph_sec_trk.push_back(p_hit);
    }
    if (vph_sec_trk.empty() || ph_trk_end != vph_sec_trk.back()) {
        if (error) *error = kEndNotFound;
        return HitPtr{};
    }

    auto dist2 = [&](HitPtr const& ph1, HitPtr const& ph2) -> double {
        return pow((ph1->PeakTime() - ph2->PeakTime()) * fTick2cm, 2)
            + pow(GetSpace(ph1->WireID()) - GetSpace(ph2->WireID()), 2);
    };
    std::sort(
        vph_sec_trk.begin(),
        vph_sec_trk.end(),
        [&](HitPtr const& ph1, HitPtr const& ph2) -> bool {
            return dist2(ph2, ph_trk_end) > dist2(ph1, ph_trk_end);
        }
    );

    HitPtrVec::iterator iph_body = std::find_if(
        vph_sec_trk.begin(),
        vph_sec_trk.end(),
        [&](HitPtr const& h) -> bool {
            return dist2(h, ph_trk_end) > fBodyDistance * fBodyDistance;
        }
    );
    if (std::distance(iph_body, vph_sec_trk.end()) < fRegN) {
        if (error) *error = kSmallBody;
        return HitPtr{};
    }

    HitPtrVec vph_reg{iph_body, iph_body+fRegN};
    auto orientation = [&](HitPtrVec const& vph) -> std::pair<double, double> {
        LinearRegression reg;
        for (HitPtr const& ph : vph) {
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

    HitPtrVec vph_near;
    for (HitPtr const& ph_ev : vph_ev) {
        if (ana::tpc2sec[geoDet][ph_ev->WireID().TPC]
            != ana::tpc2sec[geoDet][ph_trk_end->WireID().TPC]) continue;

        // check if the hit is in a track
        art::Ptr<recob::Track> pt_hit = fop_hit2trk.at(ph_ev.key());
        if (pt_hit) {
            // check if the hit is on the body of the track
            if (pt_hit.key() == p_trk.key()
                && std::find_if(
                    iph_body, vph_sec_trk.end(),
                    [&](HitPtr const& ph) -> bool {
                        return ph.key() == ph_ev.key();
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

    auto score = [&](HitPtr const& ph1, HitPtr const& ph2, double theta, double sigma) -> double {
        double dz = GetSpace(ph2->WireID()) - GetSpace(ph1->WireID());
        double dt = (ph2->PeakTime() - ph1->PeakTime()) * fTick2cm;
        double da = atan2(dt, dz) - theta;
        da = abs(da) < TMath::Pi() ? da : da - (da>0?1:-1)*2*TMath::Pi();
        double r = sqrt(dt*dt + dz*dz);

        return TMath::Gaus(da, 0, sigma) / r;
    };

    HitPtrVec vph_sec{vph_reg.begin(), vph_reg.end()};
    std::reverse(vph_sec.begin(), vph_sec.end());
    HitPtr ph_prev = ph_trk_end;
    while (vph_near.size()) {
        HitPtrVec::iterator iph_max = std::max_element(
            vph_near.begin(), vph_near.end(),
            [&](HitPtr const& ph1, HitPtr const& ph2) -> bool {
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
        // for (HitPtr const& ph : vph_reg) {
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
    HitPtr ph_max;
    for (auto iph_sec=vph_sec.begin(); iph_sec!=vph_sec.end(); iph_sec++) {
        HitPtrVec::iterator jph_sec = std::distance(iph_sec, vph_sec.end()) > trailing_radius
            ? iph_sec+trailing_radius : vph_sec.end();
        unsigned l = std::distance(iph_sec, jph_sec);

        double dQ = std::accumulate(
            iph_sec, jph_sec, 0.,
            [](double sum, HitPtr const& ph) {
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
    return ph_max;
}   



DEFINE_ART_MODULE(ana::Tagchecks)