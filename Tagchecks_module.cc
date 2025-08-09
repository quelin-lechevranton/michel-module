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
    float fMichelSpaceRadius; // in cm
    float fMichelTickRadius; // in ticks
    float fCoincidenceWindow; // in ticks
    float fCoincidenceRadius; // in cm

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

    bool TrueTagDownward;
    int TrueTagPdg;
    std::string TrueTagEndProcess;

    int TrueTagHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside };

    ana::Hits MuonHits;
    ana::Hit MuonEndHit;
    ana::Hit MuonTrueEndHit;
    bool TagEndIsInWindowT, TagEndIsInVolumeYZ;

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
};


ana::Tagchecks::Tagchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fLog(p.get<bool>("Log", false)),
    fKeepOutside(p.get<bool>("KeepOutside", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 40.F)), // in cm
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius", 20.F)), //in cm
    fCoincidenceWindow(p.get<float>("CoincidenceWindow", 1.F)), // in ticks
    fCoincidenceRadius(p.get<float>("CoincidenceRadius", 1.F)) // in cm
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
    fMichelTickRadius = fMichelSpaceRadius / fDriftVelocity / fSamplingRate;

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
        << "  Michel Space Radius: " << fMichelSpaceRadius << " cm"
        << " (" << fMichelTickRadius << " ticks)" << std::endl
        << "  Coincidence Window: " << fCoincidenceWindow << " ticks" << std::endl;

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

    tMuon->Branch("TagCathodeCrossing", &TagCathodeCrossing);
    tMuon->Branch("TagAnodeCrossing", &TagAnodeCrossing);
    tMuon->Branch("CutTrackLength", &CutTrackLength);
    tMuon->Branch("TagEndIsInWindowT", &TagEndIsInWindowT);
    tMuon->Branch("TagEndIsInVolumeYZ", &TagEndIsInVolumeYZ);

    tMuon->Branch("TrueTagPdg", &TrueTagPdg);
    tMuon->Branch("TrueTagDownward", &TrueTagDownward);
    tMuon->Branch("TrueTagEndProcess", &TrueTagEndProcess);
    tMuon->Branch("TrueTagHasMichel", &TrueTagHasMichel);


    MuonHits.SetBranches(tMuon);
    MuonEndHit.SetBranches(tMuon, "End");
    MuonTrueEndHit.SetBranches(tMuon, "TrueEnd");

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

        HitPtrVec vp_hit_muon = fmp_trk2hit.at(p_trk.key());
        ASSERT(vp_hit_muon.size() >= 2*ana::LinearRegression::nmin)

        CutTrackLength = p_trk->Length();
        TagCathodeCrossing = (
            geoLowX.ContainsPosition(p_trk->Start()) 
            && geoHighX.ContainsPosition(p_trk->End())
        ) || (
            geoHighX.ContainsPosition(p_trk->Start())
            && geoLowX.ContainsPosition(p_trk->End())
        );

        // Anode crossing: SUPPOSITION: Muon is downward
        TagIsUpright =  IsUpright(*p_trk);
        geo::Point_t Start = TagIsUpright ? p_trk->Start() : p_trk->End();
        geo::Point_t End = TagIsUpright ? p_trk->End() : p_trk->Start();
        if (geoDet == kPDVD)
            TagAnodeCrossing = geoHighX.ContainsYZ(Start.Y(), Start.Z(), 0.8);
        else if (geoDet == kPDHD)
            TagAnodeCrossing = (
                geoHighX.ContainsYZ(Start.Y(), Start.Z(), 0.8)
                || geoLowX.ContainsYZ(Start.Y(), Start.Z(), 0.8)
            );

        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        ASSERT(mcp)

        TrueTagPdg = mcp->PdgCode();
        TrueTagEndProcess = mcp->EndProcess();

        if (geoDet == kPDVD)
            TrueTagDownward = mcp->Position(0).X() > mcp->EndPosition().X();
        else if (geoDet == kPDHD)
            TrueTagDownward = mcp->Position(0).Y() > mcp->EndPosition().Y();

        resetMuon();

        HitPtrVec vp_hit_muon_sorted = GetSortedHits(
            vp_hit_muon, 
            End.Z() > Start.Z() ? 1 : -1
        );
        MuonEndHit = GetHit(vp_hit_muon_sorted.back());

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

        // fiducial cuts
        TagEndIsInWindowT = wireWindow.isInside(MuonEndHit.tick, fMichelTickRadius);
        TagEndIsInVolumeYZ = geoHighX.InFiducialY(End.Y(), fMichelSpaceRadius)
            && geoHighX.InFiducialZ(End.Z(), fMichelSpaceRadius);

        ASSERT(fKeepOutside or (TagEndIsInWindowT and TagEndIsInVolumeYZ))

        // we found a muon candidate!
        if (fLog) std::cout << "\t" "\033[1;93m" "e" << iEvent << "m" << EventNMuon << " (" << iMuon << ")" "\033[0m" << std::endl;

        EventiMuon.push_back(iMuon);

        // getting all muon hits
        for (HitPtr const& p_hit_muon : vp_hit_muon_sorted)
            if (p_hit_muon->View() == geo::kW)
                MuonHits.push_back(GetHit(p_hit_muon));

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
        reg.normalize();
    if (p_side_reg) {
        p_side_reg->clear();
        p_side_reg->push_back(side_reg[0]);
        p_side_reg->push_back(side_reg[1]);
    }
    for (int side=0; side<2; side++) {
        std::sort(
            side_hit[side].begin(),
            side_hit[side].end(),
            [&, &reg=side_reg[side]](
                HitPtr const& h1, HitPtr const& h2
            ) -> bool {
                double const s1 = reg.projection(
                    GetSpace(h1->WireID()),
                    h1->PeakTime() * fTick2cm
                );
                double const s2 = reg.projection(
                    GetSpace(h2->WireID()),
                    h2->PeakTime() * fTick2cm
                );
                return (s2 - s1) * dirz > 0;
            }
        );
    }
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

DEFINE_ART_MODULE(ana::Tagchecks)