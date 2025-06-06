////////////////////////////////////////////////////////////////////////
// Class:       Agnochecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Agnochecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"

namespace ana {
    class Agnochecks;
}


class ana::Agnochecks : public art::EDAnalyzer {
public:
    explicit Agnochecks(fhicl::ParameterSet const& p);
    Agnochecks(Agnochecks const&) = delete;
    Agnochecks(Agnochecks&&) = delete;
    Agnochecks& operator=(Agnochecks const&) = delete;
    Agnochecks& operator=(Agnochecks&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;

private:

    // Utilities
    art::ServiceHandle<art::TFileService> tfs;

    const geo::GeometryCore* asGeo;
    const geo::WireReadoutGeom* asWire;
    const detinfo::DetectorPropertiesService* asDetProp;
    const detinfo::DetectorClocksService* asDetClocks;
    geo::BoxBoundedGeo geoHighX, geoLowX;

    // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    int geoDet;
    enum EnumDet { kPDVD, kPDHD };

    // Detector Properties
    // float fADC2MeV;
    float fTick2cm; // cm/tick
    float fSamplingRate; // µs/tick
    float fDriftVelocity; // cm/µs
    float fChannelPitch;

    bounds<float> wireWindow;
    // bounds3D<float> lower_bounds, upper_bounds;
    // std::map<int,ana::bounds<unsigned>> map_tpc_ch;
    // std::map<int,float> map_ch_z;

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
    float fNearbySpaceRadius; // in cm
    float fCoincidenceWindow; // in ticks
    float fCoincidenceRadius; // in cm

    // Output Variables
    TTree* tEvent;
    bool EventIsReal;
    unsigned iEvent=0;
    unsigned EventNMuon;
    std::vector<unsigned> EventiMuon;

    ana::Hits EventHits;
    // ana::Hits EventUHits, EventVHits;


    TTree* tMuon;
    unsigned iMuon=0;

    bool MuonIsAnti;
    std::string MuonEndProcess;
    int MuonHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside };

    float MuonTrackLength;
    int MuonTrackIsNotBroken;
    enum EnumIsBroken { kBroken, kLastOfBroken, kNotBroken };
    // ana::Points MuonTrackPoints;
    ana::Point MuonEndTrackPoint;
    // ana::Point MuonTrueEndPoint;    
    // float MuonTrueEndPointT;
    // ana::Point MuonTrueEndMomentum;
    float MuonTrueEndEnergy;
    // ana::Points MuonSpacePoints;
    // ana::Point MuonEndSpacePoint;

    ana::Hits MuonHits;
    // ana::Hits MuonUHits, MuonVHits;
    ana::Hit MuonEndHit;
    // ana::Hit MuonEndUHit, MuonEndVHit;
    ana::Hit MuonTrueEndHit;
    bool MuonEndIsInWindowT, MuonEndIsInVolumeYZ;
    bool MuonEndHasGood3DAssociation;

    ana::Hits NearbyHits;
    // ana::Hits NearbyUHits, NearbyVHits;
    // ana::Points NearbySpacePoints;
    // ana::Points NearbyHitSpacePoints;
    // ana::Points NearbyHitPoints;
    // std::vector<float> NearbyHitPointQuality;
    // ana::Points NearbyReco3DPoints;

    float MichelTrackLength;

    ana::Hits MichelHits, SphereHits;

    float MichelTrueEnergy, MichelHitEnergy, SphereEnergy, TrueSphereEnergy;
    // unsigned SphereTruePositive, SphereFalsePositive;
    // float SphereEnergyTruePositive, SphereEnergyFalsePositive;


    void resetEvent();
    void resetMuon();

    bool IsUpright(recob::Track const& T);
    double GetSpace(geo::WireID);
    ana::Hit GetHit(art::Ptr<recob::Hit> const p_hit);
    art::Ptr<recob::Hit> GetDeepestHit(
        std::vector<art::Ptr<recob::Hit>> const&,
        bool increazing_z,
        geo::View_t = geo::kW
    );
    std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> GetEndsHits(
        std::vector<art::Ptr<recob::Hit>> const&,
        geo::View_t = geo::kW
    );
};


ana::Agnochecks::Agnochecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fLog(p.get<bool>("Log", false)),
    fKeepOutside(p.get<bool>("KeepOutside", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 40.F)), // in cm
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius", 20.F)), //in cm
    fNearbySpaceRadius(p.get<float>("NearbySpaceRadius", 40.F)), //in cm
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
        const std::string   process     = prod[0],
                            label       = prod[1],
                            instance    = prod[2],
                            type        = prod[3];

        const art::InputTag tag = art::InputTag{label,instance};

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
        else if (type == "recob::PFParticle")       tag_pfp = tag;
        tag_r3d = art::InputTag{"reco3d", ""};
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

    fChannelPitch = geo::WireGeo::WirePitch(
        asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 0}),
        asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 1})
    );
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

    for (unsigned t=0; t<asGeo->NTPC(); t++) {
        for (unsigned p=0; p<asWire->Nplanes(); p++) {
            geo::PlaneID pid{0, t, p};
            geo::WireGeo w0 = asWire->Wire(geo::WireID{pid, 0});
            geo::WireGeo w1 = asWire->Wire(geo::WireID{pid, 1});

            int dy = w1.GetCenter().Y() > w0.GetCenter().Y() ? 1 : -1;
            int dz = w1.GetCenter().Z() > w0.GetCenter().Z() ? 1 : -1;

            plane2axis[pid] = { dy * w0.CosThetaZ(), dz * w0.SinThetaZ() };
            plane2pitch[pid] = geo::WireGeo::WirePitch(w0, w1);
        }
    }

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector: " << std::vector<std::string>{"PDVD", "PDHD"}[geoDet]
        << "  (" << asGeo->DetectorName() << ")" << std::endl
        << "  Sampling Rate: " << fSamplingRate << " µs/tick" << std::endl
        << "  Drift Velocity: " << fDriftVelocity << " cm/µs" << std::endl
        << "  Channel Pitch: " << fChannelPitch << " cm" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX.Min() << " -> " << geoHighX.Max() << std::endl
        << "  LowX Bounds: " << geoLowX.Min() << " -> " << geoLowX.Max() << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl
        << "  Michel Space Radius: " << fMichelSpaceRadius << " cm"
        << " (" << fMichelTickRadius << " ticks)" << std::endl
        << "  Nearby Space Radius: " << fNearbySpaceRadius << " cm" << std::endl
        << "  Coincidence Window: " << fCoincidenceWindow << " ticks" << std::endl;

    tEvent = tfs->make<TTree>("event","");

    tEvent->Branch("isReal", &EventIsReal);
    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    EventHits.SetBranches(tEvent);
    // EventUHits.SetBranches(tEvent, "U");
    // EventVHits.SetBranches(tEvent, "V");

    tMuon = tfs->make<TTree>("muon","");

    tMuon->Branch("iEvent", &iEvent);
    tMuon->Branch("iMuon", &iMuon);
    tMuon->Branch("iMuonInEvent", &EventNMuon);

    tMuon->Branch("IsAnti", &MuonIsAnti);
    tMuon->Branch("EndProcess", &MuonEndProcess);
    tMuon->Branch("HasMichel", &MuonHasMichel);

    tMuon->Branch("TrackLength", &MuonTrackLength);
    tMuon->Branch("TrackIsNotBroken", &MuonTrackIsNotBroken);
    tMuon->Branch("EndIsInWindowT", &MuonEndIsInWindowT);
    tMuon->Branch("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ);
    tMuon->Branch("EndHasGood3DAssociation", &MuonEndHasGood3DAssociation);

    MuonHits.SetBranches(tMuon);
    // MuonUHits.SetBranches(tMuon, "U");
    // MuonVHits.SetBranches(tMuon, "V");
    MuonEndHit.SetBranches(tMuon, "End");
    // MuonEndUHit.SetBranches(tMuon, "EndU");
    // MuonEndVHit.SetBranches(tMuon, "EndV");
    MuonTrueEndHit.SetBranches(tMuon, "TrueEnd");
    // MuonTrackPoints.SetBranches(tMuon, "Track");
    MuonEndTrackPoint.SetBranches(tMuon, "EndTrack");
    // MuonTrueEndPoint.SetBranches(tMuon, "TrueEnd");
    // tMuon->Branch("TrueEndPointT", &MuonTrueEndPointT); // ns
    // tMuon->Branch("TrueEndMomentumX", &MuonTrueEndMomentum.x); // GeV
    // tMuon->Branch("TrueEndMomentumY", &MuonTrueEndMomentum.y); // GeV
    // tMuon->Branch("TrueEndMomentumZ", &MuonTrueEndMomentum.z); // GeV
    // tMuon->Branch("TrueEndEnergy", &MuonTrueEndEnergy); // GeV
    // MuonSpacePoints.SetBranches(tMuon, "Space");
    // MuonEndSpacePoint.SetBranches(tMuon, "EndSpace");

    tMuon->Branch("MichelTrackLength", &MichelTrackLength); // cm
    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy); // MeV
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy); // MeV
    tMuon->Branch("SphereEnergy", &SphereEnergy); // MeV
    tMuon->Branch("TrueSphereEnergy", &TrueSphereEnergy); // MeV

    // tMuon->Branch("SphereTruePositive", &SphereTruePositive);
    // tMuon->Branch("SphereFalsePositive", &SphereFalsePositive);
    // tMuon->Branch("SphereEnergyTruePositive", &SphereEnergyTruePositive);
    // tMuon->Branch("SphereEnergyFalsePositive", &SphereEnergyFalsePositive);

    MichelHits.SetBranches(tMuon, "Michel");
    SphereHits.SetBranches(tMuon, "Sphere");
    NearbyHits.SetBranches(tMuon, "Nearby");
    // NearbyUHits.SetBranches(tMuon, "NearbyU");
    // NearbyVHits.SetBranches(tMuon, "NearbyV");
    // NearbySpacePoints.SetBranches(tMuon, "NearbySpace");
    // NearbyHitSpacePoints.SetBranches(tMuon, "NearbyHitSpace");
    // NearbyHitPoints.SetBranches(tMuon, "NearbyHit");
    // tMuon->Branch("NearbyHitPointQuality", &NearbyHitPointQuality);
    // NearbyReco3DPoints.SetBranches(tMuon, "NearbyReco3D");

}

void ana::Agnochecks::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fTick2cm = fDriftVelocity * fSamplingRate;

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) return;
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) return;
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    auto const & vh_spt = e.getHandle<std::vector<recob::SpacePoint>>(tag_spt);
    if (!vh_spt.isValid()) return;
    std::vector<art::Ptr<recob::SpacePoint>> vp_spt;
    art::fill_ptr_vector(vp_spt, vh_spt);

    auto const & vh_pfp = e.getHandle<std::vector<recob::PFParticle>>(tag_pfp);
    if (!vh_pfp.isValid()) return;

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    art::FindOneP<recob::SpacePoint> fop_hit2spt(vh_hit, e, tag_spt);
    art::FindOneP<recob::Hit> fop_spt2hit(vh_spt, e, tag_spt);

    art::FindOneP<recob::PFParticle> fop_trk2pfp(vh_trk, e, tag_trk);
    art::FindManyP<recob::SpacePoint> fmp_pfp2spt(vh_pfp, e, tag_pfp);

    auto const & vh_r3d = e.getValidHandle<std::vector<recob::SpacePoint>>(tag_r3d);
    std::vector<art::Ptr<recob::SpacePoint>> vp_r3d;
    art::fill_ptr_vector(vp_r3d, vh_r3d);


    resetEvent();

    EventIsReal = e.isRealData();

    // for (recob::Hit const& hit : *vh_hit) {
    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        switch (p_hit->View()) {
            // case geo::kU: EventUHits.push_back(GetHit(p_hit)); break;
            // case geo::kV: EventVHits.push_back(GetHit(p_hit)); break;
            case geo::kW: EventHits.push_back(GetHit(p_hit)); break;
            default: break;
        }
    }

    // std::unordered_map<int, unsigned> particle_encounter;

    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        if (fLog) std::cout << "e" << iEvent << "t" << p_trk->ID() << "\r" << std::flush;

        // no short tracks
        if (!LOG(p_trk->Length() > fTrackLengthCut)) continue;

        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);

        // tracks associated to a MCTruth muon
        // if (!mcp) continue;
        // if (!LOG(abs(mcp->PdgCode()) == 13)) continue;

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());

        if (!LOG(vp_hit_muon.size())) continue;

        resetMuon();

        bool increasing_z;
        // track end point is the deepest
        if (IsUpright(*p_trk)) {
            MuonEndTrackPoint = ana::Point{p_trk->End()};
            increasing_z = p_trk->End().Z() > p_trk->Start().Z();
        } else {
            MuonEndTrackPoint = ana::Point{p_trk->Start()};
            increasing_z = p_trk->Start().Z() > p_trk->End().Z();
        }

        if (mcp) {
            art::Ptr<recob::Hit> deephit = GetDeepestHit(
                ana::mcp2hits(mcp, vp_hit, clockData, false),
                mcp->EndZ() > mcp->Vz()
            );
            if (deephit) MuonTrueEndHit = GetHit(deephit);
        } else MuonTrueEndHit = ana::Hit{};

        art::Ptr<recob::Hit> deephit = GetDeepestHit(vp_hit_muon, increasing_z);
        if (!LOG(deephit)) continue;
        MuonEndHit = GetHit(deephit);

        // !! GetDeepestHit not ready for View != W
        // deephit = GetDeepestHit(vp_hit_muon, increasing_z, geo::kU);
        // if (!LOG(deephit)) continue;
        // MuonEndUHit = GetHit(deephit);
        // deephit = GetDeepestHit(vp_hit_muon, increasing_z, geo::kV);
        // if (!LOG(deephit)) continue;
        // MuonEndVHit = GetHit(deephit);
            
        // fiducial cuts
        MuonEndIsInWindowT = wireWindow.isInside(MuonEndHit.tick, fMichelTickRadius);
        MuonEndIsInVolumeYZ = geoHighX.InFiducialY(MuonEndTrackPoint.y, fMichelSpaceRadius)
            && geoHighX.InFiducialZ(MuonEndTrackPoint.z, fMichelSpaceRadius);

        if (!LOG(fKeepOutside or (MuonEndIsInWindowT and MuonEndIsInVolumeYZ))) continue;

        // we found a muon candidate!
        if (fLog) std::cout << "\t" "\033[1;93m" "e" << iEvent << "m" << EventNMuon << " (" << iMuon << ")" "\033[0m" << std::endl;

        EventiMuon.push_back(iMuon);

        MuonIsAnti = mcp ? mcp->PdgCode() < 0 : false;
        MuonTrackLength = p_trk->Length();

        if (mcp) {
            std::vector<art::Ptr<recob::Track>> vp_trk_from_mcp = mcp2trks(mcp, vp_trk, clockData, fmp_trk2hit);
            // std::cout << "trk#" << p_trk->ID() << " mu#" << mcp->TrackId() << " encounter#" << ++particle_encounter[mcp->TrackId()] << " #mcp2trk:" << vp_trk_from_mcp.size();
            if (vp_trk_from_mcp.size()) {
                // bool isin = std::find(vp_trk_from_mcp.begin(), vp_trk_from_mcp.end(), p_trk) != vp_trk_from_mcp.end();
                // std::cout << " \033[1;9" << (isin ? 2 : 1) << "m" << (isin ? "in" : "not in") << "\033[0m";
                // std::cout << " [ ";
                // for (art::Ptr<recob::Track> p : vp_trk_from_mcp) std::cout << p->ID() << ", ";
                // std::cout << "]" << std::endl;
                if (vp_trk_from_mcp.size() == 1)
                    MuonTrackIsNotBroken = kNotBroken;
                else {
                    bool IsDeepestTrack = true;
                    for (art::Ptr<recob::Track> p_trk_from_mcp : vp_trk_from_mcp)
                        if (geoDet == kPDVD)
                            IsDeepestTrack = IsDeepestTrack
                                && (MuonEndTrackPoint.x <= p_trk_from_mcp->Start().X()
                                && MuonEndTrackPoint.x <= p_trk_from_mcp->End().X());
                        else if (geoDet == kPDHD)
                            IsDeepestTrack = IsDeepestTrack
                                && (MuonEndTrackPoint.y <= p_trk_from_mcp->Start().Y()
                                && MuonEndTrackPoint.y <= p_trk_from_mcp->End().Y());
                    if (IsDeepestTrack)
                        MuonTrackIsNotBroken = kLastOfBroken;
                    else
                        MuonTrackIsNotBroken = kBroken;
                }
            } else MuonTrackIsNotBroken = -1;
        } else MuonTrackIsNotBroken = -1;

        if (mcp) {
            MuonEndProcess = mcp->EndProcess();
            // MuonTrueEndPoint = ana::Point{mcp->EndPosition().Vect()};
            // MuonTrueEndPointT = mcp->EndT();
            // MuonTrueEndMomentum = ana::Point{mcp->EndMomentum().Vect()};
            // MuonTrueEndEnergy = mcp->EndE();
        } else {
            MuonEndProcess = "";
            // MuonTrueEndPoint = ana::Point{};
            // MuonTrueEndPointT = 0;
            // MuonTrueEndMomentum = ana::Point{};
            // MuonTrueEndEnergy = 0;
        }

        // getting all muon hits
        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            switch (p_hit_muon->View()) {
                // case geo::kU: MuonUHits.push_back(GetHit(p_hit_muon)); break;
                // case geo::kV: MuonVHits.push_back(GetHit(p_hit_muon)); break;
                case geo::kW: MuonHits.push_back(GetHit(p_hit_muon)); break;
                default: break;
            }
        }

        // and all muon track points
        // for (unsigned i_tpt=0; i_tpt<p_trk->NumberTrajectoryPoints(); i_tpt++) {
        //     if (!p_trk->HasValidPoint(i_tpt)) continue;
        //     MuonTrackPoints.push_back(p_trk->LocationAtPoint(i_tpt));
        // }

        // and all muon space points
        // MuonEndSpacePoint = ana::Point{geoUp.MaxX(), 0., 0.};
        // art::Ptr<recob::PFParticle> p_pfp = fop_trk2pfp.at(p_trk.key());
        // std::vector<art::Ptr<recob::SpacePoint>> v_spt_muon = fmp_pfp2spt.at(p_pfp.key());
        // for (art::Ptr<recob::SpacePoint> const& p_spt : v_spt_muon) {
        //     MuonSpacePoints.push_back(p_spt->position());

        //     if (p_spt->position().x() < MuonEndSpacePoint.x)
        //         MuonEndSpacePoint = ana::Point{p_spt->position()};
        // }
        // MuonEndHasGood3DAssociation = MuonEndSpacePoint.x != geoUp.MaxX()
        //     && abs(MuonEndHit.z - MuonEndSpacePoint.z) < fCoincidenceRadius;

        // a decaying muon has nu_mu, nu_e and elec as last daughters
        simb::MCParticle const* mcp_michel = nullptr;
        if (mcp) {
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
                bool isin = false;
                for (unsigned t=0; t<asGeo->NTPC(); t++) {
                    isin = asGeo->TPC(geo::TPCID{0, t}).ContainsPosition(mcp_michel->Position().Vect());
                    if (isin) break;
                }
                if (isin)
                    MuonHasMichel = kHasMichelInside; 
                else
                    MuonHasMichel = kHasMichelOutside;

                // recob::Track const * trk_michel = truthUtil.GetRecoTrackFromMCParticle(clockData, *mcp_michel, e, tag_trk.label());
                art::Ptr<recob::Track> trk_michel = ana::mcp2trk(mcp_michel, vp_trk, clockData, fmp_trk2hit);
                if (trk_michel) 
                    MichelTrackLength = trk_michel->Length();
                else
                    MichelTrackLength = 0;

                MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;

                // std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());
                std::vector<art::Ptr<recob::Hit>> vp_hit_michel = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);
                for (art::Ptr<recob::Hit> p_hit_michel : vp_hit_michel) {
                    if (p_hit_michel->View() != geo::kW) continue;

                    MichelHits.push_back(GetHit(p_hit_michel));
                }
                MichelHitEnergy = MichelHits.energy();
            }
            else MuonHasMichel = kNoMichel;
        } else {
            MuonHasMichel = -1;
            MichelTrackLength = -1;
            MichelHitEnergy = -1;
        }


        // get induction hits nearby muon end point
        // for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        //     if (p_hit->View() != geo::kU && p_hit->View() != geo::kV) continue;

        //     art::Ptr<recob::Track> p_hit_trk = fop_hit2trk.at(p_hit.key());
        //     bool from_track = p_hit_trk and p_hit_trk->Length() > fTrackLengthCut;

        //     if (from_track and (p_trk.key() != p_hit_trk.key())) continue;

        //     float pitch = plane2pitch[p_hit->WireID()];
        //     if (p_hit->View() == geo::kU) {
        //         float dz = (p_hit->Channel() - MuonEndUHit.channel) * pitch;
        //         float dt = (p_hit->PeakTime() - MuonEndUHit.tick) * fTick2cm;
        //         float dr2 = dz*dz + dt*dt;

        //         if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

        //         NearbyUHits.push_back(GetHit(p_hit));
        //     }
        //     if (p_hit->View() == geo::kV) {
        //         float dz = (p_hit->Channel() - MuonEndVHit.channel) * pitch;
        //         float dt = (p_hit->PeakTime() - MuonEndVHit.tick) * fTick2cm;
        //         float dr2 = dz*dz + dt*dt;

        //         if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

        //         NearbyVHits.push_back(GetHit(p_hit));
        //     }
        // }


        // get all hits nearby muon end point
        Hits TrueSphereHits;
        std::vector<art::Ptr<recob::Hit>> NearbyPHits;
        for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;

            art::Ptr<recob::Track> p_hit_trk = fop_hit2trk.at(p_hit.key());
            bool from_track = p_hit_trk and p_hit_trk->Length() > fTrackLengthCut;

            if (from_track and (p_trk.key() != p_hit_trk.key())) continue;

            ana::Hit hit = GetHit(p_hit);
            if (geoDet == kPDVD) {
                if (hit.slice() != MuonEndHit.slice()) continue;
            } else {
                if (hit.tpc != MuonEndHit.tpc) continue;
            }


            float dz = (hit.space - MuonEndHit.space);
            float dt = (hit.tick - MuonEndHit.tick) * fTick2cm;
            float dr2 = dz*dz + dt*dt;

            if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

            NearbyHits.push_back(hit);
            NearbyPHits.push_back(p_hit);

            // art::Ptr<recob::SpacePoint> p_hit_spt = fop_hit2spt.at(p_hit.key());
            // if (p_hit_spt)
            //     NearbyHitSpacePoints.push_back(p_hit_spt->position());
            // else 
            //     NearbyHitSpacePoints.push_back(ana::Point{});

            if (!mcp_michel) continue;
            if (from_track) continue;

            float dzt = (hit.space - MuonTrueEndHit.space);
            float dtt = (hit.tick - MuonTrueEndHit.tick) * fTick2cm;
            float dr2t = dzt*dzt + dtt*dtt;
            if (dr2t <= fMichelSpaceRadius * fMichelSpaceRadius)
                TrueSphereHits.push_back(hit);

            if (dr2 > fMichelSpaceRadius * fMichelSpaceRadius) continue;

            SphereHits.push_back(hit);

            // checking if the hit is associated to the michel MCParticle
            // std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());
            // if (std::find(v_hit_michel.begin(), v_hit_michel.end(), &*p_hit) != v_hit_michel.end()) {
            // if (MichelHits.find(*p_hit) != MichelHits.N) {
            //     SphereTruePositive++;
            //     SphereEnergyTruePositive += hit.adc;
            // } else {
            //     SphereFalsePositive++;
            //     SphereEnergyFalsePositive += hit.adc;
            // }
        } // end of loop over event hits
        SphereEnergy = SphereHits.energy();
        TrueSphereEnergy = TrueSphereHits.energy();

        /*

        // get space points nearby muon end point
        for (art::Ptr<recob::SpacePoint> const& p_spt : vp_spt) {

            art::Ptr<recob::Hit> p_spt_hit = fop_spt2hit.at(p_spt.key());
            art::Ptr<recob::Track> p_spt_trk = fop_hit2trk.at(p_spt_hit.key());

            if (p_spt_trk and p_spt_trk.key() != p_trk.key() and p_spt_trk->Length() > fTrackLengthCut) continue;

            ana::Point dpt = MuonEndSpacePoint - p_spt->position();
            dpt.x = 0;
            if (dpt.r2() > fNearbySpaceRadius * fNearbySpaceRadius) continue;

            NearbySpacePoints.push_back(p_spt->position());
        }
        
        // get reco3d space points nearby muon end point
        for (art::Ptr<recob::SpacePoint> const& p_r3d : vp_r3d) {
            ana::Point dpt = MuonEndSpacePoint - p_r3d->position();
            dpt.x = 0;
            if (dpt.r2() > fNearbySpaceRadius * fNearbySpaceRadius) continue;

            NearbyReco3DPoints.push_back(p_r3d->position());
        }

        */

        /* RECREATING SPACE POINTS FROM NEARBY HITS ATTEMPT

        // ana::Points NearbySpaceHits;
        std::cout << "mu#" << iMuon << " w/ " << NearbyPHits.size() << " nearby hits, w/ " << NearbySpacePoints.size() << " nearby space points" << std::endl;
        for (art::Ptr<recob::Hit> const& p_hit_col : NearbyPHits) {

            art::Ptr<recob::SpacePoint> p_spt = fop_hit2spt.at(p_hit_col.key());
            std::cout << "\033[1m" "hit w/ " << (p_spt ? "1" : "0") << " spt: " "\033[0m";
            if (p_spt) std::cout << " (" << p_spt->position().x() << ", " << p_spt->position().y() << ", " << p_spt->position().z() << ")";
            std::cout << std::endl;

            // assuming MuonEndHasGood3DAssociation is true
            geo::WireGeo const wiregeo_col = asWire->Wire(p_hit_col->WireID());
            float y = 0;
            float z = wiregeo_col.GetCenter().Z();
            float x;
            if (geoUp.ContainsX(MuonEndSpacePoint.x))
                x = MuonEndSpacePoint.x - (p_hit_col->PeakTime() - MuonEndHit.tick) * fTick2cm;
            else
                x = MuonEndSpacePoint.x + (p_hit_col->PeakTime() - MuonEndHit.tick) * fTick2cm;

            struct Coincidence {
                float y;
                recob::Hit const* hit;
            };
            std::vector<struct Coincidence> V_coincidences, U_coincidences;

            std::cout << "  " << (MuonEndHasGood3DAssociation ? "good3D" : "bad3D") << " EndSpt: " << MuonEndSpacePoint << "  dtick: " << p_hit_col->PeakTime() - MuonEndHit.tick << std::endl;
            geo::TPCGeo const tpcgeo = asGeo->TPC(geo::TPCID{0, p_hit_col->WireID().TPC});
            std::cout << "  TPC " << p_hit_col->WireID().TPC << ": " << tpcgeo.Min() << " -> " << tpcgeo.Max() << std::endl;

            for (recob::Hit const& hit_ind : *vh_hit) {

                if (hit_ind.WireID().TPC != p_hit_col->WireID().TPC) continue;
                if (hit_ind.View() != geo::kU and hit_ind.View() != geo::kV) continue;
                if (abs(hit_ind.PeakTime() - p_hit_col->PeakTime()) > fCoincidenceWindow) continue;

                geo::WireGeo const wiregeo_ind = asWire->Wire(hit_ind.WireID());
                float co_y = geo::WiresIntersection(wiregeo_col, wiregeo_ind).Y();

                if (!tpcgeo.ContainsYZ(co_y,z)) continue;
                struct Coincidence co = {co_y, &hit_ind};

                switch (hit_ind.View()) {
                    case geo::kU: U_coincidences.push_back(co); break;
                    case geo::kV: V_coincidences.push_back(co); break;
                    default: continue;
                }
            }
            std::cout << "  U coincidences: " << U_coincidences.size() << ", V coincidences: " << V_coincidences.size() << std::endl;
            if (U_coincidences.empty() or V_coincidences.empty()) continue;

            float min_dy = 600.F;
            std::vector<float> barys_y;
            bool has_good_coincidence = false;
            for (struct Coincidence const& U_co : U_coincidences) {
                for (struct Coincidence const& V_co : V_coincidences) {

                    float dy = abs(U_co.y - V_co.y);

                    // std::cout << "  Upt: " << U_co.pt << " w/ " << U_co.hit->Integral() << " Vpt: " << V_co.pt << " w/ " << V_co.hit->Integral() << " bary: " << bary << " w/ dy: " << dy << std::endl;
                    float bary_y = (U_co.y * U_co.hit->Integral() + V_co.y * V_co.hit->Integral()) / (U_co.hit->Integral() + V_co.hit->Integral());

                    if (has_good_coincidence) {
                        if (dy < fCoincidenceRadius) {
                            barys_y.push_back(bary_y);
                        }
                    } else if (dy < fCoincidenceRadius) {
                        has_good_coincidence = true;
                        barys_y.clear();
                        barys_y.push_back(bary_y);
                    } else if (dy < min_dy) {
                        min_dy = dy;
                        barys_y.clear();
                        barys_y.push_back(bary_y);
                    }

                    // if ((bary - MuonEndSpacePoint).r2() > fNearbySpaceRadius * fNearbySpaceRadius) continue;

                    // std::cout << "  " << bary << std::endl;
                }
            }
            // std::cout << "  best bary: " << best_bary << " w/ dy: " << min_dy << std::endl;
            float mean_y = std::accumulate(barys_y.begin(), barys_y.end(), 0.F) / barys_y.size();
            unsigned n = 0;
            for (float bary_y : barys_y) {
                if (abs(bary_y - mean_y) < 2 * fCoincidenceRadius) {
                    y += bary_y;
                    n++;
                }
            }
            if (n) y /= n;
            else {
                y = mean_y;
                has_good_coincidence = false;
            }

            ana::Point pt{x, y, z};
            std::cout << "  custom spt: " << pt << " oof " << barys_y.size() << " w/ " << (has_good_coincidence ? 0.F : min_dy) << " point quality" << std::endl;
            NearbyHitPoints.push_back(pt);
            NearbyHitPointQuality.push_back(has_good_coincidence ? 0.F : min_dy);
        }

        */


        tMuon->Fill();
        iMuon++;
        EventNMuon++;
    } // end of loop over tracks

    tEvent->Fill();
    iEvent++;
}

void ana::Agnochecks::beginJob() {}
void ana::Agnochecks::endJob() {}


void ana::Agnochecks::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
    // EventUHits.clear();
    // EventVHits.clear();
}
void ana::Agnochecks::resetMuon() {
    // MuonTrackPoints.clear();
    // MuonSpacePoints.clear();
    MuonHits.clear();
    // MuonUHits.clear();
    // MuonVHits.clear();

    NearbyHits.clear();
    // NearbyUHits.clear();
    // NearbyVHits.clear();
    // NearbySpacePoints.clear();
    // NearbyHitSpacePoints.clear();
    // NearbyHitPoints.clear();
    // NearbyHitPointQuality.clear();
    // NearbyReco3DPoints.clear();

    MichelTrackLength = 0;

    MichelTrueEnergy = 0;

    MichelHits.clear();
    MichelHitEnergy = 0;

    SphereHits.clear();
    SphereEnergy = 0;
    TrueSphereEnergy = 0;

    // SphereTruePositive = 0;
    // SphereFalsePositive = 0;
    // SphereEnergyTruePositive = 0;
    // SphereEnergyFalsePositive = 0;
}


double ana::Agnochecks::GetSpace(geo::WireID wid) {
    return plane2axis[(geo::PlaneID) wid].space(asWire->Wire(wid));
}

ana::Hit ana::Agnochecks::GetHit(art::Ptr<recob::Hit> const p_hit) {
    geo::WireID wid = p_hit->WireID();
    // if (geoDet == kPDHD)
    //     for (int t : (int[]){0, 4, 3, 7})
    //         if (wireid.TPC == t)
    //             return ana::Hit{};
    
    return ana::Hit{
        wid.TPC,
        float(GetSpace(wid)),
        p_hit->Channel(),
        p_hit->PeakTime(),
        p_hit->Integral()
    };
}

bool ana::Agnochecks::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}

art::Ptr<recob::Hit> ana::Agnochecks::GetDeepestHit(
    std::vector<art::Ptr<recob::Hit>> const& vp_hit,
    bool increasing_z,
    geo::View_t view
) {
    if (vp_hit.empty()) return art::Ptr<recob::Hit>{};

    art::Ptr<recob::Hit> DeepestHit;
    if (geoDet == kPDVD) {

        // test if the muon goes into the bottom volume
        bool in_bot = false;
        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            unsigned tpc = p_hit->WireID().TPC;
            if (tpc < 8) {
                in_bot = true;
                break;
            }
        }

        // basic linear regression on the hits that are in the last volume
        unsigned n=0;
        double mz=0, mt=0, mt2=0, mzt=0;
        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            if (in_bot && p_hit->WireID().TPC >= 8) continue;
            // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
            double z = GetSpace(p_hit->WireID());
            double t = p_hit->PeakTime() * fTick2cm;
            mz += z; mt += t; mt2 += t*t, mzt += z*t;
            n++;
        }
        mz /= n; mt /= n; mt2 /= n; mzt /= n;
        double cov = mzt - mz*mt;
        double vart = mt2 - mt*mt;

        // z ~ m*t + p
        double m = cov / vart;
        double p = mz - m*mt;

        double extrem_s = in_bot ?
            std::numeric_limits<double>::max() // search min_s
            : std::numeric_limits<double>::lowest(); // search max_s

        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            if (in_bot && p_hit->WireID().TPC >= 8) continue;
            // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
            double z = GetSpace(p_hit->WireID());
            double t = p_hit->PeakTime() * fTick2cm;

            // projection on the axis of the track
            double s = (t + m*(z-p)) / (1 + m*m);
            if (in_bot) {
                if (extrem_s > s) {
                    extrem_s = s;
                    DeepestHit = p_hit;
                }
            } else {
                if (extrem_s < s) {
                    extrem_s = s;
                    DeepestHit = p_hit;
                }
            }
        }
    } else if (geoDet == kPDHD) {

        // test if the muon crosses the cathod
        unsigned n_left=0, n_right=0;
        double mz_left=0, mz_right=0;
        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            unsigned tpc = p_hit->WireID().TPC;
            // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
            double z = GetSpace(p_hit->WireID());
            if (tpc == 1 || tpc == 5) {
                mz_left += z;
                n_left++;
            } else if (tpc == 2 || tpc == 6) {
                mz_right += z;
                n_right++;
            }
        }
        mz_left /= n_left; mz_right /= n_right;
        std::pair<unsigned, unsigned> tpcs;
        if (
            (increasing_z && mz_left > mz_right)
            || (!increasing_z && mz_left < mz_right)
        ) {
            tpcs.first = 1;
            tpcs.second = 5;
        } else {
            tpcs.first = 2;
            tpcs.second = 6;
        }

        // basic linear regression on the hits that are in the last volume
        unsigned n=0;
        double mz=0, mt=0, mz2=0, mzt=0;
        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            unsigned tpc = p_hit->WireID().TPC;
            if (tpc != tpcs.first && tpc != tpcs.second) continue;
            // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
            double z = GetSpace(p_hit->WireID());
            double t = p_hit->PeakTime() * fTick2cm;
            mz += z; mt += t; mz2 += z*z; mzt += z*t;
            n++;
        }
        mz /= n; mt /= n; mz2 /= n; mzt /= n;
        double cov = mzt - mz*mt;
        double varz = mz2 - mz*mz;

        // t ~ m*z + p
        double m = cov / varz;
        double p = mt - m*mz;

        double extrem_s = increasing_z ?
            std::numeric_limits<double>::lowest() // search max_s
            : std::numeric_limits<double>::max(); // search min_s

        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            unsigned tpc = p_hit->WireID().TPC;
            if (tpc != tpcs.first && tpc != tpcs.second) continue;
            // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
            double z = GetSpace(p_hit->WireID());
            double t = p_hit->PeakTime() * fTick2cm;

            // projection on the axis of the track
            double s = (z + m*(t-p)) / (1 + m*m);
            if (increasing_z) {
                if (extrem_s < s) {
                    extrem_s = s;
                    DeepestHit = p_hit;
                }
            } else {
                if (extrem_s > s) {
                    extrem_s = s;
                    DeepestHit = p_hit;
                }
            }
        }
    }
    return DeepestHit;
}


std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> ana::Agnochecks::GetEndsHits(
    std::vector<art::Ptr<recob::Hit>> const& vp_hit,
    geo::View_t view
) {
    using HitPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;
    if (vp_hit.empty()) return HitPair{};

    auto cathodeSide =
        geoDet == kPDVD ?
        [](geo::TPCID::TPCID_t tpc) -> int { return int(tpc >= 8); }
        : [](geo::TPCID::TPCID_t tpc) -> int { return tpc == 1 || tpc == 5 ? 1 : (tpc == 2 || tpc == 6 ? 0 : -1); };

    struct LinearRegression {
        unsigned n=0;
        double mz=0, mt=0, mt2=0, mzt=0;
        void add(double z, double t) {
            mz+=z; mt+=t; mt2+=t*t; mzt+=z*t;
        }
        void normalize() {
            mz/=n; mt/=n; mt2/=n; mzt/=n;
        }
        double cov() const { return mzt - mz*mt; }
        double vart() const { return mt2 - mt*mt; }
        double m() const { return cov() / vart(); }
        double p() const { return mz - m()*mt; }
        double projection(double z, double t) const {
            return (t + m()*(z-p())) / (1 + m()*m());
        }
    } reg_side0, reg_side1;

    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        geo::WireID wireid = p_hit->WireID();
        int side = cathodeSide(wireid.TPC);
        if (side == -1) continue;
        // double z = asWire->Wire(wireid).GetCenter().Z();
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        if (side == 0) {
            reg_side0.add(z, t);
        } else if (side == 1) {
            reg_side1.add(z, t);
        }
    }
    reg_side0.normalize();
    reg_side1.normalize();

    struct ProjectionEnds {
        HitPair hits;
        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::lowest();
        void test(double s, art::Ptr<recob::Hit> const& h) {
            if (s < min) {
                min=s;
                hits.first = h;
            } else if (s > max) {
                max=s;
                hits.second = h;
            }
        }
    } ends_side0, ends_side1;

    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        geo::WireID wireid = p_hit->WireID();
        int side = cathodeSide(wireid.TPC);
        if (side == -1) continue;
        // double z = asWire->Wire(wireid).GetCenter().Z();
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        if (side == 0) {
            double s = reg_side0.projection(z, t);
            ends_side0.test(s, p_hit);
        } else if (side == 1) {
            double s = reg_side1.projection(z, t);
            ends_side1.test(s, p_hit);
        }
    }

    auto d2 = [this](art::Ptr<recob::Hit> h1, art::Ptr<recob::Hit> h2) {
        // double z1 = asWire->Wire(h1->WireID()).GetCenter().Z();
        double z1 = GetSpace(h1->WireID());
        double t1 = h1->PeakTime() * fTick2cm;
        // double z2 = asWire->Wire(h2->WireID()).GetCenter().Z();
        double z2 = GetSpace(h2->WireID());
        double t2 = h2->PeakTime() * fTick2cm;
        return pow(z1-z2,2) + pow(t1-t2,2);
    };

    std::vector<double> distances(4, 0);
    distances[0] = d2(ends_side0.hits.first, ends_side1.hits.first);
    distances[1] = d2(ends_side0.hits.first, ends_side1.hits.second);
    distances[2] = d2(ends_side0.hits.second, ends_side1.hits.first);
    distances[3] = d2(ends_side0.hits.second, ends_side1.hits.second);

    switch(std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()))) {
        case 0: return HitPair{ends_side0.hits.second, ends_side1.hits.second};
        case 1: return HitPair{ends_side0.hits.second, ends_side1.hits.first};
        case 2: return HitPair{ends_side0.hits.first, ends_side1.hits.second};
        case 3: return HitPair{ends_side0.hits.first, ends_side1.hits.first};
        default: break;
    }
    return HitPair{};
}

// ana::Hit ana::Agnochecks::GetTrueEndHit( std::vector<art::Ptr<recob::Hit>> const& vp_hit, ana::Point end_z) {
//     double min_dz = std::numeric_limits<double>::max();
//     art::Ptr<recob::Hit> p_min_dz{};
//     for (art::Ptr<recob::Hit> p_hit : vp_hit) {
//         if (p_hit->View() != geo::kW) continue;
//         double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
//         double dz = abs(z - end_z);
//         if (min_dz > dz) {
//             min_dz = dz;
//             p_min_dz = p_hit;
//         }
//     }
//     std::cout << min_dz << std::endl;
//     return GetHit(p_min_dz);
// }

DEFINE_ART_MODULE(ana::Agnochecks)