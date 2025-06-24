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

using HitPtr = art::Ptr<recob::Hit>;
using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

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
    float fADCtoMeV;
    float fTick2cm; // cm/tick
    float fSamplingRate; // µs/tick
    float fDriftVelocity; // cm/µs
    float fChannelPitch; // cm/channel
    float fCathodeGap; // cm

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
    ana::Hit GetHit(HitPtr const p_hit);
    // HitPtr GetDeepestHit(
    //     HitPtrVec const&,
    //     bool increazing_z,
    //     geo::View_t = geo::kW
    // );
    HitPtrPair GetTrackEndsHits(
        HitPtrVec const& vp_hit,
        HitPtrPair *pp_cathode_crossing = nullptr,
        HitPtrVec *vp_tpc_crossing = nullptr,
        HitPtrVec *vp_sorted_hit = nullptr,
        geo::View_t view = geo::kW
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
    // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    fADCtoMeV = (geoDet == kPDVD ? 200 : 1000) * 23.6 * 1e-6 / 0.7;
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
    fCathodeGap = geoHighX.MinX() - geoLowX.MaxX();

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
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
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
    HitPtrVec vp_hit;
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
    for (HitPtr p_hit : vp_hit) {
        switch (p_hit->View()) {
            // case geo::kU: EventUHits.push_back(GetHit(p_hit)); break;
            // case geo::kV: EventVHits.push_back(GetHit(p_hit)); break;
            case geo::kW: EventHits.push_back(GetHit(p_hit)); break;
            default: break;
        }
    }

    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        if (fLog) std::cout << "e" << iEvent << "t" << p_trk->ID() << "\r" << std::flush;

        // no short tracks
        if (!LOG(p_trk->Length() > fTrackLengthCut)) continue;

        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);

        // tracks associated to a MCTruth muon
        LOG(mcp);
        // if (!mcp) continue;
        // if (!LOG(abs(mcp->PdgCode()) == 13)) continue;

        HitPtrVec vp_hit_muon = fmp_trk2hit.at(p_trk.key());

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
            // HitPtr deephit = getdeepesthit(
            //     ana::mcp2hits(mcp, vp_hit, clockData, false),
            //     mcp->EndZ() > mcp->Vz()
            // );
            
            // if (deephit) MuonTrueEndHit = GetHit(deephit);

            HitPtrPair ends;
            ends = GetTrackEndsHits(ana::mcp2hits(
                mcp, vp_hit, clockData, false
            ));

            if (ends.first && ends.second) {
                int dir_z = mcp->EndZ() > mcp->Vz() ? 1 : -1;
                float fz = GetSpace(ends.first->WireID());
                float sz = GetSpace(ends.second->WireID());
                MuonTrueEndHit = GetHit(
                    (sz-fz) * dir_z > 0 ? ends.second : ends.first
                );
            } else MuonTrueEndHit = ana::Hit{};
        } else MuonTrueEndHit = ana::Hit{};

        HitPtrPair ends;
        ends = GetTrackEndsHits(vp_hit_muon);

        if (!LOG(ends.first && ends.second)) continue;

        int dir_z = increasing_z ? 1 : -1;
        float fz = GetSpace(ends.first->WireID());
        float sz = GetSpace(ends.second->WireID());
        MuonEndHit = GetHit(
            (sz-fz) * dir_z > 0 ? ends.second : ends.first
        );

        // HitPtr deephit = GetDeepestHit(vp_hit_muon, increasing_z);
        // if (!LOG(deephit)) continue;
        // MuonEndHit = GetHit(deephit);

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

        /*
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
        */

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
        for (HitPtr const& p_hit_muon : vp_hit_muon) {
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
                HitPtrVec vp_hit_michel = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);
                for (HitPtr p_hit_michel : vp_hit_michel) {
                    if (p_hit_michel->View() != geo::kW) continue;

                    MichelHits.push_back(GetHit(p_hit_michel));
                }

                MichelHitEnergy = MichelHits.energy() * fADCtoMeV;
            } else MuonHasMichel = kNoMichel;

        } else {
            MuonHasMichel = -1;
            MichelTrackLength = -1;
            MichelHitEnergy = -1;
        }

        // get induction hits nearby muon end point
        // for (HitPtr const& p_hit : vp_hit) {
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
        // HitPtrVec NearbyPHits;
        for (HitPtr const& p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;

            art::Ptr<recob::Track> p_hit_trk = fop_hit2trk.at(p_hit.key());
            bool from_track = p_hit_trk and p_hit_trk->Length() > fTrackLengthCut;

            if (from_track and (p_trk.key() != p_hit_trk.key())) continue;

            ana::Hit hit = GetHit(p_hit);
            if (geoDet == kPDVD && hit.slice() != MuonEndHit.slice()) continue;
            if (geoDet == kPDHD && hit.tpc != MuonEndHit.tpc) continue;

            float dz = (hit.space - MuonEndHit.space);
            float dt = (hit.tick - MuonEndHit.tick) * fTick2cm;
            float dr2 = dz*dz + dt*dt;

            if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

            NearbyHits.push_back(hit);
            // NearbyPHits.push_back(p_hit);

            // art::Ptr<recob::SpacePoint> p_hit_spt = fop_hit2spt.at(p_hit.key());
            // if (p_hit_spt)
            //     NearbyHitSpacePoints.push_back(p_hit_spt->position());
            // else 
            //     NearbyHitSpacePoints.push_back(ana::Point{});

            if (from_track) continue;

            if (mcp_michel) {
                float dzt = (hit.space - MuonTrueEndHit.space);
                float dtt = (hit.tick - MuonTrueEndHit.tick) * fTick2cm;
                float dr2t = dzt*dzt + dtt*dtt;
                if (dr2t <= fMichelSpaceRadius * fMichelSpaceRadius)
                    TrueSphereHits.push_back(hit);
            }

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
        SphereEnergy = SphereHits.energy() * fADCtoMeV;
        TrueSphereEnergy = TrueSphereHits.energy() * fADCtoMeV;

        /*

        // get space points nearby muon end point
        for (art::Ptr<recob::SpacePoint> const& p_spt : vp_spt) {

            HitPtr p_spt_hit = fop_spt2hit.at(p_spt.key());
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
        for (HitPtr const& p_hit_col : NearbyPHits) {

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
    return plane2axis[wid].space(asWire->Wire(wid));
}

ana::Hit ana::Agnochecks::GetHit(HitPtr const p_hit) {
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

// HitPtr ana::Agnochecks::GetDeepestHit(
//     HitPtrVec const& vp_hit,
//     bool increasing_z,
//     geo::View_t view
// ) {
//     if (vp_hit.empty()) return HitPtr{};

//     HitPtr DeepestHit;
//     if (geoDet == kPDVD) {

//         // test if the muon goes into the bottom volume
//         bool in_bot = false;
//         for (HitPtr p_hit : vp_hit) {
//             if (p_hit->View() != view) continue;
//             unsigned tpc = p_hit->WireID().TPC;
//             if (tpc < 8) {
//                 in_bot = true;
//                 break;
//             }
//         }

//         // basic linear regression on the hits that are in the last volume
//         unsigned n=0;
//         double mz=0, mt=0, mt2=0, mzt=0;
//         for (HitPtr p_hit : vp_hit) {
//             if (p_hit->View() != view) continue;
//             if (in_bot && p_hit->WireID().TPC >= 8) continue;
//             // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
//             double z = GetSpace(p_hit->WireID());
//             double t = p_hit->PeakTime() * fTick2cm;
//             mz += z; mt += t; mt2 += t*t, mzt += z*t;
//             n++;
//         }
//         mz /= n; mt /= n; mt2 /= n; mzt /= n;
//         double cov = mzt - mz*mt;
//         double vart = mt2 - mt*mt;

//         // z ~ m*t + p
//         double m = cov / vart;
//         double p = mz - m*mt;

//         double extrem_s =
//             in_bot 
//             ?  std::numeric_limits<double>::max() // search min_s
//             : std::numeric_limits<double>::lowest(); // search max_s

//         for (HitPtr p_hit : vp_hit) {
//             if (p_hit->View() != view) continue;
//             if (in_bot && p_hit->WireID().TPC >= 8) continue;
//             // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
//             double z = GetSpace(p_hit->WireID());
//             double t = p_hit->PeakTime() * fTick2cm;

//             // projection on the axis of the track
//             double s = (t + m*(z-p)) / (1 + m*m);
//             if (in_bot) {
//                 if (extrem_s > s) {
//                     extrem_s = s;
//                     DeepestHit = p_hit;
//                 }
//             } else {
//                 if (extrem_s < s) {
//                     extrem_s = s;
//                     DeepestHit = p_hit;
//                 }
//             }
//         }
//     } else if (geoDet == kPDHD) {

//         // test if the muon crosses the cathod
//         unsigned n_left=0, n_right=0;
//         double mz_left=0, mz_right=0;
//         for (HitPtr p_hit : vp_hit) {
//             if (p_hit->View() != view) continue;
//             unsigned tpc = p_hit->WireID().TPC;
//             // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
//             double z = GetSpace(p_hit->WireID());
//             if (tpc == 1 || tpc == 5) {
//                 mz_left += z;
//                 n_left++;
//             } else if (tpc == 2 || tpc == 6) {
//                 mz_right += z;
//                 n_right++;
//             }
//         }
//         mz_left /= n_left; mz_right /= n_right;
//         std::pair<unsigned, unsigned> tpcs;
//         if (
//             (increasing_z && mz_left > mz_right)
//             || (!increasing_z && mz_left < mz_right)
//         ) {
//             tpcs.first = 1;
//             tpcs.second = 5;
//         } else {
//             tpcs.first = 2;
//             tpcs.second = 6;
//         }

//         // basic linear regression on the hits that are in the last volume
//         unsigned n=0;
//         double mz=0, mt=0, mz2=0, mzt=0;
//         for (HitPtr p_hit : vp_hit) {
//             if (p_hit->View() != view) continue;
//             unsigned tpc = p_hit->WireID().TPC;
//             if (tpc != tpcs.first && tpc != tpcs.second) continue;
//             // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
//             double z = GetSpace(p_hit->WireID());
//             double t = p_hit->PeakTime() * fTick2cm;
//             mz += z; mt += t; mz2 += z*z; mzt += z*t;
//             n++;
//         }
//         mz /= n; mt /= n; mz2 /= n; mzt /= n;
//         double cov = mzt - mz*mt;
//         double varz = mz2 - mz*mz;

//         // t ~ m*z + p
//         double m = cov / varz;
//         double p = mt - m*mz;

//         double extrem_s =
//             increasing_z 
//             ? std::numeric_limits<double>::lowest() // search max_s
//             : std::numeric_limits<double>::max(); // search min_s

//         for (HitPtr p_hit : vp_hit) {
//             if (p_hit->View() != view) continue;
//             unsigned tpc = p_hit->WireID().TPC;
//             if (tpc != tpcs.first && tpc != tpcs.second) continue;
//             // double z = asWire->Wire(p_hit->WireID()).GetCenter().Z();
//             double z = GetSpace(p_hit->WireID());
//             double t = p_hit->PeakTime() * fTick2cm;

//             // projection on the axis of the track
//             double s = (z + m*(t-p)) / (1 + m*m);
//             if (increasing_z) {
//                 if (extrem_s < s) {
//                     extrem_s = s;
//                     DeepestHit = p_hit;
//                 }
//             } else {
//                 if (extrem_s > s) {
//                     extrem_s = s;
//                     DeepestHit = p_hit;
//                 }
//             }
//         }
//     }
//     return DeepestHit;
// }

HitPtrPair ana::Agnochecks::GetTrackEndsHits(
    HitPtrVec const& vp_hit,
    HitPtrPair *pp_cathode_crossing,
    HitPtrVec *vp_tpc_crossing,
    HitPtrVec *vp_sorted_hit,
    geo::View_t view
) {
    // minimum number of hits to perform a linear regression
    unsigned const nmin = ana::LinearRegression::nmin;

    // split volume at de cathode
    auto cathodeSide =
        geoDet == kPDVD
        ? [](geo::TPCID::TPCID_t tpc) -> int {
                return tpc >= 8 ? 1 : 0;
            }
            // geoDet == kPDHD
        : [](geo::TPCID::TPCID_t tpc) -> int {
                return (tpc == 1 || tpc == 5)
                    ? 0
                    : ((tpc == 2 || tpc == 6) ? 1 : -1);
            };


    // linear regression on each side to have a curvilinear coordinate of each hit inside a track
    // z = m*t + p
    // struct LinearRegression {
    //     unsigned n=0;
    //     double mz=0, mt=0, mz2=0, mt2=0, mzt=0;
    //     void add(double z, double t) {
    //         mz+=z; mt+=t; mz2+=z*z; mt2+=t*t; mzt+=z*t; n++;
    //     }
    //     void normalize() {
    //         mz/=n; mt/=n; mz2/=n; mt2/=n; mzt/=n;
    //     }
    //     double cov() const { return mzt - mz*mt; }
    //     double varz() const { return mz2 - mz*mz; }
    //     double vart() const { return mt2 - mt*mt; }
    //     double m() const { return n<nmin ? 0 : cov()/vart(); }
    //     double p() const { return mz - m()*mt; }
    //     double r2() const { return n<nmin ? 0 : cov()*cov() / (varz()*vart()); }
    //     double projection(double z, double t) const {
    //         return (t + m()*(z-p())) / (1 + m()*m());
    //     }
    // };

    std::vector<ana::LinearRegression> side_reg(2);
    for (HitPtr const& p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        int side = cathodeSide(p_hit->WireID().TPC);
        if (side == -1) continue; // skip hits on the other side of the anodes
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        side_reg[side].add(z, t);
    }

    // if not enough hits on both sides, return empty pair
    if (side_reg[0].n < nmin && side_reg[1].n < nmin) return {};

    // compute average from sum
    for (ana::LinearRegression& reg : side_reg) reg.normalize();

    // find the track ends on each side of the cathode
    std::vector<HitPtrPair> side_ends(2);
    std::vector<bounds<double>> side_mimmax(2);
    for (HitPtr const& p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        int side = cathodeSide(p_hit->WireID().TPC);
        if (side == -1) continue; // skip hits on the other side of the anodes
        double s = side_reg[side].projection(
            GetSpace(p_hit->WireID()),
            p_hit->PeakTime() * fTick2cm
        );
        if (s > side_mimmax[side].max) {
            side_mimmax[side].max = s;
            side_ends[side].second = p_hit;
        }
        if (s < side_mimmax[side].min) {
            side_mimmax[side].min = s;
            side_ends[side].first = p_hit;
        }
    }

    // if hits are all on one side, and no other info is requested
    if (!vp_tpc_crossing && !vp_sorted_hit) {
        if (side_reg[0].n < nmin)
            return side_ends[1];
        else if (side_reg[1].n < nmin)
            return side_ends[0];
    }
    
    // given the ends of two pieces of track, find the closest ends
    auto closestHits = [&](
        HitPtrPair const& pph1,
        HitPtrPair const& pph2,
        double dmin,
        HitPtrPair *otherHits = nullptr
    ) -> HitPtrPair {

        // all combinations of pairs
        std::vector<HitPtrPair> pairs = {
            { pph1.first, pph2.first },
            { pph1.first, pph2.second },
            { pph1.second, pph2.second },
            { pph1.second, pph2.first }
        };

        // distance squared between all pairs
        std::vector<double> d2s(4, 0);
        for (unsigned i=0; i<4; i++) {
            double zf = GetSpace(pairs[i].first->WireID());
            double tf = pairs[i].first->PeakTime() * fTick2cm;
            double zs = GetSpace(pairs[i].second->WireID());
            double ts = pairs[i].second->PeakTime() * fTick2cm;
            d2s[i] = pow(zf-zs,2) + pow(tf-ts,2);
        }

        // find all distances under dmin threshold
        std::vector<unsigned> candidates_idx;
        std::vector<double>::iterator it = d2s.begin();
        while ((it = std::find_if(
                it,
                d2s.end(),
                [dmin](double d2) { return d2 < dmin*dmin; }
            )) != d2s.end())
            candidates_idx.push_back(std::distance(d2s.begin(), it++));
        
        // no candidates found
        if (candidates_idx.empty())
            return {};

        // get the closest pair
        unsigned closest_idx = *std::min_element(
            candidates_idx.begin(),
            candidates_idx.end(),
            [&d2s](unsigned i, unsigned j) { return d2s[i] < d2s[j]; });
        
        // if outermost hits are requested, get the outermost pair
        if (otherHits) {
            unsigned other_idx = (closest_idx+2) % 4; // opposite pair
            otherHits->first = pairs[other_idx].first;
            otherHits->second = pairs[other_idx].second;
        }
        return pairs[closest_idx];
    };

    HitPtrPair trk_ends, cathode_crossing;
    if (side_reg[0].n < nmin)
        trk_ends = side_ends[1];
    else if (side_reg[1].n < nmin)
        trk_ends = side_ends[0];
    else
        cathode_crossing = closestHits(
            side_ends[0],
            side_ends[1],
            2*fCathodeGap,
            &trk_ends
        );

    // if cathode crossing info is requested
    if (pp_cathode_crossing)
        *pp_cathode_crossing = cathode_crossing;
    
    // if no tpc crossing info is needed
    if (geoDet == kPDHD || !vp_tpc_crossing) {
        return trk_ends;
    }

    std::vector<HitPtrPair> per_sec_ends(ana::n_sec[geoDet]);
    
    // if a sorted list of hits is requested
    // if (vvp_sec_sorted_hits) {
    //     // get a sorted list of hits for each section (ie. pair of TPCs)
    //     vvp_sec_sorted_hits->clear();
    //     vvp_sec_sorted_hits->resize(ana::n_sec[geoDet]);
    //     for (HitPtr const& p_hit : vp_hit) {
    //         if (p_hit->View() != view) continue;
    //         int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
    //         if (s == -1) continue;
    //         vvp_sec_sorted_hits->at(s).push_back(p_hit);
    //     }

    //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
    //         int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
    //         std::sort(
    //             vvp_sec_sorted_hits->at(s).begin(), 
    //             vvp_sec_sorted_hits->at(s).end(),
    //             [&, &reg=side_reg[side]](
    //                 HitPtr const& h1, HitPtr const& h2
    //             ) -> bool {
    //                 double const s1 = reg.projection(
    //                     GetSpace(h1->WireID()),
    //                     h1->PeakTime() * fTick2cm
    //                 );
    //                 double const s2 = reg.projection(
    //                     GetSpace(h2->WireID()),
    //                     h2->PeakTime() * fTick2cm
    //                 );
    //                 return s1 < s2;
    //             }
    //         );
    //     }

    //     // get the track ends for each section
    //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
    //         if (vvp_sec_sorted_hits->at(s).size() < nmin) continue;
    //         per_sec_ends[s].first = vvp_sec_sorted_hits->at(s).front();
    //         per_sec_ends[s].second = vvp_sec_sorted_hits->at(s).back();
    //     }

    std::vector<HitPtrVec> vp_sec_hit(ana::n_sec[geoDet]);
    for (HitPtr const& p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
        if (s == -1) continue;
        vp_sec_hit[s].push_back(p_hit);
    }
    // THIS CAUSES A SEGFAULT FOR SOME REASON???
    // auto side_sort = [&](int side) {
    //     return [&](HitPtr const& h1, HitPtr const& h2) -> bool {
    //         double const s1 = side_reg[side].projection(
    //             GetSpace(h1->WireID()),
    //             h1->PeakTime() * fTick2cm
    //         );
    //         double const s2 = side_reg[side].projection(
    //             GetSpace(h2->WireID()),
    //             h2->PeakTime() * fTick2cm
    //         );
    //         return s1 < s2;
    //     };
    // };

    if (vp_sorted_hit) {
        vp_sorted_hit->clear();

        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            HitPtrVec& vp_sec_sorted = vp_sec_hit[s];
            if (vp_sec_sorted.size() < nmin) continue;

            int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
            std::sort(
                vp_sec_sorted.begin(), 
                vp_sec_sorted.end(),
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
                    return s1 < s2;
                }
            );

            // get the track ends for each section
            per_sec_ends[s].first = vp_sec_sorted.front();
            per_sec_ends[s].second = vp_sec_sorted.back();

            vp_sorted_hit->insert(
                vp_sorted_hit->end(),
                vp_sec_sorted.begin(), vp_sec_sorted.end()
            );
        }
    } else { // only get the minmax ends of each section
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            if (vp_sec_hit[s].size() < nmin) continue;
            int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
            auto minmax = std::minmax_element(
                vp_sec_hit[s].begin(),
                vp_sec_hit[s].end(),
                [&, &reg=side_reg[side]](HitPtr const& h1, HitPtr const& h2) -> bool {
                    double const s1 = reg.projection(
                        GetSpace(h1->WireID()),
                        h1->PeakTime() * fTick2cm
                    );
                    double const s2 = reg.projection(
                        GetSpace(h2->WireID()),
                        h2->PeakTime() * fTick2cm
                    );
                    return s1 < s2;
                }
            );
            per_sec_ends[s].first = *minmax.first;
            per_sec_ends[s].second = *minmax.second;
        }
    }


    // get the hits that are at the boundaries of two sections
    HitPtrVec tpc_crossing;
    bool prev = false;
    for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
        if (per_sec_ends[s].first.isNull()) {
            prev = false;
            continue;
        }
        if (prev) {
            HitPtrPair const pp_tpc_crossing = closestHits(
                per_sec_ends[s-1], per_sec_ends[s], 2
            );
            if (pp_tpc_crossing.first.isNonnull()) {
                tpc_crossing.push_back(pp_tpc_crossing.first);
                tpc_crossing.push_back(pp_tpc_crossing.second);
            }
        }
        if ((geoDet == kPDVD && s == 3) || (geoDet == kPDHD && s == 1))
            prev = false; // cathode crossing
        else
            prev = true;
    }

    *vp_tpc_crossing = tpc_crossing;
    return trk_ends;
}

DEFINE_ART_MODULE(ana::Agnochecks)