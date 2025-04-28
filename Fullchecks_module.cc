////////////////////////////////////////////////////////////////////////
// Class:       Fullchecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Fullchecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "pdvd_utils.h"

namespace ana {
    class Fullchecks;
}


class ana::Fullchecks : public art::EDAnalyzer {
public:
    explicit Fullchecks(fhicl::ParameterSet const& p);
    Fullchecks(Fullchecks const&) = delete;
    Fullchecks(Fullchecks&&) = delete;
    Fullchecks& operator=(Fullchecks const&) = delete;
    Fullchecks& operator=(Fullchecks&&) = delete;

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
    geo::BoxBoundedGeo geoUp, geoLow;

    // protoana::ProtoDUNETruthUtils truthUtil;
    // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // art::ServiceHandle<cheat::BackTrackerService> bt_serv;

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

    // Data Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag tag_mcp, tag_sed, tag_wir, tag_hit, tag_clu, tag_trk, tag_spt, tag_pfp, tag_r3d;

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
    unsigned iEvent=0;
    unsigned EventNMuon;
    std::vector<unsigned> EventiMuon;

    ana::Hits EventHits;
    ana::Hits EventUHits, EventVHits;


    TTree* tMuon;
    unsigned iMuon=0;

    bool MuonIsAnti;
    std::string MuonEndProcess;
    int MuonHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside };

    float MuonTrackLength;
    int MuonTrackIsNotBroken;
    enum EnumIsBroken { kBadAssociation = -1, kBroken, kLastOfBroken, kNotBroken };
    ana::Points MuonTrackPoints;
    ana::Point MuonEndTrackPoint;
    ana::Point MuonTrueEndPoint;    
    float MuonTrueEndPointT;
    ana::Point MuonTrueEndMomentum;
    float MuonTrueEndEnergy;
    ana::Points MuonSpacePoints;
    ana::Point MuonEndSpacePoint;

    ana::Hits MuonHits;
    ana::Hits MuonUHits, MuonVHits;
    ana::Hit MuonEndHit;
    ana::Hit MuonEndUHit, MuonEndVHit;
    ana::Hit MuonTrueEndHit;
    bool MuonEndIsInWindowT, MuonEndIsInVolumeYZ;
    bool MuonEndHasGood3DAssociation;

    ana::Hits NearbyHits;
    ana::Hits NearbyUHits, NearbyVHits;
    // ana::Points NearbySpacePoints;
    // ana::Points NearbyHitSpacePoints;
    // ana::Points NearbyHitPoints;
    // std::vector<float> NearbyHitPointQuality;
    // ana::Points NearbyReco3DPoints;

    float MichelTrackLength;

    ana::Hits MichelHits, SphereHits;

    float MichelTrueEnergy, MichelHitEnergy, SphereEnergy;
    unsigned SphereTruePositive, SphereFalsePositive;
    float SphereEnergyTruePositive, SphereEnergyFalsePositive;


    void resetEvent();
    void resetMuon();

    bool IsInUpperVolume(raw::ChannelID_t ch);
    bool IsUpright(recob::Track const& T);
    ana::Hit GetHit(recob::Hit const& hit);

    art::Ptr<recob::Hit> GetDeepestHit(std::vector<art::Ptr<recob::Hit>>, geo::View_t = geo::kW);
};


ana::Fullchecks::Fullchecks(fhicl::ParameterSet const& p)
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
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();


    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    for (std::vector<std::string> prod : vvsProducts) {
        const std::string   process     = prod[0],
                            label       = prod[1],
                            instance    = prod[2],
                            type        = prod[3];

        const art::InputTag tag = art::InputTag(label,instance);

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
        else if (type == "recob::PFParticle")       tag_pfp = tag;
        tag_r3d = art::InputTag("reco3d", "");
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
    geoLow = geo::BoxBoundedGeo{asGeo->TPC(geo::TPCID{0, 0}).Min(), asGeo->TPC(geo::TPCID{0, asGeo->NTPC()/2-1}).Max()};
    geoUp = geo::BoxBoundedGeo{asGeo->TPC(geo::TPCID{0, asGeo->NTPC()/2}).Min(), asGeo->TPC(geo::TPCID{0, asGeo->NTPC()-1}).Max()};

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Sampling Rate: " << fSamplingRate << " µs/tick" << std::endl
        << "  Drift Velocity: " << fDriftVelocity << " cm/µs" << std::endl
        << "  Channel Pitch: " << fChannelPitch << " cm" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  Upper Bounds: " << geoUp.Min() << " -> " << geoUp.Max() << std::endl
        << "  Lower Bounds: " << geoLow.Min() << " -> " << geoLow.Max() << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl
        << "  Michel Space Radius: " << fMichelSpaceRadius << " cm (" << fMichelTickRadius << " ticks)" << std::endl
        << "  Nearby Space Radius: " << fNearbySpaceRadius << " cm" << std::endl
        << "  Coincidence Window: " << fCoincidenceWindow << " ticks" << std::endl;

    tEvent = tfs->make<TTree>("event","");

    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    EventHits.SetBranches(tEvent);
    EventUHits.SetBranches(tEvent, "U");
    EventVHits.SetBranches(tEvent, "V");

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
    MuonUHits.SetBranches(tMuon, "U");
    MuonVHits.SetBranches(tMuon, "V");
    MuonEndHit.SetBranches(tMuon, "End");
    MuonEndUHit.SetBranches(tMuon, "EndU");
    MuonEndVHit.SetBranches(tMuon, "EndV");
    MuonTrueEndHit.SetBranches(tMuon, "TrueEnd");
    MuonTrackPoints.SetBranches(tMuon, "Track");
    MuonEndTrackPoint.SetBranches(tMuon, "EndTrack");
    MuonTrueEndPoint.SetBranches(tMuon, "TrueEnd");
    tMuon->Branch("TrueEndPointT", &MuonTrueEndPointT); // ns
    tMuon->Branch("TrueEndMomentumX", &MuonTrueEndMomentum.x); // GeV
    tMuon->Branch("TrueEndMomentumY", &MuonTrueEndMomentum.y); // GeV
    tMuon->Branch("TrueEndMomentumZ", &MuonTrueEndMomentum.z); // GeV
    tMuon->Branch("TrueEndEnergy", &MuonTrueEndEnergy); // GeV
    MuonSpacePoints.SetBranches(tMuon, "Space");
    MuonEndSpacePoint.SetBranches(tMuon, "EndSpace");

    tMuon->Branch("MichelTrackLength", &MichelTrackLength); // cm
    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy); // MeV
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy); // MeV
    tMuon->Branch("SphereEnergy", &SphereEnergy); // MeV

    tMuon->Branch("SphereTruePositive", &SphereTruePositive);
    tMuon->Branch("SphereFalsePositive", &SphereFalsePositive);
    tMuon->Branch("SphereEnergyTruePositive", &SphereEnergyTruePositive);
    tMuon->Branch("SphereEnergyFalsePositive", &SphereEnergyFalsePositive);

    MichelHits.SetBranches(tMuon, "Michel");
    SphereHits.SetBranches(tMuon, "Sphere");
    NearbyHits.SetBranches(tMuon, "Nearby");
    NearbyUHits.SetBranches(tMuon, "NearbyU");
    NearbyVHits.SetBranches(tMuon, "NearbyV");
    // NearbySpacePoints.SetBranches(tMuon, "NearbySpace");
    // NearbyHitSpacePoints.SetBranches(tMuon, "NearbyHitSpace");
    // NearbyHitPoints.SetBranches(tMuon, "NearbyHit");
    // tMuon->Branch("NearbyHitPointQuality", &NearbyHitPointQuality);
    // NearbyReco3DPoints.SetBranches(tMuon, "NearbyReco3D");

}

void ana::Fullchecks::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fTick2cm = fDriftVelocity * fSamplingRate;

    auto const & vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    auto const & vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    auto const & vh_spt = e.getValidHandle<std::vector<recob::SpacePoint>>(tag_spt);
    std::vector<art::Ptr<recob::SpacePoint>> vp_spt;
    art::fill_ptr_vector(vp_spt, vh_spt);

    auto const & vh_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(tag_pfp);

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

    for (recob::Hit const& hit : *vh_hit) {
        switch (hit.View()) {
            case geo::kU: EventUHits.push_back(GetHit(hit)); break;
            case geo::kV: EventVHits.push_back(GetHit(hit)); break;
            case geo::kW: EventHits.push_back(GetHit(hit)); break;
            default: break;
        }
    }

    std::unordered_map<int, unsigned> particle_encounter;

    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        if (fLog) printf("e%ut%u\r", iEvent, p_trk->ID()), fflush(stdout);

        // no short tracks
        if (!LOG(p_trk->Length() > fTrackLengthCut)) continue;

        // simb::MCParticle const* mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());
        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);

        // tracks associated to a MCTruth muon
        if (!mcp) continue;
        if (!LOG(abs(mcp->PdgCode()) == 13)) continue;

        art::Ptr<recob::Hit> deephit = GetDeepestHit(ana::mcp2hits(mcp, vp_hit, clockData, false));
        if (!LOG(deephit)) continue;
        MuonTrueEndHit = GetHit(*deephit);

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());

        if (!LOG(vp_hit_muon.size())) continue;

        deephit = GetDeepestHit(vp_hit_muon);
        if (!LOG(deephit)) continue;
        MuonEndHit = GetHit(*deephit);

        deephit = GetDeepestHit(vp_hit_muon, geo::kU);
        if (!LOG(deephit)) continue;
        MuonEndUHit = GetHit(*deephit);

        deephit = GetDeepestHit(vp_hit_muon, geo::kV);
        if (!LOG(deephit)) continue;
        MuonEndVHit = GetHit(*deephit);


        resetMuon();
            
        // track end point is the deepest
        if (IsUpright(*p_trk))
            MuonEndTrackPoint = ana::Point{p_trk->End()};
        else
            MuonEndTrackPoint = ana::Point{p_trk->Start()};

        // fiducial cuts
        MuonEndIsInWindowT = wireWindow.isInside(MuonEndHit.tick, fMichelTickRadius);
        MuonEndIsInVolumeYZ = geoUp.InFiducialY(MuonEndTrackPoint.y, fMichelSpaceRadius) and geoUp.InFiducialZ(MuonEndTrackPoint.z, fMichelSpaceRadius);

        if (!LOG(fKeepOutside or (MuonEndIsInWindowT and MuonEndIsInVolumeYZ))) continue;

        // we found a muon candidate!
        if (fLog) printf("\t\033[1;93m" "e%um%u (%u)" "\033[0m\n", iEvent, EventNMuon, iMuon);

        EventiMuon.push_back(iMuon);

        MuonIsAnti = mcp->PdgCode() < 0;
        MuonTrackLength = p_trk->Length();

        std::vector<art::Ptr<recob::Track>> vp_trk_from_mcp = mcp2trks(mcp, vp_trk, clockData, fmp_trk2hit);
        // std::cout << "trk#" << p_trk->ID() << " mu#" << mcp->TrackId() << " encounter#" << ++particle_encounter[mcp->TrackId()] << " #mcp2trk:" << vp_trk_from_mcp.size();
        if (vp_trk_from_mcp.size()) {
            // bool isin = std::find(vp_trk_from_mcp.begin(), vp_trk_from_mcp.end(), p_trk) != vp_trk_from_mcp.end();
            // std::cout << " \033[1;9" << (isin ? 2 : 1) << "m" << (isin ? "in" : "not in") << "\033[0m";
            // std::cout << " [ ";
            for (art::Ptr<recob::Track> p : vp_trk_from_mcp) std::cout << p->ID() << ", ";
            // std::cout << "]" << std::endl;
            if (vp_trk_from_mcp.size() == 1)
                MuonTrackIsNotBroken = kNotBroken;
            else {
                bool IsDeepestTrack = true;
                for (art::Ptr<recob::Track> p_trk_from_mcp : vp_trk_from_mcp)
                    IsDeepestTrack = IsDeepestTrack && (MuonEndTrackPoint.x <= p_trk_from_mcp->Start().X() && MuonEndTrackPoint.x <= p_trk_from_mcp->End().X());
                if (IsDeepestTrack)
                    MuonTrackIsNotBroken = kLastOfBroken;
                else
                    MuonTrackIsNotBroken = kBroken;
            }
        } else MuonTrackIsNotBroken = kBadAssociation;

        MuonEndProcess = mcp->EndProcess();
        MuonTrueEndPoint = ana::Point{mcp->EndPosition().Vect()};
        MuonTrueEndPointT = mcp->EndT();
        MuonTrueEndMomentum = ana::Point{mcp->EndMomentum().Vect()};
        MuonTrueEndEnergy = mcp->EndE();

        // getting all muon hits
        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            switch (p_hit_muon->View()) {
                case geo::kU: MuonUHits.push_back(GetHit(*p_hit_muon)); break;
                case geo::kV: MuonVHits.push_back(GetHit(*p_hit_muon)); break;
                case geo::kW: MuonHits.push_back(GetHit(*p_hit_muon)); break;
                default: break;
            }
        }

        // and all muon track points
        for (unsigned i_tpt=0; i_tpt<p_trk->NumberTrajectoryPoints(); i_tpt++) {
            if (!p_trk->HasValidPoint(i_tpt)) continue;
            MuonTrackPoints.push_back(p_trk->LocationAtPoint(i_tpt));
        }

        // and all muon space points
        MuonEndSpacePoint = ana::Point{geoUp.MaxX(), 0., 0.};
        art::Ptr<recob::PFParticle> p_pfp = fop_trk2pfp.at(p_trk.key());
        std::vector<art::Ptr<recob::SpacePoint>> v_spt_muon = fmp_pfp2spt.at(p_pfp.key());
        for (art::Ptr<recob::SpacePoint> const& p_spt : v_spt_muon) {
            MuonSpacePoints.push_back(p_spt->position());

            if (p_spt->position().x() < MuonEndSpacePoint.x)
                MuonEndSpacePoint = ana::Point{p_spt->position()};
        }
        MuonEndHasGood3DAssociation = MuonEndSpacePoint.x != geoUp.MaxX() and abs(MuonEndHit.z - MuonEndSpacePoint.z) < fCoincidenceRadius;

        // a decaying muon has nu_mu, nu_e and elec as last daughters
        bool has_numu = false, has_nue = false;
        simb::MCParticle const* mcp_michel = nullptr;
        for (int i_dau=mcp->NumberDaughters()-3; i_dau<mcp->NumberDaughters(); i_dau++) {
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    
            if (!mcp_dau) continue;

            switch (abs(mcp_dau->PdgCode())) {
                case 14: has_numu = true; break;
                case 12: has_nue = true; break;
                case 11: mcp_michel = mcp_dau; break;
            }
        }
        if (mcp_michel and has_numu and has_nue) {
            if (geoLow.ContainsPosition(mcp_michel->Position(0).Vect()) or geoUp.ContainsPosition(mcp_michel->Position(0).Vect())) 
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

                MichelHits.push_back(GetHit(*p_hit_michel));
            }
            MichelHitEnergy = MichelHits.energy();
        }
        else MuonHasMichel = kNoMichel;


        // get induction hits nearby muon end point
        float induction_pitch = 0.765; // cm/channel
        for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
            if (p_hit->View() != geo::kU && p_hit->View() != geo::kV) continue;

            art::Ptr<recob::Track> p_hit_trk = fop_hit2trk.at(p_hit.key());
            bool from_track = p_hit_trk and p_hit_trk->Length() > fTrackLengthCut;

            if (from_track and (p_trk.key() != p_hit_trk.key())) continue;

            if (p_hit->View() == geo::kU) {
                float dz = (p_hit->Channel() - MuonEndUHit.channel) * induction_pitch;
                float dt = (p_hit->PeakTime() - MuonEndUHit.tick) * fTick2cm;
                float dr2 = dz*dz + dt*dt;

                if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

                NearbyUHits.push_back(GetHit(*p_hit));
            }
            if (p_hit->View() == geo::kV) {
                float dz = (p_hit->Channel() - MuonEndVHit.channel) * induction_pitch;
                float dt = (p_hit->PeakTime() - MuonEndVHit.tick) * fTick2cm;
                float dr2 = dz*dz + dt*dt;

                if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

                NearbyVHits.push_back(GetHit(*p_hit));
            }
        }


        // get all hits nearby muon end point
        std::vector<art::Ptr<recob::Hit>> NearbyPHits;
        for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;

            art::Ptr<recob::Track> p_hit_trk = fop_hit2trk.at(p_hit.key());
            bool from_track = p_hit_trk and p_hit_trk->Length() > fTrackLengthCut;

            if (from_track and (p_trk.key() != p_hit_trk.key())) continue;

            ana::Hit hit = GetHit(*p_hit);
            if (hit.slice != MuonEndHit.slice) continue;

            float dz = (hit.z - MuonEndHit.z);
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
            if (dr2 > fMichelSpaceRadius * fMichelSpaceRadius) continue;

            SphereHits.push_back(hit);

            // checking if the hit is associated to the michel MCParticle
            // std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());
            // if (std::find(v_hit_michel.begin(), v_hit_michel.end(), &*p_hit) != v_hit_michel.end()) {
            if (MichelHits.find(*p_hit) != MichelHits.N) {
                SphereTruePositive++;
                SphereEnergyTruePositive += hit.adc;
            } else {
                SphereFalsePositive++;
                SphereEnergyFalsePositive += hit.adc;
            }
        } // end of loop over event hits
        SphereEnergy = SphereHits.energy();

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

void ana::Fullchecks::beginJob() {}
void ana::Fullchecks::endJob() {}


void ana::Fullchecks::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
    EventUHits.clear();
    EventVHits.clear();
}
void ana::Fullchecks::resetMuon() {
    MuonTrackPoints.clear();
    MuonSpacePoints.clear();
    MuonHits.clear();
    MuonUHits.clear();
    MuonVHits.clear();

    NearbyHits.clear();
    NearbyUHits.clear();
    NearbyVHits.clear();
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

    SphereTruePositive = 0;
    SphereFalsePositive = 0;
    SphereEnergyTruePositive = 0;
    SphereEnergyFalsePositive = 0;
}

ana::Hit ana::Fullchecks::GetHit(recob::Hit const& hit) {
    geo::WireID wireid = hit.WireID();
    geo::WireGeo wiregeo = asWire->Wire(wireid);

    unsigned slice;
    float z;
    switch (hit.View()) {
        case geo::kW:
            slice = 2*(wireid.TPC/4) + wireid.TPC%2;
            z = wiregeo.GetCenter().Z();
            break;
        case geo::kU:
            slice = wireid.TPC;
            z = wiregeo.GetCenter().Y() * sqrt(3) + wireid.TPC < 8 ? -wiregeo.GetCenter().Z() : +wiregeo.GetCenter().Z();
            break;
        case geo::kV:
            slice = wireid.TPC;
            z = wiregeo.GetCenter().Y() * sqrt(3) + wireid.TPC < 8 ? +wiregeo.GetCenter().Z() : -wiregeo.GetCenter().Z();
            break;
        default:
            return ana::Hit{};
    }
    return ana::Hit{
        slice,
        z,
        hit.Channel(),
        hit.PeakTime(),
        hit.Integral()
    };
}
bool ana::Fullchecks::IsInUpperVolume(raw::ChannelID_t ch) {
    if (ch == raw::InvalidChannelID) return false;
    return ch >= asWire->PlaneWireToChannel(geo::WireID{geo::PlaneID{geo::TPCID{0, 8}, geo::kU}, 0});
}
bool ana::Fullchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}



art::Ptr<recob::Hit> ana::Fullchecks::GetDeepestHit(std::vector<art::Ptr<recob::Hit>> vp_hit, geo::View_t view) {
    float TickUpMax = wireWindow.min, TickLowMin = wireWindow.max;
    art::Ptr<recob::Hit> HitUpMax, HitLowMin;

    // tacking min of ticks in lower volume and max of ticks in upper volume
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        if (p_hit->View() != view) continue;

        if (IsInUpperVolume(p_hit->Channel())) {
            if (p_hit->PeakTime() > TickUpMax) {
                TickUpMax = p_hit->PeakTime();
                HitUpMax = p_hit;
            }
        } else {
            if (p_hit->PeakTime() < TickLowMin) {
                TickLowMin = p_hit->PeakTime();
                HitLowMin = p_hit;
            }
        }
    }

    // if there is hits in lower volume, muon end is in upper volume
    // else muon end is in upper volume
    if (HitLowMin) return HitLowMin;
    else {
        if (HitUpMax) return HitUpMax;
        else return art::Ptr<recob::Hit>{};
    }
}



DEFINE_ART_MODULE(ana::Fullchecks)