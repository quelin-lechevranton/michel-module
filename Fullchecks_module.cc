////////////////////////////////////////////////////////////////////////
// Class:       Fullchecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Fullchecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "data_struct.h"

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

    protoana::ProtoDUNETruthUtils truthUtil;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Detector Properties
    // float fADC2MeV;
    float fTick2cm; // cm/tick
    float fSamplingRate; // µs/tick
    float fDriftVelocity; // cm/µs
    float fChannelPitch;

    bounds<float> tick_window;
    bounds3D<float> lower_bounds, upper_bounds;
    // std::map<int,ana::bounds<unsigned>> map_tpc_ch;
    // std::map<int,float> map_ch_z;

    // Data Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag tag_mcp, tag_sed, tag_wir, tag_hit, tag_clu, tag_trk, tag_spt, tag_pfp;

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


    TTree* tMuon;
    unsigned iMuon=0;
    std::vector<TBranch*> brMuon, brMichel;

    bool MuonIsAnti;
    std::string MuonEndProcess;
    int MuonHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside };

    float MuonTrackLength;
    ana::Points MuonTrackPoints;
    ana::Point MuonEndTrackPoint;
    ana::Points MuonSpacePoints;
    ana::Point MuonEndSpacePoint;

    ana::Hits MuonHits;
    ana::Hit MuonEndHit;
    bool MuonEndIsInWindowT, MuonEndIsInVolumeYZ;
    bool MuonEndHasGood3DAssociation;

    ana::Hits NearbyHits;
    ana::Points NearbySpacePoints;
    ana::Points NearbyHitSpacePoints;
    ana::Points NearbyHitPoints;
    std::vector<float> NearbyHitPointQuality;

    float MichelTrackLength;

    ana::Hits MichelHits, SphereHits;

    float MichelTrueEnergy, MichelHitEnergy, SphereEnergy;
    unsigned SphereTruePositive, SphereFalsePositive;
    float SphereEnergyTruePositive, SphereEnergyFalsePositive;


    void resetEvent();
    void resetMuon();
    void resetMichel();

    bool IsInUpperVolume(raw::ChannelID_t ch);
    bool IsUpright(recob::Track const& T);
    ana::Hit GetHit(recob::Hit const& hit);
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
    }

    fChannelPitch = geo::WireGeo::WirePitch(
        asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 0}),
        asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 1})
    );
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fTick2cm = fDriftVelocity * fSamplingRate;
    fMichelTickRadius = fMichelSpaceRadius / fDriftVelocity / fSamplingRate;

    tick_window.min = 0;
    tick_window.max = detProp.ReadOutWindowSize();

    for (unsigned tpc=0; tpc<asGeo->NTPC(); tpc++) {
        geo::TPCID tpcid{0, tpc};

        if (tpc >= 8U) {
            upper_bounds.x.min = upper_bounds.x.min > asGeo->TPC(tpcid).MinX() ? asGeo->TPC(tpcid).MinX() : upper_bounds.x.min;
            upper_bounds.x.max = upper_bounds.x.max < asGeo->TPC(tpcid).MaxX() ? asGeo->TPC(tpcid).MaxX() : upper_bounds.x.max;
            upper_bounds.y.min = upper_bounds.y.min > asGeo->TPC(tpcid).MinY() ? asGeo->TPC(tpcid).MinY() : upper_bounds.y.min;
            upper_bounds.y.max = upper_bounds.y.max < asGeo->TPC(tpcid).MaxY() ? asGeo->TPC(tpcid).MaxY() : upper_bounds.y.max;
            upper_bounds.z.min = upper_bounds.z.min > asGeo->TPC(tpcid).MinZ() ? asGeo->TPC(tpcid).MinZ() : upper_bounds.z.min;
            upper_bounds.z.max = upper_bounds.z.max < asGeo->TPC(tpcid).MaxZ() ? asGeo->TPC(tpcid).MaxZ() : upper_bounds.z.max;
        } else {
            lower_bounds.x.min = lower_bounds.x.min > asGeo->TPC(tpcid).MinX() ? asGeo->TPC(tpcid).MinX() : lower_bounds.x.min;
            lower_bounds.x.max = lower_bounds.x.max < asGeo->TPC(tpcid).MaxX() ? asGeo->TPC(tpcid).MaxX() : lower_bounds.x.max;
            lower_bounds.y.min = lower_bounds.y.min > asGeo->TPC(tpcid).MinY() ? asGeo->TPC(tpcid).MinY() : lower_bounds.y.min;
            lower_bounds.y.max = lower_bounds.y.max < asGeo->TPC(tpcid).MaxY() ? asGeo->TPC(tpcid).MaxY() : lower_bounds.y.max;
            lower_bounds.z.min = lower_bounds.z.min > asGeo->TPC(tpcid).MinZ() ? asGeo->TPC(tpcid).MinZ() : lower_bounds.z.min;
            lower_bounds.z.max = lower_bounds.z.max < asGeo->TPC(tpcid).MaxZ() ? asGeo->TPC(tpcid).MaxZ() : lower_bounds.z.max;
        }
    }

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Sampling Rate: " << fSamplingRate << " µs/tick" << std::endl
        << "  Drift Velocity: " << fDriftVelocity << " cm/µs" << std::endl
        << "  Channel Pitch: " << fChannelPitch << " cm" << std::endl
        << "  Tick Window: " << tick_window << std::endl
        << "  Upper Bounds: " << upper_bounds << std::endl
        << "  Lower Bounds: " << lower_bounds << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl
        << "  Michel Space Radius: " << fMichelSpaceRadius << " cm" << std::endl
        << "  Michel Tick Radius: " << fMichelTickRadius << " ticks" << std::endl
        << "  Nearby Space Radius: " << fNearbySpaceRadius << " cm" << std::endl
        << "  Coincidence Window: " << fCoincidenceWindow << " ticks" << std::endl;

    tEvent = tfs->make<TTree>("event","");

    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    HITS_BRANCHES(tEvent, "", EventHits);

    tMuon = tfs->make<TTree>("muon","");

    brMuon = {
        tMuon->Branch("iEvent", &iEvent),
        tMuon->Branch("iMuon", &iMuon),
        tMuon->Branch("iMuonInEvent", &EventNMuon),

        tMuon->Branch("IsAnti", &MuonIsAnti),
        tMuon->Branch("EndProcess", &MuonEndProcess),
        tMuon->Branch("HasMichel", &MuonHasMichel),

        tMuon->Branch("TrackLength", &MuonTrackLength),
        tMuon->Branch("EndIsInWindowT", &MuonEndIsInWindowT),
        tMuon->Branch("EndIsInVolumeYZ", &MuonEndIsInVolumeYZ),
        tMuon->Branch("EndHasGood3DAssociation", &MuonEndHasGood3DAssociation),

        HITS_BRANCHES(tMuon, "", MuonHits),
        HIT_BRANCHES(tMuon, "End", MuonEndHit),
        POINTS_BRANCHES(tMuon, "Track", MuonTrackPoints),
        POINT_BRANCHES(tMuon, "EndTrack", MuonEndTrackPoint),
        POINTS_BRANCHES(tMuon, "Space", MuonSpacePoints),
        POINT_BRANCHES(tMuon, "EndSpace", MuonEndSpacePoint)
    };

    brMichel = {
        tMuon->Branch("MichelTrackLength", &MichelTrackLength),
        tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy),
        tMuon->Branch("MichelHitEnergy", &MichelHitEnergy),
        tMuon->Branch("SphereEnergy", &SphereEnergy),

        tMuon->Branch("SphereTruePositive", &SphereTruePositive),
        tMuon->Branch("SphereFalsePositive", &SphereFalsePositive),
        tMuon->Branch("SphereEnergyTruePositive", &SphereEnergyTruePositive),
        tMuon->Branch("SphereEnergyFalsePositive", &SphereEnergyFalsePositive),

        HITS_BRANCHES(tMuon, "Michel", MichelHits),
        HITS_BRANCHES(tMuon, "Sphere", SphereHits),
        HITS_BRANCHES(tMuon, "Nearby", NearbyHits),
        POINTS_BRANCHES(tMuon, "NearbySpace", NearbySpacePoints),
        POINTS_BRANCHES(tMuon, "NearbyHitSpace", NearbyHitSpacePoints),
        POINTS_BRANCHES(tMuon, "NearbyHit", NearbyHitPoints),
        tMuon->Branch("NearbyHitPointQuality", &NearbyHitPointQuality)
    };
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

    resetEvent();

    for (recob::Hit const& hit : *vh_hit) {
        if (hit.View() != geo::kW) continue;
        EventHits.push_back(GetHit(hit));
    }


    struct EndPoint {
        ana::Hit hit; 
        simb::MCParticle const* mcp_michel;
        size_t trk_key;
        ana::Point spt;
    };
    std::vector<struct EndPoint> muon_endpoints;

    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        if (fLog) printf("e%ut%u\r", iEvent, p_trk->ID()), fflush(stdout);

        // no short tracks
        if (!LOG(p_trk->Length() > fTrackLengthCut)) continue;

        simb::MCParticle const* mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());

        // tracks associated to a MCTruth muon
        if (!mcp) continue;
        if (!LOG(abs(mcp->PdgCode()) == 13)) continue;

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());

        if (!LOG(vp_hit_muon.size())) continue;

        float TickUpMax = tick_window.min, TickLowMin = tick_window.max;
        art::Ptr<recob::Hit> HitUpMax, HitLowMin;

        // tacking min of ticks in lower volume and max of ticks in upper volume
        for (art::Ptr<recob::Hit> const& p_hit : vp_hit_muon) {
            if (p_hit->View() != geo::kW) continue;

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

        resetMuon();

        // vp_spt;

        // if there is hits in lower volume, muon end is in upper volume
        // else muon end is in upper volume
        if (HitLowMin) {
            MuonEndHit = GetHit(*HitLowMin);
        } else {
            if (HitUpMax)
                MuonEndHit = GetHit(*HitUpMax);
            else {
                LOG(!"no collection hit in volume");
                continue;
            }
        }
            
        // track end point is the deepest
        if (IsUpright(*p_trk))
            MuonEndTrackPoint = ana::Point{p_trk->End()};
        else
            MuonEndTrackPoint = ana::Point{p_trk->Start()};

        // fiducial cuts
        MuonEndIsInWindowT = tick_window.isInside(MuonEndHit.tick, fMichelTickRadius);
        MuonEndIsInVolumeYZ = upper_bounds.y.isInside(MuonEndTrackPoint.y, fMichelSpaceRadius) and upper_bounds.z.isInside(MuonEndTrackPoint.z, fMichelSpaceRadius);

        if (!LOG(fKeepOutside or (MuonEndIsInWindowT and MuonEndIsInVolumeYZ))) continue;

        // we found a muon candidate!
        if (fLog) printf("\t\033[1;93m" "e%um%u (%u)" "\033[0m\n", iEvent, EventNMuon, iMuon);

        EventiMuon.push_back(iMuon);

        MuonIsAnti = mcp->PdgCode() < 0;
        MuonTrackLength = p_trk->Length();
        MuonEndProcess = mcp->EndProcess();

        // getting all muon hits
        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            if (p_hit_muon->View() != geo::kW) continue;
            MuonHits.push_back(GetHit(*p_hit_muon));
        }

        // and all muon track points
        for (unsigned i_tpt=0; i_tpt<p_trk->NumberTrajectoryPoints(); i_tpt++) {
            if (!p_trk->HasValidPoint(i_tpt)) continue;
            MuonTrackPoints.push_back(p_trk->LocationAtPoint(i_tpt));
        }

        // and all muon space points
        MuonEndSpacePoint = ana::Point{upper_bounds.x.max, 0, 0};
        art::Ptr<recob::PFParticle> p_pfp = fop_trk2pfp.at(p_trk.key());
        std::vector<art::Ptr<recob::SpacePoint>> v_spt_muon = fmp_pfp2spt.at(p_pfp.key());
        for (art::Ptr<recob::SpacePoint> const& p_spt : v_spt_muon) {
            MuonSpacePoints.push_back(p_spt->position());

            if (p_spt->position().x() < MuonEndSpacePoint.x)
                MuonEndSpacePoint = ana::Point{p_spt->position()};
        }
        MuonEndHasGood3DAssociation = MuonEndSpacePoint.x != upper_bounds.x.max and abs(MuonEndHit.z - MuonEndSpacePoint.z) < fCoincidenceRadius;

        // a decaying muon has nu_mu, nu_e and el as last daughters
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
            if (lower_bounds.isInside(mcp_michel->Position(0)) or upper_bounds.isInside(mcp_michel->Position(0))) 
                MuonHasMichel = kHasMichelInside; 
            else
                MuonHasMichel = kHasMichelOutside;
        }
        else {
            MuonHasMichel = kNoMichel;
            mcp_michel = nullptr;
        }

        // saving end point info for further analysis
        muon_endpoints.push_back({MuonEndHit, mcp_michel, p_trk.key(), MuonEndSpacePoint});

        for (TBranch *b : brMuon) b->Fill();
        iMuon++;
        EventNMuon++;
    } // end of loop over tracks


    struct Nearby {
        // ana::Hits hits;
        std::vector<art::Ptr<recob::Hit>> vp_hit;
        ana::Hits sphere_hits;
        unsigned true_positive;
        unsigned false_positive;
        float energy_true_positive;
        float energy_false_positive;
        ana::Points spt;
        ana::Points hit_spt;
    };
    std::vector<struct Nearby> nearby(EventNMuon);

    // looping over all hits
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;

        art::Ptr<recob::Track> p_trk = fop_hit2trk.at(p_hit.key());
        bool from_track = p_trk and p_trk->Length() > fTrackLengthCut;

        // looping over muon end points
        for (unsigned m=0; m<EventNMuon; m++) {
            ana::Hit hit = GetHit(*p_hit);
            if (hit.slice != muon_endpoints.at(m).hit.slice) continue;

            float dz = (hit.z - muon_endpoints.at(m).hit.z);
            float dt = (hit.tick - muon_endpoints.at(m).hit.tick) * fTick2cm;
            float dr2 = dz*dz + dt*dt;

            bool from_another_track = from_track and muon_endpoints.at(m).trk_key != p_trk.key();

            if (from_another_track) continue;
            if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

            // nearby.at(m).hits.push_back(hit);
            nearby.at(m).vp_hit.push_back(p_hit);


            art::Ptr<recob::SpacePoint> p_spt = fop_hit2spt.at(p_hit.key());
            if (p_spt)
                nearby.at(m).hit_spt.push_back(p_spt->position());
            else 
                nearby.at(m).hit_spt.push_back(ana::Point{0.F, 0.F, 0.F});


            if (!muon_endpoints.at(m).mcp_michel) continue;
            if (from_track) continue;
            if (dr2 > fMichelSpaceRadius * fMichelSpaceRadius) continue;

            nearby.at(m).sphere_hits.push_back(hit);

            // checking if the hit is associated to the michel MCParticle
            std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *muon_endpoints.at(m).mcp_michel, e, tag_hit.label());
            if (std::find(v_hit_michel.begin(), v_hit_michel.end(), &*p_hit) != v_hit_michel.end()) {
                nearby.at(m).true_positive++;
                nearby.at(m).energy_true_positive += hit.adc;
            } else {
                nearby.at(m).false_positive++;
                nearby.at(m).energy_false_positive += hit.adc;
            }
        }
    } // end of loop over event hits


    // get space points nearby muon end point
    for (art::Ptr<recob::SpacePoint> const& p_spt : vp_spt) {

        art::Ptr<recob::Hit> p_hit = fop_spt2hit.at(p_spt.key());
        art::Ptr<recob::Track> p_trk = fop_hit2trk.at(p_hit.key());

        // std::cout << "spt: id" << p_spt->ID() << " (" << p_spt->position().x() << ", " << p_spt->position().y() << ", " << p_spt->position().z() << ")"
        //     << " w/ hit: " << char('U' + p_hit->View()) << " (" << p_hit->Channel() << ", " << p_hit->PeakTime() << ")"
        //     << std::endl;

        for (unsigned m=0; m<EventNMuon; m++) {
            ana::Point dpt = muon_endpoints.at(m).spt - p_spt->position();
            dpt.x = 0;
            if (dpt.r2() > fNearbySpaceRadius * fNearbySpaceRadius) continue;
            if (p_trk and p_trk.key() != muon_endpoints.at(m).trk_key and p_trk->Length() > fTrackLengthCut) continue;

            nearby.at(m).spt.push_back(p_spt->position());
        }
    }

    // filling the branches related to michel
    for (unsigned m=0; m<EventNMuon; m++) {

        resetMichel();

        // /* RECREATING SPACE POINTS FROM NEARBY HITS ATTEMPT

        // ana::Points NearbySpaceHits;
        std::cout << "mu#" << m << " " << nearby.at(m).vp_hit.size() << " nearby hits" << std::endl;
        std::cout << nearby.at(m).spt.size() << " nearby space points" << std::endl;
        for (art::Ptr<recob::Hit> const& p_hit_col : nearby.at(m).vp_hit) {

            art::Ptr<recob::SpacePoint> p_spt = fop_hit2spt.at(p_hit_col.key());
            std::cout << "\033[1m" "hit w/ " << (p_spt ? "1" : "0") << " spt: " "\033[0m";
            if (p_spt) std::cout << " (" << p_spt->position().x() << ", " << p_spt->position().y() << ", " << p_spt->position().z() << ")";
            std::cout << std::endl;

            // assuming MuonEndHasGood3DAssociation is true
            geo::WireGeo const wiregeo_col = asWire->Wire(p_hit_col->WireID());
            float y;
            float z = wiregeo_col.GetCenter().Z();
            float x = muon_endpoints.at(m).spt.x + (p_hit_col->PeakTime() - muon_endpoints.at(m).hit.tick) * fTick2cm;

            struct Coincidence {
                // ana::Point pt;
                float y;
                recob::Hit const* hit;
            };
            std::vector<struct Coincidence> V_coincidences, U_coincidences;


            std::cout << "   EndSpt: " << muon_endpoints.at(m).spt << "  dtick: " << p_hit_col->PeakTime() - muon_endpoints.at(m).hit.tick << std::endl;
            geo::TPCGeo const tpcgeo = asGeo->TPC(geo::TPCID{0, p_hit_col->WireID().TPC});
            printf("  TPC %u bounds: [%lf, %lf] x [%lf, %lf] x [%lf, %lf]", p_hit_col->WireID().TPC, tpcgeo.MinX(), tpcgeo.MaxX(), tpcgeo.MinY(), tpcgeo.MaxY(), tpcgeo.MinZ(), tpcgeo.MaxZ());

            for (recob::Hit const& hit_ind : *vh_hit) {
                if (!(hit_ind.View() == geo::kU or hit_ind.View() == geo::kV)) continue;
                if (abs(hit_ind.PeakTime() - p_hit_col->PeakTime()) > fCoincidenceWindow) continue;

                geo::WireGeo const wiregeo_ind = asWire->Wire(hit_ind.WireID());
                float co_y = geo::WiresIntersection(wiregeo_col, wiregeo_ind).y();
                // geo::Point_t pt = geo::WiresIntersection(wiregeo_col, wiregeo_ind);
                // // float co1_y = pt.y();
                // ana::Point wpt{pt};

                // auto const [start_ind, end_ind] = asWire->WireEndPoints(hit_ind.WireID());

                // double s = (z - start_ind.z()) / (end_ind.z() - start_ind.z());
                // float co_y = start_ind.y() + s * (end_ind.y() - start_ind.y());

                // std::cout << "\t\t wireIntersec pt: " << wpt << " vs. handmade pt: " << ana::Point(x, co_y, z) << std::endl;

                // if (!(upper_bounds.y.isInside(co_y) or lower_bounds.y.isInside(co_y))) continue;

                std::cout << " .  pt: " << ana::Point(x,co_y,z) << std::endl;

                if (!LOG(tpcgeo.ContainsPosition(geo::Point_t(x,co_y,z), 1))) continue;


                struct Coincidence co = {co_y, &hit_ind};
                // std::cout << " [" << char('U' + hit_ind.View()) << ", " << pt << "]";

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
                        if (dy < fCoincidenceRadius && abs(bary_y - muon_endpoints.at(m).spt.y) < fNearbySpaceRadius) {
                            barys_y.push_back(bary_y);
                        }
                        continue;
                    } 

                    if (dy < fCoincidenceRadius && abs(bary_y - muon_endpoints.at(m).spt.y) < fNearbySpaceRadius) {
                        has_good_coincidence = true;
                        barys_y.clear();
                        barys_y.push_back(bary_y);
                    } else {
                        if (dy < min_dy) {
                            min_dy = dy;
                            barys_y.clear();
                            barys_y.push_back(bary_y);
                        }
                    }

                    // if ((bary - muon_endpoints.at(m).spt).r2() > fNearbySpaceRadius * fNearbySpaceRadius) continue;

                    // std::cout << "  " << bary << std::endl;
                }
            }
            // std::cout << "  best bary: " << best_bary << " w/ dy: " << min_dy << std::endl;
            y = barys_y.empty() ? 0 : std::accumulate(barys_y.begin(), barys_y.end(), 0.F) / barys_y.size();
            ana::Point pt{x, y, z};
            std::cout << "  custom spt: " << pt << " oof " << barys_y.size() << " w/ " << (has_good_coincidence ? "good" : "bad") << " coincidence" << std::endl;
            NearbyHitPoints.push_back(pt);
        }
        // std::cout << no_coincidence << " hits w/o coincidence" << std::endl;

        // */

        // NearbyHits = nearby.at(m).hits;
        for (art::Ptr<recob::Hit> const& p_hit : nearby.at(m).vp_hit)
            NearbyHits.push_back(GetHit(*p_hit));
            
        NearbySpacePoints = nearby.at(m).spt;
        NearbyHitSpacePoints = nearby.at(m).hit_spt;

        if (!muon_endpoints.at(m).mcp_michel) {
            for (TBranch *b : brMichel) b->Fill();
            continue;
        }

        recob::Track const * trk_michel = truthUtil.GetRecoTrackFromMCParticle(clockData, *muon_endpoints.at(m).mcp_michel, e, tag_trk.label());
        if (trk_michel) 
            MichelTrackLength = trk_michel->Length();
        else
            MichelTrackLength = 0;

        MichelTrueEnergy = (muon_endpoints.at(m).mcp_michel->E() - muon_endpoints.at(m).mcp_michel->Mass()) * 1e3;

        std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *muon_endpoints.at(m).mcp_michel, e, tag_hit.label());
        for (const recob::Hit* hit_michel : v_hit_michel) {
            if (hit_michel->View() != geo::kW) continue;

            MichelHits.push_back(GetHit(*hit_michel));
        }
        MichelHitEnergy = MichelHits.energy();

        SphereHits = nearby.at(m).sphere_hits;
        SphereEnergy = SphereHits.energy();

        SphereTruePositive = nearby.at(m).true_positive;
        SphereFalsePositive = nearby.at(m).false_positive;
        SphereEnergyTruePositive = nearby.at(m).energy_true_positive;
        SphereEnergyFalsePositive = nearby.at(m).energy_false_positive;

        for (TBranch *b : brMichel) b->Fill();
    }

    tMuon->SetEntries(brMuon.front()->GetEntries());
    tEvent->Fill();
    iEvent++;
}

void ana::Fullchecks::beginJob() {}
void ana::Fullchecks::endJob() {}


void ana::Fullchecks::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
}
void ana::Fullchecks::resetMuon() {
    MuonTrackPoints.clear();
    MuonSpacePoints.clear();
    MuonHits.clear();
}
void ana::Fullchecks::resetMichel() {
    NearbyHits.clear();
    NearbySpacePoints.clear();
    NearbyHitSpacePoints.clear();
    NearbyHitPoints.clear();
    NearbyHitPointQuality.clear();

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
    return ana::Hit{
        unsigned(2*(hit.WireID().TPC/4) + hit.WireID().TPC%2),
        float(asWire->Wire(hit.WireID()).GetCenter().Z()),
        hit.Channel(),
        hit.PeakTime(),
        hit.Integral()
    };
}
bool ana::Fullchecks::IsInUpperVolume(raw::ChannelID_t ch) {
    if (ch == raw::InvalidChannelID) return false;
    return ch >= asWire->PlaneWireToChannel(geo::WireID{geo::PlaneID{geo::TPCID{0, 8}, geo::kW}, 0});
}
bool ana::Fullchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}


DEFINE_ART_MODULE(ana::Fullchecks)