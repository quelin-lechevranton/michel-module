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
    // float fADCtoMeV;
    float fSamplingRate; // µs/tick
    float fDriftVelocity; // cm/µs
    float fChannelPitch;

    // Data Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag tag_mcp, tag_sed, tag_wir, tag_hit, tag_clu, tag_trk, tag_spt, tag_pfp;

    // Input Parameters
    float fMichelSpaceRadius; // in cm
    float fMichelTickRadius; // in ticks
    float fNearbySpaceRadius; // in cm
    float fTrackLengthCut; // in cm
    bool fKeepOutside;

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

    ana::Hits MuonHits;
    ana::Hit MuonEndHit;
    bool MuonEndIsInWindowT, MuonEndIsInVolumeYZ;

    ana::Hits NearbyHits;
    ana::Points NearbySpacePoints;

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
};


ana::Fullchecks::Fullchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm
    fNearbySpaceRadius(p.get<float>("NearbySpaceRadius")), //in cm
    fTrackLengthCut(p.get<float>("TrackLengthCut")), // in cm
    fKeepOutside(p.get<bool>("KeepOutside"))
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
    fMichelTickRadius = fMichelSpaceRadius / fDriftVelocity / fSamplingRate;


    tEvent = tfs->make<TTree>("event","");

    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("Muon", &EventiMuon);
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

        POINTS_BRANCHES(tMuon, "Track", MuonTrackPoints),
        POINT_BRANCHES(tMuon, "EndTrack", MuonEndTrackPoint),
        HITS_BRANCHES(tMuon, "", MuonHits),
        HIT_BRANCHES(tMuon, "End", MuonEndHit),
        POINTS_BRANCHES(tMuon, "Space", MuonSpacePoints)
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
        POINTS_BRANCHES(tMuon, "NearbySpace", NearbySpacePoints)
    };
}

void ana::Fullchecks::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();

    auto const & vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    auto const & vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    auto const & vh_spt = e.getValidHandle<std::vector<recob::SpacePoint>>(tag_spt);
    auto const & vh_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(tag_pfp);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    // art::FindManyP<recob::SpacePoint> fmp_hit2spt(vh_hit, e, tag_spt);
    art::FindOneP<recob::PFParticle> fop_trk2pfp(vh_trk, e, tag_trk);
    art::FindManyP<recob::SpacePoint> fmp_pfp2spt(vh_pfp, e, tag_pfp);

    resetEvent();

    for (recob::Hit const& hit : *vh_hit) {
        if (hit.View() != geo::kW) continue;
        EventHits.push_back(ana::Hit{hit});
    }


    struct EndPoint {
        ana::Hit hit; 
        simb::MCParticle const* mcp_michel;
        size_t track_key;
        ana::Point spt;
    };
    std::vector<EndPoint> muon_endpoints;

    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        // no short tracks
        if (p_trk->Length() < fTrackLengthCut) continue;

        simb::MCParticle const* mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());

        // tracks associated to a MCTruth muon
        if (!mcp) continue;
        if (mcp->PdgCode() != 13) continue;

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());

        if (vp_hit_muon.empty()) continue;

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
            MuonEndHit = ana::Hit{*HitLowMin};
        } else {
            if (HitUpMax)
                MuonEndHit = ana::Hit{*HitUpMax};
            else
                continue;
        }
            
        // track end point is the deepest
        if (IsUpright(*p_trk))
            MuonEndTrackPoint = ana::Point{p_trk->End()};
        else
            MuonEndTrackPoint = ana::Point{p_trk->Start()};

        // fiducial cuts
        MuonEndIsInWindowT = tick_window.isInside(MuonEndHit.tick, fMichelTickRadius);
        MuonEndIsInVolumeYZ = upper_bounds.isInside(150.F, MuonEndTrackPoint.y, MuonEndTrackPoint.z, fMichelSpaceRadius);

        if (!fKeepOutside && !(MuonEndIsInWindowT && MuonEndIsInVolumeYZ)) continue;

        // we found a muon candidate!

        EventNMuon++;
        EventiMuon.push_back(iMuon);

        MuonIsAnti = mcp->PdgCode() < 0;
        MuonTrackLength = p_trk->Length();
        MuonEndProcess = mcp->EndProcess();

        // getting all muon hits
        for (art::Ptr<recob::Hit> const& p_hit_muon : vp_hit_muon) {
            if (p_hit_muon->View() != geo::kW) continue;
            MuonHits.push_back(ana::Hit{*p_hit_muon});
        }

        // and all muon track points
        for (unsigned i_tpt=0; i_tpt<p_trk->NumberTrajectoryPoints(); i_tpt++) {
            if (!p_trk->HasValidPoint(i_tpt)) continue;
            MuonTrackPoints.push_back(ana::Point{p_trk->LocationAtPoint(i_tpt)});
        }

        // and all muon space points
        double xmin = upper_bounds.x.max;
        ana::Point MuonEndSpt;
        art::Ptr<recob::PFParticle> p_pfp = fop_trk2pfp.at(p_trk.key());
        std::vector<art::Ptr<recob::SpacePoint>> v_spt_muon = fmp_pfp2spt.at(p_pfp.key());
        for (art::Ptr<recob::SpacePoint> const& p_spt : v_spt_muon) {
            MuonSpacePoints.push_back(ana::Point{p_spt->position()});

            if (p_spt->position().x() < xmin) {
                xmin = p_spt->position().x();
                MuonEndSpt = ana::Point{p_spt->position()};
            }
        }

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
        if (mcp_michel && has_numu && has_nue) {
            if (lower_bounds.isInside(mcp_michel->Position(0)) || upper_bounds.isInside(mcp_michel->Position(0))) 
                MuonHasMichel = kHasMichelInside; 
            else
                MuonHasMichel = kHasMichelOutside;
        }
        else {
            MuonHasMichel = kNoMichel;
            mcp_michel = nullptr;
        }

        // saving end point info for further analysis
        muon_endpoints.push_back({MuonEndHit, mcp_michel, p_trk.key(), MuonEndSpt});

        for (TBranch *b : brMuon) b->Fill();
    } // end of loop over tracks


    struct Nearby {
        ana::Hits hits;
        ana::Hits sphere_hits;
        unsigned true_positive;
        unsigned false_positive;
        float energy_true_positive;
        float energy_false_positive;
        ana::Points spt;
    };
    std::vector<Nearby> nearby(EventNMuon);

    // looping over all hits
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;

        art::Ptr<recob::Track> p_trk = fop_hit2trk.at(p_hit.key());
        bool from_track = p_trk && p_trk->Length() > fTrackLengthCut;

        // looping over muon end points
        for (unsigned m=0; m<EventNMuon; m++) {
            if (GetSlice(p_hit->Channel()) != muon_endpoints.at(m).hit.slice) continue;

            float dz = (map_ch_z[p_hit->Channel()] - muon_endpoints.at(m).hit.z);
            float dt = (p_hit->PeakTime() - muon_endpoints.at(m).hit.tick) * fDriftVelocity * fSamplingRate;
            float dr2 = dz*dz + dt*dt;

            bool from_another_track = from_track && muon_endpoints.at(m).track_key != p_trk->key();

            if (from_another_track) continue;
            if (dr2 > fNearbySpaceRadius * fNearbySpaceRadius) continue;

            nearby.at(m).hits.push_back(ana::Hit{*p_hit});

            if (!muon_endpoints.at(m).mcp_michel) continue;
            if (from_track) continue;
            if (dr2 > fMichelSpaceRadius * fMichelSpaceRadius) continue;

            nearby.at(m).sphere_hits.push_back(ana::Hit{*p_hit});

            // checking if the hit is associated to the michel MCParticle
            std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *muon_endpoints.at(m).mcp_michel, e, tag_hit.label());
            if (std::find(v_hit_michel.begin(), v_hit_michel.end(), &*p_hit) != v_hit_michel.end()) {
                nearby.at(m).true_positive++;
                nearby.at(m).energy_true_positive += p_hit->Integral();
            } else {
                nearby.at(m).false_positive++;
                nearby.at(m).energy_false_positive += p_hit->Integral();
            }
        }
    } // end of loop over event hits


    // get space points nearby muon end point
    for (recob::SpacePoint const& spt : *vh_spt) {
        for (unsigned m=0; m<EventNMuon; m++) {

            // should check if the space points is from an other track
            //
            //
            
            if ((muon_endpoints.at(m).spt - spt.position()).r2() > fNearbySpaceRadius * fNearbySpaceRadius) continue;

            nearby.at(m).spt.push_back(ana::Point{spt.position()});
        }
    }

    // filling the branches related to michel
    for (unsigned m=0; m<EventNMuon; m++) {
        
        resetMichel();

        NearbyHits = nearby.at(m).hits;
        NearbySpacePoints = nearby.at(m).spt;

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

            MichelHits.push_back(ana::Hit{*hit_michel});
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

void ana::Fullchecks::beginJob()
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    tick_window.min = 0;
    tick_window.max = detProp.ReadOutWindowSize();

    for (unsigned tpc=0; tpc<asGeo->NTPC(); tpc++) {
        geo::TPCID tpcid{0, tpc};

        if (tpc >= 8) {
            upper_bounds.x.min = upper_bounds.x.min < asGeo->TPC(tpcid).MinX() ? asGeo->TPC(tpcid).MinX() : upper_bounds.x.min;
            upper_bounds.x.max = upper_bounds.x.max > asGeo->TPC(tpcid).MaxX() ? asGeo->TPC(tpcid).MaxX() : upper_bounds.x.max;
            upper_bounds.y.min = upper_bounds.y.min < asGeo->TPC(tpcid).MinY() ? asGeo->TPC(tpcid).MinY() : upper_bounds.y.min;
            upper_bounds.y.max = upper_bounds.y.max > asGeo->TPC(tpcid).MaxY() ? asGeo->TPC(tpcid).MaxY() : upper_bounds.y.max;
            upper_bounds.z.min = upper_bounds.z.min < asGeo->TPC(tpcid).MinZ() ? asGeo->TPC(tpcid).MinZ() : upper_bounds.z.min;
            upper_bounds.z.max = upper_bounds.z.max > asGeo->TPC(tpcid).MaxZ() ? asGeo->TPC(tpcid).MaxZ() : upper_bounds.z.max;
        } else {
            lower_bounds.x.min = lower_bounds.x.min < asGeo->TPC(tpcid).MinX() ? asGeo->TPC(tpcid).MinX() : lower_bounds.x.min;
            lower_bounds.x.max = lower_bounds.x.max > asGeo->TPC(tpcid).MaxX() ? asGeo->TPC(tpcid).MaxX() : lower_bounds.x.max;
            lower_bounds.y.min = lower_bounds.y.min < asGeo->TPC(tpcid).MinY() ? asGeo->TPC(tpcid).MinY() : lower_bounds.y.min;
            lower_bounds.y.max = lower_bounds.y.max > asGeo->TPC(tpcid).MaxY() ? asGeo->TPC(tpcid).MaxY() : lower_bounds.y.max;
            lower_bounds.z.min = lower_bounds.z.min < asGeo->TPC(tpcid).MinZ() ? asGeo->TPC(tpcid).MinZ() : lower_bounds.z.min;
            lower_bounds.z.max = lower_bounds.z.max > asGeo->TPC(tpcid).MaxZ() ? asGeo->TPC(tpcid).MaxZ() : lower_bounds.z.max;
        }

        geo::PlaneID planeid{tpcid, 0};

        for (unsigned w=0; w<asWire->Nwires(planeid); w++) {
            geo::WireID wireid{planeid, w};
            geo::WireGeo const wiregeo = asWire->Wire(wireid);

            map_ch_z[asWire->PlaneWireToChannel(wireid)] = wiregeo.GetStart().Z();
        }

        map_tpc_ch[tpc] = {
            asWire->PlaneWireToChannel(geo::WireID{geo::PlaneID{tpcid, 0}, 0}),
            asWire->PlaneWireToChannel(geo::WireID{geo::PlaneID{tpcid, 0}, asWire->Nwires(geo::PlaneID{tpcid, 0})-1})
        };
    }
}




void ana::Fullchecks::endJob()
{
    // Implementation of optional member function here.
}


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


bool ana::Fullchecks::IsInUpperVolume(raw::ChannelID_t ch) {
    if (ch == raw::InvalidChannelID) return false;
    return ch >= map_tpc_ch[8].min;
}
bool ana::Fullchecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}

DEFINE_ART_MODULE(ana::Fullchecks)