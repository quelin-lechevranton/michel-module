////////////////////////////////////////////////////////////////////////
// Class:       Endchecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Trackchecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"

namespace ana {
    class Endchecks;
}

using HitPtr = art::Ptr<recob::Hit>;
using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

class ana::Endchecks : public art::EDAnalyzer {
public:
    explicit Endchecks(fhicl::ParameterSet const& p);
    Endchecks(Endchecks const&) = delete;
    Endchecks(Endchecks&&) = delete;
    Endchecks& operator=(Endchecks const&) = delete;
    Endchecks& operator=(Endchecks&&) = delete;

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

    enum EnumDet { kPDVD, kPDHD };
    EnumDet geoDet;

    // Detector Properties
    // float fADC2MeV;
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
        tag_spt, tag_pfp;

    // Input Parameters
    bool fLog;
    bool fKeepOutside;
    float fTrackLengthCut; // in cm
    float fMichelSpaceRadius; // in cm
    float fMichelTickRadius; // in ticks
    float fNearbySpaceRadius; // in cm
    float fBodyDistance; // in cm
    unsigned fRegNmin;
    float fBraggThreshold = 1.7; // in MIP dE/dx


    double GetSpace(geo::WireID);
    ana::Hit GetHit(HitPtr const p_hit);
    HitPtrPair GetTrackEndsHits( // return hits at track ends
        HitPtrVec const&, // all hits of the track
        HitPtrPair *pp_cathode_crossing = nullptr, // return hits at cathode crossing
        HitPtrVec *vp_tpc_crossing = nullptr, // return hits at TPC crossing (PDVD)
        HitPtrVec *vp_sorted_hit = nullptr, // return hits of the track per section and sorted
        geo::View_t = geo::kW // view to consider
    );
    std::pair<double, double> LinRegPCA(HitPtrVec const&);
    HitPtr GetBraggEnd( // return hit with max local dE/dx around track end
        HitPtrVec const&, // all hits in event
        art::Ptr<recob::Track> const, // track to check
        art::FindOneP<recob::Track> const& fop_hit2trk, // hits to track association
        HitPtr const trk_end, // track end as starting point
        std::vector<HitPtrVec> const& vvp_sec_sorted_hits, // hits of the track per section and sorted
        double *bragg_max = nullptr // return max of the local dE/dx 
    );
};


ana::Endchecks::Endchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fLog(p.get<bool>("Log", false)),
    fKeepOutside(p.get<bool>("KeepOutside", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 40.F)), // in cm
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius", 20.F)), //in cm
    fNearbySpaceRadius(p.get<float>("NearbySpaceRadius", 40.F)), //in cm
    fBodyDistance(p.get<float>("BodyDistance", 20.F)), // in cm
    fRegNmin(p.get<unsigned>("RegNmin", 6)),
    fBraggThreshold(p.get<float>("BraggThreshold", 1.7)) // in MIP dE/dx
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
        << "  Body Distance: " << fBodyDistance << " cm" << std::endl
        << "  Min. Reg. N: " << fRegNmin << std::endl;
}

void ana::Endchecks::analyze(art::Event const& e) {
    // auto const clockData = asDetClocks->DataFor(e);
    // auto const detProp = asDetProp->DataFor(e,clockData);
    // fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    // fDriftVelocity = detProp.DriftVelocity();
    // fTick2cm = fDriftVelocity * fSamplingRate;

    // auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    // if (!vh_hit.isValid()) return;
    // HitPtrVec vp_hit;
    // art::fill_ptr_vector(vp_hit, vh_hit);

    // auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    // if (!vh_trk.isValid()) return;
    // std::vector<art::Ptr<recob::Track>> vp_trk;
    // art::fill_ptr_vector(vp_trk, vh_trk);

    // art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    // art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);


    // std::vector<HitPtrVec> BraggTrackHits;

    // // loop over tracks to find muons
    // for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        // if (!LOG(p_trk->Length() > fTrackLengthCut)) continue;

        // simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);

        // // simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        // // tracks associated to a MCTruth muon
        // // if (!mcp) continue;
        // // if (!LOG(abs(mcp->PdgCode()) == 13)) continue;

        // HitPtrVec vp_hit_muon = fmp_trk2hit.at(p_trk.key());
        // if (!LOG(vp_hit_muon.size())) continue;

        // HitPtr mcp_end;
        // if (mcp) {
        //     HitPtrPair ends;
        //     ends = GetTrackEndsHits(ana::mcp2hits(
        //         mcp, vp_hit, clockData, false
        //     ));

        //     if (ends.first && ends.second) {
        //         int dir_z = mcp->EndZ() > mcp->Vz() ? 1 : -1;
        //         float fz = GetSpace(ends.first->WireID());
        //         float sz = GetSpace(ends.second->WireID());
        //         mcp_end = (sz-fz) * dir_z > 0 ? ends.second : ends.first;
        //     }
        // }

        // HitPtrPair trk_ends;
        // std::vector<HitPtrVec> vvp_sec_sorted_hits;
        // trk_ends = GetTrackEndsHits(vp_hit_muon, nullptr, nullptr, &vvp_sec_sorted_hits);

        // if (!LOG(trk_ends.first && trk_ends.second)) continue;

        // // fiducial cuts
        // std::pair<bool, bool> trk_isin;
        // trk_isin.first = wireWindow.isInside(trk_ends.first->PeakTime(), fMichelTickRadius)
        //     && geoHighX.InFiducialZ(GetSpace(trk_ends.first->WireID()), fMichelSpaceRadius);
        // trk_isin.second = wireWindow.isInside(trk_ends.second->PeakTime(), fMichelTickRadius)
        //     && geoHighX.InFiducialZ(GetSpace(trk_ends.second->WireID()), fMichelSpaceRadius);

        // if (!LOG(trk_isin.first || trk_isin.second)) continue;


        // HitPtrPair bragg_ends;
        // std::pair<double, double> braggs;
        // if (trk_isin.first)
        //     bragg_ends.first = GetBraggEnd(
        //         vp_hit,
        //         p_trk,
        //         fop_hit2trk,
        //         trk_ends.first,
        //         vvp_sec_sorted_hits,
        //         &braggs.first
        //     );
        // if (trk_isin.second)
        //     bragg_ends.second = GetBraggEnd(
        //         vp_hit,
        //         p_trk,
        //         fop_hit2trk,
        //         trk_ends.second,
        //         vvp_sec_sorted_hits,
        //         &braggs.second
        //     );


        // if (bragg_ends.first && braggs.first > fBraggThreshold) {
            
        // }



        // if (!LOG(fKeepOutside or (MuonEndIsInWindowT and MuonEndIsInVolumeYZ))) continue;

        // true test is broken
        // if (mcp) {
        //     std::vector<art::Ptr<recob::Track>> vp_trk_from_mcp = mcp2trks(mcp, vp_trk, clockData, fmp_trk2hit);
        //     // std::cout << "trk#" << p_trk->ID() << " mu#" << mcp->TrackId() << " encounter#" << ++particle_encounter[mcp->TrackId()] << " #mcp2trk:" << vp_trk_from_mcp.size();
        //     if (vp_trk_from_mcp.size()) {
        //         // bool isin = std::find(vp_trk_from_mcp.begin(), vp_trk_from_mcp.end(), p_trk) != vp_trk_from_mcp.end();
        //         // std::cout << " \033[1;9" << (isin ? 2 : 1) << "m" << (isin ? "in" : "not in") << "\033[0m";
        //         // std::cout << " [ ";
        //         // for (art::Ptr<recob::Track> p : vp_trk_from_mcp) std::cout << p->ID() << ", ";
        //         // std::cout << "]" << std::endl;
        //         if (vp_trk_from_mcp.size() == 1)
        //             MuonTrackIsNotBroken = kNotBroken;
        //         else {
        //             bool IsDeepestTrack = true;
        //             for (art::Ptr<recob::Track> p_trk_from_mcp : vp_trk_from_mcp)
        //                 if (geoDet == kPDVD)
        //                     IsDeepestTrack = IsDeepestTrack
        //                         && (MuonEndTrackPoint.x <= p_trk_from_mcp->Start().X()
        //                         && MuonEndTrackPoint.x <= p_trk_from_mcp->End().X());
        //                 else if (geoDet == kPDHD)
        //                     IsDeepestTrack = IsDeepestTrack
        //                         && (MuonEndTrackPoint.y <= p_trk_from_mcp->Start().Y()
        //                         && MuonEndTrackPoint.y <= p_trk_from_mcp->End().Y());
        //             if (IsDeepestTrack)
        //                 MuonTrackIsNotBroken = kLastOfBroken;
        //             else
        //                 MuonTrackIsNotBroken = kBroken;
        //         }
        //     } else MuonTrackIsNotBroken = -1;
        // } else MuonTrackIsNotBroken = -1;

        // if (mcp) {
        //     MuonEndProcess = mcp->EndProcess();
        // } else {
        //     MuonEndProcess = "";
        // }

        // a decaying muon has nu_mu, nu_e and elec as last daughters
        // simb::MCParticle const* mcp_michel = nullptr;
        // if (mcp) {
        //     bool has_numu = false, has_nue = false;
        //     for (int i_dau=mcp->NumberDaughters()-3; i_dau<mcp->NumberDaughters(); i_dau++) {
        //         simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    
        //         if (!mcp_dau) continue;

        //         switch (abs(mcp_dau->PdgCode())) {
        //             case 14: has_numu = true; break;
        //             case 12: has_nue = true; break;
        //             case 11: mcp_michel = mcp_dau; break;
        //             default: break;
        //         }
        //     }
        //     if (mcp_michel and has_numu and has_nue) {
        //         bool isin = false;
        //         for (unsigned t=0; t<asGeo->NTPC(); t++) {
        //             isin = asGeo->TPC(geo::TPCID{0, t}).ContainsPosition(mcp_michel->Position().Vect());
        //             if (isin) break;
        //         }
        //         if (isin)
        //             MuonHasMichel = kHasMichelInside; 
        //         else
        //             MuonHasMichel = kHasMichelOutside;

        //         // recob::Track const * trk_michel = truthUtil.GetRecoTrackFromMCParticle(clockData, *mcp_michel, e, tag_trk.label());
        //         art::Ptr<recob::Track> trk_michel = ana::mcp2trk(mcp_michel, vp_trk, clockData, fmp_trk2hit);
        //         if (trk_michel) 
        //             MichelTrackLength = trk_michel->Length();
        //         else
        //             MichelTrackLength = 0;

        //         MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3;

        //         // std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_michel, e, tag_hit.label());
        //         HitPtrVec vp_hit_michel = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);
        //         for (HitPtr p_hit_michel : vp_hit_michel) {
        //             if (p_hit_michel->View() != geo::kW) continue;

        //             MichelHits.push_back(GetHit(p_hit_michel));
        //         }
        //         MichelHitEnergy = MichelHits.energy();
        //     }
        //     else MuonHasMichel = kNoMichel;
        // } else {
        //     MuonHasMichel = -1;
        //     MichelTrackLength = -1;
        //     MichelHitEnergy = -1;
        // }
    // } // end of loop over tracks
}

void ana::Endchecks::beginJob() {}
void ana::Endchecks::endJob() {}

double ana::Endchecks::GetSpace(geo::WireID wid) {
    return plane2axis[(geo::PlaneID) wid].space(asWire->Wire(wid));
}

ana::Hit ana::Endchecks::GetHit(HitPtr const p_hit) {
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

HitPtrPair ana::Endchecks::GetTrackEndsHits(
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
    for (ana::LinearRegression& reg : side_reg) {
        if (reg.n < nmin) return {};
        reg.normalize();
    }

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
    auto side_sort = [&](int side) {
        return [&](HitPtr const& h1, HitPtr const& h2) -> bool {
            double const s1 = side_reg[side].projection(
                GetSpace(h1->WireID()),
                h1->PeakTime() * fTick2cm
            );
            double const s2 = side_reg[side].projection(
                GetSpace(h2->WireID()),
                h2->PeakTime() * fTick2cm
            );
            return s1 < s2;
        };
    };

    if (vp_sorted_hit) {
        vp_sorted_hit->clear();

        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            HitPtrVec& vp_sec_sorted = vp_sec_hit[s];
            if (vp_sec_sorted.size() < nmin) continue;

            int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
            std::sort(
                vp_sec_sorted.begin(), 
                vp_sec_sorted.end(),
                side_sort(side)
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
                side_sort(side)
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

std::pair<double, double> ana::Endchecks::LinRegPCA(HitPtrVec const& vp_hits) {
    unsigned n = vp_hits.size();
    if (!n) return {0, 0};

    double cov=0, varz=0, vart=0;
    for (HitPtr const& p_hit : vp_hits) {
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        cov += z * t;
        varz += z * z;
        vart += t * t;
    } 
    cov /= n;
    varz /= n;
    vart /= n;

    double corr = 1 - 4 * (varz*vart - cov*cov) / pow(varz + vart, 2);

    double lambda_plus = .5 * (varz + vart + sqrt(pow(varz-vart,2) + 4*cov*cov));

    int dirz = GetSpace(vp_hits.back()->WireID()) > GetSpace(vp_hits.front()->WireID()) ? 1 : -1;
    int dirt = vp_hits.back()->PeakTime() > vp_hits.front()->PeakTime() ? 1 : -1;
    double theta = atan2(
        dirt * abs(cov),
        dirz * abs(lambda_plus - vart)
    );

    return {corr, theta};
}

HitPtr ana::Endchecks::GetBraggEnd(
    HitPtrVec const& vp_hit,
    art::Ptr<recob::Track> const p_trk,
    art::FindOneP<recob::Track> const& fop_hit2trk,
    HitPtr const trk_end, 
    std::vector<HitPtrVec> const& vvp_sec_sorted_hits,
    double *bragg_max
) {
    int sec = ana::tpc2sec[geoDet][trk_end->WireID().TPC];

    HitPtrVec const &vp_prev_hits = vvp_sec_sorted_hits[sec];
    if (trk_end.key() == vvp_sec_sorted_hits[sec].front().key()) {
        // nothing to do
    } else if (trk_end.key() == vvp_sec_sorted_hits[sec].back().key()) {
        std::reverse(vvp_sec_sorted_hits[sec].begin(), vvp_sec_sorted_hits[sec].end());
    } else {
        LOG(!"Muon track start is not the first or last hit in its TPC section, this is unexpected.");
    }

    std::vector<float> distance;
    for (HitPtr const& p_hit : vp_prev_hits) {
        distance.push_back(
            sqrt(
                pow((p_hit->PeakTime() - vp_prev_hits.front()->PeakTime()) * fTick2cm, 2)
                + pow(GetSpace(p_hit->WireID()) - GetSpace(vp_prev_hits.front()->WireID()), 2)
            )
        );
    }
    // ptrdiff_t const body_start = std::distance(
    //     distance.begin(),
    //     std::find_if(
    //         distance.begin(), distance.end(),
    //         [&body_distance](float d) { return d > body_distance; }
    //     )
    // );
    // auto body_end = std::find_if(
    //     body_start, distance.end(),
    //     [&body_distance, &body_length](float d) { return d > body_distance + body_length; }
    // );

    auto dist2 = [&](HitPtr const& p_hit1, HitPtr const& p_hit2) {
        return pow((p_hit1->PeakTime() - p_hit2->PeakTime()) * fTick2cm, 2)
            + pow(GetSpace(p_hit1->WireID()) - GetSpace(p_hit2->WireID()), 2);
    };

    HitPtrVec::const_iterator body_start = std::find_if(
        vp_prev_hits.begin(), vp_prev_hits.end(),
        [&](HitPtr const& p_hit) {
            return dist2(p_hit, vp_prev_hits.front()) > fBodyDistance * fBodyDistance;
        }
    );


    if (std::distance(body_start, vp_prev_hits.end()) < fRegNmin) return;
    HitPtrVec reg_hits{body_start, body_start + fRegNmin};
    double corr, theta;
    std::tie(corr, theta) = LinRegPCA(reg_hits);
    double sigma = TMath::Pi() / 4 / corr;

    double z0 = GetSpace(vp_prev_hits.front()->WireID());
    double t0 = vp_prev_hits.front()->PeakTime() * fTick2cm;
    double sec0 = ana::tpc2sec[geoDet][vp_prev_hits.front()->WireID().TPC];

    HitPtrVec nearby_hits;
    for (HitPtr const& p_hit : vp_hit) {
        if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != sec0) continue;

        // check if the hit is in another track
        art::Ptr<recob::Track> p_hit_trk = fop_hit2trk.at(p_hit.key());
        if (
            p_hit_trk
            && p_hit_trk.key() != p_trk.key() 
            && p_hit_trk->Length() > fTrackLengthCut
        ) continue;

        // check if the hit is already treated
        if (std::find_if(
            body_start, vp_prev_hits.end(),
            [&](HitPtr const& p_prev_hit) {
                return p_prev_hit.key() == p_hit.key();
            }
        ) != vp_prev_hits.end()) continue;

        // check if the hit is close enough
        double d2 = dist2(p_hit, vp_prev_hits.front());
        if (d2 > 40) continue;

        nearby_hits.push_back(p_hit);
    }

    auto score = [](double z, double t, double z0, double t0, double theta, double sigma) -> double {
        double dz = z - z0;
        double dt = t - t0;
        double da = atan2(dt, dz) - theta;
        da = abs(da) < TMath::Pi() ? da : da - (da>0?1:-1) * 2 * TMath::Pi();
        double r = sqrt(dt*dt + dz*dz);

        return TMath::Gaus(da, 0, sigma) / r;
    };

    while (nearby_hits.size()) {
        HitPtrVec::iterator max_score_hit;
        double max_score = std::numeric_limits<double>::lowest();
        for (auto it = nearby_hits.begin(); it != nearby_hits.end(); it++) {
            double s = score(
                GetSpace((*it)->WireID()),
                (*it)->PeakTime() * fTick2cm,
                z0, t0, theta, sigma
            );
            if (s > max_score) {
                max_score = s;
                max_score_hit = it;
            }
        }

        nearby_hits.erase(max_score_hit);

        reg_hits.insert(reg_hits.begin(), *max_score_hit);
        reg_hits.pop_back();

        std::tie(corr, theta) = LinRegPCA(reg_hits);
        sigma = TMath::Pi() / 4 / corr;
    }

    std::vector<double> reg_dx;
    for (auto it = reg_hits.begin(); it != reg_hits.end(); it++) {
        if (it == reg_hits.begin()) {
            reg_dx.push_back(sqrt(dist2(*it, *(it+1))));
            continue;
        } 
        reg_dx.push_back( 0.5 * (
            sqrt(dist2(*it, *(it-1)))
            + sqrt(dist2(*it, *(it+1)))
        ));
    }

    unsigned const trailing_radius = 6;
    // std::vector<double> reg_local_dEdx;
    double max_local_dEdx = std::numeric_limits<double>::lowest();
    HitPtr max_local_dEdx_hit;
    for (unsigned i=0; i<reg_hits.size(); i++) {
        unsigned l = i+trailing_radius <= reg_hits.size() ? trailing_radius : reg_hits.size()-i;

        double local_dE = std::accumulate(reg_hits.begin()+i, reg_hits.begin()+i+l, 0.,
            [](double sum, HitPtr const& p_hit) {
                return sum + p_hit->Integral();
            }
        );
        local_dE /= l;

        double local_dx = std::accumulate(reg_dx.begin()+i, reg_dx.begin()+i+l, 0.);
        local_dx /= l;

        double local_dEdx = local_dE / local_dx;
        if (local_dEdx > max_local_dEdx) {
            max_local_dEdx = local_dEdx;
            max_local_dEdx_hit = reg_hits[i];
        }

        // reg_local_dEdx.push_back(local_dE / local_dx);
    }
    if (bragg_max) *bragg_max = max_local_dEdx;
    return max_local_dEdx_hit;
};




DEFINE_ART_MODULE(ana::Endchecks)