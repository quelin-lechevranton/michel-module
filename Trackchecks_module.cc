////////////////////////////////////////////////////////////////////////
// Class:       Trackchecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Trackchecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"
#include "event_display.h"

namespace ana {
    class Trackchecks;
}

class ana::Trackchecks : public art::EDAnalyzer {
public:
    explicit Trackchecks(fhicl::ParameterSet const& p);
    Trackchecks(Trackchecks const&) = delete;
    Trackchecks(Trackchecks&&) = delete;
    Trackchecks& operator=(Trackchecks const&) = delete;
    Trackchecks& operator=(Trackchecks&&) = delete;

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
        tag_spt, tag_pfp;

    // Input Parameters
    bool fLog;
    bool fKeepOutside;
    float fTrackLengthCut; // in cm
    float fMichelSpaceRadius; // in cm
    float fMichelTickRadius; // in ticks
    float fNearbySpaceRadius; // in cm
    float fCoincidenceWindow; // in ticks
    float fCoincidenceRadius; // in cm


    unsigned cn=0;

    double GetSpace(geo::WireID);
    ana::Hit GetHit(art::Ptr<recob::Hit> const p_hit);
    // art::Ptr<recob::Hit> GetDeepestHit(
    //     std::vector<art::Ptr<recob::Hit>> const&,
    //     bool increazing_z,
    //     geo::View_t = geo::kW
    // );
    std::vector<art::Ptr<recob::Hit>> GetEndsHits(
        std::vector<art::Ptr<recob::Hit>> const&,
        geo::View_t = geo::kW
    );
};


ana::Trackchecks::Trackchecks(fhicl::ParameterSet const& p)
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
}

void ana::Trackchecks::analyze(art::Event const& e) {
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

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);


    TCanvas *c = tfs->make<TCanvas>(
        Form("c%u", cn++),
        Form("run:%u, subrun:%u, event:%u", e.run(), e.subRun(), e.event()),
        1300,800
    );
    ana::drawFrame(c, int(geoDet), e.run(), e.subRun(), e.event());

    auto drawMarker = [this, c](TMarker* m, art::Ptr<recob::Hit> const& p_hit) -> void {
        if (geoDet == kPDVD) {
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            c->cd(s+1);
            m->DrawMarker(GetSpace(p_hit->WireID()), p_hit->PeakTime());
        } else if (geoDet == kPDHD) {
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) return;
            c->cd(s+1);
            m->DrawMarker(p_hit->PeakTime(), GetSpace(p_hit->WireID()));
        }
    };

    gStyle->SetPalette(kCividis);
    // TArrayI const& colors = TColor::GetPalette();
    // TMarker* m = new TMarker();
    // m->SetMarkerStyle(kFullCircle);
    // for (art::Ptr<recob::Hit> p_hit : vp_hit) {
    //     if (p_hit->View() != geo::kW) continue;
    //     float const x = std::min(p_hit->Integral() / (geoDet == kPDVD ? 200.F : 1000.F), 1.F);
    //     m->SetMarkerSize(2*x+0.1);
    //     m->SetMarkerColor(colors[int((colors.GetSize()-1)*x)]);

    //     drawMarker(m, p_hit);
    // }

    if (geoDet == kPDVD) {
        std::vector<TH2F*> h2(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            h2[s] = new TH2F(
                Form("h2_%u", s),
                "event hits",
                600, 0, 300,
                600, 0, 6000
            );
        }
        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            h2[s]->Fill(GetSpace(p_hit->WireID()), p_hit->PeakTime(), p_hit->Integral());
        }
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            c->cd(s+1);
            h2[s]->Draw("COL");
        }
    } else if (geoDet == kPDHD) {
        std::vector<TH2F*> h2(ana::n_sec[geoDet], nullptr);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            h2[s] = new TH2F(
                Form("h2_%u", s),
                "event hits",
                600, 0, 6000,
                600, 0, 464
            );
        }
        for (art::Ptr<recob::Hit> p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) continue;
            h2[s]->Fill(p_hit->PeakTime(), GetSpace(p_hit->WireID()), p_hit->Integral());
        }
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            c->cd(s+1);
            h2[s]->Draw("COL");
        }
    }


    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        // simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        // tracks associated to a MCTruth muon
        // if (!mcp) continue;
        // if (!LOG(abs(mcp->PdgCode()) == 13)) continue;

        std::vector<art::Ptr<recob::Hit>> vp_hit_muon = fmp_trk2hit.at(p_trk.key());
        if (!LOG(vp_hit_muon.size())) continue;

        // if (mcp) {
        //     art::Ptr<recob::Hit> deephit = GetDeepestHit(
        //         ana::mcp2hits(mcp, vp_hit, clockData, false),
        //         mcp->EndZ() > mcp->Vz()
        //     );
        //     if (deephit) MuonTrueEndHit = GetHit(deephit);
        // } else MuonTrueEndHit = ana::Hit{};

        std::vector<art::Ptr<recob::Hit>> trk_ends = GetEndsHits(vp_hit_muon);
        if (!LOG(trk_ends.size())) continue;

        TMarker* m = new TMarker();
        m->SetMarkerStyle(kFullTriangleUp);
        m->SetMarkerColor(kOrange+6);
        m->SetMarkerSize(2);
        // if (geoDet == kPDVD) {
        //     int s = ana::tpc2sec[geoDet][trk_ends.front()->WireID().TPC];
        //     c->cd(s+1);
        //     m->DrawMarker(GetSpace(trk_ends.front()->WireID()), trk_ends.front()->PeakTime());
        //     s = ana::tpc2sec[geoDet][trk_ends.back()->WireID().TPC];
        //     c->cd(s+1);
        //     m->DrawMarker(GetSpace(trk_ends.back()->WireID()), trk_ends.back()->PeakTime());
        // } else if (geoDet == kPDHD) {
        //     m->DrawMarker(trk_ends.front()->PeakTime(), GetSpace(trk_ends.front()->WireID()));
        //     m->DrawMarker(trk_ends.back()->PeakTime(), GetSpace(trk_ends.back()->WireID()));
        // }
        drawMarker(m, trk_ends.front());
        drawMarker(m, trk_ends.back());

        if (trk_ends.size() > 2) {
            m->SetMarkerStyle(kFullTriangleDown);
            m->SetMarkerColor(kOrange+2);
            m->SetMarkerSize(2);
            // if (geoDet == kPDVD) {
            //     m->DrawMarker(GetSpace(trk_ends[1]->WireID()), trk_ends[1]->PeakTime());
            //     m->DrawMarker(GetSpace(trk_ends[2]->WireID()), trk_ends[2]->PeakTime());
            // } else if (geoDet == kPDHD) {
            //     m->DrawMarker(trk_ends[1]->PeakTime(), GetSpace(trk_ends[1]->WireID()));
            //     m->DrawMarker(trk_ends[2]->PeakTime(), GetSpace(trk_ends[2]->WireID()));
            // }
            drawMarker(m, trk_ends[1]);
            drawMarker(m, trk_ends[2]);
        }

        // fiducial cuts
        // MuonEndIsInWindowT = wireWindow.isInside(MuonEndHit.tick, fMichelTickRadius);
        // MuonEndIsInVolumeYZ = geoHighX.InFiducialY(MuonEndTrackPoint.y, fMichelSpaceRadius)
        //     && geoHighX.InFiducialZ(MuonEndTrackPoint.z, fMichelSpaceRadius);

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
        //         std::vector<art::Ptr<recob::Hit>> vp_hit_michel = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);
        //         for (art::Ptr<recob::Hit> p_hit_michel : vp_hit_michel) {
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
    } // end of loop over tracks
}

void ana::Trackchecks::beginJob() {}
void ana::Trackchecks::endJob() {}

double ana::Trackchecks::GetSpace(geo::WireID wid) {
    return plane2axis[(geo::PlaneID) wid].space(asWire->Wire(wid));
}

ana::Hit ana::Trackchecks::GetHit(art::Ptr<recob::Hit> const p_hit) {
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

// art::Ptr<recob::Hit> ana::Trackchecks::GetDeepestHit(
//     std::vector<art::Ptr<recob::Hit>> const& vp_hit,
//     bool increasing_z,
//     geo::View_t view
// ) {
//     if (vp_hit.empty()) return art::Ptr<recob::Hit>{};

//     art::Ptr<recob::Hit> DeepestHit;
//     if (geoDet == kPDVD) {

//         // test if the muon goes into the bottom volume
//         bool in_bot = false;
//         for (art::Ptr<recob::Hit> p_hit : vp_hit) {
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
//         for (art::Ptr<recob::Hit> p_hit : vp_hit) {
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

//         double extrem_s = in_bot ?
//             std::numeric_limits<double>::max() // search min_s
//             : std::numeric_limits<double>::lowest(); // search max_s

//         for (art::Ptr<recob::Hit> p_hit : vp_hit) {
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
//         for (art::Ptr<recob::Hit> p_hit : vp_hit) {
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
//         for (art::Ptr<recob::Hit> p_hit : vp_hit) {
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

//         double extrem_s = increasing_z ?
//             std::numeric_limits<double>::lowest() // search max_s
//             : std::numeric_limits<double>::max(); // search min_s

//         for (art::Ptr<recob::Hit> p_hit : vp_hit) {
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


std::vector<art::Ptr<recob::Hit>> ana::Trackchecks::GetEndsHits(
    std::vector<art::Ptr<recob::Hit>> const& vp_hit,
    geo::View_t view
) {
    if (vp_hit.empty()) return {};

    std::function<int(geo::TPCID::TPCID_t)> cathodeSide =
        geoDet == kPDVD
      ? [](geo::TPCID::TPCID_t tpc) -> int { return int(tpc >= 8); }
      : [](geo::TPCID::TPCID_t tpc) -> int { return (tpc == 1 || tpc == 5) ? 1 : ((tpc == 2 || tpc == 6) ? 0 : -1); };


    // z = m*t + p
    struct LinearRegression {
        unsigned n=0;
        double mz=0, mt=0, mt2=0, mzt=0;
        void add(double z, double t) {
            mz+=z; mt+=t; mt2+=t*t; mzt+=z*t; n++;
        }
        void normalize() {
            mz/=n; mt/=n; mt2/=n; mzt/=n;
        }
        double cov() const { return mzt - mz*mt; }
        double vart() const { return mt2 - mt*mt; }
        double m() const { return n>3 ? cov()/vart() : 0; }
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
        std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> hits;
        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::lowest();
        void test(double s, art::Ptr<recob::Hit> const& h) {
            if (s < min) {
                min=s;
                hits.first = h;
            }
            if (s > max) {
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

    unsigned const nmin = 4;
    if (reg_side0.n < nmin && reg_side1.n < nmin)
        return {};
    else if (reg_side0.n < nmin)
        return { ends_side1.hits.first, ends_side1.hits.second };
    else if (reg_side1.n < nmin)
        return { ends_side0.hits.first, ends_side0.hits.second };

    auto d2 = [this](art::Ptr<recob::Hit> h1, art::Ptr<recob::Hit> h2) {
        double z1 = GetSpace(h1->WireID());
        double t1 = h1->PeakTime() * fTick2cm;
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
        case 0: return {
            ends_side0.hits.second,
            ends_side0.hits.first,
            ends_side1.hits.first,
            ends_side1.hits.second
        };
        case 1: return {
            ends_side0.hits.second,
            ends_side0.hits.first,
            ends_side1.hits.second,
            ends_side1.hits.first
        };
        case 2: return {
            ends_side0.hits.first,
            ends_side0.hits.second,
            ends_side1.hits.first,
            ends_side1.hits.second
        };
        case 3: return {
            ends_side0.hits.first,
            ends_side0.hits.second,
            ends_side1.hits.second,
            ends_side1.hits.first
        };
        default: break;
    }
    return {};
}

// ana::Hit ana::Trackchecks::GetTrueEndHit( std::vector<art::Ptr<recob::Hit>> const& vp_hit, ana::Point end_z) {
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

DEFINE_ART_MODULE(ana::Trackchecks)