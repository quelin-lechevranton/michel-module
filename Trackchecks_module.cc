////////////////////////////////////////////////////////////////////////
// Class:       Trackchecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Trackchecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

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
    float fTrackLengthCut; // in cm
    bool fSort;
    float fMichelSpaceRadius; // in cm
    float fMichelTickRadius; // in ticks


    unsigned cn=0;

    double GetSpace(geo::WireID);
    ana::Hit GetHit(art::Ptr<recob::Hit> const p_hit);
    // art::Ptr<recob::Hit> GetDeepestHit(
    //     std::vector<art::Ptr<recob::Hit>> const&,
    //     bool increazing_z,
    //     geo::View_t = geo::kW
    // );

    std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> GetTrackEndsHits(
        std::vector<art::Ptr<recob::Hit>> const&,
        std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> *pp_cathode_crossing = nullptr,
        std::vector<art::Ptr<recob::Hit>> *vp_tpc_crossing = nullptr,
        std::vector<std::vector<art::Ptr<recob::Hit>>> *vvp_sec_sorted_hits = nullptr,
        geo::View_t = geo::kW
    );

    /*
    std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> GetEndHits(
        std::vector<art::Ptr<recob::Hit>> const&,
        std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> *ppp_cathode_crossing,
        std::vector<art::Ptr<recob::Hit>> *pvp_tpc_crossing,
        geo::View_t = geo::kW
    );
    */
};


ana::Trackchecks::Trackchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fLog(p.get<bool>("Log", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 40.F)), // in cm
    fSort(p.get<bool>("Sort", false)), 
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius", 20.F)) //in cm
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
        << "  Sampling Rate: " << fSamplingRate << " µs/tick" << std::endl
        << "  Drift Velocity: " << fDriftVelocity << " cm/µs" << std::endl
        << "  Channel Pitch: " << fChannelPitch << " cm" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX.Min() << " -> " << geoHighX.Max() << std::endl
        << "  LowX Bounds: " << geoLowX.Min() << " -> " << geoLowX.Max() << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl
        << "  Michel Space Radius: " << fMichelSpaceRadius << " cm"
        << " (" << fMichelTickRadius << " ticks)" << std::endl;
        // << "  Nearby Space Radius: " << fNearbySpaceRadius << " cm" << std::endl
        // << "  Coincidence Window: " << fCoincidenceWindow << " ticks" << std::endl;
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
    ana::drawFrame(c, int(geoDet), e.run(), e.subRun(), e.event(), e.isRealData());

    auto drawMarker = [&](TMarker* m, art::Ptr<recob::Hit> const& p_hit) -> void {
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

    std::vector<TH2F*> h2s(ana::n_sec[geoDet]);
    if (geoDet == kPDVD) {
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            h2s[s] = new TH2F(
                Form("h2_%u", s),
                "event hits",
                600, 0, 300,
                600, 0, 6000
            );
        }
    } else if (geoDet == kPDHD) {
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            h2s[s] = new TH2F(
                Form("h2_%u", s),
                "event hits",
                600, 0, 6000,
                600, 0, 464
            );
        }
    }
    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;
        int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
        if (s == -1) continue;
        // int q = std::min(p_hit->Integral() / 1000???200, 1.F);
        int q = p_hit->Integral();
        if (geoDet == kPDVD)
            h2s[s]->Fill(GetSpace(p_hit->WireID()), p_hit->PeakTime(), q);
        else if (geoDet == kPDHD)
            h2s[s]->Fill(p_hit->PeakTime(), GetSpace(p_hit->WireID()), q);
    }
    double max=0;
    for (TH2F* h2 : h2s)
        max = h2->GetMaximum() > max ? h2->GetMaximum() : max;
    for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
        h2s[s]->SetMinimum(0);
        h2s[s]->SetMaximum(max);

        c->cd(s+1);
        TH2F* f = (TH2F*) gPad->FindObject(Form("f%u", s));
        f->SetMinimum(0);
        f->SetMaximum(max);
    }
    for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
        c->cd(s+1);
        if ((geoDet == kPDVD && s==3) || (geoDet == kPDHD && s==1))
            h2s[s]->Draw("same colz");
        else
            h2s[s]->Draw("same col");
    }

    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        bool tooSmall = p_trk->Length() < fTrackLengthCut;

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


        std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> cathode_crossing;
        std::vector<art::Ptr<recob::Hit>> tpc_crossing;
        std::vector<std::vector<art::Ptr<recob::Hit>>> vvp_sec_sorted;

        std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> trk_ends =
            GetTrackEndsHits(
                vp_hit_muon,
                &cathode_crossing,
                &tpc_crossing,
                fSort ? &vvp_sec_sorted : nullptr
            );

        // auto trk_ends = GetEndHits(
        //     vp_hit_muon,
        //     &cathode_crossing,
        //     &tpc_crossing,
        //     &
        // );


        if (!LOG(trk_ends.first && trk_ends.second)) {
            std::vector<TGraph*> gs(ana::n_sec[geoDet]);
            for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
                gs[s] = new TGraph();
                gs[s]->SetName(Form("g%u_%u", p_trk->ID(), s));
                gs[s]->SetTitle(Form("track %u, section %u", p_trk->ID(), s));
                gs[s]->SetLineWidth(1);
                gs[s]->SetLineColor(kGreen-8);
            }
            continue;
        }

        bool outsideFront = !wireWindow.isInside(trk_ends.first->PeakTime(), fMichelTickRadius)
            || !geoHighX.InFiducialZ(GetSpace(trk_ends.first->WireID()), fMichelSpaceRadius);
        bool outsideBack = !wireWindow.isInside(trk_ends.second->PeakTime(), fMichelTickRadius)
            || !geoHighX.InFiducialZ(GetSpace(trk_ends.second->WireID()), fMichelSpaceRadius);

        TMarker* m = new TMarker();
        if (cathode_crossing.first.isNull()) {
            m->SetMarkerColor(tooSmall ? kGray : kOrange+6);
            m->SetMarkerStyle(outsideFront ? kOpenSquare : kFullSquare);
            drawMarker(m, trk_ends.first);
            m->SetMarkerStyle(outsideBack ? kOpenSquare : kFullSquare);
            drawMarker(m, trk_ends.second);
        } else {
            m->SetMarkerColor(tooSmall ? kGray : kOrange+6);
            m->SetMarkerStyle(outsideFront ? kOpenTriangleUp : kFullTriangleUp);
            drawMarker(m, trk_ends.first);
            m->SetMarkerStyle(outsideBack ? kOpenTriangleUp : kFullTriangleUp);
            drawMarker(m, trk_ends.second);

            m->SetMarkerStyle(kFullTriangleDown);
            m->SetMarkerColor(tooSmall ? kGray : kPink-2);
            drawMarker(m, cathode_crossing.first);
            drawMarker(m, cathode_crossing.second);
        }
        if (tpc_crossing.size()) {
            m->SetMarkerStyle(kFullCircle);
            m->SetMarkerColor(tooSmall ? kGray : kViolet+6);
            for (art::Ptr<recob::Hit> p_hit : tpc_crossing)
                drawMarker(m, p_hit);
        }

        /*
        std::vector<TGraph*> gs(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            gs[s] = new TGraph();
            gs[s]->SetName(Form("g%u_%u", p_trk->ID(), s));
            gs[s]->SetTitle(Form("track %u, section %u", p_trk->ID(), s));
            gs[s]->SetLineWidth(1);
            gs[s]->SetLineColor(trk_ends.first.isNull() ?
                kGreen-8 : (tooSmall ? kGray : kOrange+6)
            );
        }
        for (art::Ptr<recob::Hit> p_hit : vp_hit_muon) {
            if (p_hit->View() != geo::kW) continue;
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) continue;
            if (geoDet == kPDVD)
                gs[s]->AddPoint(GetSpace(p_hit->WireID()), p_hit->PeakTime());
            else if (geoDet == kPDHD)
                gs[s]->AddPoint(p_hit->PeakTime(), GetSpace(p_hit->WireID()));
        }
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            c->cd(s+1);
            if (gs[s]->GetN()) gs[s]->Draw("same l");
        }
        */

        std::vector<TGraph*> gs(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            ana::LinearRegression reg;
            // LinearRegression const& reg = s >= ana::n_sec[geoDet]/2 ?
            //     per_side_reg.second : per_side_reg.first;
            // bool lin = reg.r2() > 0.5;

            TGraph* g = gs[s] = new TGraph();
            g->SetName(Form("g%u_%u", p_trk->ID(), s));
            g->SetTitle(Form("track %u, section %u", p_trk->ID(), s));
            g->SetLineWidth(1);
            g->SetLineColor(tooSmall ? kGray : kOrange+6);
            if (fSort) {
                for (art::Ptr<recob::Hit> const& p_hit : vvp_sec_sorted[s]) {
                    double const z = GetSpace(p_hit->WireID());
                    double const t = p_hit->PeakTime();
                    reg.add(z, t);
                    if (geoDet == kPDVD)
                        g->AddPoint(z, t);
                    else if (geoDet == kPDHD)
                        g->AddPoint(t, z);
                }
            } else {
                for (art::Ptr<recob::Hit> const& p_hit : vp_hit_muon) {
                    if (p_hit->View() != geo::kW) continue;
                    if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != int(s)) continue;
                    double const z = GetSpace(p_hit->WireID());
                    double const t = p_hit->PeakTime();
                    reg.add(z, t);
                    if (geoDet == kPDVD)
                        g->AddPoint(z, t);
                    else if (geoDet == kPDHD)
                        g->AddPoint(t, z);
                }
            }
            reg.normalize();

            c->cd(s+1);
            if (g->GetN()) g->Draw("same l");

            TF1* f = new TF1();
            if (geoDet == kPDVD) {
                f = new TF1(
                    Form("f%u_%u", p_trk->ID(), s),
                    "(x - [1]) / [0]",
                    reg.mz - 10,
                    reg.mz + 10
                );
            } else if (geoDet == kPDHD) {
                f = new TF1(
                    Form("f%u_%u", p_trk->ID(), s),
                    "[0]*x + [1]",
                    reg.mt - 50,
                    reg.mt + 50
                );
            }
            f->SetParameter(0, reg.m());
            f->SetParameter(1, reg.p());
            f->SetLineColor(reg.r2() > 0.5 ? kAzure-4 : kViolet+6);
            f->SetLineWidth(2);

            c->cd(s+1);
            f->Draw("same");
        }



        /***********************************
         * 
         *  Sort hits
         * 
         * **********************************

        // split volume at de cathode
        auto cathodeSide =
            geoDet == kPDVD
        ? [](geo::TPCID::TPCID_t tpc) -> int {
                return int(tpc >= 8);
            }
            // geoDet == kPDHD
        : [](geo::TPCID::TPCID_t tpc) -> int {
                return (tpc == 1 || tpc == 5) ? 0
                    : ((tpc == 2 || tpc == 6) ? 1 : -1);
            };

        using HitPtr = art::Ptr<recob::Hit>;
        using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
        using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

        std::pair<HitPtrVec, HitPtrVec> per_side_vph;
        for (HitPtr const& p_hit : vp_hit_muon) {
            if (p_hit->View() != geo::kW) continue;
            switch (cathodeSide(p_hit->WireID().TPC)) {
                case 0: per_side_vph.first.push_back(p_hit); break;
                case 1: per_side_vph.second.push_back(p_hit); break;
                default: break;
            }
        }

        // unsure there is enough hits on at least one side
        unsigned const nmin = 4;
        if (!LOG(per_side_vph.first.size() >= nmin || per_side_vph.second.size() >= nmin)) continue;


        // linear regression on each side to have a curvilinear coordinate of each hit inside a track
        // z = m*t + p
        struct LinearRegression {
            unsigned n=0;
            double mz=0, mt=0, mz2=0, mt2=0, mzt=0;
            void add(double z, double t) {
                mz+=z; mt+=t; mz2+=z*z; mt2+=t*t; mzt+=z*t; n++;
            }
            void normalize() {
                mz/=n; mt/=n; mz2/=n; mt2/=n; mzt/=n;
            }
            double cov() const { return mzt - mz*mt; }
            double varz() const { return mz2 - mz*mz; }
            double vart() const { return mt2 - mt*mt; }
            double m() const { return n<nmin ? 0 : cov()/vart(); }
            double p() const { return mz - m()*mt; }
            double r2() const { return n<nmin ? 0 : cov()*cov() / (varz()*vart()); }
            double projection(double z, double t) const {
                return (t + m()*(z-p())) / (1 + m()*m());
            }
        };
        std::pair<LinearRegression, LinearRegression> per_side_reg;

        for (HitPtr const& p_hit : per_side_vph.first) {
            double const z = GetSpace(p_hit->WireID());
            double const t = p_hit->PeakTime() * fTick2cm;
            per_side_reg.first.add(z, t);
        }
        for (HitPtr const& p_hit : per_side_vph.second) {
            double const z = GetSpace(p_hit->WireID());
            double const t = p_hit->PeakTime() * fTick2cm;
            per_side_reg.second.add(z, t);
        }
        per_side_reg.first.normalize();
        per_side_reg.second.normalize();

        LOG("regression done");

        // compare hits by their curvilinear coordinate
        auto comp = [&](LinearRegression const& reg, HitPtr const& h1, HitPtr const& h2) -> bool {
            double const s1 = reg.projection(
                GetSpace(h1->WireID()),
                h1->PeakTime() * fTick2cm
            );
            double const s2 = reg.projection(
                GetSpace(h2->WireID()),
                h2->PeakTime() * fTick2cm
            );
            return s1 < s2;
        };

        // get the track ends on each side
        auto minmax = [&](HitPtrVec const& vph, LinearRegression const& reg) -> HitPtrPair {
            if (vph.size() < nmin) return {};
            auto mm = std::minmax_element(
                vph.begin(), vph.end(),
                [&comp, &reg](HitPtr const& h1, HitPtr const& h2) {
                    return comp(reg, h1, h2);
                }
            );
            return { *mm.first, *mm.second };
        };
        std::pair<HitPtrPair, HitPtrPair> per_side_ends = {
            minmax(per_side_vph.first, per_side_reg.first),
            minmax(per_side_vph.second, per_side_reg.second)
        };

        LOG("side ends done");

        // get a sorted list of hits for each section (ie. pair of TPCs)
        std::vector<HitPtrVec> per_sec_vph(ana::n_sec[geoDet]);
        for (HitPtr const& ph : per_side_vph.first)
            per_sec_vph[ana::tpc2sec[geoDet][ph->WireID().TPC]].push_back(ph);
        for (HitPtr const& ph : per_side_vph.second)
            per_sec_vph[ana::tpc2sec[geoDet][ph->WireID().TPC]].push_back(ph);

        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            LinearRegression const& reg = s >= ana::n_sec[geoDet]/2 ? per_side_reg.second : per_side_reg.first;
            std::sort(
                per_sec_vph[s].begin(), per_sec_vph[s].end(),
                [&comp, &reg](HitPtr const& h1, HitPtr const& h2) -> bool {
                    return comp(reg, h1, h2);
                }
            );
        }

        LOG("sec sort done");

        // get the track ends for each section
        std::vector<HitPtrPair> per_sec_ends(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            if (per_sec_vph[s].size() < nmin) continue;
            per_sec_ends[s].first = per_sec_vph[s].front();
            per_sec_ends[s].second = per_sec_vph[s].back();
        }

        LOG("sec ends done");

        // given the ends of two pieces of track, find the closest ends
        auto closestHits = [&](HitPtrPair const& hp1, HitPtrPair const& hp2, double dmin, HitPtrPair *outermostHits = nullptr) -> HitPtrPair {
            std::vector<HitPtrPair> pairs = {
                { hp1.first, hp2.first },
                { hp1.first, hp2.second },
                { hp1.second, hp2.second },
                { hp1.second, hp2.first }
            };
            std::vector<double> d2s(4, 0);
            for (unsigned i=0; i<4; i++) {
                double zf = GetSpace(pairs[i].first->WireID());
                double tf = pairs[i].first->PeakTime() * fTick2cm;
                double zs = GetSpace(pairs[i].second->WireID());
                double ts = pairs[i].second->PeakTime() * fTick2cm;
                d2s[i] = pow(zf-zs,2) + pow(tf-ts,2);
            }

            // find all distances under d2min threshold
            std::vector<unsigned> candidates_idx;
            std::vector<double>::iterator it = d2s.begin();
            while ((it = std::find_if(it, d2s.end(), [dmin](double d2) { return d2 < dmin*dmin; })) != d2s.end())
                candidates_idx.push_back(std::distance(d2s.begin(), it++));
            
            // no candidates found
            if (candidates_idx.empty())
                return {};

            // get the closest pair
            unsigned closest_idx = *std::min_element(candidates_idx.begin(), candidates_idx.end(),
                [&d2s](unsigned i, unsigned j) { return d2s[i] < d2s[j]; });
            if (outermostHits) {
                unsigned outermost_idx = (closest_idx+2) % 4;
                outermostHits->first = pairs[outermost_idx].first;
                outermostHits->second = pairs[outermost_idx].second;
            }
            return pairs[closest_idx];
        };

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

        LOG("tpc crossing done");

        // get the hits that are at the cathode, if any
        HitPtrPair trk_ends, cathode_crossing;
        if (per_side_vph.first.size() < nmin)
            trk_ends = per_side_ends.second;
        else if (per_side_vph.second.size() < nmin)
            trk_ends = per_side_ends.first;
        else
            cathode_crossing = closestHits(per_side_ends.first, per_side_ends.second, 2*fCathodeGap, &trk_ends);
        
        if (!LOG(trk_ends.first.isNonnull())) continue;

        // check if the track ends are outside the fiducial volume
        bool outsideFront = !wireWindow.isInside(trk_ends.first->PeakTime(), fMichelTickRadius)
            || !geoHighX.InFiducialZ(GetSpace(trk_ends.first->WireID()), fMichelSpaceRadius);
        bool outsideBack = !wireWindow.isInside(trk_ends.second->PeakTime(), fMichelTickRadius)
            || !geoHighX.InFiducialZ(GetSpace(trk_ends.second->WireID()), fMichelSpaceRadius);


        // PLOT

        TMarker* m = new TMarker();
        if (cathode_crossing.first.isNull()) {
            m->SetMarkerColor(tooSmall ? kGray : kOrange+6);
            m->SetMarkerStyle(outsideFront ? kOpenSquare : kFullSquare);
            drawMarker(m, trk_ends.first);
            m->SetMarkerStyle(outsideBack ? kOpenSquare : kFullSquare);
            drawMarker(m, trk_ends.second);
        } else {
            m->SetMarkerColor(tooSmall ? kGray : kOrange+6);
            m->SetMarkerStyle(outsideFront ? kOpenTriangleUp : kFullTriangleUp);
            drawMarker(m, trk_ends.first);
            m->SetMarkerStyle(outsideBack ? kOpenTriangleUp : kFullTriangleUp);
            drawMarker(m, trk_ends.second);

            m->SetMarkerStyle(kFullTriangleDown);
            m->SetMarkerColor(tooSmall ? kGray : kPink-2);
            drawMarker(m, cathode_crossing.first);
            drawMarker(m, cathode_crossing.second);
        }
        if (tpc_crossing.size()) {
            m->SetMarkerStyle(kFullCircle);
            m->SetMarkerColor(tooSmall ? kGray : kPink-2);
            for (HitPtr p_hit : tpc_crossing)
                drawMarker(m, p_hit);
        }

        std::vector<TGraph*> gs(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            LinearRegression const& reg = s >= ana::n_sec[geoDet]/2 ?
                per_side_reg.second : per_side_reg.first;
            bool lin = reg.r2() > 0.5;

            TGraph* g = gs[s] = new TGraph();
            g->SetName(Form("g%u_%u", p_trk->ID(), s));
            g->SetTitle(Form("track %u, section %u", p_trk->ID(), s));
            g->SetLineWidth(1);
            g->SetLineColor(
                tooSmall ? kGray : (
                    lin ? kOrange+6 : kViolet+6
                )
            );
            for (HitPtr const& p_hit : per_sec_vph[s]) {
                if (geoDet == kPDVD)
                    g->AddPoint(GetSpace(p_hit->WireID()), p_hit->PeakTime());
                else if (geoDet == kPDHD)
                    g->AddPoint(p_hit->PeakTime(), GetSpace(p_hit->WireID()));
            }

            c->cd(s+1);
            if (g->GetN()) g->Draw("same l");

            TF1* f = new TF1();
            if (geoDet == kPDVD) {
                f = new TF1(
                    Form("f%u_%u", p_trk->ID(), s),
                    "(x - [1]) / [0]",
                    reg.mz - 10,
                    reg.mz + 10
                );
            } else if (geoDet == kPDHD) {
                f = new TF1(
                    Form("f%u_%u", p_trk->ID(), s),
                    "[0]*x + [1]",
                    (reg.mt - 10) / fTick2cm,
                    (reg.mt + 10) / fTick2cm
                );
            }
            f->SetParameter(0, reg.m() / fTick2cm);
            f->SetParameter(1, reg.p());
            f->SetLineColor(kAzure-4);
            f->SetLineWidth(2);

            c->cd(s+1);
            f->Draw("same");
        }

        **********************************
         * 
         *  Sort hits
         * 
         * ***********************************/



                    





























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

    c->Write();
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

/*
std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> ana::Trackchecks::GetEndHits(
    std::vector<art::Ptr<recob::Hit>> const& vp_hit,
    std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> *ppp_cathode_crossing,
    std::vector<art::Ptr<recob::Hit>> *pvp_tpc_crossing,
    geo::View_t view
) {
    unsigned const nmin = 4;
    if (vp_hit.size() < nmin) return {};

    std::function<int(geo::TPCID::TPCID_t)> cathodeSide =
        geoDet == kPDVD
      ? [](geo::TPCID::TPCID_t tpc) -> int { return int(tpc >= 8); }
      : [](geo::TPCID::TPCID_t tpc) -> int {
            return (tpc == 1 || tpc == 5) ? 1
                : ((tpc == 2 || tpc == 6) ? 0 : -1);
        };

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
        double m() const { return n<nmin ? 0 : cov()/vart(); }
        double p() const { return mz - m()*mt; }
        double projection(double z, double t) const {
            return (t + m()*(z-p())) / (1 + m()*m());
        }
    } side0_reg, side1_reg;

    // std::vector<std::vector<art::Ptr<recob::Hit>>> vvp_sec_hit(ana::n_sec[geoDet]);
    // for (art::Ptr<recob::Hit> p_hit : vp_hit) {
    //     if (p_hit->View() != view) continue;
    //     geo::WireID wireid = p_hit->WireID();
    //     int side = cathodeSide(wireid.TPC);
    //     if (side == -1) continue;
    //     vvp_sec_hit[ana::tpc2sec[geoDet][wireid.TPC]].push_back(p_hit);
    // }

    // std::vector<std::vector<unsigned>> vvi_sec_sort_idx(ana::n_sec[geoDet]);
    // for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
    //     std::vector<unsigned> &vi = vvi_sec_sort_idx[s];
    //     vi.resize(vvp_sec_hit[s].size());
    //     std::iota(vi.begin(), vi.end(), 0);

    //     LinearRegression const& reg = s >= ana::n_sec[geoDet]/2 ? side1_reg : side0_reg;
    //     std::sort(vi.begin(), vi.end(),
    //         [this, s, &reg, &vvp_sec_hit](unsigned i, unsigned j) {
    //             double si = reg.projection(
    //                 GetSpace(vvp_sec_hit[s][i]->WireID()),
    //                 vvp_sec_hit[s][i]->PeakTime() * fTick2cm
    //             );
    //             double sj = reg.projection(
    //                 GetSpace(vvp_sec_hit[s][j]->WireID()),
    //                 vvp_sec_hit[s][j]->PeakTime() * fTick2cm
    //             );
    //             return si < sj;
    //         }
    //     );
    // }

    std::vector<unsigned> sec_nhit(ana::n_sec[geoDet], 0);

    // unsigned nhit=0;
    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        geo::WireID wireid = p_hit->WireID();
        int side = cathodeSide(wireid.TPC);
        if (side == -1) continue;
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        if (side == 0)
            side0_reg.add(z, t);
        else if (side == 1)
            side1_reg.add(z, t);

        if (geoDet == kPDVD) // for tpc crossing
            sec_nhit[ana::tpc2sec[geoDet][wireid.TPC]]++;
    }
    side0_reg.normalize();
    side1_reg.normalize();

    if (side0_reg.n < nmin && side1_reg.n < nmin) return {};

    using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;
    struct ProjectionEnds {
        HitPtrPair hits;
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
    };
    ProjectionEnds side0_ends, side1_ends;
    std::vector<ProjectionEnds> sec_ends(ana::n_sec[geoDet]);

    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        geo::WireID wireid = p_hit->WireID();
        int side = cathodeSide(wireid.TPC);
        if (side == -1) continue;
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        if (side == 0) {
            double s = side0_reg.projection(z, t);
            side0_ends.test(s, p_hit);
            
            if (geoDet == kPDVD) // for tpc crossing
                sec_ends[ana::tpc2sec[geoDet][wireid.TPC]].test(s, p_hit);

        } else if (side == 1) {
            double s = side1_reg.projection(z, t);
            side1_ends.test(s, p_hit);

            if (geoDet == kPDVD) // for tpc crossing
                sec_ends[ana::tpc2sec[geoDet][wireid.TPC]].test(s, p_hit);
        }
    }

    // auto d2 = [&](art::Ptr<recob::Hit> h1, art::Ptr<recob::Hit> h2) -> double {
    //     double z1 = GetSpace(h1->WireID());
    //     double t1 = h1->PeakTime() * fTick2cm;
    //     double z2 = GetSpace(h2->WireID());
    //     double t2 = h2->PeakTime() * fTick2cm;
    //     return pow(z1-z2,2) + pow(t1-t2,2);
    // };

    auto closestHits = [&](HitPtrPair const& hp0, HitPtrPair const& hp1, double dmin, HitPtrPair *outermostHits = nullptr) -> HitPtrPair {
        std::vector<HitPtrPair> pairs = {
            { hp0.first, hp1.first },
            { hp0.first, hp1.second },
            { hp0.second, hp1.second },
            { hp0.second, hp1.first }
        };
        std::vector<double> d2s(4, 0);
        for (unsigned i=0; i<4; i++) {
            double zf = GetSpace(pairs[i].first->WireID());
            double tf = pairs[i].first->PeakTime() * fTick2cm;
            double zs = GetSpace(pairs[i].second->WireID());
            double ts = pairs[i].second->PeakTime() * fTick2cm;
            d2s[i] = pow(zf-zs,2) + pow(tf-ts,2);
        }

        // find all distances under d2min threshold
        std::vector<unsigned> candidates_idx;
        std::vector<double>::iterator it = d2s.begin();
        while ((it = std::find_if(it, d2s.end(), [dmin](double d2) { return d2 < dmin*dmin; })) != d2s.end())
            candidates_idx.push_back(std::distance(d2s.begin(), it++));
        
        // no candidates found
        if (candidates_idx.empty())
            return {};

        // get the closest pair
        unsigned closest_idx = *std::min_element(candidates_idx.begin(), candidates_idx.end(),
            [&d2s](unsigned i, unsigned j) { return d2s[i] < d2s[j]; });
        if (outermostHits) {
            unsigned outermost_idx = (closest_idx+2) % 4;
            outermostHits->first = pairs[outermost_idx].first;
            outermostHits->second = pairs[outermost_idx].second;
        }
        return pairs[closest_idx];

        // if (candidates_idx.size() == 1) {
        //     closest_idx = candidates_idx.front();
            // if (outermostHits) {
            //     unsigned outer_idx = (idx+2) % 4;
            //     outermostHits->first = pairs[outer_idx].first;
            //     outermostHits->second = pairs[outer_idx].second;
            // }
            // return pairs[idx];
        // }
        // if (candidates_idx.size() > 1) {
            // find the pair with the smallest distance
            // unsigned min_idx = *std::min_element(candidates_idx.begin(), candidates_idx.end(),
            //     [&d2s](unsigned i1, unsigned i2) { return d2s[i1] < d2s[i2]; });
            // if (outermostHits) {
            //     unsigned outer_idx = (min_idx+2) % 4;
            //     outermostHits->first = pairs[outer_idx].first;
            //     outermostHits->second = pairs[outer_idx].second;
            // }
            // return pairs[min_idx];
        // }
            
        // switch(std::distance(d2s.begin(), std::min_element(d2s.begin(), d2s.end()))) {
        //     case 0: 
        //         if (outermostHits) {
        //             outermostHits->first = hp0.second;
        //             outermostHits->second = hp1.second;
        //         }
        //         return { hp0.first, hp1.first };
        //     case 1:
        //         if (outermostHits) {
        //             outermostHits->first = hp0.second;
        //             outermostHits->second = hp1.first;
        //         }
        //         return { hp0.first, hp1.second };
        //     case 2:
        //         if (outermostHits) {
        //             outermostHits->first = hp0.first;
        //             outermostHits->second = hp1.second;
        //         }
        //         return { hp0.second, hp1.first };
        //     case 3:
        //         if (outermostHits) {
        //             outermostHits->first = hp0.first;
        //             outermostHits->second = hp1.first;
        //         }
        //         return { hp0.second, hp1.second };
        //     default: break;
        // }
        // return {};
    };

    // tpc crossing
    if (geoDet == kPDVD) {
        bool prev = false; // is there hits in the previous section?
        for (unsigned s=0; s<8; s++) {
            if (sec_nhit[s] < nmin) {
                prev = false;
                continue;
            }
            if (prev) {
                HitPtrPair const pp_tpc_crossing = closestHits(
                    sec_ends[s-1].hits, sec_ends[s].hits, 2
                );
                if (pp_tpc_crossing.first.isNonnull() && pvp_tpc_crossing) {
                    pvp_tpc_crossing->push_back(pp_tpc_crossing.first);
                    pvp_tpc_crossing->push_back(pp_tpc_crossing.second);
                }
            }
            if (s == 3)
                prev = false; // cathode crossing
            else
                prev = true;
        }
    }

    // no cathode crossing
    if (side0_reg.n < nmin)
        return side1_ends.hits;
    else if (side1_reg.n < nmin)
        return side0_ends.hits;

    // cathode crossing
    HitPtrPair ends;
    HitPtrPair const pp_cathode_crossing = closestHits(
        side0_ends.hits, side1_ends.hits, 2*fCathodeGap, &ends
    );

    if (pp_cathode_crossing.first.isNonnull() && ppp_cathode_crossing) {
        ppp_cathode_crossing->first = pp_cathode_crossing.first;
        ppp_cathode_crossing->second = pp_cathode_crossing.second;
    }
    return ends;

    // std::vector<double> d2s(4, 0);
    // d2s[0] = d2(side0_ends.hits.first, side1_ends.hits.first);
    // d2s[1] = d2(side0_ends.hits.first, side1_ends.hits.second);
    // d2s[2] = d2(side0_ends.hits.second, side1_ends.hits.first);
    // d2s[3] = d2(side0_ends.hits.second, side1_ends.hits.second);

    // switch(std::distance(d2s.begin(), std::min_element(d2s.begin(), d2s.end()))) {
    //     case 0: 
    //         ppp_cathode_crossing->first = side0_ends.hits.first;
    //         ppp_cathode_crossing->second = side1_ends.hits.first;
    //         return {
    //             side0_ends.hits.second,
    //             // side0_ends.hits.first,
    //             // side1_ends.hits.first,
    //             side1_ends.hits.second
    //         };
    //     case 1: 
    //         ppp_cathode_crossing->first = side0_ends.hits.first;
    //         ppp_cathode_crossing->second = side1_ends.hits.second;
    //         return {
    //             side0_ends.hits.second,
    //             // side0_ends.hits.first,
    //             // side1_ends.hits.second,
    //             side1_ends.hits.first
    //         };
    //     case 2: 
    //         ppp_cathode_crossing->first = side0_ends.hits.second;
    //         ppp_cathode_crossing->second = side1_ends.hits.first;
    //         return {
    //             side0_ends.hits.first,
    //             // side0_ends.hits.second,
    //             // side1_ends.hits.first,
    //             side1_ends.hits.second
    //         };
    //     case 3:
    //         ppp_cathode_crossing->first = side0_ends.hits.second;
    //         ppp_cathode_crossing->second = side1_ends.hits.second;
    //         return {
    //             side0_ends.hits.first,
    //             // side0_ends.hits.second,
    //             // side1_ends.hits.second,
    //             side1_ends.hits.first
    //         };
    //     default: break;
    // }
    // return {};
}
    
*/

std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> ana::Trackchecks::GetTrackEndsHits(
    std::vector<art::Ptr<recob::Hit>> const& vp_hit,
    std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>> *pp_cathode_crossing,
    std::vector<art::Ptr<recob::Hit>> *vp_tpc_crossing,
    std::vector<std::vector<art::Ptr<recob::Hit>>> *vvp_sec_sorted_hits,
    geo::View_t view
) {
    using HitPtr = art::Ptr<recob::Hit>;
    using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
    using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

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
    for (ana::LinearRegression& reg : side_reg) reg.normalize();

    // if not enough hits on both sides, return empty pair
    if (std::all_of(
            side_reg.begin(),
            side_reg.end(),
            [](const ana::LinearRegression& reg) { return reg.n < nmin; }
        )) return {};

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
    if (!vp_tpc_crossing && !vvp_sec_sorted_hits) {
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
    if (vvp_sec_sorted_hits) {

        // get a sorted list of hits for each section (ie. pair of TPCs)
        vvp_sec_sorted_hits->clear();
        vvp_sec_sorted_hits->resize(ana::n_sec[geoDet]);
        for (HitPtr const& p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) continue;
            vvp_sec_sorted_hits->at(s).push_back(p_hit);
        }

        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
            std::sort(
                vvp_sec_sorted_hits->at(s).begin(), 
                vvp_sec_sorted_hits->at(s).end(),
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
        }

        // get the track ends for each section
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            if (vvp_sec_sorted_hits->at(s).size() < nmin) continue;
            per_sec_ends[s].first = vvp_sec_sorted_hits->at(s).front();
            per_sec_ends[s].second = vvp_sec_sorted_hits->at(s).back();
        }

    } else { // only get the minmax ends of each section

        std::vector<HitPtrVec> per_sec_vph(ana::n_sec[geoDet]);
        for (HitPtr const& p_hit : vp_hit) {
            if (p_hit->View() != view) continue;
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) continue;
            per_sec_vph[s].push_back(p_hit);
        }

        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            if (per_sec_vph[s].size() < nmin) continue;
            int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
            auto minmax = std::minmax_element(
                per_sec_vph[s].begin(), per_sec_vph[s].end(),
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

DEFINE_ART_MODULE(ana::Trackchecks)