////////////////////////////////////////////////////////////////////////
// Class:       TrackDisplay
// Plugin Type: analyzer (Unknown Unknown)
// File:        TrackDisplay_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "event_display.h"

using HitPtr = art::Ptr<recob::Hit>;
using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

namespace ana {
    class TrackDisplay;
    struct MarkerStyle {
        Color_t c = kBlack;
        Style_t m = kFullCircle;
        Size_t s = 1.;
    };
    struct LineStyle {
        Color_t c = kBlack;
        Style_t l = kSolid;
        Width_t w = 1;
    };
    void inline setMarkerStyle(TAttMarker* m, MarkerStyle const& ms) {
        m->SetMarkerColor(ms.c);
        m->SetMarkerStyle(ms.m);
        m->SetMarkerSize(ms.s);
    };
    void inline setLineStyle(TAttLine* l, LineStyle const& ls) {
        l->SetLineColor(ls.c);
        l->SetLineStyle(ls.l);
        l->SetLineWidth(ls.w);
    };
}

class ana::TrackDisplay : public art::EDAnalyzer {
public:
    explicit TrackDisplay(fhicl::ParameterSet const& p);
    TrackDisplay(TrackDisplay const&) = delete;
    TrackDisplay(TrackDisplay&&) = delete;
    TrackDisplay& operator=(TrackDisplay const&) = delete;
    TrackDisplay& operator=(TrackDisplay&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:

    // Utilities
    art::ServiceHandle<art::TFileService> asFile;

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
        tag_hit, tag_clu, tag_trk, tag_shw,
        tag_spt, tag_pfp;

    // Input Parameters
    bool fLog;
    float fTrackLengthCut; // in cm

    unsigned ev=0;

    Int_t pal = kCividis;
    std::vector<Color_t> vc_pass = {kBlue, kBlue-3, kBlue+2, kAzure-2, kAzure+2, kAzure+7};
    std::vector<Color_t> vc_fail = {kRed, kRed-3, kRed+3, kPink-2, kPink-8, kPink+7};
    MarkerStyle
        ms_ev = {kBlack, 50, 0.5},
        ms_end = {kViolet+6, kFullSquare},
        ms_cc = {kViolet+6, kFullTriangleDown},
        ms_sc = {kViolet+6, kFullCircle},
        ms_pass = {vc_pass.front(), kFullCircle},
        ms_fail = {vc_fail.front(), kFullCircle},
        ms_back = {kGray, kFullCircle, 0.5};

    LineStyle
        ls_pass = {vc_pass.front(), kSolid, 2},
        ls_fail = {vc_fail.front(), kSolid, 2},
        ls_back = {kGray, kSolid, 1};

    bool IsUpright(recob::Track const& T);
    double GetSpace(geo::WireID);
    ana::Hit GetHit(HitPtr const p_hit);
    HitPtrVec GetSortedHits(
        HitPtrVec const& vp_hit,
        int dirz,
        HitPtrPair *pp_cathode_crossing = nullptr,
        HitPtrVec *vp_section_crossing = nullptr,
        std::vector<ana::LinearRegression> *p_side_reg = nullptr,
        geo::View_t view = geo::kW
    );
};


ana::TrackDisplay::TrackDisplay(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fLog(p.get<bool>("Log", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 40.F)) // in cm
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
        else if (type == "recob::Shower")           tag_shw = tag;
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
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl;
}

void ana::TrackDisplay::beginJob() {
    // TCanvas *c = asFile->make<TCanvas>(
    //     "legend", "legend",
    //     1300, 800
    // );
    // c->cd();

    // float top_margin = 0.1;
    // float indent = 0.05;
    // float title_indent = 0.04;
    // float label_indent = 0.1;
    // float line_height = 0.05;
    // float line_length = 0.03;
    // float ncol = 4;
    // Style_t font = 43;
    // Style_t font_bold = 63;
    // Float_t title_size = 20;
    // Float_t label_size = 16;

    // TText* title = new TText();
    // title->SetTextFont(font_bold);
    // title->SetTextSize(title_size);
    // TText* label = new TText();
    // label->SetTextFont(font);
    // label->SetTextSize(label_size);

    // title->DrawText(title_indent, 1 - top_margin, "Event Markers");

    // c->Write();
}

void ana::TrackDisplay::analyze(art::Event const& e) {
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

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    if (!vh_shw.isValid()) return;
    std::vector<art::Ptr<recob::Shower>> vp_shw;
    art::fill_ptr_vector(vp_shw, vh_shw);

    art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    auto drawMarker = [&](TCanvas* hc, HitPtr const& p_hit, MarkerStyle const& ms) -> void {
        TMarker *m = new TMarker();
        setMarkerStyle(m, ms);
        if (geoDet == kPDVD) {
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            hc->cd(s+1);
            m->DrawMarker(GetSpace(p_hit->WireID()), p_hit->PeakTime() * fTick2cm);
        } else if (geoDet == kPDHD) {
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) return;
            hc->cd(s+1);
            m->DrawMarker(p_hit->PeakTime() * fTick2cm, GetSpace(p_hit->WireID()));
        }
    };

    auto drawGraph = [&](TCanvas* hc, HitPtrVec const& vp_hit, char const* draw, MarkerStyle const& ms={}, LineStyle const& ls={}) -> void {
        unsigned static gn=0;
        std::vector<TGraph*> gs(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            gs[s] = new TGraph();
            gs[s]->SetEditable(kFALSE);
            gs[s]->SetName(Form("g%u_s%u", gn, s));
            setMarkerStyle(gs[s], ms);
            setLineStyle(gs[s], ls);
        }
        for (HitPtr p_hit : vp_hit) {
            if (p_hit->View() != geo::kW) continue;
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) continue;
            if (geoDet == kPDVD)
                gs[s]->AddPoint(GetSpace(p_hit->WireID()), p_hit->PeakTime() * fTick2cm);
            else if (geoDet == kPDHD)
                gs[s]->AddPoint(p_hit->PeakTime() * fTick2cm, GetSpace(p_hit->WireID()));
        }
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            if (!gs[s]->GetN()) continue;
            hc->cd(s+1);
            gs[s]->Draw(draw);
        }
        gn++;
    };

    auto drawGraph2D = [&](TCanvas* tc, art::Ptr<recob::Track> const& p_trk, MarkerStyle const& ms={}, LineStyle const& ls={}) -> void {
        TGraph2D* g = new TGraph2D();
        setMarkerStyle(g, ms);
        setLineStyle(g, ls);
        for (unsigned it=0; it<p_trk->NumberTrajectoryPoints(); it++) {
            if (!p_trk->HasValidPoint(it)) continue;
            geo::Point_t pt = p_trk->LocationAtPoint(it);
            if (geoDet == kPDVD)
                g->AddPoint(pt.Y(), pt.Z(), pt.X());
            else if (geoDet == kPDHD)
                g->AddPoint(pt.Z(), pt.X(), pt.Y());
        }
        tc->cd();
        g->Draw("same line");
    };

    std::vector<char const*> cuts = {
        "None",
        Form("TrackLength >= %.0f cm", fTrackLengthCut),
        "EndInVolume (20 cm)",
        "CathodeCrossing",
        "AnodeCrossing (20 cm)"
    };

    gStyle->SetPalette(pal);

    std::vector<TCanvas*> hcs;
    for (unsigned ihc=0; ihc<cuts.size(); ihc++) {
        TCanvas* hc = asFile->make<TCanvas>(
            Form("hc%ucut%u", ev, ihc),
            Form("Hits: run:%u, subrun:%u, event:%u, cut:%s", e.run(), e.subRun(), e.event(), cuts[ihc]),
            1300,800
        );
        ana::drawFrame(hc, int(geoDet), cuts[ihc], Form("%s R:%u SR:%u E:%u", (e.isRealData()?"Data":"Simulation"), e.run(), e.subRun(), e.event()));
        drawGraph(hc, vp_hit, "p", ms_ev);
        hcs.push_back(hc);
    }

    std::vector<TCanvas*> tcs;
    for (unsigned itc=0; itc<cuts.size(); itc++) {
        TCanvas* tc = asFile->make<TCanvas>(
            Form("tc%ucut%u", ev, itc),
            Form("Tracks: run:%u, subrun:%u, event:%u, cut:%s", e.run(), e.subRun(), e.event(), cuts[itc]),
            1300,800
        );
        TGraph2D* axes = nullptr;
        if (geoDet == kPDVD) {
            axes = new TGraph2D(2, new double[2]{-400,400}, new double[2]{-10,310}, new double[2]{-400,400});
            axes->SetTitle(";Y (cm);Z (cm);X (cm)");
        } else if (geoDet == kPDHD) {
            axes = new TGraph2D(2, new double[2]{-10,470}, new double[2]{-400,400}, new double[2]{-10,610});
            axes->SetTitle(";Z (cm);X (cm);Y (cm)");
        }
        tc->cd();
        axes->Draw("ap");
        tcs.push_back(tc);
    }
        

    unsigned im=0;
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        HitPtrVec vp_hit_muon = fmp_trk2hit.at(p_trk.key());
        ASSERT(vp_hit_muon.size())

        bool isUpright =  IsUpright(*p_trk);
        geo::Point_t Start = isUpright ? p_trk->Start() : p_trk->End();
        geo::Point_t End = isUpright ? p_trk->End() : p_trk->Start();
        HitPtrPair cathode_crossing;
        HitPtrVec section_crossing;
        HitPtrVec vp_hit_muon_sorted = GetSortedHits(
            vp_hit_muon, 
            End.Z() > Start.Z() ? 1 : -1,
            &cathode_crossing,
            &section_crossing
        );
        ASSERT(vp_hit_muon_sorted.size())

        im++;

        std::vector<TCanvas*>::iterator ihc = hcs.begin();
        std::vector<TCanvas*>::iterator itc = tcs.begin();
        auto filter = [&](bool cut) -> bool {
            if (cut) {
                drawGraph(*ihc, vp_hit_muon_sorted, "l", {}, LineStyle{vc_fail[im%vc_fail.size()], ls_fail.l, ls_fail.w});
                drawGraph2D(*itc, p_trk, MarkerStyle{vc_fail[im%vc_fail.size()], ms_fail.m, ms_fail.s});
                for (auto jhc=ihc+1; jhc<hcs.end(); jhc++)
                    drawGraph(*jhc, vp_hit_muon_sorted, "l", {}, ls_back);
                for (auto jtc=itc+1; jtc<tcs.end(); jtc++)
                    drawGraph2D(*jtc, p_trk, ms_back);
                return true;
            }
            drawGraph(*ihc, vp_hit_muon_sorted, "l", {}, LineStyle{vc_pass[im%vc_pass.size()], ls_pass.l, ls_pass.w});
            drawGraph2D(*itc, p_trk, MarkerStyle{vc_pass[im%vc_pass.size()], ms_pass.m, ms_pass.s});
            for (HitPtr const& p_hit_sc : section_crossing)
                drawMarker(*ihc, p_hit_sc, ms_sc);
            drawMarker(*ihc, vp_hit_muon_sorted.front(), ms_end);
            drawMarker(*ihc, vp_hit_muon_sorted.back(), ms_end);
            if (cathode_crossing.first) {
                drawMarker(*ihc, cathode_crossing.first, ms_cc);
                drawMarker(*ihc, cathode_crossing.second, ms_cc);
            }
            ihc++;
            itc++;
            return false;
        };

        filter(false);

        bool TagTrackLength = p_trk->Length() >= fTrackLengthCut;
        if (filter(!TagTrackLength)) continue;

        bool TagEndInVolume = geoHighX.InFiducialY(End.Y(), 20.)
            && geoHighX.InFiducialZ(End.Z(), 20.);
        if (filter(!TagEndInVolume)) continue;
        // MuonEndIsInWindowT = wireWindow.isInside(MuonEndHit.tick, fMichelTickRadius);

        bool TagCathodeCrossing = (
            geoLowX.ContainsPosition(Start)
            && geoHighX.ContainsPosition(End)
        ) || (
            geoHighX.ContainsPosition(Start)
            && geoLowX.ContainsPosition(End)
        );
        if (filter(!TagCathodeCrossing)) continue;

        bool TagAnodeCrossing = false;
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
        if (filter(!TagAnodeCrossing)) continue;
    }

    for (TCanvas* hc : hcs)
        hc->Write();
    for (TCanvas* tc : tcs)
        tc->Write();
    
    ev++;
}

void ana::TrackDisplay::endJob() {}

bool ana::TrackDisplay::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}

double ana::TrackDisplay::GetSpace(geo::WireID wid) {
    return plane2axis[wid].space(asWire->Wire(wid));
}

ana::Hit ana::TrackDisplay::GetHit(HitPtr const p_hit) {
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

HitPtrVec ana::TrackDisplay::GetSortedHits(
    HitPtrVec const& vp_hit,
    int dirz,
    HitPtrPair *pp_cathode_crossing,
    HitPtrVec *vp_section_crossing,
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
            reg.normalize();
    if (p_side_reg) p_side_reg = &side_reg;
    for (int side=0; side<2; side++)
        if (side_reg[side].n >= ana::LinearRegression::nmin)
            std::sort(
                side_hit[side].begin(),
                side_hit[side].end(),
                [&, &reg=side_reg[side]](
                    HitPtr const& h1, HitPtr const& h2
                ) -> bool {
                    // int const sec1 = ana::tpc2sec[geoDet][h1->WireID().TPC];
                    // int const sec2 = ana::tpc2sec[geoDet][h2->WireID().TPC];
                    // if (sec1 != sec2) return sec2 > sec1
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
    if (vp_section_crossing) {
        vp_section_crossing->clear();
        std::vector<HitPtrVec> vp_sec_hit(ana::n_sec[geoDet]);
        for (HitPtr const& p_hit : side_hit[side_pair.first]) {
            int sec = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (sec == -1) continue;
            vp_sec_hit[sec].push_back(p_hit);
        }
        for (HitPtr const& p_hit : side_hit[side_pair.second]) {
            int sec = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (sec == -1) continue;
            vp_sec_hit[sec].push_back(p_hit);
        }
        for (unsigned sec=0; sec<ana::n_sec[geoDet]; sec++) {
            if (vp_sec_hit[sec].front() != side_hit[ana::sec2side[geoDet][sec]].front())
                vp_section_crossing->push_back(vp_sec_hit[sec].front());
            if (vp_sec_hit[sec].back() != side_hit[ana::sec2side[geoDet][sec]].back())
                vp_section_crossing->push_back(vp_sec_hit[sec].back());
        }
        // vp_section_crossing->clear();
        // HitPtr& prev_hit = side_hit[side_pair.first].front();
        // for (HitPtr const& p_hit : side_hit[side_pair.first]) {
        //     vp_sorted_hit.push_back(p_hit);
        //     if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != ana::tpc2sec[geoDet][prev_hit->WireID().TPC]) {
        //         vp_section_crossing->push_back(prev_hit);
        //         vp_section_crossing->push_back(p_hit);
        //     }
        //     prev_hit = p_hit;
        // }
        // prev_hit = side_hit[side_pair.second].front();
        // for (HitPtr const& p_hit : side_hit[side_pair.second]) {
        //     vp_sorted_hit.push_back(p_hit);
        //     if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != ana::tpc2sec[geoDet][prev_hit->WireID().TPC]) {
        //         vp_section_crossing->push_back(prev_hit);
        //         vp_section_crossing->push_back(p_hit);
        //     }
        //     prev_hit = p_hit;
        // }
        // return vp_sorted_hit;
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

DEFINE_ART_MODULE(ana::TrackDisplay)