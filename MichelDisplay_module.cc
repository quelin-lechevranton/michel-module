#include "event_display.h"

namespace ana {
    class MichelDisplay;
}

class ana::MichelDisplay : public art::EDAnalyzer, private ana::MichelDisplayer {
public:
    explicit MichelDisplay(fhicl::ParameterSet const& p);
    MichelDisplay(MichelDisplay const&) = delete;
    MichelDisplay(MichelDisplay&&) = delete;
    MichelDisplay& operator=(MichelDisplay const&) = delete;
    MichelDisplay& operator=(MichelDisplay&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    ana::Bounds<float> wireWindow;
    ana::Bounds3D<float> geoHighX, geoLowX;
    float fCathodeGap; // cm

    // Input Parameters
    bool fLog;
    float fTrackLengthCut; // in cm
    float fNearbyRadius; // cm
    float fBodyDistance; // cm
    unsigned fRegN;
    float fBraggThreshold; // in MIP dE/dx

    unsigned ev=0;

    Int_t pal = kCividis;
    std::vector<Color_t> vc_pass = {kBlue, kBlue-3, kBlue+2, kAzure-2, kAzure+2, kAzure+7};
    std::vector<Color_t> vc_fail = {kRed, kRed-3, kRed+3, kPink-2, kPink-8, kPink+7};
    ana::MarkerStyle
        ms_ev = {kGray, kPlus, 0.5},
        ms_end = {kViolet+6, kFullSquare},
        ms_cc = {kViolet+6, kFullTriangleUp},
        ms_sc = {kViolet+6, kFullCircle},
        ms_pass = {vc_pass.front(), kFullCircle},
        ms_fail = {vc_fail.front(), kFullCircle},
        ms_back = {kGray, kFullCircle, 0.5},
        ms_bragg = {kAzure+10, kFourSquaresX},
        ms_clu = {ms_bragg.c, kMultiply, 1},
        ms_michel = {kGreen-8, kOpenDoubleDiamond},
        ms_shw = {kYellow+2, kOpenCircle, 0.5};

    ana::LineStyle
        ls_pass = {vc_pass.front(), kSolid, 2},
        ls_fail = {vc_fail.front(), kSolid, 2},
        ls_back = {kGray, kSolid, 1};

    bool IsUpright(recob::Track const& T);
    enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
};

ana::MichelDisplay::MichelDisplay(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}, MichelDisplayer{p},
    fLog(p.get<bool>("Log", true)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 40.F)), // in cm
    fNearbyRadius(p.get<float>("NearbyRadius", 40.F)), //in cm
    fBodyDistance(p.get<float>("BodyDistance", 20.F)), //in cm
    fRegN(p.get<unsigned>("RegN", 6)),
    fBraggThreshold(p.get<float>("BraggThreshold", 1.7)) // in MIP dE/dx
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    wireWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(), 
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        default: break;
    }
    fCathodeGap = geoHighX.x.min - geoLowX.x.max;

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX << std::endl
        << "  LowX Bounds: " << geoLowX << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl;
}

void ana::MichelDisplay::beginJob() {
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

void ana::MichelDisplay::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) return;
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) return;
    VecPtrTrk vpt_ev;
    art::fill_ptr_vector(vpt_ev, vh_trk);

    auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    if (!vh_shw.isValid()) return;
    VecPtrShw vps_ev;
    art::fill_ptr_vector(vps_ev, vh_shw);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    std::vector<char const*> cuts = {
        "None",
        // Form("TrackLength >= %.0f cm", fTrackLengthCut),
        // "EndInVolume (20 cm)",
        // "CathodeCrossing",
        // "AnodeCrossing (20 cm)",
        "HasMichelHits",
        "BraggEndNoError",
        Form("BraggThreshold %.1f MIP", fBraggThreshold)
    };

    gStyle->SetPalette(pal);

    std::vector<TCanvas*> hcs;
    for (unsigned ihc=0; ihc<cuts.size(); ihc++) {
        TCanvas* hc = asFile->make<TCanvas>(
            Form("hc%ucut%u", ev, ihc),
            Form("Hits: run:%u, subrun:%u, event:%u, cut:%s", e.run(), e.subRun(), e.event(), cuts[ihc]),
            1300,800
        );
        ana::DrawFrame(hc, int(geoDet), cuts[ihc], Form("%s R:%u SR:%u E:%u", (e.isRealData()?"Data":"Simulation"), e.run(), e.subRun(), e.event()));
        DrawGraph(hc, vph_ev, "p", ms_ev);

        for (PtrShw const& ps_ev : vps_ev) {
            VecPtrHit vph_shw = fmp_shw2hit.at(ps_ev.key());
            if (!vph_shw.size()) continue;
            DrawGraph(hc, vph_shw, "p", ms_shw);
        }

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
    for (PtrTrk const& pt_ev : vpt_ev) {
        VecPtrHit vph_muon = fmp_trk2hit.at(pt_ev.key());
        ASSERT(vph_muon.size())

        bool isUpright =  IsUpright(*pt_ev);
        geo::Point_t Start = isUpright ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = isUpright ? pt_ev->End() : pt_ev->Start();
        ana::SortedHits sh_mu = GetSortedHits(vph_muon, End.Z() > Start.Z() ? 1 : -1);
        ASSERT(sh_mu)


        im++;

        auto DrawFail = [&](std::vector<TCanvas*>::iterator ihc, std::vector<TCanvas*>::iterator itc) -> void {
            Color_t c_fail = vc_fail[im%vc_fail.size()];
            DrawGraph(*ihc, sh_mu.vph, "l", {}, {c_fail, ls_fail.l, ls_fail.w});
            DrawGraph2D(*itc, pt_ev, {}, {c_fail, ls_fail.l, ls_fail.w});
            for (auto jhc=ihc+1; jhc!=hcs.end(); jhc++)
                DrawGraph(*jhc, sh_mu.vph, "l", {}, ls_back);
            for (auto jtc=itc+1; jtc!=tcs.end(); jtc++)
                DrawGraph2D(*jtc, pt_ev, ms_back);
        };
        auto DrawPass = [&](std::vector<TCanvas*>::iterator ihc, std::vector<TCanvas*>::iterator itc) -> void {
            Color_t c_pass = vc_pass[im%vc_pass.size()];
            DrawGraph(*ihc, sh_mu.vph, "l", {}, {c_pass, ls_pass.l, ls_pass.w});
            DrawGraph2D(*itc, pt_ev, {}, {c_pass, ls_pass.l, ls_pass.w});
            for (PtrHit const& sc : sh_mu.sc)
                DrawMarker(*ihc, sc, ms_sc);
            DrawMarker(*ihc, sh_mu.start, ms_end);
            DrawMarker(*ihc, sh_mu.end, ms_end);
            if (sh_mu.is_cc()) {
                DrawMarker(*ihc, sh_mu.cc.first, ms_cc);
                DrawMarker(*ihc, sh_mu.cc.second, ms_cc);
            }
        };

        bool TagTrackLength = pt_ev->Length() >= fTrackLengthCut;
        bool TagEndInVolume = geoHighX.isInsideYZ(End, 20.F);
        bool TagEndInWindow = wireWindow.isInside(sh_mu.end->PeakTime(), 20.F / fTick2cm);
        bool TagCathodeCrossing = (
            geoLowX.isInside(Start)
            && geoHighX.isInside(End)
        ) || (
            geoHighX.isInside(Start)
            && geoLowX.isInside(End)
        );
        bool TagAnodeCrossing = false;
        if (geoDet == kPDVD)
            TagAnodeCrossing = geoHighX.isInsideYZ(Start, 20.);
        else if (geoDet == kPDHD)
            TagAnodeCrossing = geoHighX.isInsideYZ(Start, 20.) || geoLowX.isInsideYZ(Start, 20.);

        simb::MCParticle const* mcp_mu = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        simb::MCParticle const* mcp_mi = GetMichelMCP(mcp_mu);
        VecPtrHit vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true);
        bool TrueTagHasMichelHits = !vph_mi.empty();

        // unused variables
        LOG(TagTrackLength);
        LOG(TagEndInVolume);
        LOG(TagEndInWindow);
        LOG(TagCathodeCrossing);
        LOG(TagAnodeCrossing);
        LOG(TrueTagHasMichelHits);

        std::vector<TCanvas*>::iterator ihc = hcs.begin();
        std::vector<TCanvas*>::iterator itc = tcs.begin();

        DrawPass(ihc, itc);
        ihc++; itc++;

        if (!LOG(TrueTagHasMichelHits)) {
            DrawFail(ihc, itc);
            continue;
        }
        DrawPass(ihc, itc);
        DrawGraph(*ihc, vph_mi, "p", ms_michel);
        ihc++; itc++;

        ana::Bragg bragg = GetBragg(
            sh_mu.vph,
            sh_mu.end,
            pt_ev,
            vph_ev,
            fop_hit2trk,
            { fBodyDistance, fRegN, fTrackLengthCut, fNearbyRadius }
        );

        if (!LOG(bragg)) {
            DrawFail(ihc, itc);
            continue;
        }
        
        // Sphere around track end
        VecPtrHit vph_pandora;
        int sec_end = ana::tpc2sec[geoDet][sh_mu.end->WireID().TPC];
        for (PtrHit const& ph_ev : vph_ev) {
            if (ph_ev->View() != geo::kW) continue;
            if (ana::tpc2sec[geoDet][ph_ev->WireID().TPC] != sec_end) continue;
            PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
            // PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
            // if (ps_hit) continue;
            if (GetDistance(ph_ev, sh_mu.end) > 20.F) continue;
            if (pt_hit && pt_hit->Length() > fTrackLengthCut) continue;
            vph_pandora.push_back(ph_ev);
        }

        // Sphere around bragg end
        VecPtrHit::iterator iph_bragg = std::find_if(
            bragg.vph_clu.begin(), bragg.vph_clu.end(),
            [&](PtrHit const& h) -> bool { return h.key() == bragg.end.key(); }
        );
        if (iph_bragg != bragg.vph_clu.end()) iph_bragg++;
        VecPtrHit vph_bragg;
        for (auto iph=iph_bragg; iph!=bragg.vph_clu.end(); iph++) {
            if (GetDistance(*iph, bragg.end) > 20.F) continue;
            vph_bragg.push_back(*iph);
        }

        // Cone around bragg end
        VecPtrHit vph_cone;
        Hits bary_hits;
        Vec2 bary;
        for (auto iph=iph_bragg; iph!=bragg.vph_clu.end(); iph++) {
            if ((*iph)->View() != geo::kW) continue;
            if (GetDistance(*iph, bragg.end) > 10) continue;
            bary_hits.push_back(GetHit(*iph));
        }
        if (bary_hits.size()) {
            Hit h_end = GetHit(bragg.end);
            bary = bary_hits.barycenter(h_end.section, fTick2cm);
            Vec2 end = h_end.vec(fTick2cm);
            Vec2 end_bary = bary - end;

            // float angle = end_bary.angle();
            for (auto iph=iph_bragg; iph!=bragg.vph_clu.end(); iph++) {
                if ((*iph)->View() != geo::kW) continue;
                if (GetDistance(*iph, bragg.end) > 30) continue;

                Vec2 end_hit = GetHit(*iph).vec(fTick2cm) - end;
                float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                if (cosa > cos(30.F * TMath::DegToRad())) continue;
                vph_cone.push_back(*iph);
            }
        }


        DrawPass(ihc, itc);
        DrawGraph(*ihc, vph_mi, "p", ms_michel);
        DrawMarker(*ihc, bragg.end, ms_bragg);

        DrawGraph(*ihc, bragg.vph_clu, "l", {}, {ms_clu.c, kDashed, 1} );
        DrawGraph(*ihc, vph_bragg, "p", {kRed, kOpenCircle, .5});
        DrawGraph(*ihc, vph_pandora, "p", {kOrange, kOpenCircle, 1.5});

        TMarker *m_bary = new TMarker();
        SetMarkerStyle(m_bary, {kMagenta, kOpenFourTrianglesX, 2});
        (*ihc)->cd(ana::tpc2sec[geoDet][bragg.end->WireID().TPC]+1);
        m_bary->DrawMarker(bary.space, bary.drift);
        DrawGraph(*ihc, vph_cone, "p", {kMagenta, kOpenTriangleUp, 2});

        ihc++; itc++;

        if (!LOG(bragg.max_dQdx >= fBraggThreshold * bragg.mip_dQdx)) {
            DrawFail(ihc, itc);
            continue;
        }
        DrawPass(ihc, itc);
        DrawGraph(*ihc, vph_mi, "p", ms_michel);
        DrawMarker(*ihc, bragg.end, ms_bragg);

        DrawGraph(*ihc, bragg.vph_clu, "l", {}, {ms_clu.c, kDashed, 1} );
        DrawGraph(*ihc, vph_bragg, "p", {kRed, kOpenCircle, .5});
        DrawGraph(*ihc, vph_pandora, "p", {kOrange, kOpenCircle, 1.5});

        (*ihc)->cd(ana::tpc2sec[geoDet][bragg.end->WireID().TPC]+1);
        m_bary->DrawMarker(bary.space, bary.drift);
        DrawGraph(*ihc, vph_cone, "p", {kMagenta, kOpenTriangleUp, 2});
    }

    for (TCanvas* hc : hcs)
        hc->Write();
    for (TCanvas* tc : tcs)
        tc->Write();
    
    ev++;
}

void ana::MichelDisplay::endJob() {}

bool ana::MichelDisplay::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}

DEFINE_ART_MODULE(ana::MichelDisplay)