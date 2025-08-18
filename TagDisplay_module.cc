////////////////////////////////////////////////////////////////////////
// Class:       TagDisplay
// Plugin Type: analyzer (Unknown Unknown)
// File:        TagDisplay_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "event_display.h"

namespace ana {
    class TagDisplay;
}

class ana::TagDisplay : public art::EDAnalyzer, private ana::MichelDisplayer {
public:
    explicit TagDisplay(fhicl::ParameterSet const& p);
    TagDisplay(TagDisplay const&) = delete;
    TagDisplay(TagDisplay&&) = delete;
    TagDisplay& operator=(TagDisplay const&) = delete;
    TagDisplay& operator=(TagDisplay&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    bounds<float> wireWindow;
    bounds3D<float> geoHighX, geoLowX;
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
        ms_michel = {kGreen-8, kOpenDoubleDiamond},
        ms_shw = {kYellow+2, kOpenCircle, 0.5};

    ana::LineStyle
        ls_pass = {vc_pass.front(), kSolid, 2},
        ls_fail = {vc_fail.front(), kSolid, 2},
        ls_back = {kGray, kSolid, 1};

    bool IsUpright(recob::Track const& T);
    enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
};

ana::TagDisplay::TagDisplay(fhicl::ParameterSet const& p)
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
    wireWindow = bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoLowX = bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(), 
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoHighX = bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoLowX = bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoHighX = bounds3D<float>{
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

void ana::TagDisplay::beginJob() {
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

void ana::TagDisplay::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) return;
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) return;
    std::vector<PtrTrk> vpt_ev;
    art::fill_ptr_vector(vpt_ev, vh_trk);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    if (!vh_shw.isValid()) return;
    std::vector<art::Ptr<recob::Shower>> vps_ev;
    art::fill_ptr_vector(vps_ev, vh_shw);

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
        ana::SortedHits sh_muon = GetSortedHits(vph_muon, End.Z() > Start.Z() ? 1 : -1);
        ASSERT(sh_muon)

        ana::Bragg bragg = GetBragg(
            sh_muon.vph,
            sh_muon.end,
            pt_ev,
            vph_ev,
            fop_hit2trk,
            { fBodyDistance, fRegN, fTrackLengthCut, fNearbyRadius }
        );

        im++;

        auto DrawFail = [&](std::vector<TCanvas*>::iterator ihc, std::vector<TCanvas*>::iterator itc) -> void {
            Color_t c_fail = vc_fail[im%vc_fail.size()];
            DrawGraph(*ihc, sh_muon.vph, "l", {}, {c_fail, ls_fail.l, ls_fail.w});
            DrawGraph2D(*itc, pt_ev, {}, {c_fail, ls_fail.l, ls_fail.w});
            for (auto jhc=ihc+1; jhc!=hcs.end(); jhc++)
                DrawGraph(*jhc, sh_muon.vph, "l", {}, ls_back);
            for (auto jtc=itc+1; jtc!=tcs.end(); jtc++)
                DrawGraph2D(*jtc, pt_ev, ms_back);
        };
        auto DrawPass = [&](std::vector<TCanvas*>::iterator ihc, std::vector<TCanvas*>::iterator itc) -> void {
            Color_t c_pass = vc_pass[im%vc_pass.size()];
            DrawGraph(*ihc, sh_muon.vph, "l", {}, {c_pass, ls_pass.l, ls_pass.w});
            DrawGraph2D(*itc, pt_ev, {}, {c_pass, ls_pass.l, ls_pass.w});
            for (PtrHit const& sc : sh_muon.sc)
                DrawMarker(*ihc, sc, ms_sc);
            DrawMarker(*ihc, sh_muon.vph.front(), ms_end);
            DrawMarker(*ihc, sh_muon.vph.back(), ms_end);
            if (sh_muon.isCathodeCrossing) {
                DrawMarker(*ihc, sh_muon.cc.first, ms_cc);
                DrawMarker(*ihc, sh_muon.cc.second, ms_cc);
            }
        };

        bool TagTrackLength = pt_ev->Length() >= fTrackLengthCut;
        bool TagEndInVolume = geoHighX.isInsideYZ(End, 20.F);
        bool TagEndInWindow = wireWindow.isInside(sh_muon.end->PeakTime(), 20.F / fTick2cm);
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

        if (!LOG(bragg.error == kNoError)) {
            DrawFail(ihc, itc);
            continue;
        }
        DrawPass(ihc, itc);
        DrawGraph(*ihc, vph_mi, "p", ms_michel);
        DrawMarker(*ihc, bragg.end, ms_bragg);
        DrawGraph(*ihc, bragg.vph_muon, "p", {ms_bragg.c, kMultiply, 0.5});
        ihc++; itc++;

        if (!LOG(bragg.max_dQdx >= fBraggThreshold * bragg.mip_dQdx)) {
            DrawFail(ihc, itc);
            continue;
        }
        DrawPass(ihc, itc);
        DrawGraph(*ihc, vph_mi, "p", ms_michel);
        DrawMarker(*ihc, bragg.end, ms_bragg);
        DrawGraph(*ihc, bragg.vph_muon, "p", {ms_bragg.c, kMultiply, 0.5});
    }

    for (TCanvas* hc : hcs)
        hc->Write();
    for (TCanvas* tc : tcs)
        tc->Write();
    
    ev++;
}

void ana::TagDisplay::endJob() {}

bool ana::TagDisplay::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}

DEFINE_ART_MODULE(ana::TagDisplay)