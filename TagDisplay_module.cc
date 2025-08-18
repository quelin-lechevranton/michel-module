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
        ms_bragg = {kAzure+10, kFourSquaresX};

    ana::LineStyle
        ls_pass = {vc_pass.front(), kSolid, 2},
        ls_fail = {vc_fail.front(), kSolid, 2},
        ls_back = {kGray, kSolid, 1};

    bool IsUpright(recob::Track const& T);
    enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
    PtrHit GetBraggEnd(
        VecPtrHit const& vph_trk,
        PtrHit const& ph_trk_end,
        PtrTrk const& p_trk,
        VecPtrHit const& vph_ev,
        art::FindOneP<recob::Track> const& fop_hit2trk,
        VecPtrHit *vph_sec_bragg = nullptr,
        float *max_dQdx = nullptr,
        float *mip_dQdx = nullptr,
        int *error = nullptr
    );
    double dist2(PtrHit const& ph1, PtrHit const& ph2);
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
    std::vector<art::Ptr<recob::Shower>> vp_shw;
    art::fill_ptr_vector(vp_shw, vh_shw);

    art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    std::vector<char const*> cuts = {
        "None",
        // Form("TrackLength >= %.0f cm", fTrackLengthCut),
        // "EndInVolume (20 cm)",
        // "CathodeCrossing",
        // "AnodeCrossing (20 cm)",
        "HasMichelHits",
        "BraggEndAlogithm",
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

        int TagBraggError = -1;
        float CutdQdxMax = 0.F;
        float MIPdQdx = 0.F;
        VecPtrHit vph_bragg_muon;
        PtrHit ph_bragg = GetBraggEnd(
            sh_muon.vph,
            sh_muon.end,
            pt_ev,
            vph_ev,
            fop_hit2trk,
            &vph_bragg_muon,
            &CutdQdxMax,
            &MIPdQdx,
            &TagBraggError
        );

        im++;

        std::vector<TCanvas*>::iterator ihc = hcs.begin();
        std::vector<TCanvas*>::iterator itc = tcs.begin();
        auto DrawFilter = [&](bool tag) -> bool {
            if (!tag) {
                Color_t c_fail = vc_fail[im%vc_fail.size()];
                DrawGraph(*ihc, sh_muon.vph, "l", {}, {c_fail, ls_fail.l, ls_fail.w});
                DrawGraph2D(*itc, pt_ev, {}, {c_fail, ls_fail.l, ls_fail.w});
                for (auto jhc=ihc+1; jhc!=hcs.end(); jhc++)
                    DrawGraph(*jhc, sh_muon.vph, "l", {}, ls_back);
                for (auto jtc=itc+1; jtc!=tcs.end(); jtc++)
                    DrawGraph2D(*jtc, pt_ev, ms_back);
                return false;
            }
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
            ihc++;
            itc++;
            return true;
        };

        bool TagTrackLength = pt_ev->Length() >= fTrackLengthCut;
        bool TagEndInVolume = geoHighX.isInsideYZ(End, 20.F);
        // MuonEndIsInWindowT = wireWindow.isInside(MuonEndHit.tick, fMichelTickRadius);
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

        DrawFilter(true);
        // if (!DrawFilter(LOG(TagTrackLength))) continue;
        // if (!DrawFilter(LOG(TagEndInVolume))) continue;
        // if (!DrawFilter(LOG(TagCathodeCrossing))) continue;
        // if (!DrawFilter(LOG(TagAnodeCrossing))) continue;
        // if (!DrawFilter(LOG(TagBraggError == kNoError))) continue;
        // DrawMarker(*(ihc-1), ph_bragg, ms_bragg);
        // DrawGraph(*(ihc-1), vph_bragg_muon, "p", {ms_bragg.c, kMultiply, 0.5});
        // if (!DrawFilter(LOG(CutdQdxMax >= fBraggThreshold * MIPdQdx))) continue;
        // DrawMarker(*(ihc-1), ph_bragg, ms_bragg);
        // DrawGraph(*(ihc-1), vph_bragg_muon, "p", {ms_bragg.c, kMultiply, 0.5});
        if (!DrawFilter(LOG(TrueTagHasMichelHits))) continue;
        DrawMarker(*(ihc-1), ph_bragg, ms_bragg);
        DrawGraph(*(ihc-1), vph_bragg_muon, "p", {ms_bragg.c, kMultiply, 0.5});
        if (!DrawFilter(LOG(TagBraggError == kNoError))) continue;
        DrawMarker(*(ihc-1), ph_bragg, ms_bragg);
        DrawGraph(*(ihc-1), vph_bragg_muon, "p", {ms_bragg.c, kMultiply, 0.5});
        if (!DrawFilter(LOG(CutdQdxMax >= fBraggThreshold * MIPdQdx))) continue;
        DrawMarker(*(ihc-1), ph_bragg, ms_bragg);
        DrawGraph(*(ihc-1), vph_bragg_muon, "p", {ms_bragg.c, kMultiply, 0.5});
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

PtrHit ana::TagDisplay::GetBraggEnd(
    VecPtrHit const& vph_trk,
    PtrHit const& ph_trk_end,
    PtrTrk const& p_trk,
    VecPtrHit const& vph_ev,
    art::FindOneP<recob::Track> const& fop_hit2trk,
    VecPtrHit *vph_sec_bragg,
    float *max_dQdx,
    float *mip_dQdx,
    int *error
) {
    VecPtrHit vph_sec_trk;
    for (PtrHit const& p_hit : vph_trk)
        if (ana::tpc2sec[geoDet][p_hit->WireID().TPC]
            == ana::tpc2sec[geoDet][ph_trk_end->WireID().TPC])
            vph_sec_trk.push_back(p_hit);
    
    if (vph_sec_trk.empty()
        || std::find_if(
            vph_sec_trk.begin(), vph_sec_trk.end(),
            [key=ph_trk_end.key()](PtrHit const& ph) -> bool {
                return ph.key() == key;
            }
        ) == vph_sec_trk.end()
    ) {
        if (error) *error = kEndNotFound;
        return PtrHit{};
    }

    std::sort(
        vph_sec_trk.begin(),
        vph_sec_trk.end(),
        [&](PtrHit const& ph1, PtrHit const& ph2) -> bool {
            return GetDistance(ph2, ph_trk_end) > GetDistance(ph1, ph_trk_end);
        }
    );

    VecPtrHit::iterator iph_body = std::find_if(
        vph_sec_trk.begin(),
        vph_sec_trk.end(),
        [&](PtrHit const& h) -> bool {
            return GetDistance(h, ph_trk_end) > fBodyDistance;
        }
    );

    if (mip_dQdx) {
        VecPtrHit::iterator jph_body = std::find_if(
            vph_sec_trk.begin(),
            vph_sec_trk.end(),
            [&](PtrHit const& h) -> bool {
                return GetDistance(h, ph_trk_end) > 2*fBodyDistance;
            }
        );
        float mean_dQ = 0;
        float mean_dx = 0;
        for (auto iph=iph_body; iph!=jph_body; iph++) {
            mean_dQ += (*iph)->Integral();
            mean_dx += iph==vph_sec_trk.begin()
                ? GetDistance(*iph, *(iph+1))
                : ( iph==vph_sec_trk.end()-1
                    ? GetDistance(*(iph-1), *iph)
                    : .5*GetDistance(*(iph-1), *(iph+1))
                );
        }
        *mip_dQdx = mean_dQ / mean_dx;
    }

    if (std::distance(iph_body, vph_sec_trk.end()) < fRegN) {
        if (error) *error = kSmallBody;
        return PtrHit{};
    }

    VecPtrHit vph_near;
    for (PtrHit const& ph_ev : vph_ev) {
        if (ana::tpc2sec[geoDet][ph_ev->WireID().TPC]
            != ana::tpc2sec[geoDet][ph_trk_end->WireID().TPC]) continue;

        // check if the hit is in a track
        PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
        if (pt_hit) {
            // check if the hit is on the body of the track
            if (pt_hit.key() == p_trk.key()
                && std::find_if(
                    iph_body, vph_sec_trk.end(),
                    [key=ph_ev.key()](PtrHit const& ph) -> bool {
                        return ph.key() == key;
                    }
                ) != vph_sec_trk.end()
            ) continue;
            // check if the hit is on a long track
            if (
                pt_hit.key() != p_trk.key()
                && pt_hit->Length() > fTrackLengthCut
            ) continue;
        }

        // check if the hit is close enough
        if (GetDistance(ph_ev, ph_trk_end) > fNearbyRadius) continue;

        vph_near.push_back(ph_ev);
    }
    DEBUG(vph_near.empty())

    VecPtrHit vph_reg{iph_body, iph_body+fRegN};
    std::reverse(vph_reg.begin(), vph_reg.end());
    auto orientation = [&](VecPtrHit const& vph) -> std::pair<double, double> {
        LinearRegression reg;
        for (PtrHit const& ph : vph) {
            double z = GetSpace(ph->WireID());
            double t = ph->PeakTime() * fTick2cm;
            reg.add(z, t);
        }
        reg.compute();
        // DEBUG(reg.corr == 0)
        int dirz = GetSpace(vph.back()->WireID())
            > GetSpace(vph.front()->WireID())
            ? 1 : -1;
        int dirt = vph.back()->PeakTime()
            > vph.front()->PeakTime()
            ? 1 : -1;
        double theta = reg.theta(dirz, dirt);
        double sigma = TMath::Pi() / 4 / reg.corr;
        return std::make_pair(theta, sigma);
    };
    std::pair<double, double> reg = orientation(vph_reg);

    auto score = [&](PtrHit const& ph1, PtrHit const& ph2, double theta, double sigma) -> double {
        double dz = GetSpace(ph2->WireID()) - GetSpace(ph1->WireID());
        double dt = (ph2->PeakTime() - ph1->PeakTime()) * fTick2cm;
        double da = atan2(dt, dz) - theta;
        da = abs(da) < TMath::Pi() ? da : da - (da>0?1:-1)*2*TMath::Pi();
        double r = sqrt(dt*dt + dz*dz);

        return TMath::Gaus(da, 0, sigma) / r;
    };

    VecPtrHit vph_sec{vph_reg.begin(), vph_reg.end()};
    // std::reverse(vph_sec.begin(), vph_sec.end());
    PtrHit ph_prev = ph_trk_end;
    while (vph_near.size()) {
        VecPtrHit::iterator iph_max = std::max_element(
            vph_near.begin(), vph_near.end(),
            [&](PtrHit const& ph1, PtrHit const& ph2) -> bool {
                return score(ph_prev, ph1, reg.first, reg.second)
                    < score(ph_prev, ph2, reg.first, reg.second);
            }
        );

        vph_near.erase(iph_max);
        // vph_reg.insert(vph_reg.begin(), *iph_max);
        // vph_reg.pop_back();
        vph_reg.erase(vph_reg.begin());
        vph_reg.push_back(*iph_max);
        vph_sec.push_back(*iph_max);
        ph_prev = *iph_max;

        reg = orientation(vph_reg);
    }

    unsigned const trailing_radius = 6;
    double max = std::numeric_limits<double>::lowest();
    PtrHit ph_max;
    for (auto iph_sec=vph_sec.begin(); iph_sec!=vph_sec.end(); iph_sec++) {
        VecPtrHit::iterator jph_sec = std::distance(iph_sec, vph_sec.end()) > trailing_radius
            ? iph_sec+trailing_radius : vph_sec.end();
        unsigned l = std::distance(iph_sec, jph_sec);

        double dQ = std::accumulate(
            iph_sec, jph_sec, 0.,
            [](double sum, PtrHit const& ph) {
                return sum+ph->Integral();
            }
        );
        dQ /= l;

        double dx = iph_sec == vph_sec.begin()
            ? GetDistance(*iph_sec, *(iph_sec+1))
            : ( iph_sec == vph_sec.end()-1
                ? GetDistance(*(iph_sec-1), *iph_sec)
                : .5*GetDistance(*(iph_sec-1), *(iph_sec+1))
            );
        for (auto iph=iph_sec+1; iph!=jph_sec; iph++)
            dx += iph == vph_sec.end()-1
                ? GetDistance(*(iph-1), *iph)
                : .5*GetDistance(*(iph-1), *(iph+1));
        dx /= l;

        double dQdx = dQ / dx;
        if (dQdx > max) {
            max = dQdx;
            ph_max = *iph_sec;
        }
    }
    if (vph_sec_bragg) *vph_sec_bragg = vph_sec;
    if (max_dQdx) *max_dQdx = float(max);
    if (error) *error = kNoError;
    return ph_max;
}   

DEFINE_ART_MODULE(ana::TagDisplay)