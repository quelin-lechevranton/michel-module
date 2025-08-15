////////////////////////////////////////////////////////////////////////
// Class:       Trackchecks
// Plugin Type: analyzer (Unknown Unknown)
// File:        Trackchecks_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "event_display.h"

using HitPtr = art::Ptr<recob::Hit>;
using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

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
    bool fSort;
    float fMichelSpaceRadius; // in cm
    float fMichelTickRadius; // in ticks

    unsigned cn=0;
    unsigned gn=0;

    Int_t pal = kCividis;
    MarkerStyle
        ms_ev = {kGray, kMultiply, 0.5},
        ms_mcp_mu = {kSpring+2, kOpenDiamond},
        ms_mcp_me = {kAzure-4, kOpenDoubleDiamond},
        ms_mu_in = {kOrange+6, kFullSquare},
        ms_mu_out = {kOrange+6, kOpenSquare},
        ms_cc_in = {kOrange+6, kFullTriangleDown},
        ms_cc_out = {kOrange+6, kOpenTriangleDown},
        ms_cc = {kPink-2, kFullTriangleUp},
        ms_sc = {kPink-2, kFullCircle},
        ms_shw = {kGreen-8, kOpenCircle};

    LineStyle
        ls_mu = {kOrange+6, kSolid},
        ls_small = {kGray, kSolid},
        ls_broken = {kRed-4, kDashed},
        ls_no_end = {kGreen-8, kSolid},
        ls_reg_good = {kAzure-4, kSolid, 2},
        ls_reg_bad = {kViolet+6, kSolid, 2};


    double GetSpace(geo::WireID);
    ana::Hit GetHit(HitPtr const p_hit);
    HitPtrPair GetTrackEndsHits( // return hits at track ends
        HitPtrVec const&, // all hits of the track
        HitPtrPair *pp_cathode_crossing = nullptr, // return hits at cathode crossing
        HitPtrVec *vp_tpc_crossing = nullptr, // return hits at TPC crossing (PDVD)
        HitPtrVec *vp_sorted_hit = nullptr, // return hits of the track per section and sorted
        geo::View_t = geo::kW // view to consider
    );
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
}

void ana::Trackchecks::beginJob() {
    TCanvas *c = asFile->make<TCanvas>(
        "legend", "legend",
        1300, 800
    );
    c->cd();

    float top_margin = 0.1;
    float indent = 0.05;
    float title_indent = 0.04;
    float label_indent = 0.1;
    float line_height = 0.05;
    float line_length = 0.03;
    float ncol = 4;
    Style_t font = 43;
    Style_t font_bold = 63;
    Float_t title_size = 20;
    Float_t label_size = 16;

    TText* title = new TText();
    title->SetTextFont(font_bold);
    title->SetTextSize(title_size);
    TText* label = new TText();
    label->SetTextFont(font);
    label->SetTextSize(label_size);

    title->DrawText(title_indent, 1 - top_margin, "Event Markers");

    TMarker* m_ev = new TMarker();
    SetMarkerStyle(m_ev, ms_ev);
    m_ev->DrawMarker(indent, 1 - top_margin - line_height);
    label->DrawText(label_indent, 1 - top_margin - line_height, "Event Hits");

    TMarker* m_mcp_mu = new TMarker();
    SetMarkerStyle(m_mcp_mu, ms_mcp_mu);
    m_mcp_mu->DrawMarker(indent, 1 - top_margin - 2 * line_height);
    label->DrawText(label_indent, 1 - top_margin - 2 * line_height, "MCParticle Muon End Hit");

    TMarker* m_mcp_me = new TMarker();
    SetMarkerStyle(m_mcp_me, ms_mcp_me);
    m_mcp_me->DrawMarker(indent, 1 - top_margin - 3 * line_height);
    label->DrawText(label_indent, 1 - top_margin - 3 * line_height, "MCParticle Michel Hits");

    title->DrawText(1/ncol + title_indent, 1 - top_margin, "Tracks");

    TLine* l_small = new TLine();
    SetLineStyle(l_small, ls_small);
    l_small->DrawLine(1/ncol + indent, 1 - top_margin - line_height,
                      1/ncol + indent + line_length, 1 - top_margin - line_height);
    label->DrawText(1/ncol + label_indent, 1 - top_margin - line_height, "Small Tracks");

    TLine* l_broken = new TLine();
    SetLineStyle(l_broken, ls_broken);
    l_broken->DrawLine(1/ncol + indent, 1 - top_margin - 2 * line_height,
                       1/ncol + indent + line_length, 1 - top_margin - 2 * line_height);
    label->DrawText(1/ncol + label_indent, 1 - top_margin - 2 * line_height, "Broken Tracks");

    TLine* l_no_end = new TLine();
    SetLineStyle(l_no_end, ls_no_end);
    l_no_end->DrawLine(1/ncol + indent, 1 - top_margin - 3 * line_height,
                       1/ncol + indent + line_length, 1 - top_margin - 3 * line_height);
    label->DrawText(1/ncol + label_indent, 1 - top_margin - 3 * line_height, "End Algorithm failed");

    TLine* l_mu = new TLine();
    SetLineStyle(l_mu, ls_mu);
    l_mu->DrawLine(1/ncol + indent, 1 - top_margin - 4 * line_height,
                   1/ncol + indent + line_length, 1 - top_margin - 4 * line_height);
    label->DrawText(1/ncol + label_indent, 1 - top_margin - 4 * line_height, "Tracks");

    TLine* l_reg_good = new TLine();
    SetLineStyle(l_reg_good, ls_reg_good);
    l_reg_good->DrawLine(1/ncol + indent, 1 - top_margin - 5 * line_height,
                         1/ncol + indent + line_length, 1 - top_margin - 5 * line_height);
    label->DrawText(1/ncol + label_indent, 1 - top_margin - 5 * line_height, "Regression Good");

    TLine* l_reg_bad = new TLine();
    SetLineStyle(l_reg_bad, ls_reg_bad);
    l_reg_bad->DrawLine(1/ncol + indent, 1 - top_margin - 6 * line_height,
                         1/ncol + indent + line_length, 1 - top_margin - 6 * line_height);
    label->DrawText(1/ncol + label_indent, 1 - top_margin - 6 * line_height, "Regression Bad");

    title->DrawText(2/ncol + title_indent, 1 - top_margin, "Track Ends Markers");

    TMarker* m_mu_in = new TMarker();
    SetMarkerStyle(m_mu_in, ms_mu_in);
    m_mu_in->DrawMarker(2/ncol + indent, 1 - top_margin - line_height);
    label->DrawText(2/ncol + label_indent, 1 - top_margin - line_height, "Muon Track In");

    TMarker* m_mu_out = new TMarker();
    SetMarkerStyle(m_mu_out, ms_mu_out);
    m_mu_out->DrawMarker(2/ncol + indent, 1 - top_margin - 2 * line_height);
    label->DrawText(2/ncol + label_indent, 1 - top_margin - 2 * line_height, "Muon Track Out");

    TMarker* m_mu_cc_in = new TMarker();
    SetMarkerStyle(m_mu_cc_in, ms_cc_in);
    m_mu_cc_in->DrawMarker(2/ncol + indent, 1 - top_margin - 3 * line_height);
    label->DrawText(2/ncol + label_indent, 1 - top_margin - 3 * line_height, "Cathode Crossing Muon In");

    TMarker* m_mu_cc_out = new TMarker();
    SetMarkerStyle(m_mu_cc_out, ms_cc_out);
    m_mu_cc_out->DrawMarker(2/ncol + indent, 1 - top_margin - 4 * line_height);
    label->DrawText(2/ncol + label_indent, 1 - top_margin - 4 * line_height, "Cathode Crossing Muon Out");

    TMarker* m_sc = new TMarker();
    SetMarkerStyle(m_sc, ms_sc);
    m_sc->DrawMarker(2/ncol + indent, 1 - top_margin - 5 * line_height);
    label->DrawText(2/ncol + label_indent, 1 - top_margin - 5 * line_height, "Section Crossing");

    TMarker* m_shw = new TMarker();
    SetMarkerStyle(m_shw, ms_shw);
    m_shw->DrawMarker(2/ncol + indent, 1 - top_margin - 6 * line_height);
    label->DrawText(2/ncol + label_indent, 1 - top_margin - 6 * line_height, "Shower Hits");

    c->Write();
}

void ana::Trackchecks::analyze(art::Event const& e) {
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

    TCanvas *c = asFile->make<TCanvas>(
        Form("c%u", cn++),
        Form("run:%u, subrun:%u, event:%u", e.run(), e.subRun(), e.event()),
        1300,800
    );
    ana::DrawFrame(c, int(geoDet), "", Form("%s R:%u SR:%u E:%u", (e.isRealData()?"Data":"Simulation"), e.run(), e.subRun(), e.event()));

    auto drawMarker = [&](HitPtr const& p_hit, MarkerStyle const& ms) -> void {
        TMarker *m = new TMarker();
        SetMarkerStyle(m, ms);
        if (geoDet == kPDVD) {
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            c->cd(s+1);
            m->DrawMarker(GetSpace(p_hit->WireID()), p_hit->PeakTime() * fTick2cm);
        } else if (geoDet == kPDHD) {
            int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) return;
            c->cd(s+1);
            m->DrawMarker(p_hit->PeakTime() * fTick2cm, GetSpace(p_hit->WireID()));
        }
    };

    auto drawGraph = [&](HitPtrVec const& vp_hit, char const* draw, MarkerStyle const& ms={}, LineStyle const& ls={}) -> void {
        std::vector<TGraph*> gs(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            gs[s] = new TGraph();
            gs[s]->SetEditable(kFALSE);
            gs[s]->SetName(Form("g%u_s%u", gn++, s));
            SetMarkerStyle(gs[s], ms);
            SetLineStyle(gs[s], ls);
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
            c->cd(s+1);
            gs[s]->Draw(draw);
        }
    };

    gStyle->SetPalette(pal);

    // all event hits TGraph
    drawGraph(vp_hit, "p", ms_ev);

    // all event hits simulate TScatter
    // TArrayI const& colors = TColor::GetPalette();
    // TMarker* m = new TMarker();
    // m->SetMarkerStyle(kFullCircle);
    // for (HitPtr p_hit : vp_hit) {
    //     if (p_hit->View() != geo::kW) continue;
    //     float const x = std::min(p_hit->Integral() / (geoDet == kPDVD ? 200.F : 1000.F), 1.F);
    //     m->SetMarkerSize(2*x+0.1);
    //     m->SetMarkerColor(colors[int((colors.GetSize()-1)*x)]);

    //     drawMarker(m, p_hit);
    // }


    // all event hits TH2F
    // std::vector<TH2F*> h2s(ana::n_sec[geoDet]);
    // if (geoDet == kPDVD) {
    //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
    //         h2s[s] = new TH2F(
    //             Form("h2_%u", s),
    //             "event hits",
    //             600, 0, 300,
    //             600, 0, 6000
    //         );
    //     }
    // } else if (geoDet == kPDHD) {
    //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
    //         h2s[s] = new TH2F(
    //             Form("h2_%u", s),
    //             "event hits",
    //             600, 0, 6000,
    //             600, 0, 464
    //         );
    //     }
    // }
    // for (HitPtr p_hit : vp_hit) {
    //     if (p_hit->View() != geo::kW) continue;
    //     int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
    //     if (s == -1) continue;
    //     // int q = std::min(p_hit->Integral() / 1000???200, 1.F);
    //     int q = p_hit->Integral();
    //     if (geoDet == kPDVD)
    //         h2s[s]->Fill(GetSpace(p_hit->WireID()), p_hit->PeakTime() * fTick2cm, q);
    //     else if (geoDet == kPDHD)
    //         h2s[s]->Fill(p_hit->PeakTime() * fTick2cm, GetSpace(p_hit->WireID()), q);
    // }
    // double max=0;
    // for (TH2F* h2 : h2s)
    //     max = h2->GetMaximum() > max ? h2->GetMaximum() : max;
    // for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
    //     h2s[s]->SetMinimum(0);
    //     h2s[s]->SetMaximum(max);

    //     c->cd(s+1);
    //     TH2F* f = (TH2F*) gPad->FindObject(Form("f%u", s));
    //     f->SetMinimum(0);
    //     f->SetMaximum(max);
    // }
    // for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
    //     c->cd(s+1);
    //     if ((geoDet == kPDVD && s==3) || (geoDet == kPDHD && s==1))
    //         h2s[s]->Draw("same colz");
    //     else
    //         h2s[s]->Draw("same col");
    // }




    auto const& vh_mcp = e.getHandle<std::vector<simb::MCParticle>>(tag_mcp);
    if (vh_mcp) {
        for (simb::MCParticle const& mcp : *vh_mcp) {
            if (abs(mcp.PdgCode()) != 13) continue; 
            // TGraph* g_mcp_mu = new TGraph();
            // g_mcp_mu->SetMarkerColor(theme.c_mcp_mu);
            // g_mcp_mu->SetMarkerStyle(theme.m_mcp_mu);
            // g_mcp_mu->SetMarkerSize(theme.s_mcp_mu);
            // drawGraph(g_mcp_mu, ana::mcp2hits(&mcp, vp_hit, clockData, false), "p");

            HitPtrVec vp_mcp_hit = ana::mcp2hits(&mcp, vp_hit, clockData, false);
            HitPtrPair mcp_ends = GetTrackEndsHits(vp_mcp_hit);
            if (mcp_ends.first && mcp_ends.second) {
                int dir_z = mcp.EndZ() > mcp.Vz() ? 1 : -1;
                float fz = GetSpace(mcp_ends.first->WireID());
                float sz = GetSpace(mcp_ends.second->WireID());
                HitPtr mcp_end = (sz-fz) * dir_z > 0 ? mcp_ends.second : mcp_ends.first;

                drawMarker(mcp_end, ms_mcp_mu);
            } 

            if (mcp.NumberDaughters() >= 3) {
                simb::MCParticle const * mcp_michel = nullptr;
                bool has_numu = false, has_nue = false;
                for (int i_dau=mcp.NumberDaughters()-3; i_dau<mcp.NumberDaughters(); i_dau++) {
                    simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp.Daughter(i_dau));    
                    if (!mcp_dau) continue;
                    switch (abs(mcp_dau->PdgCode())) {
                        case 14: has_numu = true; break;
                        case 12: has_nue = true; break;
                        case 11: mcp_michel = mcp_dau; break;
                        default: break;
                    }
                }
                if (mcp_michel and has_numu and has_nue) {
                    drawGraph(ana::mcp2hits(mcp_michel, vp_hit, clockData, true), "p", ms_mcp_me);
                }
            }
        }
    }
            
    // loop over tracks to find muons
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        // bool tooSmall = p_trk->Length() < fTrackLengthCut;

        HitPtrVec vp_hit_muon = fmp_trk2hit.at(p_trk.key());
        if (!LOG(vp_hit_muon.size())) continue;

        HitPtrPair cathode_crossing;
        HitPtrVec tpc_crossing;
        HitPtrVec vp_sorted_hit;

        HitPtrPair trk_ends =
            GetTrackEndsHits(
                vp_hit_muon,
                &cathode_crossing,
                &tpc_crossing,
                fSort ? &vp_sorted_hit : nullptr
            );


        // no reconstructed track ends
        if (!LOG(trk_ends.first && trk_ends.second)) {
            std::vector<TGraph*> gs(ana::n_sec[geoDet]);
            for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
                gs[s] = new TGraph();
                gs[s]->SetName(Form("g%u_%u", p_trk->ID(), s));
                gs[s]->SetTitle(Form("track %u, section %u", p_trk->ID(), s));
                SetLineStyle(gs[s], ls_no_end);
            }
            continue;
        }

        // fiducial tests
        bool outsideFront = !wireWindow.isInside(trk_ends.first->PeakTime() * fTick2cm, fMichelTickRadius)
            || !geoHighX.InFiducialZ(GetSpace(trk_ends.first->WireID()), fMichelSpaceRadius);
        bool outsideBack = !wireWindow.isInside(trk_ends.second->PeakTime() * fTick2cm, fMichelTickRadius)
            || !geoHighX.InFiducialZ(GetSpace(trk_ends.second->WireID()), fMichelSpaceRadius);

        // broken tests
        bool isBroken = false;
        simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        if (mcp) {
            std::vector<art::Ptr<recob::Track>> vp_trk_from_mcp = mcp2trks(mcp, vp_trk, clockData, fmp_trk2hit);
            unsigned long_trk_count = 0;
            for (art::Ptr<recob::Track> const& p_trk_from_mcp : vp_trk_from_mcp)
                if (p_trk_from_mcp->Length() > fTrackLengthCut)
                    long_trk_count++;
            isBroken = long_trk_count > 1;

            // if (vp_trk_from_mcp.size()) {
            //     if (vp_trk_from_mcp.size() == 1)
                //     MuonTrackIsNotBroken = kNotBroken;
                // else {
                //     bool IsDeepestTrack = true;
                //     for (art::Ptr<recob::Track> p_trk_from_mcp : vp_trk_from_mcp)
                //         if (geoDet == kPDVD)
                //             IsDeepestTrack = IsDeepestTrack
                //                 && (MuonEndTrackPoint.x <= p_trk_from_mcp->Start().X()
                //                 && MuonEndTrackPoint.x <= p_trk_from_mcp->End().X());
                //         else if (geoDet == kPDHD)
                //             IsDeepestTrack = IsDeepestTrack
                //                 && (MuonEndTrackPoint.y <= p_trk_from_mcp->Start().Y()
                //                 && MuonEndTrackPoint.y <= p_trk_from_mcp->End().Y());
                //     if (IsDeepestTrack)
                //         MuonTrackIsNotBroken = kLastOfBroken;
                //     else
                //         MuonTrackIsNotBroken = kBroken;
                // }
            // }
        }

        // track hits and regression
        if (fSort)
            drawGraph(vp_sorted_hit, "l", {}, isBroken ? ls_broken : ls_mu);
        else
            drawGraph(vp_hit_muon, "l", {}, isBroken ? ls_broken : ls_mu);

        // std::vector<TGraph*> gs(ana::n_sec[geoDet]);
        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            ana::LinearRegression reg;
            // LinearRegression const& reg = s >= ana::n_sec[geoDet]/2 ?
            //     per_side_reg.second : per_side_reg.first;
            // bool lin = reg.r2() > 0.5;

            // TGraph* g = gs[s] = new TGraph();
            // g->SetName(Form("g%u_%u", p_trk->ID(), s));
            // g->SetTitle(Form("track %u, section %u", p_trk->ID(), s));
            // g->SetLineWidth(1);
            // g->SetLineColor(tooSmall ? theme.c_small : theme.c_mu);
            // if (fSort) {
            //     for (HitPtr const& p_hit : vvp_sec_sorted[s]) {
            //         double const z = GetSpace(p_hit->WireID());
            //         double const t = p_hit->PeakTime() * fTick2cm;
            //         reg.add(z, t);
            //         if (geoDet == kPDVD)
            //             g->AddPoint(z, t);
            //         else if (geoDet == kPDHD)
            //             g->AddPoint(t, z);
            //     }
            // } else {
            //     for (HitPtr const& p_hit : vp_hit_muon) {
            //         if (p_hit->View() != geo::kW) continue;
            //         if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != int(s)) continue;
            //         double const z = GetSpace(p_hit->WireID());
            //         double const t = p_hit->PeakTime() * fTick2cm;
            //         reg.add(z, t);
            //         if (geoDet == kPDVD)
            //             g->AddPoint(z, t);
            //         else if (geoDet == kPDHD)
            //             g->AddPoint(t, z);
            //     }
            // }
            // reg.normalize();

            // c->cd(s+1);
            // if (g->GetN()) g->Draw("l");

            for (HitPtr const& p_hit : vp_hit_muon) {
                if (p_hit->View() != geo::kW) continue;
                if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != int(s)) continue;
                reg.add(GetSpace(p_hit->WireID()), p_hit->PeakTime() * fTick2cm);
            }
            reg.compute();

            TF1* f = new TF1();
            if (geoDet == kPDVD) {
                double xmin = reg.m > 1 ? reg.mx - 10 : reg.m * (reg.my - 10) + reg.m;
                double xmax = reg.m > 1 ? reg.mx + 10 : reg.m * (reg.my + 10) + reg.m;
                f = new TF1(
                    Form("reg%u_s%u", p_trk->ID(), s),
                    "(x - [1]) / [0]",
                    xmin, xmax
                );
            } else if (geoDet == kPDHD) {
                double xmin = reg.m > 1 ? (reg.mx - 10 - reg.m) / reg.m : reg.my - 10;
                double xmax = reg.m > 1 ? (reg.mx + 10 - reg.m) / reg.m : reg.my + 10;
                f = new TF1(
                    Form("reg%u_s%u", p_trk->ID(), s),
                    "[0]*x + [1]",
                    xmin, xmax
                );
            }
            f->SetParameter(0, reg.m);
            f->SetParameter(1, reg.m);
            SetLineStyle(f, reg.r2 > 0.5 ? ls_reg_good : ls_reg_bad);

            c->cd(s+1);
            f->Draw("same");
        }

        // track ends and cathode/section crossing
        if (cathode_crossing.first.isNonnull()) {
            drawMarker(trk_ends.first, outsideFront ? ms_cc_in : ms_cc_in);
            drawMarker(trk_ends.second, outsideBack ? ms_cc_out : ms_cc_in);

            drawMarker(cathode_crossing.first, ms_cc);
            drawMarker(cathode_crossing.second, ms_cc);
        } else {
            drawMarker(trk_ends.first, outsideFront ? ms_mu_in : ms_mu_in);
            drawMarker(trk_ends.second, outsideBack ? ms_mu_out : ms_mu_in);
        }
        if (tpc_crossing.size()) {
            for (HitPtr const& p_hit : tpc_crossing)
                drawMarker(p_hit, ms_sc);
        }

        // simb::MCParticle const* mcp = ana::trk2mcp(p_trk, clockData, fmp_trk2hit);
        // if (mcp && abs(mcp->PdgCode()) == 13) {
        //     HitPtrVec vp_mcp_hit = ana::mcp2hits(
        //         mcp, vp_hit, clockData, false
        //     );
        //     HitPtrPair mcp_ends = GetTrackEndsHits(vp_mcp_hit);
        //     if (mcp_ends.first && mcp_ends.second) {
        //         int dir_z = mcp->EndZ() > mcp->Vz() ? 1 : -1;
        //         float fz = GetSpace(mcp_ends.first->WireID());
        //         float sz = GetSpace(mcp_ends.second->WireID());
        //         HitPtr mcp_end = (sz-fz) * dir_z > 0 ? mcp_ends.second : mcp_ends.first;

        //         m->SetMarkerColor(theme.c_mcp_mu);
        //         m->SetMarkerStyle(theme.m_mcp_mu);
        //         m->SetMarkerSize(theme.s_mcp_mu);
        //         drawMarker(m, mcp_end);
        //     } 

        //     TGraph* g_mcp_mu = new TGraph();
        //     g_mcp_mu->SetMarkerColor(theme.c_mcp_mu);
        //     g_mcp_mu->SetMarkerStyle(theme.m_mcp_mu);
        //     g_mcp_mu->SetMarkerSize(theme.s_mcp_mu);
        //     drawGraph(g_mcp_mu, vp_mcp_hit, "p");

        //     if (mcp->NumberDaughters() >= 3) {
        //         simb::MCParticle const * mcp_michel = nullptr;
        //         bool has_numu = false, has_nue = false;
        //         for (int i_dau=mcp->NumberDaughters()-3; i_dau<mcp->NumberDaughters(); i_dau++) {
        //             simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    
        //             if (!mcp_dau) continue;

        //             switch (abs(mcp_dau->PdgCode())) {
        //                 case 14: has_numu = true; break;
        //                 case 12: has_nue = true; break;
        //                 case 11: mcp_michel = mcp_dau; break;
        //                 default: break;
        //             }
        //         }
        //         if (mcp_michel and has_numu and has_nue) {
        //             HitPtrVec vp_mcp_me = ana::mcp2hits(
        //                 mcp_michel, vp_hit, clockData, false
        //             );

        //             TGraph* g_mcp_me = new TGraph();
        //             g_mcp_me->SetMarkerColor(theme.c_mcp_me);
        //             g_mcp_me->SetMarkerStyle(theme.m_mcp_me);
        //             g_mcp_me->SetMarkerSize(theme.s_mcp_me);
        //             drawGraph(g_mcp_me, vp_mcp_me, "p");
        //         }
        //     }
        // }



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

        using HitPtr = HitPtr;
        using HitPtrVec = HitPtrVec;
        using HitPtrPair = HitPtrPair;

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
            f->SetParameter(0, reg.m / fTick2cm);
            f->SetParameter(1, reg.m);
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
    } // end of loop over tracks

    for (art::Ptr<recob::Shower> const& p_shw : vp_shw) {
        HitPtrVec vp_hit_shw = fmp_shw2hit.at(p_shw.key());
        if (vp_hit_shw.empty()) continue;
        drawGraph(vp_hit_shw, "p", ms_shw);
    }
        
    c->Write();
}

void ana::Trackchecks::endJob() {}

double ana::Trackchecks::GetSpace(geo::WireID wid) {
    return plane2axis[wid].space(asWire->Wire(wid));
}

ana::Hit ana::Trackchecks::GetHit(HitPtr const p_hit) {
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

HitPtrPair ana::Trackchecks::GetTrackEndsHits(
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
    if (side_reg[0].n < nmin && side_reg[1].n < nmin) return {};

    // compute average from sum
    for (ana::LinearRegression& reg : side_reg) reg.compute();

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
    // THIS CAUSES A SEGFAULT FOR SOME REASON???
    // auto side_sort = [&](int side) {
    //     return [&](HitPtr const& h1, HitPtr const& h2) -> bool {
    //         double const s1 = side_reg[side].projection(
    //             GetSpace(h1->WireID()),
    //             h1->PeakTime() * fTick2cm
    //         );
    //         double const s2 = side_reg[side].projection(
    //             GetSpace(h2->WireID()),
    //             h2->PeakTime() * fTick2cm
    //         );
    //         return s1 < s2;
    //     };
    // };

    if (vp_sorted_hit) {
        vp_sorted_hit->clear();

        for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
            HitPtrVec& vp_sec_sorted = vp_sec_hit[s];
            if (vp_sec_sorted.size() < nmin) continue;

            int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
            std::sort(
                vp_sec_sorted.begin(), 
                vp_sec_sorted.end(),
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
                [&, &reg=side_reg[side]](HitPtr const& h1, HitPtr const& h2) -> bool {
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