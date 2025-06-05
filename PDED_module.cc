////////////////////////////////////////////////////////////////////////
// Class:       PDED
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDHDED_module.cc
//
// Generated at Fri Feb 21 03:35:05 2025 by Jeremy Quelin Lechevranton using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "utils.h"
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TColor.h>
#include <TText.h>

namespace ana {
    class PDED;
}

class ana::PDED : public art::EDAnalyzer {
public:
    explicit PDED(fhicl::ParameterSet const& p);
    PDED(PDED const&) = delete;
    PDED(PDED&&) = delete;
    PDED& operator=(PDED const&) = delete;
    PDED& operator=(PDED&&) = delete;

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

    // art::ServiceHandle<cheat::ParticleIn_secentoryService> pi_serv;
    // art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    int geoDet;
    enum EnumDet { kPDVD, kPDHD };

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
    std::map<unsigned, int> tpc2sec;
    std::map<int, std::pair<unsigned, unsigned>> sec2tpc;

    // Data Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag tag_mcp, tag_sed, tag_wir,
        tag_hit, tag_clu, tag_trk,
        tag_spt, tag_pfp, tag_r3d;


    unsigned cn=0;


    // std::vector<TCanvas*> cs;
    double GetSpace(geo::WireID);
    // void drawGraph(std::vector<TPad*>, std::vector<art::Ptr<recob::Hit>>, char const* draw_opt, int color = kBlack, int style = kFullCircle, float size = .3, float width = .3);
};


ana::PDED::PDED(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products"))
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

    wireWindow = bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    geoLowX = geo::BoxBoundedGeo{
        asGeo->TPC(geo::TPCID{0, 1}).Min(),
        asGeo->TPC(geo::TPCID{0, 5}).Max()
    };
    geoHighX = geo::BoxBoundedGeo{
        asGeo->TPC(geo::TPCID{0, 2}).Min(),
        asGeo->TPC(geo::TPCID{0, 6}).Max()
    };

    if (geoDet == kPDVD) {
        tpc2sec = {
            {0, 0}, {2, 0},
            {1, 1}, {3, 1},
            {4, 2}, {6, 2},
            {5, 3}, {7, 3},
            {8, 4}, {10, 4},
            {9, 5}, {11, 5},
            {12, 6}, {14, 6},
            {13, 7}, {15, 7}
        };
        sec2tpc = {
            {0, {0, 2}},
            {1, {1, 3}},
            {2, {4, 6}},
            {3, {5, 7}},
            {4, {8, 10}},
            {5, {9, 11}},
            {6, {12, 14}},
            {7, {13, 15}}
        };
    } else if (geoDet == kPDHD) {
        tpc2sec = {
            {4, -1}, {0, -1},
            {5, 0}, {1, 0},
            {6, 1}, {2, 1},
            {7, -1}, {3, -1}
        };
        sec2tpc = {
            {0, {1, 5}},
            {1, {2, 6}}
        };
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
        << "  (" << asGeo->DetectorName() << ")" << std::endl
        << "  Sampling Rate: " << fSamplingRate << " µs/tick" << std::endl
        << "  Drift Velocity: " << fDriftVelocity << " cm/µs" << std::endl
        << "  Channel Pitch: " << fChannelPitch << " cm" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX.Min() << " -> " << geoHighX.Max() << std::endl
        << "  LowX Bounds: " << geoLowX.Min() << " -> " << geoLowX.Max() << std::endl;
    
}

void ana::PDED::analyze(art::Event const& e) {
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

    auto const & vh_spt = e.getHandle<std::vector<recob::SpacePoint>>(tag_spt);
    if (!vh_spt.isValid()) return;
    std::vector<art::Ptr<recob::SpacePoint>> vp_spt;
    art::fill_ptr_vector(vp_spt, vh_spt);

    auto const & vh_pfp = e.getHandle<std::vector<recob::PFParticle>>(tag_pfp);
    if (!vh_pfp.isValid()) return;

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    art::FindOneP<recob::SpacePoint> fop_hit2spt(vh_hit, e, tag_spt);
    art::FindOneP<recob::Hit> fop_spt2hit(vh_spt, e, tag_spt);

    art::FindOneP<recob::PFParticle> fop_trk2pfp(vh_trk, e, tag_trk);
    art::FindManyP<recob::SpacePoint> fmp_pfp2spt(vh_pfp, e, tag_pfp);

    // std::cout << "e" << e.event() << "\r" << std::flush;

    // TCanvas *c = new TCanvas(
    //     Form("e%u", e.event()), 
    //     Form("e%u", e.event())
    // );

    TCanvas *c = tfs->make<TCanvas>(
        Form("c%u", cn++),
        Form("run:%u, subrun:%u, event:%u", e.run(), e.subRun(), e.event()),
        1300,800
    );
    // cs.push_back(c);


    // unsigned const n_section = 2;
    // bool const axis_label = true;
    // int const font = 43; // 80 = courrier, 3 = font size in pixel
    // float const font_size = 20;

    Style_t const font = 43;
    struct binning {
        unsigned n;
        float min, max;
        binning(unsigned n, float m, float M) : n(n), min(m), max(M) {}
        binning(float m, float M, float s) : n((M-m)/s), min(m), max(M) {}
        float step(void) const { return (max-min)/n; }
    }   b_z{(float) geoHighX.MinZ(), (float) geoHighX.MaxZ(), fChannelPitch},
        b_t{wireWindow.min, wireWindow.max, fChannelPitch / fTick2cm};

    unsigned n_sec = 0;
    std::function<TH2F(unsigned)> frame;
    Style_t font_size;
    struct { Float_t l, r, b, t; } pad_margin;
    Float_t title_offset_x, title_offset_y;
    std::string X, Y;
    if (geoDet == kPDVD) {
        c->Divide(4, 2);
        n_sec = 8;
        frame = [b_z, b_t](unsigned s) -> TH2F {
            return TH2F(Form("f%u", s), ";Z;T",
                    b_z.n, b_z.min, b_z.max,
                    b_t.n, b_t.min, b_t.max);
        };
        font_size = 12;
        pad_margin = {0.14, 0.04, 0.09, 0.06};
        title_offset_x = 1.3;
        title_offset_y = 1.7;
    } else if (geoDet == kPDHD) {
        c->Divide(2, 1);
        n_sec = 2;
        frame = [b_z, b_t](unsigned s) -> TH2F {
            return TH2F(Form("f%u", s), ";T;Z",
                    b_t.n, b_t.min, b_t.max,
                    b_z.n, b_z.min, b_z.max);
        };
        font_size = 20;
        pad_margin = {0.1, 0.04, 0.09, 0.06};
        title_offset_x = 1.5;
        title_offset_y = 1.5;
    }
        
    for (unsigned s=0; s<n_sec; s++) {
        c->cd(s+1);
        gPad->SetMargin(
            pad_margin.l, pad_margin.r,
            pad_margin.b, pad_margin.t
        );
        gPad->SetTicks(1, 1);
        TH2F* f = new TH2F(frame(s));
        f->SetStats(kFALSE);
        f->SetTitleFont(font, "xyz");
        f->SetLabelFont(font, "xyz");
        f->SetTitleSize(font_size, "xyz");
        f->SetLabelSize(font_size, "xyz");
        f->SetTitleOffset(title_offset_x, "x");
        f->SetTitleOffset(title_offset_y, "y");
        for (TAxis* ax : {f->GetXaxis(), f->GetYaxis()}) ax->CenterTitle();
        // f->Draw(s ? "rx" : "");
        f->Draw();

        TText* t = new TText(
            gPad->GetLeftMargin(),
            1-gPad->GetTopMargin()+0.01,
            Form("TPC %u & %u", sec2tpc[s].first, sec2tpc[s].second)
        );
        t->SetNDC();
        t->SetTextFont(font);
        t->SetTextSize(font_size);
        t->SetTextAlign(kHAlignLeft + kVAlignBottom);
        t->Draw();

        if ((geoDet == kPDHD && s==1) || (geoDet == kPDVD && s==3)) {
            TText* tt = new TText(
                1-gPad->GetRightMargin(),
                1-gPad->GetTopMargin()+0.01,
                Form("R:%u-SR:%u-E:%u", e.run(), e.subRun(), e.event())
            );
            tt->SetNDC();
            tt->SetTextFont(103); 
            tt->SetTextSize(font_size);
            tt->SetTextAlign(kHAlignRight + kVAlignBottom);
            tt->Draw();
        }
    }

    gStyle->SetPalette(kCividis);
    TArrayI const& colors = TColor::GetPalette();

    struct range {
        float min, max;
        range(float m, float M) : min(m), max(M) {}
        float normalize(float x) const {
            if (x <= min) return 0.F;
            if (x >= max) return 1.F;
            return (x-min) / (max-min);
        }
    } r_adc{0.F, 300.F};

    TMarker* m = new TMarker();
    m->SetMarkerStyle(kFullCircle);
    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;
        int s = tpc2sec[p_hit->WireID().TPC];
        if (s == -1) continue;
        float const x = r_adc.normalize(p_hit->Integral());
        m->SetMarkerSize(2*x+0.1);
        m->SetMarkerColor(colors[int((colors.GetSize()-1)*x)]);

        c->cd(s+1);
        if (geoDet == kPDVD)
            m->DrawMarker(GetSpace(p_hit->WireID()), p_hit->PeakTime());
        else if (geoDet == kPDHD)
            m->DrawMarker(p_hit->PeakTime(), GetSpace(p_hit->WireID()));
    }

    // unsigned g=0;
    // for (art::Ptr<recob::Track> p_trk : vp_trk) {
    //     // if (p_trk->Length() > 40) continue;
    //     std::vector<art::Ptr<recob::Hit>> vp_hit_from_trk = fmp_trk2hit.at(p_trk.key());
    //     drawGraph(ps, vp_hit_from_trk, "l", kOrange-10+g, 0, 0, 1);
    //     g++;
    // }

    c->Write();
}

void ana::PDED::beginJob() {}
void ana::PDED::endJob() {}

double ana::PDED::GetSpace(geo::WireID wid) {
    return plane2axis[(geo::PlaneID) wid].space(asWire->Wire(wid));
}

// void ana::PDED::drawGraph(std::vector<TPad*> ps, std::vector<art::Ptr<recob::Hit>> vp_hit, char const* draw_opt, int color = kBlack, int style = kFullCircle, float size = .3, float width = .3) {
//     std::map<unsigned, int> static tpc2sec = {
//         {0, -1},
//         {1, 0},
//         {2, 1},
//         {3, -1},
//         {4, -1},
//         {5, 0},
//         {6, 1},
//         {7, -1}
//     };
//     unsigned static ng = 0;
//     ng++;
//     std::vector<TGraph*> gs{ps.size()};
//     for (unsigned s=0; s<ps.size(); s++) {
//         TGraph* g = gs[s] = new TGraph();
//         g->SetName(Form("g%u_%u", ng, s));
//         g->SetEditable(kFALSE);
//         g->SetMarkerColor(color);
//         g->SetMarkerStyle(style);
//         g->SetMarkerSize(size);
//         g->SetLineColor(color);
//         g->SetLineWidth(width);
//     }
//     for (art::Ptr<recob::Hit> p_hit : vp_hit) {
//         if (p_hit->View() != geo::kW) continue;
//         int s = tpc2sec[p_hit->WireID().TPC];
//         if (s == -1) continue;
//         gs[s]->AddPoint(p_hit->PeakTime(), GetSpace(p_hit->WireID()));
//     }
//     for (unsigned s=0; s<ps.size(); s++) {
//         ps[s]->cd();
//         if (!gs[s]->GetN()) continue;
//         gs[s]->Draw(draw_opt);
//     }
// }

DEFINE_ART_MODULE(ana::PDED)