////////////////////////////////////////////////////////////////////////
// Class:       PDHDED
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

namespace ana {
    class PDHDED;
}

class ana::PDHDED : public art::EDAnalyzer {
public:
    explicit PDHDED(fhicl::ParameterSet const& p);
    PDHDED(PDHDED const&) = delete;
    PDHDED(PDHDED&&) = delete;
    PDHDED& operator=(PDHDED const&) = delete;
    PDHDED& operator=(PDHDED&&) = delete;

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
        tag_spt, tag_pfp, tag_r3d;


    unsigned cn=0;


    // std::vector<TCanvas*> cs;
    double GetSpace(geo::WireID);
};


ana::PDHDED::PDHDED(fhicl::ParameterSet const& p)
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
        tag_r3d = art::InputTag{"reco3d", ""};
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

void ana::PDHDED::analyze(art::Event const& e) {
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
        Form("run:%u, subrun:%u, event:%u", e.run(), e.subRun(), e.event())
    );
    // cs.push_back(c);

    std::map<unsigned, int> map_tpc_sec = {
        {0, -1},
        {1, 0},
        {2, 1},
        {3, -1},
        {4, -1},
        {5, 0},
        {6, 1},
        {7, -1}
    };
    unsigned const n_section = 2;
    bool const axis_label = true;
    int const font = 43; // 80 = courrier, 3 = font size in pixel
    float const font_size = 20;

    struct binning {
        unsigned n;
        float min, max;
        binning(unsigned n, float m, float M) : n(n), min(m), max(M) {}
        binning(float m, float M, float s) : n((M-m)/s), min(m), max(M) {}
        float step(void) const { return (max-min)/n; }
    }   binZ{(float) geoHighX.MinZ(), (float) geoHighX.MaxZ(), fChannelPitch},
        binT{wireWindow.min, wireWindow.max, fChannelPitch / fTick2cm};

    std::vector<TPad*> ps{n_section};
    for (unsigned s=0; s<n_section; s++) {
        TPad* p = ps[s] = new TPad(
            Form("pad_s%u", s),
            Form("pad_s%u", s),
            0.5*s, 0.,
            0.5*(s+1), 1.
        );

        // double pad_margin = axis_label ? 0.1 : 0.05;
        // p->SetMargin(pad_margin, pad_margin, pad_margin, pad_margin);
        // if (!axis_label) {
        //     if (s) p->SetLeftMargin(pad_margin / 2);
        //     else p->SetRightMargin(pad_margin / 2);
        // }
        // if (axis_label) p->SetTicks(1, 1); // ticks on top and right

        TH2F* f = new TH2F(
            Form("frame_s%u", s),
            axis_label ? ";T;Z" : "", 
            binT.n, binT.min, binT.max,
            binZ.n, binZ.min, binZ.max
        );
        f->SetStats(kFALSE);
        f->SetDirectory(nullptr);
        // f->SetTitleSize(0, "xyz"); // no title for all axis

        // if (axis_label) {
        //     f->SetTitleFont(font);
        //     f->SetLabelFont(font);
        //     f->SetTitleSize(font_size, "xy");
        //     f->SetLabelSize(font_size, "xy");
        //     f->SetTitleOffset(1.5, "xy");
        // } else f->SetTitleOffset(0.4, "xy");
        // for (TAxis *ax : {f->GetXaxis(), f->GetYaxis()}) {
        //     ax->CenterTitle();
        //     if (!axis_label) ax->SetNdivisions(0); // remove ticks and label
        // }

        c->cd();
        p->Draw();
        p->cd();
        // f->Draw(s ? "rx" : "");
        f->Draw();
    }

    Color_t const color = kGray;
    Style_t const style = kFullSquare;
    Size_t const size = 0.3;
    Width_t const width = 0.3;
    char const * draw_opt = "p";
    std::vector<TGraph*> gs{n_section};
    for (unsigned s=0; s<n_section; s++) {
        TGraph* g = gs[s] = new TGraph();
        g->SetName(Form("g%u_%u", 0, s));
        g->SetEditable(kFALSE);
        g->SetMarkerColor(color);
        g->SetMarkerStyle(style);
        g->SetMarkerSize(size);
        g->SetLineColor(color);
        g->SetLineWidth(width);
    }
    for (art::Ptr<recob::Hit> p_hit : vp_hit) {
        if (p_hit->View() != geo::kW) continue;
        int s = map_tpc_sec[p_hit->WireID().TPC];
        if (s == -1) continue;
        gs[s]->AddPoint(p_hit->PeakTime(), GetSpace(p_hit->WireID()));
    }
    for (unsigned s=0; s<n_section; s++) {
        ps[s]->cd();
        if (!gs[s]->GetN()) continue;
        gs[s]->Draw(draw_opt);
    }

    // c->Draw();
    // std::vector<TPad*> ps = drawFrame(c);
    // drawGraph(ps, EventHits, "p", kGray, kFullSquare, .3);
    // // drawHisto(ps, EventHits);

    // unsigned o=0;
    // for (unsigned mu=0; mu<EventNMuon; mu++) {
    //     muon->GetEntry(EventiMuon->at(mu));
    //     // if (!MuonHasMichel) continue;

    //     // -- Only Cathode Crossing
    //     // bool left=false, right=false;
    //     // for (Hit hit : MuonHits) {
    //     //     if (hit.tpc == 1 || hit.tpc == 5)
    //     //         left = true;
    //     //     if (hit.tpc == 2 || hit.tpc == 6)
    //     //         right = true;
    //     //     if (left && right) break;
    //     // }
    //     // if (!left || !right) continue;

    //     drawGraph(ps, MuonHits, "l", kOrange-4+o++, 0, 0., 1.);
    //     drawMarker(ps, TrueEnd, kC6Pistachio, 89, 3);
    //     drawMarker(ps, MuonEndHit, kC6Cyclamen, kFullCrossX, 2);
    //     drawGraph(ps, MichelHits, "p", kC6Blue);
    // }

    c->Write();

    // std::cout << "next event? (enter)"; 
    // std::string input;
    // std::getline(std::cin, input);
}

void ana::PDHDED::beginJob() {}
void ana::PDHDED::endJob() {}

double ana::PDHDED::GetSpace(geo::WireID wid) {
    return plane2axis[(geo::PlaneID) wid].space(asWire->Wire(wid));
}

DEFINE_ART_MODULE(ana::PDHDED)