#include "utils.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMarker.h>
#include <TLine.h>
#include <TStyle.h>
#include <TColor.h>
#include <TText.h>

#include <TGraph2D.h>

namespace ana {
    inline void DrawFrame(TCanvas* c, int geoDet, char const* left_title="", char const* right_title="") {
        Style_t const font = 43;
        unsigned n_sec=0;
        Style_t font_size=0;
        struct { Float_t l, r, b, t; } pad_margin={};
        Float_t title_offset_x=0, title_offset_y=0;

        c->SetMargin(0, 0, 0, 0);
        if (geoDet == kPDVD) {
            c->Divide(4, 2, 0, 0);
            n_sec = 8;
            font_size = 12;
            pad_margin = {0.14, 0.04, 0.09, 0.06};
            title_offset_x = 1.3;
            title_offset_y = 1.7;
        } else if (geoDet == kPDHD) {
            c->Divide(2, 1, 0, 0);
            n_sec = 2;
            font_size = 20;
            pad_margin = {0.1, 0.04, 0.09, 0.06};
            title_offset_x = 1.5;
            title_offset_y = 1.5;
        }

        TText* t = new TText();
        t->SetNDC();
        t->SetTextFont(font);
        t->SetTextSize(font_size);
        t->SetTextAlign(kHAlignLeft + kVAlignBottom);

        TText* title = new TText();
        title->SetNDC();
        title->SetTextFont(font);
        title->SetTextSize(font_size);
        title->SetTextAlign(kHAlignRight + kVAlignBottom);

        for (unsigned s=0; s<n_sec; s++) {
            c->cd(s+1);
            gPad->SetMargin(
                pad_margin.l, pad_margin.r,
                pad_margin.b, pad_margin.t
            );
            gPad->SetTicks(1, 1);
            TH2F* f = new TH2F();  // -Werror=maybe-uninitialized
            if (geoDet == kPDVD) {
                f = new TH2F(Form("f%u", s), ";Z (cm);T (cm)",
                    600, 0, 300,
                    // 600, 0, 6000
                    600, 0, 480
                );
            } else if (geoDet == kPDHD) {
                f = new TH2F(Form("f%u", s), ";T (cm);Z (cm)",
                    // 600, 0, 6000,
                    600, 0, 480,
                    600, 0, 464
                );
            }
            f->SetStats(kFALSE);
            f->SetTitleFont(font, "xyz");
            f->SetLabelFont(font, "xyz");
            f->SetTitleSize(font_size, "xyz");
            f->SetLabelSize(font_size, "xyz");
            f->SetTitleOffset(title_offset_x, "x");
            f->SetTitleOffset(title_offset_y, "y");
            for (TAxis* ax : {f->GetXaxis(), f->GetYaxis()}) ax->CenterTitle();
            f->Draw();

            t->DrawText(
                gPad->GetLeftMargin(),
                1-gPad->GetTopMargin()+0.01,
                Form("TPC %u & %u", sec2tpc[geoDet][s].first, sec2tpc[geoDet][s].second)
            );

            if (s == 0) {
                title->DrawText(
                    1-gPad->GetRightMargin(),
                    1-gPad->GetTopMargin()+0.01,
                    left_title
                );
            }
            if ((geoDet == kPDHD && s==1) || (geoDet == kPDVD && s==3)) {
                title->DrawText(
                    1-gPad->GetRightMargin(),
                    1-gPad->GetTopMargin()+0.01,
                    right_title
                );
            }
        }
    }

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
    void inline SetMarkerStyle(TAttMarker* m, MarkerStyle const& ms) {
        m->SetMarkerColor(ms.c);
        m->SetMarkerStyle(ms.m);
        m->SetMarkerSize(ms.s);
    };
    void inline SetLineStyle(TAttLine* l, LineStyle const& ls) {
        l->SetLineColor(ls.c);
        l->SetLineStyle(ls.l);
        l->SetLineWidth(ls.w);
    };
    class MichelDisplayer : public MichelAnalyzer {
    public:
        void DrawMarker(TCanvas* c, PtrHit const& ph, MarkerStyle const& ms) const;
        void DrawGraph(TCanvas* c, VecPtrHit const& vph, char const* draw, MarkerStyle const& ms={}, LineStyle const& ls={}) const;
        void DrawGraph2D(TCanvas* c, PtrTrk const& pt, MarkerStyle const& ms={}, LineStyle const& ls={}) const; 
    };

    void MichelDisplayer::DrawMarker(TCanvas* c, PtrHit const& ph, MarkerStyle const& ms) const {
        TMarker *m = new TMarker();
        SetMarkerStyle(m, ms);
        if (geoDet == kPDVD) {
            int s = tpc2sec[geoDet][ph->WireID().TPC];
            c->cd(s+1);
            m->DrawMarker(GetSpace(ph->WireID()), ph->PeakTime() * fTick2cm);
        } else if (geoDet == kPDHD) {
            int s = tpc2sec[geoDet][ph->WireID().TPC];
            if (s == -1) return;
            c->cd(s+1);
            m->DrawMarker(ph->PeakTime() * fTick2cm, GetSpace(ph->WireID()));
        }
    }
    void MichelDisplayer::DrawGraph(TCanvas* c, VecPtrHit const& vph, char const* draw, MarkerStyle const& ms, LineStyle const& ls) const {
        unsigned static gn=0;
        std::vector<TGraph*> gs(n_sec[geoDet]);
        for (unsigned s=0; s<n_sec[geoDet]; s++) {
            gs[s] = new TGraph();
            gs[s]->SetEditable(kFALSE);
            gs[s]->SetName(Form("g%u_s%u", gn, s));
            SetMarkerStyle(gs[s], ms);
            SetLineStyle(gs[s], ls);
        }
        for (HitPtr p_hit : vph) {
            if (p_hit->View() != geo::kW) continue;
            int s = tpc2sec[geoDet][p_hit->WireID().TPC];
            if (s == -1) continue;
            if (geoDet == kPDVD)
                gs[s]->AddPoint(GetSpace(p_hit->WireID()), p_hit->PeakTime() * fTick2cm);
            else if (geoDet == kPDHD)
                gs[s]->AddPoint(p_hit->PeakTime() * fTick2cm, GetSpace(p_hit->WireID()));
        }
        for (unsigned s=0; s<n_sec[geoDet]; s++) {
            if (!gs[s]->GetN()) continue;
            c->cd(s+1);
            gs[s]->Draw(draw);
        }
        gn++;
    }
    void MichelDisplayer::DrawGraph2D(TCanvas* c, PtrTrk const& pt, MarkerStyle const& ms, LineStyle const& ls) const {
        TGraph2D* g = new TGraph2D();
        SetMarkerStyle(g, ms);
        SetLineStyle(g, ls);
        for (unsigned it=0; it<pt->NumberTrajectoryPoints(); it++) {
            if (!pt->HasValidPoint(it)) continue;
            geo::Point_t p = pt->LocationAtPoint(it);
            if (geoDet == kPDVD)
                g->AddPoint(p.Y(), p.Z(), p.X());
            else if (geoDet == kPDHD)
                g->AddPoint(p.Z(), p.X(), p.Y());
        }
        c->cd();
        g->Draw("same line");
    }
}