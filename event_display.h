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
    inline void drawFrame(TCanvas* c, int geoDet, unsigned r=0, unsigned sr=0, unsigned e=0, int real=-1) {}
    inline void drawFrame(TCanvas* c, int geoDet, char const* left_title="", char const* right_title="") {
        Style_t const font = 43;
        unsigned n_sec = 0;
        Style_t font_size;
        struct { Float_t l, r, b, t; } pad_margin;
        Float_t title_offset_x, title_offset_y;

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
    // inline void drawEventHits(TCanvas* c, int geoDet, std::vector<art::Ptr<recob::Hit>> vp_hit, std::function<double(geo::WireID)> GetSpace) {
    //     float const max_adc = geoDet == kPDVD ? 1000 : 200;
    //     gStyle->SetPalette(kCividis);
    //     TArrayI const& colors = TColor::GetPalette();
    //     TMarker* m = new TMarker();
    //     m->SetMarkerStyle(kFullCircle);
    //     for (art::Ptr<recob::Hit> p_hit : vp_hit) {
    //         if (p_hit->View() != geo::kW) continue;
    //         int s = tpc2sec[geoDet][p_hit->WireID().TPC];
    //         if (s == -1) continue;
    //         float const x = p_hit->Integral() / max_adc;
    //         m->SetMarkerSize(2*x+0.1);
    //         m->SetMarkerColor(colors[int((colors.GetSize()-1)*x)]);

    //         c->cd(s+1);
    //         if (geoDet == kPDVD)
    //             m->DrawMarker(GetSpace(p_hit->WireID()), p_hit->PeakTime());
    //         else if (geoDet == kPDHD)
    //             m->DrawMarker(p_hit->PeakTime(), GetSpace(p_hit->WireID()));
    //     }
    // }
}