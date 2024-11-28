////////////////////////////////////////////////////////////////////////
// Class:       Michecks
// Plugin Type: analyzer
// File:        Michecks_module.cc
//
// Generated at Wed Nov  6 10:31:28 2024 by Jeremy Quelin Lechevranton
////////////////////////////////////////////////////////////////////////

// Art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"

#include "larcorealg/Geometry/Exceptions.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"

// ProtoDUNE includes
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

// ROOT includes
#include "TLorentzVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEllipse.h"

// std includes
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::min;
using std::max;


namespace ana {
  class Michecks;
  struct Binning {
        int n;
        double min, max;
  };
}


class ana::Michecks : public art::EDAnalyzer {
public:
    explicit Michecks(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Michecks(Michecks const&) = delete;
    Michecks(Michecks&&) = delete;
    Michecks& operator=(Michecks const&) = delete;
    Michecks& operator=(Michecks&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

private:

    // Utilities
    art::ServiceHandle<art::TFileService> tfs;

    const geo::Geometry* asGeo;
    const detinfo::DetectorPropertiesService* asDetProp;
    const detinfo::DetectorClocksService* asDetClocks;

    protoana::ProtoDUNETruthUtils truthUtil;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Conversion factors
    float fADCtoE = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    float fChannelPitch = 0.5; // cm/channel

    // Verbosity
    int iLogLevel;
    int iFlagSpecial = 1;
    int iFlagDetails = 2;

    // Products
    vector<vector<string>> vvsProducts;
    art::InputTag   tag_mcp,
                    tag_sed,
                    tag_wir,
                    tag_hit,
                    tag_clu,
                    tag_trk,
                    tag_spt;

    // Variables
    float fMichelTimeRadius, // in µs
          fMichelSpaceRadius, // in cm
          fTrackLengthCut; // in cm

    vector<vector<ana::Binning>> vvbChan;
    ana::Binning bTick, bE;

    // TGraph2D* gMuon;
    // TGraph2D* gMichel;
    // TH1F* hWS;
    // TH1F* hSignal;
    // TH1F* hWire;
    // TH2F* hWT;

    vector<vector<TH2F*>> hHT; // [plane][section]
    // vector<vector<TGraph*>> gMichelHits; // [plane][section]
    TGraph2D* gPhantom;

    vector<vector<TH2F*>> hTrueHits; // [plane][section]
    vector<vector<TH2F*>> hFullHits; // [plane][section]
    vector<vector<TH2F*>> hCylinderHits; // [plane][section]
    vector<raw::ChannelID_t> vchMuonEnd;

    vector<vector<TH1F*>> h1E; 
    vector<TH2F*> h2Corr;
    enum hEkeys {kTrue, kReco, kFull, kCylinder, kN};


    // Functions
    geo::WireID GetWireID(geo::Point_t const& P, int plane);
    raw::ChannelID_t GetChannel(geo::Point_t const& P, int plane);

    bool IsInVolume(double x, double y, double z, double eps);
    bool IsInVolume(float x, float y, float z, float eps);
    bool IsInVolume(TLorentzVector const& V, float eps);
    bool IsInVolume(geo::Point_t const& P);

    bool IsUpright(recob::Track const& T);

    int GetSection(unsigned int ch);
    int GetPlane(unsigned int ch, int sec);

    bool Logging(bool cond, int tab, string msg, string succ, string fail);

    int iCorr(int n, int i, int j);
};


ana::Michecks::Michecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel",1)),
    vvsProducts(p.get<vector<vector<string>>>("Products")),

    fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("LengthCut")) // in cm
{

    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    // Retrieving product tags
    for (vector<string> prod : vvsProducts) {

        const string    process     = prod[0],
                        label       = prod[1],
                        instance    = prod[2],
                        type        = prod[3];

        const auto  tag = art::InputTag(label,instance);

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
    }

    int section_n = asGeo->Nchannels()/12;

    // Binning of the Frame Channel-Tick for which plane and section
    vvbChan = vector<vector<Binning>>(3, vector<Binning>(4));
    for (int plane=0; plane<3; plane++) {
        for (int section=0; section<4; section++) {
            vvbChan[plane][section].n    = section_n;
            vvbChan[plane][section].min  = (3*section + plane) * section_n;
            vvbChan[plane][section].max  = (3*section + plane + 1) * section_n - 1;
        }
    }
    bTick.n = detProp.ReadOutWindowSize()/4;
    bTick.min = 0.;
    bTick.max = detProp.ReadOutWindowSize();

    bE.n = 30;
    bE.min = 0.;
    bE.max = 90.;

    gPhantom = new TGraph2D();
    gPhantom->SetTitle("Michel without hits");
    gPhantom->SetMarkerStyle(20);
    gPhantom->SetMarkerSize(0.5);
    gPhantom->SetMarkerColor(kAzure+7);

    // gMuon = new TGraph2D();
    // gMuon->SetTitle("Truly Decaying Muon Track Points");    
    // gMuon->SetMarkerStyle(20);
    // gMuon->SetMarkerSize(0.5);
    // gMuon->SetMarkerColor(kOrange-3);


    // gMichel = new TGraph2D();
    // gMichel->SetTitle("True Michel Electron MCParticle Points");
    // gMichel->SetMarkerStyle(20);
    // gMichel->SetMarkerSize(0.5);
    // gMichel->SetMarkerColor(kAzure+7);

    // hWS = new TH1F("hWS","#Wire per channel",4,0,4);
    // hSignal = new TH1F("hSignal","Signal",60,0,60);


    // Initializing TH2F Channel-Tick for each plane and section
    hHT = vector<vector<TH2F*>>(3, vector<TH2F*>(4));
    hTrueHits = vector<vector<TH2F*>>(3, vector<TH2F*>(4));
    hFullHits = vector<vector<TH2F*>>(3, vector<TH2F*>(4));
    hCylinderHits = vector<vector<TH2F*>>(3, vector<TH2F*>(4));
    for (int plane=0; plane<3; plane++) {
        for (int section=0; section<4; section++) {
            hHT[plane][section] = new TH2F(
                Form("hHT_%c%d",'U'+plane,section),"",
                bTick.n, bTick.min, bTick.max,
                vvbChan[plane][section].n, vvbChan[plane][section].min, vvbChan[plane][section].max
            );

            hTrueHits[plane][section] = new TH2F(
                Form("hTrueHits_%c%d",'U'+plane,section),"",
                bTick.n, bTick.min, bTick.max,
                vvbChan[plane][section].n, vvbChan[plane][section].min, vvbChan[plane][section].max
            );

            hFullHits[plane][section] = new TH2F(
                Form("hFullHits_%c%d",'U'+plane,section),"",
                bTick.n, bTick.min, bTick.max,
                vvbChan[plane][section].n, vvbChan[plane][section].min, vvbChan[plane][section].max
            );

            hCylinderHits[plane][section] = new TH2F(
                Form("hCylinderHits_%c%d",'U'+plane,section),"",
                bTick.n, bTick.min, bTick.max,
                vvbChan[plane][section].n, vvbChan[plane][section].min, vvbChan[plane][section].max
            );
        }
    }

    // Initializing TGraph Channel-Tick for each plane and section
    // gMichelHits = vector<vector<TGraph*>>(3, vector<TGraph*>(4));
    // for (int plane=0; plane<3; plane++) {
    //     for (int section=0; section<4; section++) {
            // gMichelHits[plane][section] = new TGraph();
    //     }
    // }


    vector<const char*> axes = {
        "True Kinetic Energy (MeV)",
        "Reconstructed Energy (~MeV)",
        "Full Reconstructed Energy (~MeV)",
        "Reconstructed Energy in Cylinder (~MeV)",
    };
    h1E = vector<vector<TH1F*>>(kN, vector<TH1F*>(2));
    for (int i=0; i<kN; i++) {
        for (int anti=0; anti<2; anti++) {
            h1E[i][anti] = new TH1F(
                Form("h1E_%d%s",i,anti ? "posi" : "elec"),
                Form(";%s;#",axes[i]),
                bE.n, bE.min, bE.max
            );
        }
    }
    h2Corr = vector<TH2F*>(kN*(kN-1)/2);
    for (int i=0; i<kN-1; i++) {
        for (int j=i+1; j<kN; j++) {
            h2Corr[iCorr(kN,i,j)] = new TH2F(
                Form("h2Corr_%d%d",i,j),
                Form("%s;%s;%s","Michel Energy Correlations",axes[i],axes[j]),
                bE.n, bE.min, bE.max,
                bE.n, bE.min, bE.max
            );
        }
    }
}

void ana::Michecks::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    // auto const detProp = asDetProp->DataFor(e,clockData);

    if (iLogLevel >= iFlagDetails) cout << "evt#" << e.id().event() << "\r" << flush;

    auto const vh_trk = e.getValidHandle<vector<recob::Track>>(tag_trk);
    vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    vector<raw::ChannelID_t> vch_mu_end;

    vector<double> vTrueE, vRecoADC, vFullRecoADC;
    vector<int> vAnti;

    if (iLogLevel >= iFlagDetails) cout << "\tlooping over " << vp_trk.size() << " tracks..." << endl;
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {

        if (iLogLevel >= iFlagDetails) cout << "\ttrk#" << p_trk->ID()+1 << "\r" << flush;

        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());
        if (Logging(
            !mcp, 
            2, "trk to mcp...", "done", "failed")
        ) continue;    

        if (Logging(
            abs(mcp->PdgCode()) != 13,
            2, "is muon...", "yes", "no")
        ) continue;
        if (iLogLevel >= iFlagSpecial) cout << "muon pdg: " << mcp->PdgCode() << endl;

        if (Logging(
            mcp->EndProcess() != "Decay",
            2, "is decaying...", "yes", "no")
        ) continue;

        if (Logging(
            mcp->NumberDaughters() == 0,
            2, "looping over muon daughters...", "", "none")
        ) continue;

        int i_dau = mcp->NumberDaughters() - 1;
        for (; i_dau >= 0; i_dau--) {
            if (iLogLevel >= iFlagDetails) cout << "\t\tdau#" << i_dau+1 << "\r" << flush;
            
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));    

            if (Logging(
                !mcp_dau,
                3, "id to mcp...", "done", "failed")
            ) continue;

            if (Logging(
                abs(mcp_dau->PdgCode()) != 11 || mcp_dau->Process() != "Decay",
                3, "is michel...", "yes", "no")
            ) continue;

            if (Logging(
                !IsInVolume(mcp_dau->Position(0), fMichelSpaceRadius),
                3, "is inside...", "yes", "no")
            ) i_dau = -1;
            break;
        } // end loop over muon daughters
        if (i_dau == -1) continue;

        simb::MCParticle const * mcp_mich = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));

        // double TrueE = mcp_mich->E()*1e3;
        double TrueE = (mcp_mich->E() - mcp_mich->Mass())*1e3;

        int anti = mcp_mich->PdgCode() > 0 ? 0 : 1;
        if (iLogLevel >= iFlagSpecial) cout << "\tMichel pdg: " << mcp_mich->PdgCode() << endl;
        vAnti.push_back(anti);

        // hTrueE[anti]->Fill(TrueE);
        h1E[kTrue][anti]->Fill(TrueE);
        vTrueE.push_back(TrueE);
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tTrueE: " << TrueE << " MeV" << endl;

        vector<const recob::Hit*> hits_michel = truthUtil.GetMCParticleHits(clockData, *mcp_mich, e, tag_hit.label(), false);

        double RecoADC = 0;
        float prev_SummedADC = 0;
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tlooping over michel's " << hits_michel.size() << " hits...";
        for (recob::Hit const* hit : hits_michel) {

            int plane = hit->View();
            int section = GetSection(hit->Channel());
            // gMichelHits[plane][section]->AddPoint(hit->PeakTime(), hit->Channel());
            hTrueHits[plane][section]->Fill(hit->PeakTime(), hit->Channel());

            if (plane != 2) continue;

            if (hit->SummedADC() == prev_SummedADC) continue;

            RecoADC += hit->SummedADC(); 
            prev_SummedADC = hit->SummedADC();
        } // end loop over michel hits
        if (iLogLevel >= iFlagDetails) cout << "\033[92m" << " done" << "\033[0m" << endl;

        if (iLogLevel >= iFlagSpecial) {
            if (RecoADC < 1000) {
                cout << "michel with low RecoADC: " << hits_michel.size() << " hits @" << geo::Point_t(mcp_mich->Position(0).Vect()) << endl;
            }
        }
        if (!hits_michel.size()) {
            gPhantom->AddPoint(mcp_mich->Vx(0), mcp_mich->Vy(0), mcp_mich->Vz(0));
        }

        // hRecoADC[anti]->Fill(RecoADC*fADCtoE);
        h1E[kReco][anti]->Fill(RecoADC*fADCtoE);
        vRecoADC.push_back(RecoADC*fADCtoE);
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tRecoADC: " << RecoADC << endl;

        // hTrueRecoCorr->Fill(TrueE, RecoADC*fADCtoE);
        h2Corr[iCorr(kN,kTrue,kReco)]->Fill(TrueE, RecoADC*fADCtoE);


        double FullRecoADC = RecoADC;
        prev_SummedADC = 0;

        if (iLogLevel >= iFlagDetails) cout << "\t\t\tlooping over michel's " << mcp_mich->NumberDaughters() << " daughters..." << endl;
        for (int i_dau=0; i_dau< mcp_mich->NumberDaughters(); i_dau++) {
            if (iLogLevel >= iFlagDetails) cout << "\t\t\tdau#" << i_dau+1 << "\r" << flush;

            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp_mich->Daughter(i_dau));

            if (Logging(
                !mcp_dau,
                4, "id to mcp...", "done", "failed")
            ) continue;

            if (Logging(
                abs(mcp_dau->PdgCode()) != 11 && abs(mcp_dau->PdgCode()) != 22,
                4, "is electron...", "yes", "no")
            ) continue;
            
            vector<const recob::Hit*> hits_dau = truthUtil.GetMCParticleHits(clockData, *mcp_dau, e, tag_hit.label(), false);

            if (iLogLevel >= iFlagDetails) cout << "\t\t\t\tlooping over electron's " << hits_dau.size() << " hits...";
            for (recob::Hit const* hit : hits_dau) {
                // if (hit->View() != 2) continue;

                int plane = hit->View();
                int section = GetSection(hit->Channel());
                hFullHits[plane][section]->Fill(hit->PeakTime(), hit->Channel());

                if (plane != 2) continue;
                if (hit->SummedADC() == prev_SummedADC) continue;

                FullRecoADC += hit->SummedADC();
                prev_SummedADC = hit->SummedADC();
            } // end loop over electron hits
            if (iLogLevel >= iFlagDetails) cout << "\033[92m" << " done" << "\033[0m" << endl;
        } // end loop over michel daughters

        // hFullRecoADC[anti]->Fill(FullRecoADC*fADCtoE);
        h1E[kFull][anti]->Fill(FullRecoADC*fADCtoE);
        vFullRecoADC.push_back(FullRecoADC*fADCtoE);
        if (iLogLevel >= iFlagDetails) cout << "\t\t\tFullRecoADC: " << FullRecoADC << endl;

        // hRecoFullRecoCorr->Fill(RecoADC*fADCtoE, FullRecoADC*fADCtoE);
        // hTrueFullRecoCorr->Fill(TrueE, FullRecoADC*fADCtoE);
        h2Corr[iCorr(kN,kReco,kFull)]->Fill(RecoADC*fADCtoE, FullRecoADC*fADCtoE);
        h2Corr[iCorr(kN,kTrue,kFull)]->Fill(TrueE, FullRecoADC*fADCtoE);

        geo::Point_t trk_end;
        if (IsUpright(*p_trk)) {
            trk_end = p_trk->End();
        } else {
            trk_end = p_trk->Start();
        }




        // if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "\t\tMuon Track End Point: " << "\033[0m" << trk_end << endl
        //                     << "\033[94m" << "\t\tMuon MCP End Point: " << "\033[0m" << "(" << mcp_mich->EndPosition().X() << ", " << mcp_mich->EndPosition().Y() << ", " << mcp_mich->EndPosition().Z() << ")" << endl;

        // geo::TPCID tpc_mu_trk_end = asGeo->FindTPCAtPosition(trk_end);
        // geo::TPCID tpc_mu_mcp_end = asGeo->FindTPCAtPosition(geo::Point_t(mcp_mich->EndPosition().Vect()));


        // if (Logging(
        //     !tpc_mu_trk_end.isValid || !tpc_mu_mcp_end.isValid,
        //     2, "mu end tpc is valid...", "yes", "no")
        // ) {vch_mu_end.push_back(raw::InvalidChannelID); continue;}

        // if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "\t\tMuon Track End TPC: " << "\033[0m" << tpc_mu_trk_end.TPC << endl
        //                     << "\033[94m" << "\t\tMuon MCP End TPC: " << "\033[0m" << tpc_mu_mcp_end.TPC << endl;
        
        // geo::WireID wire_mu_trk_end;
        // try {
        //     wire_mu_trk_end = asGeo->NearestWireID(trk_end, geo::PlaneID(tpc_mu_trk_end,2));
        // } catch (geo::InvalidWireError const& e) {
        //     if (iLogLevel >= iFlagDetails) cout << "\033[91m" << "Invalid Track Wire Error: " << e.what() << "\033[0m" << endl;
        //     wire_mu_trk_end = e.suggestedWireID();
        // }
        // geo::WireID wire_mu_mcp_end;
        // try {
        //     wire_mu_mcp_end = asGeo->NearestWireID(geo::Point_t(mcp_mich->EndPosition().Vect()), geo::PlaneID(tpc_mu_mcp_end,2));
        // } catch (geo::InvalidWireError const& e) {
        //     if (iLogLevel >= iFlagDetails) cout << "\033[91m" << "Invalid MCP Wire Error: " << e.what() << "\033[0m" << endl;
        //     wire_mu_mcp_end = e.suggestedWireID();
        // }

        // if (Logging(
        //     !wire_mu_trk_end.isValid || !wire_mu_mcp_end.isValid,
        //     2, "muon end wire is valid...", "yes", "no")
        // ) {vch_mu_end.push_back(raw::InvalidChannelID); continue;}

        // if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "\t\tMuon Track End Wire: " << "\033[0m" << wire_mu_trk_end.Wire << endl
        //                     << "\033[94m" << "\t\tMuon MCP End Wire: " << "\033[0m" << wire_mu_mcp_end.Wire << endl;

        // raw::ChannelID_t ch_mu_trk_end = asGeo->PlaneWireToChannel(wire_mu_trk_end);
        // raw::ChannelID_t ch_mu_mcp_end = asGeo->PlaneWireToChannel(wire_mu_mcp_end);


        raw::ChannelID_t ch_mu_mcp_end = GetChannel(geo::Point_t(mcp_mich->EndPosition().Vect()), 2);
        raw::ChannelID_t ch_mu_trk_end = GetChannel(trk_end, 2);


        if(Logging(
            ch_mu_trk_end == raw::InvalidChannelID || ch_mu_mcp_end == raw::InvalidChannelID,
            2, "mu end channel is valid...", "yes", "no")
        ) {vch_mu_end.push_back(raw::InvalidChannelID); continue;}

        if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "\t\tMuon Track End Channel: " << "\033[0m" << ch_mu_trk_end << endl
                            << "\033[94m" << "\t\tMichel End Channel: " << "\033[0m" << ch_mu_mcp_end << endl;

        vch_mu_end.push_back(ch_mu_trk_end);

        // for (size_t i_tpt= trk.FirstPoint(); i_tpt < trk.LastPoint(); i_tpt++) {
        //     if (!trk.HasValidPoint(i_tpt)) continue;

        //     if (iLogLevel >= iFlagDetails) cout << "\t\t\tAdding muon point#" << i_tpt+1 
        //                 << "/" << trk.NPoints() 
        //                 << ":(" << trk.LocationAtPoint(i_tpt).X() << ","
        //                 << trk.LocationAtPoint(i_tpt).Y() << ","
        //                 << trk.LocationAtPoint(i_tpt).Z() << ")"
        //                 << "\r" << flush;
        //     //             << endl;
            
        //     gMuon->AddPoint(
        //         trk.LocationAtPoint(i_tpt).X(),
        //         trk.LocationAtPoint(i_tpt).Y(),
        //         trk.LocationAtPoint(i_tpt).Z()
        //     );
        // } // end loop over track points

    } // end loop over tracks

    for (raw::ChannelID_t ch : vch_mu_end) {
        vchMuonEnd.push_back(ch);
    }

    auto const & vh_hit = e.getValidHandle<vector<recob::Hit>>(tag_hit);
    vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);

    vector<double> vCylinderRecoADC(vch_mu_end.size(), 0);

    if (iLogLevel >= iFlagDetails) cout << "\tlooping over all hits...";
    float prev_summedADC = 0;
    for (art::Ptr<recob::Hit> const& p_hit : vp_hit) {

        int plane = p_hit->View();
        int section = GetSection(p_hit->Channel());

        hHT[plane][section]->Fill(p_hit->PeakTime(), p_hit->Channel());

        // if (plane != 2) continue;

        vector<art::Ptr<recob::Track>> vp_trk = fmp_hit2trk.at(p_hit.key());

        bool BelongsToLongTrack = false;
        for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
            BelongsToLongTrack = p_trk->Length() > fTrackLengthCut;
            if (BelongsToLongTrack) break;
        }
        if (BelongsToLongTrack) continue;

        for (unsigned int i=0; i<vch_mu_end.size(); i++) {

            if (vch_mu_end[i] == raw::InvalidChannelID) continue;

            if (abs(int(p_hit->Channel() - vch_mu_end[i])) > fMichelSpaceRadius / fChannelPitch) continue;

            hCylinderHits[plane][section]->Fill(p_hit->PeakTime(), p_hit->Channel());

            if (plane != 2) continue;

            if (p_hit->SummedADC() == prev_summedADC) continue;
            vCylinderRecoADC[i] += p_hit->SummedADC() * fADCtoE;
            prev_summedADC = p_hit->SummedADC();
        } // end loop over muon ends
    } // end loop over hits
    if (iLogLevel >= iFlagDetails) cout << "\033[92m" << " done" << "\033[0m" << endl;


    if ( vTrueE.size() != vRecoADC.size() || vTrueE.size() != vFullRecoADC.size() || vTrueE.size() != vCylinderRecoADC.size() ) {
        cerr << "\033[91m" << "Size mismatch between vectors: " << "\033[0m"
             << vTrueE.size() << " " << vRecoADC.size() << " " << vFullRecoADC.size() << " " << vCylinderRecoADC.size() << endl;
        return;
    }

    for (unsigned int i=0; i<vTrueE.size(); i++) {
        // hCylinderRecoADC[vAnti[i]]->Fill(vCylinderRecoADC[i]);
        // hTrueCylinderRecoCorr->Fill(vTrueE[i], vCylinderRecoADC[i]);
        // hFullRecoCylinderRecoCorr->Fill(vFullRecoADC[i], vCylinderRecoADC[i]);
        h1E[kCylinder][vAnti[i]]->Fill(vCylinderRecoADC[i]);
        h2Corr[iCorr(kN,kTrue,kCylinder)]->Fill(vTrueE[i], vCylinderRecoADC[i]);
        h2Corr[iCorr(kN,kFull,kCylinder)]->Fill(vFullRecoADC[i], vCylinderRecoADC[i]);
    } // end loop over michel ends

} // end analyze

void ana::Michecks::beginJob()
{
    // auto const clockData = asDetClocks->DataForJob();
    // auto const detProp = asDetProp->DataForJob(clockData);

    if (iLogLevel >= iFlagDetails) {
        geo::CryostatID cryoid{0};
        cout << "\033[94m" << "Michecks::beginJob: Detector dimension =========================================" << "\033[0m" << endl
             << "Number of channels: " << asGeo->Nchannels() << endl
             << "Number of ticks: " << "???" << endl
             << "Cryostat coordinates: " << asGeo->Cryostat(cryoid).Min() << " - " << asGeo->Cryostat(cryoid).Max() << endl
             << "\tvolume: " << asGeo->Cryostat(cryoid).Width() << " x " << asGeo->Cryostat(cryoid).Height() << " x " << asGeo->Cryostat(cryoid).Length() << endl;
        geo::TPCID tpc0{cryoid, 0}, tpcN{cryoid, asGeo->NTPC(cryoid)-1};
        cout << "\tactive coordinates: " << asGeo->TPC(tpc0).Min() << " - " << asGeo->TPC(tpcN).Max() << endl
             << "\tactive volume: " << (asGeo->TPC(tpc0).Max() - asGeo->TPC(tpcN).Min()).X() << " x " << (asGeo->TPC(tpc0).Max() - asGeo->TPC(tpcN).Min()).Y() << " x " << (asGeo->TPC(tpc0).Max() - asGeo->TPC(tpcN).Min()).Z() << endl
             << "TPCs:" << endl; 
        for (unsigned int i_tpc=0; i_tpc < asGeo->NTPC(); i_tpc++) {
            geo::TPCID tpcid{cryoid, i_tpc};
            cout << "\tTPC#" << i_tpc << "\tcoordinates: " << asGeo->TPC(tpcid).Min() << " - " << asGeo->TPC(tpcid).Max() << endl
                 << "\t\tvolume: " << asGeo->TPC(tpcid).Width() << " x " << asGeo->TPC(tpcid).Height() << " x " << asGeo->TPC(tpcid).Length() << endl
                 << "\t\tactive volume: " << asGeo->TPC(tpcid).ActiveWidth() << " x " << asGeo->TPC(tpcid).ActiveHeight() << " x " << asGeo->TPC(tpcid).ActiveLength() << endl
                 << "\t\tPlanes:";
            for (unsigned int i_plane=0; i_plane < asGeo->Nplanes(); i_plane++) {
                geo::PlaneID planeid{tpcid, i_plane};
                vector<string> plane_names = {"U","V","W"};
                cout << "\t" << plane_names[i_plane] << " #Wires: " << asGeo->Nwires(planeid);
            } // end loop over Planes
            cout << endl;
        } // end loop over TPCs
        cout << "\033[94m" << "End of Michecks::beginJob ======================================================" << "\033[0m" << endl;
    }
} // end beginJob


void ana::Michecks::endJob()
{

    if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "Michecks::endJob: Plotting section =============================================" << "\033[0m" << endl;

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);




    TCanvas* c3d = tfs->make<TCanvas>("c3d","Decay Muon and Michel Electron");

    c3d->cd();
    gPhantom->Draw("p");
    // gMuon->Draw("p");
    // gMichel->Draw("same p");
    c3d->Write();

    // TCanvas* c2 = tfs->make<TCanvas>("cwc","Wire content");

    // c2->Divide(2,2);
    // c2->cd(1);
    // gPad->SetLogy();
    // hWS->Draw();
    // c2->cd(2);
    // gPad->SetLogy();
    // hSignal->Draw();
    // c2->cd(3);
    // hWire->Draw();
    // c2->Write();

    // TCanvas* c3 = tfs->make<TCanvas>("cwt","Wire Time TH2F");

    // c3->cd();
    // gPad->SetLogz();
    // hWT->Draw("colz");
    // // hWT->ProjectionX("px",i,i+1)->Draw("same");
    // c3->Write();

    if (iLogLevel >= iFlagDetails) cout 
        << "sampling rate: " << detinfo::sampling_rate(clockData) << endl
        << "number of ticks?: " << detProp.NumberTimeSamples() << endl
        << "readout window size: " << detProp.ReadOutWindowSize() << endl;


    vector<TCanvas*> cPlanes(3);
    for (int plane=0; plane<3; plane++) {
        cPlanes[plane] = new TCanvas(Form("c%c",'U'+plane),Form("%c Plane",'U'+plane));

        if (iLogLevel >= iFlagDetails) cout << "plane#" << plane << "\r" << flush;
        cPlanes[plane]->Divide(2,2);

        for (int section=0; section<4; section++) {

            if (iLogLevel >= iFlagDetails) cout << "\tsec#" << section << "\r" << flush;
            cPlanes[plane]->cd(section+1);

            gPad->DrawFrame(
                bTick.min, vvbChan[plane][section].min, 
                bTick.max, vvbChan[plane][section].max, 
                Form("Plane %c Section %d;Tick;Channel", 'U'+plane, section)
            );
            gPad->SetLogz();

            hHT[plane][section]->SetMinimum(0.1);
            hHT[plane][section]->Draw("colz same");

            if (iLogLevel >= iFlagDetails) cout << "\t\tDrawn Hit/Time" << endl;

            hTrueHits[plane][section]->SetMarkerStyle(20);
            hTrueHits[plane][section]->SetMarkerSize(0.5);
            hTrueHits[plane][section]->SetMarkerColor(kAzure+7);
            hTrueHits[plane][section]->Draw("same");

            // TGraph* g = gMichelHits[plane][section];

            // for (int i=0; i<g->GetN(); i++) {

            //     float michel_tick = g->GetPointX(i);
            //     float michel_channel = g->GetPointY(i);
            //     float tick_radius = fMichelTimeRadius / (detinfo::sampling_rate(clockData) * 1e-3);
            //     float channel_radius = fMichelSpaceRadius / 0.5; // pitch in cm WHERE TO FIND IT ?

            //     auto e = new TEllipse(
            //         michel_tick,
            //         michel_channel,
            //         tick_radius,
            //         channel_radius
            //     );
            //     e->SetFillStyle(0);
            //     e->SetLineWidth(1);
            //     e->SetLineColor(kRed);
            //     e->Draw();

            //     // if (iLogLevel >= iFlagDetails) cout << "\t\tDrawn ellipse at (" 
            //     //             << michel_tick << "," << michel_channel << ") " 
            //     //             << "of radii (" << tick_radius << "," << channel_radius << ")" << endl;
            // } // end loop over michel hits
            // if (iLogLevel >= iFlagDetails) cout << "\t\tDrawn ellipses" << endl;
        } // end loop over sections
        cPlanes[plane]->Write();
    } // end loop over planes

    for (raw::ChannelID_t ch : vchMuonEnd) {
        if (ch == raw::InvalidChannelID) continue;
        int section = GetSection(ch);
        int plane = GetPlane(ch, section);

        cPlanes[plane]->cd(section+1);
        
        TLine* lminus = new TLine(0, ch - fMichelSpaceRadius, detProp.ReadOutWindowSize(), ch - fMichelSpaceRadius);
        TLine* lplus = new TLine(0, ch + fMichelSpaceRadius, detProp.ReadOutWindowSize(), ch + fMichelSpaceRadius);

        lminus->SetLineColor(kRed);
        lplus->SetLineColor(kRed);

        lminus->Draw();
        lplus->Draw();
    }

    for (int plane=0; plane<3; plane++) {
        cPlanes[plane]->Write();
    }


    TCanvas* cE = tfs->make<TCanvas>("cE","Michel Electron Energy");
    cE->Divide(2,2);
    for (int i=0; i<kN; i++) {
        cE->cd(i+1);
        h1E[i][0]->SetLineColor(kAzure+7);
        h1E[i][0]->SetLineWidth(2);
        h1E[i][1]->SetLineColor(kViolet-3);
        h1E[i][1]->SetLineWidth(2);
        int first = h1E[i][0]->GetMaximum() > h1E[i][1]->GetMaximum() ? 0 : 1;
        h1E[i][first]->Draw();
        h1E[i][1-first]->Draw("same");
    }
    cE->Write();

    TCanvas* cCorr = tfs->make<TCanvas>("cCorr","Michel Energy Correlation");
    if (kN % 2) {
        cCorr->Divide(kN,(kN-1)/2);
    } else {
        cCorr->Divide(kN-1,kN/2);
    }
    for (unsigned int i=0; i<h2Corr.size(); i++) {
        cCorr->cd(i+1);
        h2Corr[i]->Draw("colz");
    }
    cCorr->Write();
            

    if (iLogLevel >= iFlagDetails) cout << "\033[94m" << "End of Michecks::endJob ========================================================" << "\033[0m" << endl;
} // end endJob


bool ana::Michecks::IsInVolume(double x, double y, double z, double eps) {
    geo::CryostatID cryoid{0};
    geo::TPCID tpc0{cryoid, 0}, tpcN{cryoid, asGeo->NTPC(cryoid)-1};
    return (x-eps > asGeo->TPC(tpc0).MinX() and x+eps < asGeo->TPC(tpcN).MaxX() and
            y-eps > asGeo->TPC(tpc0).MinY() and y+eps < asGeo->TPC(tpcN).MaxY() and
            z-eps > asGeo->TPC(tpc0).MinZ() and z+eps < asGeo->TPC(tpcN).MaxZ());   
}
bool ana::Michecks::IsInVolume(float x, float y, float z, float eps) {
    return IsInVolume(double(x), double(y), double(z), double(eps));
}
bool ana::Michecks::IsInVolume(TLorentzVector const& V, float eps) {
    return IsInVolume(V.X(), V.Y(), V.Z(), double(eps));
}
bool ana::Michecks::IsInVolume(geo::Point_t const& P) {
    return IsInVolume(P.X(), P.Y(), P.Z(), 0.);
}
bool ana::Michecks::IsUpright(recob::Track const& T) {
    return T.Start().X() > T.End().X();
}


int ana::Michecks::GetSection(unsigned int ch) {
    return 4.*ch / asGeo->Nchannels();
}
int ana::Michecks::GetPlane(unsigned int ch, int sec) {
    return (12.*ch / asGeo->Nchannels() - 3*sec);
}
geo::WireID ana::Michecks::GetWireID(geo::Point_t const& P, int plane) {
    geo::TPCID tpc = asGeo->FindTPCAtPosition(P);
    if (!tpc.isValid) return geo::WireID();
    geo::WireID wire;
    try {
        wire = asGeo->NearestWireID(P, geo::PlaneID(tpc,2));
    } catch (geo::InvalidWireError const& e) {
        return e.suggestedWireID();
    }
    return wire;
}
raw::ChannelID_t ana::Michecks::GetChannel(geo::Point_t const& P, int plane) {
    geo::WireID wire = GetWireID(P, plane);
    if (!wire.isValid) return raw::InvalidChannelID;
    return asGeo->PlaneWireToChannel(wire);
}


bool ana::Michecks::Logging(bool cond, int tab, string msg, string succ, string fail) {
    if (iLogLevel >= iFlagDetails) {
        for (int i=0; i<tab; i++) cout << "\t";
        cout << msg << " ";
        if (!cond) cout << "\033[92m" << succ << "\033[0m" << endl;
        else cout << "\033[91m" << fail << "\033[0m" << endl;
    }
    return cond;
}

// Upper triangle of a symmetric matrix
int ana::Michecks::iCorr(int n, int i, int j) {
    return j - 1 + i*n - (i*(i+3))/2;
}

DEFINE_ART_MODULE(ana::Michecks)
