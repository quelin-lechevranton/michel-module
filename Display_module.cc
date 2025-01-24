////////////////////////////////////////////////////////////////////////
// Class:       Display
// Plugin Type: analyzer
// File:        Display_module.cc
//
// Generated at Tue Jan 14 15:36:14 2025 by Jeremy Quelin Lechevranton
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
#include "larcore/Geometry/WireReadout.h"

#include "larcorealg/Geometry/Exceptions.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

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
#include "TH3F.h"
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
#include <map>
#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <utility>


namespace ana {
    class Display;
}


class ana::Display : public art::EDAnalyzer {
public:
    explicit Display(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Display(Display const&) = delete;
    Display(Display&&) = delete;
    Display& operator=(Display const&) = delete;
    Display& operator=(Display&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

private:

    // Utilities
    art::ServiceHandle<art::TFileService> tfs;

    const geo::GeometryCore* asGeo;
    const geo::WireReadoutGeom* asWire;
    const detinfo::DetectorPropertiesService* asDetProp;
    const detinfo::DetectorClocksService* asDetClocks;

    protoana::ProtoDUNETruthUtils truthUtil;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Verbosity
    int iLogLevel;
    enum EnumFlag { kImportant, kBasics, kInfos, kDetails };

    // Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag   tag_mcp,
                    tag_sed,
                    tag_wir,
                    tag_hit,
                    tag_clu,
                    tag_trk,
                    tag_spt;

    
    TGraph2D* g_spt;

};



ana::Display::Display(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products"))
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Display::Display: ============================================================" << "\033[0m" << std::endl;
    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    // Retrieving product tags
    for (std::vector<std::string> prod : vvsProducts) {

        const std::string   process     = prod[0],
                            label       = prod[1],
                            instance    = prod[2],
                            type        = prod[3];

        const art::InputTag tag = art::InputTag(label,instance);

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
    }


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Display::Display ======================================================" << "\033[0m" << std::endl;
}

void ana::Display::analyze(art::Event const& e)
{

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Display::analyze: Initialization evt#" << std::setw(5) << e.event() << " ====================================" << "\033[0m" << std::endl;

    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);

    // auto const & vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    // std::vector<art::Ptr<recob::Hit>> vp_hit;
    // art::fill_ptr_vector(vp_hit, vh_hit);

    // auto const & vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    // std::vector<art::Ptr<recob::Track>> vp_trk;
    // art::fill_ptr_vector(vp_trk, vh_trk);

    // art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    // art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);

    auto const & vh_spt = e.getValidHandle<std::vector<recob::SpacePoint>>(tag_spt);
    // std::vector<art::Ptr<recob::SpacePoint>> vp_spt;
    // art::fill_ptr_vector(vp_spt, vh_spt);

    for (recob::SpacePoint const& spt : *vh_spt) {
        g_spt->AddPoint(spt.position().Y(), spt.position().Z(), spt.position().X());
    }


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Display::analyze =======================================================" << "\033[0m" << std::endl;
} // end analyze

void ana::Display::beginJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Display::beginJob: ============================================================" << "\033[0m" << std::endl;

    g_spt = new TGraph2D();

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Display::beginJob ======================================================" << "\033[0m" << std::endl;
} // end beginJob


void ana::Display::endJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Display::endJob: ==============================================================" << "\033[0m" << std::endl;

    TCanvas* c_spt = tfs->make<TCanvas>("c_spt", "SpacePoints");

    TGraph2D* axes = new TGraph2D(2, new double[2]{-400,400}, new double[2]{-10,310}, new double[2]{-400,400});
    axes->SetTitle(";Y (cm);Z (cm);X (cm)");

    g_spt->SetMarkerStyle(20);
    g_spt->SetMarkerSize(0.2);

    c_spt->cd();
    axes->Draw("AP");
    g_spt->Draw("same P");

    c_spt->Write();

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Display::endJob ========================================================" << "\033[0m" << std::endl;
} // end endJob


DEFINE_ART_MODULE(ana::Display)
