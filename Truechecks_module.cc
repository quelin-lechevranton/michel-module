////////////////////////////////////////////////////////////////////////
// Class:       Truechecks
// Plugin Type: analyzer
// File:        Truechecks_module.cc
//
// Generated at Tue Feb 4 15:23:48 2025 by Jeremy Quelin Lechevranton
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
// #include "lardata/ArtDataHelper/TrackUtils.h"

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


// #include "ROOT/RVec.hxx"
// #include "RVec.hxx"

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
    class Truechecks;
    struct Binning {
        int n;
        float min, max;
    };
}

class ana::Truechecks : public art::EDAnalyzer {
public:
    explicit Truechecks(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Truechecks(Truechecks const&) = delete;
    Truechecks(Truechecks&&) = delete;
    Truechecks& operator=(Truechecks const&) = delete;
    Truechecks& operator=(Truechecks&&) = delete;

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

    // Conversion factors
    float feltoMeV = 23.6 * 1e-6 / 0.7; // 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    float fADCtoMeV = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    float fChannelPitch;
    float fDriftVelocity; // cm/µs
    float fSamplingRate; // µs/tick




    // Verbosity
    int iLogLevel;
    enum EnumFlag { kImportant, kBasics, kInfos, kDetails };


    // Diagnostic Variables
    std::map<std::string, unsigned> map_mup_endproc;
    std::map<std::string, unsigned> map_mum_endproc;
    unsigned n_mup=0, n_mum=0, n_mep=0, n_mem=0;
    unsigned n_cme_wh=0, n_cme_nh=0;
    float mean_cme_h=0;


    // Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag   tag_mcp,
                    tag_sed,
                    tag_wir,
                    tag_hit,
                    tag_clu,
                    tag_trk,
                    tag_spt;

    bool Log(bool cond, int flag, int tab, std::string msg, std::string succ, std::string fail);
    std::string GetParticleName(int pdg);
};



ana::Truechecks::Truechecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products"))
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Truechecks::Truechecks: ============================================================" << "\033[0m" << std::endl;
    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

    geo::WireGeo const wiregeo1 = asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 0});
    geo::WireGeo const wiregeo2 = asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 1});
    fChannelPitch = geo::WireGeo::WirePitch(wiregeo1, wiregeo2);

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

    // fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    // fDriftVelocity = detProp.DriftVelocity();

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Truechecks::Truechecks ======================================================" << "\033[0m" << std::endl;
}

void ana::Truechecks::analyze(art::Event const& e)
{

    if (iLogLevel >= kBasics) {
        std::cout << "\033[93m" << "Truechecks::analyze: Initialization evt#" << std::setw(5) << e.id().event() << " ====================================" << "\033[0m" << std::endl;
    }

    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();

    // auto const & vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    // std::vector<art::Ptr<recob::Hit>> vp_hit;
    // art::fill_ptr_vector(vp_hit, vh_hit);

    auto const & vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
    std::vector<art::Ptr<recob::Track>> vp_trk;
    art::fill_ptr_vector(vp_trk, vh_trk);

    // art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    // art::FindManyP<recob::Track> fmp_hit2trk(vh_hit, e, tag_trk);
    // art::FindManyP<recob::SpacePoint> fmp_hit2spt(vh_hit, e, tag_spt);

    std::map<std::string, unsigned> map_decay_count;
    std::map<std::string, unsigned> map_el_proc_count;

    std::vector<unsigned> n_dau;
    std::vector<unsigned> n_muioni;
    std::vector<std::vector<int>> dau_pdg;
    std::vector<std::vector<std::string>> dau_process;

    if (iLogLevel >= kInfos) std::cout << "looping over " << vp_trk.size() << " tracks..." << std::endl;
    for (art::Ptr<recob::Track> const& p_trk : vp_trk) {
        if (iLogLevel >= kInfos) std::cout << "trk#" << p_trk->ID() << "\r" << std::flush;

        if (p_trk->Length() < 40) continue;
        simb::MCParticle const * mcp = truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e, tag_trk.label());
        if (!mcp) continue;
        if (abs(mcp->PdgCode()) != 13) continue;

        map_decay_count[mcp->EndProcess()]++;

        if (mcp->EndProcess() == "Transportation") continue;

        if (mcp->PdgCode() > 0) {
            n_mum++;
            map_mum_endproc[mcp->EndProcess()]++;
        } else {
            n_mup++;
            map_mup_endproc[mcp->EndProcess()]++;
        }

        unsigned tp_n_dau=0;
        unsigned tp_n_muioni=0;
        std::vector<int> tp_dau_pdg;
        std::vector<std::string> tp_dau_process;
        for (int i_dau=0; i_dau<mcp->NumberDaughters(); i_dau++) {
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));
            if (!mcp_dau) continue;

            tp_n_dau++;

            if (mcp_dau->Process() == "muIoni" && mcp_dau->PdgCode() == 11) {
                tp_n_muioni++;
                continue;
            }

            tp_dau_pdg.push_back(mcp_dau->PdgCode());
            tp_dau_process.push_back(mcp_dau->Process());
        }
        n_dau.push_back(tp_n_dau);
        n_muioni.push_back(tp_n_muioni);
        dau_pdg.push_back(tp_dau_pdg);
        dau_process.push_back(tp_dau_process);
    
        if (mcp->NumberDaughters() < 3) continue;

        bool has_numu=false, has_nue=false;
        simb::MCParticle const * mcp_mich = nullptr;
        for (int i_dau=mcp->NumberDaughters()-3; i_dau<mcp->NumberDaughters(); i_dau++) {
            simb::MCParticle const * mcp_dau = pi_serv->TrackIdToParticle_P(mcp->Daughter(i_dau));
            if (!mcp_dau) continue;
            switch (abs(mcp_dau->PdgCode())) {
                case 14: has_numu = true; break;
                case 12: has_nue = true; break;
                case 11: mcp_mich = mcp_dau; break;
                default: break;
            }
        }

        if (!(has_nue and has_numu and mcp_mich)) continue;

        if (mcp->PdgCode() > 0) n_mem++;
        else n_mep++;

        map_el_proc_count[mcp_mich->Process()]++;

        std::vector<const recob::Hit*> v_hit_michel = truthUtil.GetMCParticleHits(clockData, *mcp_mich, e, tag_hit.label());
        if (mcp->EndProcess() == "muMinusCaptureAtRest") {
            if (v_hit_michel.size()) {
                n_cme_wh++;
                mean_cme_h += v_hit_michel.size();
            } else {
                n_cme_nh++;
            }
        }

        std::cout << "\tlooking at michel (" << GetParticleName(mcp_mich->PdgCode()) << ") / mcp::TrackID: " <<  mcp_mich->TrackId() << " / " << v_hit_michel.size() << " hits" << std::endl;
        float mich_ide_energy = 0;
        float mich_hit_energy = 0;

        unsigned i_hit=0;
        for (const recob::Hit* hit_michel : v_hit_michel) {
            if (hit_michel->View() != geo::kW) continue;

            mich_hit_energy += hit_michel->Integral() * fADCtoMeV;

            std::vector<sim::TrackIDE> v_tid = bt_serv->HitToTrackIDEs(clockData, *hit_michel); 
            std::cout << "\thit#" << i_hit++
                << " Integral: " << hit_michel->Integral() * fADCtoMeV
                << " MeV  / ROIADC: " << hit_michel->ROISummedADC() * fADCtoMeV
                << " MeV / "  << v_tid.size() << " trackIDEs"
                << std::endl;

            unsigned i_tid=0;
            for (const sim::TrackIDE& tid : v_tid) {
                simb::MCParticle const * mcp_tid = pi_serv->TrackIdToParticle_P(tid.trackID);
                std::cout << "\t\t\ttIDE#" << i_tid++
                    << " trackID: " << (tid.trackID == mcp_mich->TrackId() ? "\033[92m" : "\033[91m") << tid.trackID << "\033[0m"
                    << " (" << GetParticleName(mcp_tid->PdgCode()) << ")"
                    << " energy: " << tid.energy
                    << " energyFrac: " << tid.energyFrac
                    << " numElectrons: " << tid.numElectrons
                    << " (" << tid.numElectrons * feltoMeV << " MeV)"
                    << std::endl;
                if (tid.trackID == mcp_mich->TrackId()) {
                   mich_ide_energy += tid.energy; 
                }
            }
        }

        std::cout << "\tMichel IDE energy: " << mich_ide_energy << " MeV"
            << " / Michel Hit energy: " << mich_hit_energy << " MeV"
            << " / Michel True K-energy: " << (mcp_mich->E() - mcp_mich->Mass()) * 1e3 << " MeV"
            << std::endl;
    }

    // std::cout << "======= MU END PROC =======" << std::endl;

    // for (auto const& [key, val] : map_decay_count) {
    //     std::cout << key << ": " << val << std::endl;
    // }

    // std::cout << "======= MICHEL PROC =======" << std::endl;

    // for (auto const& [key, val] : map_el_proc_count) {
    //     std::cout << key << ": " << val << std::endl;
    // }

    // std::cout << "======= DAU DETAILS =======" << std::endl;

    // for (unsigned i=0; i<n_dau.size(); i++) {
    //     char ith[4];
    //     switch (i) {
    //         case 0: strcpy(ith, "st"); break;
    //         case 1: strcpy(ith, "nd"); break;
    //         case 2: strcpy(ith, "rd"); break;
    //         default: strcpy(ith, "th"); break;
    //     } 
    //     std::cout << i+1 << ith << " µ has " << n_dau[i] << " daughters" << std::endl;
    //     std::cout << "  e- muIoni (x" << n_muioni[i] << ")" << std::endl;
    //     for (unsigned j=0; j<dau_pdg[i].size(); j++) {
    //         std::cout << "  " << std::setw(6) << GetParticleName(dau_pdg[i][j]) << " " << dau_process[i][j] << std::endl;
    //     }
    // }

    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Truechecks::analyze =======================================================" << "\033[0m" << std::endl;
} // end analyze

void ana::Truechecks::beginJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Truechecks::beginJob: ============================================================" << "\033[0m" << std::endl;
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Truechecks::beginJob ======================================================" << "\033[0m" << std::endl;
} // end beginJob


void ana::Truechecks::endJob()
{
    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "Truechecks::endJob: ==============================================================" << "\033[0m" << std::endl;

    std::cout << "µ+ decay rate: " << 100.*n_mep / n_mup << "% (" << n_mup << ")" << std::endl;
    for (auto const& [key, val] : map_mup_endproc) {
        std::cout << "  " << key << ": " << val << std::endl;
    }
    std::cout << "µ- decay rate: " << 100.*n_mem / n_mum << "% (" << n_mum << ")" << std::endl;
    for (auto const& [key, val] : map_mum_endproc) {
        std::cout << "  " << key << ": " << val << std::endl;
    }
    std::cout << "µ- decaying after capture: " << n_cme_wh + n_cme_nh << "% (" << (n_cme_wh + n_cme_nh) / map_mum_endproc["muMinusCaptureAtRest"] << ")" << std::endl;
    std::cout << "  w/ hits: " << n_cme_wh << " (~" << mean_cme_h/n_cme_wh << " hits/michel)" << std::endl;
    std::cout << "  w/o hit: " << n_cme_nh << std::endl;


    if (iLogLevel >= kBasics) std::cout << "\033[93m" << "End of Truechecks::endJob ========================================================" << "\033[0m" << std::endl;
} // end endJob


bool ana::Truechecks::Log(bool cond, int flag, int tab, std::string msg, std::string succ, std::string fail) {
    if (iLogLevel >= flag) {
        std::cout << std::string(tab,'\t') << msg << " ";
        if (cond) std::cout << "\033[92m" << succ << "\033[0m" << std::endl;
        else std::cout << "\033[91m" << fail << "\033[0m" << std::endl;
    }
    return cond;
}

std::string ana::Truechecks::GetParticleName(int pdg) {

    std::vector<std::string> periodic_table = { "",
        "H",                                                                                                  "He", 
        "Li", "Be",                                                             "B",  "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg",                                                             "Al", "Si", "P",  "S",  "Cl", "Ar",
        "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe"
    };

    switch (pdg) {
        case 11: return "e-";
        case -11: return "e+";
        case 12: return "ve";
        case -12: return "-ve";
        case 13: return "µ-";
        case -13: return "µ+";
        case 14: return "vµ";
        case -14: return "-vµ";
        case 22: return "γ";
        case 2212: return "p";
        case 2112: return "n";
    }

    if (pdg > 1000000000) {
        unsigned ex = pdg % 10;
        unsigned A = (pdg / 10) % 1000;
        unsigned Z = (pdg / 10000) % 1000;
        unsigned L = (pdg / 10000000);
        if (L==100 && Z && Z < periodic_table.size()) return Form("%u%s%s", A, periodic_table[Z].c_str(), ex ? "*" : "");
    }

    return Form("%d", pdg);
}

DEFINE_ART_MODULE(ana::Truechecks)
