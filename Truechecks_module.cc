////////////////////////////////////////////////////////////////////////
// Class:       Truechecks
// Plugin Type: analyzer
// File:        Truechecks_module.cc
//
// Generated at Tue Feb 4 15:23:48 2025 by Jeremy Quelin Lechevranton
////////////////////////////////////////////////////////////////////////

#include "utils.h"

using HitPtr = art::Ptr<recob::Hit>;
using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

namespace ana {
    class Truechecks;
}

class ana::Truechecks : public art::EDAnalyzer {
public:
    explicit Truechecks(fhicl::ParameterSet const& p);
    Truechecks(Truechecks const&) = delete;
    Truechecks(Truechecks&&) = delete;
    Truechecks& operator=(Truechecks const&) = delete;
    Truechecks& operator=(Truechecks&&) = delete;

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

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Products
    std::vector<std::vector<std::string>> vvsProducts;
    art::InputTag tag_mcp, tag_sed, tag_wir,
        tag_hit, tag_clu, tag_trk, tag_shw,
        tag_spt, tag_pfp;

    int geoDet;
    enum EnumDet { kPDVD, kPDHD };

    // Detector Properties
    float fADC2el; // e-/ADC.tick
    float fADC2MeV; // MeV/ADC.tick
    float fTick2cm; // cm/tick
    float fSamplingRate; // µs/tick
    float fDriftVelocity; // cm/µs
    float fCathodeGap; // cm

    std::map<geo::PlaneID, ana::axis> plane2axis;
    std::map<geo::PlaneID, double> plane2pitch;

    // Input
    float fMichelRadius;

    // Output
    TTree* tMuon;
    bool IsAnti;
    std::string EndProcess;
    int HasMichel;
    ana::Hits Hits;

    float MichelTrueEnergy;
    float MichelHitEnergy;
    float SharedEnergy;
    std::vector<float> MichelSphereTrueEnergy;
    std::vector<float> MichelSphereEnergy; 

    std::string GetParticleName(int pdg);
    double GetSpace(geo::WireID);
    ana::Hit GetHit(HitPtr const p_hit);
    // HitPtrPair GetTrackEndsHits( // return hits at track ends
    //     HitPtrVec const&, // all hits of the track
    //     HitPtrPair *pp_cathode_crossing = nullptr, // return hits at cathode crossing
    //     HitPtrVec *vp_tpc_crossing = nullptr, // return hits at TPC crossing (PDVD)
    //     HitPtrVec *vp_sorted_hit = nullptr, // return hits of the track per section and sorted
    //     geo::View_t = geo::kW // view to consider
    // );
    HitPtrVec GetSortedHits(
        HitPtrVec const&,
        int dir_z = 1,
        geo::View_t = geo::kW
    );
};



ana::Truechecks::Truechecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")),
    fADC2el(p.get<float>("ADC2el", 0.F)), // e-/ADC.tick
    fMichelRadius(p.get<float>("MichelRadius", 20.F))
{
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

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
    // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor

    fADC2MeV = fADC2el * 23.6 * 1e-6 / 0.7;

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

    tMuon = tfs->make<TTree>("muon", "Muon Tree"); 
    tMuon->Branch("IsAnti", &IsAnti);
    tMuon->Branch("EndProcess", &EndProcess);
    // tMuon->Branch("HasMichel", &HasMichel, "HasMichel/I");
    // tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy, "MichelTrueEnergy/F");
    // tMuon->Branch("MichelHitEnergy", &MichelHitEnergy, "MichelHitEnergy/F");
    // tMuon->Branch("SharedEnergy", &SharedEnergy, "SharedEnergy/F");
    // tMuon->Branch("MichelSphereTrueEnergy", &MichelSphereTrueEnergy, "MichelSphereTrueEnergy/F");
    // tMuon->Branch("MichelSphereEnergy", &MichelSphereEnergy, "MichelSphereEnergy/F");
    tMuon->Branch("HasMichel", &HasMichel);

    Hits.SetBranches(tMuon, "");

    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy);
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy);
    tMuon->Branch("SharedEnergy", &SharedEnergy);
    tMuon->Branch("MichelSphereTrueEnergy", &MichelSphereTrueEnergy);
    tMuon->Branch("MichelSphereEnergy", &MichelSphereEnergy);
}

void ana::Truechecks::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
    fDriftVelocity = detProp.DriftVelocity();
    fTick2cm = fDriftVelocity * fSamplingRate;

    auto const& vh_mcp = e.getValidHandle<std::vector<simb::MCParticle>>(tag_mcp);

    auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    for (simb::MCParticle const& mcp : *vh_mcp) {
        if (abs(mcp.PdgCode()) != 13) continue;

        IsAnti = mcp.PdgCode() < 0;
        EndProcess = mcp.EndProcess();

        HitPtrVec vp_mcp_hit = ana::mcp2hits(&mcp, vp_hit, clockData, false);
        if (vp_mcp_hit.empty()) continue;

        HitPtrVec vp_mcp_sorted_hit = GetSortedHits(vp_mcp_hit, (mcp.EndZ() > mcp.Vz() ? 1 : -1));

        Hits.clear();
        for (HitPtr const& p_hit : vp_mcp_sorted_hit) {
            if (p_hit->View() != geo::kW) continue;
            Hits.push_back(GetHit(p_hit));
        }

        HitPtr mcp_end = vp_mcp_sorted_hit.back();

        simb::MCParticle const* mcp_michel = nullptr;
        if (mcp.NumberDaughters() >= 3) {
            bool has_numu = false, has_nue = false;
            for (int i_dau=mcp.NumberDaughters()-3; i_dau<mcp.NumberDaughters(); i_dau++) {
                simb::MCParticle const* mcp_dau = pi_serv->TrackIdToParticle_P(mcp.Daughter(i_dau));    
                if (!mcp_dau) continue;
                switch (abs(mcp_dau->PdgCode())) {
                    case 14: has_numu = true; break;
                    case 12: has_nue = true; break;
                    case 11: mcp_michel = mcp_dau; break;
                    default: break;
                }
            }
            if (has_numu && has_nue && mcp_michel) {
                bool isin = false;
                for (unsigned t=0; t<asGeo->NTPC(); t++) {
                    isin = asGeo->TPC(geo::TPCID{0, t}).ContainsPosition(mcp_michel->Position().Vect());
                    if (isin) break;
                }
                if (isin)
                    HasMichel = 2;
                else
                    HasMichel = 1;
            } else {
                mcp_michel = nullptr;
                HasMichel = 0;
            }
        }

        MichelTrueEnergy = -1.F;
        MichelHitEnergy = -1.F;
        SharedEnergy = -1.F;
        MichelSphereTrueEnergy.clear();
        MichelSphereEnergy.clear();

        if (!mcp_michel) {
            tMuon->Fill();
            continue;
        }

        MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3; // MeV

        HitPtrVec vp_michel_hit = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);

        // HitPtrVec shared_hits;
        // float shared_e = 0.F;
        for (HitPtr const& p_hit : vp_michel_hit) {
            // collection hits
            if (p_hit->View() != geo::kW) continue;

            // also from the mother muon
            if (std::find_if(
                vp_mcp_hit.begin(),
                vp_mcp_hit.end(),
                [k=p_hit.key()](HitPtr const& p) { return p.key() == k; }
            ) != vp_mcp_hit.end()) {
                // shared_hits.push_back(p_hit);
                // shared_e += p_hit->Integral() * fADC2MeV;
                SharedEnergy += p_hit->Integral();



                // shared is Me + mu energy



            }
            MichelHitEnergy += p_hit->Integral();
        }
        SharedEnergy *= fADC2MeV;
        MichelHitEnergy *= fADC2MeV;

        float Oz = GetSpace(mcp_end->WireID());
        float Ot = mcp_end->PeakTime() * fTick2cm;

        std::vector<float> radii = { 10, 20, 30, 40, 50 };
        // std::vector<float> sphere_true_e(radii.size(), 0.F);
        // std::vector<float> sphere_e(radii.size(), 0.F);
        for (unsigned i=0; i<radii.size(); i++) {
            float r2 = pow(radii[i], 2);
            
            // float r2 = fMichelRadius * fMichelRadius;

            // HitPtrVec sphere_true_hits;
            float sphere_true_e = 0.F;
            for (HitPtr const& p_hit : vp_michel_hit) {
                // collection hits
                if (p_hit->View() != geo::kW) continue;

                // same section of the detector
                if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != ana::tpc2sec[geoDet][mcp_end->WireID().TPC]) continue;

                float z = GetSpace(p_hit->WireID());
                float t = p_hit->PeakTime() * fTick2cm;
                float dr2 = pow(z-Oz, 2) + pow(t-Ot, 2);

                // at less then r cm from muon's end
                if (dr2 > r2) continue;

                // not from the mother muon
                if (std::find_if(
                    vp_mcp_hit.begin(),
                    vp_mcp_hit.end(),
                    [k=p_hit.key()](HitPtr const& p) { return p.key() == k; }
                ) != vp_mcp_hit.end()) continue;

                // sphere_true_hits.push_back(p_hit);
                sphere_true_e += p_hit->Integral() * fADC2MeV;
                // sphere_true_e[i] += p_hit->Integral() * fADC2MeV;
            }

            // HitPtrVec sphere_hits;
            float sphere_e = 0.F;
            for (HitPtr const& p_hit : vp_hit) {
                // collection hits
                if (p_hit->View() != geo::kW) continue;

                // same section of the detector
                if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != ana::tpc2sec[geoDet][mcp_end->WireID().TPC]) continue;

                float z = GetSpace(p_hit->WireID());
                float t = p_hit->PeakTime() * fTick2cm;
                float dr2 = pow(z-Oz, 2) + pow(t-Ot, 2);

                // at less then r cm from muon's end
                if (dr2 > r2) continue;

                // not from the mother muon
                if (std::find_if(
                    vp_mcp_hit.begin(),
                    vp_mcp_hit.end(),
                    [k=p_hit.key()](HitPtr const& p) { return p.key() == k; }
                ) != vp_mcp_hit.end()) continue;

                // not from other muons?
                std::vector<sim::TrackIDE> hit_ides = bt_serv->HitToTrackIDEs(clockData, *p_hit);
                if (!hit_ides.empty()) {
                    std::vector<sim::TrackIDE>::const_iterator source_ide_it = std::max_element(
                        hit_ides.begin(),
                        hit_ides.end(),
                        [](sim::TrackIDE const& a, sim::TrackIDE const& b) { return a.energy < b.energy; }
                    );
                    if (source_ide_it != hit_ides.end()) {
                        simb::MCParticle const* source_mcp = pi_serv->TrackIdToParticle_P(source_ide_it->trackID);
                        if (source_mcp && abs(source_mcp->PdgCode()) == 13)
                            continue;
                    }
                }
                
                // sphere_hits.push_back(p_hit);
                sphere_e += p_hit->Integral() * fADC2MeV;
                // sphere_e[i] += p_hit->Integral() * fADC2MeV;
            }

            MichelSphereEnergy.push_back(sphere_e);
            MichelSphereTrueEnergy.push_back(sphere_true_e);
        }

        tMuon->Fill();
    }
} // end analyze

void ana::Truechecks::beginJob() {}  
void ana::Truechecks::endJob() {}

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

double ana::Truechecks::GetSpace(geo::WireID wid) {
    return plane2axis[wid].space(asWire->Wire(wid));
}

ana::Hit ana::Truechecks::GetHit(HitPtr const p_hit) {
    geo::WireID wid = p_hit->WireID();
    return ana::Hit{
        wid.TPC,
        float(GetSpace(wid)),
        p_hit->Channel(),
        p_hit->PeakTime(),
        p_hit->Integral()
    };
}

// HitPtrPair ana::Truechecks::GetTrackEndsHits(
//     HitPtrVec const& vp_hit,
//     HitPtrPair *pp_cathode_crossing,
//     HitPtrVec *vp_tpc_crossing,
//     HitPtrVec *vp_sorted_hit,
//     geo::View_t view
// ) {
//     // minimum number of hits to perform a linear regression
//     unsigned const nmin = ana::LinearRegression::nmin;

//     // split volume at de cathode
//     auto cathodeSide =
//         geoDet == kPDVD
//         ? [](geo::TPCID::TPCID_t tpc) -> int {
//                 return tpc >= 8 ? 1 : 0;
//             }
//             // geoDet == kPDHD
//         : [](geo::TPCID::TPCID_t tpc) -> int {
//                 return (tpc == 1 || tpc == 5)
//                     ? 0
//                     : ((tpc == 2 || tpc == 6) ? 1 : -1);
//             };


//     // linear regression on each side to have a curvilinear coordinate of each hit inside a track
//     // z = m*t + p
//     // struct LinearRegression {
//     //     unsigned n=0;
//     //     double mz=0, mt=0, mz2=0, mt2=0, mzt=0;
//     //     void add(double z, double t) {
//     //         mz+=z; mt+=t; mz2+=z*z; mt2+=t*t; mzt+=z*t; n++;
//     //     }
//     //     void normalize() {
//     //         mz/=n; mt/=n; mz2/=n; mt2/=n; mzt/=n;
//     //     }
//     //     double cov() const { return mzt - mz*mt; }
//     //     double varz() const { return mz2 - mz*mz; }
//     //     double vart() const { return mt2 - mt*mt; }
//     //     double m() const { return n<nmin ? 0 : cov()/vart(); }
//     //     double p() const { return mz - m()*mt; }
//     //     double r2() const { return n<nmin ? 0 : cov()*cov() / (varz()*vart()); }
//     //     double projection(double z, double t) const {
//     //         return (t + m()*(z-p())) / (1 + m()*m());
//     //     }
//     // };

//     std::vector<ana::LinearRegression> side_reg(2);
//     for (HitPtr const& p_hit : vp_hit) {
//         if (p_hit->View() != view) continue;
//         int side = cathodeSide(p_hit->WireID().TPC);
//         if (side == -1) continue; // skip hits on the other side of the anodes
//         double z = GetSpace(p_hit->WireID());
//         double t = p_hit->PeakTime() * fTick2cm;
//         side_reg[side].add(z, t);
//     }

//     // if not enough hits on both sides, return empty pair
//     if (side_reg[0].n < nmin && side_reg[1].n < nmin) return {};

//     // compute average from sum
//     for (ana::LinearRegression& reg : side_reg) reg.normalize();

//     // find the track ends on each side of the cathode
//     std::vector<HitPtrPair> side_ends(2);
//     std::vector<bounds<double>> side_mimmax(2);
//     for (HitPtr const& p_hit : vp_hit) {
//         if (p_hit->View() != view) continue;
//         int side = cathodeSide(p_hit->WireID().TPC);
//         if (side == -1) continue; // skip hits on the other side of the anodes
//         double s = side_reg[side].projection(
//             GetSpace(p_hit->WireID()),
//             p_hit->PeakTime() * fTick2cm
//         );
//         if (s > side_mimmax[side].max) {
//             side_mimmax[side].max = s;
//             side_ends[side].second = p_hit;
//         }
//         if (s < side_mimmax[side].min) {
//             side_mimmax[side].min = s;
//             side_ends[side].first = p_hit;
//         }
//     }

//     // if hits are all on one side, and no other info is requested
//     if (!vp_tpc_crossing && !vp_sorted_hit) {
//         if (side_reg[0].n < nmin)
//             return side_ends[1];
//         else if (side_reg[1].n < nmin)
//             return side_ends[0];
//     }
    
//     // given the ends of two pieces of track, find the closest ends
//     auto closestHits = [&](
//         HitPtrPair const& pph1,
//         HitPtrPair const& pph2,
//         double dmin,
//         HitPtrPair *otherHits = nullptr
//     ) -> HitPtrPair {

//         // all combinations of pairs
//         std::vector<HitPtrPair> pairs = {
//             { pph1.first, pph2.first },
//             { pph1.first, pph2.second },
//             { pph1.second, pph2.second },
//             { pph1.second, pph2.first }
//         };

//         // distance squared between all pairs
//         std::vector<double> d2s(4, 0);
//         for (unsigned i=0; i<4; i++) {
//             double zf = GetSpace(pairs[i].first->WireID());
//             double tf = pairs[i].first->PeakTime() * fTick2cm;
//             double zs = GetSpace(pairs[i].second->WireID());
//             double ts = pairs[i].second->PeakTime() * fTick2cm;
//             d2s[i] = pow(zf-zs,2) + pow(tf-ts,2);
//         }

//         // find all distances under dmin threshold
//         std::vector<unsigned> candidates_idx;
//         std::vector<double>::iterator it = d2s.begin();
//         while ((it = std::find_if(
//                 it,
//                 d2s.end(),
//                 [dmin](double d2) { return d2 < dmin*dmin; }
//             )) != d2s.end())
//             candidates_idx.push_back(std::distance(d2s.begin(), it++));
        
//         // no candidates found
//         if (candidates_idx.empty())
//             return {};

//         // get the closest pair
//         unsigned closest_idx = *std::min_element(
//             candidates_idx.begin(),
//             candidates_idx.end(),
//             [&d2s](unsigned i, unsigned j) { return d2s[i] < d2s[j]; });
        
//         // if outermost hits are requested, get the outermost pair
//         if (otherHits) {
//             unsigned other_idx = (closest_idx+2) % 4; // opposite pair
//             otherHits->first = pairs[other_idx].first;
//             otherHits->second = pairs[other_idx].second;
//         }
//         return pairs[closest_idx];
//     };

//     HitPtrPair trk_ends, cathode_crossing;
//     if (side_reg[0].n < nmin)
//         trk_ends = side_ends[1];
//     else if (side_reg[1].n < nmin)
//         trk_ends = side_ends[0];
//     else
//         cathode_crossing = closestHits(
//             side_ends[0],
//             side_ends[1],
//             2*fCathodeGap,
//             &trk_ends
//         );

//     // if cathode crossing info is requested
//     if (pp_cathode_crossing)
//         *pp_cathode_crossing = cathode_crossing;
    
//     // if no tpc crossing info is needed
//     if (geoDet == kPDHD || !vp_tpc_crossing) {
//         return trk_ends;
//     }

//     std::vector<HitPtrPair> per_sec_ends(ana::n_sec[geoDet]);
    
//     // if a sorted list of hits is requested
//     // if (vvp_sec_sorted_hits) {
//     //     // get a sorted list of hits for each section (ie. pair of TPCs)
//     //     vvp_sec_sorted_hits->clear();
//     //     vvp_sec_sorted_hits->resize(ana::n_sec[geoDet]);
//     //     for (HitPtr const& p_hit : vp_hit) {
//     //         if (p_hit->View() != view) continue;
//     //         int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
//     //         if (s == -1) continue;
//     //         vvp_sec_sorted_hits->at(s).push_back(p_hit);
//     //     }

//     //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//     //         int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
//     //         std::sort(
//     //             vvp_sec_sorted_hits->at(s).begin(), 
//     //             vvp_sec_sorted_hits->at(s).end(),
//     //             [&, &reg=side_reg[side]](
//     //                 HitPtr const& h1, HitPtr const& h2
//     //             ) -> bool {
//     //                 double const s1 = reg.projection(
//     //                     GetSpace(h1->WireID()),
//     //                     h1->PeakTime() * fTick2cm
//     //                 );
//     //                 double const s2 = reg.projection(
//     //                     GetSpace(h2->WireID()),
//     //                     h2->PeakTime() * fTick2cm
//     //                 );
//     //                 return s1 < s2;
//     //             }
//     //         );
//     //     }

//     //     // get the track ends for each section
//     //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//     //         if (vvp_sec_sorted_hits->at(s).size() < nmin) continue;
//     //         per_sec_ends[s].first = vvp_sec_sorted_hits->at(s).front();
//     //         per_sec_ends[s].second = vvp_sec_sorted_hits->at(s).back();
//     //     }

//     std::vector<HitPtrVec> vp_sec_hit(ana::n_sec[geoDet]);
//     for (HitPtr const& p_hit : vp_hit) {
//         if (p_hit->View() != view) continue;
//         int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
//         if (s == -1) continue;
//         vp_sec_hit[s].push_back(p_hit);
//     }
//     // THIS CAUSES A SEGFAULT FOR SOME REASON???
//     // auto side_sort = [&](int side) {
//     //     return [&](HitPtr const& h1, HitPtr const& h2) -> bool {
//     //         double const s1 = side_reg[side].projection(
//     //             GetSpace(h1->WireID()),
//     //             h1->PeakTime() * fTick2cm
//     //         );
//     //         double const s2 = side_reg[side].projection(
//     //             GetSpace(h2->WireID()),
//     //             h2->PeakTime() * fTick2cm
//     //         );
//     //         return s1 < s2;
//     //     };
//     // };

//     if (vp_sorted_hit) {
//         vp_sorted_hit->clear();

//         for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//             HitPtrVec& vp_sec_sorted = vp_sec_hit[s];
//             if (vp_sec_sorted.size() < nmin) continue;

//             int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
//             std::sort(
//                 vp_sec_sorted.begin(), 
//                 vp_sec_sorted.end(),
//                 [&, &reg=side_reg[side]](
//                     HitPtr const& h1, HitPtr const& h2
//                 ) -> bool {
//                     double const s1 = reg.projection(
//                         GetSpace(h1->WireID()),
//                         h1->PeakTime() * fTick2cm
//                     );
//                     double const s2 = reg.projection(
//                         GetSpace(h2->WireID()),
//                         h2->PeakTime() * fTick2cm
//                     );
//                     return s1 < s2;
//                 }
//             );

//             // get the track ends for each section
//             per_sec_ends[s].first = vp_sec_sorted.front();
//             per_sec_ends[s].second = vp_sec_sorted.back();

//             vp_sorted_hit->insert(
//                 vp_sorted_hit->end(),
//                 vp_sec_sorted.begin(), vp_sec_sorted.end()
//             );
//         }
//     } else { // only get the minmax ends of each section
//         for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//             if (vp_sec_hit[s].size() < nmin) continue;
//             int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
//             auto minmax = std::minmax_element(
//                 vp_sec_hit[s].begin(),
//                 vp_sec_hit[s].end(),
//                 [&, &reg=side_reg[side]](HitPtr const& h1, HitPtr const& h2) -> bool {
//                     double const s1 = reg.projection(
//                         GetSpace(h1->WireID()),
//                         h1->PeakTime() * fTick2cm
//                     );
//                     double const s2 = reg.projection(
//                         GetSpace(h2->WireID()),
//                         h2->PeakTime() * fTick2cm
//                     );
//                     return s1 < s2;
//                 }
//             );
//             per_sec_ends[s].first = *minmax.first;
//             per_sec_ends[s].second = *minmax.second;
//         }
//     }


//     // get the hits that are at the boundaries of two sections
//     HitPtrVec tpc_crossing;
//     bool prev = false;
//     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//         if (per_sec_ends[s].first.isNull()) {
//             prev = false;
//             continue;
//         }
//         if (prev) {
//             HitPtrPair const pp_tpc_crossing = closestHits(
//                 per_sec_ends[s-1], per_sec_ends[s], 2
//             );
//             if (pp_tpc_crossing.first.isNonnull()) {
//                 tpc_crossing.push_back(pp_tpc_crossing.first);
//                 tpc_crossing.push_back(pp_tpc_crossing.second);
//             }
//         }
//         if ((geoDet == kPDVD && s == 3) || (geoDet == kPDHD && s == 1))
//             prev = false; // cathode crossing
//         else
//             prev = true;
//     }

//     *vp_tpc_crossing = tpc_crossing;
//     return trk_ends;
// }

HitPtrVec ana::Truechecks::GetSortedHits(
    HitPtrVec const& vp_hit,
    int dir_z,
    geo::View_t view
) {
    unsigned const nmin = ana::LinearRegression::nmin;

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

    std::vector<ana::LinearRegression> side_reg(2);
    std::vector<HitPtrVec> side_hit(2);
    for (HitPtr const& p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        int side = cathodeSide(p_hit->WireID().TPC);
        if (side == -1) continue; // skip hits on the other side of the anodes
        double z = GetSpace(p_hit->WireID());
        double t = p_hit->PeakTime() * fTick2cm;
        side_reg[side].add(z, t);
        side_hit[side].push_back(p_hit);
    }

    // if not enough hits on both sides, return empty pair
    if (side_reg[0].n < nmin && side_reg[1].n < nmin) return {};

    // compute average from sum
    for (ana::LinearRegression& reg : side_reg) reg.normalize();

    // find the track ends on each side of the cathode
    for (int side=0; side<2; side++) {
        std::sort(
            side_hit[side].begin(),
            side_hit[side].end(),
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
                return (s2 - s1) * dir_z > 0;
            }
        );
    }

    HitPtrVec sorted_hit;
    std::pair<unsigned, unsigned> side_pair = (side_reg[1].mz - side_reg[0].mz) * dir_z > 0
        ? std::make_pair(0, 1)
        : std::make_pair(1, 0);

    sorted_hit.insert(
        sorted_hit.end(),
        side_hit[side_pair.first].begin(), side_hit[side_pair.first].end()
    );
    sorted_hit.insert(
        sorted_hit.end(),
        side_hit[side_pair.second].begin(), side_hit[side_pair.second].end()
    );
    return sorted_hit;
}

DEFINE_ART_MODULE(ana::Truechecks)
