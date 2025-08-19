////////////////////////////////////////////////////////////////////////
// Class:       MichelTruth
// Plugin Type: analyzer
// File:        MichelTruth_module.cc
//
// Generated at Tue Feb 4 15:23:48 2025 by Jeremy Quelin Lechevranton
////////////////////////////////////////////////////////////////////////

#include "utils.h"

namespace ana {
    class MichelTruth;
}

class ana::MichelTruth : public art::EDAnalyzer, public ana::MichelAnalyzer {
public:
    explicit MichelTruth(fhicl::ParameterSet const& p);
    MichelTruth(MichelTruth const&) = delete;
    MichelTruth(MichelTruth&&) = delete;
    MichelTruth& operator=(MichelTruth const&) = delete;
    MichelTruth& operator=(MichelTruth&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    ana::Bounds3D<float> geoHighX, geoLowX;
    ana::Bounds<float> wireWindow;
    float fCathodeGap; // cm

    // Input Parameters
    bool fLog;
    float fMichelRadius;
    bool fKeepTransportation;

    // Output Variables
    unsigned evRun, evSubRun, evEvent;

    TTree* tEvent;
    unsigned iEvent=0;
    unsigned EventNMuon;
    std::vector<unsigned> EventiMuon;
    ana::Hits EventHits;

    TTree* tMuon;
    unsigned iMuon=0;

    bool IsAnti;
    std::string EndProcess;
    int HasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside, kHasMichelFiducial };
    ana::Hits Hits;
    std::vector<float> HitProjection;
    ana::Hit EndHit;
    ana::Point EndPoint;

    int RegDirZ;
    double RegM, RegP, RegR2;
    ana::LinearRegression MuonReg;

    float MichelTrueEnergy;
    ana::Hits MichelHits;
    std::vector<float> MichelHitDist;
    std::vector<bool> MichelHitIsShared;
    std::vector<float> MichelHitTIDEEnergy;
    std::vector<float> MichelHitEveTIDEEnergy;
    std::vector<float> MichelHitSimIDEEnergy;
    float MichelHitEnergy;
    // float MichelHitTIDEEnergy;
    // float MichelHitEveTIDEEnergy;
    // float MichelHitSimIDEEnergy;
    // float SharedEnergy;
    // std::vector<float> MichelSphereTrueEnergy;
    // std::vector<float> MichelSphereEnergy; 
    ana::Hits NearbyHits;

    std::string GetParticleName(int pdg);
};

ana::MichelTruth::MichelTruth(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}, MichelAnalyzer{p},
    fLog(p.get<bool>("Log", false)),
    fMichelRadius(p.get<float>("MichelRadius", 20.F)),
    fKeepTransportation(p.get<bool>("KeepTransportation", true))
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
    wireWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(), 
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        default: break;
    }
    fCathodeGap = geoHighX.x.min - geoLowX.x.max;

    tEvent = asFile->make<TTree>("event","");

    tEvent->Branch("eventRun", &evRun);
    tEvent->Branch("eventSubRun", &evSubRun);
    tEvent->Branch("eventEvent", &evEvent);

    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    EventHits.SetBranches(tEvent);

    tMuon = asFile->make<TTree>("muon", "Muon Tree"); 

    tMuon->Branch("iEvent", &iEvent);
    tMuon->Branch("iMuon", &iMuon);
    tMuon->Branch("iMuonInEvent", &EventNMuon);

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
    tMuon->Branch("HitProjection", &HitProjection);
    EndHit.SetBranches(tMuon, "End");
    EndPoint.SetBranches(tMuon, "End");

    tMuon->Branch("RegDirZ", &RegDirZ);
    MuonReg.SetBranches(tMuon, "Reg");
    // tMuon->Branch("RegM", &RegM);
    // tMuon->Branch("RegP", &RegP);
    // tMuon->Branch("RegR2", &RegR2);

    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy);
    MichelHits.SetBranches(tMuon, "Michel");
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy);
    tMuon->Branch("MichelHitDist", &MichelHitDist);
    tMuon->Branch("MichelHitIsShared", &MichelHitIsShared);
    tMuon->Branch("MichelHitTIDEEnergy", &MichelHitTIDEEnergy);
    tMuon->Branch("MichelHitEveTIDEEnergy", &MichelHitEveTIDEEnergy);
    tMuon->Branch("MichelHitSimIDEEnergy", &MichelHitSimIDEEnergy);
    // tMuon->Branch("SharedEnergy", &SharedEnergy);
    // tMuon->Branch("MichelSphereTrueEnergy", &MichelSphereTrueEnergy);
    // tMuon->Branch("MichelSphereEnergy", &MichelSphereEnergy);
    NearbyHits.SetBranches(tMuon, "Nearby");
}

void ana::MichelTruth::analyze(art::Event const& e)
{
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const& vh_mcp = e.getValidHandle<std::vector<simb::MCParticle>>(tag_mcp);

    auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) return;
    std::vector<art::Ptr<recob::Hit>> vp_hit;
    art::fill_ptr_vector(vp_hit, vh_hit);

    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    EventNMuon = 0;
    EventiMuon.clear();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();

    EventHits.clear();
    for (PtrHit p_hit : vp_hit)
        if (p_hit->View() == geo::kW)
            EventHits.push_back(GetHit(p_hit));

    for (simb::MCParticle const& mcp : *vh_mcp) {
        if (abs(mcp.PdgCode()) != 13) continue;

        IsAnti = mcp.PdgCode() < 0;
        EndProcess = mcp.EndProcess();

        if (!fKeepTransportation && EndProcess == "Transportation") continue;

        VecPtrHit vph_mcp = ana::mcp2hits(&mcp, vp_hit, clockData, false);
        if (vph_mcp.empty()) continue;

        RegDirZ = (mcp.EndZ() > mcp.Vz() ? 1 : -1);
        ana::SortedHits sh_muon = GetSortedHits(vph_mcp, RegDirZ);
        ASSERT(sh_muon)

        EndHit = GetHit(sh_muon.end);
        MuonReg = sh_muon.regs[ana::sec2side[geoDet][sh_muon.secs.back()]];

        Hits.clear();
        HitProjection.clear();
        for (PtrHit const& p_hit : sh_muon.vph) {
            ana::Hit hit = GetHit(p_hit);
            Hits.push_back(hit);
            HitProjection.push_back(
                sh_muon.regs[ana::tpc2side[geoDet][hit.tpc]].projection(
                    hit.space, hit.tick * fTick2cm
                )
            );
        }

        EndPoint = ana::Point(mcp.EndPosition().Vect());

        simb::MCParticle const* mcp_michel = GetMichelMCP(&mcp);
        if (mcp_michel) {
            HasMichel = (
                geoHighX.isInside(mcp_michel->Position().Vect(), 20.F)
                || geoLowX.isInside(mcp_michel->Position().Vect(), 20.F)
            ) ? kHasMichelFiducial : (
                geoHighX.isInside(mcp_michel->Position().Vect())
                || geoLowX.isInside(mcp_michel->EndPosition().Vect())
                ? kHasMichelInside
                : kHasMichelOutside
            );
        } else {
            HasMichel = 0;
        }

        EventiMuon.push_back(iMuon);

        MichelTrueEnergy = -1.F;
        MichelHits.clear();
        MichelHitDist.clear();
        MichelHitIsShared.clear();
        MichelHitEnergy = -1.F;
        MichelHitTIDEEnergy.clear();
        MichelHitEveTIDEEnergy.clear();
        MichelHitSimIDEEnergy.clear();
        NearbyHits.clear();

        if (!mcp_michel) {
            tMuon->Fill();
            iMuon++;
            EventNMuon++;
            continue;
        }

        MichelTrueEnergy = (mcp_michel->E() - mcp_michel->Mass()) * 1e3; // MeV
        MichelHitEnergy = 0.F;
        // MichelHitTIDEEnergy = 0.F;
        // MichelHitEveTIDEEnergy = 0.F;
        // MichelHitSimIDEEnergy = 0.F;

        VecPtrHit vp_michel_hit = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);

        // float Oz = EndHit.space;
        // float Ot = EndHit.tick * fTick2cm;

        // std::vector<float> radii = { 10, 20, 30, 40, 50 };
        // float r2_max = pow(radii.back(), 2);
        float r2_max = pow(60, 2);

        // VecPtrHit shared_hits;
        // float shared_e = 0.F;
        for (PtrHit const& p_hit : vp_michel_hit) {
            // collection hits
            if (p_hit->View() != geo::kW) continue;

            // also from the mother muon
            if (std::find_if(
                vph_mcp.begin(),
                vph_mcp.end(),
                [k=p_hit.key()](PtrHit const& p) { return p.key() == k; }
            ) != vph_mcp.end()) {
                // shared_hits.push_back(p_hit);
                // shared_e += p_hit->Integral() * fADC2MeV;
                // SharedEnergy += p_hit->Integral();
                MichelHitIsShared.push_back(true);


                // shared is Me + mu energy

            } else {
                MichelHitIsShared.push_back(false);
            }
            ana::Hit hit = GetHit(p_hit);
            MichelHits.push_back(hit);

            float tide_energy = 0.F;
            float eve_tide_energy = 0.F;
            float sim_ide_energy = 0.F;
            for (sim::TrackIDE const& ide : bt_serv->HitToTrackIDEs(clockData, p_hit))
                if (ide.trackID == mcp_michel->TrackId())
                    tide_energy += ide.energy; // MeV
                    // MichelHitTIDEEnergy += ide.energy; // MeV
            for (sim::TrackIDE const& ide : bt_serv->HitToEveTrackIDEs(clockData, p_hit))
                if (ide.trackID == mcp_michel->TrackId())
                    eve_tide_energy += ide.energy; // MeV
                    // MichelHitEveTIDEEnergy += ide.energy; // MeV
            for (sim::IDE const* ide : bt_serv->HitToSimIDEs_Ps(clockData, p_hit))
                if (ide->trackID == mcp_michel->TrackId())
                    sim_ide_energy += ide->energy; // MeV
                    // MichelHitSimIDEEnergy += ide->energy; // MeV
            
            MichelHitTIDEEnergy.push_back(tide_energy);
            MichelHitEveTIDEEnergy.push_back(eve_tide_energy);
            MichelHitSimIDEEnergy.push_back(sim_ide_energy);

            if (hit.section != EndHit.section) {
                MichelHitDist.push_back(-1.F);
            } else {
                // float z = GetSpace(p_hit->WireID());
                // float t = p_hit->PeakTime() * fTick2cm;
                float dr2 = pow(hit.space-EndHit.space, 2) + pow((hit.tick-EndHit.tick)*fTick2cm, 2);
                MichelHitDist.push_back(sqrt(dr2));
            }
        }
        // SharedEnergy *= fADC2MeV;
        MichelHitEnergy = MichelHits.energy();

        // MichelSphereTrueEnergy.resize(radii.size(), 0.F);
        // MichelSphereEnergy.resize(radii.size(), 0.F);
        
        // for (PtrHit const& p_hit : vp_michel_hit) {
        //     // collection hits
        //     if (p_hit->View() != geo::kW) continue;

        //     // same section of the detector
        //     if (ana::tpc2sec[geoDet][p_hit->WireID().TPC] != EndHit.section) continue;

        //     float z = GetSpace(p_hit->WireID());
        //     float t = p_hit->PeakTime() * fTick2cm;
        //     float dr2 = pow(z-Oz, 2) + pow(t-Ot, 2);

        //     // at less then r cm from muon's end
        //     if (dr2 > r2_max) continue;

        //     // not from the mother muon
        //     if (std::find_if(
        //         vph_mcp.begin(),
        //         vph_mcp.end(),
        //         [k=p_hit.key()](PtrHit const& p) { return p.key() == k; }
        //     ) != vph_mcp.end()) continue;

        //     for (int i=radii.size()-1; i>=0; i--) {
        //         float r2 = pow(radii[i], 2);
        //         if (dr2 > r2) break;
        //         MichelSphereTrueEnergy[i] += p_hit->Integral() * fADC2MeV;
        //     }
        // }

        for (PtrHit const& p_hit : vp_hit) {
            // collection hits
            if (p_hit->View() != geo::kW) continue;

            ana::Hit hit = GetHit(p_hit);

            // same section of the detector
            if (hit.section != EndHit.section) continue;

            float dr2 = pow(hit.space-EndHit.space, 2) + pow((hit.tick-EndHit.tick) * fTick2cm, 2);

            // at less then r cm from muon's end
            if (dr2 > r2_max) continue;

            // not from the mother muon
            if (std::find_if(
                vph_mcp.begin(),
                vph_mcp.end(),
                [k=p_hit.key()](PtrHit const& p) { return p.key() == k; }
            ) != vph_mcp.end()) continue;

            // not from other muons?
            std::vector<sim::TrackIDE> hit_ides = bt_serv->HitToTrackIDEs(clockData, p_hit);
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
            
            NearbyHits.push_back(hit);

            // for (int i=radii.size()-1; i>=0; i--) {
            //     float r2 = pow(radii[i], 2);
            //     if (dr2 > r2) break;
            //     MichelSphereEnergy[i] += p_hit->Integral() * fADC2MeV;
            // }
        }

        tMuon->Fill();
        iMuon++;
        EventNMuon++;
    }
    tEvent->Fill();
    iEvent++;
} // end analyze

void ana::MichelTruth::beginJob() {}  
void ana::MichelTruth::endJob() {}

std::string ana::MichelTruth::GetParticleName(int pdg) {

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
        if (L==100 && Z && Z < periodic_table.size())
            return Form("%u%s%s", A, periodic_table[Z].c_str(), ex ? "*" : "");
    }

    return Form("%d", pdg);
}

DEFINE_ART_MODULE(ana::MichelTruth)
