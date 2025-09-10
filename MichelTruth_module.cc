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
    ana::LinearRegression MuonReg;

    float MichelTrueEnergy;
    float MichelTrackLength, MichelShowerLength;
    ana::Hits MichelHits;
    std::vector<float> MichelHitDist;
    std::vector<float> MichelHitMuonAngle;
    ana::Vec2 MichelBary;
    float MichelBaryAngle, MichelBaryMuonAngle;
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

    void resetMuon();

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
    MuonReg.SetBranches(tMuon, "");

    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy);
    MichelHits.SetBranches(tMuon, "Michel");
    tMuon->Branch("MichelTrackLength", &MichelTrackLength);
    tMuon->Branch("MichelShowerLength", &MichelShowerLength);
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy);
    tMuon->Branch("MichelHitDist", &MichelHitDist);
    tMuon->Branch("MichelHitMuonAngle", &MichelHitMuonAngle);
    MichelBary.SetBranches(tMuon, "MichelBary");
    tMuon->Branch("MichelBaryAngle", &MichelBaryAngle);
    tMuon->Branch("MichelBaryMuonAngle", &MichelBaryMuonAngle);
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
    clockData = &asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e, *clockData);
    fTick2cm = detinfo::sampling_rate(*clockData) * 1e-3 * detProp.DriftVelocity();

    auto const& vh_mcp = e.getValidHandle<std::vector<simb::MCParticle>>(tag_mcp);

    auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
        return;
    }
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Track handle" "\033[0m" << std::endl;
        return;
    }
    art::fill_ptr_vector(vpt_ev, vh_trk);

    auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    if (!vh_shw.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Shower handle" "\033[0m" << std::endl;
        return;
    }
    art::fill_ptr_vector(vps_ev, vh_shw);

    // fmp_trk2hit = &art::FindManyP<recob::Hit>(vh_trk, e, tag_trk);
    // fop_hit2trk = &art::FindOneP<recob::Track>(vh_hit, e, tag_trk);
    // fmp_shw2hit = &art::FindManyP<recob::Hit>(vh_shw, e, tag_shw);

    EventNMuon = 0;
    EventiMuon.clear();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();

    EventHits.clear();
    for (PtrHit ph_ev : vph_ev)
        if (ph_ev->View() == geo::kW)
            EventHits.push_back(GetHit(ph_ev));

    for (simb::MCParticle const& mcp : *vh_mcp) {
        ASSERT(abs(mcp.PdgCode()) == 13)

        IsAnti = mcp.PdgCode() < 0;
        EndProcess = mcp.EndProcess();

        if (!fKeepTransportation && EndProcess == "Transportation") continue;

        VecPtrHit vph_mcp = mcp2hits(&mcp, false);
        ASSERT(vph_mcp.size())

        RegDirZ = (mcp.EndZ() > mcp.Vz() ? 1 : -1);
        ana::SortedHits sh_mu = GetSortedHits(vph_mcp, RegDirZ);
        ASSERT(sh_mu)

        resetMuon();

        EndHit = GetHit(sh_mu.end);
        MuonReg = sh_mu.regs[ana::sec2side[geoDet][sh_mu.secs.back()]];

        for (PtrHit const& ph_mu : sh_mu.vph) {
            ana::Hit hit = GetHit(ph_mu);
            Hits.push_back(hit);
            HitProjection.push_back(
                sh_mu.regs[ana::tpc2side[geoDet][hit.tpc]].projection(
                    hit.space, hit.tick * fTick2cm
                )
            );
        }

        EndPoint = ana::Point(mcp.EndPosition().Vect());

        simb::MCParticle const* mcp_mi = GetMichelMCP(&mcp);
        if (mcp_mi) {
            HasMichel = (
                geoHighX.isInside(mcp_mi->Position().Vect(), 20.F)
                || geoLowX.isInside(mcp_mi->Position().Vect(), 20.F)
            ) ? kHasMichelFiducial : (
                geoHighX.isInside(mcp_mi->Position().Vect())
                || geoLowX.isInside(mcp_mi->EndPosition().Vect())
                ? kHasMichelInside
                : kHasMichelOutside
            );
        }

        EventiMuon.push_back(iMuon);

        if (!mcp_mi) {
            tMuon->Fill();
            iMuon++;
            EventNMuon++;
            continue;
        }

        MichelTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3; // MeV
        MichelHitEnergy = 0.F;
        // MichelHitTIDEEnergy = 0.F;
        // MichelHitEveTIDEEnergy = 0.F;
        // MichelHitSimIDEEnergy = 0.F;

        PtrTrk pt_mi = mcp2trk(mcp_mi);
        MichelTrackLength = pt_mi ? pt_mi->Length() : -1.F;
        PtrShw ps_mi = mcp2shw(mcp_mi);
        MichelShowerLength = ps_mi ? ps_mi->Length() : -1.F;

        VecPtrHit vph_mi = mcp2hits(mcp_mi, true);

        // float Oz = EndHit.space;
        // float Ot = EndHit.tick * fTick2cm;

        // std::vector<float> radii = { 10, 20, 30, 40, 50 };
        // float r2_max = pow(radii.back(), 2);

        // VecPtrHit shared_hits;
        // float shared_e = 0.F;
        ana::Hits bary_hits;
        for (PtrHit const& ph_mi : vph_mi) {
            // collection hits
            if (ph_mi->View() != geo::kW) continue;

            // also from the mother muon
            if (std::find_if(
                vph_mcp.begin(),
                vph_mcp.end(),
                [k=ph_mi.key()](PtrHit const& p) { return p.key() == k; }
            ) != vph_mcp.end()) {
                // shared_hits.push_back(ph_mi);
                // shared_e += ph_mi->Integral() * fADC2MeV;
                // SharedEnergy += ph_mi->Integral();
                MichelHitIsShared.push_back(true);

                // shared is Me + mu energy

            } else {
                MichelHitIsShared.push_back(false);
            }
            ana::Hit hit = GetHit(ph_mi);
            MichelHits.push_back(hit);
            
            // G4 Energy Fractions
            float tide_energy = 0.F;
            float eve_tide_energy = 0.F;
            float sim_ide_energy = 0.F;
            for (sim::TrackIDE const& ide : bt_serv->HitToTrackIDEs(*clockData, ph_mi))
                if (ide.trackID == mcp_mi->TrackId())
                    tide_energy += ide.energy; // MeV
                    // MichelHitTIDEEnergy += ide.energy; // MeV
            for (sim::TrackIDE const& ide : bt_serv->HitToEveTrackIDEs(*clockData, ph_mi))
                if (ide.trackID == mcp_mi->TrackId())
                    eve_tide_energy += ide.energy; // MeV
                    // MichelHitEveTIDEEnergy += ide.energy; // MeV
            for (sim::IDE const* ide : bt_serv->HitToSimIDEs_Ps(*clockData, ph_mi))
                if (ide->trackID == mcp_mi->TrackId())
                    sim_ide_energy += ide->energy; // MeV
                    // MichelHitSimIDEEnergy += ide->energy; // MeV
            
            MichelHitTIDEEnergy.push_back(tide_energy);
            MichelHitEveTIDEEnergy.push_back(eve_tide_energy);
            MichelHitSimIDEEnergy.push_back(sim_ide_energy);

            // Distance from muon
            MichelHitDist.push_back(GetDistance(ph_mi, sh_mu.end));

            // Angle with muon
            // MichelHitVec.push_back(hit.vec(fTick2cm) - EndHit.vec(fTick2cm));
            float hit_angle = (hit.vec(fTick2cm) - EndHit.vec(fTick2cm)).angle();
            float da = hit_angle - MuonReg.theta(RegDirZ);
            da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            MichelHitMuonAngle.push_back(da);

            if (GetDistance(ph_mi, sh_mu.end) > 10) continue;
            bary_hits.push_back(hit);
        }
        // SharedEnergy *= fADC2MeV;
        MichelHitEnergy = MichelHits.energy();

        if (bary_hits.size()) {
            MichelBary = bary_hits.barycenter(fTick2cm);
            // MichelBaryVec = (MichelBary - EndHit.vec(fTick2cm));
            MichelBaryAngle = (MichelBary - EndHit.vec(fTick2cm)).angle();
            float da = MichelBaryAngle - MuonReg.theta(RegDirZ);
            da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            MichelBaryMuonAngle = da;
        }

        // MichelSphereTrueEnergy.resize(radii.size(), 0.F);
        // MichelSphereEnergy.resize(radii.size(), 0.F);
        
        for (PtrHit const& ph_ev : vph_ev) {
            if (ph_ev->View() != geo::kW) continue;
            ana::Hit hit = GetHit(ph_ev);
            if (hit.section != EndHit.section) continue;

            if (GetDistance(ph_ev, sh_mu.end) > 60) continue;

            // not from the mother muon
            if (std::find_if(
                vph_mcp.begin(),
                vph_mcp.end(),
                [k=ph_ev.key()](PtrHit const& p) { return p.key() == k; }
            ) != vph_mcp.end()) continue;

            // not from other muons?
            std::vector<sim::TrackIDE> hit_ides = bt_serv->HitToTrackIDEs(*clockData, ph_ev);
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
            //     MichelSphereEnergy[i] += ph_ev->Integral() * fADC2MeV;
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

void ana::MichelTruth::resetMuon() {
    Hits.clear();
    HitProjection.clear();
    
    // Michel
    HasMichel = kNoMichel;
    MichelTrueEnergy = -1.F;
    MichelHitEnergy = -1.F;
    MichelHits.clear();
    MichelHitDist.clear();
    MichelHitMuonAngle.clear();
    MichelBary = ana::Vec2{0.F, 0.F};
    MichelBaryAngle = 100.F;
    MichelBaryMuonAngle = 100.F;
    MichelHitIsShared.clear();
    MichelHitTIDEEnergy.clear();
    MichelHitEveTIDEEnergy.clear();
    MichelHitSimIDEEnergy.clear();
    NearbyHits.clear();
}


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
