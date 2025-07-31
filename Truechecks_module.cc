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
    bool fKeepTransportation;

    // Output
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
    ana::Hits Hits;
    ana::Hit EndHit;
    ana::Point EndPoint;

    int RegDirZ;
    double RegM, RegP, RegR2;

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
    fMichelRadius(p.get<float>("MichelRadius", 20.F)),
    fKeepTransportation(p.get<bool>("KeepTransportation", true))
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
    tEvent = tfs->make<TTree>("event","");

    tEvent->Branch("eventRun", &evRun);
    tEvent->Branch("eventSubRun", &evSubRun);
    tEvent->Branch("eventEvent", &evEvent);

    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    EventHits.SetBranches(tEvent);

    tMuon = tfs->make<TTree>("muon", "Muon Tree"); 

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
    EndHit.SetBranches(tMuon, "End");
    EndPoint.SetBranches(tMuon, "End");

    tMuon->Branch("RegDirZ", &RegDirZ);
    tMuon->Branch("RegM", &RegM);
    tMuon->Branch("RegP", &RegP);
    tMuon->Branch("RegR2", &RegR2);

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

    EventNMuon = 0;
    EventiMuon.clear();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();

    EventHits.clear();
    for (HitPtr p_hit : vp_hit)
        if (p_hit->View() == geo::kW)
            EventHits.push_back(GetHit(p_hit));

    for (simb::MCParticle const& mcp : *vh_mcp) {
        if (abs(mcp.PdgCode()) != 13) continue;

        IsAnti = mcp.PdgCode() < 0;
        EndProcess = mcp.EndProcess();

        if (!fKeepTransportation && EndProcess == "Transportation") continue;

        HitPtrVec vp_mcp_hit = ana::mcp2hits(&mcp, vp_hit, clockData, false);
        if (vp_mcp_hit.empty()) continue;

        Hits.clear();
        for (HitPtr const& p_hit : vp_mcp_hit) {
            if (p_hit->View() != geo::kW) continue;
            Hits.push_back(GetHit(p_hit));
        }

        // HitPtrVec vp_mcp_sorted_hit = GetSortedHits(vp_mcp_hit, (mcp.EndZ() > mcp.Vz() ? 1 : -1));
        // if (vp_mcp_sorted_hit.empty()) continue;

        int RegDirZ = (mcp.EndZ() > mcp.Vz() ? 1 : -1);
        std::vector<ana::LinearRegression> side_reg(2);
        std::vector<HitPtrVec> side_hit(2);
        for (HitPtr const& p_hit : vp_mcp_hit) {
            if (p_hit->View() != geo::kW) continue;
            int side = ana::tpc2side[geoDet][p_hit->WireID().TPC];
            if (side == -1) continue;
            double z = GetSpace(p_hit->WireID());
            double t = p_hit->PeakTime() * fTick2cm;
            side_reg[side].add(z, t);
            side_hit[side].push_back(p_hit);
        }
        if (side_reg[0].n < ana::LinearRegression::nmin && side_reg[1].n < ana::LinearRegression::nmin) continue;
        for (ana::LinearRegression& reg : side_reg)
            if (reg.n >= ana::LinearRegression::nmin)
                reg.normalize();
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
                    return (s2 - s1) * RegDirZ > 0;
                    // return s2 > s1;
                }
            );
        }
        HitPtrVec vp_mcp_sorted_hit;
        if (side_reg[0].n < ana::LinearRegression::nmin) {
            vp_mcp_sorted_hit = side_hit[1];
            RegM = side_reg[1].m();
            RegP = side_reg[1].p();
            RegR2 = side_reg[1].r2();
        } else if (side_reg[1].n < ana::LinearRegression::nmin) {
            vp_mcp_sorted_hit = side_hit[0];
            RegM = side_reg[0].m();
            RegP = side_reg[0].p();
            RegR2 = side_reg[0].r2();
        } else {
            std::pair<unsigned, unsigned> side_pair = (side_reg[1].mz - side_reg[0].mz) * RegDirZ > 0
                ? std::make_pair(0, 1)
                : std::make_pair(1, 0);

            vp_mcp_sorted_hit.insert(
                vp_mcp_sorted_hit.end(),
                side_hit[side_pair.first].begin(), side_hit[side_pair.first].end()
            );
            vp_mcp_sorted_hit.insert(
                vp_mcp_sorted_hit.end(),
                side_hit[side_pair.second].begin(), side_hit[side_pair.second].end()
            );
            RegM = side_reg[side_pair.second].m();
            RegP = side_reg[side_pair.second].p();
            RegR2 = side_reg[side_pair.second].r2();
        }

        EndHit = GetHit(vp_mcp_sorted_hit.back());

        EndPoint = ana::Point(mcp.EndPosition().Vect());

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

        HitPtrVec vp_michel_hit = ana::mcp2hits(mcp_michel, vp_hit, clockData, true);

        // float Oz = EndHit.space;
        // float Ot = EndHit.tick * fTick2cm;

        // std::vector<float> radii = { 10, 20, 30, 40, 50 };
        // float r2_max = pow(radii.back(), 2);
        float r2_max = pow(60, 2);

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
        MichelHitEnergy = MichelHits.energy() * fADC2MeV;


        // MichelSphereTrueEnergy.resize(radii.size(), 0.F);
        // MichelSphereEnergy.resize(radii.size(), 0.F);
        
        // for (HitPtr const& p_hit : vp_michel_hit) {
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
        //         vp_mcp_hit.begin(),
        //         vp_mcp_hit.end(),
        //         [k=p_hit.key()](HitPtr const& p) { return p.key() == k; }
        //     ) != vp_mcp_hit.end()) continue;

        //     for (int i=radii.size()-1; i>=0; i--) {
        //         float r2 = pow(radii[i], 2);
        //         if (dr2 > r2) break;
        //         MichelSphereTrueEnergy[i] += p_hit->Integral() * fADC2MeV;
        //     }
        // }

        for (HitPtr const& p_hit : vp_hit) {
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
                vp_mcp_hit.begin(),
                vp_mcp_hit.end(),
                [k=p_hit.key()](HitPtr const& p) { return p.key() == k; }
            ) != vp_mcp_hit.end()) continue;

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
        ana::tpc2sec[geoDet][wid.TPC],
        float(GetSpace(wid)),
        p_hit->Channel(),
        p_hit->PeakTime(),
        p_hit->Integral()
    };
}

HitPtrVec ana::Truechecks::GetSortedHits(
    HitPtrVec const& vp_hit,
    int dir_z,
    geo::View_t view
) {
    unsigned const static nmin = ana::LinearRegression::nmin;

    std::vector<ana::LinearRegression> side_reg(2);
    std::vector<HitPtrVec> side_hit(2);
    for (HitPtr const& p_hit : vp_hit) {
        if (p_hit->View() != view) continue;
        int side = ana::tpc2side[geoDet][p_hit->WireID().TPC];
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
                // return (s2 - s1) * dir_z > 0;
                return s2 > s1;
            }
        );
    }

    HitPtrVec sorted_hit;
    // std::pair<unsigned, unsigned> side_pair = (side_reg[1].mz - side_reg[0].mz) * dir_z > 0
    //     ? std::make_pair(0, 1)
    //     : std::make_pair(1, 0);

    // sorted_hit.insert(
    //     sorted_hit.end(),
    //     side_hit[side_pair.first].begin(), side_hit[side_pair.first].end()
    // );
    // sorted_hit.insert(
    //     sorted_hit.end(),
    //     side_hit[side_pair.second].begin(), side_hit[side_pair.second].end()
    // );
    sorted_hit.insert(
        sorted_hit.end(),
        side_hit[0].begin(), side_hit[0].end()
    );
    sorted_hit.insert(
        sorted_hit.end(),
        side_hit[1].begin(), side_hit[1].end()
    );

    return sorted_hit;
}

DEFINE_ART_MODULE(ana::Truechecks)
