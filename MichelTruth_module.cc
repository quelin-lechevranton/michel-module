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
    float fFiducialCut; // cm
    bool fKeepTransportation;
    unsigned fRegN;

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
    ana::Hit StartHit, EndHit;
    int RegDirZ;
    ana::LinearRegression MuonReg;
    ana::Point EndPoint;
    float EndEnergy;

    int TrackTag;
    ana::Point TrackStartPoint, TrackEndPoint;
    float TrackLength;
    ana::Hit TrackStartHit, TrackEndHit;
    ana::Hits TrackHits;
    int TrackHitAnodeCrossing, TrackHitCathodeCrossing;
    float TrackHitCathodeTick;
    std::vector<float> TrackHitdQdx;

    int HitAnodeCrossing, HitCathodeCrossing;
    enum EnumCathodeCrossing { kNoCC, kHitOnBothSides, kAlignedHitOnBothSides };
    float HitCathodeTick;
    ana::Hits Hits;
    std::vector<float> HitProjection;
    std::vector<float> HitdQdx;

    bool HasMichel;
    // enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside, kHasMichelFiducial };

    ana::Point MichelStartPoint, MichelEndPoint;
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

    bool IsUpright(recob::Track const& T);
    std::string GetParticleName(int pdg);
};

ana::MichelTruth::MichelTruth(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}, MichelAnalyzer{p},
    fLog(p.get<bool>("Log", false)),
    fFiducialCut(p.get<float>("FiducialCut", 10.F)),
    fKeepTransportation(p.get<bool>("KeepTransportation", true)),
    fRegN(p.get<unsigned>("RegN", 6))
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

    if (fLog) {
        std::cout << "=== Geometry ===" << std::endl;
        std::cout << "  Detector: " << (geoDet==kPDVD ? "PDVD" : (geoDet==kPDHD ? "PDHD" : "Unknown")) << std::endl;
        std::cout << "  Low X: " << geoLowX << std::endl;
        std::cout << "  High X: " << geoHighX << std::endl;
        std::cout << "  Wire window: " << wireWindow << " ticks" << std::endl;
        std::cout << "================" << std::endl;
    }



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

    Hits.SetBranches(tMuon, "");
    tMuon->Branch("HitProjection", &HitProjection);
    tMuon->Branch("HitdQdx", &HitdQdx);
    StartHit.SetBranches(tMuon, "Start");
    EndHit.SetBranches(tMuon, "End");
    EndPoint.SetBranches(tMuon, "End");
    tMuon->Branch("EndEnergy", &EndEnergy);

    tMuon->Branch("TrackTag", &TrackTag);
    TrackStartPoint.SetBranches(tMuon, "TrackStart");
    TrackEndPoint.SetBranches(tMuon, "TrackEnd");
    tMuon->Branch("TrackLength", &TrackLength);
    TrackStartHit.SetBranches(tMuon, "TrackStart");
    TrackEndHit.SetBranches(tMuon, "TrackEnd");
    TrackHits.SetBranches(tMuon, "Track");
    tMuon->Branch("TrackHitAnodeCrossing", &TrackHitAnodeCrossing);
    tMuon->Branch("TrackHitCathodeCrossing", &TrackHitCathodeCrossing);
    tMuon->Branch("TrackHitCathodeTick", &TrackHitCathodeTick);
    tMuon->Branch("TrackHitdQdx", &TrackHitdQdx);

    tMuon->Branch("HitAnodeCrossing", &HitAnodeCrossing);
    tMuon->Branch("HitCathodeCrossing", &HitCathodeCrossing);
    tMuon->Branch("HitCathodeTick", &HitCathodeTick);

    tMuon->Branch("RegDirZ", &RegDirZ);
    MuonReg.SetBranches(tMuon, "");

    tMuon->Branch("HasMichel", &HasMichel);
    MichelStartPoint.SetBranches(tMuon, "MichelStart");
    MichelEndPoint.SetBranches(tMuon, "MichelEnd");
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
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e, clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const& vh_mcp = e.getValidHandle<std::vector<simb::MCParticle>>(tag_mcp);

    auto const& vh_hit = e.getValidHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Track handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrTrk vpt_ev;
    art::fill_ptr_vector(vpt_ev, vh_trk);

    auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    if (!vh_shw.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Shower handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrShw vps_ev;
    art::fill_ptr_vector(vps_ev, vh_shw);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);

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
        if (abs(mcp.PdgCode()) != 13) continue;

        if (!fKeepTransportation && mcp.EndProcess() == "Transportation") continue;

        VecPtrHit vph_mcp = ana::mcp2hits(&mcp, vph_ev, clockData, false);
        ASSERT(vph_mcp.size())

        RegDirZ = (mcp.EndZ() > mcp.Vz() ? 1 : -1);
        ana::SortedHits sh_mu = GetSortedHits(vph_mcp, RegDirZ);
        ASSERT(sh_mu)

        resetMuon();

        EventiMuon.push_back(iMuon);
        IsAnti = mcp.PdgCode() < 0;
        EndProcess = mcp.EndProcess();
        StartHit = GetHit(sh_mu.start);
        EndHit = GetHit(sh_mu.end);
        MuonReg = sh_mu.regs[ana::sec2side[geoDet][sh_mu.secs.back()]];
        EndPoint = ana::Point(mcp.EndPosition().Vect());
        EndEnergy = (mcp.EndE() - mcp.Mass()) * 1e3; // MeV

        PtrTrk pt_mu = ana::mcp2trk(&mcp, vpt_ev, clockData, fmp_trk2hit);
        TrackTag = 0;
        if (pt_mu) {
            TrackTag++;
            if (IsUpright(*pt_mu)) {
                TrackStartPoint = ana::Point(pt_mu->Start());
                TrackEndPoint = ana::Point(pt_mu->End());
            } else {
                TrackStartPoint = ana::Point(pt_mu->End());
                TrackEndPoint = ana::Point(pt_mu->Start());
            }
            TrackLength = pt_mu->Length();
            
            VecPtrHit vph_trk = fmp_trk2hit.at(pt_mu.key());
            ana::SortedHits sh_trk = GetSortedHits(vph_trk, TrackEndPoint.z > TrackStartPoint.z ? 1 : -1);
            if (sh_trk) {
                TrackTag++;
                TrackStartHit = GetHit(sh_trk.start);
                TrackEndHit = GetHit(sh_trk.end);
                for (PtrHit const& ph_trk : sh_trk.vph) {
                    if (ph_trk->View() != geo::kW) continue;
                    TrackHits.push_back(GetHit(ph_trk));
                }

                if (geoDet == kPDVD)
                    TrackHitAnodeCrossing = TrackStartHit.section < 4 
                        && geoHighX.z.isInside(TrackStartHit.space, fFiducialCut)
                        && wireWindow.isInside(TrackStartHit.tick, fFiducialCut/fTick2cm);
                else if (geoDet == kPDHD)
                    TrackHitAnodeCrossing = false;

                if (sh_trk.is_cc()) {
                    if ((sh_trk.cc.first->PeakTime()-sh_trk.cc.second->PeakTime())*fTick2cm < 3 * fCathodeGap)
                        TrackHitCathodeCrossing = kAlignedHitOnBothSides;
                    else
                        TrackHitCathodeCrossing = kHitOnBothSides;

                    TrackHitCathodeTick = sh_trk.cc.second->PeakTime();
                } else {
                    TrackHitCathodeCrossing = kNoCC;
                    TrackHitCathodeTick = -1;
                }

                VecPtrHit vph_trk_sec;
                for (PtrHit const& ph : sh_trk.vph) {
                    ana::Hit hit = GetHit(ph);
                    Hits.push_back(hit);
                    HitProjection.push_back(
                        sh_trk.regs[ana::tpc2side[geoDet][hit.tpc]].projection(
                            hit.space, hit.tick * fTick2cm
                        )
                    );
                    if (hit.section == sh_trk.secs.back()) {
                        vph_trk_sec.push_back(ph);
                    }
                }
                if (vph_trk_sec.size() > fRegN) {
                    TrackTag++;
                    TrackHitdQdx.assign(fRegN, 0.F);
                    for (auto iph=vph_trk_sec.begin()+fRegN; iph!=vph_trk_sec.end(); iph++) {
                        VecPtrHit::iterator jph = iph-fRegN;

                        double dQ = std::accumulate(
                            jph, iph, 0.,
                            [](double sum, PtrHit const& ph) {
                                return sum+ph->Integral();
                            }
                        ) / fRegN;

                        double dx = 0;
                        for (auto kph=jph; kph!=iph; kph++) {
                            if (kph == vph_trk_sec.begin())
                                dx += GetDistance(*kph, *(kph+1));
                            else if (kph == vph_trk_sec.end()-1)
                                dx += GetDistance(*(kph-1), *kph);
                            else
                                dx += .5 * (
                                    GetDistance(*(kph-1), *kph)
                                    + GetDistance(*kph, *(kph+1))
                                );
                        }
                        dx /= fRegN;

                        TrackHitdQdx.push_back(dQ / dx);
                    }
                }
            }
        }

        if (geoDet == kPDVD)
            HitAnodeCrossing = StartHit.section < 4 
                && geoHighX.z.isInside(StartHit.space, fFiducialCut)
                && wireWindow.isInside(StartHit.tick, fFiducialCut/fTick2cm);
        else if (geoDet == kPDHD)
            HitAnodeCrossing = false;

        if (sh_mu.is_cc()) {
            if ((sh_mu.cc.first->PeakTime()-sh_mu.cc.second->PeakTime())*fTick2cm < 3 * fCathodeGap)
                HitCathodeCrossing = kAlignedHitOnBothSides;
            else
                HitCathodeCrossing = kHitOnBothSides;

            HitCathodeTick = sh_mu.cc.second->PeakTime();
        } else {
            HitCathodeCrossing = kNoCC;
            HitCathodeTick = -1;
        }

        VecPtrHit vph_mu_sec;
        for (PtrHit const& ph_mu : sh_mu.vph) {
            ana::Hit hit = GetHit(ph_mu);
            Hits.push_back(hit);
            HitProjection.push_back(
                sh_mu.regs[ana::tpc2side[geoDet][hit.tpc]].projection(
                    hit.space, hit.tick * fTick2cm
                )
            );
            if (hit.section == sh_mu.secs.back()) {
                vph_mu_sec.push_back(ph_mu);
            }
        }
        ASSERT(vph_mu_sec.size() > fRegN)
        HitdQdx.assign(fRegN, 0.F);
        for (auto iph=vph_mu_sec.begin()+fRegN; iph!=vph_mu_sec.end(); iph++) {
            VecPtrHit::iterator jph = iph-fRegN;

            double dQ = std::accumulate(
                jph, iph, 0.,
                [](double sum, PtrHit const& ph) {
                    return sum+ph->Integral();
                }
            ) / fRegN;

            double dx = 0;
            for (auto kph=jph; kph!=iph; kph++) {
                if (kph == vph_mu_sec.begin())
                    dx += GetDistance(*kph, *(kph+1));
                else if (kph == vph_mu_sec.end()-1)
                    dx += GetDistance(*(kph-1), *kph);
                else
                    dx += .5 * (
                        GetDistance(*(kph-1), *kph)
                        + GetDistance(*kph, *(kph+1))
                    );
            }
            dx /= fRegN;

            HitdQdx.push_back(dQ / dx);
        }

        simb::MCParticle const* mcp_mi = GetMichelMCP(&mcp);
        if (!mcp_mi) {
            tMuon->Fill();
            iMuon++;
            EventNMuon++;
            continue;
        }

        HasMichel = true;
        MichelStartPoint = ana::Point(mcp_mi->Position().Vect());
        MichelEndPoint = ana::Point(mcp_mi->EndPosition().Vect());
        // if (geoHighX.isInside(mcp_mi->Position().Vect(), 20.F)
        //     || geoLowX.isInside(mcp_mi->Position().Vect(), 20.F)
        // )
        //     HasMichel = kHasMichelFiducial;
        // else if (geoHighX.isInside(mcp_mi->Position().Vect())
        //     || geoLowX.isInside(mcp_mi->Position().Vect())
        // )
        //     HasMichel = kHasMichelInside;
        // else
        //     HasMichel = kHasMichelOutside;

        MichelTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3; // MeV
        MichelHitEnergy = 0.F;
        // MichelHitTIDEEnergy = 0.F;
        // MichelHitEveTIDEEnergy = 0.F;
        // MichelHitSimIDEEnergy = 0.F;

        PtrTrk pt_mi = ana::mcp2trk(mcp_mi, vpt_ev, clockData, fmp_trk2hit);
        MichelTrackLength = pt_mi ? pt_mi->Length() : -1.F;
        PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
        MichelShowerLength = ps_mi ? ps_mi->Length() : -1.F;

        VecPtrHit vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true);

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
            for (sim::TrackIDE const& ide : bt_serv->HitToTrackIDEs(clockData, ph_mi))
                if (ide.trackID == mcp_mi->TrackId())
                    tide_energy += ide.energy; // MeV
                    // MichelHitTIDEEnergy += ide.energy; // MeV
            for (sim::TrackIDE const& ide : bt_serv->HitToEveTrackIDEs(clockData, ph_mi))
                if (ide.trackID == mcp_mi->TrackId())
                    eve_tide_energy += ide.energy; // MeV
                    // MichelHitEveTIDEEnergy += ide.energy; // MeV
            for (sim::IDE const* ide : bt_serv->HitToSimIDEs_Ps(clockData, ph_mi))
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
            std::vector<sim::TrackIDE> hit_ides = bt_serv->HitToTrackIDEs(clockData, ph_ev);
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
    TrackStartPoint = ana::Point();
    TrackEndPoint = ana::Point();
    TrackLength = -1.F;
    TrackStartHit = ana::Hit();
    TrackEndHit = ana::Hit();
    TrackHits.clear();
    TrackHitdQdx.clear();

    Hits.clear();
    HitProjection.clear();
    HitdQdx.clear();
    
    // Michel
    HasMichel = false;
    MichelStartPoint = ana::Point();
    MichelEndPoint = ana::Point();
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

bool ana::MichelTruth::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
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
