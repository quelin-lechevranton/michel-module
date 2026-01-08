#include "utils.h"


/* READ ME

Naming conventions,
some variables are prefixed by their type followed by an underscore:
 - `ph_`  for `art::Ptr<recob::Hit>`
 - `vph_` for `std::vector<art::Ptr<recob::Hit>>`

 - `pt_`  for `art::Ptr<recob::Track>`
 - `vpt_` for `std::vector<art::Ptr<recob::Track>>`

 - `ps_`  for `art::Ptr<recob::Shower>`
 - `vps_` for `std::vector<art::Ptr<recob::Shower>>`

 _ `tag_` for `art::InputTag`
 - `mcp_` for `simb::MCParticle*`

 - `sh_` for custom `ana::SortedHits`, a structure containing info on the sorted list of hits of a track

*/


namespace ana {
    class MichelAnalysis;
    const float kBadPoint = -500.F;
}

class ana::MichelAnalysis : 
    public art::EDAnalyzer, 
    private ana::MichelAnalyzer 
{
public:
    explicit MichelAnalysis(fhicl::ParameterSet const& p);
    MichelAnalysis(MichelAnalysis const&) = delete;
    MichelAnalysis(MichelAnalysis&&) = delete;
    MichelAnalysis& operator=(MichelAnalysis const&) = delete;
    MichelAnalysis& operator=(MichelAnalysis&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    ana::Bounds<float> wireWindow;
    ana::Bounds3D<float> geoHighX, geoLowX;
    float fCathodeGap; // cm

    // Input Parameters
    bool        fLog;
    bool        fKeepAll;
    float       fTrackLengthCut; // in cm
    float       fFiducialLength; // in cm
    float       fBarycenterRadius; // in cm
    float       fMichelRadius; // in cm
    // float       fNearbyRadius; // in cm
    // float       fBodyDistance; // in cm
    unsigned    fRegN;
    // bool        fCone;
    // bool        fBragg;
    // unsigned    fBraggN;
    // float       fBraggThreshold; // in MIP

    // Output Variables
    TTree *evTree, *muTree;

    // Event information
    unsigned                evRun, 
                            evSubRun, 
                            evEvent,
                            evIndex=0,
                            evNMuon;
    bool                    evIsData;
    std::vector<unsigned>   evMuonIndices;
    ana::Hits               evHits;

    // Track information
    unsigned                    muIndex=0;
    float                       muLength;
    ana::Point                  muStartPoint, 
                                muEndPoint;
    // bool                        muIsUpright, // Supposition: Muon is downward
    //                             muEndInYZ;

    // Hit from track information
    ana::Hits               muHits,
                            muSphereHits;
    ana::Hit                muStartHit, 
                            muEndHit;
    float                   muEndHitX, 
                            muEndHitY,
                            muStartHitX,
                            muStartHitY;
    ana::LinearRegression   muReg;
    std::vector<float>      muTopHitdQds,
                            muBotHitdQds,
                            muSphereHitMuonAngle;
    bool                    muHitError,
                            muAnodeCrossing,
                            muCathodeCrossing,
                            muCathodeMisaligned,
                            // muEndInX,
                            // muEndInY,
                            // muEndInZ,
                            // muEndInT,
                            muBaryHasLongTrack;
    float                   muSphereMaxShowerEnergy,
                            muSphereEnergy,
                            muSphereEnergyTP;    

    // Nearby Hits Barycenter
    ana::Hits   muBaryHits;
    ana::Vec2   muBary;
    float       muBaryAngle,
                muBaryMuonAngle;

    // Cone
    // ana::Hits muConeHits;
    // float muConeEnergy;
    // float muConeEnergyTP;
    // ana::Hits muKeyholeHits;
    // float muKeyholeEnergy;
    // float muKeyholeEnergyTP;


    // Bragg information
    // int BraggError;
    // enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
    // float MIPdQds;
    // float BraggdQds;
    // ana::Hit BraggEndHit;
    // ana::Hits BraggMuonHits;

    // Sphere
    // ana::Hits BraggSphereHits;
    // std::vector<float> BraggSphereHitMuonAngle;
    // float BraggSphereEnergy;
    // float BraggSphereEnergyTP;

    // ana::Hits BraggBaryHits;
    // ana::Vec2 BraggBary;
    // float BraggBaryAngle;
    // float BraggBaryMuonAngle;
    // bool BraggBaryHasLongTrack;

    // Cone
    // ana::Hits   NearbyBaryHits;
    // ana::Vec2   NearbyBary;
    // float       NearbyBaryAngle;
    // float       NearbyBaryMuonAngle;
    // ana::Hits   BraggConeHits;
    // float       BraggConeEnergy;
    // float       BraggConeEnergyTP;
    // ana::Hits   BraggKeyholeHits;
    // float       BraggKeyholeEnergy;
    // float       BraggKeyholeEnergyTP;

    // Truth information
    int truPdg;
    std::string truEndProcess;
    ana::Point truStartPoint;
    ana::Point truEndPoint;
    float truEndEnergy;
    // bool TrueDownward;
    ana::Hit truStartHit;
    ana::Hit truEndHit;
    // int truRegDirZ;
    ana::LinearRegression truReg;

    int truHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside, kHasMichelFiducial };
    float miTrueEnergy;
    float miTrackLength;
    // float miShowerLength;
    ana::Hits miHits;
    std::vector<float> miHitEnergyFrac;
    std::vector<float> miHitMuonAngle;
    float miHitEnergy;

    unsigned miBaryNHit;
    ana::Vec2 miBary;
    float miBaryAngle;
    float miBaryMuonAngle;
    float miConeEnergy;
    float miKeyholeEnergy;

    void resetEvent(void);
    void resetMuon(void);
};

ana::MichelAnalysis::MichelAnalysis(fhicl::ParameterSet const& p) : 
    EDAnalyzer{p}, 
    MichelAnalyzer{p},
    fLog(p.get<bool>("Log", true)),
    fKeepAll(p.get<bool>("KeepAll", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 30.F)), // in cm
    fFiducialLength(p.get<float>("FiducialLength", 20.F)), // in cm
    fBarycenterRadius(p.get<float>("BarycenterRadius", 10.F)), // in cm
    fMichelRadius(p.get<float>("MichelRadius", 20.F)), //in cm
    // fNearbyRadius(p.get<float>("NearbyRadius", 40.F)), //in cm
    // fBodyDistance(p.get<float>("BodyDistance", 20.F)), //in cm
    fRegN(p.get<unsigned>("RegN", 6))
    // fCone(p.get<bool>("Cone", false)),
    // fBragg(p.get<bool>("Bragg", false)),
    // fBraggN(p.get<unsigned>("BraggN", 6)),
    // fBraggThreshold(p.get<float>("BraggThreshold", 1.7)) // in MIP dE/dx
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
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

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX << std::endl
        << "  LowX Bounds: " << geoLowX << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl
        << "  Fiducial Length: " << fFiducialLength << " cm" << std::endl
        << "  Barycenter Radius: " << fBarycenterRadius << " cm" << std::endl
        << "  Michel Space Radius: " << fMichelRadius << " cm" << std::endl;

    evTree = asFile->make<TTree>("event","");

    evTree->Branch("eventRun", &evRun);
    evTree->Branch("eventSubRun", &evSubRun);
    evTree->Branch("eventEvent", &evEvent);
    evTree->Branch("isData", &evIsData);
    evTree->Branch("iEvent", &evIndex);
    evTree->Branch("NMuon", &evNMuon);
    evTree->Branch("iMuon", &evMuonIndices);
    evHits.SetBranches(evTree);

    muTree = asFile->make<TTree>("muon","");

    // Event
    muTree->Branch("eventRun", &evRun);
    muTree->Branch("eventSubRun", &evSubRun);
    muTree->Branch("eventEvent", &evEvent);
    muTree->Branch("isData", &evIsData);
    muTree->Branch("iEvent", &evIndex);
    muTree->Branch("iMuon", &muIndex);
    muTree->Branch("iMuonInEvent", &evNMuon);

    // Track
    muTree->Branch("muLength", &muLength);
    muStartPoint.SetBranches(muTree, "Start");
    muEndPoint.SetBranches(muTree, "End");
    // muTree->Branch("TrkEndInVolumeYZ", &muEndInVolumeYZ);

    // muTree->Branch("MIPdQds", &MIPdQds);

    // Hit
    muTree->Branch("TrkHitError", &muHitError);
    muTree->Branch("TrkHitCathodeCrossing", &muCathodeCrossing);
    muTree->Branch("TrkHitAnodeCrossing", &muAnodeCrossing);
    muStartHit.SetBranches(muTree, "Start");
    muEndHit.SetBranches(muTree, "End");
    muTree->Branch("StartHitX", &muStartHitX);
    muTree->Branch("StartHitY", &muStartHitY);
    muTree->Branch("EndHitX", &muEndHitX);
    muTree->Branch("EndHitY", &muEndHitY);
    muReg.SetBranches(muTree, "");
    // muTree->Branch("TrkHitEndInVolumeX", &muEndInX);
    // muTree->Branch("TrkHitEndInWindow", &muEndInT);
    muHits.SetBranches(muTree, "");
    // muTree->Branch("EndSecHitdQds", &muEndSecHitdQds);
    muTree->Branch("TopHitdQds", &muTopHitdQds);
    muTree->Branch("BotHitdQds", &muBotHitdQds);
    muSphereHits.SetBranches(muTree, "PandoraSphere");
    muTree->Branch("PandoraSphereHitMuonAngle", &muSphereHitMuonAngle);
    // muTree->Branch("PandoraSphereHasShower", &muSphereHasShower);
    // muTree->Branch("PandoraBaryHasShower", &muBaryHasShower);
    muTree->Branch("PandoraSphereEnergy", &muSphereEnergy); // ADC
    muTree->Branch("PandoraSphereEnergyTP", &muSphereEnergyTP); // ADC

    muBaryHits.SetBranches(muTree, "PandoraBary");
    muBary.SetBranches(muTree, "PandoraBary");
    muTree->Branch("PandoraBaryAngle", &muBaryAngle);
    muTree->Branch("PandoraBaryMuonAngle", &muBaryMuonAngle);
    muTree->Branch("PandoraBaryHasLongTrack", &muBaryHasLongTrack);
    muTree->Branch("SphereMaxShowerEnergy", &muSphereMaxShowerEnergy);

    // if (fCone) {
    //     // Cone
    //     PandoraConeHits.SetBranches(muTree, "PandoraCone");
    //     muTree->Branch("PandoraConeEnergy", &PandoraConeEnergy); // ADC
    //     muTree->Branch("PandoraConeEnergyTP", &PandoraConeEnergyTP); // ADC
    //     PandoraKeyholeHits.SetBranches(muTree, "PandoraCone");
    //     muTree->Branch("PandoraKeyholeEnergy", &PandoraKeyholeEnergy); // ADC
    //     muTree->Branch("PandoraKeyholeEnergyTP", &PandoraKeyholeEnergyTP); // ADC
    // }

    // if (fBragg) {
    //     // Bragg
    //     muTree->Branch("BraggError", &BraggError);
    //     muTree->Branch("BraggdQds", &BraggdQds);
    //     BraggEndHit.SetBranches(muTree, "BraggEnd");
    //     BraggMuonHits.SetBranches(muTree, "BraggMuon");

    //     // Sphere
    //     BraggSphereHits.SetBranches(muTree, "BraggSphere");
    //     muTree->Branch("BraggSphereHitMuonAngle", &BraggSphereHitMuonAngle);
    //     muTree->Branch("BraggSphereEnergy", &BraggSphereEnergy); // ADC
    //     muTree->Branch("BraggSphereEnergyTP", &BraggSphereEnergyTP); // ADC

    //     BraggBaryHits.SetBranches(muTree, "BraggBary");
    //     BraggBary.SetBranches(muTree, "BraggBary");
    //     muTree->Branch("BraggBaryAngle", &BraggBaryAngle);
    //     muTree->Branch("BraggBaryMuonAngle", &BraggBaryMuonAngle);
    //     muTree->Branch("BraggBaryHasLongTrack", &BraggBaryHasLongTrack);

        // if (fCone) {
        //     // Cone
        //     NearbyBaryHits.SetBranches(muTree, "NearbyBary");
        //     NearbyBary.SetBranches(muTree, "NearbyBary");
        //     muTree->Branch("NearbyBaryAngle", &NearbyBaryAngle);
        //     muTree->Branch("NearbyBaryMuonAngle", &NearbyBaryMuonAngle);
        //     BraggConeHits.SetBranches(muTree, "BraggCone");
        //     muTree->Branch("BraggConeEnergy", &BraggConeEnergy); // ADC
        //     muTree->Branch("BraggConeEnergyTP", &BraggConeEnergyTP); // ADC
        //     BraggKeyholeHits.SetBranches(muTree, "BraggCone");
        //     muTree->Branch("BraggKeyholeEnergy", &BraggKeyholeEnergy); // ADC
        //     muTree->Branch("BraggKeyholeEnergyTP", &BraggKeyholeEnergyTP); // ADC
        // }
    // }

    // Truth
    muTree->Branch("TruePdg", &truPdg);
    muTree->Branch("TrueEndProcess", &truEndProcess);
    truStartPoint.SetBranches(muTree, "TrueStart");
    truEndPoint.SetBranches(muTree, "TrueEnd");
    muTree->Branch("TrueEndEnergy", &truEndEnergy);
    // muTree->Branch("TrueDownward", &TrueDownward);
    truStartHit.SetBranches(muTree, "TrueStart");
    truEndHit.SetBranches(muTree, "TrueEnd");
    // muTree->Branch("TrueRegDirZ", &truRegDirZ);
    truReg.SetBranches(muTree, "True");

    muTree->Branch("TrueHasMichel", &truHasMichel);
    muTree->Branch("MichelTrueEnergy", &miTrueEnergy); // MeV
    muTree->Branch("MichelTrackLength", &miTrackLength); // cm
    // muTree->Branch("MichelShowerLength", &miShowerLength); // cm
    miHits.SetBranches(muTree, "Michel");
    muTree->Branch("MichelHitEnergyFrac", &miHitEnergyFrac);
    muTree->Branch("MichelHitMuonAngle", &miHitMuonAngle);
    muTree->Branch("MichelHitEnergy", &miHitEnergy); // ADC

    muTree->Branch("MichelBaryNHit", &miBaryNHit);
    miBary.SetBranches(muTree, "MichelBary");
    muTree->Branch("MichelBaryAngle", &miBaryAngle); // rad
    muTree->Branch("MichelBaryMuonAngle", &miBaryMuonAngle); // rad
    // muTree->Branch("MichelConeEnergy", &miConeEnergy); // ADC
    // muTree->Branch("MichelKeyholeEnergy", &miKeyholeEnergy); // ADC
}

void ana::MichelAnalysis::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
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

    // auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    // if (!vh_shw.isValid()) {
    //     std::cout << "\033[1;91m" "No valid recob::Shower handle" "\033[0m" << std::endl;
    //     return;
    // }
    // VecPtrShw vps_ev;
    // art::fill_ptr_vector(vps_ev, vh_shw);

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    // art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    // art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    resetEvent();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    evIsData = e.isRealData();

    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            evHits.push_back(GetHit(p_hit));

    // loop over tracks to find stopping muons
    for (PtrTrk const& pt_ev : vpt_ev) {
        if (fLog) std::cout << "e" << evIndex << "t" << pt_ev->ID() << "\r" << std::flush;
        resetMuon();

        VecPtrHit vph_mu = fmp_trk2hit.at(pt_ev.key());
        std::vector<recob::TrackHitMeta const*> const& vhm_mu = fmp_trk2hit.data(pt_ev.key());
        std::map<size_t, unsigned> map_hitkey_to_metaidx;
        ASSERT(vph_mu.size())

        ASSERT(vph_mu.size() == vhm_mu.size())
        std::vector<unsigned> bad_hit_indices;
        for (unsigned i=0; i<vph_mu.size(); i++) {
            if (vph_mu[i]->View() != geo::kW) continue;
            if (!pt_ev->HasValidPoint(vhm_mu[i]->Index())) {
                bad_hit_indices.push_back(i);
            } else {
                map_hitkey_to_metaidx[vph_mu[i].key()] = i;
            }
        }
        for (int i=bad_hit_indices.size()-1; i>=0; i--)
            vph_mu.erase(vph_mu.begin() + bad_hit_indices[i]);


        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        simb::MCParticle const* mcp_mi = nullptr;
        VecPtrHit vph_mcp_mu, vph_mi;
        std::vector<float> energyFracs_mi;
        if (mcp) {
            mcp_mi = GetMichelMCP(mcp);
            vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true, &energyFracs_mi);
        }

        if (fLog) std::cout << "\t" "\033[1;93m" "e" << evIndex << "m" << evNMuon << " (" << muIndex << ")" "\033[0m" << std::endl;
        evMuonIndices.push_back(muIndex);

        // ============================
        // Dump basic track information
        muLength = pt_ev->Length();
        // TrkChi2 = pt_ev->Chi2();
        // TrkChi2PerNdof = pt_ev->Chi2PerNdof();

        LOG(muLength >= fTrackLengthCut);
        if (!fKeepAll && muLength < fTrackLengthCut) continue;

        bool is_up =  IsUpright(*pt_ev);
        geo::Point_t Start = is_up ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = is_up ? pt_ev->End() : pt_ev->Start();
        muStartPoint = ana::Point(Start);
        muEndPoint = ana::Point(End);

        // TrkEndInVolumeYZ = geoHighX.isInsideYZ(End, fFiducialLength);
        // LOG(TrkEndInVolumeYZ);
        // if (!fKeepAll && !TrkEndInVolumeYZ) continue;

        // ===================================================
        // We need some more information on the track, mainly:
        // - which is the last point? (ASSUMING DOWNWARD TRACK, for cosmic muons)
        // - does the track cross the cathode? (for t0)
        // - some quality tags...?
        // For that, we make linear regressions on each side of the cathode
        // on the hits associated to the track given by Pandora

        /*
        
            Sorting hits for decreasing X,
            in PDVD that assums downward muons

            in case of hits on both sides (cathode crossing): 
            sh_mu.cc.first -> last hit in the top volume (side1)
            sh_mu.cc.second -> first hit in the bot volume (side0)
        
        */
        // int dirZ = End.Z() > Start.Z() ? 1 : -1;
        // ana::SortedHits sh_mu = GetSortedHits(vph_mu, dirZ);
        ana::SortedHits sh_mu = geoDet == kPDVD
            ? GetSortedHits_dirX(vph_mu, -1) // PDVD: decreasing X <-> downward
            : GetSortedHits_dirX(vph_mu, End.X() > Start.X() ? 1 : -1); // PDHD
        muHitError = !sh_mu;


        // Sort Hits by their Index in the TrackHitMeta data
        // Never tested
        /*
        ana::SortedHits sh_mu;
        sh_mu.vph = vph_mu; // collection hits with valid points
        std::sort(sh_mu.vph.begin(), sh_mu.vph.end(),
            [&map_hitkey_to_metaidx](PtrHit const& a, PtrHit const& b) {
                return map_hitkey_to_metaidx.at(a.key()) < map_hitkey_to_metaidx.at(b.key());
            }
        );
        int prev_sec = -1, prev_side = -1;
        PtrHit& prev_ph = sh_mu.vph.front();
        for (PtrHit const& ph : sh_mu.vph) {
            int sec = ana::tpc2sec.at(geoDet).at(ph->WireID().TPC);
            int side = ana::tpc2side.at(geoDet).at(ph->WireID().TPC);
            if (prev_sec == -1 || sec != prev_sec) {
                sh_mu.secs.push_back(sec);
                sh_mu.sc.push_back(prev_ph);
                sh_mu.sc.push_back(ph);
            }
            prev_sec = sec;
            if (side != -1) {
                double z = GetSpace(ph->WireID());
                double t = ph->PeakTime() * fTick2cm;
                sh_mu.regs[side].add(z, t);
            }
            if (prev_side != -1 && side != prev_side) {
                sh_mu.cc = std::make_pair(prev_ph, ph);
            }
            prev_side = side;
            prev_ph = ph;
            // sh_mu.bot_index?
            // sh_mu.endsec_index?
        }
        sh_mu.regs[0].compute();
        sh_mu.regs[1].compute();

        sh_mu.start = sh_mu.vph.front();
        sh_mu.end = sh_mu.vph.back();
        */


        LOG(!muHitError);
        if (!fKeepAll && muHitError) continue;
        if (!muHitError) {
            muStartHit = GetHit(sh_mu.start);
            size_t start_track_idx = vhm_mu[map_hitkey_to_metaidx.at(sh_mu.start.key())]->Index();
            muStartHitY = pt_ev->HasValidPoint(start_track_idx)
                ? pt_ev->LocationAtPoint(start_track_idx).Y()
                : kBadPoint;
            bool start_in_Y = geoLowX.y.isInside(muStartHitY, fFiducialLength) || geoHighX.y.isInside(muStartHitY, fFiducialLength);
            bool start_in_Z = geoLowX.z.isInside(muStartHit.space, fFiducialLength) || geoHighX.z.isInside(muStartHit.space, fFiducialLength);
            bool start_in_T = wireWindow.isInside(muStartHit.tick, fFiducialLength / fTick2cm);

            muEndHit = GetHit(sh_mu.end);
            size_t end_track_idx = vhm_mu[map_hitkey_to_metaidx.at(sh_mu.end.key())]->Index();
            muEndHitY = pt_ev->HasValidPoint(end_track_idx)
                ? pt_ev->LocationAtPoint(end_track_idx).Y()
                : kBadPoint;
            bool end_in_Y = geoLowX.y.isInside(muEndHitY, fFiducialLength) || geoHighX.y.isInside(muEndHitY, fFiducialLength);
            bool end_in_Z = geoLowX.z.isInside(muEndHit.space, fFiducialLength) || geoHighX.z.isInside(muEndHit.space, fFiducialLength);
            bool end_in_T = wireWindow.isInside(muEndHit.tick, fFiducialLength / fTick2cm);

            muReg = sh_mu.end_reg(geoDet);
            LOG(muReg.r2 >= 0.4);
            if (!fKeepAll && muReg.r2 < 0.4) continue;


            muCathodeCrossing = sh_mu.is_cc();
            muCathodeMisaligned = sh_mu.is_cc()
                && abs(sh_mu.cc.first->PeakTime()-sh_mu.cc.second->PeakTime())*fTick2cm < 3 * fCathodeGap;

            switch (geoDet) {
            case kPDVD: /* ASSUMS DOWNWARD MUON */
                muAnodeCrossing = muStartHit.section < 4
                    && geoHighX.z.isInside(muStartHit.space, fFiducialLength)
                    && wireWindow.isInside(muStartHit.tick, fFiducialLength/fTick2cm);
                break;
            case kPDHD:
                muAnodeCrossing =
                    geoHighX.z.isInside(muStartHit.space, fFiducialLength)
                    && wireWindow.isInside(muStartHit.tick, fFiducialLength/fTick2cm);
                break;
            }

            LOG(muCathodeCrossing || muAnodeCrossing);
            if (!fKeepAll && (!muAnodeCrossing && !muCathodeCrossing)) continue;

            LOG(muAnodeCrossing);
            bool start_in_X=false, end_in_X=false;
            if (muCathodeCrossing) {
                // side0: neg X | side1: pos X
                int start_side = ana::tpc2side.at(geoDet).at(muStartHit.tpc);
                int end_side = ana::tpc2side.at(geoDet).at(muEndHit.tpc);

                muStartHitX = start_side == 0
                    ? -(fCathodeGap/2) - (sh_mu.cc.second->PeakTime() - muStartHit.tick) * fTick2cm
                    : +(fCathodeGap/2) + (sh_mu.cc.first->PeakTime() - muStartHit.tick) * fTick2cm;
                start_in_X = start_side == 0
                    ? geoLowX.x.isInside(muStartHitX, fFiducialLength)
                    : geoHighX.x.isInside(muStartHitX, fFiducialLength);

                muEndHitX = end_side == 0
                    ? -(fCathodeGap/2) - (sh_mu.cc.second->PeakTime() - muEndHit.tick) * fTick2cm
                    : +(fCathodeGap/2) + (sh_mu.cc.first->PeakTime() - muEndHit.tick) * fTick2cm;
                end_in_X = end_side == 0
                    ? geoLowX.x.isInside(muEndHitX, fFiducialLength)
                    : geoHighX.x.isInside(muEndHitX, fFiducialLength);
            }

            bool start_in_XYZT = start_in_X && start_in_Y && start_in_Z && start_in_T;
            bool end_in_XYZT = end_in_X && end_in_Y && end_in_Z && end_in_T;
            LOG(start_in_XYZT);
            LOG(end_in_XYZT);
            if (!fKeepAll && !(start_in_XYZT || end_in_XYZT)) continue;

            for (PtrHit const& ph_mu : sh_mu.vph) {
                if (ph_mu->View() != geo::kW) continue;
                ana::Hit hit = GetHit(ph_mu);
                muHits.push_back(hit);
            }

            VecPtrHit vph_ev_endsec;
            for (PtrHit const& ph_ev : vph_ev) {
                if (ph_ev->View() != geo::kW) continue;
                int sec = ana::tpc2sec.at(geoDet).at(ph_ev->WireID().TPC);
                if (sec != sh_mu.end_sec()) continue;
                vph_ev_endsec.push_back(ph_ev);
            }

            float mu_end_angle = muReg.theta(muEndHit.space > muStartHit.space ? 1 : -1);
            // integrate charges around muon endpoint
            muSphereEnergy = 0;
            muSphereEnergyTP = 0;
            muBaryHasLongTrack = false;
            muSphereMaxShowerEnergy = 0;
            for (PtrHit const& ph_ev : vph_ev_endsec) {
                float dist = GetDistance(ph_ev, sh_mu.end);
                if (dist > fMichelRadius) continue;

                PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                if (pt_hit && pt_hit->Length() > fTrackLengthCut) {
                    if (dist < fBarycenterRadius) muBaryHasLongTrack = true;
                    continue;
                }

                // PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
                // if (ps_hit)
                //     for (double energy : ps_hit->Energy())
                //         if (energy > muSphereMaxShowerEnergy)
                //             muSphereMaxShowerEnergy = energy;

                if (dist < fBarycenterRadius) {
                    muBaryHits.push_back(GetHit(ph_ev));
                }

                ana::Hit hit = GetHit(ph_ev);
                muSphereEnergy += ph_ev->ROISummedADC();
                muSphereHits.push_back(hit);

                float da = (hit.vec(fTick2cm) - muEndHit.vec(fTick2cm)).angle() - mu_end_angle;
                da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                muSphereHitMuonAngle.push_back(da);

                if (std::find_if(
                    vph_mi.begin(), vph_mi.end(),
                    [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                ) != vph_mi.end())
                    muSphereEnergyTP = ph_ev->ROISummedADC();
            }

            // Cone
            LOG(muBaryHits.size());
            if (muBaryHits.size()) {
                muBary = muBaryHits.barycenter(fTick2cm);
                ana::Vec2 end_bary = muBary - muEndHit.vec(fTick2cm);
                muBaryAngle = end_bary.angle();
                float da = muBaryAngle - mu_end_angle;
                da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                muBaryMuonAngle = da;

                // if (fCone) {
                //     // float angle = end_bary.angle();
                //     PandoraConeEnergy = 0;
                //     PandoraConeEnergyTP = 0;
                //     for (PtrHit const& ph_ev : vph_ev_endsec) {
                //         ana::Hit hit = GetHit(ph_ev);
                //         PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                //         if (pt_hit && pt_hit->Length() > fTrackLengthCut) continue;

                //         float dist = GetDistance(ph_ev, sh_mu.end);
                //         if (dist > 30) continue;

                //         ana::Vec2 end_hit = hit.vec(fTick2cm) - muEndHit.vec(fTick2cm);
                //         float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                //         if (dist > 5
                //             && cosa < cos(30.F * TMath::DegToRad())
                //         ) continue;

                //         bool tp = std::find_if(
                //             vph_mi.begin(), vph_mi.end(),
                //             [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                //         ) != vph_mi.end();

                //         PandoraKeyholeEnergy += ph_ev->ROISummedADC();
                //         PandoraKeyholeHits.push_back(hit);

                //         if (tp) PandoraKeyholeEnergyTP += ph_ev->ROISummedADC();

                //         if (cosa < cos(30.F * TMath::DegToRad())) continue;

                //         PandoraConeEnergy += ph_ev->ROISummedADC();
                //         PandoraConeHits.push_back(hit);

                //         if (tp) PandoraConeEnergyTP += ph_ev->ROISummedADC();
                //     }
                // }
            }

            // End dQds
            // muEndSecHitdQds = GetdQds(sh_mu.endsec_it(), sh_mu.vph.end(), fRegN);
            muTopHitdQds = GetdQds(sh_mu.vph.begin(), sh_mu.bot_it(), fRegN);
            muBotHitdQds = GetdQds(sh_mu.bot_it(), sh_mu.vph.end(), fRegN);

            // LOG(!muEndSecHitdQds.empty());
            // if (!fKeepAll && muEndSecHitdQds.empty()) continue;

            // MIPdQds = 0;
            // if (!fBragg || std::distance(sh_mu.endsec_it(), sh_mu.vph.end()) < 2 * fBraggN) {
            //     BraggError = true;
            // } else {
            //     BraggError = false;

            //     VecPtrHit vph_mu_bragg(sh_mu.endsec_it(), sh_mu.vph.end() - 2*fBraggN);
            //     VecPtrHit vph_mu_tail(sh_mu.vph.end() - 2*fBraggN, sh_mu.vph.end());

            //     int n = fBraggN;
            //     while (n--) {
            //         float min_dist = std::numeric_limits<float>::max();
            //         PtrHit closest_hit;
            //         for (PtrHit const& ph_ev : vph_ev_endsec) {
            //             // not already in muon
            //             if (std::find_if(
            //                 vph_mu_bragg.begin(), vph_mu_bragg.end(),
            //                 [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
            //             ) != vph_mu_bragg.end()) continue;
            //             if (std::find_if(
            //                 vph_mu_tail.begin(), vph_mu_tail.end(),
            //                 [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
            //             ) != vph_mu_tail.end()) continue;

            //             float dist = GetDistance(ph_ev, vph_mu_tail.back());
            //             if (dist < min_dist) {
            //                 min_dist = dist;
            //                 closest_hit = ph_ev;
            //             }
            //         }
            //         if (closest_hit) vph_mu_tail.push_back(closest_hit);
            //         else break;
            //     }

                // unsigned i_dQds_max;
                // std::vector<float> bragg_dQds = GetdQds(vph_mu_tail, fRegN, &i_dQds_max);

                // BraggdQds = bragg_dQds[i_dQds_max];
                // for (unsigned i=0; i<=i_dQds_max; i++) {
                //     vph_mu_bragg.push_back(vph_mu_tail[i]);
                // }

                // VecPtrHit::iterator it_dQ_max = std::max_element(
                //     vph_mu_tail.begin(), vph_mu_tail.end(),
                //     [&](PtrHit const& a, PtrHit const& b) -> bool {
                //         return a->ROISummedADC() < b->ROISummedADC();
                //     }
                // );

                // MIPdQds = std::accumulate(
                //     sh_mu.endsec_it(), sh_mu.vph.end() - fBraggN, 0.F,
                //     [&](float sum, PtrHit const& ph) -> float {
                //         return sum + ph->ROISummedADC();
                //     }
                // ) / (sh_mu.vph.end() - fBraggN - sh_mu.endsec_it());
                // BraggdQds = (*it_dQ_max)->ROISummedADC();
                // for (VecPtrHit::iterator it=vph_mu_tail.begin(); it!=it_dQ_max+1; ++it)
                //     vph_mu_bragg.push_back(*it);



                // BraggEndHit = GetHit(vph_mu_bragg.back());
                // for (PtrHit const& ph_mu : vph_mu_bragg)
                //     BraggMuonHits.push_back(GetHit(ph_mu));

                // BraggBaryHasLongTrack = false;
                // for (PtrHit const& ph_ev : vph_ev_endsec) {
                //     double dist = GetDistance(ph_ev, vph_mu_bragg.back());
                //     if (dist > fMichelRadius) continue;

                //     PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                //     // not from other long track
                //     if (pt_hit
                //         && pt_hit.key() != pt_ev.key()
                //         && pt_hit->Length() > fMichelRadius
                //     ) {
                //         if (dist < fBarycenterRadius) BraggBaryHasLongTrack = true;
                //         continue;
                //     }

                //     if (std::find_if(
                //         vph_mu_bragg.begin(), vph_mu_bragg.end(),
                //         [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
                //     ) != vph_mu_bragg.end()) continue;

                //     if (dist < fBarycenterRadius) {
                //         BraggBaryHits.push_back(GetHit(ph_ev));
                //     }

                //     // PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
                //     // BraggSphereHasShower = bool(ps_hit);

                //     ana::Hit hit = GetHit(ph_ev);
                //     BraggSphereEnergy += ph_ev->ROISummedADC();
                //     BraggSphereHits.push_back(hit);

                //     float da = (hit.vec(fTick2cm) - TrkEndHit.vec(fTick2cm)).angle() - TrkReg.theta(dirZ);
                //     da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                //     BraggSphereHitMuonAngle.push_back(da);

                //     if (std::find_if(
                //         vph_mi.begin(), vph_mi.end(),
                //         [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                //     ) != vph_mi.end())
                //         BraggSphereEnergyTP = ph_ev->ROISummedADC();
                // }

                // LOG(BraggBaryHits.size());
                // if (BraggBaryHits.size()) {
                //     BraggBary = BraggBaryHits.barycenter(fTick2cm);
                //     ana::Vec2 end_bary = BraggBary - BraggEndHit.vec(fTick2cm);
                //     BraggBaryAngle = end_bary.angle();
                //     float da = BraggBaryAngle - TrkReg.theta(dirZ);
                //     da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                //     BraggBaryMuonAngle = da;
                // }

                // if (fCone) {

                // }
                // // Cone
                // for (PtrHit const& ph_near : vph_near) {
                //     if (GetDistance(ph_near, bragg.end) > 10) continue;
                //     NearbyBaryHits.push_back(GetHit(ph_near));
                // }

                // LOG(NearbyBaryHits.size());
                // if (NearbyBaryHits.size()) {
                //     NearbyBary = NearbyBaryHits.barycenter(fTick2cm);
                //     ana::Vec2 end_bary = NearbyBary - BraggEndHit.vec(fTick2cm);
                //     NearbyBaryAngle = end_bary.angle();
                //     float da = NearbyBaryAngle - TrkReg.theta(dirZ);
                //     da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                //     NearbyBaryMuonAngle = da;

                //     // float angle = end_bary.angle();
                //     BraggConeEnergy = 0;
                //     BraggConeEnergyTP = 0;
                //     for (PtrHit const& ph_near : vph_near) {
                //         float dist = GetDistance(ph_near, bragg.end);
                //         if (dist > 30) continue;

                //         ana::Hit hit = GetHit(ph_near);
                //         ana::Vec2 end_hit = hit.vec(fTick2cm) - BraggEndHit.vec(fTick2cm);
                //         float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                //         if (dist > 5
                //             && cosa < cos(30.F * TMath::DegToRad())
                //         ) continue;

                //         bool tp = std::find_if(
                //             vph_mi.begin(), vph_mi.end(),
                //             [&ph_near](PtrHit const& h) -> bool { return h.key() == ph_near.key(); }
                //         ) != vph_mi.end();

                //         BraggKeyholeEnergy += ph_near->ROISummedADC();
                //         BraggKeyholeHits.push_back(hit);

                //         if (tp) BraggKeyholeEnergyTP += ph_near->ROISummedADC();

                //         if (cosa < cos(30.F * TMath::DegToRad())) continue;

                //         BraggConeEnergy += ph_near->ROISummedADC();
                //         BraggConeHits.push_back(hit);

                //         if (tp) BraggConeEnergyTP += ph_near->ROISummedADC();
                //     }
                // }
            // }
        }

        // Truth Information
        LOG(mcp);
        if (mcp) {
            truPdg = mcp->PdgCode();
            truEndProcess = mcp->EndProcess();
            truStartPoint = ana::Point(mcp->Position().Vect());
            truEndPoint = ana::Point(mcp->EndPosition().Vect());
            truEndEnergy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV

            // if (geoDet == kPDVD)
            //     TrueDownward = mcp->Position(0).X() > mcp->EndPosition().X();
            // else if (geoDet == kPDHD)
            //     TrueDownward = mcp->Position(0).Y() > mcp->EndPosition().Y();

            ana::SortedHits sh_mcp = GetSortedHits_dirX(vph_mcp_mu, mcp->EndX() > mcp->Vx() ? 1 : -1);

            LOG(sh_mcp);
            if (sh_mcp) {
                truStartHit = GetHit(sh_mcp.start);
                truEndHit = GetHit(sh_mcp.end);
                truReg = sh_mcp.end_reg(geoDet);

                LOG(mcp_mi);
                if (mcp_mi) {
                    truHasMichel = (
                        geoHighX.isInside(mcp_mi->Position().Vect(), 20.F)
                        || geoLowX.isInside(mcp_mi->Position().Vect(), 20.F)
                    ) ? kHasMichelFiducial : (
                        geoHighX.isInside(mcp_mi->Position().Vect())
                        || geoLowX.isInside(mcp_mi->EndPosition().Vect())
                        ? kHasMichelInside
                        : kHasMichelOutside
                    );
                    miTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3;


                    PtrTrk pt_mi = ana::mcp2trk(mcp_mi, vpt_ev, clockData, fmp_trk2hit);
                    miTrackLength = pt_mi ? pt_mi->Length() : -1.F;
                    // PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
                    // MichelShowerLength = ps_mi ? ps_mi->Length() : -1.F;

                    float mu_end_angle = sh_mcp.end_reg(geoDet).theta(mcp->EndZ() > mcp->Vz() ? 1 : -1);
                    Hits bary_hits;
                    // for (PtrHit const& ph_mi : vph_mi) {
                    for (size_t i=0; i<vph_mi.size(); i++) {
                        PtrHit const& ph_mi = vph_mi[i];
                        float energyFrac = energyFracs_mi[i];

                        if (ph_mi->View() != geo::kW) continue;
                        Hit hit = GetHit(ph_mi);
                        miHits.push_back(hit);
                        miHitEnergyFrac.push_back(energyFrac);

                        if (hit.section != truEndHit.section) {
                            miHitMuonAngle.push_back(100);
                        } else {
                            float da = (hit.vec(fTick2cm) - truEndHit.vec(fTick2cm)).angle() - mu_end_angle;
                            da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                            miHitMuonAngle.push_back(da);
                        }

                        if (GetDistance(ph_mi, sh_mcp.end) > fBarycenterRadius) continue;
                        bary_hits.push_back(GetHit(ph_mi));
                    }
                    miHitEnergy = miHits.energy();

                    LOG(miBaryNHit);
                    if (bary_hits.size()) {
                        miBary = bary_hits.barycenter(fTick2cm);
                        ana::Vec2 end_bary = miBary - truEndHit.vec(fTick2cm);
                        miBaryAngle = end_bary.angle();
                        float da = miBaryAngle - mu_end_angle;
                        da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                        miBaryMuonAngle = da;

                        // if (fCone) {
                        //     // float angle = end_bary.angle();
                        //     miConeEnergy = 0;
                        //     for (PtrHit const& ph_mi : vph_mi) {
                        //         if (ph_mi->View() != geo::kW) continue;
                        //         float dist = GetDistance(ph_mi, sh_mcp.end);
                        //         if (dist > 30) continue;

                        //         ana::Vec2 end_hit = GetHit(ph_mi).vec(fTick2cm) - truEndHit.vec(fTick2cm);
                        //         float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                        //         if (dist > 5
                        //             && cosa < cos(30.F * TMath::DegToRad())
                        //         ) continue;

                        //         miKeyholeEnergy += ph_mi->ROISummedADC();

                        //         if (cosa < cos(30.F * TMath::DegToRad())) continue;

                        //         miConeEnergy += ph_mi->ROISummedADC();
                        //     }
                        // }
                    }
                }
            }
        }


        muTree->Fill();
        muIndex++;
        evNMuon++;
    } // end of loop over tracks
    evTree->Fill();
    evIndex++;
}

void ana::MichelAnalysis::beginJob() {}
void ana::MichelAnalysis::endJob() {}

void ana::MichelAnalysis::resetEvent() {
    evNMuon = 0;
    evMuonIndices.clear();
    evHits.clear();
}
void ana::MichelAnalysis::resetMuon() {
    // Hit
    muCathodeCrossing = -1;
    muAnodeCrossing = false;
    muStartHit = ana::Hit{};
    muEndHit = ana::Hit{};
    muEndHitX = kBadPoint;
    muEndHitY = kBadPoint;
    muStartHitX = kBadPoint;
    muStartHitY = kBadPoint;
    // muEndInX = false;
    // muEndInY = false;
    // muEndInZ = false;
    // muEndInT = false;
    muReg = ana::LinearRegression{};
    muHits.clear();
    // muEndSecHitdQds.clear();
    muTopHitdQds.clear();
    muBotHitdQds.clear();
    muSphereHits.clear();
    muSphereHitMuonAngle.clear();
    // muSphereHasShower = false;
    // muBaryHasShower = false;
    muSphereEnergy = -1.F;
    muSphereEnergyTP = -1.F;
    muBaryHasLongTrack = false;
    muSphereMaxShowerEnergy = -1.F;

    // Cone
    muBaryHits.clear();
    muBary = ana::Vec2{0,0};
    muBaryAngle = -10.F;
    muBaryMuonAngle = -10.F;
    // muConeHits.clear();
    // muConeEnergy = -1.F;
    // muConeEnergyTP = -1.F;
    // muKeyholeHits.clear();
    // muKeyholeEnergy = -1.F;
    // muKeyholeEnergyTP = -1.F;


    // Bragg
    // MIPdQds = 0;
    // BraggdQds = -1.F;
    // BraggEndHit = ana::Hit{};
    // BraggMuonHits.clear();
    // BraggSphereHits.clear();
    // BraggSphereHitMuonAngle.clear();
    // BraggSphereEnergy = -1.F;
    // BraggSphereEnergyTP = -1.F;

    // BraggBaryHits.clear();
    // BraggBary = ana::Vec2{0,0};
    // BraggBaryAngle = -10.F;
    // BraggBaryMuonAngle = -10.F;

    // NearbyBaryHits.clear();
    // NearbyBary = ana::Vec2{0,0};
    // NearbyBaryAngle = -10.F;
    // NearbyBaryMuonAngle = -10.F;
    // BraggConeHits.clear();
    // BraggConeEnergy = -1.F;
    // BraggConeEnergyTP = -1.F;
    // BraggKeyholeHits.clear();
    // BraggKeyholeEnergy = -1.F;
    // BraggKeyholeEnergyTP = -1.F;

    // Muon Truth
    truPdg = 0;
    truEndProcess = "";
    truStartPoint = ana::Point{};
    truEndPoint = ana::Point{};
    truEndEnergy = -1.F;
    // truDownward = false;
    truStartHit = ana::Hit{};
    truEndHit = ana::Hit{};
    truReg = ana::LinearRegression{};

    // Michel Truth
    truHasMichel = kNoMichel;
    miTrueEnergy = -1.F;
    miTrackLength = -1.F;
    // miShowerLength = -1.F;
    miHits.clear();
    miHitEnergyFrac.clear();
    miHitMuonAngle.clear();
    miHitEnergy = -1.F;

    miBaryNHit = 0;
    miBary = ana::Vec2{0,0};
    miBaryAngle = -10.F;
    miBaryMuonAngle = -10.F;
    // miConeEnergy = -1.F;
    // miKeyholeEnergy = -1.F;
}

DEFINE_ART_MODULE(ana::MichelAnalysis)