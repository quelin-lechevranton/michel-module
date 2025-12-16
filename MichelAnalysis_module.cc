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
    bool fLog;
    bool fKeepAll;
    float fTrackLengthCut; // in cm
    float fFiducialLength; // in cm
    float fBarycenterRadius; // in cm
    float fMichelRadius; // in cm
    float fNearbyRadius; // in cm
    float fBodyDistance; // in cm
    unsigned fRegN;
    bool fCone;
    bool fBragg;
    unsigned fBraggN;
    float fBraggThreshold; // in MIP

    // Output Variables
    TTree *tEvent, *tMuon;

    // Event information
    unsigned evRun, evSubRun, evEvent;
    bool EventIsReal;
    unsigned iEvent=0;
    unsigned EventNMuon;
    std::vector<unsigned> EventiMuon;
    ana::Hits EventHits;

    // Track information
    unsigned iMuon=0;

    float TrkLength;
    float TrkChi2, TrkChi2PerNdof;
    bool TrkIsUpright; // Supposition: Muon is downward
    ana::Point TrkStartPoint, TrkEndPoint;
    bool TrkEndInVolumeYZ;

    // Hit from track information
    bool TrkHitError;
    int TrkHitCathodeCrossing;
    enum EnumCathodeCrossing { kNoCC, kHitOnBothSides, kAlignedHitOnBothSides };
    bool TrkHitAnodeCrossing;
    ana::Hit TrkStartHit, TrkEndHit;
    float TrkEndHitX;
    int TrkRegDirZ;
    ana::LinearRegression TrkReg;
    bool TrkHitEndInVolumeX;
    bool TrkHitEndInWindow;
    ana::Hits TrkHits;
    std::vector<float> PandoraTrkHitdQds;
    ana::Hits PandoraSphereHits;
    std::vector<float> PandoraSphereHitMuonAngle;
    // bool PandoraSphereHasShower;
    // bool PandoraBaryHasShower;
    bool PandoraBaryHasLongTrack;
    float PandoraSphereEnergy;
    float PandoraSphereEnergyTP;
    // Cone
    ana::Hits PandoraBaryHits;
    ana::Vec2 PandoraBary;
    float PandoraBaryAngle;
    float PandoraBaryMuonAngle;
    ana::Hits PandoraConeHits;
    float PandoraConeEnergy;
    float PandoraConeEnergyTP;
    ana::Hits PandoraKeyholeHits;
    float PandoraKeyholeEnergy;
    float PandoraKeyholeEnergyTP;


    // Bragg information
    int BraggError;
    enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
    float MIPdQds;
    float BraggdQds;
    ana::Hit BraggEndHit;
    ana::Hits BraggMuonHits;

    // Sphere
    ana::Hits BraggSphereHits;
    std::vector<float> BraggSphereHitMuonAngle;
    float BraggSphereEnergy;
    float BraggSphereEnergyTP;

    ana::Hits BraggBaryHits;
    ana::Vec2 BraggBary;
    float BraggBaryAngle;
    float BraggBaryMuonAngle;
    bool BraggBaryHasLongTrack;

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
    int TruePdg;
    std::string TrueEndProcess;
    ana::Point MuonTrueStartPoint;
    ana::Point MuonTrueEndPoint;
    float MuonTrueEndEnergy;
    bool TrueDownward;
    ana::Hit MuonTrueStartHit;
    ana::Hit MuonTrueEndHit;
    int MuonTrueRegDirZ;
    ana::LinearRegression MuonTrueReg;

    int TrueHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside, kHasMichelFiducial };
    float MichelTrueEnergy;
    float MichelTrackLength;
    // float MichelShowerLength;
    ana::Hits MichelHits;
    std::vector<float> MichelHitEnergyFrac;
    std::vector<float> MichelHitMuonAngle;
    float MichelHitEnergy;

    unsigned MichelBaryNHit;
    ana::Vec2 MichelBary;
    float MichelBaryAngle;
    float MichelBaryMuonAngle;
    float MichelConeEnergy;
    float MichelKeyholeEnergy;

    void resetEvent();
    void resetMuon();

    // HitPtrPair GetTrackEndsHits(
    //     HitPtrVec const& vp_hit,
    //     HitPtrPair *pp_cathode_crossing = nullptr,
    //     HitPtrVec *vp_tpc_crossing = nullptr,
    //     HitPtrVec *vp_sorted_hit = nullptr,
    //     geo::View_t view = geo::kW
    // );
};

ana::MichelAnalysis::MichelAnalysis(fhicl::ParameterSet const& p) : 
    EDAnalyzer{p}, 
    MichelAnalyzer{p},
    fLog(p.get<bool>("Log", true)),
    fKeepAll(p.get<bool>("KeepAll", false)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 20.F)), // in cm
    fFiducialLength(p.get<float>("FiducialLength", 10.F)), // in cm
    fBarycenterRadius(p.get<float>("BarycenterRadius", 10.F)), // in cm
    fMichelRadius(p.get<float>("MichelRadius", 20.F)), //in cm
    fNearbyRadius(p.get<float>("NearbyRadius", 40.F)), //in cm
    fBodyDistance(p.get<float>("BodyDistance", 20.F)), //in cm
    fRegN(p.get<unsigned>("RegN", 6)),
    fCone(p.get<bool>("Cone", false)),
    fBragg(p.get<bool>("Bragg", false)),
    fBraggN(p.get<unsigned>("BraggN", 6)),
    fBraggThreshold(p.get<float>("BraggThreshold", 1.7)) // in MIP dE/dx
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

    tEvent = asFile->make<TTree>("event","");

    tEvent->Branch("eventRun", &evRun);
    tEvent->Branch("eventSubRun", &evSubRun);
    tEvent->Branch("eventEvent", &evEvent);
    tEvent->Branch("isReal", &EventIsReal);
    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    EventHits.SetBranches(tEvent);

    tMuon = asFile->make<TTree>("muon","");

    // Event
    tMuon->Branch("eventRun", &evRun);
    tMuon->Branch("eventSubRun", &evSubRun);
    tMuon->Branch("eventEvent", &evEvent);
    tEvent->Branch("isReal", &EventIsReal);
    tMuon->Branch("iEvent", &iEvent);
    tMuon->Branch("iMuon", &iMuon);
    tMuon->Branch("iMuonInEvent", &EventNMuon);

    // Track
    tMuon->Branch("TrkLength", &TrkLength);
    tMuon->Branch("TrkChi2", &TrkChi2);
    tMuon->Branch("TrkChi2PerNdof", &TrkChi2PerNdof);
    tMuon->Branch("TrkIsUpright", &TrkIsUpright);
    TrkStartPoint.SetBranches(tMuon, "Start");
    TrkEndPoint.SetBranches(tMuon, "End");
    tMuon->Branch("TrkEndInVolumeYZ", &TrkEndInVolumeYZ);

    tMuon->Branch("MIPdQds", &MIPdQds);

    // Hit
    tMuon->Branch("TrkHitError", &TrkHitError);
    tMuon->Branch("TrkHitCathodeCrossing", &TrkHitCathodeCrossing);
    tMuon->Branch("TrkHitAnodeCrossing", &TrkHitAnodeCrossing);
    TrkStartHit.SetBranches(tMuon, "Start");
    TrkEndHit.SetBranches(tMuon, "End");
    tMuon->Branch("EndHitX", &TrkEndHitX);
    tMuon->Branch("RegDirZ", &TrkRegDirZ);
    TrkReg.SetBranches(tMuon, "");
    tMuon->Branch("TrkHitEndInVolumeX", &TrkHitEndInVolumeX);
    tMuon->Branch("TrkHitEndInWindow", &TrkHitEndInWindow);
    TrkHits.SetBranches(tMuon, "");
    tMuon->Branch("HitdQds", &PandoraTrkHitdQds);
    PandoraSphereHits.SetBranches(tMuon, "PandoraSphere");
    tMuon->Branch("PandoraSphereHitMuonAngle", &PandoraSphereHitMuonAngle);
    // tMuon->Branch("PandoraSphereHasShower", &PandoraSphereHasShower);
    // tMuon->Branch("PandoraBaryHasShower", &PandoraBaryHasShower);
    tMuon->Branch("PandoraSphereEnergy", &PandoraSphereEnergy); // ADC
    tMuon->Branch("PandoraSphereEnergyTP", &PandoraSphereEnergyTP); // ADC

    PandoraBaryHits.SetBranches(tMuon, "PandoraBary");
    PandoraBary.SetBranches(tMuon, "PandoraBary");
    tMuon->Branch("PandoraBaryAngle", &PandoraBaryAngle);
    tMuon->Branch("PandoraBaryMuonAngle", &PandoraBaryMuonAngle);
    tMuon->Branch("PandoraBaryHasLongTrack", &PandoraBaryHasLongTrack);

    if (fCone) {
        // Cone
        PandoraConeHits.SetBranches(tMuon, "PandoraCone");
        tMuon->Branch("PandoraConeEnergy", &PandoraConeEnergy); // ADC
        tMuon->Branch("PandoraConeEnergyTP", &PandoraConeEnergyTP); // ADC
        PandoraKeyholeHits.SetBranches(tMuon, "PandoraCone");
        tMuon->Branch("PandoraKeyholeEnergy", &PandoraKeyholeEnergy); // ADC
        tMuon->Branch("PandoraKeyholeEnergyTP", &PandoraKeyholeEnergyTP); // ADC
    }

    if (fBragg) {
        // Bragg
        tMuon->Branch("BraggError", &BraggError);
        tMuon->Branch("BraggdQds", &BraggdQds);
        BraggEndHit.SetBranches(tMuon, "BraggEnd");
        BraggMuonHits.SetBranches(tMuon, "BraggMuon");

        // Sphere
        BraggSphereHits.SetBranches(tMuon, "BraggSphere");
        tMuon->Branch("BraggSphereHitMuonAngle", &BraggSphereHitMuonAngle);
        tMuon->Branch("BraggSphereEnergy", &BraggSphereEnergy); // ADC
        tMuon->Branch("BraggSphereEnergyTP", &BraggSphereEnergyTP); // ADC

        BraggBaryHits.SetBranches(tMuon, "BraggBary");
        BraggBary.SetBranches(tMuon, "BraggBary");
        tMuon->Branch("BraggBaryAngle", &BraggBaryAngle);
        tMuon->Branch("BraggBaryMuonAngle", &BraggBaryMuonAngle);
        tMuon->Branch("BraggBaryHasLongTrack", &BraggBaryHasLongTrack);

        // if (fCone) {
        //     // Cone
        //     NearbyBaryHits.SetBranches(tMuon, "NearbyBary");
        //     NearbyBary.SetBranches(tMuon, "NearbyBary");
        //     tMuon->Branch("NearbyBaryAngle", &NearbyBaryAngle);
        //     tMuon->Branch("NearbyBaryMuonAngle", &NearbyBaryMuonAngle);
        //     BraggConeHits.SetBranches(tMuon, "BraggCone");
        //     tMuon->Branch("BraggConeEnergy", &BraggConeEnergy); // ADC
        //     tMuon->Branch("BraggConeEnergyTP", &BraggConeEnergyTP); // ADC
        //     BraggKeyholeHits.SetBranches(tMuon, "BraggCone");
        //     tMuon->Branch("BraggKeyholeEnergy", &BraggKeyholeEnergy); // ADC
        //     tMuon->Branch("BraggKeyholeEnergyTP", &BraggKeyholeEnergyTP); // ADC
        // }
    }

    // Truth
    tMuon->Branch("TruePdg", &TruePdg);
    tMuon->Branch("TrueEndProcess", &TrueEndProcess);
    MuonTrueStartPoint.SetBranches(tMuon, "TrueStart");
    MuonTrueEndPoint.SetBranches(tMuon, "TrueEnd");
    tMuon->Branch("TrueEndEnergy", &MuonTrueEndEnergy);
    tMuon->Branch("TrueDownward", &TrueDownward);
    MuonTrueStartHit.SetBranches(tMuon, "TrueStart");
    MuonTrueEndHit.SetBranches(tMuon, "TrueEnd");
    tMuon->Branch("TrueRegDirZ", &MuonTrueRegDirZ);
    MuonTrueReg.SetBranches(tMuon, "True");

    tMuon->Branch("TrueHasMichel", &TrueHasMichel);
    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy); // MeV
    tMuon->Branch("MichelTrackLength", &MichelTrackLength); // cm
    // tMuon->Branch("MichelShowerLength", &MichelShowerLength); // cm
    MichelHits.SetBranches(tMuon, "Michel");
    tMuon->Branch("MichelHitEnergyFrac", &MichelHitEnergyFrac);
    tMuon->Branch("MichelHitMuonAngle", &MichelHitMuonAngle);
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy); // ADC

    tMuon->Branch("MichelBaryNHit", &MichelBaryNHit);
    MichelBary.SetBranches(tMuon, "MichelBary");
    tMuon->Branch("MichelBaryAngle", &MichelBaryAngle); // rad
    tMuon->Branch("MichelBaryMuonAngle", &MichelBaryMuonAngle); // rad
    tMuon->Branch("MichelConeEnergy", &MichelConeEnergy); // ADC
    tMuon->Branch("MichelKeyholeEnergy", &MichelKeyholeEnergy); // ADC


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

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    // art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    // art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    resetEvent();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    EventIsReal = e.isRealData();

    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            EventHits.push_back(GetHit(p_hit));

    // loop over tracks to find stopping muons
    for (PtrTrk const& pt_ev : vpt_ev) {
        if (fLog) std::cout << "e" << iEvent << "t" << pt_ev->ID() << "\r" << std::flush;
        resetMuon();

        VecPtrHit vph_mu = fmp_trk2hit.at(pt_ev.key());
        ASSERT(vph_mu.size())

        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        simb::MCParticle const* mcp_mi = nullptr;
        VecPtrHit vph_mcp_mu, vph_mi;
        std::vector<float> energyFracs_mi;
        if (mcp) {
            mcp_mi = GetMichelMCP(mcp);
            vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true, &energyFracs_mi);
        }

        if (fLog) std::cout << "\t" "\033[1;93m" "e" << iEvent << "m" << EventNMuon << " (" << iMuon << ")" "\033[0m" << std::endl;
        EventiMuon.push_back(iMuon);

        // ============================
        // Dump basic track information
        TrkLength = pt_ev->Length();
        TrkChi2 = pt_ev->Chi2();
        TrkChi2PerNdof = pt_ev->Chi2PerNdof();

        LOG(TrkLength >= fTrackLengthCut);
        if (!fKeepAll && TrkLength < fTrackLengthCut) continue;

        TrkIsUpright =  IsUpright(*pt_ev);
        geo::Point_t Start = TrkIsUpright ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = TrkIsUpright ? pt_ev->End() : pt_ev->Start();
        TrkStartPoint = ana::Point(Start);
        TrkEndPoint = ana::Point(End);

        TrkEndInVolumeYZ = geoHighX.isInsideYZ(End, fFiducialLength);

        LOG(TrkEndInVolumeYZ);
        if (!fKeepAll && !TrkEndInVolumeYZ) continue;


        // ===================================================
        // We need some more information on the track, mainly:
        // - which is the last point? (ASSUMING DOWNWARD TRACK, for cosmic muons)
        // - does the track cross the cathode? (for t0)
        // - some quality tags...?
        // For that, we make linear regressions on each side of the cathode
        // on the hits associated to the track given by Pandora

        // SUPPOSITION: Muon is downward
        TrkRegDirZ = End.Z() > Start.Z() ? 1 : -1;
        ana::SortedHits sh_mu = GetSortedHits(vph_mu, TrkRegDirZ);
        TrkHitError = !sh_mu;


        LOG(!TrkHitError);
        if (!fKeepAll && TrkHitError) continue;
        if (!TrkHitError) {
            TrkStartHit = GetHit(sh_mu.start);
            TrkEndHit = GetHit(sh_mu.end);
            TrkReg = sh_mu.end_reg(geoDet);
            TrkHitEndInWindow = wireWindow.isInside(TrkEndHit.tick, fFiducialLength / fTick2cm);

            LOG(TrkReg.r2 >= 0.4);
            if (!fKeepAll && TrkReg.r2 < 0.4) continue;
            LOG(TrkHitEndInWindow);
            if (!fKeepAll && !TrkHitEndInWindow) continue;

            if (sh_mu.is_cc()) {
                if (abs(sh_mu.cc.first->PeakTime()-sh_mu.cc.second->PeakTime())*fTick2cm < 3 * fCathodeGap)
                    TrkHitCathodeCrossing = kAlignedHitOnBothSides;
                else
                    TrkHitCathodeCrossing = kHitOnBothSides;
            } else
                TrkHitCathodeCrossing = kNoCC;

            /* ASSUMS DOWNWARD MUONS FOR PDVD */
            if (geoDet == kPDVD)
                TrkHitAnodeCrossing = TrkStartHit.section < 4
                    && geoHighX.z.isInside(TrkStartHit.space, fFiducialLength)
                    && wireWindow.isInside(TrkStartHit.tick, fFiducialLength/fTick2cm);
            else if (geoDet == kPDHD)
                TrkHitAnodeCrossing =
                    geoHighX.z.isInside(TrkStartHit.space, fFiducialLength)
                    && wireWindow.isInside(TrkStartHit.tick, fFiducialLength/fTick2cm);

            LOG(TrkHitCathodeCrossing != kNoCC || TrkHitAnodeCrossing);
            if (!fKeepAll && (!TrkHitAnodeCrossing && !TrkHitCathodeCrossing)) continue;

            LOG(TrkHitCathodeCrossing != kAlignedHitOnBothSides);
            LOG(TrkHitAnodeCrossing);
            if (TrkHitCathodeCrossing == kAlignedHitOnBothSides) {
                if (geoDet == kPDVD) {
                    TrkEndHitX = geoLowX.x.max - abs(TrkEndHit.tick - sh_mu.cc.second->PeakTime()) * fTick2cm;
                    TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
                } else if (geoDet == kPDHD) {
                    int cc_sec = ana::tpc2sec.at(geoDet).at(sh_mu.cc.second->WireID().TPC);
                    if (cc_sec == 0) {
                        TrkEndHitX = geoLowX.x.max - abs(TrkEndHit.tick - sh_mu.cc.second->PeakTime()) * fTick2cm;
                        TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
                    } else if (cc_sec == 1) {
                        TrkEndHitX = geoHighX.x.min + abs(TrkEndHit.tick - sh_mu.cc.second->PeakTime()) * fTick2cm;
                        TrkHitEndInVolumeX = geoHighX.x.isInside(TrkEndHitX, fFiducialLength);
                    }
                }
            } else if (TrkHitCathodeCrossing == kNoCC && TrkHitAnodeCrossing) {
                if (geoDet == kPDVD) {
                    TrkEndHitX = geoHighX.x.max - abs(TrkEndHit.tick - TrkStartHit.tick) * fTick2cm;
                    TrkHitEndInVolumeX = geoHighX.x.isInside(TrkEndHitX, fFiducialLength);
                } else if (geoDet == kPDHD) {
                    if (TrkStartHit.section == 0) {
                        TrkEndHitX = geoLowX.x.min + abs(TrkEndHit.tick - TrkStartHit.tick) * fTick2cm;
                        TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
                    } else if (TrkStartHit.section == 1) {
                        TrkEndHitX = geoHighX.x.max - abs(TrkEndHit.tick - TrkStartHit.tick) * fTick2cm;
                        TrkHitEndInVolumeX = geoHighX.x.isInside(TrkEndHitX, fFiducialLength);
                    }
                }
            }

            LOG(TrkHitEndInVolumeX);
            if (!fKeepAll && !TrkHitEndInVolumeX) continue;

            for (PtrHit const& ph_mu : sh_mu.vph) {
                if (ph_mu->View() != geo::kW) continue;
                ana::Hit hit = GetHit(ph_mu);
                TrkHits.push_back(hit);
            }

            /* COMPARE TO AGNOCHECKS
            HitPtrPair ends = GetTrackEndsHits(vph_mu);
            if (!LOG(ends.first && ends.second)) continue;
            bool increasing_z = IsUpright(*pt_ev)
                ? pt_ev->End().Z() > pt_ev->Start().Z()
                : pt_ev->Start().Z() > pt_ev->End().Z();
            int dir_z = increasing_z ? 1 : -1;
            float fz = GetSpace(ends.first->WireID());
            float sz = GetSpace(ends.second->WireID());
            HitPtr end = (sz-fz) * dir_z > 0 ? ends.second : ends.first;
            TrkEndHit = GetHit(end);
            */

            VecPtrHit vph_ev_endsec;
            for (PtrHit const& ph_ev : vph_ev) {
                if (ph_ev->View() != geo::kW) continue;
                int sec = ana::tpc2sec.at(geoDet).at(ph_ev->WireID().TPC);
                if (sec != sh_mu.end_sec()) continue;
                vph_ev_endsec.push_back(ph_ev);
            }

            // integrate charges around muon endpoint
            PandoraSphereEnergy = 0;
            PandoraSphereEnergyTP = 0;
            PandoraBaryHasLongTrack = false;
            for (PtrHit const& ph_ev : vph_ev_endsec) {
                float dist = GetDistance(ph_ev, sh_mu.end);
                if (dist > fMichelRadius) continue;

                PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                if (pt_hit && pt_hit->Length() > fTrackLengthCut) {
                    if (dist < fBarycenterRadius) PandoraBaryHasLongTrack = true;
                    continue;
                }

                // PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
                // PandoraSphereHasShower = bool(ps_hit);

                if (dist < fBarycenterRadius) {
                    PandoraBaryHits.push_back(GetHit(ph_ev));
                }

                ana::Hit hit = GetHit(ph_ev);
                PandoraSphereEnergy += ph_ev->ROISummedADC();
                PandoraSphereHits.push_back(hit);

                float da = (hit.vec(fTick2cm) - TrkEndHit.vec(fTick2cm)).angle() - TrkReg.theta(TrkRegDirZ);
                da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                PandoraSphereHitMuonAngle.push_back(da);

                if (std::find_if(
                    vph_mi.begin(), vph_mi.end(),
                    [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                ) != vph_mi.end())
                    PandoraSphereEnergyTP = ph_ev->ROISummedADC();
            }

            // Cone
            LOG(PandoraBaryHits.size());
            if (PandoraBaryHits.size()) {
                PandoraBary = PandoraBaryHits.barycenter(fTick2cm);
                ana::Vec2 end_bary = PandoraBary - TrkEndHit.vec(fTick2cm);
                PandoraBaryAngle = end_bary.angle();
                float da = PandoraBaryAngle - TrkReg.theta(TrkRegDirZ);
                da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                PandoraBaryMuonAngle = da;

                if (fCone) {
                    // float angle = end_bary.angle();
                    PandoraConeEnergy = 0;
                    PandoraConeEnergyTP = 0;
                    for (PtrHit const& ph_ev : vph_ev_endsec) {
                        ana::Hit hit = GetHit(ph_ev);
                        PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                        if (pt_hit && pt_hit->Length() > fTrackLengthCut) continue;

                        float dist = GetDistance(ph_ev, sh_mu.end);
                        if (dist > 30) continue;

                        ana::Vec2 end_hit = hit.vec(fTick2cm) - TrkEndHit.vec(fTick2cm);
                        float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                        if (dist > 5
                            && cosa < cos(30.F * TMath::DegToRad())
                        ) continue;

                        bool tp = std::find_if(
                            vph_mi.begin(), vph_mi.end(),
                            [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                        ) != vph_mi.end();

                        PandoraKeyholeEnergy += ph_ev->ROISummedADC();
                        PandoraKeyholeHits.push_back(hit);

                        if (tp) PandoraKeyholeEnergyTP += ph_ev->ROISummedADC();

                        if (cosa < cos(30.F * TMath::DegToRad())) continue;

                        PandoraConeEnergy += ph_ev->ROISummedADC();
                        PandoraConeHits.push_back(hit);

                        if (tp) PandoraConeEnergyTP += ph_ev->ROISummedADC();
                    }
                }
            }

            // End dQds
            PandoraTrkHitdQds = GetdQds(sh_mu.endsec_it(), sh_mu.vph.end(), fRegN);

            LOG(!PandoraTrkHitdQds.empty());
            if (!fKeepAll && PandoraTrkHitdQds.empty()) continue;



            MIPdQds = 0;


            if (!fBragg || std::distance(sh_mu.endsec_it(), sh_mu.vph.end()) < 2 * fBraggN) {
                BraggError = true;
            } else {
                BraggError = false;

                VecPtrHit vph_mu_bragg(sh_mu.endsec_it(), sh_mu.vph.end() - 2*fBraggN);
                VecPtrHit vph_mu_tail(sh_mu.vph.end() - 2*fBraggN, sh_mu.vph.end());

                int n = fBraggN;
                while (n--) {
                    float min_dist = std::numeric_limits<float>::max();
                    PtrHit closest_hit;
                    for (PtrHit const& ph_ev : vph_ev_endsec) {
                        // not already in muon
                        if (std::find_if(
                            vph_mu_bragg.begin(), vph_mu_bragg.end(),
                            [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
                        ) != vph_mu_bragg.end()) continue;
                        if (std::find_if(
                            vph_mu_tail.begin(), vph_mu_tail.end(),
                            [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
                        ) != vph_mu_tail.end()) continue;

                        float dist = GetDistance(ph_ev, vph_mu_tail.back());
                        if (dist < min_dist) {
                            min_dist = dist;
                            closest_hit = ph_ev;
                        }
                    }
                    if (closest_hit) vph_mu_tail.push_back(closest_hit);
                    else break;
                }

                // unsigned i_dQds_max;
                // std::vector<float> bragg_dQds = GetdQds(vph_mu_tail, fRegN, &i_dQds_max);

                // BraggdQds = bragg_dQds[i_dQds_max];
                // for (unsigned i=0; i<=i_dQds_max; i++) {
                //     vph_mu_bragg.push_back(vph_mu_tail[i]);
                // }

                VecPtrHit::iterator it_dQ_max = std::max_element(
                    vph_mu_tail.begin(), vph_mu_tail.end(),
                    [&](PtrHit const& a, PtrHit const& b) -> bool {
                        return a->ROISummedADC() < b->ROISummedADC();
                    }
                );

                MIPdQds = std::accumulate(
                    sh_mu.endsec_it(), sh_mu.vph.end() - fBraggN, 0.F,
                    [&](float sum, PtrHit const& ph) -> float {
                        return sum + ph->ROISummedADC();
                    }
                ) / (sh_mu.vph.end() - fBraggN - sh_mu.endsec_it());
                BraggdQds = (*it_dQ_max)->ROISummedADC();
                for (VecPtrHit::iterator it=vph_mu_tail.begin(); it!=it_dQ_max+1; ++it)
                    vph_mu_bragg.push_back(*it);



                BraggEndHit = GetHit(vph_mu_bragg.back());
                for (PtrHit const& ph_mu : vph_mu_bragg)
                    BraggMuonHits.push_back(GetHit(ph_mu));

                BraggBaryHasLongTrack = false;
                for (PtrHit const& ph_ev : vph_ev_endsec) {
                    double dist = GetDistance(ph_ev, vph_mu_bragg.back());
                    if (dist > fMichelRadius) continue;

                    PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                    // not from other long track
                    if (pt_hit
                        && pt_hit.key() != pt_ev.key()
                        && pt_hit->Length() > fMichelRadius
                    ) {
                        if (dist < fBarycenterRadius) BraggBaryHasLongTrack = true;
                        continue;
                    }

                    if (std::find_if(
                        vph_mu_bragg.begin(), vph_mu_bragg.end(),
                        [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
                    ) != vph_mu_bragg.end()) continue;

                    if (dist < fBarycenterRadius) {
                        BraggBaryHits.push_back(GetHit(ph_ev));
                    }

                    // PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
                    // BraggSphereHasShower = bool(ps_hit);

                    ana::Hit hit = GetHit(ph_ev);
                    BraggSphereEnergy += ph_ev->ROISummedADC();
                    BraggSphereHits.push_back(hit);

                    float da = (hit.vec(fTick2cm) - TrkEndHit.vec(fTick2cm)).angle() - TrkReg.theta(TrkRegDirZ);
                    da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                    BraggSphereHitMuonAngle.push_back(da);

                    if (std::find_if(
                        vph_mi.begin(), vph_mi.end(),
                        [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                    ) != vph_mi.end())
                        BraggSphereEnergyTP = ph_ev->ROISummedADC();
                }

                LOG(BraggBaryHits.size());
                if (BraggBaryHits.size()) {
                    BraggBary = BraggBaryHits.barycenter(fTick2cm);
                    ana::Vec2 end_bary = BraggBary - BraggEndHit.vec(fTick2cm);
                    BraggBaryAngle = end_bary.angle();
                    float da = BraggBaryAngle - TrkReg.theta(TrkRegDirZ);
                    da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                    BraggBaryMuonAngle = da;
                }

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
                //     float da = NearbyBaryAngle - TrkReg.theta(TrkRegDirZ);
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
            }
        }

        // Truth Information
        LOG(mcp);
        if (mcp) {
            TruePdg = mcp->PdgCode();
            TrueEndProcess = mcp->EndProcess();
            MuonTrueStartPoint = ana::Point(mcp->Position().Vect());
            MuonTrueEndPoint = ana::Point(mcp->EndPosition().Vect());
            MuonTrueEndEnergy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV

            if (geoDet == kPDVD)
                TrueDownward = mcp->Position(0).X() > mcp->EndPosition().X();
            else if (geoDet == kPDHD)
                TrueDownward = mcp->Position(0).Y() > mcp->EndPosition().Y();

            MuonTrueRegDirZ = mcp->EndZ() > mcp->Vz() ? 1 : -1;
            ana::SortedHits sh_mcp = GetSortedHits(vph_mcp_mu, MuonTrueRegDirZ);

            LOG(sh_mcp);
            if (sh_mcp) {
                MuonTrueStartHit = GetHit(sh_mcp.start);
                MuonTrueEndHit = GetHit(sh_mcp.end);
                MuonTrueReg = sh_mcp.end_reg(geoDet);

                LOG(mcp_mi);
                if (mcp_mi) {
                    TrueHasMichel = (
                        geoHighX.isInside(mcp_mi->Position().Vect(), 20.F)
                        || geoLowX.isInside(mcp_mi->Position().Vect(), 20.F)
                    ) ? kHasMichelFiducial : (
                        geoHighX.isInside(mcp_mi->Position().Vect())
                        || geoLowX.isInside(mcp_mi->EndPosition().Vect())
                        ? kHasMichelInside
                        : kHasMichelOutside
                    );
                    MichelTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3;


                    PtrTrk pt_mi = ana::mcp2trk(mcp_mi, vpt_ev, clockData, fmp_trk2hit);
                    MichelTrackLength = pt_mi ? pt_mi->Length() : -1.F;
                    // PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
                    // MichelShowerLength = ps_mi ? ps_mi->Length() : -1.F;

                    float mu_end_angle = sh_mcp.end_reg(geoDet).theta(MuonTrueRegDirZ);
                    // for (PtrHit const& ph_mi : vph_mi) {
                    for (size_t i=0; i<vph_mi.size(); i++) {
                        PtrHit const& ph_mi = vph_mi[i];
                        float energyFrac = energyFracs_mi[i];

                        if (ph_mi->View() != geo::kW) continue;
                        Hit hit = GetHit(ph_mi);
                        MichelHits.push_back(hit);
                        MichelHitEnergyFrac.push_back(energyFrac);

                        if (hit.section != MuonTrueEndHit.section) {
                            MichelHitMuonAngle.push_back(100);
                        } else {
                            float mu_hit_angle = (hit.vec(fTick2cm) - MuonTrueEndHit.vec(fTick2cm)).angle() - mu_end_angle;
                            mu_hit_angle = abs(mu_hit_angle) > TMath::Pi()
                                ? mu_hit_angle - (mu_hit_angle>0 ? 1 : -1) * 2 * TMath::Pi()
                                : mu_hit_angle;
                            MichelHitMuonAngle.push_back(mu_hit_angle);
                        }
                    }
                    MichelHitEnergy = MichelHits.energy();

                    // Cone
                    Hits bary_hits;
                    for (PtrHit const& ph_mi : vph_mi) {
                        if (ph_mi->View() != geo::kW) continue;
                        if (GetDistance(ph_mi, sh_mcp.end) > 10) continue;
                        bary_hits.push_back(GetHit(ph_mi));
                    }
                    MichelBaryNHit = bary_hits.size();

                    LOG(MichelBaryNHit);
                    if (bary_hits.size()) {
                        MichelBary = bary_hits.barycenter(fTick2cm);
                        ana::Vec2 end_bary = MichelBary - MuonTrueEndHit.vec(fTick2cm);
                        MichelBaryAngle = end_bary.angle();
                        float da = MichelBaryAngle - mu_end_angle;
                        da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                        MichelBaryMuonAngle = da;

                        if (fCone) {
                            // float angle = end_bary.angle();
                            MichelConeEnergy = 0;
                            for (PtrHit const& ph_mi : vph_mi) {
                                if (ph_mi->View() != geo::kW) continue;
                                float dist = GetDistance(ph_mi, sh_mcp.end);
                                if (dist > 30) continue;

                                ana::Vec2 end_hit = GetHit(ph_mi).vec(fTick2cm) - MuonTrueEndHit.vec(fTick2cm);
                                float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                                if (dist > 5
                                    && cosa < cos(30.F * TMath::DegToRad())
                                ) continue;

                                MichelKeyholeEnergy += ph_mi->ROISummedADC();

                                if (cosa < cos(30.F * TMath::DegToRad())) continue;

                                MichelConeEnergy += ph_mi->ROISummedADC();
                            }
                        }
                    }
                }
            }
        }


        tMuon->Fill();
        iMuon++;
        EventNMuon++;
    } // end of loop over tracks
    tEvent->Fill();
    iEvent++;
}

void ana::MichelAnalysis::beginJob() {}
void ana::MichelAnalysis::endJob() {}

void ana::MichelAnalysis::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
}
void ana::MichelAnalysis::resetMuon() {
    // Hit
    TrkHitCathodeCrossing = -1;
    TrkHitAnodeCrossing = false;
    TrkStartHit = ana::Hit{};
    TrkEndHit = ana::Hit{};
    TrkEndHitX = -500.F;
    TrkHitEndInVolumeX = false;
    TrkHitEndInWindow = false;
    TrkRegDirZ = 0;
    TrkReg = ana::LinearRegression{};
    TrkHits.clear();
    PandoraTrkHitdQds.clear();
    PandoraSphereHits.clear();
    PandoraSphereHitMuonAngle.clear();
    // PandoraSphereHasShower = false;
    // PandoraBaryHasShower = false;
    PandoraSphereEnergy = -1.F;
    PandoraSphereEnergyTP = -1.F;

    // Cone
    PandoraBaryHits.clear();
    PandoraBary = ana::Vec2{0,0};
    PandoraBaryAngle = -10.F;
    PandoraBaryMuonAngle = -10.F;
    PandoraConeHits.clear();
    PandoraConeEnergy = -1.F;
    PandoraConeEnergyTP = -1.F;
    PandoraKeyholeHits.clear();
    PandoraKeyholeEnergy = -1.F;
    PandoraKeyholeEnergyTP = -1.F;


    // Bragg
    MIPdQds = 0;
    BraggdQds = -1.F;
    BraggEndHit = ana::Hit{};
    BraggMuonHits.clear();
    BraggSphereHits.clear();
    BraggSphereHitMuonAngle.clear();
    BraggSphereEnergy = -1.F;
    BraggSphereEnergyTP = -1.F;

    BraggBaryHits.clear();
    BraggBary = ana::Vec2{0,0};
    BraggBaryAngle = -10.F;
    BraggBaryMuonAngle = -10.F;

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
    TruePdg = 0;
    TrueEndProcess = "";
    MuonTrueStartPoint = ana::Point{};
    MuonTrueEndPoint = ana::Point{};
    MuonTrueEndEnergy = -1.F;
    TrueDownward = false;
    MuonTrueStartHit = ana::Hit{};
    MuonTrueEndHit = ana::Hit{};
    MuonTrueRegDirZ = 0;
    MuonTrueReg = ana::LinearRegression{};

    // Michel Truth
    TrueHasMichel = kNoMichel;
    MichelTrueEnergy = -1.F;
    MichelTrackLength = -1.F;
    // MichelShowerLength = -1.F;
    MichelHits.clear();
    MichelHitEnergyFrac.clear();
    MichelHitMuonAngle.clear();
    MichelHitEnergy = -1.F;

    MichelBaryNHit = 0;
    MichelBary = ana::Vec2{0,0};
    MichelBaryAngle = -10.F;
    MichelBaryMuonAngle = -10.F;
    MichelConeEnergy = -1.F;
    MichelKeyholeEnergy = -1.F;
}

DEFINE_ART_MODULE(ana::MichelAnalysis)
