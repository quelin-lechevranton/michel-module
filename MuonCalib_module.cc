#include "utils.h"

namespace ana {
    class MuonCalib;
}

// using HitPtr = art::Ptr<recob::Hit>;
// using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
// using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

class ana::MuonCalib : public art::EDAnalyzer, private ana::MichelAnalyzer {
public:
    explicit MuonCalib(fhicl::ParameterSet const& p);
    MuonCalib(MuonCalib const&) = delete;
    MuonCalib(MuonCalib&&) = delete;
    MuonCalib& operator=(MuonCalib const&) = delete;
    MuonCalib& operator=(MuonCalib&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    ana::Bounds<float> wireWindow;
    ana::Bounds3D<float> geoHighX, geoLowX;
    float fCathodeGap; // cm

    bool fLog;
    float fFiducialLength;

    std::unordered_map<std::string, unsigned> track_content;

    TTree* tEvent;
    TTree* tTrack;

    unsigned evRun, evSubRun, evEvent;
    unsigned iEvent=0;
    unsigned EventNTrack;
    std::vector<unsigned> EventiTrack;
    ana::Hits EventHits;

    // Track information
    unsigned iTrack=0;

    float TrkLength;
    bool TrkIsUpright; // Supposition: Muon is downward
    ana::Point TrkStartPoint, TrkEndPoint;
    bool TrkCathodeCrossing;
    bool TrkAnodeCrossing;

    ana::Hits TrkHits;
    ana::Hit TrkStartHit, TrkEndHit;
    float TrkEndHitX;
    int TrkRegDirZ;
    ana::LinearRegression TrkReg;
    float TrkChi2;
    float TrkChi2PerNdof;
    bool TrkHitEndInVolumeX;

    std::vector<float> EndSecHitdQdx;
    std::vector<float> TopHitdQdx;
    std::vector<float> BotHitdQdx;

    int TruePdg;
    std::string TrueEndProcess;
    ana::Point MuonTrueStartPoint;
    ana::Point MuonTrueEndPoint;
    float MuonTrueEndEnergy;
    bool TrueDownward;

    void resetEvent();
    void resetTrack();
};

ana::MuonCalib::MuonCalib(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, MichelAnalyzer{p},
  fLog(p.get<bool>("Log", true)),
  fFiducialLength(p.get<float>("FiducialLength", 10.F)) // in cm
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
    
    track_content["total"] = 0;
    track_content["with hits"] = 0;
    track_content[">4 hits in one section"] = 0;
    track_content["r2 > 0.4"] = 0;
    track_content["a. anode crossing"] = 0;
    track_content["b. cathode crossing"] = 0;
    track_content["c. cathode + bot anode crossing"] = 0;
    track_content["d. anode and cathode crossing"] = 0;

    tEvent = asFile->make<TTree>("event", "");

    tEvent->Branch("Run", &evRun);
    tEvent->Branch("SubRun", &evSubRun);
    tEvent->Branch("Event", &evEvent);
    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NTrack", &EventNTrack);
    tEvent->Branch("iTrack", &EventiTrack);
    EventHits.SetBranches(tEvent, "Event");

    tTrack = asFile->make<TTree>("track", "");

    tTrack->Branch("Run", &evRun);
    tTrack->Branch("SubRun", &evSubRun);
    tTrack->Branch("Event", &evEvent);
    tTrack->Branch("iEvent", &iEvent);
    tTrack->Branch("iTrack", &iTrack);
    tTrack->Branch("iTrackInEvent", &EventNTrack);

    tTrack->Branch("TrkLength", &TrkLength);
    tTrack->Branch("TrkIsUpright", &TrkIsUpright);
    TrkStartPoint.SetBranches(tTrack, "Start");
    TrkEndPoint.SetBranches(tTrack, "End");
    tTrack->Branch("TrkCathodeCrossing", &TrkCathodeCrossing);
    tTrack->Branch("TrkAnodeCrossing", &TrkAnodeCrossing);

    TrkHits.SetBranches(tTrack, "");
    TrkStartHit.SetBranches(tTrack, "Start");
    TrkEndHit.SetBranches(tTrack, "End");
    tTrack->Branch("TrkEndHitX", &TrkEndHitX);
    tTrack->Branch("TrkRegDirZ", &TrkRegDirZ);
    TrkReg.SetBranches(tTrack, "TrkReg");
    tTrack->Branch("TrkChi2", &TrkChi2);
    tTrack->Branch("TrkChi2PerNdof", &TrkChi2PerNdof);
    tTrack->Branch("TrkHitEndInVolumeX", &TrkHitEndInVolumeX);

    tTrack->Branch("EndSecHitdQdx", &EndSecHitdQdx);
    tTrack->Branch("TopHitdQdx", &TopHitdQdx);
    tTrack->Branch("BotHitdQdx", &BotHitdQdx);

    tTrack->Branch("TruePdg", &TruePdg);
    tTrack->Branch("TrueEndProcess", &TrueEndProcess);
    MuonTrueStartPoint.SetBranches(tTrack, "MuonTrueStart");
    MuonTrueEndPoint.SetBranches(tTrack, "MuonTrueEnd");
    tTrack->Branch("MuonTrueEndEnergy", &MuonTrueEndEnergy);
    tTrack->Branch("TrueDownward", &TrueDownward);
}

void ana::MuonCalib::analyze(art::Event const& e) {
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

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);

    resetEvent();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();

    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            EventHits.push_back(GetHit(p_hit));

    track_content.clear();

    for (PtrTrk pt_ev : vpt_ev) {
        if (fLog) std::cout << "e" << iEvent << "t" << pt_ev->ID() << "\r" << std::flush;

        track_content["total"]++;
        resetTrack();

        VecPtrHit vph_trk = fmp_trk2hit.at(pt_ev.key());
        ASSERT(vph_trk.size())
        track_content["with hits"]++;

        TrkIsUpright =  IsUpright(*pt_ev);
        geo::Point_t Start = TrkIsUpright ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = TrkIsUpright ? pt_ev->End() : pt_ev->Start();

        TrkRegDirZ = End.Z() > Start.Z() ? 1 : -1;
        ana::SortedHits sh_trk = GetSortedHits(vph_trk, TrkRegDirZ);
        ASSERT(sh_trk)
        track_content[">4 hits in one section"]++;

        TrkStartHit = GetHit(sh_trk.start);
        TrkEndHit = GetHit(sh_trk.end);
        TrkReg = sh_trk.end_reg(geoDet);
        TrkChi2 = pt_ev->Chi2();
        TrkChi2PerNdof = pt_ev->Chi2PerNdof();

        ASSERT(TrkReg.r2 > 0.4)
        track_content["r2 > 0.4"]++;
        EventiTrack.push_back(iTrack);

        VecPtrHit vph_trk_endsec;
        for (PtrHit const& ph_trk : sh_trk.vph) {
            if (ph_trk->View() != geo::kW) continue;
            ana::Hit hit = GetHit(ph_trk);
            TrkHits.push_back(hit);

            if (hit.section != sh_trk.end_sec()) continue;
            vph_trk_endsec.push_back(ph_trk);
        }
        EndSecHitdQdx = GetdQdx(vph_trk_endsec, 6);

        TopHitdQdx = GetdQdx(sh_trk.vph.begin(), sh_trk.bot_it(), 6);
        BotHitdQdx = GetdQdx(sh_trk.bot_it(), sh_trk.vph.end(), 6);

        TrkCathodeCrossing = sh_trk.is_cc()
            && abs(sh_trk.cc.first->PeakTime()-sh_trk.cc.second->PeakTime())*fTick2cm < 2 * fCathodeGap;

        /* ASSUMS DOWNWARD MUONS FOR PDVD */
        if (geoDet == kPDVD)
            TrkAnodeCrossing = TrkStartHit.section < 4 
                && geoHighX.z.isInside(TrkStartHit.space, fFiducialLength)
                && wireWindow.isInside(TrkStartHit.tick, fFiducialLength/fTick2cm);
        else if (geoDet == kPDHD)
            TrkAnodeCrossing = 
                geoHighX.z.isInside(TrkStartHit.space, fFiducialLength)
                && wireWindow.isInside(TrkStartHit.tick, fFiducialLength/fTick2cm);

        if (TrkCathodeCrossing) {
            if (geoDet == kPDVD) {
                TrkEndHitX = geoLowX.x.max - abs(TrkEndHit.tick - sh_trk.cc.second->PeakTime()) * fTick2cm;
                TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
            } else if (geoDet == kPDHD) {
                int cc_sec = ana::tpc2sec[geoDet][sh_trk.cc.second->WireID().TPC];
                if (cc_sec == 0) {
                    TrkEndHitX = geoLowX.x.max - abs(TrkEndHit.tick - sh_trk.cc.second->PeakTime()) * fTick2cm;
                    TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
                } else if (cc_sec == 1) {
                    TrkEndHitX = geoHighX.x.min + abs(TrkEndHit.tick - sh_trk.cc.second->PeakTime()) * fTick2cm;
                    TrkHitEndInVolumeX = geoHighX.x.isInside(TrkEndHitX, fFiducialLength);
                }
            }
        } else {
            TrkHitEndInVolumeX = false;
        }

        if (TrkAnodeCrossing)
            track_content["a. anode crossing"]++;
        if (TrkCathodeCrossing)
            track_content["b. cathode crossing"]++;
        if (TrkCathodeCrossing && !TrkHitEndInVolumeX)
            track_content["c. cathode + bot anode crossing"]++;
        if (TrkAnodeCrossing && TrkCathodeCrossing)
            track_content["d. anode and cathode crossing"]++;


        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        if (mcp) {
            TruePdg = mcp->PdgCode();
            TrueEndProcess = mcp->EndProcess();
            MuonTrueStartPoint = ana::Point(mcp->Position().Vect());
            MuonTrueEndPoint = ana::Point(mcp->EndPosition().Vect());
            MuonTrueEndEnergy = (mcp->EndE() - mcp->Mass()) * 1e3;

            if (geoDet == kPDVD)
                TrueDownward = MuonTrueStartPoint.x > MuonTrueEndPoint.x;
            else if (geoDet == kPDHD)
                TrueDownward = MuonTrueStartPoint.y > MuonTrueEndPoint.y;
        }

        tTrack->Fill();
        iTrack++;
        EventNTrack++;
    }
    tEvent->Fill();
    iEvent++;
}

void ana::MuonCalib::beginJob() {}
void ana::MuonCalib::endJob() {
    std::cout << "\033[1;93m" "Track Selection Summary:" "\033[0m" << std::endl;
    std::cout << "  Total Tracks: " << track_content["total"] << std::endl;
    std::cout << "  With Hits: " << track_content["with hits"] << std::endl;
    std::cout << "  >4 Hits in One Section: " << track_content[">4 hits in one section"] << std::endl;
    std::cout << "  r2 > 0.4: " << track_content["r2 > 0.4"] << std::endl;
    std::cout << "  Anode Crossing: " << track_content["a. anode crossing"] << std::endl;
    std::cout << "  Cathode Crossing: " << track_content["b. cathode crossing"] << std::endl;
    std::cout << "  Cathode + Bottom Anode Crossing: " << track_content["c. cathode + bot anode crossing"] << std::endl;
    std::cout << "  Anode and Cathode Crossing: " << track_content["d. anode and cathode crossing"] << std::endl;
}

void ana::MuonCalib::resetEvent() {
    EventNTrack = 0;
    EventiTrack.clear();
    EventHits.clear();
}

void ana::MuonCalib::resetTrack() {
    TrkLength = 0.F;
    TrkIsUpright = false;

    TrkStartPoint = ana::Point();
    TrkEndPoint = ana::Point();
    TrkCathodeCrossing = false;
    TrkAnodeCrossing = false;

    TrkHits.clear();
    TrkStartHit = ana::Hit();
    TrkEndHit = ana::Hit();
    TrkEndHitX = -500.F;
    TrkRegDirZ = 0;
    TrkReg = ana::LinearRegression();
    TrkChi2 = 0.F;
    TrkChi2PerNdof = 0.F;
    TrkHitEndInVolumeX = false;

    EndSecHitdQdx.clear();
    TopHitdQdx.clear();
    BotHitdQdx.clear(); 

    TruePdg = 0;
    TrueEndProcess = "";
    MuonTrueStartPoint = ana::Point();
    MuonTrueEndPoint = ana::Point();
    MuonTrueEndEnergy = -1.F;
    TrueDownward = false;
}

DEFINE_ART_MODULE(ana::MuonCalib)