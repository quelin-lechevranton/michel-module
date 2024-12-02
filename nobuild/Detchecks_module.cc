////////////////////////////////////////////////////////////////////////
// Class:       Detchecks
// Plugin Type: analyzer
// File:        Detchecks_module.cc
//
// Generated at Wed Nov  6 10:31:28 2024 by Jeremy Quelin Lechevranton
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
#include "lardata/ArtDataHelper/TrackUtils.h"

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

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::min;
using std::max;


namespace ana {
  class Detchecks;
  struct Binning {
        int n;
        double min, max;
  };
}


class ana::Detchecks : public art::EDAnalyzer {
public:
    explicit Detchecks(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    Detchecks(Detchecks const&) = delete;
    Detchecks(Detchecks&&) = delete;
    Detchecks& operator=(Detchecks const&) = delete;
    Detchecks& operator=(Detchecks&&) = delete;

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

    // Conversion factors
    float fADCtoE = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    float fChannelPitch = 0.5; // cm/channel

};


ana::Detchecks::Detchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    iLogLevel(p.get<int>("LogLevel")),
    vvsProducts(p.get<vector<vector<string>>>("Products")),

    fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm

    fTrackLengthCut(p.get<float>("LengthCut")) // in cm
{

    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

}

void ana::Detchecks::analyze(art::Event const& e)
{

} // end analyze

void ana::Detchecks::beginJob()
{
    // auto const clockData = asDetClocks->DataForJob();
    // auto const detProp = asDetProp->DataForJob(clockData);

    if (iLogLevel >= iFlagDetails) {
        geo::CryostatID cryoid{0};
        cout << "\033[94m" << "Detchecks::beginJob: Detector dimension =========================================" << "\033[0m" << endl
             << "Number of channels: " << asWire->Nchannels() << endl
             << "Number of ticks: " << "???" << endl
             << "Cryostat coordinates: " << asGeo->Cryostat(cryoid).Min() << " - " << asGeo->Cryostat(cryoid).Max() << endl
             << "\tvolume: " << asGeo->Cryostat(cryoid).Width() << " x " << asGeo->Cryostat(cryoid).Height() << " x " << asGeo->Cryostat(cryoid).Length() << endl;
        geo::TPCID tpcid0{cryoid, 0}, tpcidN{cryoid, asGeo->NTPC(cryoid)-1};
        cout << "\tactive coordinates: " << asGeo->TPC(tpcid0).Min() << " - " << asGeo->TPC(tpcidN).Max() << endl
             << "\tactive volume: " << (asGeo->TPC(tpcid0).Max() - asGeo->TPC(tpcidN).Min()).X() << " x " << (asGeo->TPC(tpcid0).Max() - asGeo->TPC(tpcidN).Min()).Y() << " x " << (asGeo->TPC(tpcid0).Max() - asGeo->TPC(tpcidN).Min()).Z() << endl
             << "TPCs:" << endl; 
        for (unsigned int tpc=0; tpc < asGeo->NTPC(); tpc++) {
            geo::TPCID tpcid{cryoid, tpc};
            cout << "\tTPC#" << tpc << "\tcoordinates: " << asGeo->TPC(tpcid).Min() << " - " << asGeo->TPC(tpcid).Max() << endl
                 << "\t\tvolume: " << asGeo->TPC(tpcid).Width() << " x " << asGeo->TPC(tpcid).Height() << " x " << asGeo->TPC(tpcid).Length() << endl
                 << "\t\tactive volume: " << asGeo->TPC(tpcid).ActiveWidth() << " x " << asGeo->TPC(tpcid).ActiveHeight() << " x " << asGeo->TPC(tpcid).ActiveLength() << endl
                 << "\t\tPlanes:";
            for (unsigned int plane=0; plane < asWire->Nplanes(); plane++) {
                geo::PlaneID planeid{tpcid, plane};
                cout << "\t" << 'U'+plane << " #Wires: " << asWire->Nwires(planeid);
            } // end loop over Planes
            cout << endl;
        } // end loop over TPCs
        cout << "\033[94m" << "End of Detchecks::beginJob ======================================================" << "\033[0m" << endl;
    }
} // end beginJob


void ana::Detchecks::endJob()
{

} // end endJob


DEFINE_ART_MODULE(ana::Detchecks)
