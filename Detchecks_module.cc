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

// #include "nusimdata/SimulationBase/MCTruth.h"
// #include "nusimdata/SimulationBase/MCParticle.h"
// #include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
// #include "lardataobj/RecoBase/Wire.h"
// #include "lardataobj/RecoBase/Track.h"
// #include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "larcorealg/Geometry/Exceptions.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ProtoDUNE includes
// #include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
// #include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
// #include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

// ROOT includes


// std includes
#include <iostream>
#include <sstream>
// #include <vector>
// #include <string>
// #include <iterator>

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

    geo::WireID GetWireID(geo::Point_t const& P, geo::View_t plane);
    raw::ChannelID_t GetChannel(geo::Point_t const& P, geo::View_t plane);

};


ana::Detchecks::Detchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}
    // iLogLevel(p.get<int>("LogLevel")),

    // fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    // fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm
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
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);

    geo::CryostatID cryoid{0};
    std::cout << "\033[93m" << "Detchecks::beginJob: Detector dimension =========================================" << "\033[0m" << std::endl
        << "Detector name: " << asGeo->DetectorName() << std::endl
        << "Number of channels: " << asWire->Nchannels() << std::endl
        << "Number of ticks: " << detProp.ReadOutWindowSize() << std::endl
        << "Channel pitch: " << "(0.5)" << " cm/channel" << std::endl
        << "Sampling rate: " << detinfo::sampling_rate(clockData) << " ns/tick" << std::endl
        << "Drift velocity: " << detProp.DriftVelocity() << " cm/µs" << std::endl;

    std::cout << "Cryostat:" << std::endl
        << "\t" << "coordinates: " << asGeo->Cryostat(cryoid).Min() 
                                    << " - " << asGeo->Cryostat(cryoid).Max() 
                                    << std::endl
        << "\t" << "volume: " << asGeo->Cryostat(cryoid).Width() << " x " << asGeo->Cryostat(cryoid).Height() << " x " << asGeo->Cryostat(cryoid).Length() << std::endl;

    geo::TPCID tpcid0{cryoid, 0};
    geo::TPCID tpcidN{cryoid, asGeo->NTPC(cryoid)-1};

    geo::TPCGeo tpcgeo0 = asGeo->TPC(tpcid0);
    geo::TPCGeo tpcgeoN = asGeo->TPC(tpcidN);
    std::cout << "\t" << "active coordinates: " << tpcgeo0.Min() << " - " << tpcgeoN.Max() << std::endl;
    geo::Vector_t diag = tpcgeoN.Max() - tpcgeo0.Min();
    std::cout << "\t" << "active volume: "  << diag.X() << " x " << diag.Y() << " x " << diag.Z() << std::endl;

    std::cout << "TPCs:" << std::endl; 
    for (unsigned int tpc=0; tpc < asGeo->NTPC(); tpc++) {
        geo::TPCID tpcid{cryoid, tpc};
        std::cout << "\t" << "TPC#" << tpc << "\tcoordinates: " << asGeo->TPC(tpcid).Min() << " - " << asGeo->TPC(tpcid).Max() << std::endl
                << "\t\t" << "volume: " << asGeo->TPC(tpcid).Width() << " x " << asGeo->TPC(tpcid).Height() << " x " << asGeo->TPC(tpcid).Length() << std::endl
                << "\t\t" << "active volume: " << asGeo->TPC(tpcid).ActiveWidth() << " x " << asGeo->TPC(tpcid).ActiveHeight() << " x " << asGeo->TPC(tpcid).ActiveLength() << std::endl
                << "\t\t" << "Planes:";
        for (unsigned int plane=0; plane < asWire->Nplanes(); plane++) {
            geo::PlaneID planeid{tpcid, plane};
            std::cout << "\t" << char('U'+plane) << " #Wires: " << asWire->Nwires(planeid);
        } // end loop over Planes
        std::cout << std::endl;

        geo::PlaneID collectionid{tpcid, geo::kW};
        unsigned int minchan = asWire->Nchannels();
        unsigned int maxchan = 0;
        for (unsigned int wire=0; wire < asWire->Nwires(collectionid); wire++) {
            geo::WireID wireid{collectionid, wire};
            if (asWire->PlaneWireToChannel(wireid) == raw::InvalidChannelID) continue;
            unsigned int chan = asWire->PlaneWireToChannel(wireid);
            if (chan < minchan) minchan = chan;
            if (chan > maxchan) maxchan = chan;
        } // end loop over Wires
        std::cout << "\t\t" << "W plane Channel range: " << minchan << " - " << maxchan << std::endl;
        std::cout << std::endl;
    } // end loop over TPCs

    std::cout << "----------------------TEST----------------------" << std::endl;

    int p = -1;
    for (unsigned int chan=0; chan < asWire->Nchannels(); chan++) {
        int plane = asWire->View(raw::ChannelID_t(chan));
        if (plane == p) continue;
        std::cout << "Channel#" << chan << " onward is from " << char('U'+plane) << " plane" << std::endl;
        p = plane;
    }

    std::cout << "----------------------TEST----------------------" << std::endl;

    geo::WireID wireid{geo::PlaneID{tpcid0, geo::kW}, 0};
    geo::WireID wireid1{geo::PlaneID{tpcid0, geo::kW}, 1};
    geo::WireGeo const wiregeo = asWire->Wire(wireid);
    geo::WireGeo const wiregeo1 = asWire->Wire(wireid1);

    std::cout << "channel: " << asWire->PlaneWireToChannel(wireid) << " " << wiregeo.GetStart() << " -> " << wiregeo.GetEnd() << std::endl;
    std::cout << "\t" << "pitch: " << geo::WireGeo::WirePitch(wiregeo, wiregeo1) << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << wiregeo.WireInfo() << std::endl;

    std::vector<double> pitch;
    for (unsigned int tpc=0; tpc<asGeo->NTPC(); tpc++) {
        geo::TPCID tpcid{cryoid, tpc};
        geo::PlaneID planeid{tpcid, geo::kW};

        for (unsigned int wire=0; wire<asWire->Nwires(planeid); wire++) {
            geo::WireID wireid{planeid, wire};
            geo::WireID wireid1{planeid, wire+1};
            if (asWire->PlaneWireToChannel(wireid) == raw::InvalidChannelID || asWire->PlaneWireToChannel(wireid1) == raw::InvalidChannelID) {
                std::cerr << "\033[91m" << "ERROR: Invalid channel for collection wire " << wire << " in TPC " << tpc << "\033[0m" << std::endl;
                continue;
            }

            geo::WireGeo const wiregeo = asWire->Wire(wireid);
	        std::cout << "    {" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetStart().Z() << "}," << std::endl;

            if (wire == asWire->Nwires(planeid)-1) continue;
            geo::WireGeo const wiregeo1 = asWire->Wire(wireid1);
            pitch.push_back(geo::WireGeo::WirePitch(wiregeo, wiregeo1));
        }
    }

    std::cout << "----------------------TEST----------------------" << std::endl;

    std::sort(pitch.begin(), pitch.end());
    std::map<double,unsigned int> mpp;
    for (double p : pitch) {
	mpp[p]++;
    }
    for (std::pair<double,unsigned int> pp : mpp) {
	std::cout << "pitch: " << pp.first << " count: " << pp.second << std::endl;
    }

    std::cout << "----------------------TEST----------------------" << std::endl;

    // geo::BoxBoundedGeo::Coord_t y0 = tpcgeo0.Min().Y();
    // geo::BoxBoundedGeo::Coord_t yN = tpcgeoN.Max().Y();
    // geo::BoxBoundedGeo::Coord_t z0 = tpcgeo0.Min().Z();
    // geo::BoxBoundedGeo::Coord_t zN = tpcgeoN.Max().Z();

    // const unsigned int n = 10;
    // const double dy = (yN - y0) / n;
    // const double dz = (zN - z0) / n;

    // for (unsigned int i=0; i<=n; i++) {
        // for (unsigned int j=0; j<=n; j++) {
            // double y = y0 + j*dy;
            // double z = z0 + i*dz;

            // geo::Point_t pt_low = {-50, y, z};
            // geo::Point_t pt_upp = {50, y, z};

            // raw::ChannelID_t ch_low = GetChannel(pt_low, geo::kW);
            // raw::ChannelID_t ch_upp = GetChannel(pt_upp, geo::kW);

            // std::cout << "(Y,Z)=(" << y << "," << z << ")" << "\r" << std::flush;
            // std::cout << std::string(4,'\t') << "LOW: ch: " << ch_low << "\r" << std::flush;
            // std::cout << std::string(7,'\t') << "UPP: ch: " << ch_upp << std::endl;
        // }
    // }

    std::cout << "\033[93m" << "End of Detchecks::beginJob ======================================================" << "\033[0m" << std::endl;
} // end beginJob


void ana::Detchecks::endJob()
{

} // end endJob

geo::WireID ana::Detchecks::GetWireID(geo::Point_t const& P, geo::View_t plane) {
    geo::TPCID tpcid = asGeo->FindTPCAtPosition(P);
    if (!tpcid.isValid) return geo::WireID();
    geo::PlaneGeo const& planegeo = asWire->Plane(tpcid,plane);
    geo::WireID wireid;
    try {
        wireid = planegeo.NearestWireID(P);
    } catch (geo::InvalidWireError const& e) {
        return e.suggestedWireID();
    }
    return wireid;
}
raw::ChannelID_t ana::Detchecks::GetChannel(geo::Point_t const& P, geo::View_t plane) {
    geo::WireID wireid = GetWireID(P, plane);
    if (!wireid.isValid) return raw::InvalidChannelID;
    return asWire->PlaneWireToChannel(wireid);
}

DEFINE_ART_MODULE(ana::Detchecks)
