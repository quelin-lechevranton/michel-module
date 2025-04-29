////////////////////////////////////////////////////////////////////////
// Class:       Detchecks
// Plugin Type: analyzer
// File:        Detchecks_module.cc
//
// Generated at Wed Nov  6 10:31:28 2024 by Jeremy Quelin Lechevranton
////////////////////////////////////////////////////////////////////////

#include "pdvd_utils.h"

namespace ana {
  class Detchecks;
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
    const detinfo::LArPropertiesService* asLarProp;

    // Input Parameters
    bool pWireGeo;
    bool pWireEnds;
    bool pTPCBounds;
    bool pTPCChannels;
    bool pChannelPitch;
    bool pCollectionZ;
    bool pViewUCoord;
    bool pViewVCoord;

    // Conversion factors
    float fADCtoE = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    float fChannelPitch = 0.5; // cm/channel

    // geo::WireID GetWireID(geo::Point_t const& P, geo::View_t plane);
    // raw::ChannelID_t GetChannel(geo::Point_t const& P, geo::View_t plane);
};


ana::Detchecks::Detchecks(fhicl::ParameterSet const& p)
    : EDAnalyzer{p},
    pWireGeo(p.get<bool>("WireGeo", false)),
    pWireEnds(p.get<bool>("WireEnds", false)),
    pTPCBounds(p.get<bool>("TPCBounds", false)),
    pTPCChannels(p.get<bool>("TPCChannels", false)),
    pChannelPitch(p.get<bool>("ChannelPitch", false)),
    pCollectionZ(p.get<bool>("CollectionZ", false)),
    pViewUCoord(p.get<bool>("ViewUCoord", false)),
    pViewVCoord(p.get<bool>("ViewVCoord", false))


    // fMichelTimeRadius(p.get<float>("MichelTimeRadius")), //in µs
    // fMichelSpaceRadius(p.get<float>("MichelSpaceRadius")), //in cm
{
    // Basic Utilities
    asGeo = &*art::ServiceHandle<geo::Geometry>();
    asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();
    asLarProp = &*art::ServiceHandle<detinfo::LArPropertiesService>();
}

void ana::Detchecks::analyze(art::Event const& e)
{

} // end analyze

void ana::Detchecks::beginJob()
{
    detinfo::DetectorClocksData const clockData = asDetClocks->DataForJob();
    detinfo::DetectorPropertiesData const detProp = asDetProp->DataForJob(clockData);
    detinfo::LArProperties const * larProp = asLarProp->provider();

    geo::CryostatID cryoid{0};
    std::cout << "\033[93m" << "Detchecks::beginJob: ============================================================" << "\033[0m" << std::endl
        << "--- Geometry ---" << std::endl
        << "Detector name: " << asGeo->DetectorName() << std::endl
        << "Number of cryostats: " << asGeo->Ncryostats() << std::endl
        << "Number of TPCs: " << asGeo->NTPC() << std::endl
        << "--- WireReadout ---" << std::endl
        << "Number of channels: " << asWire->Nchannels() << std::endl
        << "--- DetectorProperties ---" << std::endl
        << "Number of ticks: " << detProp.ReadOutWindowSize() << std::endl
        << "Drift velocity: " << detProp.DriftVelocity() << " cm/µs" << std::endl
        << "--- DetectorClocks ---" << std::endl
        << "Sampling rate: " << detinfo::sampling_rate(clockData) << " ns/tick" << std::endl
        << "Time to tick: " << clockData.Time2Tick(1.) << " tick/?s" << std::endl
        // << "Tick to X" << detProp.ConvertTicksToX(??) << " cm/tick" << std::endl
        <<  "--- LArProperties ---" << std::endl
        << "Atomic mass of the liquid " << larProp->AtomicMass() << " g/mol" << std::endl;

    // std::cout << "Cryostat:" << std::endl
    //     << "\t" << "coordinates: " << asGeo->Cryostat(cryoid).Min() 
    //                                 << " - " << asGeo->Cryostat(cryoid).Max() 
    //                                 << std::endl
    //     << "\t" << "volume: " << asGeo->Cryostat(cryoid).Width() << " x " << asGeo->Cryostat(cryoid).Height() << " x " << asGeo->Cryostat(cryoid).Length() << std::endl;

    geo::TPCID tpcid0{cryoid, 0};
    geo::TPCID tpcidN{cryoid, asGeo->NTPC(cryoid)-1};

    // geo::TPCGeo tpcgeo0 = asGeo->TPC(tpcid0);
    // geo::TPCGeo tpcgeoN = asGeo->TPC(tpcidN);
    // std::cout << "\t" << "active coordinates: " << tpcgeo0.Min() << " - " << tpcgeoN.Max() << std::endl;
    // geo::Vector_t diag = tpcgeoN.Max() - tpcgeo0.Min();
    // std::cout << "\t" << "active volume: "  << diag.X() << " x " << diag.Y() << " x " << diag.Z() << std::endl;

    // std::cout << "TPCs:" << std::endl; 
    // for (unsigned int tpc=0; tpc < asGeo->NTPC(); tpc++) {
    //     geo::TPCID tpcid{cryoid, tpc};
    //     std::cout << "\t" << "TPC#" << tpc << "\tcoordinates: " << asGeo->TPC(tpcid).Min() << " - " << asGeo->TPC(tpcid).Max() << std::endl
    //             << "\t\t" << "volume: " << asGeo->TPC(tpcid).Width() << " x " << asGeo->TPC(tpcid).Height() << " x " << asGeo->TPC(tpcid).Length() << std::endl
    //             << "\t\t" << "active volume: " << asGeo->TPC(tpcid).ActiveWidth() << " x " << asGeo->TPC(tpcid).ActiveHeight() << " x " << asGeo->TPC(tpcid).ActiveLength() << std::endl
    //             << "\t\t" << "Planes:";
    //     for (unsigned int plane=0; plane < asWire->Nplanes(); plane++) {
    //         geo::PlaneID planeid{tpcid, plane};
    //         std::cout << "\t" << char('U'+plane) << " #Wires: " << asWire->Nwires(planeid);
    //     } // end loop over Planes
    //     std::cout << std::endl;

    //     geo::PlaneID collectionid{tpcid, geo::kW};
    //     unsigned int minchan = asWire->Nchannels();
    //     unsigned int maxchan = 0;
    //     for (unsigned int wire=0; wire < asWire->Nwires(collectionid); wire++) {
    //         geo::WireID wireid{collectionid, wire};
    //         if (asWire->PlaneWireToChannel(wireid) == raw::InvalidChannelID) continue;
    //         unsigned int chan = asWire->PlaneWireToChannel(wireid);
    //         if (chan < minchan) minchan = chan;
    //         if (chan > maxchan) maxchan = chan;
    //     } // end loop over Wires
    //     std::cout << "\t\t" << "W plane Channel range: " << minchan << " - " << maxchan << std::endl;
    //     std::cout << std::endl;
    // } // end loop over TPCs

    // std::cout << "----------------------TEST----------------------" << std::endl;

    // int p = -1;
    // for (unsigned int chan=0; chan < asWire->Nchannels(); chan++) {
    //     int plane = asWire->View(raw::ChannelID_t(chan));
    //     if (plane == p) continue;
    //     std::cout << "Channel#" << chan << " onward is from " << char('U'+plane) << " plane" << std::endl;
    //     p = plane;
    // }

    // std::cout << "----------------------TEST----------------------" << std::endl;

    std::cout << "\033[93m" "Cryostat Bounds: " "\033[0m" << std::endl;
    for (unsigned c=0; c<asGeo->Ncryostats(); c++) {
        geo::CryostatGeo cryogeo = asGeo->Cryostat(geo::CryostatID{c});
        std::cout << cryogeo.Min() << " -> " << cryogeo.Max() << std::endl;

        // std::cout << "\t{" << c << ", {"
        //     << cryogeo.MinX() << ", " << cryogeo.MaxX() << ", "
        //     << cryogeo.MinY() << ", " << cryogeo.MaxY() << ", "
        //     << cryogeo.MinZ() << ", " << cryogeo.MaxZ() << " }"
        //     << "}," << std::endl;
        // }
    }

    if (pWireGeo) {
        std::cout << "\033[93m" "Wire Geometry:" "\033[0m" << std::endl;
        for (unsigned t=0; t<asGeo->NTPC(); t++) {
            geo::TPCID tpcid{cryoid, t};
            geo::TPCGeo tpcgeo = asGeo->TPC(tpcid);

            std::cout << "TPC#" << t << ": " << tpcgeo.Min() << " -> " << tpcgeo.Max() << std::endl;

            geo::PlaneID Uid{tpcid, geo::kU};
            std::cout << "  U plane" << std::endl;
            std::cout << "\t" << "number of wires: " << asWire->Nwires(Uid) << std::endl;
            geo::WireGeo Uw0 = asWire->Wire(geo::WireID{Uid, 0});
            geo::WireGeo Uw1 = asWire->Wire(geo::WireID{Uid, 1});
            std::cout << "\t" << "wire #0: " << Uw0.GetStart() << " -> " << Uw0.GetEnd() << std::endl;
            std::cout << "\t" << "wire #1: " << Uw1.GetStart() << " -> " << Uw1.GetEnd() << std::endl;
            std::cout << "\t" << "wires Z->Y angle: " << Uw0.ThetaZ(true) << "º" << std::endl;
            std::cout << "\t" << "channel pitch: " << geo::WireGeo::WirePitch(Uw0, Uw1) << " cm" << std::endl;
            
            geo::PlaneID Vid{tpcid, geo::kV};
            std::cout << "  V plane" << std::endl;
            std::cout << "\t" << "number of wires: " << asWire->Nwires(Vid) << std::endl;
            geo::WireGeo Vw0 = asWire->Wire(geo::WireID{Vid, 0});
            geo::WireGeo Vw1 = asWire->Wire(geo::WireID{Vid, 1});
            std::cout << "\t" << "wire #0: " << Vw0.GetStart() << " -> " << Vw0.GetEnd() << std::endl;
            std::cout << "\t" << "wire #1: " << Vw1.GetStart() << " -> " << Vw1.GetEnd() << std::endl;
            std::cout << "\t" << "wires Z->Y angle: " << Vw0.ThetaZ(true) << "º" << std::endl;
            std::cout << "\t" << "channel pitch: " << geo::WireGeo::WirePitch(Vw0, Vw1) << " cm" << std::endl;

            geo::PlaneID Wid{tpcid, geo::kW};
            std::cout << "  W plane" << std::endl;
            std::cout << "\t" << "number of wires: " << asWire->Nwires(Wid) << std::endl;
            geo::WireGeo Ww0 = asWire->Wire(geo::WireID{Wid, 0});
            geo::WireGeo Ww1 = asWire->Wire(geo::WireID{Wid, 1});
            std::cout << "\t" << "wire #0: " << Ww0.GetStart() << " -> " << Ww0.GetEnd() << std::endl;
            std::cout << "\t" << "wire #1: " << Ww1.GetStart() << " -> " << Ww1.GetEnd() << std::endl;
            std::cout << "\t" << "wires Z->Y angle: " << Ww0.ThetaZ(true) << "º" << std::endl;
            std::cout << "\t" << "channel pitch: " << geo::WireGeo::WirePitch(Ww0, Ww1) << " cm" << std::endl;
        }
    }

    if (pWireEnds) {
        std::cout << "\033[93m" "Wire Ends:" "\033[0m" << std::endl;
        for (unsigned int tpc=0; tpc<asGeo->NTPC(); tpc++) {
            geo::TPCID tpcid{cryoid, tpc};
            geo::TPCGeo tpcgeo = asGeo->TPC(tpcid);
            geo::PlaneID planeidU{tpcid, geo::kU};
            geo::PlaneID planeidV{tpcid, geo::kV};
            geo::PlaneID planeidW{tpcid, geo::kW};

            std::cout << "TPC#" << tpc << ": " << tpcgeo.Min() << " -> " << tpcgeo.Max() << std::endl;
            std::cout << "  U plane, " << asWire->Nwires(planeidU) << " wires:" << std::endl;
            for (unsigned int wire=0; wire<asWire->Nwires(planeidU); wire++) {
                geo::WireID wireid{planeidU, wire};
                geo::WireGeo const wiregeo = asWire->Wire(wireid);
                std::cout << "\t" << wire << ": " << wiregeo.GetStart() << " -> " << wiregeo.GetEnd() << " ch: " << asWire->PlaneWireToChannel(wireid) << std::endl;
                // std::cout << Form("<path data-tpc=\"%u\" data-plane=\"U\" data-wire=\"%u\" data-channel=\"%u\" d=\"M %f %f L %f %f\" stroke=\"green\" stroke-width=\"0.1\" />", tpc, wire, asWire->PlaneWireToChannel(wireid), wiregeo.GetStart().Y(), wiregeo.GetStart().Z(), wiregeo.GetEnd().Y(), wiregeo.GetEnd().Z()) << std::endl;
            }
            std::cout << "  V plane, " << asWire->Nwires(planeidV) << " wires:" << std::endl;
            for (unsigned int wire=0; wire<asWire->Nwires(planeidV); wire++) {
                geo::WireID wireid{planeidV, wire};
                geo::WireGeo const wiregeo = asWire->Wire(wireid);
                std::cout << "\t" << wire << ": " << wiregeo.GetStart() << " -> " << wiregeo.GetEnd() << " ch: " << asWire->PlaneWireToChannel(wireid) << std::endl;
                // std::cout << Form("<path data-tpc=\"%u\" data-plane=\"V\" data-wire=\"%u\" data-channel=\"%u\" d=\"M %f %f L %f %f\" stroke=\"blue\" stroke-width=\"0.1\" />", tpc, wire, asWire->PlaneWireToChannel(wireid), wiregeo.GetStart().Y(), wiregeo.GetStart().Z(), wiregeo.GetEnd().Y(), wiregeo.GetEnd().Z()) << std::endl;
            }
            std::cout << "  W plane, " << asWire->Nwires(planeidW) << " wires:" << std::endl;
            for (unsigned int wire=0; wire<asWire->Nwires(planeidW); wire++) {
                geo::WireID wireid{planeidW, wire};
                geo::WireGeo const wiregeo = asWire->Wire(wireid);
                std::cout << "\t" << wire << ": " << wiregeo.GetStart() << " -> " << wiregeo.GetEnd() << " ch: " << asWire->PlaneWireToChannel(wireid) << std::endl;
                // std::cout << Form("<path data-tpc=\"%u\" data-plane=\"W\" data-wire=\"%u\" data-channel=\"%u\" d=\"M %f %f L %f %f\" stroke=\"red\" stroke-width=\"0.1\" />", tpc, wire, asWire->PlaneWireToChannel(wireid), wiregeo.GetStart().Y(), wiregeo.GetStart().Z(), wiregeo.GetEnd().Y(), wiregeo.GetEnd().Z()) << std::endl;
            }
        }
    }


    if (pTPCBounds) {
        std::cout << "\033[93m" "TPC Bounds: {TPC, {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax} }" "\033[0m" << std::endl;
        for (unsigned tpc=0; tpc<asGeo->NTPC(); tpc++) {
            geo::TPCID tpcid{cryoid, tpc};
            geo::TPCGeo tpcgeo = asGeo->TPC(tpcid);

            std::cout << "\t{" << tpc << ", {"
                << tpcgeo.MinX() << ", " << tpcgeo.MaxX() << ", "
                << tpcgeo.MinY() << ", " << tpcgeo.MaxY() << ", "
                << tpcgeo.MinZ() << ", " << tpcgeo.MaxZ() << " }"
                << "}," << std::endl;
        }
    }

    if (pTPCChannels) {
        std::cout << "\033[93m" "TPC Channels: {TPC, {{view, {ch min, ch max}}, ...}}" "\033[0m" << std::endl;
        for (unsigned tpc=0; tpc<asGeo->NTPC(); tpc++) {
            geo::TPCID tpcid{cryoid, tpc};
            std::cout << "\t{" << tpc << ", {\r" << std::flush;
            for (unsigned view=0; view<3; view++) {
                geo::PlaneID planeid{tpcid, view};

                unsigned min=13000, max=0;
                for (unsigned wire=0; wire<asWire->Nwires(planeid); wire++) {
                    geo::WireID wireid{planeid, wire};
                    raw::ChannelID_t ch = asWire->PlaneWireToChannel(wireid);
                    if (ch == raw::InvalidChannelID) continue;
                    min = ch < min ? ch : min; 
                    max = ch > max ? ch : max;
                }
                std::cout << "\t\t{" << view << ", {" << min << ", " << max << "} }" << std::endl;
            }
            std::cout << "\t}}," << std::endl;
        }
    }

    if (pCollectionZ) {
        std::cout << "\033[93m" "Collection Z: {channel, Z}" "\033[0m" << std::endl;
        for (unsigned tpc=0; tpc<asGeo->NTPC(); tpc++) {
            std::cout << "TPC#" << tpc << std::endl;
            geo::TPCID tpcid{cryoid, tpc};
            geo::PlaneID planeid{tpcid, geo::kW};
            for (unsigned wire=0; wire<asWire->Nwires(planeid); wire++) {
                geo::WireID wireid{planeid, wire};
                geo::WireGeo const wiregeo = asWire->Wire(wireid);
                std::cout << "\t{" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetCenter().Z() << "}," << std::endl;
            }
        }
    }

    if (pViewUCoord) {
        std::cout << "\033[93m" "View U Coordinates: {channel, sqrt(3)*Y-Z}" "\033[0m" << std::endl;
        for (unsigned tpc=0; tpc<asGeo->NTPC(); tpc++) {
            std::cout << "TPC#" << tpc << std::endl;
            geo::TPCID tpcid{cryoid, tpc};
            geo::PlaneID planeid{tpcid, geo::kU};
            for (unsigned wire=0; wire<asWire->Nwires(planeid); wire++) {
                geo::WireID wireid{planeid, wire};
                geo::WireGeo const wiregeo = asWire->Wire(wireid);
                std::cout << "\t{" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetCenter().Y() * sqrt(3) - wiregeo.GetCenter().Z() << "}," << std::endl;
                std::cout << "\t{" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetStart().Y() * sqrt(3) - wiregeo.GetStart().Z() << "}," << std::endl;
                std::cout << "\t{" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetEnd().Y() * sqrt(3) - wiregeo.GetEnd().Z() << "}," << std::endl;
            }
        }
    }

    if (pViewVCoord) {
        std::cout << "\033[93m" "View V Coordinates: {channel, sqrt(3)*Y+Z}" "\033[0m" << std::endl;
        for (unsigned tpc=0; tpc<asGeo->NTPC(); tpc++) {
            std::cout << "TPC#" << tpc << std::endl;
            geo::TPCID tpcid{cryoid, tpc};
            geo::PlaneID planeid{tpcid, geo::kV};
            for (unsigned wire=0; wire<asWire->Nwires(planeid); wire++) {
                geo::WireID wireid{planeid, wire};
                geo::WireGeo const wiregeo = asWire->Wire(wireid);
                std::cout << "\t{" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetCenter().Y() * sqrt(3) + wiregeo.GetCenter().Z() << "}," << std::endl;
                std::cout << "\t{" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetStart().Y() * sqrt(3) + wiregeo.GetStart().Z() << "}," << std::endl;
                std::cout << "\t{" << asWire->PlaneWireToChannel(wireid) << ", " << wiregeo.GetEnd().Y() * sqrt(3) + wiregeo.GetEnd().Z() << "}," << std::endl;
            }
        }
    }

    std::cout << "\033[93m" << "End of Detchecks::beginJob ======================================================" << "\033[0m" << std::endl;
} // end beginJob


void ana::Detchecks::endJob()
{

} // end endJob

DEFINE_ART_MODULE(ana::Detchecks)
