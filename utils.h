#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include <TTree.h>
#include <TBranch.h>

#include <cstdio>
#include <algorithm>

#define LOG(x) (fLog ? printf("\tLOG " #x ": " "\033[1;9%dm" "%s" "\033[0m\n", x?2:1, x?"true":"false") : 0, x)
#define ASSERT(x)  if (!(fLog ? printf("\tAST " #x ": " "\033[1;9%dm" "%s" "\033[0m\n", x?2:1, x?"true":"false") : 0, x)) continue; 
#define DEBUG(x) if ((fLog ? printf("\tDBG " #x ": " "\033[1;9%dm" "%s" "\033[0m\n", x?1:2, x?"true":"false") : 0, x)) exit(1);


using PtrHit = art::Ptr<recob::Hit>;
using VecPtrHit = std::vector<art::Ptr<recob::Hit>>;

namespace ana {
    enum EnumDet { kPDVD, kPDHD };
    std::vector<unsigned> n_sec = {
        8, // PDVD
        2  // PDHD
    };
    std::vector<std::map<unsigned, int>> tpc2sec = {
        { // PDVD
            {0, 4}, {2, 4},
            {1, 5}, {3, 5},
            {4, 6}, {6, 6},
            {5, 7}, {7, 7},
            {8, 0}, {10, 0},
            {9, 1}, {11, 1},
            {12, 2}, {14, 2},
            {13, 3}, {15, 3}
        }, { // PDHD
            {4, -1}, {0, -1},
            {5, 0}, {1, 0},
            {6, 1}, {2, 1},
            {7, -1}, {3, -1}
        }
    };
    std::vector<std::map<int, std::pair<unsigned, unsigned>>> sec2tpc = {
        { // PDVD
            {0, {8, 10}},
            {1, {9, 11}},
            {2, {12, 14}},
            {3, {13, 15}},
            {4, {0, 2}},
            {5, {1, 3}},
            {6, {4, 6}},
            {7, {5, 7}}
        }, { // PDHD
            {0, {1, 5}},
            {1, {2, 6}}
        }
    };
    std::vector<std::map<unsigned, int>> tpc2side = {
        { // PDVD
            { 0, 0 }, { 1, 0 }, { 2, 0 }, { 3, 0 },
            { 4, 0 }, { 5, 0 }, { 6, 0 }, { 7, 0 },
            { 8, 1 }, { 9, 1 }, {10, 1 }, {11, 1 },
            {12, 1 }, {13, 1 }, {14, 1 }, {15, 1 }
        }, { // PDHD
            { 0, -1 }, { 1, 0 }, { 2, 1 }, { 3, -1 },
            { 4, -1 }, { 5, 0 }, { 6, 1 }, { 7, -1 }
        }
    };
    // std::vector<std::map<int, std::vector<int>>> side2sec = {
    //     { // PDVD
    //         { 0, {4, 5, 6, 7} }, { 1, {0, 1, 2, 3} }
    //     }, { // PDHD
    //         { 0, {0} }, { 1, {1} }
    //     }
    // };
    std::vector<std::map<int, int>> sec2side = {
        { // PDVD
            {0, 1}, {1, 1}, {2, 1}, {3, 1},
            {4, 0}, {5, 0}, {6, 0}, {7, 0},
        }, { // PDHD
            {0, 0}, {1, 1} 
        }
    };
    std::map<int, int> side2dirx = {
        { 0, -1 }, { 1, 1 }
    };

    // y = m*x + p
    struct LinearRegression {
        static constexpr unsigned nmin = 4;
        unsigned n=0;
        double mx=0, my=0, mx2=0, my2=0, mxy=0;
        double varx=0, vary=0, cov=0;
        double m=0, p=0, r2=0;
        double lp=0, corr=0;
        void add(double x, double y) {
            mx+=x; my+=y; mx2+=x*x; my2+=y*y; mxy+=x*y; n++;
        }
        void compute() {
            mx/=n; my/=n; mx2/=n; my2/=n; mxy/=n;
            varx=mx2-mx*mx;
            vary=my2-my*my;
            cov=mxy-mx*my;
            m= n<nmin ? 0 : cov/varx;
            p= my - m*mx;
            r2= n<nmin ? 0 : cov*cov / (varx*vary);
            lp= .5*(varx+vary + sqrt(pow(varx+vary, 2) + 4*cov*cov));
            corr= 1 - 4 * (varx*vary - cov*cov) / pow(varx+vary, 2);
        }
        double projection(double x, double y) const {
            return (x + m*(y-p)) / (1 + m*m);
        }
        double theta(int dirx) {
            return atan2(dirx * abs(cov), dirx* abs(lp-vary));
        }
        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sRegM", pre), &m);
            t->Branch(Form("%sRegP", pre), &p);
            t->Branch(Form("%sRegR2", pre), &r2);
        }
    };

    template<typename T>
    struct bounds {
        T min, max;
        bounds() : min(std::numeric_limits<T>::max()), max(std::numeric_limits<T>::lowest()) {}
        bounds(T m, T M) : min(m), max(M) {}
        bool isInside(T x, float r=0) const { 
            return min+r <= x && x <= max-r;
        }

        friend std::ostream& operator<<(std::ostream& os, const bounds& b) {
            return os << "[" << b.min << ", " << b.max << "]";
        }
    };
    template<typename T>
    struct bounds3D {
        bounds<T> x, y, z;
        bounds3D() : x(), y(), z() {}
        bounds3D(geo::Point_t const& min, geo::Point_t const& max) :
            x{T(min.x()), T(max.x())}, y{T(min.y()), T(max.y())}, z{T(min.z()), T(max.z())} {}
        // bounds3D(geo::BoxBoundedGeo const& bb) :
        //     x{bb.MinX(), bb.MaxX()}, y{bb.MinY(), bb.MaxY()}, z{bb.MinZ(), bb.MaxZ()} {}
        bool isInside(geo::Point_t const& p, float r=0) const {
            return x.isInside(p.x(), r) && y.isInside(p.y(), r) && z.isInside(p.z(), r);
        }
        bool isInside(TVector3 const& p, float r=0) const {
            return x.isInside(p.x(), r) && y.isInside(p.y(), r) && z.isInside(p.z(), r);
        }
        bool isInsideYZ(geo::Point_t const& p, float r=0) const {
            return y.isInside(p.y(), r) && z.isInside(p.z(), r);
        }
        geo::Point_t min() const {
            return geo::Point_t{x.min, y.min, z.min};
        }
        geo::Point_t max() const {
            return geo::Point_t{x.max, y.max, z.max};
        }

        friend std::ostream& operator<<(std::ostream& os, const bounds3D& b) {
            return os << b.min() << " -> " << b.max();
        }
    };

    struct axis {
        double ay, az;
        double space(geo::WireGeo wiregeo) {
            return ay * wiregeo.GetCenter().Y() + az * wiregeo.GetCenter().Z();
        }
    };

    // class Hit: public recob::Hit {};

    struct Hit {
        unsigned tpc;
        int section;
        float space;
        unsigned channel;
        float tick;
        float adc;
        Hit() : tpc(0), section(0), space(0), channel(0), tick(0), adc(0) {}
        Hit(unsigned T, int S, float s, unsigned c, float t, float a) :
            tpc(T), section(S), space(s), channel(c), tick(t), adc(a) {}

        // unsigned slice() { return 2*(tpc/4) + tpc%2; }
        friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
            return os << "tpc:" << hit.tpc << " space:" << hit.space << " ch:" << hit.channel << " tick:" << hit.tick << " ADC:" << hit.adc;
        }
        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sHitTPC", pre), &tpc);
            t->Branch(Form("%sHitSection", pre), &section);
            t->Branch(Form("%sHitSpace", pre), &space);
            t->Branch(Form("%sHitChannel", pre), &channel);
            t->Branch(Form("%sHitTick", pre), &tick);
            t->Branch(Form("%sHitADC", pre), &adc);
        }
    };
    struct Hits {
        unsigned N;
        std::vector<unsigned> tpc;
        std::vector<int> section;
        std::vector<float> space;
        std::vector<unsigned> channel;
        std::vector<float> tick;
        std::vector<float> adc;
        Hits() : N(0), tpc(), section(), space(), channel(), tick(), adc() {}
        void push_back(Hit const& hit) {
            N++;
            tpc.push_back(hit.tpc);
            section.push_back(hit.section);
            space.push_back(hit.space);
            channel.push_back(hit.channel);
            tick.push_back(hit.tick);
            adc.push_back(hit.adc);
        }
        void clear() {
            N = 0;
            tpc.clear();
            section.clear();
            space.clear();
            channel.clear();
            tick.clear();
            adc.clear();
        }
        unsigned size() const { return N; }
        bool empty() const { return !N; }
        float energy() const {
            float e = 0;
            for (unsigned i=0; i<N; i++) e += adc[i];
            return e;
        }

        void SetBranches(TTree *t, const char* pre="") {
            t->Branch(Form("%sNHit", pre), &N);
            t->Branch(Form("%sHitTPC", pre), &tpc);
            t->Branch(Form("%sHitSection", pre), &section);
            t->Branch(Form("%sHitSpace", pre), &space);
            t->Branch(Form("%sHitChannel", pre), &channel);
            t->Branch(Form("%sHitTick", pre), &tick);
            t->Branch(Form("%sHitADC", pre), &adc);
        }

        Hit at(unsigned i) const { return Hit{tpc[i], section[i], space[i], channel[i], tick[i], adc[i]}; }
        struct iterator {
            const Hits* hits;
            unsigned i;
            iterator(Hits const* h, unsigned i=0) : hits(h), i(i) {}
            iterator& operator++() { i++; return *this; }
            bool operator!=(iterator const& h) const { return i != h.i; }
            Hit operator*() const { return hits->at(i); }
        };
        iterator begin() const { return iterator(this); }
        iterator end() const { return iterator(this, N); }
        
        bool find(recob::Hit hit) {
            auto it_chan = std::find(channel.begin(), channel.end(), hit.Channel());
            auto it_tick = std::find(tick.begin(), tick.end(), hit.PeakTime());
            unsigned i = std::distance(channel.begin(), it_chan);
            if (i == std::distance(tick.begin(), it_tick)) return i;
            return N;
        }
    };
    struct Point {
        float x, y, z;
        Point() : x(0), y(0), z(0) {}
        Point(float x, float y, float z) : x(x), y(y), z(z) {}
        Point(double x, double y, double z) : x(x), y(y), z(z) {}
        Point(geo::Point_t const& p) : x(p.x()), y(p.y()), z(p.z()) {}
        Point(TVector3 const& v) : x(v.x()), y(v.y()), z(v.z()) {}
        float r2() const { return x*x + y*y + z*z; }

        Point operator+(Point const& p) const { return Point{x+p.x, y+p.y, z+p.z}; }
        Point operator+(geo::Point_t const& p) const { return Point{x+(float)p.x(), y+(float)p.y(), z+(float)p.z()}; }
        Point operator+=(Point const& p) { x += p.x; y += p.y; z += p.z; return *this; }
        Point operator-(Point const& p) const { return Point{x-p.x, y-p.y, z-p.z}; }
        Point operator-(geo::Point_t const& p) const { return Point{x-(float)p.x(), y-(float)p.y(), z-(float)p.z()}; }
        Point operator-=(Point const& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
        Point operator*(float f) const { return Point{f*x, f*y, f*z}; }

        operator bool() const { return !(x==0 && y==0 && z==0); }

        friend std::ostream& operator<<(std::ostream& os, const Point& p) {
            return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
        }

        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sPointX", pre), &x);
            t->Branch(Form("%sPointY", pre), &y);
            t->Branch(Form("%sPointZ", pre), &z);
        }
    };
    struct Points {
        unsigned N;
        std::vector<float> x, y, z;
        Points() : N(0), x(), y(), z() {}
        void push_back(Point const& p) {
            N++;
            x.push_back(p.x);
            y.push_back(p.y);
            z.push_back(p.z);
        }
        void push_back(geo::Point_t const& p) { push_back(Point{p}); }
        void clear() {
            N = 0;
            x.clear();
            y.clear();
            z.clear();
        }
        unsigned size() const { return N; }
        bool empty() const { return !N; }
        Point barycenter() const {
            Point b;
            for (unsigned i=0; i<N; i++) b += at(i);
            return b * (1.F / N);
        }


        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sNPoint", pre), &N);
            t->Branch(Form("%sPointX", pre), &x);
            t->Branch(Form("%sPointY", pre), &y);
            t->Branch(Form("%sPointZ", pre), &z);
        }

        Point at(unsigned i) const { return Point{x[i], y[i], z[i]}; }
        struct iterator {
            const Points* points;
            unsigned i;
            iterator(Points const* p, unsigned i=0) : points(p), i(i) {}
            iterator& operator++() { i++; return *this; }
            bool operator!=(iterator const& p) { return i != p.i; }
            Point operator*() { return points->at(i); }
        };
        iterator begin() const { return iterator(this); }
        iterator end() const { return iterator(this, N); }
    };



    // reco to truth functions
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    simb::MCParticle const* trk2mcp(art::Ptr<recob::Track> const& p_trk, detinfo::DetectorClocksData const& clockData, art::FindManyP<recob::Hit> const& fmp_trk2hit) {
        std::unordered_map<int, float> map_tid_ene;
        for (art::Ptr<recob::Hit> const& p_hit : fmp_trk2hit.at(p_trk.key()))
            for (sim::TrackIDE ide : bt_serv->HitToTrackIDEs(clockData, p_hit))
                map_tid_ene[ide.trackID] += ide.energy;

        float max_ene = -1;
        int tid_max = 0;
        for (std::pair<int, float> p : map_tid_ene)
            if (p.second > max_ene)
                max_ene = p.second, tid_max = p.first;
        return max_ene == -1 ? nullptr : pi_serv->TrackIdToParticle_P(tid_max);
    }

    art::Ptr<recob::Track> mcp2trk(simb::MCParticle const* mcp, std::vector<art::Ptr<recob::Track>> const& vp_trk, detinfo::DetectorClocksData const& clockData, art::FindManyP<recob::Hit> const& fmp_trk2hit) {
        std::unordered_map<art::Ptr<recob::Track>, unsigned> map_trk_nhit;
        for (art::Ptr<recob::Track> p_trk : vp_trk)
            map_trk_nhit[p_trk] += bt_serv->TrackIdToHits_Ps(clockData, mcp->TrackId(), fmp_trk2hit.at(p_trk.key())).size();

        unsigned max = 0;
        art::Ptr<recob::Track> p_trk_from_mcp;
        for (std::pair<art::Ptr<recob::Track>, unsigned> p : map_trk_nhit)
            if (p.second > max)
                max = p.second, p_trk_from_mcp = p.first;
        return p_trk_from_mcp;
    }

    std::vector<art::Ptr<recob::Track>> mcp2trks(simb::MCParticle const* mcp, std::vector<art::Ptr<recob::Track>> const& vp_trk, detinfo::DetectorClocksData const& clockData, art::FindManyP<recob::Hit> const& fmp_trk2hit) {
        std::vector<art::Ptr<recob::Track>> vp_trk_from_mcp;
        for (art::Ptr<recob::Track> p_trk : vp_trk) {
            simb::MCParticle const* mcp_from_trk = trk2mcp(p_trk, clockData, fmp_trk2hit);
            if (mcp_from_trk && mcp_from_trk->TrackId() == mcp->TrackId())
                vp_trk_from_mcp.push_back(p_trk);
        }
        return vp_trk_from_mcp;
    }

    std::vector<art::Ptr<recob::Hit>> mcp2hits(simb::MCParticle const* mcp, std::vector<art::Ptr<recob::Hit>> const& vp_hit, detinfo::DetectorClocksData const& clockData, bool use_eve) {
        std::vector<art::Ptr<recob::Hit>> vp_hit_from_mcp;
        for (art::Ptr<recob::Hit> p_hit : vp_hit)
            for (sim::TrackIDE ide : (use_eve ? bt_serv->HitToEveTrackIDEs(clockData, p_hit) : bt_serv->HitToTrackIDEs(clockData, p_hit)))
                if (ide.trackID == mcp->TrackId())
                    vp_hit_from_mcp.push_back(p_hit);
        return vp_hit_from_mcp;
    }
































    struct SortedHits {
        VecPtrHit vph;
        unsigned i_cathode_crossing;
        std::vector<unsigned> vi_section_crossing;
        std::vector<LinearRegression> regs;
        SortedHits() : vph(), i_cathode_crossing(0), vi_section_crossing(), regs(2) {}
        bool isCathodeCrossing() const {
            return i_cathode_crossing != vph.size()
                && i_cathode_crossing != 0;
        }
        operator bool() const {
            return !vph.empty();
        }
        PtrHit lastHit(int dirz) const {
            int dirx = (regs[1].mx - regs[0].mx) * dirz > 0 ? 1 : -1;
            if (isCathodeCrossing()) {
                if (dirx > 0) {
                    if (dirz > 0)
                        return vph.back();
                    else 
                        return vph[i_cathode_crossing];
                } else {
                    if (dirz > 0)
                        return vph[i_cathode_crossing-1];
                    else
                        return vph.front();
                }
            } else {
                if (dirz > 0)
                    return vph.back();
                else
                    return vph.front();
            }
        }
    };

    class MichelAnalyzer {
    public:
        // Utilities
        art::ServiceHandle<art::TFileService> asFile;
        art::ServiceHandle<cheat::ParticleInventoryService> asPartInv;
        art::ServiceHandle<cheat::BackTrackerService> asBackTrack;

        const geo::GeometryCore* asGeo;
        const geo::WireReadoutGeom* asWire;
        const detinfo::DetectorPropertiesService* asDetProp;
        const detinfo::DetectorClocksService* asDetClocks;

        int geoDet;
        enum EnumDet { kPDVD, kPDHD };

        float fTick2cm;

        axis GetAxis(geo::PlaneID) const;
        double GetSpace(geo::WireID) const;
        Hit GetHit(art::Ptr<recob::Hit> const&) const;
        SortedHits GetSortedHits(
            VecPtrHit const& vph_unsorted,
            geo::View_t view = geo::kW
        ) const;
    };

    axis MichelAnalyzer::GetAxis(geo::PlaneID pid) const {
        geo::WireGeo w0 = asWire->Wire(geo::WireID{pid, 0});
        geo::WireGeo w1 = asWire->Wire(geo::WireID{pid, 1});
        int dy = w1.GetCenter().Y() > w0.GetCenter().Y() ? 1 : -1;
        int dz = w1.GetCenter().Z() > w0.GetCenter().Z() ? 1 : -1;
        return { dy * w0.CosThetaZ(), dz * w0.SinThetaZ() };
    }
    double MichelAnalyzer::GetSpace(geo::WireID wid) const {
        return GetAxis(wid).space(asWire->Wire(wid));
    }
    Hit MichelAnalyzer::GetHit(art::Ptr<recob::Hit> const& ph) const {
        geo::WireID wid = ph->WireID();
        return {
            wid.TPC,
            tpc2sec[geoDet][wid.TPC],
            float(GetSpace(wid)),
            ph->Channel(),
            ph->PeakTime(),
            ph->Integral()
        };
    }
    SortedHits MichelAnalyzer::GetSortedHits(
        VecPtrHit const& vph_unsorted,
        geo::View_t view
    ) const {
        SortedHits sh;
        for (PtrHit const& ph : vph_unsorted) {
            if (ph->View() != view) continue;
            int side = tpc2side[geoDet][ph->WireID().TPC];
            if (side == -1) continue;
            double z = GetSpace(ph->WireID());
            double t = ph->PeakTime() * fTick2cm;
            sh.regs[side].add(z, t);
            if (side) // highX
                sh.vph.push_back(ph);
            else {  // lowX
                sh.vph.insert(sh.vph.begin(), ph);
                sh.i_cathode_crossing++; 
            }
        }
        std::vector<unsigned> ns(2);
        ns[0] = sh.i_cathode_crossing;
        ns[1] = sh.vph.size() - sh.i_cathode_crossing;
        if ((
            0 < ns[0] && ns[0] <= LinearRegression::nmin
        ) || (
            0 < ns[1] && ns[1] <= LinearRegression::nmin
        )) {
            return SortedHits{};
        }
        sh.regs[0].compute();
        sh.regs[1].compute();
        auto comp = [&](LinearRegression const& reg, PtrHit const& ph1, PtrHit const& ph2) {
            double s1 = reg.projection(GetSpace(ph1->WireID()), ph1->PeakTime() * fTick2cm);
            double s2 = reg.projection(GetSpace(ph2->WireID()), ph2->PeakTime() * fTick2cm);
            return s1 < s2; // z increasing
        };
        std::sort(
            sh.vph.begin(), sh.vph.begin()+sh.i_cathode_crossing,
            [&](PtrHit const& ph1, PtrHit const& ph2) {
                return comp(sh.regs[0], ph1, ph2);
            }
        );
        std::sort(
            sh.vph.begin()+sh.i_cathode_crossing, sh.vph.end(),
            [&](PtrHit const& ph1, PtrHit const& ph2) {
                return comp(sh.regs[1], ph1, ph2);
            }
        );
        int prev_sec = tpc2sec[geoDet][sh.vph.front()->WireID().TPC];
        for (unsigned i=1; i<sh.vph.size(); i++) {
            int sec = tpc2sec[geoDet][(sh.vph[i])->WireID().TPC];
            if (sec != prev_sec) {
                sh.vi_section_crossing.push_back(i-1);
                sh.vi_section_crossing.push_back(i);
                prev_sec = sec;
            }
        }
        return sh;
    }
}
