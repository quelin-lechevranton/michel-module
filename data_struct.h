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


#define HIT_BRANCHES(t, pre, h) \
    t->Branch(pre "HitSlice", &h.slice), \
    t->Branch(pre "HitZ", &h.z), \
    t->Branch(pre "HitChannel", &h.channel), \
    t->Branch(pre "HitTick", &h.tick), \
    t->Branch(pre "HitADC", &h.adc)

#define HITS_BRANCHES(t, pre, h) \
    t->Branch(pre "NHit", &h.N), \
    t->Branch(pre "HitSlice", &h.slice), \
    t->Branch(pre "HitZ", &h.z), \
    t->Branch(pre "HitChannel", &h.channel), \
    t->Branch(pre "HitTick", &h.tick), \
    t->Branch(pre "HitADC", &h.adc)

#define POINT_BRANCHES(t, pre, p) \
    t->Branch(pre "PointX", &p.x), \
    t->Branch(pre "PointY", &p.y), \
    t->Branch(pre "PointZ", &p.z)

#define POINTS_BRANCHES(t, pre, p) \
    t->Branch(pre "NPoint", &p.N), \
    t->Branch(pre "PointX", &p.x), \
    t->Branch(pre "PointY", &p.y), \
    t->Branch(pre "PointZ", &p.z)

#define LOG(x) (fLog ? std::cout << '\t' << #x << ": " << (x ? "\e[1;92m" : "\e[1;91m") << (x) << "\e[0m" << std::endl : void(), x)


namespace ana {
    template<typename T>
    struct bounds {
        T min, max;
        bounds() : min(std::numeric_limits<T>::max()), max(std::numeric_limits<T>::lowest()) {}
        bounds(T m, T M) : min(m), max(M) {}
        bool isInside(T x, float r=0) { 
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
        bounds3D(T xm, T xM, T ym, T yM, T zm, T zM) : x(xm, xM), y(ym, yM), z(zm, zM) {}
        bool isInside(T x, T y, T z, float r=0) {
            return this->x.isInside(x, r) && this->y.isInside(y, r) && this->z.isInside(z, r);
        }
        bool isInside(TLorentzVector const& v, float r=0) {
            return isInside(v.X(), v.Y(), v.Z(), r);
        }


        friend std::ostream& operator<<(std::ostream& os, const bounds3D& b) {
            return os << b.x << "x" << b.y << "x" << b.z;
        }
    };


    // some detector specific data
    bounds<float> tick_window;
    bounds3D<float> lower_bounds, upper_bounds;

    // std::map<raw::ChannelID_t, float> map_ch_z;
    // std::map<int,bounds<unsigned>> map_tpc_ch;
    std::map<int,std::pair<int,int>> map_sl_tpc = {
        {0, {0, 2}},
        {1, {1, 3}},
        {2, {4, 6}},
        {3, {5, 7}},
        {4, {8, 10}},
        {5, {9, 11}},
        {6, {12, 14}},
        {7, {13, 15}}
    };

    float fADCtoMeV = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
 


    // detector specific conversions
    unsigned GetSlice(raw::ChannelID_t ch, std::map<int,bounds<unsigned>>& map_tpc_ch) {
        bounds<unsigned> ch_bounds;
        for (auto const& [sl, tpcs] : map_sl_tpc) {
            ch_bounds = map_tpc_ch[tpcs.first];
            if (ch_bounds.min <= ch && ch <= ch_bounds.max) return sl;

            ch_bounds = map_tpc_ch[tpcs.second];
            if (ch_bounds.min <= ch && ch <= ch_bounds.max) return sl;
        }
        return -1;
    };
    float GetZ(raw::ChannelID_t ch, std::map<int,float>& map_ch_z) {
        return map_ch_z[ch];
    }

    struct Hit {
        // unsigned view;
        unsigned slice;
        float z;
        int channel;
        float tick;
        float adc;
        Hit() : slice(0), z(0), channel(0), tick(0), adc(0) {}
        Hit(unsigned s, float z, int c, float t, float a) :
            slice(s), z(z), channel(c), tick(t), adc(a) {}
        Hit(recob::Hit const& hit, std::map<int,bounds<unsigned>>& map_tpc_ch, std::map<int,float>& map_ch_z) :
            slice(GetSlice(hit.Channel(), map_tpc_ch)),
            z(GetZ(hit.Channel(), map_ch_z)),
            channel(hit.Channel()),
            tick(hit.PeakTime()),
            adc(hit.Integral()) {}


        friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
            return os << "sl:" << hit.slice << " z:" << hit.z << " ch:" << hit.channel << " tick:" << hit.tick << " ADC:" << hit.adc;
        }
    };
    struct Hits {
        unsigned N;
        // std::vector<unsigned> view;
        std::vector<unsigned> slice;
        std::vector<float> z;
        std::vector<int> channel;
        std::vector<float> tick;
        std::vector<float> adc;
        Hits() : N(0), slice(), z(), channel(), tick(), adc() {}
        void push_back(Hit const& hit) {
            N++;
            slice.push_back(hit.slice);
            z.push_back(hit.z);
            channel.push_back(hit.channel);
            tick.push_back(hit.tick);
            adc.push_back(hit.adc);
        }
        void clear() {
            N = 0;
            slice.clear();
            z.clear();
            channel.clear();
            tick.clear();
            adc.clear();
        }
        float energy() const {
            float e = 0;
            for (unsigned i=0; i<N; i++) e += adc[i];
            return e * fADCtoMeV;
        }

        // for range-loop
        Hit at(unsigned i) const { return Hit{slice[i], z[i], channel[i], tick[i], adc[i]}; }
        struct Hits_iterator {
            const Hits* hits;
            unsigned i;
            Hits_iterator(Hits const* h, unsigned i) : hits(h), i(i) {}
            Hits_iterator& operator++() { i++; return *this; }
            bool operator!=(Hits_iterator const& h) { return i != h.i; }
            Hit operator*() { return hits->at(i); }
        }
        Hits_iterator begin() const { return Hits_iterator(this, 0); }
        Hits_iterator end() const { return Hits_iterator(this, N); }
    };
    struct Point {
        float x, y, z;
        Point() : x(0), y(0), z(0) {}
        Point(float x, float y, float z) : x(x), y(y), z(z) {}
        Point(geo::Point_t const& p) : x(p.X()), y(p.Y()), z(p.Z()) {}
        Point operator-(Point const& p) { return Point{x-p.x, y-p.y, z-p.z}; }
        Point operator-(geo::Point_t const& p) { return Point{x-(float)p.x(), y-(float)p.y(), z-(float)p.z()}; }
        float r2() { return x*x + y*y + z*z; }


        friend std::ostream& operator<<(std::ostream& os, const Point& p) {
            return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
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
        void clear() {
            N = 0;
            x.clear();
            y.clear();
            z.clear();
        }

        // for range-loop
        Point at(unsigned i) const { return Point{x[i], y[i], z[i]}; }
        struct Points_iterator {
            const Points* points;
            unsigned i;
            Points_iterator(Points const* p, unsigned i) : points(p), i(i) {}
            Points_iterator& operator++() { i++; return *this; }
            bool operator!=(Points_iterator const& p) { return i != p.i; }
            Point operator*() { return points->at(i); }
        }
        Points_iterator begin() const { return Points_iterator(this, 0); }
        Points_iterator end() const { return Points_iterator(this, N); }
    };
}