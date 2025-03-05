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

#include <cstdio>

#define LOG(x) (fLog ? printf("\t" #x ": " "\033[1;9%dm" "%s" "\033[0m\n", x?2:1, x?"true":"false") : 0, x)


namespace ana {
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


    float fADCtoMeV = 200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor

    struct Hit {
        // unsigned view;
        unsigned slice;
        float z;
        unsigned channel;
        float tick;
        float adc;
        Hit() : slice(0), z(0), channel(0), tick(0), adc(0) {}
        Hit(unsigned s, float z, unsigned c, float t, float a) :
            slice(s), z(z), channel(c), tick(t), adc(a) {}

        friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
            return os << "sl:" << hit.slice << " z:" << hit.z << " ch:" << hit.channel << " tick:" << hit.tick << " ADC:" << hit.adc;
        }
        void SetBranches(TTree* t, const char* pre) const {
            t->Branch(Form("%sHitSlice", pre), &slice);
            t->Branch(Form("%sHitZ", pre), &z);
            t->Branch(Form("%sHitChannel", pre), &channel);
            t->Branch(Form("%sHitTick", pre), &tick);
            t->Branch(Form("%sHitADC", pre), &adc);
        }
    };
    struct Hits {
        unsigned N;
        // std::vector<unsigned> view;
        std::vector<unsigned> slice;
        std::vector<float> z;
        std::vector<unsigned> channel;
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
        unsigned size() const { return N; }
        bool empty() const { return !N; }
        float energy() const {
            float e = 0;
            for (unsigned i=0; i<N; i++) e += adc[i];
            return e * fADCtoMeV;
        }

        void SetBranches(TTree *t, const char* pre) const {
            t->Branch(Form("%sNHit", pre), &N);
            t->Branch(Form("%sHitSlice", pre), &slice);
            t->Branch(Form("%sHitZ", pre), &z);
            t->Branch(Form("%sHitChannel", pre), &channel);
            t->Branch(Form("%sHitTick", pre), &tick);
            t->Branch(Form("%sHitADC", pre), &adc);
        }



        Hit at(unsigned i) const { return Hit{slice[i], z[i], channel[i], tick[i], adc[i]}; }
        struct iterator {
            const Hits* hits;
            unsigned i;
            iterator(Hits const* h, unsigned i) : hits(h), i(i) {}
            iterator& operator++() { i++; return *this; }
            bool operator!=(iterator const& h) const { return i != h.i; }
            Hit operator*() const { return hits->at(i); }
        };
        iterator begin() const { return iterator(this, 0); }
        iterator end() const { return iterator(this, N); }
    };
    struct Point {
        float x, y, z;
        Point() : x(0), y(0), z(0) {}
        Point(float x, float y, float z) : x(x), y(y), z(z) {}
        Point(double x, double y, double z) : x(x), y(y), z(z) {}
        Point(geo::Point_t const& p) : x(p.x()), y(p.y()), z(p.z()) {}
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

        void SetBranches(TTree* t, const char* pre) const {
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


        void SetBranches(TTree* t, const char* pre) const {
            t->Branch(Form("%sNPoint", pre), &N);
            t->Branch(Form("%sPointX", pre), &x);
            t->Branch(Form("%sPointY", pre), &y);
            t->Branch(Form("%sPointZ", pre), &z);
        }

        Point at(unsigned i) const { return Point{x[i], y[i], z[i]}; }
        struct iterator {
            const Points* points;
            unsigned i;
            iterator(Points const* p, unsigned i) : points(p), i(i) {}
            iterator& operator++() { i++; return *this; }
            bool operator!=(iterator const& p) { return i != p.i; }
            Point operator*() { return points->at(i); }
        };
        iterator begin() const { return iterator(this, 0); }
        iterator end() const { return iterator(this, N); }
    };
}