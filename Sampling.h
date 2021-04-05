#ifndef INC_LRAY_SAMPLING_H_
#define INC_LRAY_SAMPLING_H_
/**
@file Sampling.h
@author t-sakai
@date 2019/09/15 create

[1] Eric Heitz, Laurent Belcour, Victor Ostromoukhov, David Coeurjolly, and Jean-Claude Iehl,
"A Low-Discrepancy Sampler that Distributes Monte Carlo Errors as a Blue Noise in Screen Space",
ACM SIGGRAPH Talk 2019
https://belcour.github.io/blog/research/2019/06/17/sampling-bluenoise.html

[2] Hammersley Points on the Hemisphere
http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html

[3] Per Christensen, Andrew Kensler, Charlie Kilpatrick,
"Progressive Multi-Jittered Sample Sequences",
EGSR 2018,
https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/paper.pdf

[4] Matt Pharr,
"Efficient Generation of Points that Satisfy Two-Dimensional Elementary Intervals",
JCGT 2019,
http://jcgt.org/published/0008/01/04/

[4]
http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/

[5] Eric Heitz, Laurent Belcour, Victor Ostromoukhov, David Coeurjolly, and Jean-Claude Iehl,
"A Low-Discrepancy Sampler that Distributes Monte Carlo Errors as a Blue Noise in Screen Space",
ACM Siggraph Talk 2019,
https://belcour.github.io/blog/research/2019/06/17/sampling-bluenoise.html

[6] A. Ahmed and H. Perrier and D. Coeurjolly and V. Ostromoukhov and J. Guo and D. Yan and H. Huang and O. Deussen,
"Low-Discrepancy Blue Noise Sampling",
ACM Transactions on Graphics 2016,
http://graphics.uni-konstanz.de/publikationen/Ahmed2016LowdiscrepancyBlue/index.html
http://liris.cnrs.fr/ldbn/

[7]
"Blue Noise through Optimal Transport",
ACM Siggraph Asia 2012,
http://www.geometry.caltech.edu/BlueNoise/

[8]
"Sequences with Low-Discrepancy Blue-Noise 2-D Projections"
https://perso.liris.cnrs.fr/dcoeurjo/publications/perrier18eg.html

[9] Alexander Keller and Wolfgang Heidrich,
"Interleaved Sampling",
Proceedings of the Eurographics Workshop on Rendering 2001
http://www.cs.ubc.ca/nest/imager/tr/2001/keller2001a/keller.2001a.pdf

[10]
"Image Synthesis by Rank-1 Lattices",
https://link.springer.com/chapter/10.1007/978-3-540-74496-2_12

[11]
"Rank-1 Lattices in Computer Graphics",

[13]
https://www.slideshare.net/SAMSI_Info/program-on-quasimonte-carlo-and-highdimensional-sampling-methods-for-applied-mathematics-opening-workshop-lattice-rules-for-quasimonte-carlos-pierre-lecuyer-aug-28-2017
*/

#define LRAY_STANDALONE (1)

#if LRAY_STANDALONE
#include <cstdint>
#include <cassert>
#include <array>

#define LRAY_NULL nullptr
#define LRAY_ASSERT(exp) assert(exp)
#define LRAY_MALLOC(size) malloc(size)
#define LRAY_MALLOC_TYPE(type, size) reinterpret_cast<type*>(malloc(size))
#define LRAY_FREE(ptr) free(ptr);(ptr)=LRAY_NULL
#endif

namespace lray
{
#if LRAY_STANDALONE
    typedef int8_t s8;
    typedef int16_t s16;
    typedef int32_t s32;
    typedef int64_t s64;

    typedef uint8_t u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    typedef float f32;
    typedef double f64;

    using std::move;

    template<class T>
    void swap(T& x0, T& x1)
    {
        T tmp = x0;
        x0 = x1;
        x1 = tmp;
    }

    u32 bitreverse(u32 x);

    /**
    @brief [0 1)
    */
    f32 rand();
    void srand(u32 seed);
    u32 cryptRandom();
    
    s32 gcd(s32 a, s32 b);
    void extended_euclid(s32& d, s32& x, s32& y, s32 a, s32 b);
    void gaussianBasisReduction(s64& b1x, s64& b1y, s64& b2x, s64& b2y, s32 gx, s32 gy, s32 n);
    void searchExhausitiveR1Lattice(s32& x, s32& y, s32 n);
#endif

    struct Sample2
    {
        f32 x_;
        f32 y_;
    };

    //---------------------------------------------------------
    //---
    //--- Point sequences
    //---
    //---------------------------------------------------------

    //---------------------------------------------------------
    //--- Halton
    //---------------------------------------------------------
    template<s32 D>
    class Halton
    {
    public:
        static const s32 Dimension = D;

        /**
        @param generators ... generating vector whose all components are a prime number
        */
        Halton(const std::array<s32, Dimension>& generators);

        /**
        @param generators ... generating vector whose all components are a prime number
        @param start .. start index
        */
        Halton(const std::array<s32, Dimension>& generators, s32 start);

        void next(std::array<f32, Dimension>& x);
    private:
        f32 next(f32 prev, f32 invPrime) const;
        f32 at(s32 index, s32 prime) const;

        std::array<f32, Dimension> invPrimes_;
        std::array<f32, Dimension> previous_;
    };

    template<s32 D>
    Halton<D>::Halton(const std::array<s32, Dimension>& generators)
        :Halton(generators, 0)
    {}

    template<s32 D>
    Halton<D>::Halton(const std::array<s32, Dimension>& generators, s32 start)
    {
        for(s32 i=0; i<Dimension; ++i){
            invPrimes_[i] = 1.0f/generators[i];
            previous_[i] = at(start, generators[i]);
        }
    }

    template<s32 D>
    void Halton<D>::next(std::array<f32, Dimension>& x)
    {
        for(s32 i=0; i<Dimension; ++i){
            x[i] = previous_[i];
            previous_[i] = next(previous_[i], invPrimes_[i]);
        }
    }

    template<s32 D>
    f32 Halton<D>::next(f32 prev, f32 invPrime) const
    {
        float r = 1.0f - prev - 0.000001f;
        float f = invPrime;
        if(f < r) {
            return prev + f;
        } else {
            float h = f;
            float hh;
            do {
                hh = h;
                h *= f;
            } while(h >= r);
            return prev + hh + h - 1.0f;
        }
    }

    template<s32 D>
    f32 Halton<D>::at(s32 index, s32 prime) const
    {
        s32 i0 = index;
        f32 f = 1.0f/prime;
        f32 result = 0.0f;
        while(0<i0){
            s32 i1 = (i0%prime);
            f32 r = static_cast<f32>(i0)-i1*prime;
            result += f*r;
            i0 = i1;
        }
        return result;
    }

    //---------------------------------------------------------
    //--- R2
    //---------------------------------------------------------
    class R2
    {
    public:
        static f32 generate1(s32 index);
        static Sample2 generate2(s32 index);
    private:
        friend class JitteredR2;

        static const f64 G[3];
    };

    //---------------------------------------------------------
    //--- JitteredR2
    //---------------------------------------------------------
    class JitteredR2
    {
    public:
        explicit JitteredR2(s32 maxSamples);
        ~JitteredR2();

        Sample2 generate2(s32 index, f32 lamda=1.0f);
        static Sample2 generateRandom2(s32 index, f32 lamda=1.0f);
    private:
        JitteredR2(const JitteredR2&) = delete;
        JitteredR2(JitteredR2&&) = delete;
        JitteredR2& operator=(const JitteredR2&) = delete;
        JitteredR2& operator=(JitteredR2&&) = delete;

        f32 fractionalPowerX(s32 x, s32 y);
        Sample2 getU(s32 index);

        s32* a_;
        s32* s_;
        s32* c_;
    };

    //---------------------------------------------------------
    //---
    //--- Point sets
    //---
    //---------------------------------------------------------
    //---------------------------------------------------------
    //--- Rank1Lattice
    //---------------------------------------------------------
    template<s32 D>
    class Rank1Lattice
    {
    public:
        static const s32 Dimension = D;

        /**
        @param numPoints ... number of points to generate
        @param generators ... generating vector whose all components coprime to numPoints
        */
        Rank1Lattice(s32 numPoints, const std::array<s32, Dimension>& generators);

        /**
        @param numPoints ... number of points to generate
        @param generators ... generating vector whose all components coprime to numPoints
        @param shifts ... shift vector
        */
        Rank1Lattice(s32 numPoints, const std::array<s32, Dimension>& generators, const std::array<f32, Dimension>& shifts);

        void next(std::array<f32, Dimension>& x);
    private:

        s32 index_;
        s32 numPoints_;
        f32 inv_;
        std::array<s32, Dimension> generators_;
        std::array<f32, Dimension> shifts_;
    };

    template<s32 D>
    Rank1Lattice<D>::Rank1Lattice(s32 numPoints, const std::array<s32, Dimension>& generators)
        :index_(0)
        ,numPoints_(numPoints)
        ,inv_(1.0f/numPoints_)
        ,generators_(generators)
        ,shifts_{}
    {}

    template<s32 D>
    Rank1Lattice<D>::Rank1Lattice(s32 numPoints, const std::array<s32, Dimension>& generators, const std::array<f32, Dimension>& shifts)
        :Rank1Lattice(numPoints, generators)
        ,shifts_(shifts)
    {}

    template<s32 D>
    void Rank1Lattice<D>::next(std::array<f32, Dimension>& x)
    {
        f32 r = inv_*index_;
        for(s32 i=0; i<Dimension; ++i){
            f32 n = r * generators_[i];
            x[i] = n - ::floorf(n);
            LRAY_ASSERT(0.0f<=x[i] && x[i]<1.0f);
        }
        ++index_;
        if(numPoints_<=index_){
            index_ = 0;
        }
    }

    //---------------------------------------------------------
    //--- hammersley
    //---------------------------------------------------------
    void hammersley(f32& x0, f32& x1, s32 index, s32 numPoints);

    //---------------------------------------------------------
    //--- ProgressiveJittered
    //---------------------------------------------------------
    class ProgressiveJittered
    {
    public:
        static void generate(s32 M, Sample2* samples);
    private:
        static void extendSequence(s32 N, Sample2* samples);
        static Sample2 generateSamplePoint(s32 i, s32 j, s32 xhalf, s32 yhalf, f32 invn);
    };

    //---------------------------------------------------------
    //--- Stratified
    //---------------------------------------------------------
    class Stratified
    {
    public:
        Stratified(s32 maxSamples);
        ~Stratified();

        s32 getResolution(s32 N) const;
        void generate(s32 N, Sample2* samples);
        void shuffle(s32 resolution);
    private:
        Stratified(const Stratified&) = delete;
        Stratified(Stratified&&) = delete;
        Stratified& operator=(const Stratified&) = delete;
        Stratified& operator=(Stratified&&) = delete;

        s32 maxSamples_;
        Sample2* samples_;
    };

#if 0
    //---------------------------------------------------------
    //--- ProgressiveMultiJittered
    //---------------------------------------------------------
    /**
    @warn This is still a work-in-progress, so never use.
    */
    class ProgressiveMultiJittered
    {
    public:
        explicit ProgressiveMultiJittered(s32 maxSamples);
        ~ProgressiveMultiJittered();

        void generate(s32 M, Sample2* samples);
    private:
        ProgressiveMultiJittered(const ProgressiveMultiJittered&) = delete;
        ProgressiveMultiJittered(ProgressiveMultiJittered&&) = delete;
        ProgressiveMultiJittered& operator=(const ProgressiveMultiJittered&) = delete;
        ProgressiveMultiJittered& operator=(ProgressiveMultiJittered&&) = delete;

        void extendSequenceEven(s32 N, Sample2* samples);
        void extendSequenceOdd(s32 N, Sample2* samples);
        void markOccupiedStrata(s32 N, Sample2* samples);

        Sample2 generateSamplePoint(f32 i, f32 j, f32 xhalf, f32 yhalf, s32 n, s32 N, const Sample2* samples);
        static f32 minDist(const Sample2& candidate, s32 numSamples, const Sample2* samples);

        s32 numSamples_;
        f32* xhalves_;
        f32* yhalves_;
        bool* occupied1Dx_;
        bool* occupied1Dy_;
    };

    //---------------------------------------------------------
    //--- CCVT
    //---------------------------------------------------------
    /**
    @brief Capacity-Constrained Voronoi Tessellation
    */
    class CCVT
    {
    public:
        struct Site
        {
            Sample2 central_;
            s32 top_;
        };

        struct Point
        {
            f32 x_;
            f32 y_;
        };

        struct Criterion
        {
            s32 index_;
            f32 value_;
        };

        static f32 sqrDistance(const Sample2& x0, const Sample2& x1);
        static f32 sqrDistance(const Point& x0, const Sample2& x1);

        /**
        @param N ... number of sites
        @param M ... number of samples
        @param sites ... sites of voronoi
        */
        static void generate(s32 N, s32 M, Sample2* samples);
    private:
        static void generateInitialSites(s32 N, s32 M, Site* sites, Point* points);
        static Sample2 updateCenter(s32 n, s32 M, Point* points);

        s32 numSites_;
        s32 numSamples_;
        f32* h0_;
        f32* h1_;
    };
#endif

    //---------------------------------------------------------
    //---
    //---------------------------------------------------------
    f32 radicalInverseVanDerCorput(u32 index);

#if 0
    f32 vanDerCorput(u32 n, u32 base)
    {
        f32 vdc = 0.0f;
        f32 inv = 1.0f/base;
        f32 factor = inv;

        while(n){
            vdc += static_cast<f32>(n%base) * factor;
            n /= base;
            factor *= inv;
        }
        return vdc;
    }

    f32 radicalInverseSobol(u32 i, u32 scramble)
    {
        for(u32 v=1U<<31; i; i>>=1, v ^= v>>1){
            if(i&1){
                scramble ^= v;
            }
        }
        return static_cast<f32>(scramble) / static_cast<f32>(0x100000000L);
    }

    f32 radicalInverseLarcherPillichshammer(u32 i, u32 scramble)
    {
        for(u32 v=1U<<31; i; i>>=1, v |= v>>1){
            if(i&1){
                scramble ^= v;
            }
        }
        return static_cast<f32>(scramble) / static_cast<f32>(0x100000000L);
    }

    inline void sobol02(f32& v0, f32& v1, u32 i, u32 scramble0, u32 scramble1)
    {
        v0 = radicalInverseVanDerCorput(i, scramble0);
        v1 = radicalInverseSobol(i, scramble1);
    }
#endif
}
#endif  //INC_LRAY_SAMPLING_H_
