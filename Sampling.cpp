/**
@file Sampling.cpp
@author t-sakai
@date 2019/09/15 create
*/
#include "Sampling.h"
#if LRAY_STANDALONE
#include <Windows.h>
#include "Sort.h"
#endif

namespace lray
{
#if LRAY_STANDALONE
    u32 bitreverse(u32 x)
    {
#if defined(__GNUC__)
        x = __builtin_bswap32(x);
#elif defined(__clang__)
        x = __builtin_bswap32(x);
#else
        x = (x << 16) | (x >> 16);
        x = ((x & 0x00FF00FFU) << 8) | ((x & 0xFF00FF00U) >> 8);
#endif
        x = ((x & 0x0F0F0F0FU) << 4) | ((x & 0xF0F0F0F0U) >> 4);
        x = ((x & 0x33333333U) << 2) | ((x & 0xCCCCCCCCU) >> 2);
        x = ((x & 0x55555555U) << 1) | ((x & 0xAAAAAAAAU) >> 1);
        return x;
    }

    namespace
    {
        static const u32 WELL_N = 16;

        u32 WELL_state[WELL_N];
        u32 WELL_index = 0;
    }

    f32 rand()
    {
        u32 a, b, c, d;

        a = WELL_state[WELL_index];
        c = WELL_state[(WELL_index + 13) & 15];
        b = a ^ c ^ (a << 16) ^ (c << 15);
        c = WELL_state[(WELL_index + 9) & 15];
        c ^= c >> 11;
        a = WELL_state[WELL_index] = b ^ c;
        d = a ^ ((a << 5) & 0xDA442D24UL);
        WELL_index = (WELL_index + 15) & 15;
        a = WELL_state[WELL_index];
        WELL_state[WELL_index] = a ^ b ^ d ^ (a << 2) ^ (b << 18) ^ (c << 28);

        u32 x = WELL_state[WELL_index];

        static const u32 m0 = 0x3F800000U;
        static const u32 m1 = 0x007FFFFFU;
        x = m0 | (x & m1);
        return (*(f32*)&x) - 1.000000000f;
    }

    void srand(u32 seed)
    {
        WELL_state[0] = seed;
        for(u32 i=1; i<WELL_N; ++i){
            WELL_state[i] = (1812433253 * (WELL_state[i-1]^(WELL_state[i-1] >> 30)) + i); 
        }
    }

    u32 cryptRandom()
    {
#if _WIN32
        u32 buffer = 0;
        HMODULE handle = LoadLibraryA("Advapi32.dll");
        if(NULL != handle){
            FARPROC procAddress = GetProcAddress(handle, "SystemFunction036");
            BOOLEAN result = (*(BOOLEAN (*)(PVOID, ULONG))procAddress)(&buffer, sizeof(u32));
            FreeLibrary(handle);
            if(!result){
                return 1234567U;
            }
        }

#else
        s32 fd = open("/dev/random", O_RDONLY);
        if(fd<0){
            return 1234567U;
        }
        u32 buffer = 0;
        read(fd, &buffer, sizeof(u32));
        close(fd);
#endif
        return buffer;
    }
#endif

    s32 gcd(s32 a, s32 b)
    {
        while(0!=b){
            s32 tmp = a%b;
            a = b;
            b = tmp;
        }
        return a;
    }

    void extended_euclid(s32& d, s32& x, s32& y, s32 a, s32 b)
    {
        if(0 == b){
            d = a;
            x = 1;
            y = 0;
        }else{
            s32 td, tx,ty;
            extended_euclid(td, tx, ty, b, a%b);
            d = td;
            x = ty;
            y = tx - (a/b) * ty;
        }
    }

    void gaussianBasisReduction(s64& b1x, s64& b1y, s64& b2x, s64& b2y, s32 gx, s32 gy, s32 n)
    {
        if(1==gx){
            b1x = gx;
            b1y = gy;
            b2x = 0;
            b2y = n;
        }else if(1==gy){
            b1x = gx;
            b1y = gy;
            b2x = n;
            b2y = 0;
        }else{
            s32 egx0, egx1, egx2;
            extended_euclid(egx0, egx1, egx2, n, gx);
            s32 egy0, egy1, egy2;
            extended_euclid(egy0, egy1, egy2, egx0, gy);

            s32 b1 = (n/egx0) * egy2;
            s32 b2 = (-gx/egx0) * egy2;
            b1x = egx2 * gx + egx1*n;
            b1y = egx2 * gy;
            b2x = b1 * gx + b2 * n;
            b2y = b1 * gy + egy1 * n;
        }

        s64 ib1 = b1x*b1x + b1y*b1y;
        s64 ib2 = b2x*b2x + b2y*b2y;

        for(;;){
            if(ib2<ib1){
                swap(b1x, b2x);
                swap(b1y, b2y);
            }
            f32 t = static_cast<f32>(b1x*b2x + b1y*b2y)/(b1x*b1x + b1y*b1y);
            s32 mu = static_cast<s32>(floorf(t+0.5f));
            b2x -= mu*b1x;
            b2y -= mu*b1y;
            ib1 = b1x*b1x + b1y*b1y;
            ib2 = b2x*b2x + b2y*b2y;
            if(ib1<=ib2){
                break;
            }
        }
    }

    void searchExhausitiveR1Lattice(s32& x, s32& y, s32 n)
    {
        s64 maxMin = 0;
        x=1;
        y=1;
        for(s32 j=1; j<n; ++j){
            for(s32 i=1; i<j; ++i){
                if(1 == gcd(gcd(i,j), n)){
                    s64 b1x,b1y,b2x,b2y;
                    gaussianBasisReduction(b1x, b1y, b2x, b2y, i, j, n);
                    s64 distance = b1x*b1x + b1y*b1y;
                    if(maxMin<distance){
                        maxMin = distance;
                        x=i;
                        y=j;
                    }
                }
            }
        }
    }

    //---------------------------------------------------------
    //--- R2
    //---------------------------------------------------------
    const f64 R2::G[3] = {
        1.0/1.61803398874989484820458683436563,
        1.0/1.32471795724474602596090885447809,
        1.0/1.22074408460575947536168534910883,
    };

    f32 R2::generate1(s32 index)
    {
        f32 a1 = static_cast<f32>(G[0]);
        f32 x = (0.5f + a1*index);
        return x-floorf(x);
    }

    Sample2 R2::generate2(s32 index)
    {
        f32 a1 = static_cast<f32>(G[1]);
        f32 a2 = static_cast<f32>(G[1]*G[1]);
        f32 x = (0.5f + a1*index);
        f32 y = (0.5f + a2*index);
        return {x-floorf(x), y-floorf(y)};
    }

    //---------------------------------------------------------
    //--- JitteredR2
    //---------------------------------------------------------
    JitteredR2::JitteredR2(s32 maxSamples)
        :a_(LRAY_NULL)
        ,s_(LRAY_NULL)
        ,c_(LRAY_NULL)
    {
        size_t maxSize = sizeof(s32) * maxSamples * 3;
        a_ = reinterpret_cast<s32*>(LRAY_MALLOC(maxSize));
        s_ = reinterpret_cast<s32*>(LRAY_MALLOC(maxSize));
        c_ = reinterpret_cast<s32*>(LRAY_MALLOC(maxSize));
    }

    JitteredR2::~JitteredR2()
    {
        LRAY_FREE(c_);
        LRAY_FREE(s_);
        LRAY_FREE(a_);
    }

    f32 JitteredR2::fractionalPowerX(s32 x, s32 y)
    {
        s32 n = x * y;
        size_t size = sizeof(s32) * n;

        memset(a_, 0, size);
        memset(s_, 0, size);
        a_[0] = 1;
        for(s32 j = 0; j < y; ++j){
            memset(c_, 0, size);
            s_[0] = a_[0];
            for(s32 i = 0; i < (n - 1); ++i){
                f32 z = static_cast<f32>(a_[i + 1] + a_[i] + c_[i]);
                c_[i + 1] = static_cast<s32>(z / x);
                s_[i + 1] = static_cast<s32>(z - c_[i + 1] * x); //z%x
            }
            memcpy(a_, s_, size);
        }
        f32 f = 0.0f;
        for(s32 i = 0; i < y; ++i){
            f += a_[i] * powf(x, i - y);
        }
        return f;
    }

    Sample2 JitteredR2::getU(s32 index)
    {
        static const f32 pi = 3.14159265359f;
        f32 v[2] = {fractionalPowerX(2, index + 1), fractionalPowerX(3, index + 1)};
        return {sqrtf(v[0]) * cosf(2.0f * pi * v[1]), sqrtf(v[0]) * sinf(2.0f * pi * v[1])};
    }

    Sample2 JitteredR2::generate2(s32 index, f32 lamda)
    {
        static const f32 delta0 = 0.76f;
        static const f32 i0 = 0.7f;
        static const f32 sqrt_pi = 1.77245385091f;
        Sample2 p = R2::generate2(index);
        Sample2 u = getU(index);

        f32 k = lamda * delta0 * sqrt_pi / (4.0f * sqrtf(static_cast<f32>(index) + i0));
        p.x_ += k*u.x_;
        p.y_ += k*u.y_;
        p.x_ -= floorf(p.x_);
        p.y_ -= floorf(p.y_);
        return p;
    }

    Sample2 JitteredR2::generateRandom2(s32 index, f32 lamda)
    {
        static const f32 delta0 = 0.76f;
        static const f32 i0 = 0.7f;
        static const f32 sqrt_pi = 1.77245385091f;
        Sample2 p = R2::generate2(index);
        f32 ux = lray::rand();
        f32 uy = lray::rand();

        f32 k = lamda * delta0 * sqrt_pi / (4.0f * sqrtf(static_cast<f32>(index) + i0));
        p.x_ += k*ux;
        p.y_ += k*uy;
        p.x_ -= floorf(p.x_);
        p.y_ -= floorf(p.y_);
        return p;
    }

    void hammersley(f32& x0, f32& x1, s32 index, s32 numPoints)
    {
        LRAY_ASSERT(0<=index && index<numPoints);
        x0 = static_cast<f32>(index)/numPoints;
        x1 = radicalInverseVanDerCorput(index);
    }

    //---------------------------------------------------------
    //--- ProgressiveJittered
    //---------------------------------------------------------
    void ProgressiveJittered::generate(s32 M, Sample2* samples)
    {
        LRAY_ASSERT(0==(M%4));
        LRAY_ASSERT(LRAY_NULL != samples);

        samples[0].x_ = lray::rand();
        samples[0].y_ = lray::rand();
        s32 N = 1;
        while(N<M){
            extendSequence(N, samples);
            N <<= 2;
        }
    }

    void ProgressiveJittered::extendSequence(s32 N, Sample2* samples)
    {
        f32 n = sqrtf(static_cast<f32>(N));
        f32 invn = 1.0f/n;
        for(s32 s=0; s<N; ++s){
            Sample2 old = samples[s];
            s32 i = static_cast<s32>(n * old.x_);
            s32 j = static_cast<s32>(n * old.y_);
            s32 xhalf = static_cast<s32>(2.0f * (n*old.x_ - i));
            s32 yhalf = static_cast<s32>(2.0f * (n*old.y_ - j));
            xhalf = 1-xhalf;
            yhalf = 1-yhalf;
            samples[N+s] = generateSamplePoint(i, j, xhalf, yhalf, invn);
            if(0.5f<rand()){
                xhalf = 1-xhalf;
            }else{
                yhalf = 1-yhalf;
            }
            samples[2*N+s] = generateSamplePoint(i, j, xhalf, yhalf, invn);
            xhalf = 1-xhalf;
            yhalf = 1-yhalf;
            samples[3*N+s] = generateSamplePoint(i, j, xhalf, yhalf, invn);
        }
    }

    Sample2 ProgressiveJittered::generateSamplePoint(s32 i, s32 j, s32 xhalf, s32 yhalf, f32 invn)
    {
        f32 x = (i + 0.5f * (xhalf + rand())) * invn;
        f32 y = (j + 0.5f * (yhalf + rand())) * invn;
        return {x, y};
    }

    //---------------------------------------------------------
    //--- Stratified
    //---------------------------------------------------------
    Stratified::Stratified(s32 maxSamples)
        :maxSamples_(maxSamples)
    {
        LRAY_ASSERT(0<=maxSamples_);
        s32 resolution = static_cast<s32>(ceilf(sqrtf(maxSamples_)));
        maxSamples_ = resolution*resolution;
        samples_ = LRAY_MALLOC_TYPE(Sample2, sizeof(Sample2)*resolution*resolution);
    }

    Stratified::~Stratified()
    {
        LRAY_FREE(samples_);
    }

    s32 Stratified::getResolution(s32 N) const
    {
        LRAY_ASSERT(0<=N && N<=maxSamples_);
        return static_cast<s32>(ceilf(sqrtf(N)));
    }

    void Stratified::generate(s32 N, Sample2* samples)
    {
        LRAY_ASSERT(LRAY_NULL != samples);
        s32 resolution = getResolution(N);
        f32 span = 1.0f/resolution;
        for(s32 i=0; i<resolution; ++i){
            f32 sy = static_cast<f32>(i)*span;
            for(s32 j=0; j<resolution; ++j){
                f32 x = static_cast<f32>(j)*span + span*rand();
                f32 y = sy + span*rand();
                samples_[resolution*i + j] = {x, y};
            }
        }
        shuffle(resolution);

        for(s32 i=0; i<N; ++i){
            samples[i] = samples_[i];
        }
    }

    void Stratified::shuffle(s32 resolution)
    {
        for(s32 i=0; i<resolution; ++i){
            s32 r = i*resolution;
            for(s32 j=2; j<resolution; ++j){
                s32 k = static_cast<s32>(rand()*j);
                swap(samples_[r+j], samples_[r+k]);
            }
        }

        //swap along rows
        for(s32 j = 0; j < resolution; ++j){
            for(s32 i = 2; i < resolution; ++i){
                s32 i0 = i*resolution + j;
                s32 i1 = static_cast<s32>(rand() * i) * resolution + j;
                swap(samples_[i0], samples_[i1]);
            }
        }
    }

#if 0
    //---------------------------------------------------------
    //--- ProgressiveMultiJittered
    //---------------------------------------------------------
    ProgressiveMultiJittered::ProgressiveMultiJittered(s32 maxSamples)
        :numSamples_(0)
        ,xhalves_(LRAY_NULL)
        ,yhalves_(LRAY_NULL)
        ,occupied1Dx_(LRAY_NULL)
        ,occupied1Dy_(LRAY_NULL)
    {
        xhalves_ = reinterpret_cast<f32*>(LRAY_MALLOC(sizeof(f32)*maxSamples));
        yhalves_ = reinterpret_cast<f32*>(LRAY_MALLOC(sizeof(f32)*maxSamples));
        occupied1Dx_ = reinterpret_cast<bool*>(LRAY_MALLOC(sizeof(bool)*maxSamples*2));
        occupied1Dy_ = reinterpret_cast<bool*>(LRAY_MALLOC(sizeof(bool)*maxSamples*2));
    }

    ProgressiveMultiJittered::~ProgressiveMultiJittered()
    {
        LRAY_FREE(occupied1Dy_);
        LRAY_FREE(occupied1Dy_);
        LRAY_FREE(yhalves_);
        LRAY_FREE(xhalves_);
    }

    void ProgressiveMultiJittered::generate(s32 M, Sample2* samples)
    {
        LRAY_ASSERT(0==(M%4));
        LRAY_ASSERT(LRAY_NULL != samples);

        numSamples_ = 0;
        samples[0] = {rand(), rand()};
        s32 N = 1;
        while(N<M){
            extendSequenceEven(N, samples);
            extendSequenceOdd(2*N, samples);
            N <<= 2;
        }
    }

    void ProgressiveMultiJittered::extendSequenceEven(s32 N, Sample2* samples)
    {
        s32 n = static_cast<s32>(sqrtf(N));
        markOccupiedStrata(N, samples);
        for(s32 s=0; s<N; ++s){
            Sample2 oldpt = samples[s];
            f32 i = floorf(n*oldpt.x_);
            f32 j = floorf(n*oldpt.y_);
            f32 xhalf = 1.0f - floorf(2.0f*(n*oldpt.x_-i));
            f32 yhalf = 1.0f - floorf(2.0f*(n*oldpt.y_-j));
            samples[numSamples_++] = generateSamplePoint(i, j, xhalf, yhalf, n, N, samples);
        }
    }

    void ProgressiveMultiJittered::extendSequenceOdd(s32 N, Sample2* samples)
    {
        s32 n = static_cast<s32>(sqrtf(N));
        markOccupiedStrata(N, samples);

        s32 halfN = N>>1;
        for(s32 s=0; s<halfN; ++s){
            Sample2 oldpt = samples[s];
            f32 i = floorf(n*oldpt.x_);
            f32 j = floorf(n*oldpt.y_);
            f32 xhalf = floorf(2.0f*(n*oldpt.x_-i));
            f32 yhalf = floorf(2.0f*(n*oldpt.y_-j));
            if(0.5f<rand()){
                xhalf = 1.0f-xhalf;
            }else{
                yhalf = 1.0f-yhalf;
            }
            xhalves_[s] = xhalf;
            yhalves_[s] = yhalf;
            samples[numSamples_++] = generateSamplePoint(i, j, xhalf, yhalf, n, N, samples);
        }

        for(s32 s=0; s<halfN; ++s){
            Sample2 oldpt = samples[s];
            f32 i = floorf(n*oldpt.x_);
            f32 j = floorf(n*oldpt.y_);
            f32 xhalf = 1.0f - xhalves_[s];
            f32 yhalf = 1.0f - yhalves_[s];
            samples[numSamples_++] = generateSamplePoint(i, j, xhalf, yhalf, n, N, samples);
        }
    }

    void ProgressiveMultiJittered::markOccupiedStrata(s32 N, Sample2* samples)
    {
        s32 NN = N<<1;
        for(s32 i=0; i<NN; ++i){
            occupied1Dx_[i] = false;
            occupied1Dy_[i] = false;
        }
        for(s32 s=0; s<N; ++s){
            s32 xstratum = static_cast<s32>(NN*samples[s].x_);
            s32 ystratum = static_cast<s32>(NN*samples[s].y_);
            occupied1Dx_[xstratum] = true;
            occupied1Dy_[ystratum] = true;
        }
    }

    Sample2 ProgressiveMultiJittered::generateSamplePoint(f32 i, f32 j, f32 xhalf, f32 yhalf, s32 n, s32 N, const Sample2* samples)
    {
        s32 NN = N<<1;
        f32 bestDist = 0.0f;
        s32 numCand = 10;

        s32 xstratum;
        s32 ystratum;

        Sample2 pt = {};
        Sample2 candpt;

        for(s32 t=1; t<numCand; ++t){
            do{
                candpt.x_ = (i+0.5f * (xhalf + rand()))/n;
                xstratum = static_cast<s32>(NN*candpt.x_);
            }while(occupied1Dx_[xstratum]);

            do{
                candpt.y_ = (j+0.5f * (yhalf + rand()))/n;
                ystratum = static_cast<s32>(NN*candpt.y_);
            }while(occupied1Dy_[ystratum]);

            f32 d = minDist(candpt, numSamples_, samples);
            if(bestDist<d){
                bestDist = d;
                pt = candpt;
            }
        }
        xstratum = static_cast<s32>(NN*pt.x_);
        ystratum = static_cast<s32>(NN*pt.y_);
        occupied1Dx_[xstratum] = true;
        occupied1Dy_[ystratum] = true;
        return pt;
    }

    f32 ProgressiveMultiJittered::minDist(const Sample2& candidate, s32 numSamples, const Sample2* samples)
    {
        f32 x = candidate.x_ - samples[0].x_;
        f32 y = candidate.y_ - samples[0].y_;
        f32 sqrDistance = x*x + y*y;
        for(s32 i=1; i<numSamples; ++i){
            f32 tx = candidate.x_ - samples[i].x_;
            f32 ty = candidate.y_ - samples[i].y_;
            f32 d = tx*tx + ty*ty;
            if(sqrDistance<d){
                sqrDistance = d;
            }
        }
        return sqrtf(sqrDistance);
    }

    //---------------------------------------------------------
    //--- CCVT
    //---------------------------------------------------------
    f32 CCVT::sqrDistance(const Sample2& x0, const Sample2& x1)
    {
        f32 dx = x0.x_ - x1.x_;
        f32 dy = x0.y_ - x1.y_;
        return dx*dx + dy*dy;
    }

    f32 CCVT::sqrDistance(const Point& x0, const Sample2& x1)
    {
        f32 dx = x0.x_ - x1.x_;
        f32 dy = x0.y_ - x1.y_;
        return dx*dx + dy*dy;
    }

    void CCVT::generate(s32 N, s32 M, Sample2* samples)
    {
        //create a point set
        Point* points = reinterpret_cast<Point*>(LRAY_MALLOC(sizeof(Point)*N*M));
        Site* sites = reinterpret_cast<Site*>(LRAY_MALLOC(sizeof(Site) * N));
        generateInitialSites(N, M, sites, points);

        Criterion* h0 = reinterpret_cast<Criterion*>(LRAY_MALLOC(sizeof(Criterion)*M));
        Criterion* h1 = reinterpret_cast<Criterion*>(LRAY_MALLOC(sizeof(Criterion)*M));
        bool stable;
        do{
            stable = true;
            for(s32 j = 1; j < N; ++j){
                Site& site1 = sites[j];
                Point* x1 = points + site1.top_;
                for(s32 i = 0; i < j; ++i){
                    Site& site0 = sites[i];
                    Point* x0 = points + site0.top_;
                    for(s32 k = 0; k < M; ++k){
                        h0[k].index_ = k;
                        h0[k].value_ = sqrDistance(x0[k], site0.central_) - sqrDistance(x0[k], site1.central_);
                        h1[k].index_ = k;
                        h1[k].value_ = sqrDistance(x1[k], site1.central_) - sqrDistance(x1[k], site0.central_);
                    }
                    introsort(M, h0, [](const Criterion& x0, const Criterion& x1){return x0.value_<x1.value_;});
                    introsort(M, h1, [](const Criterion& x0, const Criterion& x1){return x0.value_<x1.value_;});
                    s32 m = M-1;
                    while(0<=m && 0.0f<(h0[m].value_ + h1[m].value_)){
                        s32 k = h0[m].index_;
                        s32 l = h1[m].index_;
                        swap(x0[k], x1[l]);
                        --m;
                        stable = false;
                    }
                    if(!stable){
                        site0.central_ = updateCenter(i, M, points);
                        site1.central_ = updateCenter(j, M, points);
                    }
                }
            }
        } while(!stable);

        for(s32 i=0; i<N; ++i){
            samples[i] = sites[i].central_;
        }
        LRAY_FREE(h1);
        LRAY_FREE(h0);
        LRAY_FREE(points);
        LRAY_FREE(sites);
    }

    void CCVT::generateInitialSites(s32 N, s32 M, Site* sites, Point* points)
    {
        s32 top = 0;
        for(s32 i=0; i<N; ++i){
            sites[i].top_ = top;
            Point* s = points + top;
            for(s32 j=0; j<M; ++j){
                s[j].x_ = rand();
                s[j].y_ = rand();
            }
            sites[i].central_ = updateCenter(i, M,  points);
            top += M;
        }
    }

    Sample2 CCVT::updateCenter(s32 n, s32 M, Point* points)
    {
        s32 top = n * M;
        Point* s = points + top;
        Sample2 total = {};
        for(s32 i = 0; i < M; ++i){
            total.x_ += s[i].x_;
            total.y_ += s[i].y_;
        }
        f32 inv = 1.0f/M;
        Sample2 center = {total.x_*inv, total.y_*inv};
        return center;
    }

    namespace
    {
        class BitArray
        {
        public:
            static const s32 Max = 1024;

            BitArray();
            ~BitArray();

            void clear();
            bool get(s32 index) const;
            void set(s32 index);
            void reset(s32 index);

        private:
            u32 bits_[32];
        };
        BitArray::BitArray()
            :bits_{}
        {}

        BitArray::~BitArray()
        {}

        void BitArray::clear()
        {
            memset(bits_, 0, sizeof(u32)*32);
        }

        bool BitArray::get(s32 index) const
        {
            LRAY_ASSERT(0<=index && index<Max);
            s32 n = index>>5;
            s32 m = index - (n<<5);
            return (bits_[n] & (0x01U<<m)) == m;
        }

        void BitArray::set(s32 index)
        {
            LRAY_ASSERT(0<=index && index<Max);
            s32 n = index>>5;
            s32 m = index - (n<<5);
            bits_[n] &= ~(0x01U<<m);
        }

        void BitArray::reset(s32 index)
        {
            LRAY_ASSERT(0<=index && index<Max);
            s32 n = index>>5;
            s32 m = index - (n<<5);
            bits_[n] |= (0x01U<<m);
        }

        class StratifiedTree
        {
        public:
            static const s32 MaxBitCount = 10;
            struct Node
            {
                //BitArray occupied_;
                bool occupied_;
            };
            StratifiedTree(s32 bitCount);
            ~StratifiedTree();

            void add(f32 x);

            bool getValidOffsets(Sample2& offset, s32 x, s32 y, const Sample2& sample);
        private:
            StratifiedTree(const StratifiedTree&) = delete;
            StratifiedTree(StratifiedTree&&) = delete;
            StratifiedTree& operator=(const StratifiedTree&) = delete;
            StratifiedTree& operator=(StratifiedTree&&) = delete;

            void add(s32 node, f32 x, s32 nx);

            s32 bitCount_;
            s32 numNodes_;
            s32 leavesStart_;
            Node* nodes_;
        };

        StratifiedTree::StratifiedTree(s32 bitCount)
            :bitCount_(bitCount)
            ,numNodes_(0)
            ,leavesStart_(0)
            ,nodes_(LRAY_NULL)
        {
            LRAY_ASSERT(0<bitCount_ && bitCount_<=MaxBitCount);

            s32 start=(bitCount_>>1);
            for(s32 i=start; i<bitCount_; ++i){
                numNodes_ += 0x01U<<(i-start);
            }
            leavesStart_ = numNodes_;
            numNodes_ += 0x01U<<(bitCount_-start);

            size_t size = numNodes_ * sizeof(Node);
            nodes_ = LRAY_MALLOC_TYPE(Node, size);
            memset(nodes_, 0, size);
        }

        StratifiedTree::~StratifiedTree()
        {
            LRAY_FREE(nodes_);
        }

        void StratifiedTree::add(f32 x)
        {
            s32 xbits = (bitCount_ >> 1);

            s32 wx = 0x01<<xbits;
            x = wx * x;
            s32 nx = static_cast<s32>(x);
            add(0, x, nx);
        }

        void StratifiedTree::add(s32 node, f32 x, s32 nx)
        {
            if(leavesStart_<=node){
                if(static_cast<s32>(x) == nx){
                    nodes_[node].occupied_ = true;
                }
            }else{
                s32 left = (node+1)<<1;
                s32 right = left + 1;
                add(left, 2*x, 2*nx);
                add(right, 2*x+1, 2*nx);
                nodes_[node].occupied_ = nodes_[left].occupied_ || nodes_[right].occupied_;
            }
        }
    }
#endif
    //---------------------------------------------------------
    //---
    //---------------------------------------------------------
    f32 radicalInverseVanDerCorput(u32 index)
    {
        index = bitreverse(index);
        return static_cast<f32>(index) / static_cast<f32>(0x100000000L);
    }
}
