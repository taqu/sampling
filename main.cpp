#include "Sampling.h"
#include <stdio.h>

using namespace lray;
static const s32 TRY = 10;

typedef f64 (*FuncRMSE)(s32, Sample2*);

f64 gaussianRMSE(s32 N, Sample2* samples)
{
    f64 total = 0.0;
    for(s32 i=0; i<N; ++i){
        f64 x = samples[i].x_ - 0.5;
        f64 y = samples[i].y_ - 0.5;
        f64 d = exp(-(x*x + y*y));
        f64 v = d;
        total += v;
    }
    total /= N;
    return sqrt(::fabs(total-0.85112));//ground truth: 0.85112
}

void testRandom(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        for(s32 j = 0; j < TRY; ++j){
            for(s32 k = 0; k < i; ++k){
                samples[k].x_ = lray::rand();
                samples[k].y_ = lray::rand();
            }
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawRandom(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=0; i<N; ++i){
        samples[i].x_ = lray::rand();
        samples[i].y_ = lray::rand();
        fprintf(file, "%f,%f\n", samples[i].x_, samples[i].y_);
    }
    fclose(file);
}

void testHalton(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    Halton<2> halton({2, 3});
    std::array<f32, 2> sample;
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        for(s32 j = 0; j < TRY; ++j){
            for(s32 k = 0; k < i; ++k){
                halton.next(sample);
                samples[k].x_ = sample[0];
                samples[k].y_ = sample[1];
            }
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawHalton(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    Halton<2> halton({2, 3});
    std::array<f32, 2> sample;
    for(s32 i=0; i<N; ++i){
        halton.next(sample);
        fprintf(file, "%f,%f\n", sample[0], sample[1]);
    }
    fclose(file);
}

void testStratified(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    Stratified stratified(MaxN);
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        for(s32 j = 0; j < TRY; ++j){
            stratified.generate(i, samples);
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawStratified(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    Stratified stratified(N);
    stratified.generate(N, samples);
    for(s32 i=0; i<N; ++i){
        fprintf(file, "%f,%f\n", samples[i].x_, samples[i].y_);
    }
    fclose(file);
}

void testHammersley(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        for(s32 j = 0; j < TRY; ++j){
            for(s32 k = 0; k < i; ++k){
                hammersley(samples[k].x_, samples[k].y_, k, i);
            }
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawHammersley(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=0; i<N; ++i){
        hammersley(samples[i].x_, samples[i].y_, i, N);
        fprintf(file, "%f,%f\n", samples[i].x_, samples[i].y_);
    }
    fclose(file);
}

void testRank1Lattice(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    std::array<f32, 2> sample;
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        //Search generators
        s32 gx,gy;
        searchExhausitiveR1Lattice(gx, gy, i);
        for(s32 j = 0; j < TRY; ++j){
            Rank1Lattice<2> rank1Lattice(i, {gx,gy});
            for(s32 k = 0; k < i; ++k){
                rank1Lattice.next(sample);
                samples[k] = {sample[0], sample[1]};
            }
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawRank1Lattice(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    //Search generators
    s32 x, y;
    searchExhausitiveR1Lattice(x, y, N);
    Rank1Lattice<2> rank1Lattice(N, {x,y});
    std::array<f32, 2> sample;
    for(s32 i=0; i<N; ++i){
        rank1Lattice.next(sample);
        fprintf(file, "%f,%f\n", sample[0], sample[1]);
    }
    fclose(file);
}

void testProgressiveJittered(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        s32 N = (i+3)&~3L;
        for(s32 j = 0; j < TRY; ++j){
            ProgressiveJittered::generate(N, samples);
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawProgressiveJittered(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    ProgressiveJittered::generate(N, samples);
    for(s32 i=0; i<N; ++i){
        fprintf(file, "%f,%f\n", samples[i].x_, samples[i].y_);
    }
    fclose(file);
}

void testR2(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        for(s32 j = 0; j < TRY; ++j){
            for(s32 k = 0; k < i; ++k){
                samples[k] = R2::generate2(k);
            }
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawR2(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=0; i<N; ++i){
        samples[i] = R2::generate2(i);
        fprintf(file, "%f,%f\n", samples[i].x_, samples[i].y_);
    }
    fclose(file);
}

void testJitteredR2(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    JitteredR2 jitteredR2(MaxN);
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        for(s32 j = 0; j < TRY; ++j){
            for(s32 k = 0; k < i; ++k){
                samples[k] = jitteredR2.generate2(k);
            }
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawJitteredR2(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    JitteredR2 jitteredR2(N);
    for(s32 i=0; i<N; ++i){
        samples[i] = jitteredR2.generate2(i);
        fprintf(file, "%f,%f\n", samples[i].x_, samples[i].y_);
    }
    fclose(file);
}

void testRandomJitteredR2(s32 MaxN, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    for(s32 i=10; i<MaxN; i+=10){
        f64 rmse = 0.0;
        for(s32 j = 0; j < TRY; ++j){
            for(s32 k = 0; k < i; ++k){
                samples[k] = JitteredR2::generateRandom2(k);
            }
            rmse += gaussianRMSE(i, samples);
        }
        rmse /= TRY;
        f64 x = log10(i);
        fprintf(file, "%d, %lf, %lf\n", i, x, rmse);
    }
    fclose(file);
}

void drawRandomJitteredR2(s32 N, Sample2* samples, const char* filename)
{
    FILE* file = fopen(filename, "wb");
    JitteredR2 jitteredR2(N);
    for(s32 i=0; i<N; ++i){
        samples[i] = JitteredR2::generateRandom2(i);
        fprintf(file, "%f,%f\n", samples[i].x_, samples[i].y_);
    }
    fclose(file);
}

int main(int, char**)
{
    lray::srand(lray::cryptRandom());
    //
    static const s32 MaxSamples = 10000;
    Sample2* samples = new Sample2[MaxSamples*4];
    //testRandom(MaxSamples, samples, "random.csv");
    //testHalton(MaxSamples, samples, "halton.csv");
    //testStratified(MaxSamples, samples, "stratified.csv");
    //testHammersley(MaxSamples, samples, "hammersley.csv");
    testRank1Lattice(MaxSamples, samples, "rank1Lattice.csv");
    //testProgressiveJittered(MaxSamples, samples, "progressiveJittered.csv");
    //testR2(MaxSamples, samples, "R2.csv");
    ////testJitteredR2(MaxSamples, samples, "jitteredR2.csv");
    //testRandomJitteredR2(MaxSamples, samples, "randomJitteredR2.csv");

    //
    static const s32 DrawSamples = 1024;
    //drawRandom(DrawSamples, samples, "random_points.csv");
    //drawStratified(DrawSamples, samples, "stratified_points.csv");
    //drawHalton(DrawSamples, samples, "halton_points.csv");
    //drawHammersley(DrawSamples, samples, "hammersley_points.csv");
    drawRank1Lattice(DrawSamples, samples, "rank1Lattice_points.csv");
    //drawProgressiveJittered(DrawSamples, samples, "progressiveJittered_points.csv");
    //drawR2(DrawSamples, samples, "R2_points.csv");
    //drawJitteredR2(DrawSamples, samples, "jitteredR2_points.csv");
    //drawRandomJitteredR2(DrawSamples, samples, "randomJitteredR2_points.csv");
    delete[] samples;
    return 0;
}
