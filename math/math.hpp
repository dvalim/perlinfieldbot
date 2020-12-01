#include <math.h>
#include <algorithm>
#include "../FastNoise.h"


#pragma once

#define M_PI 3.14159265358979323846

namespace aural {
    namespace math {
        FastNoise noiseEngine;

        double constrain(double a, double b, double c) {
            if(a < b) return b;
            else if(a > c) return c;
            else return a;
        }

        double map(double val, double low1, double hi1, double low2, double hi2) { 
            return (val - low1) / (hi1 - low1) * (hi2 - low2) + low2;
        }	

        void initNoise(unsigned long long seed) {
            noiseEngine.SetNoiseType(FastNoise::PerlinFractal);
            noiseEngine.SetFractalOctaves(3);
            noiseEngine.SetSeed(seed);
        }

        double noise(double x) {
            return noiseEngine.GetNoise(x, 0, 0);
        }

        double noise(double x, double y) {
            return noiseEngine.GetNoise(x, y, 0);
        }

        double noise(double x, double y, double z) {
            return noiseEngine.GetNoise(x, y, z);
        }

        double interpolate(double val1, double val2, double percent) {
            double result = val1 + percent * (val2 - val1);
            return result;
        }

        struct vec {
            double x, y;
            vec() : x(0), y(0) {}
            vec(double x, double y) : x(x), y(y) {}
            double mag() {
                return sqrt(x*x+y*y);
            }
            double mag2() {
                return x*x+y*y;
            }
            double arctan() {
                return atan2(y, x);
            }
        };
        vec operator+(vec const& lhs, vec const& rhs) {
            return {lhs.x+rhs.x, lhs.y+rhs.y};
        }
        vec operator-(vec const& lhs, vec const& rhs) {
            return {lhs.x-rhs.x, lhs.y-rhs.y};
        }
        vec operator*(vec const& lhs, vec const& rhs) {
            return {lhs.x*rhs.x, lhs.y*rhs.y};
        }
        vec operator/(vec const& lhs, vec const& rhs) {
            return {rhs.x == 0 ? 0 : lhs.x/rhs.x, rhs.y == 0 ? 0 : lhs.y/rhs.y};
        }

        vec operator+(vec const& lhs, double const& rhs) {
            return {lhs.x+rhs, lhs.y+rhs};
        }
        vec operator-(vec const& lhs, double const& rhs) {
            return {lhs.x-rhs, lhs.y-rhs};
        }
        vec operator*(vec const& lhs, double const& rhs) {
            return {lhs.x*rhs, lhs.y*rhs};
        }
        vec operator/(vec const& lhs, double const& rhs) {
            return {rhs == 0 ? 0 : lhs.x/rhs, rhs == 0 ? 0 : lhs.y/rhs};
        }

        vec operator+(double const& lhs, vec const& rhs) {
            return {lhs+rhs.x, lhs+rhs.y};
        }
        vec operator-(double const& lhs, vec const& rhs) {
            return {lhs-rhs.x, lhs-rhs.y};
        }
        vec operator*(double const& lhs, vec const& rhs) {
            return {lhs*rhs.x, lhs*rhs.y};
        }
        vec operator/(double const& lhs, vec const& rhs) {
            return {lhs == 0 ? 0 : rhs.x/lhs, lhs == 0 ? 0 : rhs.y/lhs};
        }


    }
}
