#include <random>
#include <chrono>

#pragma once

namespace aural {
    namespace rnd {
        std::mt19937 engine;

        unsigned long long seedEngine() {
	        unsigned long long seed = std::chrono::system_clock::now().time_since_epoch().count();
            engine.seed(seed);
            return seed;
        }

        unsigned long long seedEngine(unsigned int seed) {
            engine.seed(seed);
            return seed;
        }

        int random(int lower, int upper) {
            std::uniform_int_distribution<int> uid(lower, upper);
            return uid(engine);
        }

        int random(int upper) {
            return random(0, upper);
        }

        double random(double lower, double upper) {
            std::uniform_real_distribution<double> urd(lower, upper);
            return urd(engine);
        }

        double random(double upper) {
            return random(0.0, upper);
        }

        double gaussian(double mean, double deviation) {
            std::normal_distribution<double> nd(mean, deviation);
            return nd(engine);
        }
    }
}