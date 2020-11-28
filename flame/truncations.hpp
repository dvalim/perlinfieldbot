#include "../math/math.hpp"
#include <vector>
#include <numeric>

#pragma once

namespace aural {
namespace flame {

double xpowy(math::vec v) {
    return pow(v.x+3, v.y+3);
}

double xmody(math::vec v) {
    return fmod(v.x, v.y);
}

double xdivy(math::vec v) {
    return v.y == 0 ? 0 : v.x/v.y;
}

double xmuly(math::vec v) {
    return v.x*v.y;
}

double sqrtfix(double n) {
    return sqrt(n+3);
}

double arctangent(math::vec v) {
    return atan2(v.x, v.y)*10;
}

double discard_chance;

double discard1(math::vec v) {
    return discard_chance < 50 ? v.x : v.y;
}

double dotproduct(math::vec v) {
    std::vector<double> v1 = {v.x, v.y}, v2 = {v.y, v.x};
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

double perlin(math::vec v) {
    return math::noise(v.x, v.y);
}

}
}