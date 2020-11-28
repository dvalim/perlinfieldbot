#include "../math/math.hpp"
#include "../math/rnd.hpp"
#include <vector>

#pragma once

#define R sqrt(v.x*v.x+v.y*v.y)
#define R2 (v.x*v.x+v.y*v.y)
#define O atan2(v.x, v.y)
#define PHI atan2(v.y, v.x)
#define PSY rnd::random(1.0)

namespace aural {
namespace flame {

extern double discard_chance;
double A, B, C, D, E, F;

math::vec sinusoidal(math::vec v, double weight) {
    return math::vec(weight*sin(v.x), weight*sin(v.y));
}

math::vec hyperbolic(math::vec v, double weight) {
    double r = R + 1.0e-10;
    double theta = atan2(v.x, v.y);
    double x = weight*sin(theta)/r;
    double y = weight*cos(theta)*r;
    return math::vec(x, y);
}

math::vec pdj(math::vec v, double weight) {
    return math::vec(weight*(sin(A*v.y)-cos(B*v.x)), weight*(sin(C*v.x)-cos(D*v.y)));
}

math::vec d_pdj(math::vec v, double weight) {
    double h = 0.1;
    double sqrth = sqrt(h);
    math::vec v1 = pdj(v, weight);
    math::vec v2 = pdj(math::vec(v.x+h, v.y+h), weight);
    return math::vec((v2.x-v1.x)/sqrth, (v2.y-v1.y)/sqrth);
}

math::vec julia(math::vec v, double weight) {
    double r = weight*sqrt(R);
    double theta = 0.5*atan2(v.x, v.y)+(int)(2.0*rnd::random(0.0,1.0))*M_PI;
    double x = r*cos(theta);
    double y = r*sin(theta);
    return math::vec(x, y);
}

double cosh(double x) { return 0.5*(exp(x)+exp(-x)); }
double sinh(double x) { return 0.5*(exp(x)-exp(-x)); }

math::vec sech(math::vec v, double weight) {
    double d = cos(2.0*v.y)+cosh(2.0*v.x);
    if(d) d = weight * 2.0 / d;
    return math::vec(d*cos(v.y)*cosh(v.x), -d*sin(v.y)*sinh(v.x));
}

math::vec linear(math::vec v, double weight) {
    return math::vec(v.x*weight, v.y*weight);
}

math::vec spherical(math::vec v, double weight) {
    return math::vec(weight*(1.0/(R*R))*v.x, weight*(1.0/(R*R))*v.y);
}

math::vec swirl(math::vec v, double weight) {
    auto xx = v.x*sin(R2)-v.y*cos(R2);
    auto yy = v.x*cos(R2)+v.y*sin(R2);
    return math::vec(weight*xx, weight*yy);
}

math::vec polar(math::vec v, double weight) {
    auto xx = O/M_PI;
    auto yy = R-1;
    return math::vec(weight*xx, weight*yy);
}

math::vec disc(math::vec v, double weight) {
    auto xx = O/M_PI*(sin(M_PI*R));
    auto yy = O/M_PI*(cos(M_PI*R));
    return math::vec(weight*xx, weight*yy);
}

math::vec heart(math::vec v, double weight) {
    auto xx = R*(sin(O*R));
    auto yy = R*(-cos(O*R));
    return math::vec(weight*xx, weight*yy);
}

math::vec spiral(math::vec v, double weight) {
    auto xx = (1.0/R)*(cos(O)+sin(R));
    auto yy = (1.0/R)*(sin(O)-cos(R));
    return math::vec(weight*xx, weight*yy);
}

math::vec ex(math::vec v, double weight) {
    auto p0 = sin(O+R), p1 = cos(O-R);
    auto xx = R*(p0*p0*p0+p1*p1*p1);
    auto yy = R*(p0*p0*p0-p1*p1*p1);
    return math::vec(weight*xx, weight*yy);
}

math::vec bent(math::vec v, double weight) {
    double xx, yy;
    if(v.x>=0 && v.y >= 0) xx = v.x, yy = v.y;
    else if(v.x<0&&v.y>=0) xx=2*v.x, yy = v.y;
    else if(v.x>=0&&v.y<0) xx = v.x, yy = v.y/2;
    else xx = 2*v.x, yy = v.y/2;
    return math::vec(weight*xx, weight*yy);
}

math::vec waves(math::vec v, double weight) {
    auto xx = v.x+B*sin(v.y/(C*C));
    auto yy = v.y+E*sin(v.x/(F*F));
    return math::vec(weight*xx, weight*yy);
}

math::vec fisheye(math::vec v, double weight) {
    auto xx = (2/(R+1))*v.y;
    auto yy = (2/(R+1))*v.x;
    return math::vec(weight*xx, weight*yy);
}

math::vec popcorn(math::vec v, double weight) {
    auto xx = v.x+C*sin(tan(3*v.y));
    auto yy = v.y+F*sin(tan(3*v.x));
    return math::vec(weight*xx, weight*yy);
}

math::vec exponential(math::vec v, double weight) {
    auto xx = exp(v.x-1)*cos(M_PI*v.y);
    auto yy = exp(v.x-1)*sin(M_PI*v.y);
    return math::vec(weight*xx, weight*yy);
}

math::vec cosine(math::vec v, double weight) {
    auto xx = cos(M_PI*v.x)*cosh(v.y);
    auto yy = -sin(M_PI*v.x)*sinh(v.y);
    return math::vec(weight*xx, weight*yy);
}

math::vec rings(math::vec v, double weight) {
    double mul = fmod((R+C*C),(2*C*C))-C*C+R*(1-C*C);
    auto xx = mul*cos(O);
    auto yy = mul*sin(O);
    return math::vec(weight*xx, weight*yy);
}

double fan_t;

math::vec fan(math::vec v, double weight) {
    double xx, yy;
    if(fmod(O+F, fan_t) > fan_t/2) {
        xx = R*(cos(O-fan_t/2));
        yy = R*sin(O-fan_t/2);
    } else {
        xx = R*(cos(O+fan_t/2));
        yy = R*sin(O+fan_t/2);
    }
    return math::vec(weight*xx, weight*yy);
}

math::vec cylinder(math::vec v, double weight) {
    auto xx = sin(v.x);
    auto yy = v.y;
    return math::vec(weight*xx, weight*yy);
}

math::vec gaussian(math::vec v, double weight) {
    double mul = PSY+PSY+PSY+PSY-2;
    double p5 = PSY;
    auto xx = mul*cos(2*M_PI*p5);
    auto yy = mul*sin(2*M_PI*p5);
    return math::vec(weight*xx, weight*yy);
}

double curl_p1, curl_p2;

math::vec curl(math::vec v, double weight) {
    double t1 = 1+curl_p1*v.x+curl_p2*(v.x*v.x-v.y*v.y), t2 = curl_p1*v.y+2*curl_p2*v.x*v.y;
    auto xx = (1/(t1*t1+t2*t2))*(v.x*t1+v.y*t2);
    auto yy = (1/(t1*t1+t2*t2))*(v.y*t1-v.x*t2);
    return math::vec(weight*xx, weight*yy);
}

math::vec rectangles(math::vec v, double weight) {
    double p1 = C, p2 = C;
    auto xx = (2*floor(v.x/p1)+1)*p1-v.x;
    auto yy = (2*floor(v.y/p2)+1)*p2-v.y;
    return math::vec(weight*xx, weight*yy);
}

math::vec horseshoe(math::vec v, double weight) {
    auto xx = 1/R * (v.x-v.y)*(v.x+v.y);
    auto yy = 1/R * 2*v.x*v.y;
    return math::vec(weight*xx, weight*yy);
}

math::vec handkerchief(math::vec v, double weight) {
    auto xx = R*(sin(O+R));
    auto yy = R*(cos(O-R));
    return math::vec(weight*xx, weight*yy);
}

math::vec diamond(math::vec v, double weight) {
    auto xx = sin(O)*cos(R);
    auto yy = cos(O)*sin(R);
    return math::vec(weight*xx, weight*yy);
}

math::vec power(math::vec v, double weight) {
    auto r = pow(R, sin(O));
    auto xx = r*cos(O);
    auto yy = r*sin(O);
    return math::vec(weight*xx, weight*yy);
}

math::vec bubble(math::vec v, double weight) {
    auto r = 4/(R*R+4);
    auto xx = r*v.x;
    auto yy = r*v.y;
    return math::vec(weight*xx, weight*yy);
}

math::vec tangent(math::vec v, double weight) {
    auto xx = sin(v.x)/cos(v.y);
    auto yy = tan(v.y);
    return math::vec(weight*xx, weight*yy);
}

math::vec square(math::vec v, double weight) {
    auto xx = PSY-0.5;
    auto yy = PSY-0.5;
    return math::vec(weight*xx, weight*yy);
}

math::vec cross(math::vec v, double weight) {
    auto r = sqrt(1/pow((v.x*v.x-v.y*v.y),2));
    auto xx = r*v.x;
    auto yy = r*v.y;
    return math::vec(weight*xx, weight*yy);
}

math::vec discard2(math::vec v, double weight) {
    double xx = v.x, yy = v.y;
    if(discard_chance < 50) xx = yy;
    else yy = xx;
    return math::vec(weight*xx, weight*yy);
}

math::vec blob(math::vec v, double weight) {
    double p1 = math::map(A, -1.5, 1.5, 0.8, 1.2), p2 = math::map(B, -1.5, 1.5, 0.8, 1.2), p3 = math::map(C, -1.5, 1.5, 0.8, 1.2);
    double xx = cos(O), yy = sin(O);
    double mult = p2+(p1-p2)*0.5*(sin(p3*O)+1);
    xx *= R * mult;
    yy *= R * mult;
    return math::vec(weight*xx, weight*yy);
}

math::vec ngon(math::vec v, double weight) {
    double p1 = 2, p2 = 2*M_PI/9, p3 = 5, p4 = 1;
    double t3 = PHI-p2*floor(PHI/p2);
    double t4 = (t3 > p2/2 ? t3 : t3 - p2);
    double k = p3*(1/cos(t4)-1)+p4;
    k /= pow(R, p1);
    double xx = v.x*k, yy = v.y*k;
    return math::vec(weight*xx, weight*yy);
}

math::vec secant(math::vec v, double weight) {
    return math::vec(weight*v.x, weight*(1/cos(R)));
}

math::vec hyperpdj(math::vec v, double weight) {
    return pdj(v, weight)+hyperbolic(v, weight);
}

math::vec hyperjulia(math::vec v, double weight) {
    return julia(v, weight)-hyperbolic(v, weight);
}

math::vec sechpdj(math::vec v, double weight) {
    return sech(v, weight)*pdj(v, weight);
}

math::vec julofhyp(math::vec v, double weight) {
    return julia(hyperbolic(v, weight), weight);
}

math::vec hyperten(math::vec v, double weight) {
    for(int i = 1; i <= 10; i++) v = hyperbolic(v, weight);
    return v;
}

math::vec spherewaves(math::vec v, double weight) {
    return spherical(waves(v, weight), weight);
}

math::vec juspiral(math::vec v, double weight) {
    return julia(spiral(v, weight), weight);
}

math::vec spherecorn(math::vec v, double weight) {
    return spherical(popcorn(v, weight), weight);
}

math::vec nspiral(math::vec v, double weight) {
    return ngon(spiral(v, weight), weight);
}

math::vec wavecant(math::vec v, double weight) {
    return waves(v, weight) - secant(v, weight);
}

#undef R
#undef R2
#undef O
#undef PSY

}
}