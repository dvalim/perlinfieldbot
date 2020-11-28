#include "../math/math.hpp"
#include "../math/rnd.hpp"

#pragma once

#define R sqrt(v.x*v.x+v.y*v.y)
#define R2 (v.x*v.x+v.y*v.y)
#define O atan2(v.y, v.x)
#define PSY rnd::random(1.0)

namespace aural {
namespace flame {

math::vec cissoid(double n, double weight) {
    double sinn2 = 2*sin(n)*sin(n);
    double xt = sinn2;
    double yt = sinn2*tan(n);
    return math::vec(weight*xt, weight*yt);
}

math::vec kampyle(double n, double weight) {
    double sec = 1/sin(n);
    double xt = sec;
    double yt = tan(n)*sec; 
    return math::vec(weight*xt, weight*yt);
}

math::vec astroid(double n, double weight) {
    double sinn = sin(n);
    double cosn = cos(n);
    double xt = sinn*sinn*sinn;
    double yt = cosn*cosn*cosn;
    return math::vec(weight*xt, weight*yt);
}

math::vec hyperbole(double n, double weight) {
    double xt = 1/sin(n);
    double yt = tan(n); 
    return math::vec(weight*xt, weight*yt);
}

math::vec superformula(double n, double weight) {
    double a = 1;
    double b = 1;
    double m = 6;
    double n1 = 1;
    double n2 = 7;
    double n3 = 8;
    
    double f1 = pow(abs(cos(m*n/4)/a),n2);
    double f2 = pow(abs(sin(m*n/4)/b),n3);
    double fr = pow(f1+f2,-1/n1);
    
    double xt = cos(n)*fr;
    double yt = sin(n)*fr;
    return math::vec(weight*xt, weight*yt);
}

math::vec cartesianXY(double n, double weight) {
    double offx = rnd::random(0, 100) < 50 ? -1 : 1, offy = rnd::random(0, 100) < 50 ? -1 : 1;
    return math::vec(n*weight*offx, n*weight*offy);
}

math::vec polarXY(double n, double weight) {
    double xt = cos(M_PI*2*n);
    double yt = sin(M_PI*2*n); 
    return math::vec(weight*xt, weight*yt);
}

math::vec copy(double n, double weight) {
    return math::vec(n*weight, n*weight);
}

math::vec archispiral(double n, double weight) {
    double xt = n*cos(n), yt = n*sin(n);
    return math::vec(weight*xt, weight*yt);
}

math::vec cycloid(double n, double weight) {
    double xt = cos(n)*2*cos(2*n)+1, yt = sin(n)*2*sin(2*n)+1;
    return math::vec(weight*xt, weight*yt);
}

math::vec deltoid(double n, double weight) {
    double xt = (2*cos(n))/3+(cos(2*n)/3), yt = (2*sin(n))/3-(sin(2*n)/3);
    return math::vec(weight*xt, weight*yt);
}

math::vec folium(double n, double weight) {
    double xt = (3*n)/(pow(n, 3)+1), yt = (3*n*n)/(pow(n, 3)+1);
    return math::vec(weight*xt, weight*yt);
}

double trochoid_b, trochoid_h;

math::vec epitrochoid(double n, double weight) {
    double a = 16, b = trochoid_b, h = trochoid_h;
    double xt = (a+b)*cos(n)-h*cos((n*(a+b))/b);
    double yt = (a+b)*sin(n)-h*sin((n*(a+b))/b);
    return math::vec(weight*xt, weight*yt);
}

math::vec hypotrochoid(double n, double weight) {
    double a = 16, b = trochoid_b, h = trochoid_h;
    double xt = (a-b)*cos(n)+h*cos((n*(a-b))/b);
    double yt = (a-b)*sin(n)-h*sin((n*(a-b))/b);
    return math::vec(weight*xt, weight*yt);
}

double sinc(double x) { return sin(x)/x; }

math::vec quadratrix(double n, double weight) {
    double xt = (2*cos(n))/(M_PI*sinc(n));
    double yt = (2*sin(n))/(M_PI*sinc(n));
    return math::vec(weight*xt, weight*yt);
}

math::vec kilroy(double n, double weight) {
    double xt = n;
    double yt = log(abs(sinc(n)));
    return math::vec(weight*xt, weight*yt);
}

math::vec lituus(double n, double weight) {
    double xt = cos(n)/sqrt(n);
    double yt = sin(n)/sqrt(n);
    return math::vec(weight*xt, weight*yt);
}

math::vec maltesecross(double n, double weight) {
    double xt = (2*cos(n))/(sqrt(sin(4*n)));
    double yt = (2*sin(n))/(sqrt(sin(4*n)));
    return math::vec(weight*xt, weight*yt);
}

math::vec neoid(double n, double weight) {
    double xt = cos(n)*(n+1);
    double yt = sin(n)*(n+1);
    return math::vec(weight*xt, weight*yt);
}

math::vec sextic(double n, double weight) {
    double xt = pow(cos(n/3), 3)*cos(n);
    double yt = sin(n)*pow(cos(n/3), 3);
    return math::vec(weight*xt, weight*yt);
}

math::vec eightcurve(double n, double weight) {
    double xt = sin(n);
    double yt = sin(n)*cos(n);
    return math::vec(weight*xt, weight*yt);
}

double ecycloid_a, ecycloid_b, hcycloid_a, hcycloid_b;

math::vec epicycloid(double n, double weight) {
    double a = ecycloid_a, b = ecycloid_b;
    double xt = (a+b)*cos(n)-b*cos(n*(a+b)/b);
    double yt = (a+b)*sin(n)-b*sin(n*(a+b)/b);
    return math::vec(weight*xt, weight*yt);
}

math::vec hypocycloid(double n, double weight) {
    double a = hcycloid_a, b = hcycloid_b;
    double xt = (a-b)*cos(n)+b*cos(n*(a-b)/b);
    double yt = (a-b)*sin(n)-b*sin(n*(a-b)/b);
    return math::vec(weight*xt, weight*yt);
}

math::vec trisectrix(double n, double weight) {
    double xt = (n*n-3)/(n*n+1);
    double yt = n*(n*n-3)/(n*n+1);
    return math::vec(weight*xt, weight*yt);
}

math::vec quadrifolium(double n, double weight) {
    double xt = sin(2*n)*cos(n);
    double yt = sin(n)*sin(2*n);
    return math::vec(weight*xt, weight*yt);
}

math::vec tschirnhausen(double n, double weight) {
    double xt = 1-3*n*n;
    double yt = n*(3-n*n);
    return math::vec(weight*xt, weight*yt);
}

math::vec cornoid(double n, double weight) {
    auto xx = (1-2*sin(n)*sin(n))*cos(n);
    auto yy = sin(n)*(2*cos(n)*cos(n)+1);
    return math::vec(xx*weight, yy*weight);
}

math::vec fishcurve(double n, double weight) {
    auto xx = cos(n)-(sin(n)*sin(n))/sqrt(2);
    auto yy = sin(n)*cos(n);
    return math::vec(xx*weight, yy*weight);
}

math::vec ranunculoid(double n, double weight) {
    auto xx = 6*cos(n)-cos(6*n);
    auto yy = 6*sin(n)-sin(6*n);
    return math::vec(xx*weight, yy*weight);
}

#undef R
#undef R2
#undef O
#undef PSY

}
}