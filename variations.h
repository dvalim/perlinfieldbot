#include "sfmlpp.h"
#include "FastNoise.h"
#define R sqrt(v.x*v.x+v.y*v.y)
#define R2 (v.x*v.x+v.y*v.y)
#define O atan2(v.y, v.x)
#define PHI atan2(v.x, v.y)
#define OM rd(0, 100) < 50 ? 0 : M_PI
#define TRI rd(0, 100) < 50 ? -1 : 1
#define PSY rd(0, 1)

double A, B, C, D, E, F;

sf::Vector2f sinusoidal(sf::Vector2f v, double amount) {
    return sf::Vector2f(amount*sin(v.x), amount*sin(v.y));
}

sf::Vector2f hyperbolic(sf::Vector2f v, double amount) {
    double r = R + 1.0e-10;
    double theta = atan2(v.x, v.y);
    double x = amount*sin(theta)/r;
    double y = amount*cos(theta)*r;
    return sf::Vector2f(x, y);
}

sf::Vector2f pdj(sf::Vector2f v, double amount) {
    return sf::Vector2f(amount*(sin(A*v.y)-cos(B*v.x)), amount*(sin(C*v.x)-cos(D*v.y)));
}

sf::Vector2f d_pdj(sf::Vector2f v, double amount) {
    double h = 0.1;
    double sqrth = sqrt(h);
    sf::Vector2f v1 = pdj(v, amount);
    sf::Vector2f v2 = pdj(sf::Vector2f(v.x+h, v.y+h), amount);
    return sf::Vector2f((v2.x-v1.x)/sqrth, (v2.y-v1.y)/sqrth);
}

sf::Vector2f julia(sf::Vector2f v, double amount) {
    double r = amount*sqrt(R);
    double theta = 0.5*atan2(v.x, v.y)+(int)(2.0*rd(0,1))*M_PI;
    double x = r*cos(theta);
    double y = r*sin(theta);
    return sf::Vector2f(x, y);
}

double cosh_(double x) { return 0.5*(exp(x)+exp(-x)); }
double sinh_(double x) { return 0.5*(exp(x)-exp(-x)); }

sf::Vector2f sech(sf::Vector2f v, double amount) {
    double d = cos(2.0*v.y)+cosh_(2.0*v.x);
    if(d) d = amount * 2.0 / d;
    return sf::Vector2f(d*cos(v.y)*cosh_(v.x), -d*sin(v.y)*sinh_(v.x));
}

sf::Vector2f addF(sf::Vector2f v1, sf::Vector2f v2) { return sf::Vector2f(v1.x+v2.x, v1.y+v2.y); }
sf::Vector2f subF(sf::Vector2f v1, sf::Vector2f v2) { return sf::Vector2f(v1.x-v2.x, v1.y-v2.y); }
sf::Vector2f mulF(sf::Vector2f v1, sf::Vector2f v2) { return sf::Vector2f(v1.x*v2.x, v1.y*v2.y); }
sf::Vector2f divF(sf::Vector2f v1, sf::Vector2f v2) { return sf::Vector2f(v2.x == 0 ? 0 : v1.x/v2.x, v2.y == 0 ? 0 : v1.y/v2.y); }

typedef sf::Vector2f(*variation)(sf::Vector2f, double);
typedef sf::Vector2f(*extrapolation)(double, double);
typedef double(*truncation)(sf::Vector2f);
typedef double(*transformation)(double);

map<string, variation> variations;
map<string, extrapolation> extrapolations;
map<string, truncation> truncations;
map<string, transformation> transformations;

string randVariation() {
    vector<string> t;
    for(auto i : variations)
        t.push_back(i.first);
    return t[rdint(0, t.size()-1)];
}

string randExtrapolation() {
    vector<string> t;
    for(auto i : extrapolations)
        t.push_back(i.first);
    return t[rdint(0, t.size()-1)];
}

char randOperation() {
    vector<char> t = {'+', '-', '*', '/'};
    return t[rdint(0, t.size()-1)];
}

sf::Vector2f linear(sf::Vector2f v, double amount) {
    return sf::Vector2f(v.x*amount, v.y*amount);
}

sf::Vector2f spherical(sf::Vector2f v, double amount) {
    return sf::Vector2f(amount*(1.0/(R*R))*v.x, amount*(1.0/(R*R))*v.y);
}

sf::Vector2f swirl(sf::Vector2f v, double amount) {
    auto xx = v.x*sin(R2)-v.y*cos(R2);
    auto yy = v.x*cos(R2)+v.y*sin(R2);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f polar(sf::Vector2f v, double amount) {
    auto xx = O/M_PI;
    auto yy = R-1;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f disc(sf::Vector2f v, double amount) {
    auto xx = O/M_PI*(sin(M_PI*R));
    auto yy = O/M_PI*(cos(M_PI*R));
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f heart(sf::Vector2f v, double amount) {
    auto xx = R*(sin(O*R));
    auto yy = R*(-cos(O*R));
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f spiral(sf::Vector2f v, double amount) {
    auto xx = (1.0/R)*(cos(O)+sin(R));
    auto yy = (1.0/R)*(sin(O)-cos(R));
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f ex(sf::Vector2f v, double amount) {
    auto p0 = sin(O+R), p1 = cos(O-R);
    auto xx = R*(p0*p0*p0+p1*p1*p1);
    auto yy = R*(p0*p0*p0-p1*p1*p1);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f bent(sf::Vector2f v, double amount) {
    double xx, yy;
    if(v.x>=0 && v.y >= 0) xx = v.x, yy = v.y;
    else if(v.x<0&&v.y>=0) xx=2*v.x, yy = v.y;
    else if(v.x>=0&&v.y<0) xx = v.x, yy = v.y/2;
    else xx = 2*v.x, yy = v.y/2;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f waves(sf::Vector2f v, double amount) {
    auto xx = v.x+B*sin(v.y/(C*C));
    auto yy = v.y+E*sin(v.x/(F*F));
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f fisheye(sf::Vector2f v, double amount) {
    auto xx = (2/(R+1))*v.y;
    auto yy = (2/(R+1))*v.x;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f popcorn(sf::Vector2f v, double amount) {
    auto xx = v.x+C*sin(tan(3*v.y));
    auto yy = v.y+F*sin(tan(3*v.x));
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f exponential(sf::Vector2f v, double amount) {
    auto xx = exp(v.x-1)*cos(M_PI*v.y);
    auto yy = exp(v.x-1)*sin(M_PI*v.y);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f cosine(sf::Vector2f v, double amount) {
    auto xx = cos(M_PI*v.x)*cosh_(v.y);
    auto yy = -sin(M_PI*v.x)*sinh_(v.y);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f rings(sf::Vector2f v, double amount) {
    double mul = fmod((R+C*C),(2*C*C))-C*C+R*(1-C*C);
    auto xx = mul*cos(O);
    auto yy = mul*sin(O);
    return sf::Vector2f(amount*xx, amount*yy);
}

double fan_t;

sf::Vector2f fan(sf::Vector2f v, double amount) {
    double xx, yy;
    if(fmod(O+F, fan_t) > fan_t/2) {
        xx = R*(cos(O-fan_t/2));
        yy = R*sin(O-fan_t/2);
    } else {
        xx = R*(cos(O+fan_t/2));
        yy = R*sin(O+fan_t/2);
    }
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f cylinder(sf::Vector2f v, double amount) {
    auto xx = sin(v.x);
    auto yy = v.y;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f gaussian(sf::Vector2f v, double amount) {
    double mul = PSY+PSY+PSY+PSY-2;
    double p5 = PSY;
    auto xx = mul*cos(2*M_PI*p5);
    auto yy = mul*sin(2*M_PI*p5);
    return sf::Vector2f(amount*xx, amount*yy);
}

double curl_p1, curl_p2;

sf::Vector2f curl(sf::Vector2f v, double amount) {
    double t1 = 1+curl_p1*v.x+curl_p2*(v.x*v.x-v.y*v.y), t2 = curl_p1*v.y+2*curl_p2*v.x*v.y;
    auto xx = (1/(t1*t1+t2*t2))*(v.x*t1+v.y*t2);
    auto yy = (1/(t1*t1+t2*t2))*(v.y*t1-v.x*t2);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f rectangles(sf::Vector2f v, double amount) {
    double p1 = C, p2 = C;
    auto xx = (2*floor(v.x/p1)+1)*p1-v.x;
    auto yy = (2*floor(v.y/p2)+1)*p2-v.y;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f horseshoe(sf::Vector2f v, double amount) {
    auto xx = 1/R * (v.x-v.y)*(v.x+v.y);
    auto yy = 1/R * 2*v.x*v.y;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f handkerchief(sf::Vector2f v, double amount) {
    auto xx = R*(sin(O+R));
    auto yy = R*(cos(O-R));
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f diamond(sf::Vector2f v, double amount) {
    auto xx = sin(O)*cos(R);
    auto yy = cos(O)*sin(R);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f power(sf::Vector2f v, double amount) {
    auto r = pow(R, sin(O));
    auto xx = r*cos(O);
    auto yy = r*sin(O);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f bubble(sf::Vector2f v, double amount) {
    auto r = 4/(R*R+4);
    auto xx = r*v.x;
    auto yy = r*v.y;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f tangent(sf::Vector2f v, double amount) {
    auto xx = sin(v.x)/cos(v.y);
    auto yy = tan(v.y);
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f square(sf::Vector2f v, double amount) {
    auto xx = PSY-0.5;
    auto yy = PSY-0.5;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f cross(sf::Vector2f v, double amount) {
    auto r = sqrt(1/pow((v.x*v.x-v.y*v.y),2));
    auto xx = r*v.x;
    auto yy = r*v.y;
    return sf::Vector2f(amount*xx, amount*yy);
}

sf::Vector2f cissoid(double n, double amount) {
    double sinn2 = 2*sin(n)*sin(n);
    double xt = sinn2;
    double yt = sinn2*tan(n);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f kampyle(double n, double amount) {
    double sec = 1/sin(n);
    double xt = sec;
    double yt = tan(n)*sec; 
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f astroid(double n, double amount) {
    double sinn = sin(n);
    double cosn = cos(n);
    double xt = sinn*sinn*sinn;
    double yt = cosn*cosn*cosn;
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f hyperbole(double n, double amount) {
    double xt = 1/sin(n);
    double yt = tan(n); 
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f superformula(double n, double amount) {
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
    return sf::Vector2f(amount*xt, amount*yt);
}

#ifdef NOISE
extern FastNoise pn;
extern double noise_scale, coord_scale;
//#define noise_scale noise_magnitudes[rdint(0, noise_magnitudes.size()-1)]
//#define coord_scale coord_magnitudes[rdint(0, coord_magnitudes.size()-1)]
extern vector<double> noise_magnitudes, coord_magnitudes;

double perlin(sf::Vector2f v) {
    //auto temp = rd(0.5, 3);
    return noise_scale*pn.GetNoise(v.x*coord_scale, v.y*coord_scale);
}

double noise(double n) {
    return noise_scale*pn.GetNoise(n*coord_scale, n*coord_scale);
}

sf::Vector2f noise3(double n, double amount) {
    vector<string> t;
    for(auto i : extrapolations) if(i.first != "noise3") t.push_back(i.first);
    auto v1 = extrapolations[t[rdint(0, t.size()-1)]](n, amount);
    auto v2 = extrapolations[t[rdint(0, t.size()-1)]](n, amount);
    auto n1 = noise_scale*pn.GetNoise(v1.x*coord_scale, v1.y*coord_scale);
    auto n2 = noise_scale*pn.GetNoise(v2.x*coord_scale, v2.y*coord_scale);
    return sf::Vector2f(amount*n1, amount*n2);
}
double length(sf::Vector2f v) {
    return noise_scale*sqrt(v.x*v.x*coord_scale + v.y*v.y*coord_scale);
}
#endif



sf::Vector2f cartesianXY(double n, double amount) {
    double offx = rd(0, 100) < 50 ? -1 : 1, offy = rd(0, 100) < 50 ? -1 : 1;
    return sf::Vector2f(n*amount*offx, n*amount*offy);
}

sf::Vector2f polarXY(double n, double amount) {
    double xt = cos(M_PI*2*n);
    double yt = sin(M_PI*2*n); 
    return sf::Vector2f(amount*xt, amount*yt);
}

double arctangent(sf::Vector2f v) {
    return atan2(v.x, v.y)*10;
}

double discard_chance;

double discard1(sf::Vector2f v) {
    return discard_chance < 50 ? v.x : v.y;
}

sf::Vector2f discard2(sf::Vector2f v, double amount) {
    double xx = v.x, yy = v.y;
    if(discard_chance < 50) xx = yy;
    else yy = xx;
    return sf::Vector2f(amount*xx, amount*yy);
}

double dotproduct(sf::Vector2f v) {
    vector<double> v1 = {v.x, v.y}, v2 = {v.y, v.x};
    return inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

sf::Vector2f copy(double n, double amount) {
    return sf::Vector2f(n*amount, n*amount);
}

sf::Vector2f archispiral(double n, double amount) {
    double xt = n*cos(n), yt = n*sin(n);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f cycloid(double n, double amount) {
    double xt = cos(n)*2*cos(2*n)+1, yt = sin(n)*2*sin(2*n)+1;
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f deltoid(double n, double amount) {
    double xt = (2*cos(n))/3+(cos(2*n)/3), yt = (2*sin(n))/3-(sin(2*n)/3);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f folium(double n, double amount) {
    double xt = (3*n)/(pow(n, 3)+1), yt = (3*n*n)/(pow(n, 3)+1);
    return sf::Vector2f(amount*xt, amount*yt);
}

double trochoid_b, trochoid_h;

sf::Vector2f epitrochoid(double n, double amount) {
    double a = 16, b = trochoid_b, h = trochoid_h;
    double xt = (a+b)*cos(n)-h*cos((n*(a+b))/b);
    double yt = (a+b)*sin(n)-h*sin((n*(a+b))/b);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f hypotrochoid(double n, double amount) {
    double a = 16, b = trochoid_b, h = trochoid_h;
    double xt = (a-b)*cos(n)+h*cos((n*(a-b))/b);
    double yt = (a-b)*sin(n)-h*sin((n*(a-b))/b);
    return sf::Vector2f(amount*xt, amount*yt);
}

double sinc(double x) { return sin(x)/x; }

sf::Vector2f quadratrix(double n, double amount) {
    double xt = (2*cos(n))/(M_PI*sinc(n));
    double yt = (2*sin(n))/(M_PI*sinc(n));
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f kilroy(double n, double amount) {
    double xt = n;
    double yt = log(abs(sinc(n)));
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f lituus(double n, double amount) {
    double xt = cos(n)/sqrt(n);
    double yt = sin(n)/sqrt(n);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f maltesecross(double n, double amount) {
    double xt = (2*cos(n))/(sqrt(sin(4*n)));
    double yt = (2*sin(n))/(sqrt(sin(4*n)));
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f neoid(double n, double amount) {
    double xt = cos(n)*(n+1);
    double yt = sin(n)*(n+1);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f sextic(double n, double amount) {
    double xt = pow(cos(n/3), 3)*cos(n);
    double yt = sin(n)*pow(cos(n/3), 3);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f eightcurve(double n, double amount) {
    double xt = sin(n);
    double yt = sin(n)*cos(n);
    return sf::Vector2f(amount*xt, amount*yt);
}

double ecycloid_a, ecycloid_b, hcycloid_a, hcycloid_b;

sf::Vector2f epicycloid(double n, double amount) {
    double a = ecycloid_a, b = ecycloid_b;
    double xt = (a+b)*cos(n)-b*cos(n*(a+b)/b);
    double yt = (a+b)*sin(n)-b*sin(n*(a+b)/b);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f hypocycloid(double n, double amount) {
    double a = hcycloid_a, b = hcycloid_b;
    double xt = (a-b)*cos(n)+b*cos(n*(a-b)/b);
    double yt = (a-b)*sin(n)-b*sin(n*(a-b)/b);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f trisectrix(double n, double amount) {
    double xt = (n*n-3)/(n*n+1);
    double yt = n*(n*n-3)/(n*n+1);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f quadrifolium(double n, double amount) {
    double xt = sin(2*n)*cos(n);
    double yt = sin(n)*sin(2*n);
    return sf::Vector2f(amount*xt, amount*yt);
}

sf::Vector2f tschirnhausen(double n, double amount) {
    double xt = 1-3*n*n;
    double yt = n*(3-n*n);
    return sf::Vector2f(amount*xt, amount*yt);
}

double xpowy(sf::Vector2f v) {
    return pow(v.x+3, v.y+3);
}

double xmody(sf::Vector2f v) {
    return fmod(v.x, v.y);
}

double xdivy(sf::Vector2f v) {
    return v.y == 0 ? 0 : v.x/v.y;
}

double xmuly(sf::Vector2f v) {
    return v.x*v.y;
}

double sqrtfix(double n) {
    return sqrt(n+3);
}

sf::Vector2f cornoid(double n, double amount) {
    auto xx = (1-2*sin(n)*sin(n))*cos(n);
    auto yy = sin(n)*(2*cos(n)*cos(n)+1);
    return sf::Vector2f(xx*amount, yy*amount);
}

sf::Vector2f fishcurve(double n, double amount) {
    auto xx = cos(n)-(sin(n)*sin(n))/sqrt(2);
    auto yy = sin(n)*cos(n);
    return sf::Vector2f(xx*amount, yy*amount);
}

sf::Vector2f ranunculoid(double n, double amount) {
    auto xx = 6*cos(n)-cos(6*n);
    auto yy = 6*sin(n)-sin(6*n);
    return sf::Vector2f(xx*amount, yy*amount);
}

void initialise_variations() {
#ifdef NOISE
    transformations["noise"] = noise;
    truncations["perlin"] = perlin;
    //extrapolations["noise3"] = noise3;
    truncations["length"] = length;
#endif
    discard_chance = rd(0, 100);
    trochoid_b = rd(1, M_PI);
    trochoid_h = rd(1, 6);
    A = rd(-5, 5);
    B = rd(-5, 5);
    C = rdnormal(0.5, 0.2)*(rd(0, 100) < 50 ? 1 : -1);
    D = rd(-5, 5);
    E = rd(-5, 5);
    F = rd(-5, 5);
    fan_t = M_PI*C*C;
    curl_p1 = rd(-5, 5), curl_p2 = rd(-5, 5);
    ecycloid_a = rd(1, 4), ecycloid_b = rd(1, 3);
    hcycloid_a = rd(sqrt(2), 5), hcycloid_b = rd(1/3, 3);

    transformations["sin"] = sin;
    transformations["cos"] = cos;
    transformations["tan"] = tan;
    transformations["exp"] = exp;
    transformations["log"] = log;
    transformations["sqrt"] = sqrtfix;

    truncations["arctangent"] = arctangent;
    truncations["xpowy"] = xpowy;
    truncations["xmody"] = xmody;
    truncations["xdivy"] = xdivy;
    truncations["xmuly"] = xmuly;

    extrapolations["ranunculoid"] = ranunculoid;
    extrapolations["fishcurve"] = fishcurve;
    extrapolations["cornoid"] = cornoid;
    extrapolations["neoid"] = neoid;
    //extrapolations["maltesecross"] = maltesecross;
    extrapolations["lituus"] = lituus;
    //extrapolations["kilroy"] = kilroy;
    extrapolations["quadratrix"] = quadratrix;
    extrapolations["epitrochoid"] = epitrochoid;
    extrapolations["folium"] = folium;
    extrapolations["cycloid"] = cycloid;
    extrapolations["archispiral"] = archispiral;
    extrapolations["cissoid"] = cissoid;
    extrapolations["kampyle"] = kampyle;
    extrapolations["astroid"] = astroid;
    extrapolations["superformula"] = superformula;
    extrapolations["cartesianXY"] = cartesianXY;
    extrapolations["polarXY"] = polarXY;
    extrapolations["hyperbole"] = hyperbole;
    extrapolations["copy"] = copy;
    extrapolations["sextic"] = sextic;
    extrapolations["eightcurve"] = eightcurve;
    extrapolations["epicycloid"] = epicycloid;
    extrapolations["hypocycloid"] = hypocycloid;
    extrapolations["hypotrochoid"] = hypotrochoid;
    extrapolations["trisectrix"] = trisectrix;
    extrapolations["quadrifolium"] = quadrifolium;
    extrapolations["tschirnhausen"] = tschirnhausen;

    variations["rectangles"] = rectangles;
    variations["curl"] = curl;
    //variations["gaussian"] = gaussian;
    variations["fan"] = fan;
    variations["rings"] = rings;
    variations["cosine"] = cosine;
    variations["julia"] = julia;
    variations["sinusoidal"] = sinusoidal;
    variations["sech"] = sech;
    variations["hyperbolic"] = hyperbolic;
    variations["pdj"] = pdj;
    variations["spherical"] = spherical;
    variations["linear"] = linear;
    variations["swirl"] = swirl;
    variations["polar"] = polar;
    variations["heart"] = heart;
    variations["disc"] = disc;
    variations["spiral"] = spiral;
    variations["ex"] = ex;
    variations["bent"] = bent;
    variations["waves"] = waves;
    variations["fisheye"] = fisheye;
    variations["popcorn"] = popcorn;
    variations["exponential"] = exponential;
    variations["horseshoe"] = horseshoe;
    variations["handkerchief"] = handkerchief;
    variations["diamond"] = diamond;
    variations["power"] = power;
    variations["bubble"] = bubble;
    variations["tangent"] = tangent;
    //variations["square"] = square;
    variations["cross"] = cross;
}

struct node {
    int type;
    char operation;
    string var;
    node() {};
    node(int t, char o, string v1) {
        type = t;
        operation = o;
        var = v1;
    }
};

vector<int> tree[1000];
node nodes[1000];

sf::Vector2f resolveFoldTree(int u, sf::Vector2f v, double amount) {
    node U = nodes[u];
    if(U.type == 1) return variations[U.var](v, amount);
    else if(U.type == 3) {
        if(tree[u].size()) return variations[U.var](resolveFoldTree(tree[u][0], v, amount), amount);
        else return variations[U.var](v, amount);
    }
    else {
        switch(U.operation) {
            case '+':
                return addF(resolveFoldTree(tree[u][0], v, amount), resolveFoldTree(tree[u][1], v, amount));
                break;
            case '-':
                return subF(resolveFoldTree(tree[u][0], v, amount), resolveFoldTree(tree[u][1], v, amount));
                break;
            case '*':
                return mulF(resolveFoldTree(tree[u][0], v, amount), resolveFoldTree(tree[u][1], v, amount));
                break;
            case '/':
                return divF(resolveFoldTree(tree[u][0], v, amount), resolveFoldTree(tree[u][1], v, amount));
                break;
        }
    }
}

void printFoldTree(int u) {
    node U = nodes[u];
    if(U.type == 1) cout << U.var;
    else if(U.type == 3) {
        cout << U.var;
        cout << "(";
        printFoldTree(tree[u][0]);
        cout << ")";
    }
    else {
        switch(U.operation) {
            case '+':
                cout << "(";
                printFoldTree(tree[u][0]);
                cout << ") + (";
                printFoldTree(tree[u][1]);
                cout << ")";
                break;
            case '-':
                cout << "(";
                printFoldTree(tree[u][0]);
                cout << ") - (";
                printFoldTree(tree[u][1]);
                cout << ")";
                break;
            case '*':
                cout << "(";
                printFoldTree(tree[u][0]);
                cout << ") * (";
                printFoldTree(tree[u][1]);
                cout << ")";
                break;
            case '/':
                cout << "(";
                printFoldTree(tree[u][0]);
                cout << ") / (";
                printFoldTree(tree[u][1]);
                cout << ")";
                break;
        }
    }
}

int ncount = 1, max_depth, max_ncount = 4;

void addFoldNode(int u, int depth) {
    //cout << depth << endl;
    if(depth > max_depth || nodes[u].type == 1) return;
    for(int i = 1; i <= (nodes[u].type == 2 ? 2 : 1); i++) {
        tree[u].push_back(++ncount);
        nodes[ncount] = node((depth == max_depth || ncount >= max_ncount) ? 1 : rdint(1, 3), randOperation(), randVariation());
        addFoldNode(ncount, depth+1);
    }
}

void createFoldTree() {
    node n((max_ncount == 1? 1 : 2), randOperation(),randVariation());
    //cout << n.type << endl;
    nodes[ncount] = n;
    max_depth = 3;//rd(2, 4);
    addFoldNode(1, 0);
}

string randTruncation() {
    vector<string> t;
    for(auto i : truncations)
        t.push_back(i.first);
    return t[rdint(0, t.size()-1)];
}

string randTransformation() {
    vector<string> t;
    for(auto i : transformations)
        t.push_back(i.first);
    return t[rdint(0, t.size()-1)];
}

void addFieldNode(int u, int depth) {
    if(depth >= max_depth) return;
    tree[u].push_back(++ncount);
    if(nodes[u].type == 1 || nodes[u].type == 0) {
        if(rd(0, 100) < 70 || depth == max_depth-1) nodes[ncount] = node(3, '0', randExtrapolation());
        else nodes[ncount] = node(0, '0', randTransformation());
    } else {
        int type = 2;
        if((depth != max_depth-1 && rd(0, 100) < 50)||(rd(0, 100) < 50 && depth == max_depth-2)) type = 1;
        nodes[ncount] = node(type, '0', type == 1 ? randTruncation() : randVariation());
    }
    addFieldNode(ncount, depth+1);
}

void createFieldTree() {
    int type = rd(0, 100) < 10 ? 1 : 2;
    node n(type, '0', type == 1 ? randTruncation() : randVariation());
    nodes[ncount] = n;
    max_depth = rd(2, 4);
    addFieldNode(1, 0);
}

void printFieldTree(int u, string &s) {
    s += nodes[u].var;
    if(tree[u].size()) s += " → ", printFieldTree(tree[u][0], s);
}

sf::Vector2f resolveFieldTree(int u, sf::Vector2f v) {
    if(nodes[u].type == 0) {
        auto res = transformations[nodes[u].var](v.x);
        v.x = res;
        return resolveFieldTree(tree[u][0], v);
    } else if(nodes[u].type == 1) {
        auto res = truncations[nodes[u].var](v);
        v.x = res;
        return resolveFieldTree(tree[u][0], v);
    } else if(nodes[u].type == 2) {
        auto res = variations[nodes[u].var](v, 1);
        if(tree[u].size()) return resolveFieldTree(tree[u][0], res);
        else return res;
    } else if(nodes[u].type == 3) {
        auto res = extrapolations[nodes[u].var](v.x, 1);
        if(tree[u].size()) return resolveFieldTree(tree[u][0], res);
        else return res;
    }
}

bool stringToFieldTree(string s) {
    string current = "";
    int prevtype = 2;
    for(int i = 0; i < s.size(); i++) {
        while(s[i] != '>' && i < s.size()) current += s[i], i++;
        int type = -1;
        if(transformations.count(current) && (prevtype == 1 || prevtype == 0) && i < s.size()) type = 0;
        else if(truncations.count(current) && (prevtype == 2 || prevtype == 3) && i < s.size()) type = 1;
        else if(variations.count(current) && (prevtype == 2 || prevtype == 3)) type = 2;
        else if(extrapolations.count(current) && (prevtype == 0 || prevtype == 1)) type = 3;
        else return 0;
        //cout << current << ' ' << type << endl;
        nodes[ncount] = node(type, '0', current);
        if(i < s.size()) tree[ncount].push_back(ncount+1);
        ncount++;
        current = "";
        prevtype = type;
    }
    return 1;
}

set<int> freenodes;

sf::Vector2f resolveFieldTree2(int u, sf::Vector2f v, double amount, vector<vector<int> > &tree, vector<node> &nodes) { //0 - trunc, 1- var, 2 - ext, 3 - trans, 4 - op
    node U = nodes[u];
    if(!tree[u].size()) {
        if(U.type == 0) return sf::Vector2f(truncations[U.var](v), 0);
        else return variations[U.var](v, amount);
    } else if(U.type < 4) {
        switch(U.type) {
            case 0:
                return sf::Vector2f(truncations[U.var](resolveFieldTree2(tree[u][0], v, amount, tree, nodes)), 0);
                break;
            case 1:
                return variations[U.var](resolveFieldTree2(tree[u][0], v, amount, tree, nodes), amount);
                break;
            case 2:
                return extrapolations[U.var](resolveFieldTree2(tree[u][0], v, amount, tree, nodes).x, amount);
                break;
            case 3:
                return sf::Vector2f(transformations[U.var](resolveFieldTree2(tree[u][0], v, amount, tree, nodes).x), 0);
                break;
        }
    }
    else {
        switch(U.operation) {
            case '+':
                return addF(resolveFieldTree2(tree[u][0], v, amount, tree, nodes), resolveFieldTree2(tree[u][1], v, amount, tree, nodes));
                break;
            case '-':
                return subF(resolveFieldTree2(tree[u][0], v, amount, tree, nodes), resolveFieldTree2(tree[u][1], v, amount, tree, nodes));
                break;
            case '*':
                return mulF(resolveFieldTree2(tree[u][0], v, amount, tree, nodes), resolveFieldTree2(tree[u][1], v, amount, tree, nodes));
                break;
            case '/':
                return divF(resolveFieldTree2(tree[u][0], v, amount, tree, nodes), resolveFieldTree2(tree[u][1], v, amount, tree, nodes));
                break;
        }
    }
}

string printFieldTree2(int u, string &s, vector<vector<int> > &tree, vector<node> &nodes) {
    node U = nodes[u];
    if(!tree[u].size()) s += U.var;
    else if(U.type < 4) {
        s += U.var;
        s += "(";
        printFieldTree2(tree[u][0], s, tree, nodes);
        s += ")";
    } else {
        switch(U.operation) {
            case '+':
                s += "(";
                printFieldTree2(tree[u][0], s, tree, nodes);
                s += ") + (";
                printFieldTree2(tree[u][1], s, tree, nodes);
                s += ")";
                break;
            case '-':
                s += "(";
                printFieldTree2(tree[u][0], s, tree, nodes);
                s += ") - (";
                printFieldTree2(tree[u][1], s, tree, nodes);
                s += ")";
                break;
            case '*':
                s += "(";
                printFieldTree2(tree[u][0], s, tree, nodes);
                s += ") * (";
                printFieldTree2(tree[u][1], s, tree, nodes);
                s += ")";
                break;
            case '/':
                s += "(";
                printFieldTree2(tree[u][0], s, tree, nodes);
                s += ") / (";
                printFieldTree2(tree[u][1], s, tree, nodes);
                s += ")";
                break;
        }
    }
    return s;
}

int randFromVector(vector<int> v) {
    return v[rdint(0, v.size()-1)];
}

int varchance = 20;

void addFieldNode2(int u, int depth, vector<vector<int> > &tree, vector<node> &nodes) { //0 - trunc, 1- var, 2 - ext, 3 - trans, 4 - op
    //cout << depth << endl;
    if(depth > max_depth || (nodes[u].type <= 1 && (rd(0, 100) < mathmap(depth, 0, max_depth, 0, 100)))) return;
    for(int i = 1; i <= (nodes[u].type == 4 ? 2 : 1); i++) {
        tree[u].push_back(++ncount);
        int ntype;
        if(depth >= max_depth || ncount >= max_ncount) {
            if(nodes[u].type <= 1 || nodes[u].type == 4) ntype = ((rd(0, 100) < varchance || nodes[u].type == 4) ? 1 : 2);
            else ntype = 0;
        } else {
            if(nodes[u].type <= 1 || nodes[u].type == 4) ntype = (rd(0, 100) < 20 ? ((rd(0, 100) < varchance || nodes[u].type == 4) ? 1 : 2) : 4);
            else ntype = randFromVector({0, 3});
        }
        nodes[ncount] = node(ntype, randOperation(), ntype == 0 ? (depth == max_depth || ncount >= max_ncount) ? "perlin" : randTruncation() : ntype == 1 ? randVariation() : ntype == 2 ? randExtrapolation() : randTransformation());
        addFieldNode2(ncount, depth+1, tree, nodes);
    }
}

void createFieldTree2(int maxcount, int maxdepth, vector<vector<int> > &tree, vector<node> &nodes, int vc) {
    varchance = vc;
    ncount = 1;
    for(int i = 0; i < 100; i++) tree.push_back({}), nodes.push_back(node());
    max_ncount = maxcount;
    max_depth = maxdepth;
    int type = (maxdepth == -1 ? 1 : randFromVector({1, 4}));
    node n(type, randOperation(), type == 1 ? randVariation() : randExtrapolation());
    nodes[ncount] = n;
    addFieldNode2(1, 0, tree, nodes);
}