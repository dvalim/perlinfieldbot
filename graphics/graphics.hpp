#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include "../math/math.hpp"

#pragma once

namespace aural {
    sf::RenderTexture * renderer;

    struct Rgb {
        double r, g, b, a;
        Rgb(double r, double g, double b, double a) : r(r), g(g), b(b), a(a) {}
        Rgb(double r, double g, double b) : r(r), g(g), b(b), a(255) {}
        Rgb() : r(0), g(0), b(0), a(255) {}
    };

    sf::Color convert(Rgb RGB) {
        sf::Color color(RGB.r, RGB.g, RGB.b, RGB.a);
        return color;
    }

    struct Hsv {
        double h, s, v;
        Hsv(double h, double s, double v) : h(h), s(s), v(v) {}
        Hsv() : h(0), s(0), v(0) {}
    }; 

    Rgb hsv2rgb(Hsv in) {
        double      hh, p, q, t, ff;
        long        i;
        Rgb         out;

        if(in.s <= 0.0) {
            out.r = in.v;
            out.g = in.v;
            out.b = in.v;
            return out;
        }

        hh = in.h;
        if(hh >= 360.0) hh = 0.0;
        hh /= 60.0;
        i = (long)hh;
        ff = hh - i;
        p = in.v * (1.0 - in.s);
        q = in.v * (1.0 - (in.s * ff));
        t = in.v * (1.0 - (in.s * (1.0 - ff)));

        switch(i) {
        case 0:
            out.r = in.v;
            out.g = t;
            out.b = p;
            break;
        case 1:
            out.r = q;
            out.g = in.v;
            out.b = p;
            break;
        case 2:
            out.r = p;
            out.g = in.v;
            out.b = t;
            break;

        case 3:
            out.r = p;
            out.g = q;
            out.b = in.v;
            break;
        case 4:
            out.r = t;
            out.g = p;
            out.b = in.v;
            break;
        //case 5:
        //	break;
        default:
            out.r = in.v;
            out.g = p;
            out.b = q;
            break;
        }
        out.r *= 255, out.g *= 255, out.b *= 255;
        return out;     
    }

    Hsv rgb2hsv(Rgb in) {
        in.r /= 255, in.g /= 255, in.b /= 255;
        Hsv         out;
        double      min, max, delta;

        min = in.r < in.g ? in.r : in.g;
        min = min  < in.b ? min  : in.b;

        max = in.r > in.g ? in.r : in.g;
        max = max  > in.b ? max  : in.b;

        out.v = max;                                // v
        delta = max - min;
        if (delta < 0.00001)
        {
            out.s = 0;
            out.h = 0; // undefined, maybe nan?
            return out;
        }
        if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
            out.s = (delta / max);                  // s
        } else {
            // if max is 0, then r = g = b = 0              
            // s = 0, h is undefined
            out.s = 0.0;
            out.h = NAN;                            // its now undefined
            return out;
        }
        if( in.r >= max )                           // > is bogus, just keeps compilor happy
            out.h = ( in.g - in.b ) / delta;        // between yellow & magenta
        else
        if( in.g >= max )
            out.h = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
        else
            out.h = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

        out.h *= 60.0;                              // degrees

        if( out.h < 0.0 )
            out.h += 360.0;

        return out;
    }

    sf::Color convert(Hsv HSV) {
        return convert(hsv2rgb(HSV));
    }

    void clear(sf::Color col) {
        renderer->clear(col);
    }

    sf::Color interpolate(sf::Color color1, sf::Color color2, double percent) {
        double resultRed = color1.r + percent * (color2.r - color1.r);
        double resultGreen = color1.g + percent * (color2.g - color1.g);
        double resultBlue = color1.b + percent * (color2.b - color1.b);
        double resultAlpha = color1.a + percent * (color2.a - color1.a);
        return sf::Color(resultRed, resultGreen, resultBlue, resultAlpha);
    }

    void rect(double x, double y, double w, double h) {
        sf::RectangleShape rectangle(sf::Vector2f(w, h));
        rectangle.setPosition(x, y);
        rectangle.setFillColor(sf::Color::Black);
        renderer->draw(rectangle);
    }

    void rect(double x, double y, double w, double h, sf::Color col) {
        sf::RectangleShape rectangle(sf::Vector2f(w, h));
        rectangle.setPosition(x, y);
        rectangle.setFillColor(col);
        renderer->draw(rectangle);
    }

    void rect(double x, double y, double w, double h, sf::Color col, sf::BlendMode blend) {
        sf::RectangleShape rectangle(sf::Vector2f(w, h));
        rectangle.setPosition(x, y);
        rectangle.setFillColor(col);
        renderer->draw(rectangle, blend);
    }

    void display() {
        renderer->display();
    }

    class Particle {
        public:
            math::vec pos, vel;
            double angle;
            sf::Color color;
            Particle() : pos(math::vec(0, 0)), vel(math::vec(0, 0)), angle(0), color(sf::Color::White) {}
            Particle(math::vec pos, sf::Color c) : pos(pos), vel(math::vec(0, 0)), angle(0), color(c) {}
            void draw() {
                rect(pos.x, pos.y, 1, 1, color);
            }
            void move() {
                pos = pos + vel;
            }
    };
}

#undef sf::Color