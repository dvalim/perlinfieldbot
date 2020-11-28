#define WINDOW
#define OUTPUT

#include "aural.hpp"
#include <iostream>

using namespace au;
using namespace au::app;
using namespace au::math;
using namespace au::rnd;
using namespace au::flame;

double ssin(double x, double p) {
    if(p == 0) return asin(sin(x));
    else return pow(sin(x), p);
}

double ccos(double x, double p) {
    if(p == 0) return acos(sin(x));
    else return pow(cos(x), p);
}

class baseApp : public App {
    using App::App;
    public:
        double a[7], f[7], v, p, q, x = 0, y = 0, t = 0;
        sf::Color colors[3];
        double minx = 10, miny = 10, maxx = -10, maxy = -10, color_limit = 1;
        bool colored = 0;

        struct attractor {
            double a[7], f[7], v, p, q, x = 0, y = 0, t = 0;
            double centerx, centery, radius;
            sf::Color colors[3];
            double minx = 10, miny = 10, maxx = -10, maxy = -10, color_limit = 1;

            void create() {
                for(int i = 1; i <= 6; i++)
                    a[i] = random(0.7, 1.2) * (random(0, 100) > 50 ? -1: 1),
                    f[i] = random(0.7, 1.2) * (random(0, 100) > 50 ? -1: 1);
                v = random(0.001, 0.5);
                p = random(1, 5);
                q = random(1, 5);
            }

            void init(bool colored) {
                for(double x = -1; x <= 1; x += 0.02)
                    for(double y = -1; y <= 1; y += 0.02) {
                        double xx = x*sqrt(1-y*y/2);
                        double yy = y*sqrt(1-x*x/2);
                        sf::Color color1 = interpolate(colors[0], colors[1], map(x, -1, 1, 0, 1));
                        sf::Color color2 = interpolate(colors[0], colors[1], map(y, -1, 1, 0, 1));
                        sf::Color color3 = interpolate(color1, color2, 0.5);
                        color3.a = 80 * map(radius, 80, 800, 0.2, 1.2);
                        vec f = vec(xx, yy) - vec(0, 0);
                        color3.a *= f.mag2();
                        rect(xx*radius+centerx, yy*radius+centery, 2, 2, color3/*(colored ? sf::Color(255, 255, 255, 60) : sf::Color(0, 0, 0, 60))*/);
                    }
            }

            void iterate(int screen_width, int screen_height, bool colored) {
                for(int i = 1; i <= 10000*map(radius, 80, 800, 0.2, 1); i++) {
                    double xx = x, yy = y;
                    xx = a[1]*ssin(f[1]*x, p) + a[2]*ccos(f[2]*y, q) + a[3]*ssin(f[3]*t, p);
                    yy = a[4]*ccos(f[4]*x, q) + a[5]*ssin(f[5]*y, p) + a[6]*ssin(f[6]*t, q);
                    t += v;

                    vec step(xx-x, yy-y);
                    step = step*10;
                    double tt = constrain(map(step.mag2(), 0, M_PI*M_PI*2, 0, 1), -color_limit, color_limit);
                    sf::Color color = interpolate(colors[0], colors[1], tt);
                    color.a = 10 * gaussian(1, 0.1);
                    //std::cout << step.mag2()*pow(10, 5) << ' ' << tt << "\n\n";
                    x = xx;
                    y = yy;

                    minx = std::min(minx, x);
                    miny = std::min(miny, y);
                    maxx = std::max(maxx, x);
                    maxy = std::max(maxy, y);

                    xx = map(x, minx, maxx, -1, 1);
                    yy = map(y, miny, maxy, -1, 1);

                    //xx -= screen_width/2;
                    //yy -= screen_height/2;

                    double xxx = xx * sqrt(1-(yy*yy)/2);//+screen_width/2;
                    double yyy = yy * sqrt(1-(xx*xx)/2);//+screen_height/2;

                    xxx *= radius;
                    yyy *= radius;

                    xxx += centerx;
                    yyy += centery;

                    vec f = vec(xxx, yyy) - vec(centerx, centery);
                    color.a *= map(f.mag(), 0, radius, 0, 1);

                    if(i > 2000) {
                        if(colored) rect(xxx, yyy, 1, 1, color, sf::BlendAdd);
                        else rect(xxx, yyy, 1, 1, color);
                    }
                }
            }

            attractor(double cx, double cy, double r) : centerx(cx), centery(cy), radius(r) {}
            attractor() {

            }
        };

        std::vector<attractor> attractors;
        int attractor_number = 10, timer = 0;
        attractor bg;

        virtual void init() {
            #ifdef WINDOW
            initWindow();
            #endif
            texture.create(screen_width, screen_height);
            texture.setSmooth(true);
            clock.restart();
            seed = std::chrono::system_clock::now().time_since_epoch().count();
            //seed = 1545768988279824;
        }

        void setup() {
            #ifdef WINDOW
            sf::ContextSettings settings;
            settings.antialiasingLevel = 8;
            window.create(sf::VideoMode(screen_width, screen_height), window_name, sf::Style::Close, settings);
            window.setVerticalSyncEnabled(true);
            window.setFramerateLimit(60);
            #endif

            timeLimit = 30;
            save_path = "./";

            std::cout << seed << '\n';

            /*for(int i = 1; i <= 6; i++)
                std::cout << "a" << i << " = " << a[i] << '\n';
            for(int i = 1; i <= 6; i++)
                std::cout << "f" << i << " = " << f[i] << '\n';

            std::cout << "v = " << v << "\np = " << p << "\nq = " << q << '\n';*/


            double sum_angle = 0;
            for(int i = 1; i <= attractor_number; i++) {
                double radius = 600 * map(i, 1, attractor_number, 1, 0.125);//(random(0, 100) < 50 ? random(500.0, 800.0) : random(60.0, 250.0));

                double angle = fmod(sum_angle + random(M_PI*0.8, M_PI*1.2), 2*M_PI);
                if(sum_angle) {
                    sum_angle += angle;
                    sum_angle /= 2;
                } else {
                    sum_angle += angle;
                }
                double center_radius = constrain(gaussian(800, 400), 0, 1000);
                double centerx = cos(angle)*center_radius+screen_width/2;
                double centery = sin(angle)*center_radius+screen_height/2;

                attractors.push_back(attractor(centerx, centery, radius));
                attractors.back().create();
            }

            bg = attractor(screen_width/2, screen_height/2, 1200);
            bg.create();

            int hue = random(0, 360);
            color_limit = random(1.0, 7.0);
            if(random(0, 100) < 70) {
                colored = 1;
                sf::Vertex background[] = {
                    sf::Vertex(sf::Vector2f(0, 0), convert(Hsv((hue+random(100, 260))%360, 0.4, 1.4*random(0.03, 0.04)))),
                    sf::Vertex(sf::Vector2f(screen_width, 0), convert(Hsv((hue+random(100, 260))%360, 0.4, 1.4*random(0.03, 0.04)))),
                    sf::Vertex(sf::Vector2f(screen_width, screen_height), convert(Hsv((hue+random(100, 260))%360, 0.4, 1.4*random(0.03, 0.04)))),
                    sf::Vertex(sf::Vector2f(0, screen_height), convert(Hsv((hue+random(100, 260))%360, 0.4, 1.4*random(0.03, 0.04))))
                };

                texture.draw(background, 4, sf::Quads);

                for(int i = 0; i < attractor_number; i++)
                    attractors[i].colors[0] = convert(Hsv((hue+random(0, 20))%360, 0.5, 0.8)),
                    attractors[i].colors[1] = convert(Hsv((hue+random(100, 260))%360, 0.5, 0.8)),
                    attractors[i].colors[0].a = 10,
                    attractors[i].colors[1].a = 10;

                bg.colors[0] = convert(Hsv((hue+random(0, 20))%360, 0.5, 0.6)),
                bg.colors[1] = convert(Hsv((hue+random(100, 260))%360, 0.5, 0.6));
            } else {
                texture.clear(sf::Color(200, 200, 200));
            
                for(int i = 0; i < attractor_number; i++)
                    attractors[i].colors[0] = convert(Hsv((hue+random(0, 20))%360, 0.4, 0.2)),
                    attractors[i].colors[1] = convert(Hsv((hue+random(100, 260))%360, 0.4, 0.2)),
                    attractors[i].colors[0].a = 10,
                    attractors[i].colors[1].a = 10;

                bg.colors[0] = convert(Hsv((hue+random(0, 20))%360, 0.4, 0.1)),
                bg.colors[1] = convert(Hsv((hue+random(100, 260))%360, 0.4, 0.1));
                bg.colors[0].a = 10;
                bg.colors[1].a = 10;
            }

            for(auto &a : attractors)
                a.init(colored);
            
            //std::cout << seed << '\n';
        }
        void loop() {
            #ifdef WINDOW
            checkForEvents();
            #endif

            //std::cout << x << ' ' << y << '\n';

            for(auto &a : attractors)
                a.iterate(screen_width, screen_height, colored);

            if(timer++ < 50) bg.iterate(screen_width, screen_height, colored);

            #ifdef WINDOW
            drawTextureToWindow();
            #endif
        }
        void close() {
            #ifdef WINDOW
            if(window.isOpen()) window.close();
            #endif
             #ifdef OUTPUT
            std::string info = saveFile() + "#method=fujii attractor";
            info+="#seed="+std::to_string(seed);
            std::cout << info;
            #else
            saveFile();
            #endif
        }
} testApp(2048, 2048);

AURAL_APP(testApp)