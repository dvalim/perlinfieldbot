#define WINDOW
#define OUTPUT

#include "aural.hpp"
#include <iostream>

using namespace au;
using namespace au::app;
using namespace au::math;
using namespace au::rnd;
using namespace au::flame;

class baseApp : public App {
    using App::App;
    public:
        double a[7], f[7], v, p, q, x = 0, y = 0, t = 0;
        sf::Color colors[3];
        double minx = 10, miny = 10, maxx = -10, maxy = -10, color_limit = 1;
        bool colored = 0;

        double ssin(double x, double p) {
            if(p == 0) return asin(sin(x));
            else return pow(sin(x), p);
        }

        double ccos(double x, double p) {
            if(p == 0) return acos(sin(x));
            else return pow(cos(x), p);
        }

        virtual void init() {
            #ifdef WINDOW
            initWindow();
            #endif
            texture.create(screen_width, screen_height);
            texture.setSmooth(true);
            clock.restart();
            seed = std::chrono::system_clock::now().time_since_epoch().count();
            //seed = 1542315600081299;
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

            for(int i = 1; i <= 6; i++)
                a[i] = random(0.7, 1.2) * (random(0, 100) > 50 ? -1: 1),
                f[i] = random(0.7, 1.2) * (random(0, 100) > 50 ? -1: 1);
            v = random(0.001, 0.5);
            p = random(1, 5);
            q = random(1, 5);

            /*for(int i = 1; i <= 6; i++)
                std::cout << "a" << i << " = " << a[i] << '\n';
            for(int i = 1; i <= 6; i++)
                std::cout << "f" << i << " = " << f[i] << '\n';

            std::cout << "v = " << v << "\np = " << p << "\nq = " << q << '\n';*/

            int hue = random(0, 360);
            color_limit = random(1.0, 7.0);
            if(random(0, 100) < 70) {
                colored = 1;
                sf::Vertex background[] = {
                    sf::Vertex(sf::Vector2f(0, 0), convert(Hsv((hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06)))),
                    sf::Vertex(sf::Vector2f(screen_width, 0), convert(Hsv((hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06)))),
                    sf::Vertex(sf::Vector2f(screen_width, screen_height), convert(Hsv((hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06)))),
                    sf::Vertex(sf::Vector2f(0, screen_height), convert(Hsv((hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06))))
                };

                texture.draw(background, 4, sf::Quads);

                colors[0] = convert(Hsv(hue, 0.5, 0.8));
                colors[1] = convert(Hsv((hue+random(100, 260))%360, 0.5, 0.8));
            } else {
                texture.clear(sf::Color(200, 200, 200));
            
                colors[0] = convert(Hsv(hue, 0.4, 0.2));
                colors[1] = convert(Hsv((hue+random(100, 260))%360, 0.4, 0.2));
            }
            
            //std::cout << seed << '\n';
        }
        void loop() {
            #ifdef WINDOW
            checkForEvents();
            #endif

            for(int i = 1; i <= 20000; i++) {
                double xx = x, yy = y;
                xx = a[1]*ssin(f[1]*x, p) + a[2]*ccos(f[2]*y, q) + a[3]*ssin(f[3]*t, p);
                yy = a[4]*ccos(f[4]*x, q) + a[5]*ssin(f[5]*y, p) + a[6]*ssin(f[6]*t, q);
                t += v;

                vec step(xx-x, yy-y);
                step = step*10;
                double tt = constrain(map(step.mag2(), 0, M_PI*M_PI*2, 0, 1), -color_limit, color_limit);
                sf::Color color = interpolate(colors[0], colors[1], tt);
                color.a = (colored ? 10 : 15);
                //std::cout << step.mag2()*pow(10, 5) << ' ' << tt << "\n\n";
                x = xx;
                y = yy;

                minx = std::min(minx, x);
                miny = std::min(miny, y);
                maxx = std::max(maxx, x);
                maxy = std::max(maxy, y);

                xx = map(x, minx, maxx, 50, screen_width-50);
                yy = map(y, miny, maxy, 50, screen_height-50);

                if(clock.getElapsedTime().asSeconds() > 1.5) {
                    if(colored) rect(xx, yy, 1, 1, color, sf::BlendAdd);
                    else rect(xx, yy, 1, 1, color);
                }
            }

            //std::cout << x << ' ' << y << '\n';

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