#define WINDOW
#define OUTPUT

#include "aural.hpp"
#include <iostream>

using namespace au;
using namespace au::app;
using namespace au::math;
using namespace au::rnd;
using namespace au::flame;

#define phi (sqrt(5.0)+1.0)/2.0

class baseApp : public App {
    using App::App;
    public:
        double a[7], f[7], a1[7], f1[7], n[7], r[7], v, p[7], x = 0, y = 0, t = 0, T = 0;
        sf::Color colors[3];
        double minx = 10, miny = 10, maxx = -10, maxy = -10;
        double curminx = 0, curminy = 0, curmaxx = 0, curmaxy = 0;
        bool colored = 0;
        double hue;

        double ssin(double x, double p) {
            if(p == 0) return asin(sin(x));
            else return pow(sin(x), p);
        }

        double ccos(double x, double p) {
            if(p == 0) return acos(cos(x));
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

            timeLimit = 10000;
            save_path = "./";

            for(int i = 1; i <= 6; i++)
                a[i] = random(0.7, 1.2) * (random(0, 100) > 50 ? -1: 1),
                f[i] = random(0.7, 1.2) * (random(0, 100) > 50 ? -1: 1),
                n[i] = random(0.0, 100000.0),
                r[i] = random(0.04, 0.2),
                p[i] = random(1, 5);
            v = random(0.001, 0.5);
            //p = random(1, 5);
            //q = random(1, 5);

            /*for(int i = 1; i <= 6; i++)
                std::cout << "a" << i << " = " << a[i] << '\n';
            for(int i = 1; i <= 6; i++)
                std::cout << "f" << i << " = " << f[i] << '\n';

            std::cout << "v = " << v << "\np = " << p << "\nq = " << q << '\n';*/

            hue = random(0, 360);

            if(1 || random(0, 100) < 30) {
                colored = 1;
                sf::Vertex background[] = {
                    sf::Vertex(sf::Vector2f(0, 0), convert(Hsv(((int)hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06)))),
                    sf::Vertex(sf::Vector2f(screen_width, 0), convert(Hsv(((int)hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06)))),
                    sf::Vertex(sf::Vector2f(screen_width, screen_height), convert(Hsv(((int)hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06)))),
                    sf::Vertex(sf::Vector2f(0, screen_height), convert(Hsv(((int)hue+random(100, 260))%360, 0.2, 1.4*random(0.04, 0.06))))
                };

                texture.draw(background, 4, sf::Quads);

                colors[0] = convert(Hsv(hue, 0.4, 0.7));
                colors[1] = convert(Hsv(((int)hue+random(100, 260))%360, 0.4, 0.7));
            } else {
                texture.clear(sf::Color(200, 200, 200));
            
                colors[0] = convert(Hsv(hue, 0.4, 0.15));
                colors[1] = convert(Hsv(((int)hue+random(100, 260))%360, 0.4, 0.15));
            }
            
            //std::cout << seed << '\n';
        }
        void loop() {
            #ifdef WINDOW
            checkForEvents();
            #endif

            texture.clear(sf::Color(10, 10, 10));

            x = y = t = 0;
            minx = miny = 10000;
            maxx = maxy = -10000;

            hue+=1;
            if(hue > 360) hue = 0;

            colors[0] = convert(Hsv(hue, 0.4, 0.8));
            colors[1] = convert(Hsv(((int)hue+100+(int)(noise(hue)*100))%360, 0.4, 0.8));

            for(int i = 1; i <= 6; i++) {
                a1[i] = a[i];
                f1[i] = f[i];
                n[i] += r[i];
                a1[i] = constrain(a[i]*noise(n[i])*10, -1.2, 1.2);
                f1[i] = constrain(f[i]*noise(n[i])*10, -1.2, 1.2);
            }

            for(int i = 1; i <= 30000; i++) {
                double xx = x, yy = y;
                xx = a1[1]*ssin(f1[1]*x, p[1]) + a1[2]*ccos(f1[2]*y, p[2]) + a1[3]*ssin(f1[3]*t, p[3]);// + 10*pow(noise(x*10, y*10), q);
                yy = a1[4]*ccos(f1[4]*x, p[4]) + a1[5]*ssin(f1[5]*y, p[5]) + a1[6]*ssin(f1[6]*t, p[6]);// + 10*pow(noise(y*10, x*10), p);
                t += v;

                vec step(xx-x, yy-y);
                step = step*10;
                double tt = constrain(map(step.mag2(), 0, M_PI*M_PI*2, 0, 1), 0, 1);
                sf::Color color = interpolate(colors[0], colors[1], tt);
                color.a = 200;//(colored ? 10 : 15);
                //std::cout << step.mag2()*pow(10, 5) << ' ' << tt << "\n\n";
                x = xx;
                y = yy;

                minx = std::min(minx, x);
                miny = std::min(miny, y);
                maxx = std::max(maxx, x);
                maxy = std::max(maxy, y);

                /*xx = map(x, -M_PI, M_PI, 50, screen_width-50);
                yy = map(y, -M_PI, M_PI, 50, screen_height-50);*/
                xx = map(x, curminx, curmaxx, 50, screen_width-50);
                yy = map(y, curminy, curmaxy, 50, screen_height-50);

                if(1 || clock.getElapsedTime().asSeconds() > 1.5) {
                    if(colored) rect(xx, yy, 1, 1, color, sf::BlendAdd);
                    else rect(xx, yy, 2, 2, color);
                }
            }

            //std::cout << x << ' ' << y << '\n';

            double minp = 0.1, maxp = 1.0-minp;

            curminx = maxp*curminx+minp*minx;
            curminy = maxp*curminy+minp*miny;
            curmaxx = maxp*curmaxx+minp*maxx;
            curmaxy = maxp*curmaxy+minp*maxy;

            //std::cout << x << ' ' << y << '\n';

            // Declare and load a font
            sf::Font font;
            font.loadFromFile("font.ttf");
            // Create a text
            std::string s;
            sf::Text text("", font);
            for(int i = 1; i <= 6; i++)
                s += std::to_string(a1[i]) + " ";
            s += "\n";
            for(int i = 1; i <= 6; i++)
                s += std::to_string(f1[i]) + " ";
            s += "\n";
            for(int i = 1; i <= 6; i++)
                s += std::to_string(p[i]) + " ";
            text.setCharacterSize(25);
            text.setString(s);
            //text.setStyle(sf::Text::Bold);
            text.setFillColor(sf::Color::White);
            text.setPosition(10, screen_height/2);

            // Draw it
            texture.draw(text);

            texture.display();

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
} testApp(1600, 1600);

AURAL_APP(testApp)