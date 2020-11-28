#define WINDOW
#define OUTPUT

#include "aural.hpp"
#include <iostream>

using namespace au;
using namespace au::app;
using namespace au::math;
using namespace au::rnd;
using namespace au::flame;

const int color_number = 5;

class baseApp : public App {
    using App::App;
    public:
        struct particle {
            //randomizable
            double radius, coord_scale, centerx, centery;

            vec pos, vel;
            int age, id, maxage;
            double mood;
            sf::Color color;
            void reset(double t) {
                age = maxage = random(500, 4000);
                double angle = random(0.0, 2*M_PI);
                //mood = abs(map(angle, 0, 2*M_PI, -1, 1));
                //double rad = gaussian(radius, 0.01);// * (random(0, 100) < 80 ? map(gaussian(1, 0.2), 0, 2, 0, 1.2) : 1);
                double rad = 2*gaussian(0, 0.1) + (radius)*(1.0-pow(random(1.0),7.0));
                pos = vec(cos(angle)*rad+centerx, sin(angle)*rad+centery);
                vel = vec(0, 0);
                calcMood(t);
            }
            particle(int i, sf::Color c, double t, double r, double co, double cx, double cy) {
                id = i;
                color = c;
                radius = r;
                coord_scale = co;
                centerx = cx, centery = cy;
                //c.a = 10;
                reset(t);
            }
            void update(double t) {
                pos = pos + vel;
                calcMood(t);
                if(vec(pos.x-centerx, pos.y-centery).mag() < 100) age = maxage = 5;
                if(age-- < 0 /*|| pos.x < 0 || pos.x > 1800 || pos.y < 0 || pos.y > 1800*/) reset(t);
            }
            void calcMood(double t) {
                //mood = sin(waves(pos*10, 1).x);
                vec v = pos*coord_scale;
                mood = sin(noise(v.x, v.y, t)*2*M_PI);
            }
            double moodSimilarity(particle p) {
                return 1-abs(p.mood-mood);
            }
        };

        struct colony {
            double radius, proximity, prox1, prox2, coord_scale, time_step, noise_scale, centerx, centery;
            double power, multiplier;
            sf::Color colors[color_number+1];
            double t = 0;
            int N, id;
            std::vector<particle> particles;
            colony(double r, double p, double p1, double p2, double c, double n, double ts, double cx, double cy, int NN, int i, double pow, double mult) :
                radius(r), proximity(p), prox1(p1), prox2(p2), coord_scale(c), noise_scale(n), time_step(ts), centerx(cx), centery(cy), N(NN), id(i), power(pow), multiplier(mult) {}
            
            void addParticles() {
                for(int i = 1; i <= N; i++) {
                    particles.push_back(particle(i, colors[random(0, 5)], t, radius, coord_scale, centerx, centery));
                    particles.back().color.a = 200;
                }
            }

            void update(std::vector<colony> colonies) {
                for(auto &p : particles) {
                    int close = 0;
                    double lovex = 0, lovey = 0;
                    for(auto &c : colonies) {
                        if(c.id != id) continue;
                        for(auto &o : c.particles) {
                            if(p.id == o.id) continue;
                            vec v = p.pos-o.pos;
                            double dis = v.mag();
                            double angle = v.arctan();
                            double love = pow(1.0/std::max(1.0, dis), power)*multiplier;
                            vec f = p.pos-vec(centerx, centery);
                            if(dis < proximity) love *= map(f.mag(), 0, radius, prox1, prox2);
                            //std::cout << dis << '\n';
                            //love *= map(f.mag(), 0, radius, 0.5, 2);
                            if(dis < 50) {
                                close++;
                                //p.color.a = map(dis, 0, 100, 0, 6);
                            } //else p.color.a = 130;
                            love *= p.moodSimilarity(o);
                            love *= 20;
                            if(c.id != id) love *= 0.1;
                            //std::cout << love << '\n';
                            lovex += -cos(angle)*love;
                            lovey += -sin(angle)*love;
                        }
                    }
                    //std::cout << lovex << ' ' << lovey << '\n';
                    p.vel = vec(lovex, lovey);
                    
                    vec pp = p.pos*coord_scale;
                    double n = noise(pp.x, pp.y)*noise_scale;
                    vec nv(cos(n), sin(n));
                    //std::cout << nv.x << ' ' << nv.y << '\n';
                    vec f = p.pos-vec(centerx, centery);
                    double mod = map(noise(p.pos.x, p.pos.y, t), -1, 1, 0, 1);
                    //p.vel = p.vel * nv * 5;

                    close = std::min(close, 30);
                    p.color.a = map(close, 0, 30, 220, 30);
                    p.color.a *= map(p.age, 0, p.maxage, 1, 0.1);
                    p.color.a *= constrain(map(f.mag(), 0, radius*2, 1, 0), 0, 1);
                    if(f.mag() > radius) p.color.a *= 0.9;
                    p.update(t);
                }

            t += time_step;
            }

            void draw() {
                for(auto &p : particles) {
                    rect(p.pos.x, p.pos.y, 2, 2, p.color, sf::BlendAdd);
                }
            }
        };

        

        sf::Color colors[color_number+1];
        std::vector<particle> particles;
        double t = 0, angle = 0;
        int N = 700;
        double radius = 700, proximity = 2, prox1 = -1, prox2 = -1, coord_scale = 5, time_step = 0.01, noise_scale = 7, centerx = 900, centery = 900;
        std::vector<colony> colonies;

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

            timeLimit = 60;
            save_path = "./";

            int hue = random(0, 360);

            int colony_number = random(1, 5);

            double sum_angle = 0;

            for(int j = 1; j <= colony_number; j++) {
                colors[0] = convert(Hsv((hue+random(0, 20))%360, random(0.1, 0.65), random(0.6, 0.95)));
                for(int i = 1; i <= color_number; i++)
                    colors[i] = convert(Hsv((hue+random(150, 200))%360, random(0.1, 0.65), random(0.6, 0.95)));
                //colony(double r, double p, double p1, double p2, double c, double n, double ts, double cx, double cy, int NN)
                radius = random(300.0, 900.0) * map(colony_number, 1, 5, 2, 0.5);

                double angle = fmod(sum_angle + random(M_PI*0.8, M_PI*1.2), 2*M_PI);
                if(sum_angle) {
                    sum_angle += angle;
                    sum_angle /= 2;
                } else {
                    sum_angle += angle;
                }
                double center_radius = constrain(gaussian(900, 500), 0, 1000);
                centerx = cos(angle)*center_radius+screen_width/2;
                centery = sin(angle)*center_radius+screen_height/2;

                //centerx = random(0.0, (double)screen_width);
                //centery = random(0.0, (double)screen_height);

                proximity = random(1.0, 100.0);
                prox1 = random(0.1, 10.0) * (random(0, 100) < 50 ? -1 : 1); 
                prox2 = random(0.1, 10.0) * (random(0, 100) < 50 ? -1 : 1); 
                coord_scale = random(0.5, 10.0);
                noise_scale = random(0.5, 10.0);
                time_step = random(0.0001, 0.5);
                double power = random(1.0, 2.0);
                double mult = pow(map(power, 1, 2, 1, 10), 2)/2;
                colonies.push_back(colony(radius, proximity, prox1, prox2, coord_scale, noise_scale, time_step, centerx, centery, random(400, 600), j, power, mult));
                
                for(int k = 0; k <= color_number; k++) colonies.back().colors[k] = colors[k];
                colonies.back().addParticles();
            }

            sf::Vertex background[] = {
                sf::Vertex(sf::Vector2f(0, 0), convert(Hsv((hue+random(10, 50))%360, 0.6, random(0.03, 0.06)))),
                sf::Vertex(sf::Vector2f(screen_width, 0), convert(Hsv((hue+random(10, 50))%360, 0.6, random(0.03, 0.06)))),
                sf::Vertex(sf::Vector2f(screen_width, screen_height), convert(Hsv((hue+random(10, 50))%360, 0.6, random(0.03, 0.06)))),
                sf::Vertex(sf::Vector2f(0, screen_height), convert(Hsv((hue+random(10, 50))%360, 0.6, random(0.03, 0.06))))
            };

            texture.draw(background, 4, sf::Quads);

            //texture.clear(sf::Color(220, 220, 220));
        }
        void loop() {
            #ifdef WINDOW
            checkForEvents();
            #endif

            //std::cout << colonies.back().centerx << ' ' << colonies.back().centery << '\n';

            for(auto &c : colonies) {
                c.update(colonies);
                c.draw();
            }

            #ifdef WINDOW
            texture.display();
            drawTextureToWindow();
            #endif
        }

        void close() {
            #ifdef WINDOW
            if(window.isOpen()) window.close();
            #endif
            #ifdef OUTPUT
            std::string info = saveFile() + "#method=organic growth";
            info+="#seed="+std::to_string(seed);
            std::cout << info;
            #else
            saveFile();
            #endif
        }
} testApp(1800, 1800);

AURAL_APP(testApp)