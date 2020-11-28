#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Audio.hpp>
#include <string>
#include "../flame/flame.hpp"
#include "../graphics/graphics.hpp"

#ifdef OUTPUT
#include <iostream>
#endif

#pragma once

namespace aural {
namespace app {
    class App {
        public:
            int screen_width, screen_height;
            sf::RenderTexture texture;
            unsigned long long seed;
            sf::Clock clock;
            std::string window_name, save_path;
            sf::RenderWindow window;
            double timeLimit;

            virtual void initWindow() {
                sf::ContextSettings settings;
                settings.antialiasingLevel = 8;
                window.create(sf::VideoMode(screen_width, screen_height), window_name, sf::Style::Close, settings);
                window.setVerticalSyncEnabled(true);
                window.setFramerateLimit(60);
            }

            virtual void init() {
                #ifdef WINDOW
                initWindow();
                #endif
                texture.create(screen_width, screen_height);
                texture.setSmooth(true);
                clock.restart();
                seed = std::chrono::system_clock::now().time_since_epoch().count();
            }

            virtual void setup() {

            }

            virtual void checkForEvents() {
                sf::Event event; 
                    while (window.pollEvent(event))
                        if (event.type == sf::Event::Closed) 
                            window.close();
            }

            virtual void drawTextureToWindow() {
                window.clear();
                window.draw(sf::Sprite(texture.getTexture()));
                window.display();
            }

            virtual void loop() {
                #ifdef WINDOW
                checkForEvents();
                drawTextureToWindow();
                #endif
            }

            virtual std::string saveFile() {
                sf::Image image = texture.getTexture().copyToImage();
                std::string path = save_path + std::to_string(seed) + ".jpg";
                image.saveToFile(path);
                return path;
            }

            virtual void close() {
                #ifdef WINDOW
                if(window.isOpen()) window.close();
                #endif
                #ifdef OUTPUT
                std::string info = saveFile() + "#method=method";
                info+="#seed="+std::to_string(seed);
                std::cout << info;
                #else
                saveFile();
                #endif
            }
            App() : screen_width(800), screen_height(600), window_name("aural application"), window(), save_path("./"), timeLimit(60) {}
            App(int w, int h) : screen_width(w), screen_height(h), window_name("aural application"), window(), save_path("./"), timeLimit(60) {}
        private:
    };


}
}



#define AURAL_APP(App) int main() { App.init(); renderer = &App.texture; aural::rnd::seedEngine(App.seed); aural::math::initNoise(App.seed); aural::flame::initFlame(); App.setup(); while(App.clock.getElapsedTime().asSeconds() < App.timeLimit) { App.loop(); } App.close(); }