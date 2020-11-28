#define NOISE
#include "sfmlpp.h"
#include "variations.h"
#include "FastNoise.h"
#include <iostream>

#define WINDOW

string NAME = "perlin2";
sf::RenderWindow window;
std::default_random_engine gen;
int WIDTH = 2048, HEIGHT = 2048;
bool vis = 0, is_colored = 1, circle = 0;
sf::Color bg1, bg2;
vector<sf::Color> pal;
FastNoise pn;
double xlimit = 3, ylimit = 3, timeLimit;
sf::RenderTexture renderTexture;
double weight, smoothing = 0, hue;
int octaves;

struct particle {
    double x, y;
    sf::Vector2f v;
    sf::Vertex vertex;
    particle(double x1, double y1, sf::Color c) {
        x = x1, y = y1, vertex.color = c;
    }
    void update(sf::Vector2f v1, double scale) {
        v.x = smoothing*v.x+(1-smoothing)*v1.x, v.y = smoothing*v.y+(1-smoothing)*v1.y;
        x += scale*v.x, y += scale*v.y;
    }
    
    void draw() {
        double xx = mathmap(x, -xlimit, xlimit, 10, WIDTH-10);
        double yy = mathmap(y, -ylimit, ylimit, 10, HEIGHT-10);
        vertex.position = sf::Vector2f(xx, yy);
        renderTexture.draw(&vertex, 1, sf::Points, sf::BlendAdd);
    }
};

vector<particle> points;
double step = 0.035, vector_scale = 0.001, noise_scale = 10, coord_scale = 50, contrast;
bool color_choice;
vector<vector<int> > tree1;
vector<node> nodes1;

void setup() {
    //if(rd(0, 100) < 20) circle = 1;
    if(rd(0, 100) < 50) {
        HEIGHT = 1152;
        xlimit = 5.33;
    }
    if(rd(0, 100) < 50) color_choice = 1;
    hue = rd(0, 360);
    contrast = 0;//rd(0, 1);
    octaves = 3;
    smoothing = constrain(rd(0, 0.99), 0, 1);
    for(double i = -xlimit; i <= xlimit; i += step)
        for(double j = -ylimit; j <= ylimit; j += step) {
            particle p(i+0.003*rdnormal(0, 1), j+0.003*rdnormal(0, 1), sf::Color(0, 0, 0, 30));
            points.push_back(p);
        }
    initialise_variations();

    createFieldTree2(rdint(4, 6), rdint(3, 4), tree1, nodes1, 2);
    vector<double> coord_magnitudes = {10, 50, 100, 500, 1000};
    vector<double> noise_magnitudes = {5, 10, 20, 100};
    coord_scale = coord_magnitudes[rdint(0, coord_magnitudes.size()-1)];
    noise_scale = noise_magnitudes[rdint(0, noise_magnitudes.size()-1)];
    //ylimit = 3;//rd(2.8, 3.7);
    timeLimit = rd(15, 30);
    weight = 1;
    /*string s;
    printFieldTree2(1, s, tree1, nodes1);
    cout << s << endl;*/
}

sf::Clock timerClock;

void draw() {
    int index = 0;
    for(auto &i : points) {
        /*if(is_colored) {
            int cn = (int)(100*pal.size()*mathmap(pn.GetNoise(index, 0), -1, 1, 0, 1))%pal.size();
            i.vertex.color = pal[cn];
        }*/
        
        sf::Vector2f v(i.x, i.y);
        v = mulF(resolveFieldTree2(1, v, weight, tree1, nodes1), sf::Vector2f(2, 2));
        double col = mathmap(sin(sqrtf(v.x*v.x+v.y*v.y)), -1, 1, 0, 1);
        i.vertex.color = interpolate(interpolate(pal[2], pal[3], mathmap(sin(atan2(v.y, v.x)), -1, 1, 0, 1)), interpolate(pal[0], pal[1], mathmap(sin(atan2(v.x, v.y)), -1, 1, 0, 1)), mathmap(i.x, -xlimit, xlimit, 0, 1));
        i.vertex.color.a = constrain(pow(sqrtf(v.x*v.x+v.y*v.y), 3)*4, 0, 80);
        //cout << v.x << ' ' << v.y << endl;
        i.update(v, vector_scale);
        i.draw();
        //i.edgecheck();
        index++;
    }
}

int main() {
    timerClock.restart();
    auto seed = seedgen();
    setup();
    #ifdef WINDOW
    init(window, "flowfield", 0, 60, WIDTH, HEIGHT);
    #endif
    renderTexture.create(WIDTH, HEIGHT);
    //init(renderTexture, "perlin2", false, 60, WIDTH, HEIGHT);
    //if(!vis) renderTexture.setVisible(0);
    pn.SetNoiseType(FastNoise::PerlinFractal);
    pn.SetFractalOctaves(octaves);
    pn.SetSeed(seed);

    if(is_colored) {
        pal = /*{sf::Color(230, 255, 220, 30), sf::Color(170, 180, 150, 30)};//*/randomPalette(4, 50, rd(0.5, 0.9), rd(0.2, 0.6), hue, 2);
        //for(auto i : pal) cout << (int)i.r << ' ' << (int)i.g << ' ' << (int)i.b << ' ' << (int)i.a << endl;
        int brightness = mathmap(contrast, 0, 1, 30, 220);
        //cout << brightness << endl;
        bg1 = /*sf::Color(30, 10, 20);//*/sf::Color(constrain(rdnormal(brightness, brightness/20), 0, 255), constrain(rdnormal(brightness, brightness/20), 0, 255), constrain(rdnormal(brightness, brightness/20), 0, 255));
        bg2 = /*sf::Color(10, 7, 20);//*/sf::Color(rdnormal(brightness/1.2, brightness/50), rdnormal(brightness/1.2, brightness/50), rdnormal(brightness/1.2, brightness/50));
    } else {
        bg1 = sf::Color(rdnormal(220, 1), rdnormal(220, 1), rdnormal(220, 1));
        bg2 = sf::Color(rdnormal(180, 1), rdnormal(180, 1), rdnormal(180, 1));
    }

    /*renderTexture.clear(bg1);
    vector<sf::Color> pal1 = {randomPalette(1, 250, 0.98, 0.2)[0], randomPalette(1, 250, 0.98, 0.8)[0]};
    for(int i = 0; i <= WIDTH; i++) {
        for(int j = 0; j <= HEIGHT; j++) {
            sf::Vertex point;
            point.color = interpolate(pal1[0], pal1[1], mathmap(pn.GetNoise(i/20.0, j/20.0), -1, 1, 0, 1));
            point.position = sf::Vector2f(i, j);
            //point.position = sf::Vector2f(i+0.003*rdnormal(0, 1), j+0.003*rdnormal(0, 1));
            renderTexture.draw(&point, 1, sf::Points);
        }
    }*/

    renderTexture.clear(bg1);
    double invert_hue = (hue + 180 > 360 ? hue - 180 : hue + 180);
    sf::Vertex rectangle[] =
    {
        sf::Vertex(sf::Vector2f(0, 0), randomPalette(1, rd(30, 80), 0.98, rd(0.2, 0.4), rdnormal(invert_hue, 10))[0]),
        sf::Vertex(sf::Vector2f(WIDTH, 0), randomPalette(1, rd(30, 80), 0.98, rd(0.2, 0.4), rdnormal(invert_hue, 10))[0]),
        sf::Vertex(sf::Vector2f(WIDTH, HEIGHT), randomPalette(1, rd(30, 80), 0.98, rd(0.2, 0.4), rdnormal(invert_hue, 10))[0]),
        sf::Vertex(sf::Vector2f(0, HEIGHT), randomPalette(1, rd(30, 80), 0.98, rd(0.2, 0.4), rdnormal(invert_hue, 10))[0])
    };

    renderTexture.draw(rectangle, 4, sf::Quads);

    auto image = renderTexture.getTexture().copyToImage();
    for(int i = 1; i < HEIGHT; i++) {
        for(int j = 1; j < WIDTH; j++) {
            auto pixelColor = image.getPixel(j, i);
            pixelColor = sf::Color(constrain(rdnormal(pixelColor.r, 2), pixelColor.r/2, min(255, pixelColor.r*2)), constrain(rdnormal(pixelColor.g, 2), pixelColor.g/2, min(255, pixelColor.g*2)), constrain(rdnormal(pixelColor.b, 2), pixelColor.b/2, min(255, pixelColor.b*2)));
            sf::Vertex newPixel;
            newPixel.position = sf::Vector2f(j, i);
            newPixel.color = pixelColor;
            renderTexture.draw(&newPixel, 1, sf::Points);
        }
    }
    
    //if(!vis) renderTexture.setVisible(0);
	while (timerClock.getElapsedTime().asSeconds() < timeLimit) {
		#ifdef WINDOW
		sf::Event event;
        while (window.pollEvent(event)) {
			if (event.type == sf::Event::Closed)
				window.close();
        }
        #endif

        sf::Vertex point;
        point.color = sf::Color(0, 0, 0, 0);
        point.position = sf::Vector2f(rd(0, WIDTH), rd(0, HEIGHT));
        renderTexture.draw(&point, 1, sf::Points);

        draw();

        /*for(int i = 0; i < pal.size(); i++) {
            sf::RectangleShape rec;
            rec.setSize(sf::Vector2f(200, 200));
            rec.setFillColor(pal[i]);
            rec.setPosition(WIDTH/2+i*250, HEIGHT/2);
            renderTexture.draw(rec);
        }*/

		renderTexture.display();

        #ifdef WINDOW
        window.clear();
        window.draw(sf::Sprite(renderTexture.getTexture()));
        window.display();
        #endif
	}

    #ifdef WINDOW
    window.close();
    #endif

    sf::Image screenshot = renderTexture.getTexture().copyToImage();
    string path = "./" + to_string(seed) + ".jpg";
    screenshot.saveToFile(path);
    string info = path + "#method=noise flowfield#seed="+to_string(seed);
    cout << info;
    //cout << endl << timerClock.getElapsedTime().asSeconds();
}
