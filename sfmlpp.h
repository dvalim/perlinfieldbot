#ifndef SFMLPP_H
#define SFMLPP_H

#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Audio.hpp>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <set>
#include "PerlinNoise.h"
#include <chrono>

#define xMouse(window) sf::Mouse::getPosition(window).x
#define yMouse(window) sf::Mouse::getPosition(window).y
#define isMousePressed(button) sf::Mouse::isButtonPressed(sf::Mouse::Button::button)
#define isKeyPressed(key) sf::Keyboard::isKeyPressed(sf::Keyboard::key)

using namespace std;

extern int WIDTH, HEIGHT;
extern sf::RenderWindow window2;
extern std::default_random_engine gen;
extern sf::Font font;
typedef std::chrono::system_clock myclock;

#define PI 3.14159265
#define M_PI 3.14159265

unsigned long long seedgen() { //zasiewa generator liczb losowych, tylko potrzebne do init()
	myclock::time_point beginning = myclock::now();
	unsigned long long seed = beginning.time_since_epoch().count();
	gen.seed(seed);
	return seed;
}
unsigned long long seedgen(unsigned long long seed) { //zasiewa generator liczb losowych, tylko potrzebne do init()
	gen.seed(seed);
	return seed;
}

void getScreenSize() {
	sf::VideoMode vid = sf::VideoMode::getDesktopMode();
	WIDTH = vid.width, HEIGHT = vid.height;
}

void init(sf::RenderWindow &window, string NAME, bool visible, int framerate, bool fullscreen) { //pierwszy bool ustawia widocznosc konsoli, int ilosc klatek, drugi bool czy bedzie fullscreen
	//if(!visible) ShowWindow(GetConsoleWindow(), SW_HIDE);
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;

	sf::VideoMode vid = sf::VideoMode::getFullscreenModes()[0];

	//sf::VideoMode vid = sf::VideoMode::getDesktopMode();
	WIDTH = vid.width, HEIGHT = vid.height;
	if(!fullscreen) {
		window.create(sf::VideoMode(vid.width, vid.height), NAME, sf::Style::Default, settings);
	}
	else {
		window.create(vid, NAME, sf::Style::None, settings);
		//WIDTH = window.getSize().x;
		//HEIGHT = window.getSize().y;
	}
	window.setVerticalSyncEnabled(true);
	window.setFramerateLimit(framerate);

	seedgen();
}

void init(sf::RenderWindow &window, string NAME, bool visible, int framerate, int w, int h) { //zamiast boola fullscreen mozna podac wymiary
	//if(!visible) ShowWindow(GetConsoleWindow(), SW_HIDE);
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	window.create(sf::VideoMode(w, h), NAME, sf::Style::Close, settings);
	window.setVerticalSyncEnabled(true);
	window.setFramerateLimit(framerate);

	//seedgen();

	WIDTH = w, HEIGHT = h;
}

void init(sf::RenderWindow &window, string NAME, bool visible, int framerate, int w, int h, bool none) { 
	//if(!visible) ShowWindow(GetConsoleWindow(), SW_HIDE);
	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	if(none) window.create(sf::VideoMode(w, h), NAME, sf::Style::None, settings);
	else window.create(sf::VideoMode(w, h), NAME, sf::Style::Close, settings);;
	window.setVerticalSyncEnabled(true);
	window.setFramerateLimit(framerate);

	seedgen();
}

sf::Color interpolate(sf::Color color1, sf::Color color2, double percent) {
	double resultRed = color1.r + percent * (color2.r - color1.r);
	double resultGreen = color1.g + percent * (color2.g - color1.g);
	double resultBlue = color1.b + percent * (color2.b - color1.b);
	double resultAlpha = color1.a + percent * (color2.a - color1.a);
	return sf::Color(resultRed, resultGreen, resultBlue, resultAlpha);
}

struct rgb {
	double r, g, b;
};

struct hsv {
	double h, s, v;
	hsv(double h1, double s1, double v1) {
		h = h1, s = s1, v = v1;
	}
};

rgb hsv2rgb(hsv in)
{
    double      hh, p, q, t, ff;
    long        i;
    rgb         out;

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
    return out;     
}

sf::Color toColor(rgb col, int a) {
	return sf::Color(col.r, col.g, col.b, a);
}

double mathmap(double a, double b, double c, double d, double e) { //tlumaczy wartosc z jednej wielkosci na druga
	return (a - b) / (c - b) * (e - d) + d;						   //np 50 jest pomiedzy 0 a 100, a jesli chcesz te wartosc miedzy 0 a 1000
}																   //to mathmap(50, 0, 100, 0, 1000) zwroci 500

double rd(double a, double b) { //zwraca losowa liczbe pomiedzy a i b
	std::uniform_real_distribution<double> rd1(a, b);
	return rd1(gen);
}

double rdnormal(double a, double b) {
	std::normal_distribution<double> rd1(a, b);
	return rd1(gen);
}

int rdint(int a, int b) { //zwraca losowa liczbe pomiedzy a i b
	std::uniform_int_distribution<int> rd2(a, b);
	return rd2(gen);
}

double constrain(double a, double b, double c) { //zamyka liczbe a pomiedzy b i c
	if(a < b) return b;
	else if(a > c) return c;
	else return a;
}

vector<sf::Color> randomPalette(int len, int alpha, double sat, double vib, double h) {
	vector<sf::Color> res;
	double golden_ratio_conjugate = 0.618033988749895;
	int core = 1;//rdint(1, len/2+len%2);
	vector<double> hues;
	for(int i = 1; i <= core; i++) {
		h = fmod((h/360.0 + golden_ratio_conjugate), 1) * 360.0;
		hues.push_back(h);
		hsv hs(h, sat, vib);
		auto rg = hsv2rgb(hs);
		sf::Color c(255*rg.r, 255*rg.g, 255*rg.b, alpha);
		res.push_back(c);
	}
	for(int i = 1; i <= len-core; i++) {
		auto hue = hues[rdint(0, core-1)];
		h = rdnormal(hue, 60);
		hsv hs(h, constrain(rdnormal(sat*0.95, 0.3), 0, 1), constrain(rdnormal(vib*0.95, 0.2), 0, 1));
		auto rg = hsv2rgb(hs);
		sf::Color c(255*rg.r, 255*rg.g, 255*rg.b, alpha);
		res.push_back(c);
	}
	for(int i = 0; i < len; i++)
		swap(res[i], res[rdint(0, len)]);
	return res;
}

vector<sf::Color> randomPalette(int len, int alpha, double sat, double vib, double h, int core) {
	vector<sf::Color> res;
	double golden_ratio_conjugate = 0.618033988749895;
	vector<double> hues;
	for(int i = 1; i <= core; i++) {
		h = fmod((h/360.0 + golden_ratio_conjugate), 1) * 360.0;
		hues.push_back(h);
		hsv hs(h, sat, vib);
		auto rg = hsv2rgb(hs);
		//cout << alpha << ' ';
		sf::Color c(255*rg.r, 255*rg.g, 255*rg.b, alpha);
		//cout << (int)c.a << endl;
		res.push_back(c);
	}
	for(int i = 1; i <= len-core; i++) {
		auto hue = hues[rdint(0, core-1)];
		h = rdnormal(hue, 60);
		hsv hs(h, constrain(rdnormal(sat*0.95, 0.3), 0, 1), constrain(rdnormal(vib*0.95, 0.2), 0, 1));
		auto rg = hsv2rgb(hs);
		sf::Color c(255*rg.r, 255*rg.g, 255*rg.b, alpha);
		res.push_back(c);
	}
	/*for(int i = 0; i < len; i++)
		swap(res[i], res[rdint(0, len)]);*/
	//cout << (int)res[0].a << ' ' << (int)res[1].a << endl;
	return res;
}

sf::Vector2f VectorFromAngle(double angle) { //zwraca wektor z kata w stopniach
	sf::Vector2f v;
	double radrotation = angle / 360 * acos(-1.0) * 2;
	v.x = sin(radrotation);
	v.y = -cos(radrotation);
	return v;
}

sf::Vector2f VectorFromRadianAngle(double angle, double mag) { //zwraca wektor z kata w stopniach
	sf::Vector2f v;
	v.x = sin(angle)*mag;
	v.y = -cos(angle)*mag;
	return v;
}

sf::Vector2f VectorFromAngle(double angle, double mag) { //mozna podac magnitude o ktore pomnozy wektor
	sf::Vector2f v;
	double radrotation = angle / 360 * acos(-1.0) * 2;
	v.x = sin(radrotation) * mag;
	v.y = -cos(radrotation) * mag;
	return v;
}

double dist(double x1, double y1, double x2, double y2) { //podaje dystans pomiedzy dwoma punktami
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

double dist(sf::Transformable s1, sf::Transformable s2) { //podaje dystans pomiedzy dwoma ksztaltami
	return sqrt((s1.getPosition().x - s2.getPosition().x)*(s1.getPosition().x -
	 s2.getPosition().x) + (s1.getPosition().y - s2.getPosition().y)*(s1.getPosition().y -
	  s2.getPosition().y));
}

double dist(sf::Vertex * s1, sf::Transformable s2) { //podaje dystans miedzy tablica sf::Vertex a ksztaltem
	return sqrt((s1[0].position.x - s2.getPosition().x)*(s1[0].position.x -
	 s2.getPosition().x) + (s1[0].position.y - s2.getPosition().y)*(s1[0].position.y -
	  s2.getPosition().y));
}

struct rainbowParticle { //tecza. nie ma funkcji zwracajacej kolor, po prostu RainbowParticle rp, a potem rp.color
	sf::Color color = sf::Color::Black;
	double cred, cgreen, cblue, balance, alpha = 255;
	bool yellow, green, cyan, blue, pink;
	rainbowParticle() { //default 
		balance = 6;
		yellow = green = cyan = blue = pink = 0;
		cred = cgreen = cblue = 0;
	}
	rainbowParticle(double b) { //z podana liczba, o ktora zmienia sie kolor
		balance = b;
		yellow = green = cyan = blue = pink = 0;
		cred = cgreen = cblue = 0;
	}
	rainbowParticle(double b, double a) { //z podana liczba i przezroczystoscia
		balance = b;
		yellow = green = cyan = blue = pink = 0;
		cred = cgreen = cblue = 0;
		alpha = a;
		color = sf::Color(cred, cgreen, cblue, alpha);
	}
	void update() {
		if (!yellow) {
			cgreen += balance;
			if (cgreen >= 255) cgreen = 255, yellow = 1;
		}
		else if (!green) {
			cred -= balance;
			if (cred <= 0) cred = 0, green = 1;
		}
		else if (!cyan) {
			cblue += balance;
			if (cblue >= 255) cblue = 255, cyan = 1;
		}
		else if (!blue) {
			cgreen -= balance;
			if (cgreen <= 0) cgreen = 0, blue = 1;
		}
		else if (!pink) {
			cred += balance;
			if (cred >= 255) cred = 255, pink = 1;
		}
		else {
			cblue -= balance;
			if (cblue <= 0) {
				cblue = 0;
				yellow = green = cyan = blue = pink = 0;
			}
		}
		color = sf::Color(cred, cgreen, cblue, alpha);
	}
};

sf::Color offsetColor(sf::Color color, int offset) { //zmienia kolor o offset
	color.r = constrain(color.r + offset, 0, 255);
	color.g = constrain(color.g + offset, 0, 255);
	color.b = constrain(color.b + offset, 0, 255);
	return color;
}

sf::Color offsetColor(sf::Color color, int offset, int alpha) { //zmienia kolor o offset i przezroczystosc o alpha
	color.r = constrain(color.r + offset, 0, 255);
	color.g = constrain(color.g + offset, 0, 255);
	color.b = constrain(color.b + offset, 0, 255);
	color.a = constrain(color.a + alpha, 0, 255);
	return color;
}

sf::Color offsetColor(sf::Color color, int offset1, int offset2, int offset3) { //zmienia kazda czesc koloru o inny offset
	color.r = constrain(color.r + offset1, 0, 255);
	color.g = constrain(color.g + offset2, 0, 255);
	color.b = constrain(color.b + offset3, 0, 255);
	return color;
}

sf::Color offsetColor(sf::Color color, int offset1, int offset2, int offset3, int alpha) { //j/w z przezroczystoscia
	color.r = constrain(color.r + offset1, 0, 255);
	color.g = constrain(color.g + offset2, 0, 255);
	color.b = constrain(color.b + offset3, 0, 255);
	color.a = constrain(color.a + alpha, 0, 255);
	return color;
}

struct Slider { //suwak
	sf::RectangleShape bg, rect;
	sf::CircleShape grip, circ1, circ2;
	double value, minval, maxval;
	bool setting = 0;
	string name;

	Slider(string s, double minval1, double maxval1, double x, double y, double w, sf::Color bgcolor, sf::Color rectcolor, sf::Color gripcolor) {
		name = s;					//Slider slider(nazwa, najmniejsza wartosc, najwieksza wartosc, x, y, dlugosc, kolor tla, kolor suwaka, kolor koleczka);
		double h = w * 0.2;			//trzeba zaladowac czcionke o nazwie font
		bg.setPosition(x, y);		//tutaj akurat jest ustawione ze pojawiaja sie tylko w oknie window2, jak chcesz to zmien
		bg.setSize(sf::Vector2f(w, h));
		bg.setFillColor(bgcolor);

		rect.setSize(sf::Vector2f(w*0.75, h*0.05));
		rect.setOrigin(rect.getSize().x / 2, rect.getSize().y / 2);
		rect.setFillColor(rectcolor);
		rect.setPosition(x + w / 2, y + h / 2);

		circ1.setRadius(rect.getSize().y / 2);
		circ1.setOrigin(circ1.getRadius(), circ1.getRadius());
		circ1.setFillColor(rectcolor);
		circ1.setPosition(rect.getPosition().x - rect.getSize().x/2, rect.getPosition().y);

		circ2 = circ1;
		circ2.setPosition(rect.getPosition().x + rect.getSize().x / 2, rect.getPosition().y);

		grip.setRadius(rect.getSize().y*1.5);
		grip.setOrigin(grip.getRadius(), grip.getRadius());
		grip.setFillColor(gripcolor);
		grip.setPosition(rect.getPosition());

		minval = minval1;
		maxval = maxval1;
		value = (minval + maxval) / 2;
	}

	void setValue(double val) {
		val = min(maxval, max(minval, val));
		value = val;
		grip.setPosition(sf::Vector2f(mathmap(value, minval, maxval, rect.getPosition().x - rect.getSize().x / 2, rect.getPosition().x + rect.getSize().x / 2), grip.getPosition().y));
	}

	void show() {
		window2.draw(bg);
		window2.draw(rect);
		window2.draw(circ1);
		window2.draw(circ2);
		window2.draw(grip);

		sf::Text text;
		string s = to_string(value);
		text.setFont(font);
		text.setFillColor(rect.getFillColor());
		text.setString(s);
		text.setCharacterSize(grip.getRadius() * 3);
		double y = mathmap(2.5, 1, 10, bg.getPosition().y, rect.getPosition().y);
		text.setPosition(rect.getPosition().x - text.getCharacterSize() * s.size() / 4, y);
		window2.draw(text);

		text.setString(name);
		y = mathmap(2.5, 1, 10, rect.getPosition().y, bg.getPosition().y + bg.getSize().y);
		text.setPosition(rect.getPosition().x - text.getCharacterSize() * name.size() / 4.5, y);
		window2.draw(text);
	}

	void update() { //update i show trzeba caly czas wlaczac
		if (isMousePressed(Left) &&
			dist(xMouse(window2), yMouse(window2), grip.getPosition().x, grip.getPosition().y) <= grip.getRadius()) setting = 1;
		if (setting && !isMousePressed(Left)) setting = 0;
		if (setting) {
			double x = max((double)rect.getPosition().x - (double)rect.getSize().x / 2, min((double)rect.getPosition().x + (double)rect.getSize().x / 2, (double)xMouse(window2)));
			double y = grip.getPosition().y;
			grip.setPosition(sf::Vector2f(x, y));

			value = mathmap(grip.getPosition().x, rect.getPosition().x - rect.getSize().x / 2, rect.getPosition().x + rect.getSize().x / 2, minval, maxval);
		}
	}

	double getValue() {
		return value;
	}
};

struct Button { //przycisk
	sf::RectangleShape rect;
	sf::Text text;
	bool value = 0, activated = 0, pressed = 0;
	double alpha1 = 255, alpha2 = 0;
	Button(double x, double y, double size, string s, sf::Color textColor, sf::Color bgColor) {
		text.setFont(font);
		text.setFillColor(textColor);
		text.setString(s);
		text.setCharacterSize(size);
		text.setOrigin(size * s.size() / 3.75, size / 1.5);

		rect.setOutlineColor(textColor);
		rect.setOutlineThickness(size / 20);
		rect.setFillColor(bgColor);
		rect.setPosition(x, y);
		rect.setSize(sf::Vector2f(s.size() * size * 0.7, size * 1.5));

		text.setPosition(x + rect.getSize().x / 2, y + rect.getSize().y / 2);
	}

	void update() {
		value = 0;

		if (xMouse(window2) > rect.getPosition().x && xMouse(window2) < rect.getPosition().x + rect.getSize().x &&
			yMouse(window2) > rect.getPosition().y && yMouse(window2) < rect.getPosition().y + rect.getSize().y) {
			activated = 1;
		}
		else activated = 0;

		if (activated && isMousePressed(Left)) pressed = 1;

		double offset = 50;

		if (activated && alpha1) {
			alpha1-=offset;
			alpha1 = constrain(alpha1, 0, 255);
			alpha2+=offset;
			alpha2 = constrain(alpha2, 0, 255);
		}
		else if (!activated && alpha2) {
			alpha2-=offset;
			alpha2 = constrain(alpha2, 0, 255);
			alpha1+=offset;
			alpha1 = constrain(alpha1, 0, 255);
		}

		if (pressed && !isMousePressed(Left)) {
			pressed = 0;
			value = 1;
		}
	}

	void show() {
		auto color1 = rect.getFillColor();
		color1.a = alpha1;
		auto rect2 = rect;
		rect2.setFillColor(color1);
		window2.draw(rect2);

		color1 = text.getFillColor();
		color1.a = alpha2;
		rect2.setFillColor(color1);
		window2.draw(rect2);

		color1 = text.getFillColor();
		color1.a = alpha1;
		auto text1 = text;
		text1.setFillColor(color1);
		window2.draw(text1);

		color1 = rect.getFillColor();
		color1.a = alpha2;
		text1.setFillColor(color1);
		window2.draw(text1);
	}

	bool getValue() {
		return value;
	}
};

void checkEvents(sf::RenderWindow &window) {
	sf::Event event;
		while (window.pollEvent(event)) {
			if (event.type == sf::Event::Closed)
				window.close();
		}
}

sf::RectangleShape line(double x1, double y1, double x2, double y2, double thickness, sf::Color color) {
	double dis = dist(x1, y1, x2, y2), angle = atan2(y1 - y2, x1 - x2) * 180 / PI;
	sf::RectangleShape ret(sf::Vector2f(dis, thickness));
	ret.setOrigin(0, thickness/2);
	ret.setPosition(x2, y2);
	ret.rotate(angle);
	ret.setFillColor(color);
	return ret;
}

#endif