#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <iomanip>

struct Position {
    double x = 0;
    double y = 0;
    double vx = 0;
    double vy = 0;
};

struct Planet {
    double x = 0;
    double y = 0;
    double K = 0;
};

const double GM = 398600;
const double R_Earth = 6371;


Position get_next(const Position cur, const std::vector<Planet>& planets) {
    double ax = 0;
    double ay = 0;
    for (size_t i = 0; i != planets.size(); ++i) {
        double K = planets[i].K;
        double dirx = (planets[i].x - cur.x);
        double diry = (planets[i].y - cur.y);
        double dist = dirx * dirx + diry * diry;
        ax += (K * GM / dist) * (dirx / sqrt(dist));
        ay += (K * GM / dist) * (diry / sqrt(dist));
    }
    Position new_pos;

    new_pos.x = cur.x + cur.vx / 23 + ax / 1058;
    new_pos.y = cur.y + cur.vy / 23 + ay / 1058;
    new_pos.vx = cur.vx + ax / 23;
    new_pos.vy = cur.vy + ay / 23;
    return new_pos;
}


std::vector<std::pair<double, double>> get_way(size_t seconds, Position pos, const std::vector<Planet>& planets) {
    seconds *= 23;
    std::vector<std::pair<double, double>> res;
    for (size_t i = 0; i != seconds; ++i) {
        if (i % 23 == 0) {
           res.emplace_back(pos.x, pos.y);
        }
        pos = get_next(pos, planets);
    }
    res.emplace_back(pos.x, pos.y);
    return res;
}

void print_planets(std::ofstream &out, const std::vector<Planet>& planets) {
    out << planets.size() << "\n";
    for (size_t i = 0; i != planets.size(); ++i) {
        out << planets[i].K << " " <<
        planets[i].x << " " << planets[i].y << "\n";
    }
}

void print_sol(std::ofstream &out, const std::vector<std::pair<double, double>>& sol) {
    for (size_t j = 0; j != sol.size(); ++j) {
        out << std::fixed << std::setprecision(6) << sol[j].first << " "
        << std::fixed << std::setprecision(6) << sol[j].second << "\n";
    }
}

void task_a() {
    std::ofstream out;
    out.open("task_a.txt");
    size_t seconds = 10000;
    Planet Earth = {0, 0, 1};
    std::vector<Planet> planets;
    planets.push_back(Earth);
    print_planets(out, planets);
    std::vector<double> v = {5, 9, 12};
    out << v.size() << " " << seconds << "\n";
    for (size_t i = 0; i != v.size(); ++i) {
        Position init_pos = {0, R_Earth, v[i], 0};
        auto sol = get_way(seconds, init_pos, planets);
        print_sol(out, sol);
    }
}


// считаем радиус через секунду и смотрим на близость к исходному.
double task_b() {
    Planet Earth = {0, 0, 1};
    std::vector<Planet> planets;
    planets.push_back(Earth);
    double min_difference = 10000; // просто большое число для начала
    double diff = 10000;
    double v = 5; // знаем из a, что это меньше первой космической
    while (diff <= min_difference) {
        min_difference = diff;
        v += 0.0001;
        Position init_pos = {0, R_Earth, v, 0};
        auto res = get_way(1, init_pos, planets);
        double radius = sqrt(res.back().first * res.back().first +
                 res.back().second * res.back().second);
        diff = std::abs(R_Earth - radius);
    }
    v -= 0.0001;
    // и нарисуем
    std::ofstream out;
    out.open("task_b.txt");
    print_planets(out, planets);
    double seconds = 10000;
    out << "1 " << seconds << "\n";
    Position init_pos = {0, R_Earth, v, 0};
    auto sol = get_way(seconds, init_pos, planets);
    print_sol(out, sol);
    return v;
}

void task_c() {
    Planet Earth_1 = {-20000, 0, 1};
    Planet Earth_2 = {20000, 0, 1};
    std::vector<Planet> planets;
    planets.push_back(Earth_2);
    planets.push_back(Earth_1);
    size_t seconds = 20000;
    std::ofstream out;
    out.open("task_c.txt");
    print_planets(out, planets);
    out << "9 " << seconds << "\n";
    for (double v = 1; v <= 5; v += 0.5) {
        double vx = v / sqrt(2);
        double vy = v / sqrt(2);
        Position init_pos = {0, 0, vx, vy};
        auto sol = get_way(seconds, init_pos, planets);
        print_sol(out, sol);
    }
}

int main() {
    task_c();
}
