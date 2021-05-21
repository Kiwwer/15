#include "math.h"
#include "algorithm"
#include "fstream"
#include "set"
#include "iostream"
#include "iomanip"
#include "vector"
#include "string"
#include "map"
#include "queue"
#include "stack"
#include "unordered_map"
#include "chrono"
#include "thread"

#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

using namespace std;
using namespace cv;

#define ll long long
#define ld long double
#define elif else if
#define br break
#define flag bool
#define con continue
#define INF 1000000000
#define MOD 1000000007


struct Vect {
    ll x, y;
    Vect(ll _x, ll _y) : x(_x), y(_y) {}
    Vect() {}
    bool operator < (const Vect& m) const {
        return x < m.x || (x == m.x && y < m.y);
    }
    bool operator > (const Vect& m) const {
        return x > m.x || (x == m.x && y > m.y);
    }
    bool operator == (const Vect& m) const {
        return (x == m.x && y == m.y);
    }
    bool operator != (const Vect& m) const {
        return (x != m.x || y != m.y);
    }
};

struct Triangle {
    Vect a, b, c;
    Triangle(Vect _a, Vect _b, Vect _c) : a(_a), b(_b), c(_c) {}
    Triangle() {}
    bool in(Vect p) {
        ll pl1, pl2, pl3;
        pl1 = (a.x - p.x) * (b.y - a.y) - (b.x - a.x) * (a.y - p.y);
        pl2 = (b.x - p.x) * (c.y - b.y) - (c.x - b.x) * (b.y - p.y);
        pl3 = (c.x - p.x) * (a.y - c.y) - (a.x - c.x) * (c.y - p.y);
        if ((pl1 >= 0 && pl2 >= 0 && pl3 >= 0) || (pl1 <= 0 && pl2 <= 0 && pl3 <= 0))
        {
            return true;
        }
        return false;
    }
};

map <pair <Vect, Vect>, pair <ll, ll>> faces;
vector <Triangle> delaune;
vector <Vect> pts;

pair <ll, ll> findTriangle(Vect p) {
    pair <ll, ll> res = {-1, -1};
    for (ll i = 0; i < delaune.size(); ++i) {
        if (delaune[i].in(p)) {
            if (res.first == -1)
                res.first = i;
            elif(res.second == -1)
                res.second = i;
            else
                return { -1, -1 };
        }
    }
    return res;
}

pair <Vect, Vect> getedge(pair <Vect, Vect> semiedge) {
    if (semiedge.first > semiedge.second)
        return semiedge;
    return { semiedge.second, semiedge.first };
}


// Некрасивая функция проверки на плохое ребро
bool isinCircle(Vect p, Triangle tr) {
    ld xp = p.x;
    ld yp = p.y;
    ld x1 = tr.a.x;
    ld x2 = tr.b.x;
    ld x3 = tr.c.x;
    ld y1 = tr.a.y;
    ld y2 = tr.b.y;
    ld y3 = tr.c.y;
    ld m1, m2, mx1, mx2, my1, my2, dx, dy, rsqr, drsqr, xc, yc, r;

    if (abs(y1 - y2) < 1e-6 && abs(y2 - y3) < 1e-6)
        return 0;
    if (abs(y2 - y1) < 1e-6) {
        m2 = -(x3 - x2) / (y3 - y2);
        mx2 = (x2 + x3) / 2.0;
        my2 = (y2 + y3) / 2.0;
        xc = (x2 + x1) / 2.0;
        yc = m2 * (xc - mx2) + my2;
    }
    else if (abs(y3 - y2) < 1e-6) {
        m1 = -(x2 - x1) / (y2 - y1);
        mx1 = (x1 + x2) / 2.0;
        my1 = (y1 + y2) / 2.0;
        xc = (x3 + x2) / 2.0;
        yc = m1 * (xc - mx1) + my1;
    }
    else {
        m1 = -(x2 - x1) / (y2 - y1);
        m2 = -(x3 - x2) / (y3 - y2);
        mx1 = (x1 + x2) / 2.0;
        mx2 = (x2 + x3) / 2.0;
        my1 = (y1 + y2) / 2.0;
        my2 = (y2 + y3) / 2.0;
        xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
        yc = m1 * (xc - mx1) + my1;
    }
    dx = x2 - xc;
    dy = y2 - yc;
    rsqr = dx * dx + dy * dy;
    r = sqrt(rsqr);
    dx = xp - xc;
    dy = yp - yc;
    drsqr = dx * dx + dy * dy;
    return (drsqr <= rsqr) ? 1 : 0;
}

// Чуть менее некрасивая функция флипа
void flip(pair <Vect, Vect> edge) {
    pair <ll, ll> fc = faces.at(edge);
    if (fc.first == -1 || fc.second == -1)
        return;
    bool needed = 0;
    Vect lp1, lp2;
    if (edge.first == delaune[fc.first].a && edge.second == delaune[fc.first].b ||
        edge.first == delaune[fc.first].b && edge.second == delaune[fc.first].a) {
        lp1 = delaune[fc.first].c;
        needed |= isinCircle(lp1, delaune[fc.second]);
    } elif (edge.first == delaune[fc.first].b && edge.second == delaune[fc.first].c ||
        edge.first == delaune[fc.first].c && edge.second == delaune[fc.first].b) {
        lp1 = delaune[fc.first].a;
        needed |= isinCircle(lp1, delaune[fc.second]);
    } else {
        lp1 = delaune[fc.first].b;
        needed |= isinCircle(lp1, delaune[fc.second]);
    }

    if (edge.first == delaune[fc.second].a && edge.second == delaune[fc.second].b ||
        edge.first == delaune[fc.second].b && edge.second == delaune[fc.second].a) {
        lp2 = delaune[fc.second].c;
        needed |= isinCircle(lp2, delaune[fc.first]);
    } elif(edge.first == delaune[fc.second].b && edge.second == delaune[fc.second].c ||
        edge.first == delaune[fc.second].c && edge.second == delaune[fc.second].b) {
        lp2 = delaune[fc.second].a;
        needed |= isinCircle(lp2, delaune[fc.first]);
    } else {
        lp2 = delaune[fc.second].b;
        needed |= isinCircle(lp2, delaune[fc.first]);
    }

    if (needed) {
        faces.erase(edge);
        delaune[fc.first] = { lp1, lp2, edge.first };
        delaune[fc.second] = { lp1, lp2, edge.second };
        faces[getedge({ lp1, lp2 })] = { fc.first, fc.second };
        pair <ll, ll> cura = faces[getedge({ edge.first, lp2 })];
        pair <ll, ll> curb = faces[getedge({ edge.second, lp1 })];
        if (cura.first == fc.second) {
            cura.first = fc.first;
        } else {
            cura.second = fc.first;
        }
        if (curb.first == fc.first) {
            curb.first = fc.second;
        } else {
            curb.second = fc.second;
        }
        faces[getedge({ edge.first, lp2 })] = cura;
        faces[getedge({ edge.second, lp1 })] = curb;
    }
}

// Фейковый гигантский треугольник
vector <Vect> bannedpts({ {-1000000000, -1000000000}, {1000000000, -1000000000}, {0, 1000000000} });

// O(n) (гарантированно) локализация точки
// O(log n) (гарантированно) обновление
// -
// Итого вставка: O(n)
void addpoint(Vect p) {
    pair <ll, ll> loc = findTriangle(p);
    if (loc.first == -1 && loc.second == -1)
        return;
    vector <pair <Vect, Vect>> flips;
    if (loc.second == -1) {
        pair <Vect, Vect>
            a = getedge({ delaune[loc.first].a, delaune[loc.first].b }),
            b = getedge({ delaune[loc.first].b, delaune[loc.first].c }),
            c = getedge({ delaune[loc.first].c, delaune[loc.first].a });
        // Обновление
        {
            delaune.push_back({ delaune[loc.first].a, delaune[loc.first].b, p });
            delaune.push_back({ delaune[loc.first].b, delaune[loc.first].c, p });
            faces[getedge({ delaune[loc.first].a, p })] = { delaune.size() - 2, loc.first };
            faces[getedge({ delaune[loc.first].b, p })] = { delaune.size() - 2, delaune.size() - 1 };
            faces[getedge({ delaune[loc.first].c, p })] = { delaune.size() - 1, loc.first };
            delaune[loc.first] = { delaune[loc.first].c, delaune[loc.first].a, p };
        }
        {
            pair <ll, ll> af = faces[a];
            if (af.second == loc.first)
                swap(af.first, af.second);
            af.first = delaune.size() - 2;
            faces[a] = af;
        }
        {
            pair <ll, ll> bf = faces[b];
            if (bf.second == loc.first)
                swap(bf.first, bf.second);
            bf.first = delaune.size() - 1;
            faces[b] = bf;
        }
        {
            pair <ll, ll> cf = faces[c];
            if (cf.second == loc.first)
                swap(cf.first, cf.second);
            faces[c] = cf;
        }
        flips.push_back(a);
        flips.push_back(b);
        flips.push_back(c);
    } else {
        Vect l1, l2, r1, r2;
        pair <Vect, Vect>
            a = getedge({ delaune[loc.first].a, delaune[loc.first].b }),
            b = getedge({ delaune[loc.first].b, delaune[loc.first].c }),
            c = getedge({ delaune[loc.first].c, delaune[loc.first].a });
        r1 = delaune[loc.first].a;
        r2 = delaune[loc.first].b;
        pair <Vect, Vect>
            d = getedge({ delaune[loc.second].a, delaune[loc.second].b }),
            e = getedge({ delaune[loc.second].b, delaune[loc.second].c }),
            f = getedge({ delaune[loc.second].c, delaune[loc.second].a });
        if (faces[b] == loc || faces[b] == pair <ll, ll>({ loc.second, loc.first })) {
            r1 = delaune[loc.first].b;
            r2 = delaune[loc.first].c;
            swap(a, b);
        } elif(faces[c] == loc || faces[c] == pair <ll, ll>({ loc.second, loc.first })) {
            r1 = delaune[loc.first].c;
            r2 = delaune[loc.first].a;
            swap(a, c);
        }
        if (faces[e] == loc || faces[e] == pair <ll, ll>({ loc.second, loc.first })) {
            swap(d, e);
        } elif(faces[f] == loc || faces[f] == pair <ll, ll>({ loc.second, loc.first })) {
            swap(d, f);
        }
        faces.erase(a);
        if (getedge({ delaune[loc.first].a, delaune[loc.first].b }) == a) {
            l1 = delaune[loc.first].c;
        } elif(getedge({ delaune[loc.first].b, delaune[loc.first].c }) == a) {
            l1 = delaune[loc.first].a;
        } else {
            l1 = delaune[loc.first].b;
        }
        if (getedge({ delaune[loc.second].a, delaune[loc.second].b }) == a) {
            l2 = delaune[loc.second].c;
        } elif(getedge({ delaune[loc.second].b, delaune[loc.second].c }) == a) {
            l2 = delaune[loc.second].a;
        } else {
            l2 = delaune[loc.second].b;
        }

        {
            delaune.push_back({ l1, r2, p });
            delaune.push_back({ l2, r2, p });
            faces[getedge({ l1, p })] = { delaune.size() - 2, loc.first };
            faces[getedge({ r2, p })] = { delaune.size() - 2, delaune.size() - 1 };
            faces[getedge({ l2, p })] = { delaune.size() - 1, loc.second };
            faces[getedge({ r1, p })] = { loc.first, loc.second };
            delaune[loc.first] = { l1, r1, p };
            delaune[loc.second] = { l2, r1, p };
        }
        {
            pair <ll, ll> af = faces[getedge({ l1, r2 })];
            if (af.second == loc.first)
                swap(af.first, af.second);
            af.first = delaune.size() - 2;
            faces[getedge({ l1, r2 })] = af;
        }
        {
            pair <ll, ll> bf = faces[getedge({ l2, r2 })];
            if (bf.second == loc.second)
                swap(bf.first, bf.second);
            bf.first = delaune.size() - 1;
            faces[getedge({ l2, r2 })] = bf;
        }

        flips.push_back(b);
        flips.push_back(c);
        flips.push_back(e);
        flips.push_back(f);
    }
    for (auto edge : flips) {
        flip(edge);
    }
}

string numtostr(ll num) {
    string res;
    while (num != 0) {
        res += num % 10 + '0';
        num /= 10;
    }
    reverse(res.begin(), res.end());
    return res;
}

ll SZX;
ll SZY;
// Включение отображения фейковых линий в бесконечно удалённые точки
// #define IMGDEBUG

// Включение статических размеров итоговой картинки
// #define USESTATICSCALEX 1000
// #define USESTATICSCALEY 1000
void savecur(ll num) {
#ifdef USESTATICSCALEX
    SZX = USESTATICSCALEX;
#endif
#ifdef USESTATICSCALEY
    SZY = USESTATICSCALEY;
#endif
    Mat img(SZY, SZX, CV_8UC3);
#ifdef IMGDEBUG
    for (auto tr : delaune) {
        line(img, Point(tr.a.x + SZX / 2, tr.a.y + SZY / 2), Point(tr.b.x + SZX / 2, tr.b.y + SZY / 2), Scalar(0, 0, 255), 1);
        line(img, Point(tr.a.x + SZX / 2, tr.a.y + SZY / 2), Point(tr.c.x + SZX / 2, tr.c.y + SZY / 2), Scalar(0, 0, 255), 1);
        line(img, Point(tr.b.x + SZX / 2, tr.b.y + SZY / 2), Point(tr.c.x + SZX / 2, tr.c.y + SZY / 2), Scalar(0, 0, 255), 1);
    }
#else
    for (auto tr : delaune) {
        if (bannedpts[0] != tr.a && bannedpts[1] != tr.a && bannedpts[2] != tr.a) {
            if (bannedpts[0] != tr.b && bannedpts[1] != tr.b && bannedpts[2] != tr.b) {
                line(img, Point(tr.a.x + SZX / 2, tr.a.y + SZY / 2), Point(tr.b.x + SZX / 2, tr.b.y + SZY / 2), Scalar(0, 0, 255), 1);
            }
            if (bannedpts[0] != tr.c && bannedpts[1] != tr.c && bannedpts[2] != tr.c) {
                line(img, Point(tr.a.x + SZX / 2, tr.a.y + SZY / 2), Point(tr.c.x + SZX / 2, tr.c.y + SZY / 2), Scalar(0, 0, 255), 1);
            }
        }
        if (bannedpts[0] != tr.b && bannedpts[1] != tr.b && bannedpts[2] != tr.b && bannedpts[0] != tr.c && bannedpts[1] != tr.c && bannedpts[2] != tr.c) {
            line(img, Point(tr.b.x + SZX / 2, tr.b.y + SZY / 2), Point(tr.c.x + SZX / 2, tr.c.y + SZY / 2), Scalar(0, 0, 255), 1);
        }
    }
#endif
// Сохранение картинки в формате triang###.png
    {
        string filename = "triang" + numtostr(num);
        imwrite(filename + ".png", img);
    }
}

signed main() {
    delaune.push_back({ bannedpts[0], bannedpts[1], bannedpts[2] });
    faces[getedge({ bannedpts[0], bannedpts[1] })] = { 0, -1 };
    faces[getedge({ bannedpts[1], bannedpts[2] })] = { 0, -1 };
    faces[getedge({ bannedpts[2], bannedpts[0] })] = { 0, -1 };
    Vect a, b, c;
    cin >> a.x >> a.y >> b.x >> b.y >> c.x >> c.y;
    ll maxy = -INF, miny = INF, maxx = -INF, minx = INF;
    maxx = max(a.x, max(b.x, c.x));
    minx = min(a.x, min(b.x, c.x));
    maxy = max(a.y, max(b.y, c.y));
    miny = min(a.y, min(b.y, c.y));
    addpoint(a);
    addpoint(b);
    addpoint(c);
    ll i = 1;

    SZX = 50 + 2 * max(abs(maxx), abs(minx));
    SZY = 50 + 2 * max(abs(maxy), abs(miny));
    savecur(i++);
    while (cin) {
        cin >> a.x >> a.y;
        maxx = max(maxx, a.x);
        maxy = max(maxy, a.y);
        minx = min(minx, a.x);
        miny = min(miny, a.y);

        SZX = 50 + 2 * max(abs(maxx), abs(minx));
        SZY = 50 + 2 * max(abs(maxy), abs(miny));

        addpoint(a);
        savecur(i++);
    }
    return 0;
}
