// Задачный шаблон
#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <random>
#include <vector>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#define ll long long
#define ld long double
#define elif else if
using namespace std;

// Я художник, я так вижу, мой #define - моя крепость
#define forall(i, j, breaker) for (ll i = 0, breaker = 0; i < 4 && !breaker; ++i) for (ll j = 0; j < 4 && !breaker; ++j)

class State {
public:
    vector <vector <ll>> st;
    ll SpacePos = -1;
    // Состояние хранит 16 вместо 0 для эстетики и удобства математики
    State(vector <vector <ll>> _State = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 0} }) : st(_State) {
        forall(i, j, breaker) {
            if (st[i][j] == 0) {
                SpacePos = i * 4 + j + 1;
                breaker = 1;
            }
        }
        st[(SpacePos - 1) / 4][(SpacePos - 1) % 4] = 16;
    }
    bool operator == (const State& other) const {
        return st == other.st;
    }
    bool operator != (const State& other) const {
        return st != other.st;
    }
    bool operator < (const State& other) const {
        return st < other.st;
    }
    bool operator > (const State& other) const {
        return st > other.st;
    }

    vector <ll>& operator [] (ll i) {
        return st[i];
    }
    const vector <ll>& operator [] (ll i) const {
        return st[i];
    }
    // Функция проверки состояния на решаемость
    bool checkPossibility() {
        ll InversionCount = 0;
        forall(i, j, breaker) {
            if (st[i][j] == 16)
                continue;
            for (ll ii = i; ii < 4; ++ii) {
                if (ii == i) {
                    for (ll jj = j + 1; jj < 4; ++jj) {
                        if (st[ii][jj] == 16)
                            continue;
                        if (st[i][j] > st[ii][jj])
                            ++InversionCount;
                    }
                } else {
                    for (ll jj = 0; jj < 4; ++jj) {
                        if (st[ii][jj] == 16)
                            continue;
                        if (st[i][j] > st[ii][jj])
                            ++InversionCount;
                    }
                }
            }
        }
        return !static_cast <bool>((InversionCount + (SpacePos - 1) / 4 + 1) % 2);
    }
};

// Состояние решённости
State endsy({ {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 0} });

ll func(State v) {
// Коэффициенты значимости слагаемого в оценке
#define MC 2
#define LIC 1
#define CIC 1
    ll HeuristicSum = 0;

    ll ManhattanSum = 0; // Сумма манхэттэнских расстояний до корректных позиций
    {
        forall(i, j, breaker) {
            ManhattanSum += abs((v[i][j] - 1) / 4 - i) + abs((v[i][j] - 1) % 4 - j);
        }
    }

    ll LinearInversionSum = 0; // Количество инверсий на линиях и в столбцах
    {
        for (ll i = 0; i < 4; ++i) {
            for (ll j = 0; j < 4; ++j) {
                for (ll jj = j; jj < 4; ++jj) {
                    if ((v[i][j] - 1) / 4 == i && (v[i][jj] - 1) / 4 == i && v[i][j] > v[i][jj])
                        LinearInversionSum += 2;
                    if ((v[j][i] - 1) % 4 == i && (v[jj][i] - 1) % 4 == i && v[j][i] > v[jj][i])
                        LinearInversionSum += 2;
                }
            }
        }
    }

    ll CornerInversionSum = 0; // Дополнительная сложность извлечения из углов
    {
        if (v[0][0] != 1) {
            if (v[0][1] == 2 || v[1][0] == 5)
                CornerInversionSum += 2;
        }
        if (v[0][3] != 4) {
            if (v[0][2] == 3 || v[1][3] == 6)
                CornerInversionSum += 2;
        }
        if (v[3][0] != 13) {
            if (v[2][0] == 9 || v[3][1] == 14)
                CornerInversionSum += 2;
        }
    }

    HeuristicSum = MC * ManhattanSum + LIC * LinearInversionSum + CIC * CornerInversionSum;
    return HeuristicSum;
}

// Работать будем A* в одну сторону, ибо в две влом
// Код функции свой, древний (отсюда стиль названий переменных, менять влом), но рабочий
// Комментами попытался объяснить основные моменты

// Полученная глубина состояния
map <State, ll> d;
// Прямой предок состояния на пути
map <State, State> p;
bool AStar(State start, State goal) {
    // Состояния на рассмотрении
    map <State, ll> online;
    // "Очередь" перебора с сортировкой по минимуму счёта
    set <pair <ll, State>> q;
    q.insert({ func(start), start });
    online.insert({ start, func(start) });
    d[start] = 0;
    while (q.size()) {
        pair <ll, State> v = *q.begin();
        if (v.second == goal)
            return 1;
        q.erase(v);
        for (ll i = 0; i < 4; i++) {
            State u = v.second;
            // Применение хода
            if (i == 0) {
                if ((v.second.SpacePos - 1) / 4 == 0)
                    continue;
                swap(u[(v.second.SpacePos - 1) / 4][(v.second.SpacePos - 1) % 4], u[(v.second.SpacePos - 1) / 4 - 1][(v.second.SpacePos - 1) % 4]);
                u.SpacePos -= 4;
            } elif (i == 1) {
                if ((v.second.SpacePos - 1) % 4 == 0)
                    continue;
                swap(u[(v.second.SpacePos - 1) / 4][(v.second.SpacePos - 1) % 4], u[(v.second.SpacePos - 1) / 4][(v.second.SpacePos - 1) % 4 - 1]);
                u.SpacePos -= 1;
            } elif (i == 2) {
                if ((v.second.SpacePos - 1) / 4 == 3)
                    continue;
                swap(u[(v.second.SpacePos - 1) / 4][(v.second.SpacePos - 1) % 4], u[(v.second.SpacePos - 1) / 4 + 1][(v.second.SpacePos - 1) % 4]);
                u.SpacePos += 4;
            } else {
                if ((v.second.SpacePos - 1) % 4 == 3)
                    continue;
                swap(u[(v.second.SpacePos - 1) / 4][(v.second.SpacePos - 1) % 4], u[(v.second.SpacePos - 1) / 4][(v.second.SpacePos - 1) % 4 + 1]);
                u.SpacePos += 1;
            }
            // Оценка счёта
            ll newsc = d[v.second] + 1 + func(u);
            // Постобработка, применение результатов
            if ((online.find(u) != online.end()) && newsc >= online[u])
                continue;
            if (online.find(u) == online.end() || online[u] > newsc) {
                if (online.find(u) != online.end())
                    q.erase({ online[u], u });
                p[u] = v.second;
                d[u] = d[v.second] + 1;
                online[u] = newsc;
                q.insert({ newsc, u });
            }
        }
    }
}

int main() {
    vector <vector <ll>> input(4, vector <ll>(4));
    for (ll i = 0; i < 4; ++i) {
        for (ll j = 0; j < 4; ++j) {
            cin >> input[i][j];
        }
    }
    State init(input);
    if (!init.checkPossibility()) {
        cout << "Unsolvable\n";
        return 0;
    }
    AStar(init, endsy);
    State cur = endsy, prev;
    stack <State> reverser;
    while (cur != init) {
        reverser.push(cur);
        cur = p[cur];
    }
    prev = init;
    cout << "Found solution in " << reverser.size() << endl;
    while (reverser.size()) {
        cur = reverser.top();
        reverser.pop();
// Вывод всех состояний подряд, нечитабельно и много
#ifdef SHOWALLSTATES
        forall(i, j, breaker) {
            if (cur[i][j] == 16)
                cout << 0 << " ";
            else
                cout << cur[i][j] << " ";
            if (j == 3)
                cout << endl;
    }
#endif // SHOWALLSTATES
// Вывод алгоритма для решения
#define SHOWALGORITHM
#ifdef SHOWALGORITHM
        if (cur.SpacePos == prev.SpacePos + 4)
            cout << "U";
        elif(cur.SpacePos == prev.SpacePos - 4)
            cout << "D";
        elif(cur.SpacePos == prev.SpacePos + 1)
            cout << "L";
        elif(cur.SpacePos == prev.SpacePos - 1)
            cout << "R";
        prev = cur;
#endif // SHOWALGORITHM
    }
    return 0;
}
