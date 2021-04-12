// Укреплённый задачный шаблон
#define _CRT_SECURE_NO_WARNINGS
#include "algorithm"
#include "cmath"
#include "cstddef"
#include "iomanip"
#include "iostream"
#include "list"
#include "map"
#include "memory"
#include "optional"
#include "set"
#include "stdexcept"
#include "string"
#include "queue"
#include "unordered_map"
#include "unordered_set"
#include "utility"
#include "vector"
using namespace std;
#define int long long
#define ll long long
#define INF 1e18
#define elif else if

struct STNode {
	vector <STNode*> children;
	ll termcnt;
	bool terminal;
	unordered_set <ll> nums;
	STNode() {
		children.resize(26);
		terminal = 0;
		termcnt = 0;
	}
};
// Алфавит "abcdefghijklmnopqrstuvwxyz" морально взят за константу для удобства и упрощения кода
// Можно map'ом легко расширить на более общий случай, но... Зачем?
class STrie {
public:
	// Просто техническая переменная для поиска
	ll strsize;

	STNode* root;
	STrie() {
		root = new STNode();
	}
	// Добавление строки в бор
	void add(const string& newstring) {
		strsize = max(strsize, static_cast <ll> (newstring.size()));
		STNode* cur = root;
		ll num = cur->termcnt++;
		for (ll i = 0; i < newstring.size(); ++i) {
			if (cur->children[newstring[i] - 'a'] == nullptr) {
				cur->children[newstring[i] - 'a'] = new STNode();
			}
			cur = cur->children[newstring[i] - 'a'];
			cur->nums.insert(num + i);
			cur->termcnt++;
		}
		cur->terminal = true;
	}

	// Возвращает множество подходящих первых элементов
	// (Куча костылей и разборов случаев, читать на свой страх и риск, писал после ночи без сна)
	unordered_set <ll> tryfind(const string& regex) {
		STNode* cur = root;
		unordered_set <ll> possibles;
		possibles.insert(-1);
		for (ll i = 0; i < regex.size(); ++i) {
			if (possibles.empty()) {
				break;
			} elif(regex[i] == '?' && (*possibles.begin()) == -1) {
				continue;
			} elif(regex[i] == '?' && (*possibles.begin()) != -2 && cur != root) {
				unordered_set <ll> newpossibles;
				if (cur->nums.size() < possibles.size()) {
					for (const auto& it : cur->nums) {
						if (possibles.find(it - i + 1) != possibles.end()) {
							newpossibles.insert(it - i + 1);
						}
					}
				} else {
					for (const auto& it : possibles) {
						if (cur->nums.find(it + i) != cur->nums.end()) {
							newpossibles.insert(it);
						}
					}
				}
				possibles = newpossibles;
				cur = root;
			} elif(regex[i] == '?' && (*possibles.begin()) != -2) {
				continue;
			} elif(regex[i] == '?') {
				possibles.clear();
				for (const auto& it : cur->nums) {
					if (it - i + 1 >= 0)
						possibles.insert(it - i + 1);
				}
				cur = root;
			} elif (*possibles.begin() == -1) {
				possibles.clear();
				possibles.insert(-2);
				cur = cur->children[regex[i] - 'a'];
				if (cur == nullptr) {
					possibles.clear();
					break;
				}
			} else {
				cur = cur->children[regex[i] - 'a'];
				if (cur == nullptr) {
					possibles.clear();
					break;
				}
			}
		}
		if (!cur || possibles.empty()) {
			return possibles;
		}
		if (cur != root && *possibles.begin() >= 0) {
			unordered_set <ll> newpossibles;
			if (cur->nums.size() != possibles.size()) {
				for (const auto& it : cur->nums) {
					if (possibles.find(it - regex.size() + 1) != possibles.end()) {
						newpossibles.insert(it - regex.size() + 1);
					}
				}
			} else {
				for (const auto& it : possibles) {
					if (cur->nums.find(it + regex.size() - 1) != cur->nums.end()) {
						newpossibles.insert(it);
					}
				}
			}
			possibles = newpossibles;
		} elif(cur != root) {
			possibles.clear();
			for (const auto& it : cur->nums) {
				if (it - regex.size() + 1 >= 0)
					possibles.insert(it - regex.size() + 1);
			}
		} elif (*possibles.begin() == -1) {
			possibles.clear();
			for (ll i = 0; i < strsize - regex.size() + 1; ++i)
				possibles.insert(i);
		}
		return possibles;
	}
};

// Ассимптотика O(n^2) гарантированная на предподсчёт, O(p + k) ожидаемая на запрос
// (или что-то ассимптотически близкое, точнее сказать не могу).
// Вместо суфбора в теории можно было бы использовать суфдерево, Укконеном строится за O(n),
// но не было времени довести до ума и реализовать.

signed main() {
	STrie trie;
	string input;
	cin >> input;
	for (ll i = 0; i < input.size(); ++i) {
		trie.add(input.substr(i));
	}
	string regex;
	while (1) {
		cin >> regex;
		unordered_set <ll> ans = trie.tryfind(regex);
		for (const auto& it : ans) {
			if (it + regex.size() - 1 < input.size())
				cout << it << " ";
		}
		cout << endl;
	}
	return 0;
}
