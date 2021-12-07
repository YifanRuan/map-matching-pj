#pragma GCC optimize("Ofast")

#include <bits/extc++.h>
#include <bits/stdc++.h>

using namespace std;
using namespace __gnu_pbds;

using ll = long long;

constexpr int kN = 121296, kM = 1024;

int N, M;

struct Point {
    double x, y;
};

struct Line {};

struct Road {
    Line line;
    string way_type;
};

struct Trajectory {
    vector<ll> timestamps;
    vector<Point> points;
} trajectory[kM];

int main() {
    cin.tie(nullptr)->sync_with_stdio(false);
#ifndef ONLINE_JUDGE
    freopen("sample.in", "r", stdin);
    freopen("my.out", "w", stdout);
#endif

    // Input

    cin >> N;

    for (int i = 0; i < N; ++i) {
        // TODO: Input Road Sections
    }

    cin >> M;

    for (int i = 0; i < M; ++i) {
        while (true) {
            ll timestamp;
            cin >> timestamp;
            if (timestamp == i) {
                break;
            }
            trajectory[i].timestamps.push_back(timestamp);
            double x, y;
            cin >> x >> y;
            trajectory[i].points.push_back({x, y});
        }
    }

    for (int i = 1; i <= M; ++i) {
        ;
    }
    return 0;
}
