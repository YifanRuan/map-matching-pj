#pragma GCC optimize("Ofast")

#include <bits/stdc++.h>

#define pb push_back

using namespace std;

using ll = long long;
using ld = double;
using pii = pair<int, int>;
using pid = pair<int, ld>;
using pdd = pair<ld, ld>;

// Definitions of constants

constexpr ld kPI =
    3.14159265358979323846264338327950288419716939937510582097494459230;
constexpr ld kEps = 1e-9;
constexpr ld kEarthRadius = 6371000.0;
constexpr ld kR = 50; // r: candidates within radius r, corelated to kGap
constexpr int kMaxTrajectoryLength = 2080;
constexpr int kMaxVertexNum = 58000;
constexpr int kMaxN = 130000;
constexpr ld kMaxLatitude = 31.4777782, kMaxLongitude = 122.0933554,
             kMinLatitude = 30.6405552, kMinLongitude = 121.0037503;
constexpr ld kInf = 1e18;
// constexpr ld kMaxSpeed[8] = {22, 23, 24, 25, 26, 27, 28, 29};
constexpr ld kMaxSpeed[8] = {5, 6, 7, 8, 9, 10, 11, 12};
constexpr ld kLengthPerRad = 111226.29021121707545;
constexpr ld kBeta[31] = {0,           0.49037673,  0.82918373,  1.24364564,
                          1.67079581,  2.00719298,  2.42513007,  2.81248831,
                          3.15745473,  3.52645392,  4.09511775,  4.67319795,
                          5.41088180,  6.47666590,  6.29010734,  7.80752112,
                          8.09074504,  8.08550528,  9.09405065,  11.09090603,
                          11.87752824, 12.55107715, 15.82820829, 17.69496773,
                          18.07655652, 19.63438911, 25.40832185, 23.76001877,
                          28.43289797, 32.21683062, 34.56991141};
constexpr ld kSigmaZ = 12.00;
constexpr ld kSpeed = 60;
vector<int> kGridIndex[2200][2200];
ld kDist[kMaxVertexNum];
vector<set<int>> kEdgeList;
map<pii, pdd> kShortestPath;
int kMaxLongitudeGrid, kMaxLatitudeGrid;
ld kGridSize;
vector<pid> kAdj[kMaxVertexNum];

// util functions

ld ToRed(ld deg) { return deg * 2 * kPI / 360; }

ld CalcGPSDistance(ld lon1, ld lat1, ld lon2, ld lat2) {
    ld dLat = ToRed(lat2 - lat1);
    ld dLon = ToRed(lon2 - lon1);
    ld a = sin(dLat / 2) * sin(dLat / 2) +
           cos(ToRed(lat1)) * cos(ToRed(lat2)) * sin(dLon / 2) * sin(dLon / 2);
    ld c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return kEarthRadius * c;
}

ld CalcGreatCircleDistance(ld lon1, ld lat1, ld lon2, ld lat2) {
    ld dLat = lat1 - lat2;
    ld dLon = (lon2 - lon1) * cos(lat1 * kPI / 180);
    return kLengthPerRad * sqrt(dLat * dLat + dLon * dLon);
}

ld Distance(ld lon1, ld lat1, ld lon2, ld lat2) {
    // return CalcGreatCircleDistance(lon1, lat1, lon2, lat2);
    return CalcGPSDistance(lon1, lat1, lon2, lat2);
}

ld GetBeta(int sample_period) { return kBeta[min(sample_period, 30)]; }

// forward declaration

struct Trajectory;
void FindMatchedSequence(Trajectory &traj);

// Data structures

struct Point {
    ld x, y;
    ll timestamp;

    Point() : x(0), y(0), timestamp(0) {}
    Point(ld x, ld y, ll timestamp) : x(x), y(y), timestamp(timestamp) {}
    Point(ld x, ld y) : x(x), y(y), timestamp(0) {}

    ld CalculateDistanceInMeters(const Point &other) const {
        return Distance(x, y, other.x, other.y);
    }

    ld CalcDistanceToSegment(const Point &a, const Point &b, ld &px,
                             ld &py) const {
        ld d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
        ld t = (x - a.x) * (b.x - a.x) + (y - a.y) * (b.y - a.y) / d;
        if (t < 0) {
            t = 0;
        } else if (t > 1) {
            t = 1;
        }
        px = a.x + t * (b.x - a.x);
        py = a.y + t * (b.y - a.y);
        return CalculateDistanceInMeters(Point(px, py));
    }
};

struct Edge {
    int id, start, end, way_type, c;
    string way_string;
    vector<Point> p;
    vector<ld> dist_from_start;

    ld length() { return dist_from_start[c - 1]; }
} kEdge[kMaxN];

struct Trajectory {
    int n;
    Point p[kMaxTrajectoryLength];

    Trajectory() : n(0) {}

    void insert(Point poi) { p[n++] = poi; }

    void GetCandidates() {
        kEdgeList.clear();
        for (int i = 0; i < n; ++i) {
            int row = (p[i].y - kMinLatitude) / kGridSize;
            int col = (p[i].x - kMinLongitude) / kGridSize;
            set<int> res;
            for (int j = row - 1; j <= row + 1; ++j) {
                for (int k = col - 1; k <= col + 1; ++k) {
                    for (auto it : kGridIndex[j][k]) {
                        res.insert(it);
                    }
                }
            }
            kEdgeList.pb(res);
        }
    }

    void STMatching() {
        // Precompute
        GetCandidates();

        FindMatchedSequence(*this);
    }
};

struct Probability {
    int eid, pre;
    ld score, dist_from_start;
    Probability() {}
    Probability(int eid, ld score, int pre, ld dist_from_start)
        : eid(eid), score(score), pre(pre), dist_from_start(dist_from_start) {}
};

void AddEdge(int from, int to, ld cost) { kAdj[from].pb({to, cost}); }

ld Distance(Point z, Edge e, ld &point_dist_from_start) {
    ld mn = kInf;
    for (int i = 1; i < e.c; i++) {
        Point nex;
        ld tmp = z.CalcDistanceToSegment(e.p[i - 1], e.p[i], nex.x, nex.y);
        if (tmp < mn) {
            mn = tmp;
            point_dist_from_start = e.dist_from_start[i - 1] +
                                    e.p[i - 1].CalculateDistanceInMeters(nex);
        }
    }
    return mn;
}

ld N(ld dis) {
    ld tmp = (dis / kSigmaZ);
    return exp(-0.5 * tmp * tmp) / (sqrt(2 * kPI) * kSigmaZ);
}

ld N(ld dis, Edge e) {
    ld tmp = (dis / kSigmaZ);
    return exp(-0.5 * tmp * tmp) / (sqrt(2 * kPI) * kSigmaZ) *
           kMaxSpeed[e.way_type];
}

ld Dijkstra(int s, int t, ld delta_t) {
    using Node = pair<ld, int>;
    for (int i = 0; i < kMaxVertexNum; ++i) {
        kDist[i] = kInf;
    }
    kDist[s] = 0;
    priority_queue<Node, vector<Node>, greater<Node>> pq;
    pq.push(make_pair(0, s));
    while (!pq.empty()) {
        auto u = pq.top();
        pq.pop();
        if (u.first > delta_t * kSpeed)
            return kInf;
        if (u.first > kDist[u.second])
            continue;
        if (u.second == t)
            break;
        for (auto e : kAdj[u.second]) {
            if (kDist[e.first] > u.first + e.second) {
                kDist[e.first] = u.first + e.second;
                pq.push(make_pair(kDist[e.first], e.first));
            }
        }
    }
    return kDist[t];
}

// 获得这一轮概率转移分数最大值的下标
int GetMaxProbIndex(vector<Probability> &row) {
    int ret = -1;
    ld possibility = -1;
    for (int i = 0; i < row.size(); i++) {
        if (possibility < row[i].score) {
            possibility = row.at(i).score;
            ret = i;
        }
    }
    return ret;
}

// 对 pt1, pt2 两个点构成的线段标记网格，id 为边序号
void GetCrossCells(int id, vector<Point>::iterator pt1,
                   vector<Point>::iterator pt2) {
    ld x1 = pt1->x - kMinLongitude, y1 = pt1->y - kMinLatitude,
       x2 = pt2->x - kMinLongitude, y2 = pt2->y - kMinLatitude;
    int row = y1 / kGridSize, nex_row = y2 / kGridSize, col = x1 / kGridSize,
        nex_col = x2 / kGridSize;

    // 同一行
    if (row == nex_row) {
        for (int i = min(col, nex_col); i <= max(col, nex_col); ++i) {
            if (row >= kMaxLatitudeGrid || row < 0 ||
                col >= kMaxLongitudeGrid || col < 0 ||
                (kGridIndex[row][col].size() > 0 &&
                 kGridIndex[row][col].back() == id))
                continue;
            kGridIndex[row][i].pb(id);
        }
        return;
    }

    // 同一列
    if (col == nex_col) {
        for (int i = min(row, nex_row); i < max(row, nex_row); ++i) {
            if (row >= kMaxLatitudeGrid || row < 0 ||
                col >= kMaxLongitudeGrid || col < 0 ||
                (kGridIndex[row][col].size() > 0 &&
                 kGridIndex[row][col].back() == id))
                continue;
            kGridIndex[i][col].pb(id);
        }
        return;
    }
    ld dy = y2 - y1, dx = x2 - x1, c = x1 * y2 - x2 * y1;
    // 加与坐标轴相交的点
    vector<pdd> points;
    for (int i = min(row, nex_row); i <= max(row, nex_row); ++i) {
        points.pb(make_pair((c + dx * i * kGridSize) / dy, i * kGridSize));
    }
    for (int i = min(col, nex_col); i <= max(col, nex_col); ++i) {
        points.pb(make_pair(i * kGridSize, (-c + dy * i * kGridSize) / dx));
    }
    points.pb(make_pair(x1, y1));
    points.pb(make_pair(x1, y1));
    sort(points.begin(), points.end());

    // 枚举点间直线
    for (auto it = ++points.begin(); it != points.end(); ++it) {
        auto pre_it = --it;
        ++it;
        int row = min((int)(pre_it->second / kGridSize + kEps),
                      (int)(it->second / kGridSize + kEps));
        int col = min((int)(pre_it->first / kGridSize + kEps),
                      (int)(it->first / kGridSize + kEps));
        if (row >= kMaxLatitudeGrid || row < 0 || col >= kMaxLongitudeGrid ||
            col < 0 ||
            (kGridIndex[row][col].size() > 0 &&
             kGridIndex[row][col].back() == id))
            continue;
        kGridIndex[row][col].pb(id);
    }
}

void GenerateGridIndex(int n) {
    // 初始化网格索引
    kMaxLongitudeGrid =
        int(Distance(kMaxLongitude, kMinLatitude, kMinLongitude, kMinLatitude) /
            kR) +
        2;
    kMaxLatitudeGrid =
        int((kMaxLatitude - kMinLatitude) / (kMaxLongitude - kMinLongitude) *
            ld(kMaxLongitudeGrid)) +
        2;
    kGridSize = (kMaxLongitude - kMinLongitude) / ld(kMaxLongitudeGrid);

    // 枚举线段
    for (int i = 0; i < n; ++i) {
        auto nex_p = kEdge[i].p.begin();
        auto p = nex_p++;
        while (nex_p != kEdge[i].p.end()) {
            GetCrossCells(i, p, nex_p);
            p = nex_p++; // TODO
        }
    }
}

void FindMatchedSequence(Trajectory &traj) {
    int sample_period =
        (traj.p[traj.n - 1].timestamp - traj.p[0].timestamp) / (traj.n - 1);
    ld beta = GetBeta(sample_period);
    vector<vector<Probability>> probability_matrix;
    bool can_propagate = false;
    int ct = traj.n; // last transferred, or traj.n if can not transfer

    // iterate traj points p[i]
    for (int i = 0; i < traj.n; ++i) {
        ld dist_between_two_traj_points, delta_t = -1;
        if (ct != traj.n) {
            delta_t = traj.p[i].timestamp - traj.p[ct].timestamp;
            dist_between_two_traj_points =
                traj.p[i].CalculateDistanceInMeters(traj.p[ct]);
        }
        ld mx_prob = -kInf;
        vector<Probability> probs;
        ld emission_prob;

        // 枚举候选边
        for (auto cand : kEdgeList[i]) {
            int pre = -1;
            ld dist_from_start = 0,
               dist_from_point =
                   Distance(traj.p[i], kEdge[cand], dist_from_start);
            emission_prob = N(dist_from_point);

            // 能从上一轮转移
            if (can_propagate) {
                ld mx = -kInf;
                int c = 0; // 当前候选转移点（上一轮）

                // 枚举上一轮的概率
                for (auto sp : probability_matrix.back()) {
                    ld pre_dist_from_start = sp.dist_from_start;
                    ld pre_dist_to_end =
                        kEdge[sp.eid].length() - pre_dist_from_start;
                    ld road_distance, on_graph_distance;
                    // 这里进行了一个细节处理，如果两个点在同一条边上，则直接计算他们到边起点的距离差
                    // 而且并没有对反向的情况进行处理，这样得分会更高
                    if (cand == sp.eid) {
                        on_graph_distance =
                            fabs(dist_from_start - pre_dist_from_start);
                    } else {
                        pii ind =
                            make_pair(kEdge[sp.eid].end, kEdge[cand].start);
                        // 剪枝
                        if (kShortestPath.find(ind) != kShortestPath.end() &&
                            kShortestPath[ind].first < kInf) {
                            kShortestPath[ind].first <= kSpeed *delta_t
                                ? road_distance = kShortestPath[ind].first
                                : road_distance = kInf;
                        } else {
                            if (kShortestPath.find(ind) !=
                                    kShortestPath.end() &&
                                delta_t + kEps <= kShortestPath[ind].second)
                                road_distance = kInf;
                            else {
                                road_distance =
                                    Dijkstra(kEdge[sp.eid].end,
                                             kEdge[cand].start, delta_t);
                                kShortestPath[ind] =
                                    make_pair(road_distance, delta_t);
                            }
                        }
                        // 计算距离
                        on_graph_distance =
                            road_distance + dist_from_start + pre_dist_to_end;
                    }

                    // 计算转移概率
                    ld trans_prob = exp(-fabs(dist_between_two_traj_points -
                                              on_graph_distance) /
                                        beta) /
                                    beta;
                    ld alt = sp.score * trans_prob;
                    if (mx < alt) {
                        mx = alt;
                        pre = c;
                    }
                    ++c;
                }
                // end of iterate previous probs

                // 计算得分
                emission_prob *= mx;
            }
            // end of can propagate

            probs.pb(Probability(cand, emission_prob, pre, dist_from_start));
            if (mx_prob < emission_prob)
                mx_prob = emission_prob;
        }
        // end of iterrate candidate edges

        ct = i;

        // 归一化
        for (int j = 0; j < probs.size(); ++j)
            probs[j].score /= mx_prob;
        probability_matrix.pb(probs);

        // 可能候选边是零
        if (probs.size() == 0) {
            can_propagate = false;
            ct = traj.n;
        } else {
            can_propagate = true;
        }
    }
    // end of iterate traj points p[i]

    // 反向获取路径
    vector<int> r_list;
    int c = 0, cid = GetMaxProbIndex(probability_matrix.back());
    for (int j = probability_matrix.size() - 1; j >= 0; j--) {
        if (cid != -1) {
            c = probability_matrix[j][cid].eid;
            r_list.pb(c);
            cid = probability_matrix[j][cid].pre;
        } else {
            r_list.pb(c);
            if (j)
                cid = GetMaxProbIndex(probability_matrix[j - 1]);
        }
    }
    reverse(r_list.begin(), r_list.end());
    for (auto it : r_list) {
        cout << it << " ";
    }
    cout << "\n";
}

int main() {
    // Deal with input parameters and MACROs
    cin.tie(nullptr)->sync_with_stdio(false);
#ifndef ONLINE_JUDGE
    freopen("sample.in", "r", stdin);
    freopen("my2.out", "w", stdout);
#endif
    int kN, kM;
    cin >> kN;

    // Preprocessing
    for (int i = 0; i < kN; ++i) {
        cin >> kEdge[i].id >> kEdge[i].start >> kEdge[i].end >>
            kEdge[i].way_string >> kEdge[i].way_type >> kEdge[i].c;
        for (int j = 0; j < kEdge[i].c; ++j) {
            ld x, y;
            cin >> y >> x;
            kEdge[i].p.pb(Point(x, y));
            if (j)
                kEdge[i].dist_from_start.pb(
                    kEdge[i].dist_from_start[j - 1] +
                    kEdge[i].p[j - 1].CalculateDistanceInMeters(kEdge[i].p[j]));
            else
                kEdge[i].dist_from_start.pb(0);
        }
        AddEdge(kEdge[i].start, kEdge[i].end, kEdge[i].length());
    }
    GenerateGridIndex(kN);

    // Main
    cin >> kM;
    cout << kM << "\n";
    for (int i = 0; i < kM; ++i) {
        Trajectory traj;
        while (true) {
            ll timestamp;
            cin >> timestamp;
            if (timestamp == i)
                break;
            ld x, y;
            cin >> y >> x;
            traj.insert(Point(x, y, timestamp));
        }
        traj.STMatching();
    }
    return 0;
}
