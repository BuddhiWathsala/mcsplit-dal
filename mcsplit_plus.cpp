
#include "graph.h"
#include "nausparse.h"

#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>

#include <argp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define VSCORE
//#define DEBUG
#define Best
using std::vector;
using std::cout;
using std::endl;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}


enum Heuristic {
    min_max, min_product
};

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format\vHEURISTIC can be min_max or min_product";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
static struct argp_option options[] = {
        {"quiet",                'q', 0,         0, "Quiet output"},
        {"verbose",              'v', 0,         0, "Verbose output"},
        {"dimacs",               'd', 0,         0, "Read DIMACS format"},
        {"lad",                  'l', 0,         0, "Read LAD format"},
        {"connected",            'c', 0,         0, "Solve max common CONNECTED subgraph problem"},
        {"directed",             'i', 0,         0, "Use directed graphs"},
        {"labelled",             'a', 0,         0, "Use edge and vertex labels"},
        {"vertex-labelled-only", 'x', 0,         0, "Use vertex labels, but not edge labels"},
        {"big-first",            'b', 0,         0, "First try to find an induced subgraph isomorphism, then decrement the target size"},
        {"timeout",              't', "timeout", 0, "Specify a timeout (seconds)"},
        {0}
};

static struct {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool connected;
    bool directed;
    bool edge_labelled;
    bool vertex_labelled;
    bool big_first;
    Heuristic heuristic;
    char *filename1;
    char *filename2;
    int timeout;
    int arg_num;
} arguments;

static std::atomic<bool> abort_due_to_timeout;
std::chrono::time_point<std::chrono::steady_clock> START;
vector<VtxPair> incumbent;

void set_default_arguments() {
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.dimacs = false;
    arguments.lad = false;
    arguments.connected = false;
    arguments.directed = false;
    arguments.edge_labelled = false;
    arguments.vertex_labelled = false;
    arguments.big_first = false;
    arguments.filename1 = NULL;
    arguments.filename2 = NULL;
    arguments.timeout = 0;
    arguments.arg_num = 0;
}

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'd':
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            break;
        case 'l':
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            break;
        case 'q':
            arguments.quiet = true;
            break;
        case 'v':
            arguments.verbose = true;
            break;
        case 'c':
            if (arguments.directed)
                fail("The connected and directed options can't be used together.");
            arguments.connected = true;
            break;
        case 'i':
            if (arguments.connected)
                fail("The connected and directed options can't be used together.");
            arguments.directed = true;
            break;
        case 'a':
            if (arguments.vertex_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.edge_labelled = true;
            arguments.vertex_labelled = true;
            break;
        case 'x':
            if (arguments.edge_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.vertex_labelled = true;
            break;
        case 'b':
            arguments.big_first = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                if (std::string(arg) == "min_max")
                    arguments.heuristic = min_max;
                else if (std::string(arg) == "min_product")
                    arguments.heuristic = min_product;
                else
                    fail("Unknown heuristic (try min_max or min_product)");
            } else if (arguments.arg_num == 1) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 2) {
                arguments.filename2 = arg;
            } else {
                argp_usage(state);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

/*******************************************************************************
                                     Stats
*******************************************************************************/

unsigned long long nodes_mcsplit = 0;
unsigned long long nodes_mcsolver = 0;
unsigned long long cut_branches_mcsplit = 0;
unsigned long long cut_branches_mcsolver = 0;
unsigned long long best_count = 0;
unsigned long long calls_for_optimal = 0;
unsigned long long conflicts = 0;
std::chrono::duration<double> duration;
size_t min_graph_size;
int dl = 0;
int Maxnum = 0;
int M = 1;
int num = 0;
clock_t start;

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

using gtype = double;
const int short_memory_threshold = 1e5;
const int long_memory_threshold = 1e9;

struct Bidomain {
    int l, r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;

    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent) :
            l(l),
            r(r),
            left_len(left_len),
            right_len(right_len),
            is_adjacent(is_adjacent) {};
};

void show(const vector<VtxPair> &current, const vector<Bidomain> &domains,
          const vector<int> &left, const vector<int> &right) {
    cout << "Nodes: " << nodes_mcsplit << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i = 0; i < current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    for (unsigned int i = 0; i < domains.size(); i++) {
        struct Bidomain bd = domains[i];
        cout << "Left  ";
        for (int j = 0; j < bd.left_len; j++)
            cout << left[bd.l + j] << " ";
        cout << std::endl;
        cout << "Right  ";
        for (int j = 0; j < bd.right_len; j++)
            cout << right[bd.r + j] << " ";
        cout << std::endl;
    }
    cout << "\n" << std::endl;
}

bool check_sol(const Graph &g0, const Graph &g1, const vector<VtxPair> &solution) {
    return true;
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i = 0; i < solution.size(); i++) {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j = i + 1; j < solution.size(); j++) {
            struct VtxPair p1 = solution[j];
            if (g0.adjmat[p0.v][p1.v] != g1.adjmat[p0.w][p1.w])
                return false;
        }
    }
    return true;
}

int calc_bound(const vector<Bidomain> &domains) {
    int bound = 0;
    for (const Bidomain &bd: domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

int selectV_index(const vector<int> &arr, const vector<gtype> &lgrade, int start_idx, int len) {
    int idx = -1;
    gtype max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for (int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        if (lgrade[vtx] > max_g) {
            idx = i;
            best_vtx = vtx;
            max_g = lgrade[vtx];
        } else if (lgrade[vtx] == max_g) {
            if (vtx < best_vtx) {
                idx = i;
                best_vtx = vtx;
            }
        }
    }
    return idx;
}

int select_bidomain(const vector<Bidomain> &domains, const vector<int> &left, const vector<gtype> &lgrade,
                    int current_matching_size) {
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int tie_breaker;
    unsigned int i;
    int len;
    int best = -1;
    for (i = 0; i < domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size > 0 && !bd.is_adjacent) continue;
        len = arguments.heuristic == min_max ?
              std::max(bd.left_len, bd.right_len) :
              bd.left_len * bd.right_len;//最多有这么多种连线
        if (len < min_size) {//优先计算可能连线少的情况,当domain中可能配对数最小时，选择domain中含有最大度的点那个domain
            min_size = len;//此时进行修改将度最大的点进行优化
            min_tie_breaker = left[bd.l + selectV_index(left, lgrade, bd.l, bd.left_len)];//序号越小度越大
            best = i;
        } else if (len == min_size) {
            tie_breaker = left[bd.l + selectV_index(left, lgrade, bd.l, bd.left_len)];
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

// Returns length of left half of array
int partition(vector<int> &all_vv, int start, int len, const vector<unsigned int> &adjrow) {
    int i = 0;
    for (int j = 0; j < len; j++) {
        if (adjrow[all_vv[start + j]]) {
            std::swap(all_vv[start + i], all_vv[start + j]);
            i++;
        }
    }
    return i;
}

int remove_matched_vertex(vector<int> &arr, int start, int len, const vector<int> &matched) {
    int p = 0;
    for (int i = 0; i < len; i++) {
        if (!matched[arr[start + i]]) {
            std::swap(arr[start + i], arr[start + p]);
            p++;
        }
    }
    return p;
}

vector<Bidomain> rewardfeed_RL(const vector<Bidomain> &d, int bd_idx, vector<VtxPair> &current, vector<int> &left,
                               vector<int> &right, vector<gtype> &lgrade, vector<gtype> &rgrade, const Graph &g0,
                               const Graph &g1, int v, int w,
                               bool multiway) {//每个domain均分成2个new_domain分别是与当前对应点相连或者不相连
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    current.push_back(VtxPair(v, w));
    //unsigned int old_bound = current.size() + calc_bound(d)+1;//这里的domain是已经去掉选了的点，但是该点还没纳入current所以+1
    //unsigned int new_bound(0);
    int l, r = -1;
    unsigned int i;
    int temp = 0, total = 0;
    for (const Bidomain &old_bd: d) {
        l = old_bd.l;
        r = old_bd.r;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);//将与选出来的顶点相连顶点数返回，
        int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);//并且将相连顶点依次交换到数组前面
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        //这里传递v,w选取的下标到bd_idx，j用来计算当前所分domain的下标，这样就能将bound更精确地计算
        temp = std::min(old_bd.left_len, old_bd.right_len) - std::min(left_len, right_len) -
               std::min(left_len_noedge, right_len_noedge);
        total += temp;

#ifdef DEBUG

        printf("adj=%d ,noadj=%d ,old=%d ,temp=%d \n",std::min(left_len,right_len),
                  std::min(left_len_noedge,right_len_noedge),std::min(old_bd.left_len,old_bd.right_len),temp);
           cout<<"j="<<j<< "  idx="<<bd_idx<<endl;
           cout<<"gl="<<lgrade[v]<<" gr="<<rgrade[w]<<endl;
#endif
        if (left_len_noedge && right_len_noedge)//new_domain存在的条件是同一domain内需要存在与该对应点同时相连，或者同时不相连的点的点
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        //Question:when the value of 'left_len_noedge' will be 0,
        //what does it mean if it happens?
        if (multiway && left_len && right_len) {//不与顶点相连的顶点
            auto &adjrow_v = g0.adjmat[v];
            auto &adjrow_w = g1.adjmat[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin + left_len, [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(r_begin, r_begin + right_len, [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l < l_top && r < r_top) {
                unsigned int left_label = adjrow_v[left[l]];
                unsigned int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l < l_top && adjrow_v[left[l]] == left_label);
                    do { r++; } while (r < r_top && adjrow_w[right[r]] == left_label);
                    new_d.push_back({lmin, rmin, l - lmin, r - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});//与顶点相连的顶点 标志为Ture
        }
    }
    if (total > 0) {
        conflicts++;
        lgrade[v] += total;
        if (lgrade[v] > short_memory_threshold) {
            for (i = 0; i < lgrade.size(); i++)
                lgrade[i] = lgrade[i] / 2;
        }
        rgrade[w] += total;
        if (rgrade[w] > short_memory_threshold) {
            for (i = 0; i < rgrade.size(); i++)
                rgrade[i] = rgrade[i] / 2;
        }
    }
#ifdef DEBUG
    cout<<"new domains are "<<endl;
     for (const Bidomain &testd : new_d){
         l = testd.l;
         r = testd.r;
         for(j=0;j<testd.left_len;j++)
             cout<<left[l+j] <<" ";
         cout<<" ; ";
         for(j=0;j<testd.right_len;j++)
             cout<<right[r+j]<<" " ;
         cout<<endl;
     }
#endif
    return new_d;
}


// multiway is for directed and/or labelled graphs
vector<Bidomain> rewardfeed(const vector<Bidomain> &d, int bd_idx, vector<VtxPair> &current, vector<int> &g0_matched,
                            vector<int> &g1_matched,
                            vector<int> &left, vector<int> &right, vector<gtype> &lgrade, vector<gtype> &rgrade,
                            const Graph &g0, const Graph &g1, int v, int w,
                            bool multiway) {//每个domain均分成2个new_domain分别是与当前对应点相连或者不相连
//    assert(g0_matched[v] == 0);
//    assert(g1_matched[w] == 0);
    current.push_back(VtxPair(v, w));
    g0_matched[v] = 1;
    g1_matched[w] = 1;

    int leaves_match_size = 0, v_leaf, w_leaf;
    for (unsigned int i = 0, j = 0; i < g0.leaves[v].size() && j < g1.leaves[w].size();) {
        if (g0.leaves[v][i].first < g1.leaves[w][j].first) i++;
        else if (g0.leaves[v][i].first > g1.leaves[w][j].first) j++;
        else {
            const vector<int> &leaf0 = g0.leaves[v][i].second;
            const vector<int> &leaf1 = g1.leaves[w][j].second;
            for (unsigned int p = 0, q = 0; p < leaf0.size() && q < leaf1.size();) {
                if (g0_matched[leaf0[p]]) p++;
                else if (g1_matched[leaf1[q]]) q++;
                else {
                    v_leaf = leaf0[p], w_leaf = leaf1[q];
                    p++, q++;
                    current.push_back(VtxPair(v_leaf, w_leaf));
                    g0_matched[v_leaf] = 1;
                    g1_matched[w_leaf] = 1;
                    leaves_match_size++;
                }
            }
            i++, j++;
        }
    }


    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    //unsigned int old_bound = current.size() + calc_bound(d)+1;//这里的domain是已经去掉选了的点，但是该点还没纳入current所以+1
    //unsigned int new_bound(0);
    int l, r = -1;
    int temp = 0, total = 0;
    int unmatched_left_len, unmatched_right_len;
    for (const Bidomain &old_bd: d) {
        l = old_bd.l;
        r = old_bd.r;
        if (leaves_match_size > 0 && old_bd.is_adjacent == false) {
            unmatched_left_len = remove_matched_vertex(left, l, old_bd.left_len, g0_matched);
            unmatched_right_len = remove_matched_vertex(right, r, old_bd.right_len, g1_matched);
        } else {
            unmatched_left_len = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, unmatched_left_len, g0.adjmat[v]);//将与选出来的顶点相连顶点数返回，
        int right_len = partition(right, r, unmatched_right_len, g1.adjmat[w]);//并且将相连顶点依次交换到数组前面
        int left_len_noedge = unmatched_left_len - left_len;
        int right_len_noedge = unmatched_right_len - right_len;
//这里传递v,w选取的下标到bd_idx，j用来计算当前所分domain的下标，这样就能将bound更精确地计算
        temp = std::min(old_bd.left_len, old_bd.right_len) - std::min(left_len, right_len) -
               std::min(left_len_noedge, right_len_noedge);
        total += temp;

#ifdef DEBUG

        printf("adj=%d ,noadj=%d ,old=%d ,temp=%d \n",std::min(left_len,right_len),
               std::min(left_len_noedge,right_len_noedge),std::min(old_bd.left_len,old_bd.right_len),temp);
        cout<<"j="<<j<< "  idx="<<bd_idx<<endl;
        cout<<"gl="<<lgrade[v]<<" gr="<<rgrade[w]<<endl;
#endif
        if (left_len_noedge && right_len_noedge)//new_domain存在的条件是同一domain内需要存在与该对应点同时相连，或者同时不相连的点的点
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {//不与顶点相连的顶点
            auto &adjrow_v = g0.adjmat[v];
            auto &adjrow_w = g1.adjmat[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin + left_len, [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(r_begin, r_begin + right_len, [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l < l_top && r < r_top) {
                unsigned int left_label = adjrow_v[left[l]];
                unsigned int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l < l_top && adjrow_v[left[l]] == left_label);
                    do { r++; } while (r < r_top && adjrow_w[right[r]] == left_label);
                    new_d.push_back({lmin, rmin, l - lmin, r - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});//与顶点相连的顶点 标志为Ture
        }
    }
    int domgrade = new_d.size();
    total = total + domgrade;
    if (total > 0) {
        conflicts++;

        lgrade[v] += total;
        float UB = 0.2;
        lgrade[v] = 0.97 * lgrade[v] + UB; // CHANGE
        rgrade[w] += total;

        if (lgrade[v] > short_memory_threshold)
            for (int i = 0; i < g0.n; i++)
                lgrade[i] = lgrade[i] / 2;
        if (rgrade[w] > long_memory_threshold)
            for (int i = 0; i < g1.n; i++)
                rgrade[i] = rgrade[i] / 2;
    }
#ifdef DEBUG
    cout<<"new domains are "<<endl;
  for (const Bidomain &testd : new_d){
      l = testd.l;
      r = testd.r;
      for(j=0;j<testd.left_len;j++)
          cout<<left[l+j] <<" ";
      cout<<" ; ";
      for(j=0;j<testd.right_len;j++)
          cout<<right[r+j]<<" " ;
      cout<<endl;
  }
#endif
    return new_d;
}

int selectW_index(const vector<int> &arr, const vector<gtype> &rgrade, int start_idx, int len,
                  const vector<int> &wselected) {
    int idx = -1;
    gtype max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for (int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        if (wselected[vtx] == 0) {
            if (rgrade[vtx] > max_g) {
                idx = i;
                best_vtx = vtx;
                max_g = rgrade[vtx];
            } else if (rgrade[vtx] == max_g) {
                if (vtx < best_vtx) {
                    idx = i;
                    best_vtx = vtx;
                }
            }
        }
    }
    return idx;
}

void remove_vtx_from_array(vector<int> &arr, int start_idx, int &len, int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

void remove_bidomain(vector<Bidomain> &domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

void solve(const Graph &g0, const Graph &g1, vector<gtype> &V, vector<gtype> &lgrade,
           vector<gtype> &rgrade, vector<vector<gtype>> &Q, vector<VtxPair> &incumbent,
           vector<VtxPair> &current, vector<int> &g0_matched, vector<int> &g1_matched,
           vector<Bidomain> &domains, vector<int> &left, vector<int> &right, unsigned int matching_size_goal) {
//    if (arguments.timeout && double(clock() - start) / CLOCKS_PER_SEC > arguments.timeout) {
//        cout << "time out" << endl;
//        exit(0);
//    }
    if (abort_due_to_timeout) {
        return;
    }
    nodes_mcsplit++;
    //if (arguments.verbose) show(current, domains, left, right);

    if (current.size() > incumbent.size()) {//incumbent 现任的
        num = 0;
        incumbent = current;
        calls_for_optimal = nodes_mcsplit + nodes_mcsolver;
        duration = std::chrono::steady_clock::now() - START;
//        if (!arguments.quiet) cout << "Incumbent size: " << incumbent.size() << endl;
    }

    unsigned int bound = current.size() + calc_bound(domains);//计算相连和不相连同构数的最大可能加上当前已经同构的点数
    if (bound <= incumbent.size() || bound < matching_size_goal) {//剪枝
        //num--;
        cut_branches_mcsplit++;
        return;
    }
    if (arguments.big_first && incumbent.size() == matching_size_goal)
        return;

    int bd_idx = -1;
    num++;
    if (num > Maxnum) {
        //cout<<M<<endl;
        M = M % 2 + 1;
        num = 0;
    }
    switch (M) {
        case 1: {
            bd_idx = select_bidomain(domains, left, lgrade, current.size());//选出domain中可能情况最少的下标
            break;
        }
        case 2: {
            bd_idx = select_bidomain(domains, left, V, current.size());//选出domain中可能情况最少的下标
            break;
        }
        default:
            break;
    }


    if (bd_idx == -1) {  // In the MCCS case, there may be nothing we can branch on
        //  cout<<endl;
        //  cout<<endl;
        //  cout<<"big"<<endl;
        return;
    }
    Bidomain &bd = domains[bd_idx];

    int v, w;
    int tmp_idx = -1;
    switch (M) {
        case 1: {
            tmp_idx = selectV_index(left, lgrade, bd.l, bd.left_len);
            break;
        }
        case 2: {
            tmp_idx = selectV_index(left, V, bd.l, bd.left_len);
            break;
        }
        default:
            break;
    }

    v = left[bd.l + tmp_idx];
    remove_vtx_from_array(left, bd.l, bd.left_len, tmp_idx);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    vector<int> wselected(g1.n, 0);
    bd.right_len--;//
    for (int i = 0; i <= bd.right_len; i++) {
        if (i != 0) {
            num++;
            if (num > Maxnum) {
                M = M % 2 + 1;
                num = 0;
            }
        }
        switch (M) {
            case 1: {
                tmp_idx = selectW_index(right, rgrade, bd.r, bd.right_len + 1, wselected);
                break;
            }
            case 2: {
                tmp_idx = selectW_index(right, Q[v], bd.r, bd.right_len + 1, wselected);
                break;
            }
            default:
                break;
        }
        w = right[bd.r + tmp_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + tmp_idx], right[bd.r + bd.right_len]);
#ifdef DEBUG
        cout<<"v= "<<v<<" w= "<<w<<endl;
     unsigned int m;
     for(m=0;m<left.size();m++)
         cout<<left[m]<<" ";
     cout<<endl;
     for(m=0;m<right.size();m++)
         cout<<right[m]<<" ";
     cout<<endl;
#endif
        unsigned int cur_len = current.size();
        vector<Bidomain> new_domains;
        switch (M) {
            case 1: {
                new_domains = rewardfeed_RL(domains, bd_idx, current, left, right, lgrade, rgrade, g0, g1, v, w,
                                            arguments.directed || arguments.edge_labelled);
                break;
            }
            case 2: {
                new_domains = rewardfeed(domains, bd_idx, current, g0_matched, g1_matched, left, right, V, Q[v], g0, g1,
                                         v, w,
                                         arguments.directed || arguments.edge_labelled);
                break;
            }
            default:
                break;
        }

        dl++;
        solve(g0, g1, V, lgrade, rgrade, Q, incumbent, current, g0_matched, g1_matched, new_domains, left, right,
              matching_size_goal);
        while (current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g0_matched[pr.v] = 0;
            g1_matched[pr.w] = 0;
        }
    }
    bd.right_len++;
    if (bd.left_len == 0)
        remove_bidomain(domains, bd_idx);
    solve(g0, g1, V, lgrade, rgrade, Q, incumbent, current, g0_matched, g1_matched, domains, left, right,
          matching_size_goal);
}

/**
 *
 * @param d - Current domains.
 * @param left - Left graph nodes.
 * @param right - Right graph nodes.
 * @param g0 - Left graph.
 * @param g1 - Right graph.
 * @param left_selection - Nodes of left graph.
 * @param right_selection - Nodes of right graph.
 * @param multiway - Multiway is for directed and/or labelled graphs
 * @return - New bi-domains.
 */
vector<Bidomain>
filter_domains_for_many(const vector<Bidomain> &d, vector<int> &left, vector<int> &right, const Graph &g0,
                        const Graph &g1, vector<int> left_selection, vector<int> right_selection, bool multiway) {
    vector<Bidomain> new_d;
    vector<Bidomain> current_domains = d;
    new_d.reserve(d.size());
    int n = left_selection.size();
    for (int i = 0; i < n; ++i) {
        int v = left_selection[i];
        int w = right_selection[i];
        for (const Bidomain &old_bd: current_domains) {
            int l = old_bd.l;
            int r = old_bd.r;
            // After these two partitions, left_len and right_len are the lengths of the
            // arrays of vertices with edges from v or w (int the directed case, edges
            // either from or to v or w)
            int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
            int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
            int left_len_noedge = old_bd.left_len - left_len;
            int right_len_noedge = old_bd.right_len - right_len;
            if (left_len_noedge && right_len_noedge)
                new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
            if (multiway && left_len && right_len) {
                auto &adjrow_v = g0.adjmat[v];
                auto &adjrow_w = g1.adjmat[w];
                auto l_begin = std::begin(left) + l;
                auto r_begin = std::begin(right) + r;
                std::sort(l_begin, l_begin + left_len, [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; });
                std::sort(r_begin, r_begin + right_len, [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; });
                int l_top = l + left_len;
                int r_top = r + right_len;
                while (l < l_top && r < r_top) {
                    unsigned int left_label = adjrow_v[left[l]];
                    unsigned int right_label = adjrow_w[right[r]];
                    if (left_label < right_label) {
                        l++;
                    } else if (left_label > right_label) {
                        r++;
                    } else {
                        int lmin = l;
                        int rmin = r;
                        do { l++; } while (l < l_top && adjrow_v[left[l]] == left_label);
                        do { r++; } while (r < r_top && adjrow_w[right[r]] == left_label);
                        new_d.push_back({lmin, rmin, l - lmin, r - rmin, true});
                    }
                }
            } else if (left_len && right_len) {
                new_d.push_back({l, r, left_len, right_len, true});
            }
        }
        current_domains = new_d;
        new_d.clear();
    }
    return current_domains;
}

/**
 * Returns colours for an isomorphic type with 3 nodes. Update the colours vector accordingly.
 *
 * @param g - Given graph.
 * @param selection - Selected nodes (size must be 3).
 * @param colours - Colours vector to update.
 */
int colours_for_tri_nodes(const Graph &g, vector<int> &selection, vector<int> &colours) {
    int code1 = g.adjmat[selection[0]][selection[1]];
    int code2 = g.adjmat[selection[0]][selection[2]];
    int code3 = g.adjmat[selection[1]][selection[2]];
    vector<int> new_selection;
    int total = code1 + code2 + code3;
    if (total == 3) {
        colours.push_back(0);
        colours.push_back(0);
        colours.push_back(0);
        return 3;
    }
    if (code1 == 0) {
        colours.push_back(0);
        colours.push_back(0);
        colours.push_back(1);
        new_selection.push_back(selection[0]);
        new_selection.push_back(selection[1]);
        new_selection.push_back(selection[2]);
        selection = new_selection;
        return 2;
    }
    if (code2 == 0) {
        colours.push_back(0);
        colours.push_back(0);
        colours.push_back(1);
        new_selection.push_back(selection[0]);
        new_selection.push_back(selection[2]);
        new_selection.push_back(selection[1]);
        selection = new_selection;
        return 2;
    }
    if (code3 == 0) {
        colours.push_back(0);
        colours.push_back(0);
        colours.push_back(1);
        new_selection.push_back(selection[1]);
        new_selection.push_back(selection[2]);
        new_selection.push_back(selection[0]);
        selection = new_selection;
        return 2;
    }
    return 2;
}

/**
 * Returns colours for an isomorphic type with 3 nodes. Update the colours vector accordingly.
 *
 * @param g - Given graph.
 * @param selection - Selected nodes (size must be 3).
 * @param colours - Colours vector to update.
 */
bool colours_for_tri_nodes_pair(const Graph &g1, vector<int> &right_selection,
                                const vector<int> &colours) {
    vector<int> new_right_selection;
    int g0_code1 = (colours[0] == 0);
    int g0_code2 = (colours[1] == 0);
    int g0_code3 = (colours[2] == 0);

    int g1_code1 = g1.adjmat[right_selection[0]][right_selection[1]];
    int g1_code2 = g1.adjmat[right_selection[0]][right_selection[2]];
    int g1_code3 = g1.adjmat[right_selection[1]][right_selection[2]];
    int g0_total = g0_code1 + g0_code2 + g0_code3;
    int g1_total = g1_code1 + g1_code2 + g1_code3;
    if (g0_total != g1_total) {
        return false;
    }
    if (g0_total == 3) {
        return true;
    }
    // Right graph (colours are same as above)
    if (g1_code1 == 0) {
        new_right_selection.push_back(right_selection[0]);
        new_right_selection.push_back(right_selection[1]);
        new_right_selection.push_back(right_selection[2]);
        right_selection = new_right_selection;
    }
    if (g1_code2 == 0) {
        new_right_selection.push_back(right_selection[0]);
        new_right_selection.push_back(right_selection[2]);
        new_right_selection.push_back(right_selection[1]);
        right_selection = new_right_selection;
    }
    if (g1_code3 == 0) {
        new_right_selection.push_back(right_selection[1]);
        new_right_selection.push_back(right_selection[2]);
        new_right_selection.push_back(right_selection[0]);
        right_selection = new_right_selection;
    }
    return true;
}

/**
 * Find non-selected nodes and update the neigh vector.
 *
 * @param selected_indicator - Vector to indicate selected/visited nodes.
 * @param neigh - Neighbourhood.
 */
void neighbours(vector<bool> &selected_indicator, vector<int> &neigh) {
    for (size_t i = 0; i < selected_indicator.size(); ++i) {
        if (!selected_indicator[i]) {
            neigh.push_back(i);
        }
    }
}

/**
 * Convert a vector to a string with a seperator as comma.
 * @param vec - Input vector.
 * @return - Vector as a string.
 */
std::string vector_to_string(vector<int> vec) {
    if (vec.size() > 2) {
        std::sort(vec.begin(), vec.end());
    }
    std::string s = "";
    for (int i: vec) {
        s = s + std::to_string(i) + ",";
    }
    return s;
}

/**
 * Returns a hash value for a given string.
 *
 * @param s - Input string.
 * @return - Hash.
 */
size_t string_to_hash(std::string s) {
    return std::hash<std::string>{}(s);
}

/**
 * Colour the neighbourhood subgraph of selected nodes w.r.t the colours of the selected nodes.
 *
 * @param g - Input graph.
 * @param neigh - Neighbourhood.
 * @param selection - Selected nodes.
 * @param current_colours - Colours of the selected nodes.
 * @param neigh_colours - Neighbourhood colours, which given by the connectivity to selected nodes.
 * @param neigh_unique_colours - Final colour of the neighbourhood nodes.
 * @param colour_map - Already used colour map.
 * @param last_colour - Last used colour (i.e, a natural number).
 */
void neigh_colouring(const Graph &g, vector<int> &neigh, vector<int> &selection, vector<int> &current_colours,
                     vector<vector<int> > &neigh_colours, vector<int> &neigh_unique_colours,
                     std::unordered_map<size_t, int> &colour_map, int &last_colour) {
    for (size_t i = 0; i < neigh.size(); ++i) {
        vector<int> colour;
        for (size_t j = 0; j < selection.size(); ++j) {
            if (g.adjmat[neigh[i]][selection[j]]) {
                colour.push_back(current_colours[j]);
            }
        }
        neigh_colours.push_back(colour);
    }
    for (size_t i = 0; i < neigh_colours.size(); ++i) {
        if (neigh_colours[i].size() == 0) {
            neigh_unique_colours.push_back(-1);
        } else {
            size_t hash = string_to_hash(vector_to_string(neigh_colours[i]));
            if (colour_map.find(hash) == colour_map.end()) {
                colour_map[hash] = ++last_colour;
                neigh_unique_colours.push_back(last_colour);
            } else {
                neigh_unique_colours.push_back(colour_map[hash]);
            }
        }
    }
}

struct SymmetricDomains {
    vector<Bidomain> new_domains;
    bool is_isomorphic;
    vector<int> left_selection;
    vector<int> right_selection;
    vector<int> colours;
};

/**
 * Check isomorphism of the two ordered set of nodes.
 * @param g0 - Left graph.
 * @param g1 - Right graph.
 * @param left - Ordered left selection.
 * @param right - Ordered right selection.
 * @return - Isomorphic or not.
 */
bool is_isomorphic(const Graph &g0, const Graph &g1, vector<int> left, vector<int> right) {
    size_t size = left.size();
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {
            if (g0.adjmat[left[i]][left[j]] != g1.adjmat[right[i]][right[j]]) {
                return false;
            }
        }
    }
    return true;
}

bool is_canon_check_needed(const Graph &g, vector<int> &left_selection, vector<int> &right_selection,
                           vector<int> &prev_aut_colours, std::unordered_map<int, int> &colour_map) {
    for (int colour: prev_aut_colours) {
        if (colour_map.find(colour) == colour_map.end()) {
            colour_map[colour] = 1;
        } else {
            int val = ++(colour_map[colour]);
            colour_map[colour] = val;
        }
    }
    int size = left_selection.size();
    int v = left_selection[size - 1];
    bool danger = false;
    bool danger3 = false;
    int danger_colour = -1;
    int danger_colour_connections = 0;
    for (int i = 0; i < (size - 1); ++i) {
        int v_dash = left_selection[i];
        if (g.adjmat[v][v_dash] == 1) {
            int colour = prev_aut_colours[i];
            if (colour_map.at(colour) > 1) {
                if (colour_map.at(colour) > 3) {
                    if (danger3) {
                        return true;
                    }
                    danger3 = true;
                }
                if (danger_colour == -1) {
                    danger_colour = colour;
                    danger = true;
                    ++danger_colour_connections;
                    continue;
                }
                if (colour != danger_colour) {
                    return true;
                }
                ++danger_colour_connections;
            }
        }
    }
    if (danger && danger_colour_connections < colour_map.at(danger_colour)) {
        return true;
    }
    return false;
}

struct CanonLabeling {
    sparsegraph canon_graph;
    vector<int> aut_colours;
    vector<int> lab_g;
    std::unordered_map<int, int> prev_colour_map;
    unsigned int n_symmetries = 1;
    bool is_canon_needed;
};

CanonLabeling canon_labeling_left(const Graph &g0, vector<int> &left_selection, vector<int> &right_selection,
                                  vector<int> &prev_aut_colours) {
    vector<int> aut_colours;
    size_t subgraph_size = left_selection.size();
    CanonLabeling canon_labeling;
    std::unordered_map<int, int> prev_colour_map;
    if (subgraph_size == 1) {
        aut_colours.push_back(0);
        canon_labeling.is_canon_needed = false;
    } else if (subgraph_size == 2) {
        aut_colours.push_back(0);
        aut_colours.push_back(0);
        canon_labeling.is_canon_needed = false;
        canon_labeling.n_symmetries = 2;
    } else if (subgraph_size == 3) {
        canon_labeling.n_symmetries = colours_for_tri_nodes(g0, left_selection, aut_colours);
        canon_labeling.is_canon_needed = false;
    } else {
        canon_labeling.is_canon_needed = is_canon_check_needed(g0, left_selection, right_selection,
                                                               prev_aut_colours, prev_colour_map);
        DYNALLSTAT(int, ptn_g, ptn1_size);
        DYNALLOC1(int, ptn_g, ptn1_size, subgraph_size, "malloc");
        int orbits_g[subgraph_size];
        int lab_g[subgraph_size];

        static DEFAULTOPTIONS_SPARSEGRAPH(canon_options);
        statsblk stats_g;
        sparsegraph subgraph_g;
        SG_INIT(subgraph_g);

        std::vector<unsigned int> edges_g;
        int degrees_g[subgraph_size];
        size_t v_g[subgraph_size];
        size_t count_g = 0;
        for (size_t i = 0; i < subgraph_size; ++i) {
            unsigned int node1 = left_selection[i];
            int degree = 0;
            v_g[i] = count_g;
            for (size_t j = 0; j < subgraph_size; ++j) {
                unsigned int node2 = left_selection[j];
                if (g0.adjmat[node1][node2] == 1) {
                    ++degree;
                    edges_g.push_back(j);
                }
            }
            degrees_g[i] = degree;
            count_g += degree;
        }

        int nde_g = edges_g.size();
        SG_ALLOC(subgraph_g, subgraph_size, nde_g, "malloc");

        subgraph_g.d = degrees_g;
        subgraph_g.v = v_g;
        subgraph_g.nv = subgraph_size;
        subgraph_g.nde = nde_g;
        for (int i = 0; i < nde_g; ++i) {
            subgraph_g.e[i] = edges_g[i];
        }

        SG_DECL(cg_g);
        canon_options.getcanon = TRUE;
        sparsenauty(&subgraph_g, lab_g, ptn_g, orbits_g, &canon_options, &stats_g, &cg_g);

        int last_colour = 0;
        unsigned int n_symmetries = 0;
        aut_colours = vector<int>(subgraph_size, -1);
        std::unordered_map<int, int> orbit_colours_map = std::unordered_map<int, int>();
        for (size_t i = 0; i < subgraph_size; ++i) {
            int index = lab_g[i];
            int given_colour = orbits_g[index];
            if (orbit_colours_map.find(given_colour) == orbit_colours_map.end()) {
                aut_colours[index] = last_colour;
                orbit_colours_map[given_colour] = last_colour;
                ++last_colour;
            } else {
                ++n_symmetries;
                aut_colours[index] = orbit_colours_map[given_colour];
            }
        }
        canon_labeling.canon_graph = cg_g;
        vector<int> lab_g_vec(lab_g, lab_g + subgraph_size);
        canon_labeling.lab_g = lab_g_vec;
        canon_labeling.n_symmetries = std::max(uint(1), n_symmetries);
    }
    canon_labeling.aut_colours = aut_colours;
    canon_labeling.prev_colour_map = prev_colour_map;
    return canon_labeling;
}

struct Automorphism {
    vector<int> left;
    vector<int> right;
    vector<int> aut_colours;
    bool is_isomorphic;

};

/**
 * Find symmetric (automorphism) colours.
 *
 * @param g - Given graph.
 * @param selected_nodes - Nodes in the selected subgraph.
 * @return - Symmetric colours.
 */
bool is_isomorphic_with_canon_right(const Graph &g1, vector<int> &right_selection, vector<int> &left_selection,
                                    sparsegraph &cg_g, vector<int> lab_g,
                                    vector<int> &aut_colours, vector<int> &prev_aut_colours,
                                    std::unordered_map<int, int> &colour_map, bool is_canon_needed) {
    size_t subgraph_size = right_selection.size();
    if (subgraph_size == 1) {
        return true;
    } else if (subgraph_size == 2) {
        return true;
    } else if (subgraph_size == 3) {
        return colours_for_tri_nodes_pair(g1, right_selection, aut_colours);
    } else {
        if (!is_canon_needed) {
            return true;
        }
        DYNALLSTAT(int, ptn_h, ptn2_size);
        DYNALLOC1(int, ptn_h, ptn2_size, subgraph_size, "malloc");
        int orbits_h[subgraph_size];
        int lab_h[subgraph_size];

        static DEFAULTOPTIONS_SPARSEGRAPH(canon_options);
        statsblk stats_h;
        sparsegraph subgraph_h;
        SG_INIT(subgraph_h);

        std::vector<unsigned int> edges_h;
        int degrees_h[subgraph_size];
        size_t v_h[subgraph_size];
        size_t count_h = 0;
        for (size_t i = 0; i < subgraph_size; ++i) {
            unsigned int node1 = right_selection[i];
            int degree = 0;
            v_h[i] = count_h;
            for (size_t j = 0; j < subgraph_size; ++j) {
                unsigned int node2 = right_selection[j];
                if (g1.adjmat[node1][node2] == 1) {
                    ++degree;
                    edges_h.push_back(j);
                }
            }
            degrees_h[i] = degree;
            count_h += degree;
        }

        int nde_h = edges_h.size();
        SG_ALLOC(subgraph_h, subgraph_size, nde_h, "malloc");

        subgraph_h.d = degrees_h;
        subgraph_h.v = v_h;
        subgraph_h.nv = subgraph_size;
        subgraph_h.nde = nde_h;
        for (int i = 0; i < nde_h; ++i) {
            subgraph_h.e[i] = edges_h[i];
        }

        SG_DECL(cg_h);
        canon_options.getcanon = TRUE;
        sparsenauty(&subgraph_h, lab_h, ptn_h, orbits_h, &canon_options, &stats_h, &cg_h);
        bool is_isomorphic = aresame_sg(&cg_g, &cg_h);

        if (!(is_isomorphic)) {
            return false;
        }
        vector<int> new_left_selection;
        vector<int> new_right_selection;
        for (size_t i = 0; i < subgraph_size; ++i) {
            new_left_selection.push_back(left_selection[lab_g[i]]);
            new_right_selection.push_back(right_selection[lab_h[i]]);
        }
        right_selection = new_right_selection;
        left_selection = new_left_selection;
        int last_colour = 0;
        vector<int> new_aut_colour;
        std::unordered_map<int, int> orbit_colours_map = std::unordered_map<int, int>();
        for (size_t i = 0; i < subgraph_size; ++i) {
            int a = orbits_h[lab_g[i]];
            if (orbit_colours_map.find(a) == orbit_colours_map.end()) {
                new_aut_colour.push_back(last_colour);
                orbit_colours_map[a] = last_colour;
                ++last_colour;
            } else {
                new_aut_colour.push_back(orbit_colours_map[a]);
            }
        }
        aut_colours = new_aut_colour;
        return true;
    }
}

/**
 * Recreate the partitions for a given node pair (v, w).
 * @param d - Current domains.
 * @param left - Left graph nodes.
 * @param right - Right graph nodes.
 * @param g0 - Left graph.
 * @param g1 - Right graph.
 * @param v - Left graph node.
 * @param w - Right graph node.
 * @param multiway - Multiway is for directed and/or labelled graphs
 * @return - New bi-domains.
 */
vector<Bidomain>
filter_domains(const vector<Bidomain> &d, vector<int> &left, vector<int> &right, const Graph &g0, const Graph &g1,
               int v, int w, bool multiway) {
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd: d) {
        int l = old_bd.l;
        int r = old_bd.r;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
        int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            auto &adjrow_v = g0.adjmat[v];
            auto &adjrow_w = g1.adjmat[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin + left_len, [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(r_begin, r_begin + right_len, [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l < l_top && r < r_top) {
                unsigned int left_label = adjrow_v[left[l]];
                unsigned int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l < l_top && adjrow_v[left[l]] == left_label);
                    do { r++; } while (r < r_top && adjrow_w[right[r]] == left_label);
                    new_d.push_back({lmin, rmin, l - lmin, r - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }
    return new_d;
}

/**
 * Filter domains based on the symmetric colours.
 *
 * @param d - Current domains.
 * @param left - Left graph nodes.
 * @param right - Right graph nodes.
 * @param g0 - Left graph.
 * @param g1 - Right graph.
 * @param current - Current matching node pairs.
 * @param left_watched - Already watched left nodes.
 * @param right_watched - Already watched right nodes.
 * @param left_selection - Nodes of left graph.
 * @param right_selection - Nodes of right graph.
 * @param aut_colours - Symmetric colours.
 * @param multiway - Multiway is for directed and/or labelled graphs
 * @return - New bi-domains.
 */
SymmetricDomains
filter_domains_with_sym(const vector<Bidomain> &d, vector<int> &left, vector<int> &right, const Graph &g0,
                        const Graph &g1, vector<VtxPair> &current, vector<bool> &left_watched,
                        vector<bool> &right_watched, vector<int> &left_selection, vector<int> &right_selection,
                        vector<int> &aut_colours, bool multiway) {
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    size_t subgraph_size = left_selection.size();
    SymmetricDomains sd;
    if (subgraph_size == 1) {
        int v = left_selection[0];
        int w = right_selection[0];
        new_d = filter_domains(d, left, right, g0, g1, v, w, multiway);
        current.push_back(VtxPair(v, w));
        sd.new_domains = new_d;
        sd.is_isomorphic = true;
        sd.left_selection = left_selection;
        sd.right_selection = right_selection;
        sd.colours = aut_colours;
        return sd;
    }
    vector<bool> left_selected_indicator = vector<bool>(g0.n, false);
    vector<bool> right_selected_indicator = vector<bool>(g1.n, false);
    vector<int> left_neigh;
    vector<int> right_neigh;
    vector<vector<int> > left_neigh_colours;
    vector<vector<int> > right_neigh_colours;
    vector<int> left_neigh_unique_colours;
    vector<int> right_neigh_unique_colours;
    int last_colour = 0;
    std::unordered_map<size_t, int> colour_map = std::unordered_map<size_t, int>();

    for (size_t i = 0; i < subgraph_size; ++i) {
        left_selected_indicator[left_selection[i]] = true;
        right_selected_indicator[right_selection[i]] = true;
    }

    neighbours(left_selected_indicator, left_neigh);
    neighbours(right_selected_indicator, right_neigh);
    neigh_colouring(g0, left_neigh, left_selection, aut_colours, left_neigh_colours, left_neigh_unique_colours,
                    colour_map, last_colour);
    neigh_colouring(g1, right_neigh, right_selection, aut_colours, right_neigh_colours, right_neigh_unique_colours,
                    colour_map, last_colour);

    int left_n_size = left_neigh.size();
    int right_n_size = right_neigh.size();
    vector<size_t> left_neigh_indices(left_n_size);
    vector<size_t> right_neigh_indices(right_n_size);
    vector<int> new_left_neigh(left_n_size);
    vector<int> new_right_neigh(right_neigh.size());
    std::iota(left_neigh_indices.begin(), left_neigh_indices.end(), 0);
    std::iota(right_neigh_indices.begin(), right_neigh_indices.end(), 0);

    std::sort(left_neigh_indices.begin(), left_neigh_indices.end(),
              [&](size_t i, size_t j) { return left_neigh_unique_colours[i] < left_neigh_unique_colours[j]; });
    std::sort(right_neigh_indices.begin(), right_neigh_indices.end(),
              [&](size_t i, size_t j) { return right_neigh_unique_colours[i] < right_neigh_unique_colours[j]; });

    for (int i = 0; i < left_n_size; ++i) {
        new_left_neigh[i] = left_neigh[left_neigh_indices[i]];
    }
    for (int i = 0; i < right_n_size; ++i) {
        new_right_neigh[i] = right_neigh[right_neigh_indices[i]];
    }
    std::sort(left_neigh_unique_colours.begin(), left_neigh_unique_colours.end());
    std::sort(right_neigh_unique_colours.begin(), right_neigh_unique_colours.end());

    int left_index = 0;
    int right_index = 0;
    int left_partition_len = 0;
    int right_partition_len = 0;
    while (left_index < left_n_size && right_index < right_n_size) {
        int left_colour = left_neigh_unique_colours[left_index];
        int right_colour = right_neigh_unique_colours[right_index];
        if (left_colour == right_colour) {
            left_partition_len = 0;
            right_partition_len = 0;
            int left_start_index = left_index;
            int right_start_index = right_index;
            while (left_neigh_unique_colours[left_index] == left_colour && left_index < left_n_size) {
                ++left_index;
                ++left_partition_len;
            }
            while (right_neigh_unique_colours[right_index] == right_colour && right_index < right_n_size) {
                ++right_index;
                ++right_partition_len;
            }
            if (left_partition_len > 0 && right_partition_len > 0) {
                if (left_colour == -1) {
                    new_d.push_back(
                            {left_start_index, right_start_index, left_partition_len, right_partition_len, false});
                } else {
                    new_d.push_back(
                            {left_start_index, right_start_index, left_partition_len, right_partition_len, true});
                }
            }

            if (left_index >= left_n_size || right_index >= right_n_size) {
                break;
            }
        } else {
            if (left_colour < right_colour) {
                while (left_neigh_unique_colours[left_index] == left_colour && left_index < left_n_size) {
                    ++left_index;
                }
            } else {
                while (right_neigh_unique_colours[right_index] == right_colour && right_index < right_n_size) {
                    ++right_index;
                }
            }
        }
    }
    for (size_t i = 0; i < current.size(); ++i) {
        new_left_neigh.push_back(current[i].v);
        new_right_neigh.push_back(current[i].w);
    }
    left = new_left_neigh;
    right = new_right_neigh;
    current.clear();
    for (size_t i = 0; i < subgraph_size; ++i) {
        current.push_back(VtxPair(left_selection[i], right_selection[i]));
    }
    sd.new_domains = new_d;
    sd.is_isomorphic = true;
    sd.left_selection = left_selection;
    sd.right_selection = right_selection;
    sd.colours = aut_colours;
    return sd;
}

/**
 * Find factorial of n.
 * @param n - Positive integer.
 * @return - Factorial of n.
 */
int factorial(const int n) {
    int result = 1;
    for (int i = 1; i <= n; ++i)
        result *= i;
    return result;
}

/**
 * Create all possible permutations of vector.
 * @param vec - Input vector.
 * @return - All permutations of the input vector.
 */
vector<vector<int> > permutations(vector<int> vec) {
    vector<vector<int> > all_permutations;
    std::sort(vec.begin(), vec.end());
    do {
        vector<int> permute;
        for (size_t i = 0; i < vec.size(); ++i) {
            permute.push_back(vec[i]);
        }
        all_permutations.push_back(vec);
    } while (std::next_permutation(vec.begin(), vec.end()));
    return all_permutations;
}

vector<vector<vector<int> > > permutations_for_vectors(vector<vector<int>> all_vectors) {
    vector<vector<vector<int> > > all_permutations;
    for (vector<int> vec: all_vectors) {
        all_permutations.push_back(permutations(vec));
    }
    return all_permutations;
}

/**
 * Find all possible permutations of singletons and symmetric colours.
 * @param all_vectors - Vector of vectors that contains non-repeating and repeating vectors.
 * @return - All permutations.
 */
vector<vector<vector<int> > > automorphism(vector<vector<int> > all_vectors) {
    int all_count = 1;
    vector<int> sizes;
    for (vector<int> v: all_vectors) {
        int f = factorial(v.size());
        all_count = all_count * f;
        sizes.push_back(f);
    }

    int n_vectors = all_vectors.size();
    vector<vector<vector<int> > > all = permutations_for_vectors(all_vectors);
    vector<vector<vector<int> > > permutations;
    for (int i = 0; i < all_count; ++i) {
        int carry = i;
        vector<vector<int> > entry;
        for (int j = 0; j < n_vectors; ++j) {
            entry.push_back(all[j][carry % sizes[j]]);
            carry = carry / sizes[j];
        }
        permutations.push_back(entry);
    }
    return permutations;
}

struct Partitions {
    vector<int> l_nodes;
    vector<vector<int>> r_automorphisms;
    vector<int> r_singletons;
    vector<int> sorted_colours;
};

/**
 * Create the compact vector from repeated permutation (instances of repeats) + non-repeating ones (singletons).
 * @param automorphisms - Instances of repeats.
 * @param singletons - Non-repeating ones (singletons).
 * @return - Compacted vector.
 */
vector<int> compact(vector<vector<int> > automorphisms, vector<int> singletons) {
    vector<int> selected_nodes;
    for (vector<int> vec: automorphisms) {
        selected_nodes.insert(selected_nodes.end(), vec.begin(), vec.end());
    }
    for (int num: singletons) {
        selected_nodes.push_back(num);
    }
    return selected_nodes;
}

/**
 * Permute the right vector based on the colours. No need to permute left and right both.
 *
 * @param left - Left selected nodes.
 * @param right - Right selected nodes.
 * @param colours - Automorphism colours.
 * @return - Permutations of right.
 */
Partitions auto_permutations(vector<int> left, vector<int> right, vector<int> colours) {
    int n = left.size();
    vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) { return colours[i] < colours[j]; });

    vector<int> sorted_colours;
    for (int i = 0; i < n; ++i) {
        sorted_colours.push_back(colours[indices[i]]);
    }
    int index = 0;
    vector<vector<int> > l_partitions;
    vector<vector<int> > r_partitions;
    while (index < n) {
        int current_colour = sorted_colours[index];
        vector<int> l_cell;
        vector<int> r_cell;
        if (current_colour == sorted_colours[index]) {
            while (current_colour == sorted_colours[index] && index < n) {
                l_cell.push_back(left[indices[index]]);
                r_cell.push_back(right[indices[index]]);
                ++index;
            }
            l_partitions.push_back(l_cell);
            r_partitions.push_back(r_cell);
        } else {
            ++index;
        }
    }

    vector<vector<int> > l_automorphisms;
    vector<vector<int> > r_automorphisms;
    vector<int> l_singletons;
    vector<int> r_singletons;
    vector<int> aut_colours;
    vector<int> singleton_colours;
    for (size_t k = 0; k < r_partitions.size(); ++k) {
        if (r_partitions[k].size() == 1) {
            l_singletons.push_back(l_partitions[k][0]);
            r_singletons.push_back(r_partitions[k][0]);
            singleton_colours.push_back(sorted_colours[k]);
        } else {
            l_automorphisms.push_back(l_partitions[k]);
            r_automorphisms.push_back(r_partitions[k]);
            for (size_t i = 0; i < r_partitions[k].size(); ++i) {
                aut_colours.push_back(sorted_colours[k]);
            }
        }
    }
    aut_colours.insert(aut_colours.end(), singleton_colours.begin(), singleton_colours.end());
    vector<int> left_nodes = compact(l_automorphisms, l_singletons);
    return {left_nodes, r_automorphisms, r_singletons, aut_colours};
}

/**
 * Find a vector has repeated elements or not.
 * @param vec - Input vector.
 * @return - True if it does not has repeats. False otherwise.
 */
bool is_unique(vector<int> vec) {
    std::sort(vec.begin(), vec.end());
    return ((vec.size() > 1) && std::unique(vec.begin(), vec.end()) == vec.end());
}

//int debug_index = 0;

int calculate_mcsolver_bound(const Graph &g, const vector<bool> &indicator) {
    int bound = g.n;
    for (int i = 0; i < g.n; ++i) {
        if (indicator[i]) {
            bound = bound - 1;
        }
    }
    return bound;
}

/**
 * Node selection heuristic. Find the highest degree node (i.e., node with minimum number is the highest degree).
 *
 * @param arr - Partition (i.e., same coloured node set).
 * @param start_idx - Starting index of the partition.
 * @param len - Length of the partition.
 * @return - Selected node.
 */
int find_min_value(const vector<int> &arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++) {
        if (arr[start_idx + i] < min_v) {
            min_v = arr[start_idx + i];
        }
    }
    return min_v;
}

int find_min_value1(const vector<int> &arr, int start_idx, int len, vector<bool> watched) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++) {
        int idx = start_idx + i;
        if (!watched[arr[idx]] && arr[idx] < min_v) {
            min_v = arr[idx];
        }
    }
    return min_v;
}

/**
 * Remove vertex from a partition after selection.
 *
 * @param left - Partition.
 * @param bd - Bidomain.
 * @param v - Vertex.
 */
void remove_vtx_from_left_domain(vector<int> &left, Bidomain &bd, int v) {
    int i = 0;
    while (left[bd.l + i] != v) i++;
    std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
    bd.left_len--;
}

int debug_index = 0;
/**
 * Recursive solver of MCSolver with symmetry.
 *
 * @param g0 - Left graph.
 * @param g1 - Right graph.
 * @param current - Current solution.
 * @param domains - Current domains.
 * @param left - Left graph nodes.
 * @param right - Right graph nodes.
 * @param matching_size_goal - Matching goal.
 * @param left_watched - Already watched left graph nodes.
 * @param right_watched - Already watched right graph nodes.
 * @param symmetric_bound - Bound to terminate symmetric colouring and move to MCSPlit (i.e., \tau).
 */
void mcsp_sym_without_bfs_bound(const Graph &g0, const Graph &g1, vector<VtxPair> current, vector<Bidomain> &domains,
                                vector<int> left, vector<int> right, unsigned int matching_size_goal,
                                vector<bool> left_watched, vector<bool> right_watched, vector<int> &prev_aut_colours,
                                int symmetric_bound, vector<gtype> &lgrade, vector<gtype> &rgrade,
                                vector<vector<gtype>> &Q, vector<gtype> &V,
                                vector<int> &g0_matched, vector<int> &g1_matched) {
    ++debug_index;
    int current_index = debug_index;
    if (current.size() > incumbent.size()) {
        duration = std::chrono::steady_clock::now() - START;
        incumbent = current;
        calls_for_optimal = nodes_mcsplit + nodes_mcsolver;
    }

    if (min_graph_size == incumbent.size()) {
        return;
    }
    if (abort_due_to_timeout)
        return;
    nodes_mcsolver++;

    int current_size = current.size();
    unsigned int left_bound = calculate_mcsolver_bound(g0, left_watched);
    unsigned int right_bound = calculate_mcsolver_bound(g1, right_watched);
    unsigned int bound = std::min(left_bound, right_bound);
    if (bound <= incumbent.size() || (bound - current_size) == 0 || bound < matching_size_goal) {
        cut_branches_mcsolver++;
        return;
    }

    if (arguments.big_first && incumbent.size() == matching_size_goal)
        return;

    int bd_idx = select_bidomain(domains, left, lgrade, current.size());
    if (bd_idx == -1) {
        return;
    }
    vector<int> original_prev_aut_colours = prev_aut_colours;
    Bidomain &bd = domains[bd_idx];
    int v = find_min_value1(left, bd.l, bd.left_len, left_watched);
    if (v < g0.n) {
        remove_vtx_from_left_domain(left, domains[bd_idx], v);
        int w = -1;
        bd.right_len--;
        vector<int> left_selection;
        vector<int> right_selection_init;
        for (VtxPair pair: current) {
            left_selection.push_back(pair.v);
            right_selection_init.push_back(pair.w);
        }
        left_selection.push_back(v);
        vector<bool> right_watched_copy(g1.n, false);
        int new_size = left_selection.size();
        int incumbent_size = incumbent.size();
        CanonLabeling left_canon_labeling = canon_labeling_left(g0, left_selection, right_selection_init,
                                                                prev_aut_colours);
        bool is_discrete = is_unique(left_canon_labeling.aut_colours);
        vector<int> wselected(g1.n, 0);
        for (int i = 0; i <= bd.right_len; i++) {
//            int idx = index_of_next_smallest(right, bd.r, bd.right_len + 1, w);
            int tmp_idx = -1;
            switch (M) {
                case 1: {
                    tmp_idx = selectW_index(right, rgrade, bd.r, bd.right_len + 1, wselected);
                    break;
                }
                case 2: {
                    tmp_idx = selectW_index(right, Q[v], bd.r, bd.right_len + 1, wselected);
                    break;
                }
                default:
                    break;
            }
            w = right[bd.r + tmp_idx];
            wselected[w] = 1;
            w = right[bd.r + tmp_idx];
            if (right_watched[w] == true) {
                continue;
            }
            if (v == 0) {
                symmetric_bound = symmetric_bound - 1;
            }

            bound = calculate_mcsolver_bound(g1, right_watched);
            if (bound <= incumbent.size()) {
                cut_branches_mcsolver++;
                return;
            }

            // swap w to the end of its colour class
            if (current.size() == 0) {
                right[bd.r + tmp_idx] = right[bd.r + bd.right_len];
                right[bd.r + bd.right_len] = w;
            }

            vector<int> right_selection = right_selection_init;
            vector<VtxPair> current_copy = current;
            right_selection.push_back(w);
            vector<int> left_selection_copy = left_selection;
            bool is_isomorphic_selection = is_isomorphic_with_canon_right(g1, right_selection,
                                                                          left_selection_copy,
                                                                          left_canon_labeling.canon_graph,
                                                                          left_canon_labeling.lab_g,
                                                                          left_canon_labeling.aut_colours,
                                                                          prev_aut_colours,
                                                                          left_canon_labeling.prev_colour_map,
                                                                          left_canon_labeling.is_canon_needed);
            if (is_isomorphic_selection) {
                vector<int> left_copy = left;
                vector<int> right_copy = right;
//            cout << "MCSolver " << v << ", " << W_indicator[w] << endl;
                SymmetricDomains sd = filter_domains_with_sym(domains, left_copy, right_copy, g0, g1, current_copy,
                                                              left_watched, right_watched_copy, left_selection_copy,
                                                              right_selection,
                                                              left_canon_labeling.aut_colours,
                                                              arguments.directed || arguments.edge_labelled);
                int switch_bound = (new_size + calc_bound(sd.new_domains)) / left_canon_labeling.n_symmetries;
                if (true) {
                    mcsp_sym_without_bfs_bound(g0, g1, current_copy, sd.new_domains, left_copy, right_copy,
                                               matching_size_goal, left_watched, right_watched_copy,
                                               left_canon_labeling.aut_colours, symmetric_bound - 1, lgrade, rgrade,
                                               Q, V, g0_matched, g1_matched);
                    if (abort_due_to_timeout)
                        return;
                } else {
                    Partitions p = auto_permutations(sd.left_selection, sd.right_selection,
                                                     left_canon_labeling.aut_colours);
                    vector<vector<vector<int> > > r_permutations = automorphism(p.r_automorphisms);
                    int select_len = p.l_nodes.size();
                    vector<bool> l_select_indicator(g0.n, false);
                    vector<bool> r_select_indicator(g1.n, false);
                    for (int j = 0; j < select_len; ++j) {
                        l_select_indicator[sd.left_selection[j]] = true;
                        r_select_indicator[sd.right_selection[j]] = true;
                    }
                    left_copy.clear();
                    right_copy.clear();
                    for (int j = 0; j < g0.n; ++j) {
                        if (!l_select_indicator[j] && !left_watched[j]) {
                            left_copy.push_back(j);
                        }
                    }
                    for (int j = 0; j < g1.n; ++j) {
                        if (!r_select_indicator[j]) {
                            right_copy.push_back(j);
                        }
                    }
                    int left_watched_len = 0;
                    for (int j = 0; j < g0.n; ++j) {
                        if (left_watched[j] && !l_select_indicator[j]) {
                            left_copy.push_back(j);
                            ++left_watched_len;
                        }
                    }
                    for (int j = 0; j < select_len; ++j) {
                        left_copy.push_back(sd.left_selection[j]);
                        right_copy.push_back(sd.right_selection[j]);
                    }
                    auto bd_init = vector<Bidomain>{};
                    bd_init.push_back({0, 0, g0.n - select_len - left_watched_len, g1.n - select_len, false});
                    for (vector<vector<int> > r_perm: r_permutations) {
                        vector<int> right_permutation = compact(r_perm, p.r_singletons);
                        if (is_isomorphic(g0, g1, p.l_nodes, right_permutation)) {
                            vector<int> left_copy1 = left_copy;
                            vector<int> right_copy1 = right_copy;
                            auto new_d = filter_domains_for_many(bd_init, left_copy1, right_copy1, g0, g1, p.l_nodes,
                                                                 right_permutation, false);
                            vector<VtxPair> new_current;
                            for (size_t k = 0; k < p.l_nodes.size(); ++k) {
                                new_current.push_back(VtxPair{p.l_nodes[k], right_permutation[k]});
                            }
                            vector<VtxPair> new_current1 = new_current;
                            vector<Bidomain> new_d1 = new_d;
//                            solve(g0, g1, new_current1, new_d1, left_copy1, right_copy1, matching_size_goal, bound);
                            solve(g0, g1, V, lgrade, rgrade, Q, incumbent, current,
                                  g0_matched, g1_matched, new_d1, left, right,
                                  matching_size_goal);
                        }
                        if (abort_due_to_timeout)
                            return;
                        // Just return is incorrect since it might miss some branches.
                    }
                }
                right_watched_copy[w] = true;
            }
        }
        bd.right_len++;
        // We need this removal only in the implementation to avoid incorrect partition selection, which leads to incorrect
        // memory access. But we don't do it in conceptual level since we assume SelectPartition heuristic will handle that.
        if (bd.left_len == 0)
            remove_bidomain(domains, bd_idx);
        left_watched[v] = true;
    } else {
        remove_bidomain(domains, bd_idx);
    }
    mcsp_sym_without_bfs_bound(g0, g1, current, domains, left, right, matching_size_goal, left_watched, right_watched,
                               original_prev_aut_colours, symmetric_bound - 1, lgrade, rgrade,
                               Q, V, g0_matched, g1_matched);
}


vector<VtxPair> mcs(const Graph &g0, const Graph &g1) {
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions

    vector<int> g0_matched(g0.n, 0);
    vector<int> g1_matched(g1.n, 0);
    //*******
    Maxnum = 2 * std::min(g0.n, g1.n);
    //******
    vector<gtype> lgrade(g0.n, 0);
    vector<gtype> V(g0.n, 0);
    vector<gtype> rgrade(g1.n, 0);
    vector<vector<gtype>> Q(g0.n, vector<gtype>(g1.n, 0));


    auto domains = vector<Bidomain>{};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label: g0.label) left_labels.insert(label);
    for (unsigned int label: g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label: labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (int i = 0; i < g0.n; i++)
            if (g0.label[i] == label)
                left.push_back(i);
        for (int i = 0; i < g1.n; i++)
            if (g1.label[i] == label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    vector<VtxPair> incumbent;

    if (arguments.big_first) {
        for (int k = 0; k < g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<VtxPair> current;
            solve(g0, g1, V, lgrade, rgrade, Q, incumbent, current, g0_matched, g1_matched, domains_copy, left_copy,
                  right_copy, goal);
            if (incumbent.size() == goal || abort_due_to_timeout) break;
            if (!arguments.quiet) cout << "Upper bound: " << goal - 1 << std::endl;
        }

    } else {
        vector<VtxPair> current;
        int symmetry_bound_from_graphs = std::min(g0.n, g1.n);
        int symmetric_bound = std::min((symmetry_bound_from_graphs / 2), 6);
        vector<bool> left_watched(g0.n, false);
        vector<bool> right_watched(g1.n, false);
        vector<int> prev_aut_colours;
        mcsp_sym_without_bfs_bound(g0, g1, current, domains, left, right, 1, left_watched, right_watched,
                                   prev_aut_colours, symmetric_bound, lgrade, rgrade,
                                   Q, V, g0_matched, g1_matched);
    }

    return incumbent;
}

vector<int> calculate_degrees(const Graph &g) {
    vector<int> degree(g.n, 0);
    for (int v = 0; v < g.n; v++) {
        for (int w = 0; w < g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjmat[v][w] & mask) degree[v]++;
            if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
}

int sum(const vector<int> &vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

void pack_leaves(Graph &g) {
    vector<int> deg(g.n, 0);

    for (int i = 0; i < g.n; i++)
        for (int j = 0; j < g.n; j++)
            if (i != j && g.adjmat[i][j])
                deg[i]++;

    for (int u = 0; u < g.n; u++) {
        for (int v = 0; v < g.n; v++)
            if (g.adjmat[u][v] && u != v && deg[v] == 1) {
                std::pair<unsigned int, unsigned int> labels(g.adjmat[u][v], g.label[v]);
                int pos = -1;
                for (int k = 0;; k++) {
                    if (k == int(g.leaves[u].size())) {
                        g.leaves[u].push_back(std::make_pair(labels, vector<int>()));
                    }
                    if (g.leaves[u][k].first == labels) {
                        pos = k;
                        break;
                    }
                }
//            assert(pos != -1);
                g.leaves[u][pos].second.push_back(v);
            }
        sort(g.leaves[u].begin(), g.leaves[u].end());
    }
}

int main(int argc, char **argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    struct Graph g0 = readGraph(arguments.filename1, format, arguments.directed,
                                arguments.edge_labelled, arguments.vertex_labelled);
    struct Graph g1 = readGraph(arguments.filename2, format, arguments.directed,
                                arguments.edge_labelled, arguments.vertex_labelled);

    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;
    START = std::chrono::steady_clock::now();
    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
            auto abort_time = std::chrono::steady_clock::now() + std::chrono::seconds(arguments.timeout);
            {
                /* Sleep until either we've reached the time limit,
                 * or we've finished all the work. */
                std::unique_lock<std::mutex> guard(timeout_mutex);
                while (!abort_due_to_timeout.load()) {
                    if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                        /* We've woken up, and it's due to a timeout. */
                        aborted = true;
                        break;
                    }
                }
            }
            abort_due_to_timeout.store(true);
        });
    }
    //  auto start = std::chrono::steady_clock::now();
    start = clock();

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);

    // As implemented here, g1_dense and g0_dense are false for all instances
    // in the Experimental Evaluation section of the paper.  Thus,
    // we always sort the vertices in descending order of degree (or total degree,
    // in the case of directed graphs.  Improvements could be made here: it would
    // be nice if the program explored exactly the same search tree if both
    // input graphs were complemented.
//  #ifdef VSCORE
//    vector<int> lgrade(g0.n,0);
//    vector<int> rgrade(g1.n,0);
//  #endif
    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) > g1.n * (g1.n - 1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a] < g0_deg[b]) : (g0_deg[a] > g0_deg[b]);
    });

    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) > g0.n * (g0.n - 1);
    std::stable_sort(std::begin(vv1), std::end(vv1),
                     [&](int a, int b) {//????????????????????????????????????????????????????
                         return g0_dense ? (g1_deg[a] < g1_deg[b]) : (g1_deg[a] > g1_deg[b]);
                     });

#if 0
    int idx;
    for(idx=0;idx<vv0.size();idx++)
       cout<<vv0[idx]<<"  ";
    cout<<endl;
    for(idx=0;idx<g0_deg.size();idx++)
       cout<<g0_deg[idx]<<"  ";
    cout<<endl;
    for(idx=0;idx<vv1.size();idx++)
       cout<<vv1[idx]<<"  ";
    cout<<endl;
    for(idx=0;idx<g1_deg.size();idx++)
       cout<<"("<<idx<<","<<g1_deg[idx]<<")  ";
    cout<<endl;

#endif
    struct Graph g0_sorted = induced_subgraph(g0, vv0);
    struct Graph g1_sorted = induced_subgraph(g1, vv1);
#if 0
    int idx;
    for(idx=0;idx<g0.n;idx++)
        cout<<g0.label[idx]<<" ";
    cout<<endl;
    for(idx=0;idx<g0_sorted.n;idx++)
        cout<<g0_sorted.label[idx]<<" ";
    cout<<endl;
#endif

    pack_leaves(g0_sorted);
    pack_leaves(g1_sorted);
    vector<VtxPair> solution = mcs(g0_sorted, g1_sorted);

    // Convert to indices from original, unsorted graphs
    for (auto &vtx_pair: solution) {
        vtx_pair.v = vv0[vtx_pair.v];
        vtx_pair.w = vv1[vtx_pair.w];
    }

    // auto stop = std::chrono::steady_clock::now();
    // auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
//    clock_t time_elapsed = clock() - start;
//    clock_t time_find = bestfind - start;
    /* Clean up the timeout thread */
    std::chrono::duration<double> time_elapsed = std::chrono::steady_clock::now() - START;

    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }
//    if (!check_sol(g0, g1, solution))
//        fail("*** Error: Invalid solution\n");
//
//    cout << "Solution size " << solution.size() << std::endl;
//    for (int i=0; i<g0.n; i++)
//        for (unsigned int j=0; j<solution.size(); j++)
//            if (solution[j].v == i)
//                cout << "(" << solution[j].v << " -> " << solution[j].w << ") ";
//    cout << std::endl;
//
//    cout<<"Nodes:                      " << nodes << endl;
//    cout<<"Cut branches:               "<<cutbranches<<endl;
//    cout<<"Conflicts:                    " <<conflicts<< endl;
//    printf("CPU time (ms):              %15ld\n", time_elapsed * 1000 / CLOCKS_PER_SEC);
//    printf("FindBest time (ms):              %15ld\n", time_find * 1000 / CLOCKS_PER_SEC);
//#ifdef Best
//    cout<<"Best nodes:                 "<<bestnodes<<endl;
//    cout<<"Best count:                 "<<bestcount<<endl;
//#endif
//    if (aborted)
//        cout << "TIMEOUT" << endl;

    cout << incumbent.size() << ", " << check_sol(g0_sorted, g1_sorted, incumbent) << ", "
         << duration.count() << ", " << time_elapsed.count() << ", " << nodes_mcsolver << ", " << nodes_mcsplit
         << ", " << calls_for_optimal << ", " << cut_branches_mcsolver << ", " << cut_branches_mcsplit << ", " <<
         aborted << endl;

}

