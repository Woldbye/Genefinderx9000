// #TLE solution for addresssearch
// Brute force search through genes
#include <bits/stdc++.h>
#include <functional>
#include <iterator>
#include <experimental/iterator>
using namespace std;

#define ALL(o) (o).begin(), (o).end() 
#define FOR(iter, star, finish) for(iter = star; iter < finish; ++iter)
#define ODD(i) (bool) (i & 0x1)
#define EVEN(i) !ODD(i)
#define vec vector
#define PROGRAM_SUCCESS 0
#define os_joiner(del) experimental::make_ostream_joiner(cout, del)

typedef pair<int, int> pi;
typedef unsigned long long int ui64;

class Dictionary
{
private:
  vec<string> src;

public:
  void init(vec<string> &genes) { src = genes; };

  set<int> find_all(const string &needle) {
    set<int> ret; ui64 j;
    for (int i = 0; i < (int)src.size(); i++){
      if ((j = src[i].find(needle)) != string::npos)
        ret.insert(i);
    }
    return ret;
  };
};

Dictionary dict;
set<int> ans;
string s, needle;
int n, N, // count of dictionary words 
    q, Q; // count of queries

void print_ans() {
  cout << ans.size();
  if (ans.size() > 0) {
    cout << " ";
    copy(ALL(ans), os_joiner(" "));
  }
  cout << "\n";
}

int main() {
  cin >> N;
  vec<string> genes(N);
  
  FOR(n,0,N) cin >> genes[n];
  dict.init(genes);
  
  cin >> Q;
  FOR(q,0,Q) {
    cin >> needle;
    ans = dict.find_all(needle);
    print_ans();
  }

  return PROGRAM_SUCCESS;
}
