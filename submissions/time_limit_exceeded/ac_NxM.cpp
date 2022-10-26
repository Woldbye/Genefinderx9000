/**
 * TLE: O(N×M)
 * Aho-Corasick algorithm without exit link optimization
 */
#include <bits/stdc++.h>
#include <functional>
#include <iterator>
#include <experimental/iterator>
#include <chrono>
using namespace std;
using namespace std::chrono;

#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,fma")
#pragma GCC optimize("unroll-loops")
#define fast_io() ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);

#define beg() begin()
#define ALL(o) (o).beg(), (o).end() 
#define FOR(iter, star, finish) for(iter = star; iter < finish; ++iter)
#define vec vector
#define PROGRAM_SUCCESS 0
#define os_joiner(del) experimental::make_ostream_joiner(cout, del)
#define FAIL -1
#define PB push_back

typedef unsigned long long int ui64;
typedef pair<int, int> pi;
typedef vec<vec<int>> vvi;

template <int MAX_N>
class MultiString
{
private:
  string s;
  int snum[MAX_N];
  int N;     // count of substrings in s

public:
  MultiString() : s(""), N(0) {
    memset(snum, -1, sizeof(snum));
  }

  inline string& str()   { return s; }
  inline char &operator[](size_t i) { return s[i]; }
  inline int length() { return s.length(); }   // Count of characters
  inline int size()   { return N; }  // Count of substrings
  
  void clear() {
    s = "";
    memset(snum, -1, sizeof(snum)*length());
    N = 0;
  }

  // Append word to end of this string
  void append(string &word)
  {
    if (word.empty()) return;
    for (char c : word){
      snum[length()] = N;
      s.PB(c);
    }
    ++N;
  }

  /**
   * @brief returns the number of the string enclosing the index 'i' 
   */
  inline int strnum(int i) { return snum[i];}

  friend ostream &operator<<(ostream &os, MultiString &mstr)
  {
    for (int k = 0, i = 0; k < mstr.size(); ++k) {
      os << "\n[" << k << "]: \"";
      while (mstr.strnum(i) == k) os << mstr[i++];
      os << "\"";
    }
    return os;
  }
};

/**
 * @brief Aho-Corasick Automaton for answering multiple containment queries in parralel
 * @warning will cause stack-overflows in-case of large sample of queries or
 *          large arrays alphabets.
 * @tparam MAX_N Maximum length of text to search
 * @tparam MAX_M Maximum length of query
 * @tparam A_LEN Length of alphabet
 * @tparam (*cindx)(char) A mapping function to map any given char 'c' to an index from [0..(A_LEN-1)]
 */
template<int MAX_N, int MAX_M, int A_LEN, int (*cindx)(char)>
class AC_Matcher
{
private:
  vec<string> queries;
  int M,                // Number of nodes in trie
      fail[MAX_M],      // Suffix links : fail[m] : Link to the Longest-common-suffix of node m in the trie
      go[MAX_M][A_LEN], // Go-to Trie   : go[m][cindx(c)] - unencountered chars in alphabet are set to point to root
      out[MAX_M];       // out[m]       : contains index of all queries that ends at node 'm'
  
  void clear(int m = FAIL) {
    if (m == FAIL) m = M;
    memset(fail, FAIL, sizeof(fail[0])*m);
    memset(go, FAIL, sizeof(go[0][0])*A_LEN*m); 
    memset(out, FAIL, sizeof(out[0])*m);
    // memset(_exit, FAIL, sizeof(_exit[0])*m);
    M = 1;
  }

  /** Utility function for finding next node to visit from current pos */
  int findNode(int pos, const char c) {
    int ci = cindx(c);
    while(go[pos][ci] == FAIL) pos = fail[pos];
    return go[pos][ci];
  }

public:
  /** @brief returns count of nodes in Aho-Corasick automaton */
  int size() { return M; }
  constexpr int root() { return 0; }

  AC_Matcher() { clear(MAX_M); }

  void build(vec<string>& keywords) {
    clear(M);
    queries = keywords;
    queue<int> que; // queue for bfs of trie to init suflinks
    int K = queries.size();
    // Build trie, use 1 based k indexing to use 0 as false out
    for (int m = root(), k = 1; k <= K; ++k, m=root()) {
      for (const char &c : queries[k-1]){
        int ci = cindx(c);
        if (go[m][ci] == FAIL) {
          go[m][ci] = M++; // Add a new node
        }
        m = go[m][ci];
      }
      out[m] = k;
    }
    // Loop all edges of root node
    for (int ci = 0; ci < A_LEN; ++ci) {
      int& kid = go[root()][ci];
      // Add edge to root node for missing chars in queries
      if (kid == FAIL) kid = root();  
      else {
        fail[kid] = root(); // Link depth 1 nodes to root
        que.push(kid);    
      }
    } 

    // Initalizes suffix links in O(M) by BFS 
    //! note lemma: Suffix links point from longer strings to shorter strings
    for (int m; !que.empty(); que.pop())
    {
      m = que.front();
      // Find suffix links for all outgoing edges of nd 'm'
      for (int ci = 0, kid, f; ci < A_LEN; ++ci) {
        if ((kid = go[m][ci]) == FAIL) continue;
        // Follow failure link of m
        for (f = fail[m];
             go[f][ci] == FAIL;
             f = fail[f]);

        fail[kid] = go[f][ci];
        que.push(kid);
      }
    }
  };

  /**
   * @brief Returns a vector of all matches to queries in s in O(N+M+z) where 
   *        z = count of matches.
   *        M = Total length of queries.
   *        N = Length of source string.
   *        Note: Algorithm is output-sensitive
   * @return vec<vec<int>> of same size as 'queries', 
   *         The k'th index contain a list of all matches to query 'k' 
   */
  vvi find(MultiString<MAX_N>&s)
  {
    const int N = s.length();
    vvi ans(queries.size(), vec<int>());
    for (int n = 0, m = root(); n < N; ++n) {
      if (n && s.strnum(n) != s.strnum(n-1)) m = root();
      m = findNode(m, s[n]);
      // Follow suffix links : O(M)
      for (int f = m;
           f != root();
           f = fail[f]) //! should use exit links to achieve O(Z) 
      {
        if (out[f] > root()) {
          ans[out[f] - 1].PB(n - queries[out[f] - 1].size() + 1);
        }
      } 
    }
    return ans;
  }
};

constexpr int MAX_N = 1000001; // Max length of dictionary
constexpr int MAX_M = 1000001; // MAX_M = |query1| + ... + |queryQ|
constexpr int ALPHA_SIZE = 4; // 'A', 'C', 'G', 'T', 

// For now lets try with T = 1
int ALPHA_INDX_MAPPER(const char c) {  // Alphabet mapper
  switch (c) {
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 0; // A
  }
}

string s; 
MultiString<MAX_N> ms; // text to search
AC_Matcher<MAX_N, MAX_M, ALPHA_SIZE, ALPHA_INDX_MAPPER> matcher;

int main()
{
  fast_io();
  int I, K; // I=MS.size, K=Count of queries
  cin >> I;
  for (int i = 0; i < I; ++i) {
    cin >> s;
    ms.append(s);
  }
  
  cin >> K;
  vec<string> queries(K, "");
  for (int k = 0; k < K; ++k) cin >> queries[k];
  // Answer queries
  matcher.build(queries);
  vec<vec<int>> matches = matcher.find(ms);
  auto joiner = os_joiner(" ");
  for (int k = 0; k < K; ++k) {
    set<int> ans;
    for (const int match : matches[k]) 
      ans.insert(ms.strnum(match));
    
    cout << ans.size();
    if (!ans.empty()) {
      cout << " ";
      copy(ALL(ans), joiner);
    }
    cout << "\n";
  }

  return PROGRAM_SUCCESS;
}