// TLE: Suffix Array utilized as a Suffix Trie 
// O(N×Log^2(N)) construction 
// O(M) to answer all queries 
#include <bits/stdc++.h>
#include <functional>
#include <iterator>
#include <experimental/iterator>
using namespace std;

#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,fma")
#pragma GCC optimize("unroll-loops")
#define fast_io() ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);

#define ALL(o) (o).beg(), (o).end() 
#define FOR(iter, star, finish) for(iter = star; iter < finish; ++iter)
#define vec vector
#define PROGRAM_SUCCESS 0
#define PB push_back
#define beg begin
#define os_joiner(del) experimental::make_ostream_joiner(cout, del)

typedef unsigned long long int ui64;
typedef long long int i64;
typedef pair<int, int> pi;

static const char SEPERATOR = '~';
static const int ALPHA_SIZE = '4'; // 'A', 'G', 'C', 'T', '~'
constexpr int MAX_N = 1000001; // Max length of dictionary

pi NONE_PI(-1, -1);


/**
 * @brief A string composed of several strings.
 *        Useful for extending algorithms that are made to operate on
 *        a single string to multiple strings (i.e. suffix arrays etc.).
 */
class MultiString {
private:
  string s;
  vec<int> snum; // snum[n] : denotes the string number of the string enclosing char n numbered from 0.._size-1
  int _size;     // count of substrings in s

  /**
   * @brief Find bound of the multistring starting from pos.
   * @param i index of starting point
   * @param forward A boolean indicating whether to search forwards or backwards from pos
   * @return if there was no seperator found from i to 0
   *            -1
   *         if there was no seperator found from i to length()-1
   *            length() 
   *         otherwise
   *            index of SEPERATOR 
   */
  int find_bound(int i, bool forward) {
    if (i >= (int) s.length()) return -1;
    int inc = forward ? 1 : -1;
    for (char c = s[i];
        c != SEPERATOR && i > 0 && i < length();
        i+=inc, c=s[i]);
    return i;
  }

public:
  MultiString() : _size(0) {}

  inline const char *c_str() { return s.c_str(); }
  inline char &operator[](size_t i) { return s[i]; }
  inline int length() { return s.length(); }   // Count of characters
  inline int size() { return _size; }  // Count of substrings

  // Append word to end of this string
  void append(string &word)
  {
    if (word.empty()) return;
    s += SEPERATOR;
    snum.push_back(_size);
    for (char c : word)
    {
      s.push_back(c);
      snum.push_back(_size);
    }
    // s += SEPERATOR;
    ++_size;
  }

  string get(size_t i) {
    string out = "",
           pre = prefix(i),
           suf = suffix(i);
    if (!pre.empty() && pre[0] != SEPERATOR) out+=pre;
    if (!suf.empty() && suf[0] != SEPERATOR) out+=suf;
    return out;
  }

  /**
   * @brief returns the number of the string enclosing the index 'i' 
   */
  inline int strnum(int i) { return snum[i];}

  string prefix(size_t i) {
    int lo = max(0,find_bound(i, false)+1);
    string key = s.substr(lo, i - lo+1);
    if (key.empty()) return string(1, SEPERATOR);
    return key;
  }

  string suffix(size_t i) {
    int hi = find_bound(++i, true);
    string key = s.substr(i, hi-i);
    if (key.empty()) return string(1,SEPERATOR);
    else return key;
  }

  friend ostream& operator<<(ostream& os, MultiString& mstr) {
    os << "\n[" << 0 << "]: \"";
    for (int i =1; i < mstr.length(); i++) {
      char c = mstr[i];
      if (c == SEPERATOR)
        os << "\"\n[" << mstr.strnum(i) << "]: \"";
      else
        os << c;
    }
    return os;
  }
};

class LongestCommonPrefixArray
{
private:
  vec<int> src; // src[n] : length of the longest common prefix between sa[n], sa[n-1]

  // init longest-common-prefix array in O(N) : Kasai's Algorithm
  void init_kasai(vec<int> &sa, string &s)
  {
    const int N = (int)s.size();
    int n, len = 0;

    // Stores the inverse index of the suffix array
    // rank[n] = the position of suffix[n] in sa
    vec<int> rank(N, 0);

    FOR(n, 0, N) rank[sa[n]] = n;

    FOR(n, 0, N)
    {
      if (rank[n] > 0)
      {                                 // suffix n has 'prv' in the suffix array.
        for (int prv = sa[rank[n] - 1]; // prv = preceding suffix
             s[n + len] == s[prv + len];
             len++)
          ;
        src[rank[n]] = len;
        if (len > 0)
          len--; // if len > 0 -1 else 0
      }
    }
  }

  // O(N)
  void init_updown()
  {
    stack<int> st;
    int n, prv_i = -1, N = src.size();
    st.push(0);
    FOR(n, 0, N+1)
    {
      while (src[n] < src[st.top()])
      {
        prv_i = st.top();
        st.pop();
        if (!st.empty() && src[prv_i] != src[st.top()] && src[n] <= src[st.top()])
          down[st.top()] = prv_i;
      }

      if (src[n] >= src[st.top()])
      {
        if (prv_i != -1)
        {
          up[n] = prv_i;
          prv_i = -1;
        }
        st.push(n);
      }
    }
  }

  // O(N)
  void init_next()
  {
    stack<int> st;
    int n, prv_i, N = src.size();
    st.push(0);
    for (n = 1; n < N; st.push(n), n++)
    {
      while (!st.empty() && src[n] < src[st.top()])
        st.pop();

      if (!st.empty() && src[n] == src[st.top()]){
        prv_i = st.top();
        st.pop();
        next[prv_i] = n;
      }
    }
  }

public:
  vec<int> up,   // up[n]   : max q of all q ϵ {0..n-1}   | where S[sa[q]..].contains(S[sa[n]..])
           down, // down[n] : min q of all q ϵ {n+1..N-1} | where S[sa[q]..].contains(S[sa[n]..]))
           next; // next[n] : min q of all q ϵ {n+1..N-1} | nearest entry,  where lowest-common-prefix(S[sa[n]..], S[next[]]) == 0

  LongestCommonPrefixArray(){};

  ~LongestCommonPrefixArray()
  {
    up.clear();
    down.clear();
    next.clear();
    src.clear();
  };

  void build(vec<int> &sa, string s)
  {
    int N = (int)sa.size();
    src.resize(N, 0);
    up.resize(N, -1);
    next.resize(N, -1);
    down.resize(N, -1);
    init_kasai(sa, s);
    init_updown();
    init_next();
  }
  
  int &operator[](size_t i) { return src[i]; }

  /**
   * @brief Returns the count of common prefixes in λ[lo..hi] 
   */
  int &range(int lo, int hi)
  {
    return (lo < up[hi + 1] && up[hi+1] < hi+1) ? src[up[hi + 1]] : src[down[lo]];
  }
};

/**
 * @brief Represents the suffix node λ[lo..hi]
 */
class SuffixNode
{
public:
  int lo,  // lowest index
      hi;  // highest index

  SuffixNode(int lo, int hi) : lo(lo), hi(hi)
  {}

  bool isLeaf() { return lo == hi; }
  
  friend ostream& operator<<(ostream& os, const SuffixNode nd) {
    return os << "λ[" << nd.lo << ".." << nd.hi << "]";
  }
};

class SuffixArray {
  private:
    // subtree roots ordered from first to last 
    // contains O(α) entries in the worst case
    vec<SuffixNode> _root; 

    /**
     * @brief Constructs the suffix array from the input string in O(N*log^2(N))
     */
    void build() {
      struct suffix {
        int rank[2]; // alphabetical rank of prv:0 and cur:1
        int i;       // p_idx;
      };
      
      auto cmp = [](struct suffix &a, struct suffix &b) {
        return a.rank[0] == b.rank[0] ? 
          (a.rank[1] < b.rank[1] ? 1 : 0) : 
          (a.rank[0] < b.rank[0] ? 1 : 0);
      };

      int i;
      const int N = (int)s.size();
      vec<suffix> suf(N, suffix{});       // for storing suffixes
      vec<vec<int>> P(1, vec<int>(N, 0)); // init length 1

      // if only lowercase letters, change to s[i]-'a'+1
      FOR(i,0,(int)N) P[0][i] = s[i]; // step 0

      // P[step][i] stores the position of the i-th longest suffix
      for (int step = 1, len = 1; 
          len < (int) (2*N);        // log(N) iterations
          step++, len <<= 1)  
      {                     
        P.PB(vec<int>(N,0));

        FOR(i,0,N) {          // init entries 
          suf[i].rank[0] = P[step - 1][i];  
          suf[i].rank[1] = (i + len) < (int) N ? P[step - 1][i + len] : -1; 
          suf[i].i = i;
        }

        sort(ALL(suf), cmp);   // replace with radix sort to improve running time to O(N*log(N))

        FOR(i,1,N) {           // Update P[step]
          if (suf[i].rank[0] == suf[i-1].rank[0] 
          && suf[i].rank[1] == suf[i-1].rank[1])
            P[step][suf[i].i] = P[step][suf[i - 1].i];
          else     
            P[step][suf[i].i] = i;
        }    
      }

      // init suffix array
      FOR(i,0,N) sa.push_back(suf[i].i);  
    }

    vec<SuffixNode> find_kids(int lo, int hi) {
      vec<SuffixNode> ret;
      SuffixNode nd(lo, hi);
      if (nd.isLeaf()) { return ret; }
      // If up is within interval start at up, else start at down
      int i = (lo < lcp.up[hi+1] && lcp.up[hi+1] <= hi) ?
              lcp.up[hi+1] : lcp.down[lo];
      ret.emplace_back(lo, i-1);
      for (int j = lcp.next[i];
           j != -1 && i != j;
           i = j, j = lcp.next[i])
        ret.emplace_back(i, j-1);
      ret.emplace_back(i, hi);
      return ret;
    }

    // O(α) init intervals of suffix array
    void init_root() {
      int i, j, N = sa.size();
      for (i = 0, j = lcp.next[i];
            j != -1;
            i = j, j = lcp.next[i])
        _root.emplace_back(i, j - 1);
      _root.emplace_back(i, N-1);
    };

  public:
    string s;     // src string
    vec<int> sa;  // suffix array - consider changing to std::array for performance
    LongestCommonPrefixArray lcp;

    SuffixArray(string& s) : s(s), lcp() {
      this->build();
      lcp.build(sa, s);
      init_root();
    }
    ~SuffixArray() {
      sa.clear();
    }

    /**
     * @brief find all matches to s2 in suffix array in O(M*α+k) 
     *        where k = number of matches
     *              α = size of string alphabet 
     *              M = size of suffix array i.e. length of src string
     * @param s2 needle to search for
     * @return pair {lo,hi} containing indices denoting the bounds of all matches to needle
     *         occurrences of a particular substring will appear as a contiguous block in the suffix array, 
               so the width of this block is the number of occurrences.
    */
    pi cfind(const string &s2) {
      // O(hi-lo) auxiliary method to compare needle[lo..hi] with S[su[λ[lo]+lo..su[λ[lo]]+hi]
      auto str_rng_eq = [this, &s2](SuffixNode &node, int lo, int hi) -> bool
      {
        auto S1 = s.c_str(), S2 = s2.c_str();
        const int k = hi - lo;
        return equal(S1 + sa[node.lo] + lo, S1 + sa[node.lo] + lo + k,
                     S2 + lo, S2 + lo + k);
      };
      
      vec<SuffixNode> kids(_root);
      vec<SuffixNode>::iterator kid = kids.beg();
      int M = (int)s2.length(), m = 0;
      
      while (m < M && kid != kids.end()) 
      {
        // if (*kid->hi)
        // cout << "Visiting: " << *kid << endl;
        int hi = kid->isLeaf() ? M : min(lcp.range(kid->lo, kid->hi), M);
        // if range s1[sa[lo] + m.. sa[lo]+hi] equals s2[m..hi]
        if (str_rng_eq(*kid, m, hi)) {
          if (kid->isLeaf()) return pi{kid->lo, kid->hi};
          // iterate down the suffix tree
          m = hi; 
          if (m < M) {
            kids = find_kids(kid->lo, kid->hi);
            kid = kids.beg();
          }
        } else {
          ++kid;
        }
      }
      
      if (kid == kids.end()) return NONE_PI;
      return pi{kid->lo, kid->hi};
    }

    void print() {
      printf("Printing suffix array #%d\n ", size());
      printf("n  | sa[n] | lcp[n] | up[n] | down[n] | next[n] | suffix[n]\n ");
      for (int i = 0; i < (int) s.size(); ++i)
        printf("%-3d   %-3d      %-3d      %-3d      %-3d      %-3d     %s\n ", 
        i, sa[i], lcp[i], lcp.up[i], lcp.down[i], lcp.next[i], this->get(sa[i]).c_str());
      cout.flush();
      // int i, j, k;
      cout << "\n";
      auto print_q = [this](vec<SuffixNode>::iterator qbeg, 
                            vec<SuffixNode>::iterator qend, 
                            int lo, int hi)
      {
        cout << SuffixNode(lo, hi) << "="
             << lcp.range(lo,hi) << " { ";
        copy(qbeg, qend, os_joiner(", "));
        cout << " }\n";
      };

      printf("Printing λ-suffix-array tree #%d\n", size());
      print_q(ALL(_root), 0, (int)s.size() - 1);
      for(auto r : _root) {
        auto kids = find_kids(r.lo, r.hi);
        print_q(ALL(kids), r.lo, r.hi);
        for(auto kid : kids) {
          auto kids2 = find_kids(kid.lo,kid.hi);
          print_q(ALL(kids2), kid.lo,kid.hi);
        }
      }
    }

    int size() { return s.size(); }
    int& operator[](size_t i) { return sa[i]; }
    string get(size_t i) { return s.substr(i); }
};

MultiString mstr;

int main() {
  fast_io();

  int N, // Count of entries
      Q, // Count of queries
      n, q; 

  string s, needle;
  cin >> N;
  
  FOR(n,0,N) {
    cin >> s;
    mstr.append(s);
  }

  s = mstr.c_str();
  SuffixArray sa(s);
  cin >> Q;

  FOR(q,0,Q) {
    cin >> needle;
    set<int> ans; // ordered set ascending
    auto [lo, hi] = sa.cfind(needle);
    if (lo == -1 || hi == -1) // if not found
    {
      cout << 0 << endl;
      continue;
    }
    for (int i = lo; i <= hi; i++) {
      ans.insert(mstr.strnum(sa[i]));
    }
    cout << ans.size() << " ";
    copy(ALL(ans), os_joiner(" "));
    cout << "\n";
  }
  return PROGRAM_SUCCESS;
}