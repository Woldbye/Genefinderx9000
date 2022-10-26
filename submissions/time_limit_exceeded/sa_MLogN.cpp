// #TLE solution for genefinder
// Dictionary containing N strings build in O(|N|*log|N|)
//            where |N| denotes the combined length of the N strings
// Queries of length M answered in O(M*logN)
#include <bits/stdc++.h>
#include <functional>
#include <experimental/iterator>
using namespace std;

#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,fma")
#pragma GCC optimize("unroll-loops")
#define fast_io() ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);

#define ALL(o) (o).begin(), (o).end() 
#define FOR(iter, star, finish) for(iter = star; iter < finish; ++iter)
#define ODD(i) (bool) (i & 0x1)
#define EVEN(i) !ODD(i)
#define vec vector
#define PROGRAM_SUCCESS 0
#define PB push_back
#define os_joiner(del) experimental::make_ostream_joiner(cout, del)
#define PROGRAM_SUCCESS 0

typedef unsigned long long int ui64;
typedef long long int i64;
typedef pair<int, int> pi;

static const char SEPERATOR = '~';
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

// https://web.stanford.edu/class/cs97si/suffix-array.pdf
// https://github.itu.dk/algorithms/aps-2021/blob/master/week05-string-algorithms/lecture_notes.md
class SuffixArray {
  private:
    const string s;  // src string
    void build_lcp(); // O(N) : build longest-common-prefix array
    void build_sa();  // O(N*log(N)^2) build suffix array

  public:
    vec<int> sa;     // suffix array
    vec<int> lcpa;   // lcpa[n] : lcpa[n] = length of the longest-common-prefix between strings (sa[n], sa[n-1])

    SuffixArray(string& s) : s(s) {
      this->build_sa();
      this->build_lcp();
    }
    ~SuffixArray() {
      sa.clear();
      lcpa.clear();
    }

    string get(size_t i) { return s.substr(i); } 

    // O(M*log2(N)) : Find the index of lower bound match to needle
    //! change to knuth-morris-pratt algorithm to achieve O(N+M)
    int bfind(const string& s2) {
      const int M = (int) s2.size(), N = (int) sa.size(); 
      
      // Do simple binary search for the pat in txt using the
      // built suffix array
      int lo = 0, hi = N-1;  // Initialize left and right indexes
      while (lo <= hi)
      {
          // See if 'pat' is prefix of middle suffix in suffix array
          int mid = lo + (hi - lo)/2;
          int res = strncmp(s2.c_str(), s.c_str()+sa[mid], M);
  
          // If match found at the middle, print it and return
          if (res == 0) return mid;
          if (res < 0) hi = mid - 1;
          else lo = mid + 1;
      }
      return -1;
    }

    pi find_bounds(int i, const string& needle) {
      // printf("find_bounds(src#%s i#%d needle#%s)\n", s.c_str(), i, needle.c_str());
      const int M = (int)needle.size(),
                N = (int)sa.size();
      int hi,lo;
      for (hi = i+1; hi < N  && lcpa[hi]   >= M; ++hi);
      for (lo = i-1; lo >= 0 && lcpa[lo+1] >= M; --lo);
      return {lo+1, hi-1};
    }
    
    void print() {
      printf("Printing suffix array #%d\n ", size());
      printf("n  | sa[n] | lcp[n] | suffix[n]\n ");
      for (int i = 0; i < (int) s.size(); ++i)
        printf("%-3d   %-3d      %-3d     %s\n ", i, sa[i], lcpa[i], this->get(sa[i]).c_str());
      cout.flush();
    }

    int size() { return s.size(); }
    int& lcp(size_t i) { return lcpa[i]; }
    int& operator[](size_t i) { return sa[i]; }
};

// Build longest-common-prefix array in O(N) : Kasai's Algorithm
// https://www.codingninjas.com/codestudio/library/longest-common-prefix-from-suffix-array
// https://stackoverflow.com/questions/57748988/kasai-algorithm-for-constructing-lcp-array-practical-example
void SuffixArray::build_lcp()
{
  const int N = (int)s.size();
  int n, len = 0;
  
  // Stores the inverse index of the suffix array 
  // rank[n] = the position of suffix[n] in sa
  vec<int> rank(N, 0);  
  lcpa.resize(N, 0);
  
  FOR(n,0,N) rank[sa[n]] = n;

  FOR(n,0,N) {
    if (rank[n] > 0) { // suffix n has 'prv' in the suffix array.
      for (int prv = sa[rank[n] - 1]; // prv = preceding suffix
           s[n + len] == s[prv + len];
           len++);
      lcpa[rank[n]] = len;
      if (len > 0) len--;   // if len > 0 -1 else 0
    }
  }
}

/**
 * @brief Constructs the suffix array from the input string in O(N*log^2(N))
 */
void SuffixArray::build_sa() {
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
  int step = 1, len = 1;
  for (; 
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

    FOR(i,1,N) {         // Update P[step]
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

class Dictionary
{
private:
  SuffixArray* sa;
  MultiString mstr;

public:
  void init(vec<string> &genes) {
    string str;
    for (auto &s : genes)
      mstr.append(s);
    str = mstr.c_str();
    sa = new SuffixArray(str);
  };

  set<int> find_all(const string &needle) {
    set<int> ret; ui64 j;
    int i = sa->bfind(needle);
    if (i != -1) {
      auto [lo, hi] = sa->find_bounds(i, needle);

      // occurrences of a particular substring will appear as a 
      // contiguous block in the suffix array, 
      // so the width of this block is the number of occurrences.
      for (int i = lo; i <= hi; i++) {
        ret.insert(mstr.strnum((*sa)[i]));
      }
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