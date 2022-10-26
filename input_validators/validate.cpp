#include <bits/stdc++.h>
#include <experimental/iterator>
using namespace std;

#define vec vector
#define FOR(it, start, stop) for (it = start; it < stop; ++it)

#define SUCCESS_EXIT 42

const int MAX_LENGTH = 10 ^ 5;
const int MAX_N = 10 ^ 5;
const int MAX_M = 10 ^ 6;

inline bool isGene(const char c)
{
  switch (c) {
    case 'A': return true;
    case 'G': return true;
    case 'T': return true;
    case 'C': return true;
    default: return false;
  }
}

inline bool isValidEntry(const string& s) {
  if (s.length() == 0) // No empty lines
    return false;  
  for (const char c : s) {
    if (!isGene(c)) return false;
  }
  return true;
}

int main(int argc, char** argv) {
  string line;
  int K;
  assert(cin.peek() != '0');
  for (int i = 0, len; i < 2; ++i)
  {
    assert(cin >> K);
    assert(K > 0);
    int l = 0;
    cin.ignore(); // ignore \n
    set<string> ss;
    for (int k = 0; k < K; ++k) {
      getline(cin, line);
      assert(isValidEntry(line));
      ss.insert(line);
      l += line.length();
    }
    assert(ss.size() == K); // Assert no duplicates
    assert(l <= (!i) ? MAX_N : MAX_M);
    ++i;
  }
  string err_msg = "\"";
  if (!cin)
  {
    string tmp;
    while (cin >> tmp) {
      cin >> tmp;
      err_msg += tmp;
    }
  }
  err_msg += "\"";
  assert(err_msg.length() == 2);
  return SUCCESS_EXIT;
}