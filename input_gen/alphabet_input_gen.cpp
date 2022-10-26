/**
 * @brief Generates random input of desired length for genematcherx9000.
 *        The input is written to stdout. 
 *        Bash usage example:
 *        ./in_gen 1000 5 > sample.in 
 * @warning can take quite a long time to generate for large numbers
 */
#include <bits/stdc++.h>
#include <functional>
#include <experimental/iterator>
#include <chrono>
using namespace std;
using namespace chrono;

#define ALL(o) (o).begin(), (o).end()
#define beg() begin()
#define vec vector
#define os_joiner(del) experimental::make_ostream_joiner(cout, del)

typedef unordered_set<string> hmap;

const int ARGS = 5,
          MIN_ARGS = 2;
vec<char> alphabet = {'A', 'C', 'G', 'T'};

int fail()
{
  printf("Expected %d non-negative integer arguments: N Q\n", ARGS-1);
  printf("1: N   - Maximum length of the combined animal database\n");
  printf("2: M   - Maximum length of query strings, expected to be less than nL\n");
  printf("3: RND - 3 if only N random \n 2 if only M random \n 1 if both random\n0 if minimum length\nDefault: RND=0 \n");
  return 1;
}

// Generate a database of length N where all entries are random sized
void generate_random(int N) {
  hmap db;
  auto rnd_char = [](vec<char>& alpha) { return alpha[rand() % alpha.size()]; };
  
  int n = 0, fail = 0; // Use failure counter: we might be stuck in a situation where it can only construct length 1, but all possibilities are already pushed
  while (n < N && fail < 10) {
    string animal;
    int len = rand() % (N-n) + 1; 
    for (int i = 0; i < len; ++i) animal.push_back(rnd_char(alphabet));
    if (!db.count(animal)) {
      db.emplace(animal);
      n += animal.length();
      fail = 0;
    } else {
      ++fail;
    }
  }

  cout << db.size() << "\n";
  for (auto& animal : db) cout << animal << "\n";
}

// Generates a random database of maximum length N. 
// The database will contain the animals in the database will each be
// minimum length, such that no animals are equal.
// Therefore if an animal of length k is in the database, all permutations of animals 
// for length  [1..k-1] is also in the DB
void generate_minimum(int N) {
  vec<string> out;
  vec<vec<string>> db(1, vec<string>(1,""));
  int n = 0;
  auto db_add = [&db, &n](int l, string &s)
  {
    db[l].push_back(s);
    n += s.length();
  };

  for (int l = 1; n + l < N; ++l){
    db.push_back(vec<string>());
    
    for(string& s : db[l-1]) {
      // Push back all permutations
      for (size_t ci = 0; ci <  alphabet.size() && (n+l < N); ci++) {
        char c = alphabet[ci];
        string cur = s + c;
        db_add(l, cur);
      }
    }    
  }

  // Flatten db into out
  for (size_t l = 1; l < db.size(); ++l) out.insert(out.end(), ALL(db[l]));
  random_shuffle(ALL(out)); // Randomize order

  cout << out.size() << "\n";
  for (auto& animal : out) cout << animal << "\n";
}

const int RND = -1, MINIMUM = -2;

/**
 * @brief Generate and print a random database of strings consisting of the alphabet
 *        defined and of combined length N.
 *        if I=RND          -> initializes random length entries
 *        else if I=MINIMUM -> initializes minimum length entries (note combined length might not be exactly N)
 **/
void generate(bool isRandomLength, int N) {
  if (isRandomLength) generate_random(N);
  else {
    generate_minimum(N);
  } 
}

int main(int argc, char** argv)
{
  srand(time(NULL)); // Generate random seed

  if (argc < MIN_ARGS+1) return fail();

  int N = stoi(argv[1]),
      M = stoi(argv[2]),
      RND = (argv[3] != NULL) ? stoi(argv[3]) : 0;
  bool isRnd = RND != 0;
  // auto start = high_resolution_clock::now();
  generate(isRnd && RND != 2, M);
  generate(isRnd && RND != 3, N);
  // auto stop = high_resolution_clock::now();
  // cout << "Time taken: " << duration_cast<microseconds>(stop - start).count() << "ms"<< endl;
  return 0;
}