/*
 Hierarchically clusters phrases.
 Running time: O(N*C^2).
 
 We want to cluster the phrases so that the pairwise mutual information between
 clusters is maximized.  This mutual information is a sum over terms between
 each pair of clusters: q2[a, b] for clusters a and b.  The trick is to compute
 quickly the loss of mutual information when two clusters a and b are merged.
 
 The four structures p1, p2, q2, L2 allow this quick computation.
 p1[a] = probability of of cluster a.
 p2[a, b] = probability of cluster a followed by cluster b.
 q2[a, b] = contribution to the mutual information from clusters a and b (computed from p2[a, b]).
 L2[a, b] = the loss of mutual information if clusters a and b were merged.
 
 Changes:
 * Removed hash tables for efficiency.
 * Notation: a is an phrase (sequence of words), c is a cluster, s is a slot.
 
 To cut down memory usage:
 * Change double to float.
 Ideas:
 * Hashing vectors is really slow.
 * Find intuition behind algorithm based on simple cases
 * Test clustering algorithm on artificial generated data.  Generate a text
 with a class-based ngram model.
 */
#include <Rcpp.h>

#include "brown/basic/std.h"
#include "brown/basic/stl-basic.h"
#include "brown/basic/stl-utils.h"
#include "brown/basic/str.h"
#include "brown/basic/strdb.h"
#include "brown/basic/union-set.h"
#include "brown/basic/mem-tracker.h"
//#include "brown/basic/opt.h"
#ifndef _WIN32
#include <unistd.h>
#endif
#include <condition_variable>
#include <mutex>
#include <thread>

#ifdef _WIN32
#include <windows.h>
int getpid()
{
  return GetCurrentProcessId();
}
#endif

class semaphore
{
private:
  mutex _m;
  condition_variable _cv;
  int _count;
  
public:
  semaphore(int count = 0) : _count(count) {}
  
  inline void notify()
  {
    unique_lock<mutex> lock(_m);
    ++_count;
    _cv.notify_one();
  }
  
  inline void wait()
  {
    unique_lock<mutex> lock(_m);
    _cv.wait(lock, [this]() { return _count > 0; });
    --_count;
  }
};

class BrownClustering{
  
public:
};

const char *delim_str = "$#$";

typedef IntPair _;

StrDB db;                // word database
IntVec phrase_freqs;     // phrase a < N -> number of times a appears in the text
IntVecVec left_phrases;  // phrase a < N -> list of phrases that appear to left of a in the text
IntVecVec right_phrases; // phrase a < N -> list of phrases that appear to right of a in the text
IntIntPairMap cluster_tree; // cluster c -> the 2 sub-clusters that merged to create c
int delim_word;

IntVec freq_order_phrases; // List of phrases in decreasing order of frequency.

// Allows for very quick (inverse Ackermann) lookup of clusters and merging
// of clusters.  Each phrase points to an arbitrary representative phrase of
// the cluster.
UnionSet phrase2rep;   // phrase a -> the rep phrase in the same cluster as a
IntIntMap rep2cluster; // rep phrase a -> the cluster that contains a
IntIntMap cluster2rep; // cluster a -> the rep phrase in cluster a

// Store all the phrases efficiently.  Just for printing out.
// For each phrase length, we store a flattened list of words.
IntVecVec phrases; // length of phrase -> flattened list of words

// Each cluster will occupy a slot.  There will always be two extra slots
// as intermediate scratch space.
IntVec slot2cluster; // slot index -> cluster (-1 if none exists)
IntIntMap cluster2slot; // cluster -> slot index
int free_slot1, free_slot2; // two free slots
int nslots;

// Partial results that allow quick computation and update of mutual information.
// Mutual information is the sum of all the q2 terms.
// Update p1, p2, q2 for 0..N-1, but L2 only for 0..initC-1.
DoubleVec p1;        // slot s (containing cluster a) -> probability Pr(a)
DoubleVecVec p2;     // slots s, t (containing clusters a, b) -> probability Pr(a, b)
DoubleVecVec q2;     // slots s, t (contianing clusters a, b) -> contribution to mutual information
DoubleVecVec L2;     // slots s, t (containing clusters a, b) -> loss of mutual information if merge a and b

int curr_cluster_id; // ID to assign to a new cluster
int stage2_cluster_offset; // start of the IDs of clusters created in stage 2

double curr_minfo; // Mutual info, should be sum of all q2's

// Map phrase to the KL divergence to its cluster
DoubleVec kl_map[2];

// Variables used to control the thread pool
vector<semaphore> thread_idle;
vector<semaphore> thread_start;
thread * threads;
struct Compute_L2_Job {
  int s;
  int t;
  int u;
  bool is_type_a;
};
Compute_L2_Job the_job;
bool all_done = false;

#define FOR_SLOT(s)                        \
for(int s = 0; s < len(slot2cluster); s++) \
  for(bool _tmp = true; slot2cluster[s] != -1 && _tmp; _tmp = false)

// We store only L2[s, t] for which the cluster ID in slot s is smaller
// than the one in slot t.
#define ORDER_VALID(s, t) (slot2cluster[s] < slot2cluster[t])

#define num_phrases(l) (len(phrases[l])/(l))

int N; // number of phrases
int T; // length of text

// Output a phrase.
struct Phrase { Phrase(int a) : a(a) { } int a; };
// For pretty-printing of clusters.
struct Cluster { Cluster(int c) : c(c) { } int c; };
std::string PhraseText(const Phrase &phrase) {
  std::string out;
  // Decode the phrase ID into the length and the offset in phrases.
  int a = phrase.a;
  int l; for(l = 1; a >= num_phrases(l); a -= num_phrases(l), l++);
  
  foridx(i, l) {
    if(i > 0) out += ' ';
    out += db[phrases[l][a*l+i]];
  }
  return out;
}
std::string ClusterText(const Cluster &cluster){
  std::string out;
  int c = cluster.c;
  out += std::to_string(c);
  
  int a;
  bool more;
  if(c < N)
    a = c, more = false;
  else {
    assert(contains(cluster2rep, c));
    a = cluster2rep[c], more = true;
  }
  out += '(';
  out += PhraseText(Phrase(a));
  if(more) out += "|...";
  out += ')';
  return out;
}

#define Slot(s) Cluster(slot2cluster[s])



////////////////////////////////////////////////////////////

// p2[s, t] + p2[t, s].
inline double bi_p2(int s, int t) {
  if(s == t) return p2[s][s];
  return p2[s][t] + p2[t][s];
}

// q2[s, t] + q2[t, s].
inline double bi_q2(int s, int t) {
  if(s == t) return q2[s][s];
  return q2[s][t] + q2[t][s];
}

// Hypothetical p1[st] = p1[s] + p1[t].
inline double hyp_p1(int s, int t) {
  return p1[s] + p1[t];
}

//// hyp_p2

// Hypothetical p2[st, u] = p2[s, u] + p2[t, u].
inline double hyp_p2(const IntPair &st, int u) {
  return p2[st.first][u] + p2[st.second][u];
}

// Hypothetical p2[u, st] = p2[u, s] + p2[u, t].
inline double hyp_p2(int u, const IntPair &st) {
  return p2[u][st.first] + p2[u][st.second];
}

inline double bi_hyp_p2(const IntPair &st, int u) {
  return hyp_p2(st, u) + hyp_p2(u, st);
}

// Hypothetical p2[st, st] = p2[s, s] + p2[s, t] + p2[t, s] + p2[t, t].
inline double hyp_p2(const IntPair &st) {
  return p2[st.first][st.first] + p2[st.first][st.second] +
    p2[st.second][st.first] + p2[st.second][st.second];
}

//// hyp_q2

inline double p2q(double pst, double ps, double pt) {
  if(feq(pst, 0.0)) return 0.0;
  return pst * log2(pst / (ps*pt));
}

// Hypothetical q2[st, u].
inline double hyp_q2(const IntPair &st, int u) {
  return p2q(hyp_p2(st, u), hyp_p1(st.first, st.second), p1[u]);
}

// Hypothetical q2[u, st].
inline double hyp_q2(int u, const IntPair &st) {
  return p2q(hyp_p2(u, st), hyp_p1(st.first, st.second), p1[u]);
}

inline double bi_hyp_q2(const IntPair &st, int u) {
  return hyp_q2(st, u) + hyp_q2(u, st);
}

// Hypothetical q2[st, st].
inline double hyp_q2(const IntPair &st) {
  double p = hyp_p2(make_pair(st.first, st.second)); // p2[st,st]
  double P = hyp_p1(st.first, st.second);
  return p2q(p, P, P);
}

////////////////////////////////////////////////////////////

// Return slot.
void put_cluster_in_slot(int a, int s) {
  cluster2slot[a] = s;
  slot2cluster[s] = a;
}
inline int put_cluster_in_free_slot(int a) {
  int s = -1;
  
  // Find available slot.
  if(free_slot1 != -1)      s = free_slot1, free_slot1 = -1;
  else if(free_slot2 != -1) s = free_slot2, free_slot2 = -1;
  assert(s != -1);
  
  put_cluster_in_slot(a, s);
  return s;
}

inline void free_up_slots(int s, int t) {
  free_slot1 = s;
  free_slot2 = t;
  cluster2slot.erase(slot2cluster[s]);
  cluster2slot.erase(slot2cluster[t]);
  slot2cluster[s] = slot2cluster[t] = -1;
}

void init_slot(int s) {
  // Clear any entries that relates to s.
  // The p1 and L2 will be filled in densely, so they
  // will be overwritten anyway.
  FOR_SLOT(t)
  p2[s][t] = q2[s][t] = p2[t][s] = q2[t][s] = 0;
}

void add_to_set(const IntVec &phrases, IntIntMap &phrase_counts, int offset) {
  forvec(_, int, a, phrases)
  phrase_counts[a+offset]++;
}

bool is_good_phrase(const IntVec &phrase) {
  if(len(phrase) == 1) return phrase[0] != delim_word && phrase[0] != -1; // Can't be delimiter or an invalid word
  
  // HACK HACK HACK - pick out some phrases
  // Can't be too many delim words.
  int di = index_of(phrase, delim_word, 1);
  if(di > 0 && di < len(phrase)-1) return false; // Delimiter must occur at the ends
  if(phrase[0] == delim_word && phrase[len(phrase)-1] == delim_word) return false; // Only one delimiter allowed
  
  // If every word is capitalized with the exception of some function
  // words which must go in the middle
  forvec(i, int, a, phrase) {
    bool at_end = i == 0 || i == len(phrase)-1;
    const string &word = db[a];
    bool is_upper = isupper(word[0]);
    
    if(at_end && !is_upper) return false; // Ends must be uppercase
    if(is_upper) continue; // Ok
    if(word[0] == '\'' || word == "of" || word == "and") continue; // Ok
    return false;
  }
  return true;
}


void read_text(bool paths2map, int plen, int min_occur, int initC, IntVec text) {
  IntVecIntMap vec2phrase;
  T = len(text);
  delim_word = db.lookup(delim_str, false, -1);
  //if(!paths2map) db.destroy_s2i(); // Conserve memory.
  
  // Count the phrases that we care about so we can map them all to integers.
  track_block("Counting phrases", "", false) {
    phrases.resize(plen+1);
    for(int l = 1; l <= plen; l++) {
      // Count.
      IntVecIntMap freqs; // phrase vector -> number of occurrences
      for(int i = 0; i < T-l+1; i++) {
        IntVec a_vec = subvector(text, i, i+l);
        if(!is_good_phrase(a_vec)) continue;
        freqs[a_vec]++;
      }
      
      forcmap(const IntVec &, a_vec, int, count, IntVecIntMap, freqs) {
        if(count < min_occur) continue;
        
        int a = len(phrase_freqs);
        phrase_freqs.push_back(count);
        vec2phrase[a_vec] = a;
        forvec(_, int, w, a_vec) phrases[l].push_back(w);
      }
      
      logs(len(freqs) << " distinct phrases of length " << l << ", keeping " << num_phrases(l) << " which occur at least " << min_occur << " times");
    }
  }
  
  N = len(phrase_freqs); // number of phrases
  
  track_block("Finding left/right phrases", "", false) {
    left_phrases.resize(N);
    right_phrases.resize(N);
    for(int l = 1; l <= plen; l++) {
      for(int i = 0; i < T-l+1; i++) {
        IntVec a_vec = subvector(text, i, i+l);
        if(!contains(vec2phrase, a_vec)) continue;
        int a = vec2phrase[a_vec];
        
        // Left
        for(int ll = 1; ll <= plen && i-ll >= 0; ll++) {
          IntVec aa_vec = subvector(text, i-ll, i);
          if(!contains(vec2phrase, aa_vec)) continue;
          int aa = vec2phrase[aa_vec];
          left_phrases[a].push_back(aa);
          //logs(i << ' ' << Cluster(a) << " L");
        }
        
        // Right
        for(int ll = 1; ll <= plen && i+l+ll <= T; ll++) {
          IntVec aa_vec = subvector(text, i+l, i+l+ll);
          if(!contains(vec2phrase, aa_vec)) continue;
          int aa = vec2phrase[aa_vec];
          right_phrases[a].push_back(aa);
          //logs(i << ' ' << Cluster(a) << " R");
        }
      }
    }
  }
  
  //destroy(text);
  initC = min(initC, N);
  
  logs("Text length: " << T << ", " << N << " phrases, " << len(db) << " words");
}

// O(C) time.
double compute_s1(int s) { // compute s1[s]
  double q = 0.0;
  
  for(int t = 0; t < len(slot2cluster); t++) {
    if (slot2cluster[t] == -1) continue;
    q += bi_q2(s, t);
  }
  
  return q;
}

// O(C) time.
double compute_L2(int s, int t) { // compute L2[s, t]
  assert(ORDER_VALID(s, t));
  // st is the hypothetical new cluster that combines s and t
  
  // Lose old associations with s and t
  double l = 0.0;
  for (int w = 0; w < len(slot2cluster); w++) {
    if ( slot2cluster[w] == -1) continue;
    l += q2[s][w] + q2[w][s];
    l += q2[t][w] + q2[w][t];
  }
  l -= q2[s][s] + q2[t][t];
  l -= bi_q2(s, t);
  
  // Form new associations with st
  FOR_SLOT(u) {
    if(u == s || u == t) continue;
    l -= bi_hyp_q2(make_pair(s, t), u);
  }
  l -= hyp_q2(make_pair(s, t)); // q2[st, st]
  return l;
}

void repcheck(bool chk) {
  if(!chk) return;
  double sum;
  
  assert_eq(len(rep2cluster), len(cluster2rep));
  assert_eq(len(rep2cluster), len(cluster2slot));
  
  assert(free_slot1 == -1 || slot2cluster[free_slot1] == -1);
  assert(free_slot2 == -1 || slot2cluster[free_slot2] == -1);
  FOR_SLOT(s) {
    assert(contains(cluster2slot, slot2cluster[s]));
    assert(cluster2slot[slot2cluster[s]] == s);
  }
  
  sum = 0.0;
  FOR_SLOT(s) FOR_SLOT(t) {
    double q = q2[s][t];
    //logs(s << ' ' << t << ' ' << p2[s][t] << ' ' << p1[s] << ' ' << p1[t]);
    assert_feq(q, p2q(p2[s][t], p1[s], p1[t]));
    sum += q;
  }
  assert_feq(sum, curr_minfo);
  
  FOR_SLOT(s) FOR_SLOT(t) {
    if(!ORDER_VALID(s, t)) continue;
    double l = L2[s][t];
    assert(l + TOL >= 0);
    assert_feq(l, compute_L2(s, t));
  }
}



// c is new cluster that has been just formed from a and b
// Want to compute L2[d, e]
// O(1) time.
double compute_L2_using_old(int s, int t, int u, int v, int w) {
  assert(ORDER_VALID(v, w));
  assert(v != u && w != u);
  
  double l = L2[v][w];
  
  // Remove old associations between v and w with s and t
  l -= bi_q2(v, s) + bi_q2(w, s) + bi_q2(v, t) + bi_q2(w, t);
  l += bi_hyp_q2(make_pair(v, w), s) + bi_hyp_q2(make_pair(v, w), t);
  
  // Add new associations between v and w with u (ab)
  l += bi_q2(v, u) + bi_q2(w, u);
  l -= bi_hyp_q2(make_pair(v, w), u);
  
  return l;
}

// return q2
double set_p2_q2_from_count(int s, int t, int count) {
  double pst = (double)count / (T-1); // p2[s,t]
  double ps = p1[s];
  double pt = p1[t];
  double qst = p2q(pst, ps, pt); // q2[s,t]
  p2[s][t] = pst;
  q2[s][t] = qst;
  return qst;
}

// O(N lg N) time.
// Sort the phrases by decreasing frequency and then set the initC most frequent
// phrases to be in the initial cluster.
bool phrase_freq_greater(int a, int b) {
  return phrase_freqs[a] > phrase_freqs[b];
}
void create_initial_clusters(int initC) {
  // track("create_initial_clusters()", "", true);
  // 
  // freq_order_phrases.resize(N);
  // foridx(a, N) freq_order_phrases[a] = a;
  // 
  // logs("Sorting " << N << " phrases by frequency");
  // sort(freq_order_phrases.begin(), freq_order_phrases.end(), phrase_freq_greater);
  // 
  // if(initC > N){
  //   initC = N;
  // }
  
  // Initialize slots
  logs("Selecting top " << initC << " phrases to be initial clusters");
  nslots = initC+2;
  slot2cluster.resize(nslots);
  free_up_slots(initC, initC+1);
  
  // Create the inital clusters.
  phrase2rep.Init(N); // Init union-set: each phrase starts out in its own cluster
  curr_minfo = 0.0;
  foridx(s, initC) {
    int a = freq_order_phrases[s];
    put_cluster_in_slot(a, s);
    
    rep2cluster[a] = a;
    cluster2rep[a] = a;
  }
  
  // Allocate memory
  p1.resize(nslots);
  matrix_resize(p2, nslots, nslots);
  matrix_resize(q2, nslots, nslots);
  matrix_resize(L2, nslots, nslots);
  
  FOR_SLOT(s) init_slot(s);
  
  // Compute p1
  FOR_SLOT(s) {
    int a = slot2cluster[s];
    p1[s] = (double)phrase_freqs[a] / T;
  }
  
  // Compute p2, q2, curr_minfo
  FOR_SLOT(s) {
    int a = slot2cluster[s];
    IntIntMap right_phrase_freqs;
    
    // Find collocations of (a, b), where both are clusters.
    forvec(_, int, b, right_phrases[a])
      if(contains(cluster2slot, b))
        right_phrase_freqs[b]++;
      
      forcmap(int, b, int, count, IntIntMap, right_phrase_freqs) {
        int t = cluster2slot[b];
        curr_minfo += set_p2_q2_from_count(s, t, count);
      }
  }
}

// Output the ncollocs bigrams that have the highest mutual information.
void output_best_collocations(string collocs_file, int ncollocs) {
  if(collocs_file.empty()) return;
  logs("Writing to " << collocs_file);
  
  vector< pair<double, IntPair> > collocs;
  FOR_SLOT(s) FOR_SLOT(t) {
    collocs.push_back(pair<double, IntPair>(q2[s][t], make_pair(slot2cluster[s], slot2cluster[t])));
  }
  ncollocs = min(ncollocs, len(collocs));
  partial_sort(collocs.begin(), collocs.begin()+ncollocs, collocs.end(), greater< pair<double, IntPair> >());
  
  ofstream out(collocs_file.c_str());
  assert(out);
  for(int i = 0; i < ncollocs; i++) {
    const IntPair &ab = collocs[i].second;
    out << collocs[i].first << '\t' << PhraseText(Phrase(ab.first)) << '\t' << PhraseText(Phrase(ab.second)) << endl;
  }
}

// O(C^3) time.
void compute_L2() {
  track("compute_L2()", "", true);
  
  track_block("Computing L2", "", false)
    FOR_SLOT(s) {
      track_block("L2", "L2[" << ClusterText(Slot(s)) << ", *]", false)
      FOR_SLOT(t) {
        if(!ORDER_VALID(s, t)) continue;
        double l = L2[s][t] = compute_L2(s, t);
        logs("L2[" << ClusterText(Slot(s)) << "," << ClusterText(Slot(t)) << "] = " << l << ", resulting minfo = " << curr_minfo-l);
      }
    }
}

// Add new phrase as a cluster.
// Compute its L2 between a and all existing clusters.
// O(C^2) time, O(T) time over all calls.
void incorporate_new_phrase(int a, int num_threads) {
  track("incorporate_new_phrase()", ClusterText(Cluster(a)), false);
  
  int s = put_cluster_in_free_slot(a);
  init_slot(s);
  cluster2rep[a] = a;
  rep2cluster[a] = a;
  
  // Compute p1
  p1[s] = (double)phrase_freqs[a] / T;
  
  // Overall all calls: O(T)
  // Compute p2, q2 between a and everything in clusters
  IntIntMap freqs;
  freqs.clear(); // right bigrams
  forvec(_, int, b, right_phrases[a]) {
    b = phrase2rep.GetRoot(b);
    if(!contains(rep2cluster, b)) continue;
    b = rep2cluster[b];
    if(!contains(cluster2slot, b)) continue;
    freqs[b]++;
  }
  forcmap(int, b, int, count, IntIntMap, freqs) {
    curr_minfo += set_p2_q2_from_count(cluster2slot[a], cluster2slot[b], count);
    logs(ClusterText(Cluster(a)) << ' ' << ClusterText(Cluster(b)) << ' ' << count << ' ' << set_p2_q2_from_count(cluster2slot[a], cluster2slot[b], count));
  }
  
  freqs.clear(); // left bigrams
  forvec(_, int, b, left_phrases[a]) {
    b = phrase2rep.GetRoot(b);
    if(!contains(rep2cluster, b)) continue;
    b = rep2cluster[b];
    if(!contains(cluster2slot, b)) continue;
    freqs[b]++;
  }
  forcmap(int, b, int, count, IntIntMap, freqs) {
    curr_minfo += set_p2_q2_from_count(cluster2slot[b], cluster2slot[a], count);
    logs(ClusterText(Cluster(b)) << ' ' << ClusterText(Cluster(a)) << ' ' << count << ' ' << set_p2_q2_from_count(cluster2slot[b], cluster2slot[a], count));
  }
  
  curr_minfo -= q2[s][s]; // q2[s, s] was double-counted
  
  // Update L2: O(C^2)
  track_block("Update L2", "", false) {
    
    the_job.s = s;
    the_job.is_type_a = true;
    // start the jobs
    for (int ii=0; ii<num_threads; ii++) {
      thread_start[ii].notify(); // the thread waits for this lock to begin
    }
    // wait for them to be done
    for (int ii=0; ii<num_threads; ii++) {
      thread_idle[ii].wait();  // the thread releases the lock to finish
    }
  }
}


void update_L2(int thread_id, int num_threads) {
  
  while (true) {
    
    // wait for mutex to unlock to begin the job
    thread_start[thread_id].wait();
    if ( all_done ) break;  // mechanism to close the threads
    int num_clusters = len(slot2cluster);
    
    if (the_job.is_type_a) {
      int s = the_job.s;
      
      for(int t=thread_id; t < num_clusters; t += num_threads) { // L2[s, *], L2[*, s]
        if (slot2cluster[t] == -1) continue;
        if (s == t) continue;
        int S, T;
        if(ORDER_VALID(s, t)) S = s, T = t;
        else                  S = t, T = s;
        L2[S][T] = compute_L2(S, T);
      }
      
      for(int t=thread_id; t < num_clusters; t += num_threads) {
        if (slot2cluster[t] == -1) continue;
        if (t == s) continue;
        FOR_SLOT(u) {
          if(u == s) continue;
          if(!ORDER_VALID(t, u)) continue;
          L2[t][u] += bi_q2(t, s) + bi_q2(u, s) - bi_hyp_q2(make_pair(t, u), s);
        }
      }
      
    } else {       // this is a type B job
      int s = the_job.s;
      int t = the_job.t;
      int u = the_job.u;
      
      for (int v = thread_id; v < num_clusters; v += num_threads) {
        if ( slot2cluster[v] == -1) continue;
        for ( int w = 0; w < num_clusters; w++) {
          if ( slot2cluster[w] == -1) continue;
          if(!ORDER_VALID(v, w)) continue;
          
          if(v == u || w == u)
            L2[v][w] = compute_L2(v, w);
          else
            L2[v][w] = compute_L2_using_old(s, t, u, v, w);
        }
      }
    }
    
    // signal that the thread is done by unlocking the mutex
    thread_idle[thread_id].notify();
  }
}

// O(C^2) time.
// Merge clusters a (in slot s) and b (in slot t) into c (in slot u).
void merge_clusters(int s, int t, int num_threads) {
  assert(ORDER_VALID(s, t));
  int a = slot2cluster[s];
  int b = slot2cluster[t];
  int c = curr_cluster_id++;
  int u = put_cluster_in_free_slot(c);
  
  free_up_slots(s, t);
  
  // Record merge in the cluster tree
  cluster_tree[c] = make_pair(a, b);
  curr_minfo -= L2[s][t];
  
  // Update relationship between clusters and rep phrases
  int A = cluster2rep[a];
  int B = cluster2rep[b];
  phrase2rep.Join(A, B);
  int C = phrase2rep.GetRoot(A); // New rep phrase of cluster c (merged a and b)
  
  track("Merging clusters", ClusterText(Cluster(a)) << " and " << ClusterText(Cluster(b)) << " into " << c << ", lost " << L2[s][t], false);
  
  cluster2rep.erase(a);
  cluster2rep.erase(b);
  rep2cluster.erase(A);
  rep2cluster.erase(B);
  cluster2rep[c] = C;
  rep2cluster[C] = c;
  
  // Compute p1: O(1)
  p1[u] = p1[s] + p1[t];
  
  // Compute p2: O(C)
  p2[u][u] = hyp_p2(make_pair(s, t));
  FOR_SLOT(v) {
    if(v == u) continue;
    p2[u][v] = hyp_p2(make_pair(s, t), v);
    p2[v][u] = hyp_p2(v, make_pair(s, t));
  }
  
  // Compute q2: O(C)
  q2[u][u] = hyp_q2(make_pair(s, t));
  FOR_SLOT(v) {
    if(v == u) continue;
    q2[u][v] = hyp_q2(make_pair(s, t), v);
    q2[v][u] = hyp_q2(v, make_pair(s, t));
  }
  
  // Compute L2: O(C^2)
  track_block("Compute L2", "", false) {
    the_job.s = s;
    the_job.t = t;
    the_job.u = u;
    the_job.is_type_a = false;
    
    // start the jobs
    for (int ii=0; ii<num_threads; ii++) {
      thread_start[ii].notify(); // the thread waits for this lock to begin
    }
    // wait for them to be done
    for (int ii=0; ii<num_threads; ii++) {
      thread_idle[ii].wait();  // the thread releases the lock to finish
    }
  }
}
void merge_clusters(const IntPair &st, int num_threads) { merge_clusters(st.first, st.second, num_threads); }

// MAKE SURE THIS IS NOT DEFINED FOR EFFICIENCY!
//#define PRINT_RANKED

// Merge the optimal pair of clusters that result in the least amount of lost
// mutual information.
// Return the slots.
// O(C^2) time.
IntPair find_opt_clusters_to_merge() {
  track("find_opt_clusters_to_merge()", "", false);
  int best_s = -1, best_t = -1;
  double min_l = 1e30;
  
  // Pick two clusters to merge
  FOR_SLOT(s) {
    FOR_SLOT(t) {
      if(!ORDER_VALID(s, t)) continue;
      // Consider merging clusters in slots s and t.
      double l = L2[s][t];
#ifndef PRINT_RANKED
      logs("If merge clusters " << ClusterText(Slot(s)) << " and " << ClusterText(Slot(t)) << ", lose " << l << ", resulting minfo = " << curr_minfo-l);
#endif
      if(l < min_l) {
        min_l = l;
        best_s = s;
        best_t = t;
      }
    }
  }
  
#ifdef PRINT_RANKED
  vector< pair<double, IntPair> > merges;
  FOR_SLOT(s) {
    FOR_SLOT(t) {
      if(!ORDER_VALID(s, t)) continue;
      merges.push_back(pair<double, IntPair>(L2[s][t], make_pair(s, t)));
    }
  }
  sort(merges.begin(), merges.end());
  for(int i = 0; i < len(merges); i++) {
    const IntPair &st = merges[i].second;
    int s = st.first;
    int t = st.second;
    double l = merges[i].first;
    logs("If merge clusters " << ClusterText(Slot(s)) << " and " << ClusterText(Slot(t)) << ", lose " << l << ", resulting minfo = " << curr_minfo-l);
  }
#endif
  
  return IntPair(best_s, best_t);
}

int phrase2cluster(int a) {
  a = phrase2rep.GetRoot(a);
  assert2(contains(rep2cluster, a), a);
  return rep2cluster[a];
}

double kl_divergence(const IntIntMap &a_count2, int a_count1, const IntPairIntMap &count2,
                   const IntIntMap &count1, int ca, bool right) {
  double kl = 0;
  forcmap(int, cb, int, count, IntIntMap, a_count2) {
    double p = (double)count/a_count1; // P(cb | a)
    IntPair cab = right ? IntPair(ca, cb) : IntPair(cb, ca);
    double q = (double)(count2.find(cab)->second)/count1.find(ca)->second; // P(cb | ca)
    kl += p * log(p/q);
  }
  return kl;
}

// Motivation: each word has it's own identity (characterized by a
// distribution of its context).  The cluster has a distribution over
// contexts.  We can define an assignment of a word to a cluster by comparing
// this similarity.
// For each cluster, compute the cluster distributions.
void compute_cluster_distribs() {
  track("compute_cluster_distribs()", "", true);
  
  IntPairIntMap count2; // (cluster a, cluster b) -> number of times a-b appears
  IntIntMap count1; // cluster a -> number of times a appears
  
  // Compute cluster distributions
  foridx(a, N) {
    int ca = phrase2cluster(a);
    forvec(_, int, b, right_phrases[a]) {
      int cb = phrase2cluster(b);
      count2[IntPair(ca, cb)]++;
      count1[ca]++;
      count1[cb]++;
    }
  }
  
  // For each word (phrase), compute its distribution
  kl_map[0].resize(N);
  kl_map[1].resize(N);
  foridx(a, N) {
    int ca = phrase2cluster(a);
    IntIntMap a_count2;
    int a_count1 = 0;
    //double kl;
    
    // Left distribution
    a_count2.clear(), a_count1 = 0;
    forvec(_, int, b, left_phrases[a]) {
      int cb = phrase2cluster(b);
      a_count2[cb]++;
      a_count1++;
    }
    //kl = kl_map[0][a] = kl_divergence(a_count2, a_count1, count2, count1, ca, false);
    kl_map[0][a] = kl_divergence(a_count2, a_count1, count2, count1, ca, false);
    //logs("Left-KL(" << Phrase(a) << " | " << Cluster(ca) << ") = " << kl);
    
    // Right distribution
    a_count2.clear(), a_count1 = 0;
    forvec(_, int, b, right_phrases[a]) {
      int cb = phrase2cluster(b);
      a_count2[cb]++;
      a_count1++;
    }
    //kl = kl_map[1][a] = kl_divergence(a_count2, a_count1, count2, count1, ca, true);
    kl_map[1][a] = kl_divergence(a_count2, a_count1, count2, count1, ca, true);
    //logs("Right-KL(" << Phrase(a) << " | " << Cluster(ca) << ") = " << kl);
  }
}


void do_clustering(bool chk, int num_threads, int initC) {
  track("do_clustering()", "", true);
  
  compute_L2();
  repcheck(chk);
  
  // start the threads
  { vector<semaphore> i(num_threads); thread_start.swap(i); }
  { vector<semaphore> i(num_threads); thread_idle.swap(i); }
  threads = new thread[num_threads];
  for (int ii=0; ii<num_threads; ii++) {
    thread_start[ii].notify();
    thread_idle[ii].notify();
    threads[ii] = thread(update_L2, ii, num_threads);
  }
  
  curr_cluster_id = N; // New cluster ids will start at N, after all the phrases.
  
  // Stage 1: Maintain initC clusters.  For each of the phrases initC..N-1, make
  // it into a new cluster.  Then merge the optimal pair among the initC+1
  // clusters.
  // O(N*C^2) time.
  track_block("Stage 1", "", false) {
    mem_tracker.report_mem_usage();
    for(int i = initC; i < len(freq_order_phrases); i++) { // Merge phrase new_a
      int new_a = freq_order_phrases[i];
      track("Merging phrase", i << '/' << N << ": " << ClusterText(Cluster(new_a)), true);
      logs("Mutual info: " << curr_minfo);
      incorporate_new_phrase(new_a, num_threads);
      repcheck(chk);
      IntPair cl = find_opt_clusters_to_merge();
      merge_clusters(find_opt_clusters_to_merge(), num_threads);
      repcheck(chk);
    }
  }
  
  compute_cluster_distribs();
  
  stage2_cluster_offset = curr_cluster_id;
  
  // Stage 2: Merge the initC clusters in an hierarchical manner.
  // O(C^3) time.
  track_block("Stage 2", "", false) {
    mem_tracker.report_mem_usage();
    track_foridx(i, initC-1, "Clustering", true) {
      logs("Mutual info of " << len(cluster2slot) << " clusters: " << curr_minfo);
      merge_clusters(find_opt_clusters_to_merge(), num_threads);
      repcheck(chk);
    }
  }
  
  // finish the threads
  all_done = true;
  for (int ii=0; ii<num_threads; ii++) {
    thread_start[ii].notify(); // thread will grab this to start
    threads[ii].join();
  }
  { vector<semaphore> i; thread_start.swap(i); }
  { vector<semaphore> i; thread_idle.swap(i); }
  delete [] threads;
  
  logs("Done: 1 cluster left: mutual info = " << curr_minfo);
  mem_tracker.report_mem_usage();
  //assert(feq(curr_minfo, 0.0));
}

struct StackItem {
  StackItem(int a, int path_i, char ch) : a(a), path_i(path_i), ch(ch) { }
  int a;
  int path_i;
  char ch;
};



Rcpp::DataFrame brown_collocations(int ncollocs) {
  vector< pair<double, IntPair> > collocs;
  FOR_SLOT(s) FOR_SLOT(t) {
    collocs.push_back(pair<double, IntPair>(q2[s][t], make_pair(slot2cluster[s], slot2cluster[t])));
  }
  ncollocs = len(collocs);
  //ncollocs = min(ncollocs, len(collocs));
  partial_sort(collocs.begin(), collocs.begin()+ncollocs, collocs.end(), greater< pair<double, IntPair> >());
  
  std::vector<std::string> term1;
  std::vector<std::string> term2;
  std::vector<double> collocation;
  for(int i = 0; i < ncollocs; i++) {
    const IntPair &ab = collocs[i].second;
    term1.push_back(PhraseText(Phrase(ab.first)));
    term2.push_back(PhraseText(Phrase(ab.second)));
    collocation.push_back(collocs[i].first);
  }
  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("term1") = term1,
    Rcpp::Named("term2") = term2,
    Rcpp::Named("collocation") = collocation,
    Rcpp::Named("stringsAsFactors") = false
  );
  return out;
}


Rcpp::DataFrame brown_clusterassignment() {
  std::vector<std::string> word;
  std::vector<std::string> word_path;
  std::vector<int> freq;
  std::vector<double> kl_left;
  std::vector<double> kl_right;
  // The cluster tree is composed of the top part, which consists
  // of Stage 2 merges, and the bottom part, which consists of stage 1 merges.
  // Print out paths from the root only through the stage 2 merges.
  char path[16384];
  vector<StackItem> stack;
  stack.push_back(StackItem(cluster2slot.begin()->first, 0, '\0'));
  
  while(!stack.empty()) {
    // Take off a stack item (a node in the tree).
    StackItem item = stack.back();
    int a = item.a;
    int path_i = item.path_i;
    if(item.ch)
      path[path_i-1] = item.ch;
    stack.pop_back();
    
    // Look at the node's children (if any).
    IntIntPairMap::const_iterator it = cluster_tree.find(a);
    if(it == cluster_tree.end()) {
      path[path_i] = '\0';
      std::string p(path);
      word.push_back(PhraseText(Phrase(a)));
      word_path.push_back(p);
      freq.push_back(phrase_freqs[a]);
      kl_left.push_back(kl_map[0][a]);
      kl_right.push_back(kl_map[1][a]);
      // if(out_paths) paths_out << path << '\t' << Phrase(a) << '\t' << phrase_freqs[a] << endl;
      // if(out_map) map_out << Phrase(a) << '\t'
      //                     << path << "-L " << kl_map[0][a] << '\t'
      //                     << path << "-R " << kl_map[1][a] << '\t'
      //                     << path << "-freq " << phrase_freqs[a] << endl;
    }
    else {
      const IntPair &children = it->second;
      // Only print out paths through the part of the tree constructed in stage 2.
      bool extend = a >= stage2_cluster_offset;
      int new_path_i = path_i + extend;
      stack.push_back(StackItem(children.second, new_path_i, extend ? '1' : '\0'));
      stack.push_back(StackItem(children.first, new_path_i, extend ? '0' : '\0'));
    }
  }
  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("word") = word,
    Rcpp::Named("path") = word_path,
    Rcpp::Named("freq") = freq,
    Rcpp::Named("kl_left") = kl_left,
    Rcpp::Named("kl_right") = kl_right,
    Rcpp::Named("stringsAsFactors") = false
  );
  return out;
}




// [[Rcpp::export]]
Rcpp::List cluster_brown(std::vector<std::string> x, 
                         std::vector<std::string> vocabulary, 
                  int min_occur = 5, int initC = 100,
                  int ncollocs = 500,
                  int plen = 1,
                  int num_threads = 1,
                  bool chk = false,
                  bool print_stats = false,
                  bool paths2map = false) {
  string featvec_file = "";
  //init_log;
  
  track_mem(db);
  track_mem(phrase_freqs);
  track_mem(left_phrases);
  track_mem(right_phrases);
  track_mem(cluster_tree);
  track_mem(freq_order_phrases);
  track_mem(phrase2rep);
  track_mem(rep2cluster);
  track_mem(cluster2rep);
  track_mem(phrases);
  track_mem(slot2cluster);
  track_mem(cluster2slot);
  track_mem(p1);
  track_mem(p2);
  track_mem(q2);
  track_mem(L2);

  bool use_restrict = vocabulary.size() > 0;
  if(use_restrict){
    for(unsigned int i = 0; i < vocabulary.size(); i++){
      std::string s = vocabulary[i];
      db.lookup(s.c_str(), true, -1); // use only this vocabulary
    }
  }
  IntVec text;
  for(unsigned int i = 0; i < x.size(); i++){
    std::string s = x[i];
    int w = db.lookup(s.c_str(), !use_restrict, -1); // incorporate words either from the vocabulary set or all
    text.push_back(w);
  }
  if(db.size() == 0){
    Rcpp::stop("No data provided");
  }

  //read_restrict_text(restrict_file);
  read_text(paths2map, plen, min_occur, initC, text);
  Rcpp::DataFrame collocations;
  Rcpp::DataFrame clusterassignment;
  
  track("create_initial_clusters()", "", true);
  
  freq_order_phrases.resize(N);
  foridx(a, N) freq_order_phrases[a] = a;
  
  logs("Sorting " << N << " phrases by frequency");
  sort(freq_order_phrases.begin(), freq_order_phrases.end(), phrase_freq_greater);
  
  if(initC > N){
    initC = N;
  }
  create_initial_clusters(initC);
  collocations = brown_collocations(ncollocs);
  do_clustering(chk, num_threads, initC);
  clusterassignment = brown_clusterassignment();
  
  db.clear();                // word database
  phrase_freqs.clear();     // phrase a < N -> number of times a appears in the text
  left_phrases.clear();  // phrase a < N -> list of phrases that appear to left of a in the text
  right_phrases.clear(); // phrase a < N -> list of phrases that appear to right of a in the text
  cluster_tree.clear(); // cluster c -> the 2 sub-clusters that merged to create c
  freq_order_phrases.clear(); // List of phrases in decreasing order of frequency.
  //UnionSet phrase2rep;   // phrase a -> the rep phrase in the same cluster as a
  rep2cluster.clear(); // rep phrase a -> the cluster that contains a
  cluster2rep.clear(); // cluster a -> the rep phrase in cluster a
  phrases.clear(); // length of phrase -> flattened list of words
  slot2cluster.clear(); // slot index -> cluster (-1 if none exists)
  cluster2slot.clear(); // cluster -> slot index
  p1.clear();        // slot s (containing cluster a) -> probability Pr(a)
  p2.clear();     // slots s, t (containing clusters a, b) -> probability Pr(a, b)
  q2.clear();     // slots s, t (contianing clusters a, b) -> contribution to mutual information
  L2.clear();     // slots s, t (containing clusters a, b) -> loss of mutual information if merge a and b
  kl_map[2].clear();
  
  /*
   // Variables used to control the thread pool
   vector<semaphore> thread_idle;
   vector<semaphore> thread_start;
   thread * threads;
   */
  //Compute_L2_Job the_job;
  all_done = false;

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("initC") = initC,
    Rcpp::Named("min_occur") = min_occur,
    Rcpp::Named("clusters") = clusterassignment,
    Rcpp::Named("collocations") = collocations
  );
  out.attr("class") = "brown";
  return out;
}


