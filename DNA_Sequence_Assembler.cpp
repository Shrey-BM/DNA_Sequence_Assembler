/*
   DNA Sequence Assembler - DSA Project 2
   Group 7: Shreyansh, Pradnesh

   The program takes short DNA reads and tries to join them together using
   the overlaps between reads. It can use entered reads, sample data, or
   generate its own reads for testing.

   We use an overlap graph to find possible joins, dynamic programming
   to score overlaps with errors, and a greedy approach to build sequences.
*/
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
using namespace std;

// Settings can be changed from the menu.
int MIN_OVERLAP = 8;          // shortest overlap that is accepted
int MAX_ERROR_PERCENT = 20;   // max errors allowed in an overlap (% of its length)

// Input and working data
vector<string> reads;         // reads exactly as entered
vector<string> cleaned;       // reads after removing duplicates and contained reads
string originalGenome = "";   // known only for sample / simulated data (to check the result)

// Overlap graph
struct Edge {
    int from, to;   // read 'from' is followed by read 'to'
    int len;        // overlap length (how many bases at the start of 'to' are shared)
    int errors;     // mismatches + insertions + deletions inside the overlap
    int score;      // len - 2 * errors   (higher = better overlap)

    bool operator<(const Edge& o) const {    // used to order edges in the max-heap
        if (score != o.score) return score < o.score;
        if (len != o.len) return len < o.len;
        if (from != o.from) return from > o.from;
        return to > o.to;
    }
};

vector<vector<Edge>> graph;   // graph[i] stores the overlaps that start from read i
int edgeCount = 0;

// Helper functions
bool validDNA(string& s) {    // converts to capitals, accepts only A, C, G, T
    for (char& c : s) {
        c = toupper(c);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return !s.empty();
}

// Remove duplicate reads before doing the overlap checks.
vector<string> uniqueReads() {
    // The set makes duplicate checking quick.
    unordered_set<string> seen;
    vector<string> unique;
    for (string& r : reads)
        if (seen.insert(r).second) unique.push_back(r);
    return unique;
}

// Find and score the overlap between two reads using DP.
// We compare the end of a with the start of b, allowing a few errors.
// D[i][j] stores the minimum edits needed to match j bases of b with a part of a
// ending at i. The match can start anywhere in a.
// We also use the table to check whether b is contained in a.
Edge scoreOverlap(const vector<string>& r, int ai, int bi, bool& bInsideA) {
    const string& a = r[ai];
    const string& b = r[bi];
    int n = a.size(), m = b.size();

    vector<vector<int>> D(n + 1, vector<int>(m + 1));
    for (int i = 0; i <= n; i++) D[i][0] = 0;     // the overlap can start at any point in a
    for (int j = 1; j <= m; j++) D[0][j] = j;     // bases of b that are not matched cost one

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            D[i][j] = min({ D[i - 1][j - 1] + cost,   // match or mismatch
                            D[i - 1][j] + 1,          // deletion
                            D[i][j - 1] + 1 });       // insertion
        }
    }

    // Check if b is mostly contained somewhere inside a.
    // We use a slightly relaxed limit here so small errors do not remove a read.
    int bestInside = m;
    for (int i = 0; i <= n; i++) bestInside = min(bestInside, D[i][m]);
    bInsideA = (bestInside * 200 <= m * MAX_ERROR_PERCENT);

    // For an actual overlap, it has to reach the end of a, so use the last row.
    Edge best = {ai, bi, 0, 0, 0};
    for (int j = MIN_OVERLAP; j < m && j <= n; j++) {   // b should contribute at least some new bases
        int e = D[n][j];
        if (e * 100 <= j * MAX_ERROR_PERCENT) {         // within the allowed error limit
            int score = j - 2 * e;
            if (score > best.score) best = {ai, bi, j, e, score};
        }
    }
    return best;   // A length of zero means no overlap was found.
}

// Build the graph using the overlaps we find.
void buildGraph() {
    vector<string> u = uniqueReads();
    int n = u.size();

    // Index reads by short k-mers to avoid checking every pair.
    // Reads sharing a k-mer are treated as possible matches, so DP is only run on
    // these candidate pairs instead of all pairs.
    int K = min(MIN_OVERLAP, 12);
    unordered_map<string, vector<int>> kmerIndex;
    for (int i = 0; i < n; i++) {
        unordered_set<string> own;                     // to avoid adding the same k-mer twice for one read
        for (int p = 0; p + K <= (int)u[i].size(); p++) {
            string kmer = u[i].substr(p, K);
            if (own.insert(kmer).second) kmerIndex[kmer].push_back(i);
        }
    }
    vector<vector<char>> candidate(n, vector<char>(n, 0));
    for (auto& entry : kmerIndex)
        for (int x : entry.second)
            for (int y : entry.second)
                if (x != y) candidate[x][y] = 1;

    // Check candidate pairs and keep track of overlaps and contained reads.
    vector<Edge> found;
    vector<bool> removed(n, false);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!candidate[i][j]) continue;
            bool jInsideI;
            Edge e = scoreOverlap(u, i, j, jInsideI);
            // If one read is contained in another, we can skip the shorter one.
            if (jInsideI && (u[i].size() > u[j].size() || (u[i].size() == u[j].size() && i < j)))
                removed[j] = true;
            if (e.len > 0) found.push_back(e);
        }
    }

    // Keep the reads we need and renumber them.
    vector<int> newId(n, -1);
    cleaned.clear();
    for (int i = 0; i < n; i++) {
        if (removed[i]) continue;
        newId[i] = cleaned.size();
        cleaned.push_back(u[i]);
    }

    // Store the remaining overlaps in an adjacency list.
    graph.assign(cleaned.size(), vector<Edge>());
    edgeCount = 0;
    for (Edge e : found) {
        if (removed[e.from] || removed[e.to]) continue;
        e.from = newId[e.from];
        e.to = newId[e.to];
        graph[e.from].push_back(e);
        edgeCount++;
    }
}

// Edit distance for two complete sequences. Only two rows are stored to save space.
// This is used to compare the result with the original sequence when available.
int editDistance(const string& a, const string& b) {
    vector<int> prev(b.size() + 1), cur(b.size() + 1);
    for (int j = 0; j <= (int)b.size(); j++) prev[j] = j;
    for (int i = 1; i <= (int)a.size(); i++) {
        cur[0] = i;
        for (int j = 1; j <= (int)b.size(); j++) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            cur[j] = min({ prev[j - 1] + cost, prev[j] + 1, cur[j - 1] + 1 });
        }
        swap(prev, cur);
    }
    return prev[b.size()];
}

// Join reads using the greedy approach.
void assemble() {
    if (reads.empty()) { cout << "  No reads loaded\n"; return; }
    buildGraph();
    int n = cleaned.size();

    // Put all overlaps in a max-heap so the best one is considered first.
    priority_queue<Edge> heap;
    for (int i = 0; i < n; i++)
        for (Edge& e : graph[i]) heap.push(e);

    // These arrays keep track of which reads are connected in each chain.
    vector<int> next(n, -1), prev(n, -1), overlapLen(n, 0), otherEnd(n);
    for (int i = 0; i < n; i++) otherEnd[i] = i;   // initially, each read is a separate chain

    cout << "\n  Reads: " << reads.size() << " entered, " << n
         << " after removing duplicates / contained reads\n";
    cout << "  Overlap graph: " << n << " nodes, " << edgeCount << " edges\n\n";
    cout << "  Greedy merges (best overlap first):\n";

    int merges = 0;
    while (!heap.empty() && merges < n - 1) {
        Edge e = heap.top();
        heap.pop();
        if (next[e.from] != -1 || prev[e.to] != -1) continue;   // each read can have only one next and one previous read
        if (otherEnd[e.from] == e.to) continue;                 // avoid forming a cycle

        next[e.from] = e.to;
        prev[e.to] = e.from;
        overlapLen[e.from] = e.len;

        int head = otherEnd[e.from];    // start of the left chain
        int tail = otherEnd[e.to];      // end of the right chain
        otherEnd[head] = tail;
        otherEnd[tail] = head;
        merges++;

        if (merges <= 15)
            cout << "   " << merges << ". R" << e.from + 1 << " -> R" << e.to + 1
                 << "   overlap " << e.len << ", errors " << e.errors << ", score " << e.score << "\n";
    }
    if (merges > 15) cout << "   ... " << merges - 15 << " more merges\n";

    // Follow each chain and build the final contigs.
    vector<string> contigs;
    int leftover = 0;                                // reads that could not be joined
    for (int i = 0; i < n; i++) {
        if (prev[i] != -1) continue;                 // this read is already part of another chain
        if (next[i] == -1 && n > 1) { leftover++; continue; }
        string contig = cleaned[i];
        for (int cur = i; next[cur] != -1; cur = next[cur])
            contig += cleaned[next[cur]].substr(overlapLen[cur]); // skip the bases already covered by the overlap
        contigs.push_back(contig);
    }
    if (contigs.empty()) { cout << "\n  No overlaps found - reads could not be joined\n"; return; }
    sort(contigs.begin(), contigs.end(),
         [](const string& a, const string& b) { return a.size() > b.size(); });

    cout << "\n  Result: " << contigs.size() << " contig(s)";
    if (leftover > 0) cout << ", " << leftover << " leftover read(s) that matched nothing";
    cout << "\n";
    for (int i = 0; i < (int)contigs.size(); i++) {
        cout << "   Contig " << i + 1 << " (" << contigs[i].size() << " bp): ";
        if (contigs[i].size() <= 100) cout << contigs[i] << "\n";
        else cout << contigs[i].substr(0, 60) << "...\n";
    }

    // If the original sequence is known, compare it with our longest contig.
    if (!originalGenome.empty()) {
        int d = editDistance(contigs[0], originalGenome);
        double identity = 100.0 * (1.0 - (double)d / max(contigs[0].size(), originalGenome.size()));
        cout << "\n  Original genome: " << originalGenome.size() << " bp\n";
        cout << "  Longest contig vs original: edit distance " << d
             << ", identity " << identity << "%\n";
        if (d == 0) cout << "  >> Genome rebuilt EXACTLY <<\n";
    }
}

// Display the overlap graph.
void showGraph() {
    if (reads.empty()) { cout << "  No reads loaded\n"; return; }
    buildGraph();
    cout << "\n  Nodes (reads after cleaning):\n";
    for (int i = 0; i < (int)cleaned.size() && i < 30; i++)
        cout << "   R" << i + 1 << ": " << cleaned[i] << "\n";

    cout << "\n  Edges (adjacency list):\n";
    int shown = 0;
    for (int i = 0; i < (int)graph.size(); i++) {
        for (Edge& e : graph[i]) {
            if (shown++ >= 40) break;
            cout << "   R" << e.from + 1 << " -> R" << e.to + 1 << "   overlap " << e.len
                 << ", errors " << e.errors << ", score " << e.score << "\n";
        }
    }
    cout << "  Total: " << cleaned.size() << " nodes, " << edgeCount << " edges\n";
}

// Input options
void enterReads() {
    int k;
    cout << "How many reads: "; cin >> k;
    for (int i = 0; i < k; i++) {
        string r;
        cout << "Read " << i + 1 << ": "; cin >> r;
        if (validDNA(r)) reads.push_back(r);
        else { cout << "  Only A, C, G, T allowed - read skipped\n"; }
    }
    originalGenome = "";
    cout << "  Total reads: " << reads.size() << "\n";
}

void loadSample() {
    // Sample reads from an 80 bp sequence, with overlaps and shuffled order.
    // The sample also has an error, a duplicate, and a read contained in another.
    originalGenome = "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTT";
    reads = { "AATTCTCATTGTGTTTCGGA", "GTTATTATTTGT",         "GTAGAGGTCGAGGTTATTAT",
              "ATGCTGGGTAGAGGTCGAGG", "ATTGTGTTTCGGAACTTGCG", "TCGCGGTTATTATTTGTTAC",
              "TGTTTCGGAACTTGCGTTTT", "CAAAGATATGCTGGGTAGAG", "AGACTTTCAAAGATATGCTG",
              "TGTTACCAATTCTCATTGTG", "ATGCTGGGTAGAGGTCGAGG", "TATTATTTGTTACCAATTCT" };
    cout << "  Loaded " << reads.size() << " sample reads (80 bp genome)\n";
}

void simulate() {
    int G, L, N, errPct, seed;
    cout << "Genome length: ";             cin >> G;
    cout << "Read length: ";               cin >> L;
    cout << "Number of reads: ";           cin >> N;
    cout << "Error rate per base (%): ";   cin >> errPct;
    cout << "Random seed: ";               cin >> seed;
    if (G < 10 || L < 5 || L > G || N < 2 || errPct < 0) { cout << "  Invalid values\n"; return; }

    srand(seed);
    const char bases[] = "ACGT";
    originalGenome = "";
    for (int i = 0; i < G; i++) originalGenome += bases[rand() % 4];

    reads.clear();
    for (int k = 0; k < N; k++) {
        int start = (k == 0) ? 0 : (k == 1) ? G - L : rand() % (G - L + 1);   // make the first two reads cover the beginning and end of the genome
        string exact = originalGenome.substr(start, L), r = "";
        for (char c : exact) {                       // randomly introduce sequencing errors
            if (rand() % 1000 < errPct * 10) {
                int type = rand() % 3;
                if (type == 0) r += bases[rand() % 4];              // substitution
                else if (type == 1) { r += c; r += bases[rand() % 4]; }  // insertion
                // deletion (the base is skipped)
            } else r += c;
        }
        reads.push_back(r);
    }
    cout << "  Simulated " << N << " reads of ~" << L << " bp from a " << G
         << " bp genome (coverage " << (double)N * L / G << "x)\n";
}

void showReads() {
    if (reads.empty()) { cout << "  No reads loaded\n"; return; }
    for (int i = 0; i < (int)reads.size() && i < 40; i++)
        cout << "   " << i + 1 << ". " << reads[i] << "  (" << reads[i].size() << " bp)\n";
    if (reads.size() > 40) cout << "   ... " << reads.size() - 40 << " more\n";
}

// Main menu
int main() {
    cout.setf(ios::fixed);
    cout.precision(2);
    int choice;
    bool running = true;

    while (running) {
        cout << "\nDNA SEQUENCE ASSEMBLER\n"
             << "1. Enter reads\n" 
             << "2. Load sample reads\n"
             << "3. Simulate sequencing (random genome)\n"
             << "4. Show reads\n"
             << "5. Show overlap graph\n"
             << "6. Assemble genome\n"
             << "7. Settings (min overlap = " << MIN_OVERLAP << ", max error = " << MAX_ERROR_PERCENT << "%)\n"
             << "8. Clear reads\n"
             << "9. Exit\n"
             << "Choice: ";
        if (!(cin >> choice)) break;

        switch (choice) {
        case 1: enterReads();  break;
        case 2: loadSample();  break;
        case 3: simulate();    break;
        case 4: showReads();   break;
        case 5: showGraph();   break;
        case 6: assemble();    break;
        case 7:
            cout << "Min overlap: ";      cin >> MIN_OVERLAP;
            cout << "Max error (%): ";    cin >> MAX_ERROR_PERCENT;
            if (MIN_OVERLAP < 1) MIN_OVERLAP = 1;
            if (MAX_ERROR_PERCENT < 0) MAX_ERROR_PERCENT = 0;
            break;
        case 8:
            reads.clear();
            originalGenome = "";
            cout << "  Cleared\n";
            break;
        case 9: running = false; break;
        default: cout << "  Invalid choice\n";
        }
    }
    cout << "Goodbye!\n";
    return 0;
}
