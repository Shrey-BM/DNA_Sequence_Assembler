# DNA Sequence Assembler

A C++ program that rebuilds a DNA sequence from short, overlapping fragments (reads), built for our **Algorithms and Data Structures Lab (Project 2)**.

A DNA sequencer can't read a whole genome in one go. It gives thousands of short reads in random order, and some of them contain errors. This program finds how the reads overlap and joins them back together into the original sequence.

**Method:** Overlap graph → Dynamic programming overlap scoring → Greedy assembly

---

## Features

- Overlaps are scored with **edit distance (DP)**, so reads with wrong, extra or missing bases can still be joined
- **k-mer hash index** so DP only runs on read pairs that share a short word (about 7× faster on 600 reads)
- **Greedy assembly** that always merges the best-overlapping pair first, using a max-heap
- Removes **duplicate reads** and reads that sit **inside other reads**
- **O(1) loop check** when joining reads into chains
- Built-in **sample data** and a **sequencing simulator** (random genome + errors) that checks the result against the real genome

---

## Data Structures and Algorithms Used

| Data structure / algorithm | Where it is used |
|---|---|
| Graph (adjacency list) | Overlap graph: node = read, edge A → B = end of A overlaps start of B |
| Dynamic programming | `scoreOverlap()`: edit distance table for overlaps and containment |
| Greedy algorithm | `assemble()`: merges the best overlap first |
| Priority queue (max-heap) | Keeps all overlap edges, highest score on top |
| Hash table | `unordered_set` removes duplicates, `unordered_map` stores the k-mer index |
| Linked list | `next[]` / `prev[]` arrays chain the merged reads into contigs |
| Vector / array | Reads, DP table, contigs |

---

## How It Works

1. **Load reads**: typed in, from the sample, or from the simulator
2. **Remove duplicates** with a hash set
3. **Build a k-mer index** (hash table) to find candidate read pairs
4. **Score overlaps with DP** (edit distance) for each candidate pair
5. **Remove reads** that lie inside other reads
6. **Build the overlap graph** (adjacency list)
7. **Greedy merge**: pop the best edge from the max-heap and join the two reads if neither is already joined on that side and it doesn't form a loop
8. **Walk each chain** to build the final contig(s) and compare with the original genome

### DP recurrence

```
D[i][0] = 0                     (overlap can start anywhere in read A)
D[0][j] = j
D[i][j] = min( D[i-1][j-1] + (A[i] != B[j]),
               D[i-1][j] + 1,
               D[i][j-1] + 1 )

An overlap of length j is kept if D[n][j] <= 20% of j
score = j - 2 x errors
```

---

## How to Run

### Compile
```bash
g++ -std=c++17 DNA_Sequence_Assembler.cpp -o DNA_Sequence_Assembler
```

### Run
```bash
# Windows
.\DNA_Sequence_Assembler.exe

# Linux / macOS
./DNA_Sequence_Assembler
```


### Menu
```
1. Enter reads
2. Load sample reads
3. Simulate sequencing (random genome)
4. Show reads
5. Show overlap graph
6. Assemble genome
7. Settings (min overlap, max error %)
8. Clear reads
9. Exit
```

### Quick demo
- Type `2` → `6` to assemble the built-in sample (80 bp genome, rebuilt exactly)
- Type `8` → `3` and enter `5000`, `100`, `600`, `1`, `9` → then `6` for a simulated 5,000 bp genome with 1% errors

---

## Results

| Genome | Reads | Errors | Contigs | Match |
|---|---|---|---|---|
| 80 bp (sample) | 12 × 20 bp | 1 base | 1 | exact |
| 500 bp | 100 × 50 bp | 0% | 1 | exact |
| 1,000 bp | 200 × 60 bp | 2% | 1 | 96.9% |
| 2,000 bp | 300 × 80 bp | 1% | 1 | 98.6% |
| 5,000 bp | 600 × 100 bp | 1% | 1 | 98.9% |
| 10,000 bp | 1,000 × 100 bp | 0% | 1 | exact |

- Every error-free genome was rebuilt exactly
- With 1–2% errors, the contig matches about 97–99% of the real genome
- Simulator results can vary slightly between computers, because `rand()` produces different sequences on different compilers

---

## Time Complexity

| Step | Time |
|---|---|
| Remove duplicates | O(N·L) |
| Build k-mer index | O(N·L) |
| DP overlap scoring | O(C·L²) |
| Greedy merging | O(E log E) |
| Build contigs | O(G) |

N = number of reads, L = read length, C = candidate pairs, E = overlap edges, G = genome length

---

## Limitations

- Repeats longer than a read can confuse the greedy step
- Reads from the opposite DNA strand (reverse complement) are not handled
- Above about 3% read errors, the result can split into several contigs
- No consensus step, so a few read errors can end up in the final contig

---

## Team (Group 7)

| Name | Roll No. |
|---|---|
| Shreyansh Nayak | BM25BTECH11025 |
| Pradnesh Suvvada | BM25BTECH11020 |
