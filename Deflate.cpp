/*
 * File:   Deflate.cpp
 *
 * Greedy LZ77 with hash chains and dynamic Huffman blocks, roughly the speed and the ratio
 * of the fast levels of zlib.
 */

#include "Deflate.h"
#include <algorithm>
#include <cstring>
#include <queue>

using namespace std;

namespace {
	const int WINDOW = 32768;
	const int HASH_BITS = 15;
	const int MAX_CHAIN = 16;
	const int MIN_MATCH = 4;
	const int MAX_MATCH = 258;
	const size_t BLOCK_SYMBOLS = 1 << 16;

	const int lenBase[29] = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258 };
	const int lenExtra[29] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0 };
	const int distBase[30] = { 1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577 };
	const int distExtra[30] = { 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13 };
	// order of the code length code lengths in the block header
	const int clOrder[19] = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };

	struct Tables {
		uint8_t lenCode[MAX_MATCH + 1];
		uint8_t distCode[WINDOW + 1];
		Tables() {
			for (int c = 0; c < 29; c++) {
				int end = c < 28 ? lenBase[c + 1] : MAX_MATCH + 1;
				for (int l = lenBase[c]; l < end; l++) {
					lenCode[l] = (uint8_t)c;
				}
			}
			for (int c = 0; c < 30; c++) {
				int end = c < 29 ? distBase[c + 1] : WINDOW + 1;
				for (int d = distBase[c]; d < end; d++) {
					distCode[d] = (uint8_t)c;
				}
			}
		}
	};

	const Tables& tables() {
		static Tables t;
		return t;
	}

	class BitWriter {
	public:
		vector<uint8_t>& out;
		uint64_t acc = 0;
		int nb = 0;
		BitWriter(vector<uint8_t>& out) : out(out) {}
		void put(uint32_t v, int n) {
			acc |= (uint64_t)v << nb;
			nb += n;
			while (nb >= 8) {
				out.push_back((uint8_t)acc);
				acc >>= 8;
				nb -= 8;
			}
		}
		void align() {
			if (nb > 0) {
				out.push_back((uint8_t)acc);
			}
			acc = 0;
			nb = 0;
		}
	};

	// literal (dist == 0) or match
	struct Symbol {
		uint16_t litLen;
		uint16_t dist;
	};

	// Huffman code lengths limited to maxBits, the limit is reached by flattening the frequencies
	void huffmanLengths(vector<uint32_t> freq, int maxBits, vector<uint8_t>& len) {
		int n = (int)freq.size();
		len.assign(n, 0);
		while (true) {
			vector<int> sym;
			for (int s = 0; s < n; s++) {
				if (freq[s] > 0) {
					sym.push_back(s);
				}
			}
			int m = (int)sym.size();
			if (m == 1) {
				len[sym[0]] = 1;
				return;
			}
			// leaves 0..m-1, inner nodes m..2m-2, a parent has always a higher index than its children
			vector<int> parent(2 * m - 1, -1);
			typedef pair<uint64_t, int> Node;
			priority_queue<Node, vector<Node>, greater<Node> > q;
			for (int i = 0; i < m; i++) {
				q.push(Node(freq[sym[i]], i));
			}
			int next = m;
			while (q.size() > 1) {
				Node a = q.top();
				q.pop();
				Node b = q.top();
				q.pop();
				parent[a.second] = next;
				parent[b.second] = next;
				q.push(Node(a.first + b.first, next++));
			}
			vector<int> depth(2 * m - 1, 0);
			int maxDepth = 0;
			for (int i = 2 * m - 3; i >= 0; i--) {
				depth[i] = depth[parent[i]] + 1;
				maxDepth = max(maxDepth, depth[i]);
			}
			if (maxDepth <= maxBits) {
				for (int i = 0; i < m; i++) {
					len[sym[i]] = (uint8_t)depth[i];
				}
				return;
			}
			for (int s = 0; s < n; s++) {
				if (freq[s] > 0) {
					freq[s] = (freq[s] >> 1) | 1;
				}
			}
		}
	}

	// canonical codes, bit reversed as deflate writes them from the least significant bit
	void canonicalCodes(const vector<uint8_t>& len, vector<uint16_t>& code) {
		int count[16] = { 0 }, next[16] = { 0 };
		for (uint8_t l : len) {
			count[l]++;
		}
		count[0] = 0;
		int c = 0;
		for (int b = 1; b < 16; b++) {
			c = (c + count[b - 1]) << 1;
			next[b] = c;
		}
		code.assign(len.size(), 0);
		for (size_t s = 0; s < len.size(); s++) {
			int l = len[s];
			if (l > 0) {
				int v = next[l]++, r = 0;
				for (int b = 0; b < l; b++) {
					r = (r << 1) | ((v >> b) & 1);
				}
				code[s] = (uint16_t)r;
			}
		}
	}

	// at least two used symbols keep the codes complete for every inflater
	void ensureTwo(vector<uint32_t>& freq) {
		int used = 0;
		for (uint32_t f : freq) {
			used += f > 0;
		}
		for (size_t s = 0; used < 2 && s < freq.size(); s++) {
			if (freq[s] == 0) {
				freq[s] = 1;
				used++;
			}
		}
	}

	void writeBlock(BitWriter& bw, const vector<Symbol>& syms, bool final) {
		const Tables& t = tables();
		vector<uint32_t> litFreq(286, 0), distFreq(30, 0);
		for (const Symbol& s : syms) {
			if (s.dist == 0) {
				litFreq[s.litLen]++;
			}
			else {
				litFreq[257 + t.lenCode[s.litLen]]++;
				distFreq[t.distCode[s.dist]]++;
			}
		}
		litFreq[256] = 1;
		ensureTwo(litFreq);
		ensureTwo(distFreq);
		vector<uint8_t> litLen, distLen;
		huffmanLengths(litFreq, 15, litLen);
		huffmanLengths(distFreq, 15, distLen);
		int hlit = 286, hdist = 30;
		while (hlit > 257 && litLen[hlit - 1] == 0) {
			hlit--;
		}
		while (hdist > 1 && distLen[hdist - 1] == 0) {
			hdist--;
		}

		// run length coding of the code lengths: (symbol, extra bits value)
		vector<uint8_t> all(litLen.begin(), litLen.begin() + hlit);
		all.insert(all.end(), distLen.begin(), distLen.begin() + hdist);
		vector<pair<int, int> > rle;
		for (size_t i = 0; i < all.size();) {
			int cur = all[i];
			size_t r = 1;
			while (i + r < all.size() && all[i + r] == cur) {
				r++;
			}
			i += r;
			if (cur == 0) {
				while (r >= 11) {
					size_t k = min<size_t>(r, 138);
					rle.push_back(make_pair(18, (int)k - 11));
					r -= k;
				}
				if (r >= 3) {
					rle.push_back(make_pair(17, (int)r - 3));
					r = 0;
				}
			}
			else {
				rle.push_back(make_pair(cur, 0));
				r--;
				while (r >= 3) {
					size_t k = min<size_t>(r, 6);
					rle.push_back(make_pair(16, (int)k - 3));
					r -= k;
				}
			}
			for (; r > 0; r--) {
				rle.push_back(make_pair(cur, 0));
			}
		}
		vector<uint32_t> clFreq(19, 0);
		for (auto& p : rle) {
			clFreq[p.first]++;
		}
		ensureTwo(clFreq);
		vector<uint8_t> clLen;
		huffmanLengths(clFreq, 7, clLen);
		int hclen = 19;
		while (hclen > 4 && clLen[clOrder[hclen - 1]] == 0) {
			hclen--;
		}
		vector<uint16_t> litCode, distCode, clCode;
		canonicalCodes(litLen, litCode);
		canonicalCodes(distLen, distCode);
		canonicalCodes(clLen, clCode);

		bw.put(final ? 1 : 0, 1);
		bw.put(2, 2);
		bw.put(hlit - 257, 5);
		bw.put(hdist - 1, 5);
		bw.put(hclen - 4, 4);
		for (int i = 0; i < hclen; i++) {
			bw.put(clLen[clOrder[i]], 3);
		}
		for (auto& p : rle) {
			bw.put(clCode[p.first], clLen[p.first]);
			if (p.first == 16) {
				bw.put(p.second, 2);
			}
			else if (p.first == 17) {
				bw.put(p.second, 3);
			}
			else if (p.first == 18) {
				bw.put(p.second, 7);
			}
		}
		for (const Symbol& s : syms) {
			if (s.dist == 0) {
				bw.put(litCode[s.litLen], litLen[s.litLen]);
			}
			else {
				int lc = t.lenCode[s.litLen], dc = t.distCode[s.dist];
				bw.put(litCode[257 + lc], litLen[257 + lc]);
				bw.put(s.litLen - lenBase[lc], lenExtra[lc]);
				bw.put(distCode[dc], distLen[dc]);
				bw.put(s.dist - distBase[dc], distExtra[dc]);
			}
		}
		bw.put(litCode[256], litLen[256]);
	}

	inline uint32_t hash4(const uint8_t* p) {
		uint32_t v;
		memcpy(&v, p, 4);
		return (v * 2654435761u) >> (32 - HASH_BITS);
	}

	uint32_t adler32(const uint8_t* data, size_t n) {
		uint32_t a = 1, b = 0;
		while (n > 0) {
			size_t k = min<size_t>(n, 5552);
			n -= k;
			for (size_t i = 0; i < k; i++) {
				a += data[i];
				b += a;
			}
			data += k;
			a %= 65521;
			b %= 65521;
		}
		return (b << 16) | a;
	}
}

void deflateChunk(const uint8_t* data, size_t n, bool last, vector<uint8_t>& out) {
	BitWriter bw(out);
	vector<int32_t> head(1 << HASH_BITS, -1), prev(WINDOW, -1);
	vector<Symbol> syms;
	syms.reserve(BLOCK_SYMBOLS);
	auto insert = [&](size_t p) {
		if (p + 4 <= n) {
			uint32_t h = hash4(data + p);
			prev[p & (WINDOW - 1)] = head[h];
			head[h] = (int32_t)p;
		}
	};
	size_t i = 0;
	while (i < n) {
		int bestLen = 0, bestDist = 0;
		if (i + MIN_MATCH <= n) {
			int maxLen = (int)min<size_t>(MAX_MATCH, n - i);
			int32_t cand = head[hash4(data + i)];
			for (int chain = 0; chain < MAX_CHAIN && cand >= 0 && i - cand <= WINDOW; chain++) {
				const uint8_t* a = data + i;
				const uint8_t* b = data + cand;
				if (b[bestLen] == a[bestLen]) {
					int l = 0;
					while (l < maxLen && a[l] == b[l]) {
						l++;
					}
					if (l > bestLen) {
						bestLen = l;
						bestDist = (int)(i - cand);
						if (l == maxLen) {
							break;
						}
					}
				}
				int32_t nx = prev[cand & (WINDOW - 1)];
				if (nx >= cand) {
					break;
				}
				cand = nx;
			}
		}
		if (bestLen >= MIN_MATCH) {
			syms.push_back({ (uint16_t)bestLen, (uint16_t)bestDist });
			for (int k = 0; k < bestLen; k++) {
				insert(i + k);
			}
			i += bestLen;
		}
		else {
			syms.push_back({ data[i], 0 });
			insert(i);
			i++;
		}
		if (syms.size() >= BLOCK_SYMBOLS) {
			writeBlock(bw, syms, last && i == n);
			syms.clear();
			if (last && i == n) {
				bw.align();
				return;
			}
		}
	}
	if (last) {
		writeBlock(bw, syms, true);
	}
	else {
		if (!syms.empty()) {
			writeBlock(bw, syms, false);
		}
		// empty stored block aligns the output to a byte (sync flush)
		bw.put(0, 3);
		bw.align();
		const uint8_t sync[4] = { 0, 0, 0xFF, 0xFF };
		out.insert(out.end(), sync, sync + 4);
	}
	bw.align();
}

vector<uint8_t> zlibCompress(const uint8_t* data, size_t n) {
	vector<uint8_t> out = { 0x78, 0x01 };
	deflateChunk(data, n, true, out);
	uint32_t a = adler32(data, n);
	for (int s = 24; s >= 0; s -= 8) {
		out.push_back((uint8_t)(a >> s));
	}
	return out;
}

vector<uint8_t> deflateParallel(const uint8_t* data, size_t n) {
	const size_t chunk = 1 << 20;
	long long nChunk = max<long long>(1, (long long)((n + chunk - 1) / chunk));
	vector<vector<uint8_t> > parts(nChunk);
#pragma omp parallel for schedule(dynamic)
	for (long long c = 0; c < nChunk; c++) {
		size_t start = c * chunk;
		deflateChunk(data + start, min(chunk, n - start), c == nChunk - 1, parts[c]);
	}
	vector<uint8_t> out;
	for (auto& p : parts) {
		out.insert(out.end(), p.begin(), p.end());
	}
	return out;
}
