/*
 * File:   Deflate.h
 *
 * Small deflate (RFC 1951) and zlib (RFC 1950) compressor for the result files,
 * the output is readable by zlib (ParaView, .NET DeflateStream).
 */

#ifndef DEFLATE_H
#define DEFLATE_H

#include <cstddef>
#include <cstdint>
#include <vector>

// raw deflate of one chunk; the output ends on a byte boundary (a final block or an empty
// stored block), so independently compressed chunks concatenate into one deflate stream
void deflateChunk(const uint8_t* data, size_t n, bool last, std::vector<uint8_t>& out);

// complete zlib stream (header, deflate, adler32)
std::vector<uint8_t> zlibCompress(const uint8_t* data, size_t n);

// one raw deflate stream compressed in parallel chunks
std::vector<uint8_t> deflateParallel(const uint8_t* data, size_t n);

#endif /* DEFLATE_H */
