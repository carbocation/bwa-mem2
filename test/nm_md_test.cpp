// Synthetic unit tests for reference ambiguity in NM/MD. No data files needed.
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include "bwa.h"

static int nt(char c)
{
    const char *bases = "ACGT";
    const char *p = strchr(bases, toupper((unsigned char)c));
    return p? p - bases : 4;
}

static void check(const std::string& reference, const std::string& packed,
                  int beg, int end, bool reverse, bool query_n)
{
    bntseq_t bns = {};
    bns.l_pac = reference.size();
    std::vector<bntamb1_t> holes;
    std::vector<uint8_t> pac((packed.size() + 3) / 4, 0);
    for (size_t i = 0; i < reference.size(); ++i) {
        assert(nt(packed[i]) < 4);
        pac[i >> 2] |= nt(packed[i]) << ((~i & 3) << 1);
        if (nt(reference[i]) == 4) {
            if (!holes.empty() && holes.back().offset + holes.back().len == (int64_t)i
                && holes.back().amb == reference[i]) ++holes.back().len;
            else {
                bntamb1_t hole = {};
                hole.offset = i; hole.len = 1; hole.amb = reference[i];
                holes.push_back(hole);
            }
        }
    }
    bns.n_holes = holes.size();
    bns.ambs = holes.empty()? NULL : &holes[0];
    std::string read = packed.substr(beg, end - beg);
    if (query_n) read[read.size()/2] = 'N';
    std::string expected_md;
    int expected_nm = 0, matches = 0;
    for (int i = 0; i < (int)read.size(); ++i) {
        char r = toupper((unsigned char)reference[beg+i]);
        if (nt(r) < 4 && nt(read[i]) < 4 && r == read[i]) ++matches;
        else {
            expected_md += std::to_string(matches) + r;
            ++expected_nm; matches = 0;
        }
    }
    expected_md += std::to_string(matches);
    std::vector<uint8_t> query(read.size());
    for (size_t i = 0; i < read.size(); ++i) {
        int base = nt(read[reverse? read.size()-1-i : i]);
        query[i] = reverse && base < 4? 3-base : base;
    }
    std::vector<uint8_t> original_query = query;
    int8_t mat[25];
    bwa_fill_scmat(1, 4, mat);
    int expected_score = 0;
    for (int i = 0; i < (int)read.size(); ++i)
        expected_score += mat[nt(packed[beg+i])*5 + nt(read[i])];
    int64_t rb = reverse? 2*bns.l_pac-end : beg;
    int64_t re = reverse? 2*bns.l_pac-beg : end;
    // Exercise both the no-gap shortcut and banded DP, plus the wrapper API.
    for (int width = 0; width <= 10; width += 10) {
        for (int wrapper = 0; wrapper < 2; ++wrapper) {
            int score = 0, n_cigar = 0, nm = -1;
            uint32_t *cigar = wrapper?
                bwa_gen_cigar(mat, 6, 1, width, &bns, &pac[0], query.size(), &query[0],
                              rb, re, &score, &n_cigar, &nm) :
                bwa_gen_cigar2(mat, 6, 1, 6, 1, width, &bns, &pac[0], query.size(), &query[0],
                               rb, re, &score, &n_cigar, &nm);
            assert(n_cigar == 1 && cigar[0] == (query.size()<<4));
            assert(nm == expected_nm);
            assert(expected_md == (char*)(cigar + n_cigar));
            assert(score == expected_score && query == original_query);
            free(cigar);
        }
        int score = 0;
        assert(bwa_gen_cigar2(mat, 6, 1, 6, 1, width, &bns, &pac[0], query.size(), &query[0],
                              rb, re, &score, NULL, NULL) == NULL);
        assert(score == expected_score && query == original_query);
    }
}

int main()
{
    bwa_verbose = 0;
    const std::string packed = "ACGTTGCACTGGATCCGATGCTAGCATGACCTAGGTACGTCAGATCGTAGCTAGCATGCTA";
    std::vector<std::string> references(1, packed);
    std::string r = packed;
    r[9] = r[10] = r[29] = r[30] = 'N';
    references.push_back(r); // holes immediately outside/at both interval ends
    r = packed;
    r.replace(8, 6, 6, 'N');
    r.replace(27, 8, 8, 'n');
    references.push_back(r); // partially intersecting runs, including lowercase
    r = packed;
    r.replace(10, 16, "NRYSWKMBDHVXnrys");
    references.push_back(r); // adjacent distinct ambiguity annotations
    references.push_back(std::string(packed.size(), 'N')); // one spanning hole
    r = packed;
    r[0] = r[packed.size()-1] = 'N';
    references.push_back(r); // holes exist but none overlaps this interval
    for (size_t i = 0; i < references.size(); ++i)
        for (int reverse = 0; reverse < 2; ++reverse)
            for (int query_n = 0; query_n < 2; ++query_n)
                check(references[i], packed, 10, 30, reverse, query_n);
    puts("NM/MD unit tests passed (96 tag calculations plus 48 score-only calls)");
    return 0;
}
