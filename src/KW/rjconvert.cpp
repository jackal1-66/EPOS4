//
//  This file will hopefullt be part of EPOS4
//  Copyright (C) 2026 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

// Build-time converter for the rj cross-section tables: for every rj*.i in
// the input directory, parse the trailing 360800 array values (fhss, fhgg,
// fhqg, fhgq, fhqq - fixed common-block dimensions) with strtof - the same
// libc conversion both gfortran's list-directed READ and the fastreadf
// runtime reader use, so values are bit-identical - and write them as a
// magic-tagged raw-float sidecar rj*.ib. The ASCII header (whose token
// count may vary) stays in the .i file, which the Fortran loader keeps
// reading; the sidecar holds only the arrays.
//
// usage: rjconvert <dir-with-rj*.i> <output-dir>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <string>
#include <vector>
#include <fstream>

static const char rjb_magic[16] = "EPOS4RJTAB-BIN1"; // 15 chars + NUL
static const size_t NARR = 360800;                   // 4*88000 + 8800

static bool convert(const std::string &in, const std::string &out)
{
    std::string buf;
    {
        std::ifstream fin(in, std::ios::binary);
        if (!fin)
            return false;
        fin.seekg(0, std::ios::end);
        std::streamoff len = fin.tellg();
        if (len <= 0)
            return false;
        buf.resize((size_t)len);
        fin.seekg(0, std::ios::beg);
        fin.read(&buf[0], len);
    }
    // tokenize: arrays are the last NARR whitespace-separated tokens
    std::vector<const char *> tok;
    tok.reserve(400000);
    const char *p = buf.data(), *end = p + buf.size();
    while (p < end)
    {
        while (p < end && (*p == ' ' || *p == '\n' || *p == '\t' || *p == '\r' || *p == '\v' || *p == '\f'))
            ++p;
        if (p >= end)
            break;
        tok.push_back(p);
        while (p < end && !(*p == ' ' || *p == '\n' || *p == '\t' || *p == '\r' || *p == '\v' || *p == '\f'))
            ++p;
    }
    if (tok.size() < NARR)
    {
        std::fprintf(stderr, "rjconvert: %s: only %zu tokens\n", in.c_str(), tok.size());
        return false;
    }

    std::vector<float> a(NARR);
    const size_t first = tok.size() - NARR;
    for (size_t i = 0; i < NARR; i++)
    {
        char *q;
        a[i] = strtof(tok[first + i], &q);
        if (q == tok[first + i])
        {
            std::fprintf(stderr, "rjconvert: %s: parse failure at token %zu\n", in.c_str(), first + i);
            return false;
        }
    }

    std::FILE *f = std::fopen(out.c_str(), "wb");
    if (!f)
        return false;
    int n = (int)NARR;
    bool ok = std::fwrite(rjb_magic, 1, 16, f) == 16 && std::fwrite(&n, sizeof(int), 1, f) == 1 && std::fwrite(a.data(), sizeof(float), NARR, f) == NARR && std::fclose(f) == 0;
    if (!ok)
        std::fprintf(stderr, "rjconvert: %s: write failed\n", out.c_str());
    return ok;
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::fprintf(stderr, "usage: rjconvert <indir> <outdir>\n");
        return 2;
    }
    DIR *d = opendir(argv[1]);
    if (!d)
    {
        std::fprintf(stderr, "rjconvert: cannot open %s\n", argv[1]);
        return 1;
    }
    int nconv = 0;
    struct dirent *e;
    while ((e = readdir(d)))
    {
        std::string name = e->d_name;
        if (name.size() < 5 || name.compare(0, 2, "rj") != 0 || name.compare(name.size() - 2, 2, ".i") != 0)
            continue;
        std::string in = std::string(argv[1]) + "/" + name;
        std::string out = std::string(argv[2]) + "/" + name.substr(0, name.size() - 2) + ".ib";
        if (!convert(in, out))
        {
            closedir(d);
            return 1;
        }
        nconv++;
    }
    closedir(d);
    std::printf("rjconvert: converted %d rj tables\n", nconv);
    return nconv > 0 ? 0 : 1;
}