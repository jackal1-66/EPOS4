//
//  This file will hopefully be part of EPOS4
//  Copyright (C) 2026 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

// Build-time converter: parses an ASCII eos3f table (the gzipped
// eos4f.eos) exactly the way EoS3f's constructor does, and writes
// it as the binary format EoS3f loads directly. Values are
// stored raw as parsed
//
// usage: eosconvert <ascii-in> <binary-out>

#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

static const char eosb_magic[16] = "EPOS4EOS3F-BIN1"; // 15 chars + NUL

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::fprintf(stderr, "usage: eosconvert <ascii-in> <binary-out>\n");
        return 2;
    }

    std::string buf;
    {
        std::ifstream fin(argv[1], std::ios::binary);
        if (!fin)
        {
            std::fprintf(stderr, "eosconvert: cannot open %s\n", argv[1]);
            return 1;
        }
        fin.seekg(0, std::ios::end);
        std::streamoff len = fin.tellg();
        if (len < 0)
            len = 0;
        buf.resize((size_t)len);
        fin.seekg(0, std::ios::beg);
        if (len > 0)
            fin.read(&buf[0], len);
    }
    const char *sp = buf.data();
    const char *se = sp + buf.size();
    bool okread = true;
    auto nextd = [&](double &v)
    {
        while (sp < se && (*sp == ' ' || *sp == '\n' || *sp == '\t' || *sp == '\r' || *sp == '\v' || *sp == '\f'))
            ++sp;
        std::from_chars_result r = std::from_chars(sp, se, v);
        if (r.ec != std::errc())
        {
            v = 0;
            okread = false;
            return;
        }
        sp = r.ptr;
    };
    auto nexti = [&](int &v)
    {
        while (sp < se && (*sp == ' ' || *sp == '\n' || *sp == '\t' || *sp == '\r' || *sp == '\v' || *sp == '\f'))
            ++sp;
        std::from_chars_result r = std::from_chars(sp, se, v);
        if (r.ec != std::errc())
        {
            v = 0;
            okread = false;
            return;
        }
        sp = r.ptr;
    };

    double B, volex0, delta0, aaa, bbb, emax, e0;
    int ne, nn;
    nextd(B);
    nextd(volex0);
    nextd(delta0);
    nextd(aaa);
    nextd(bbb);
    nextd(emax);
    nextd(e0);
    nexti(ne);
    nexti(nn);
    if (!okread || ne <= 0 || nn <= 0)
    {
        std::fprintf(stderr, "eosconvert: bad header\n");
        return 1;
    }

    const size_t npts = (size_t)ne * nn * nn * nn;
    std::vector<double> egrid((size_t)ne), ngrid((size_t)ne * nn);
    std::vector<double> pre(npts), T(npts), mub(npts), muq(npts), mus(npts);

    for (int ixe = 0; ixe < ne; ixe++)
        nextd(egrid[ixe]);
    for (size_t i = 0; i < (size_t)ne * nn; i++)
        nextd(ngrid[i]);

    // same element order as EoS3f: index = ie + ne*ib + ne*nn*iq + ne*nn*nn*is
    for (int ie = 0; ie < ne; ie++)
        for (int inb = 0; inb < nn; inb++)
            for (int inq = 0; inq < nn; inq++)
                for (int ins = 0; ins < nn; ins++)
                {
                    const size_t idx = (size_t)ie + (size_t)ne * inb + (size_t)ne * nn * inq + (size_t)ne * nn * nn * ins;
                    nextd(pre[idx]);
                    nextd(T[idx]);
                    nextd(mub[idx]);
                    nextd(muq[idx]);
                    nextd(mus[idx]);
                    if (!okread)
                    {
                        std::fprintf(stderr, "eosconvert: parse failed near point %d/%d/%d/%d\n", ie, inb, inq, ins);
                        return 1;
                    }
                }

    std::FILE *out = std::fopen(argv[2], "wb");
    if (!out)
    {
        std::fprintf(stderr, "eosconvert: cannot open %s\n", argv[2]);
        return 1;
    }
    bool ok = true;
    ok = ok && std::fwrite(eosb_magic, 1, 16, out) == 16;
    int dims[2] = {ne, nn};
    double hdr[7] = {B, volex0, delta0, aaa, bbb, emax, e0};
    ok = ok && std::fwrite(dims, sizeof(int), 2, out) == 2;
    ok = ok && std::fwrite(hdr, sizeof(double), 7, out) == 7;
    ok = ok && std::fwrite(egrid.data(), sizeof(double), egrid.size(), out) == egrid.size();
    ok = ok && std::fwrite(ngrid.data(), sizeof(double), ngrid.size(), out) == ngrid.size();
    ok = ok && std::fwrite(pre.data(), sizeof(double), npts, out) == npts;
    ok = ok && std::fwrite(T.data(), sizeof(double), npts, out) == npts;
    ok = ok && std::fwrite(mub.data(), sizeof(double), npts, out) == npts;
    ok = ok && std::fwrite(muq.data(), sizeof(double), npts, out) == npts;
    ok = ok && std::fwrite(mus.data(), sizeof(double), npts, out) == npts;
    ok = ok && std::fclose(out) == 0;
    if (!ok)
    {
        std::fprintf(stderr, "eosconvert: write failed\n");
        return 1;
    }
    std::printf("eosconvert: %s -> %s (ne=%d nn=%d, %zu points)\n", argv[1], argv[2], ne, nn, npts);
    return 0;
}