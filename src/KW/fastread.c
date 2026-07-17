//
//  This file will be hopefully part of EPOS4
//  Copyright (C) 2026 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

/* Bulk-parse *n whitespace-separated single-precision values from fname
   starting at byte offset *ioff into a[0..*n-1]; on success advance *ioff
   past the last parsed token.  Conversion uses strtof - the same libc
   routine gfortran's list-directed READ calls internally  */

void fastreadf_(const char *fname, int *ioff, float *a, const int *n,
                int *ierr)
{
    *ierr = 1;
    int fd = open(fname, O_RDONLY);
    if (fd < 0)
        return;
    struct stat st;
    if (fstat(fd, &st) != 0 || *ioff < 0 || *ioff > st.st_size)
    {
        close(fd);
        return;
    }
    size_t len = (size_t)st.st_size;
    char *map = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (map == MAP_FAILED)
        return;

    const char *p = map + *ioff;
    const char *end = map + len;
    /* records end with a newline, so strtof cannot scan past the map */
    if (len == 0 || !isspace((unsigned char)end[-1]))
    {
        munmap(map, len);
        return;
    }

    *ierr = 2;
    for (int i = 0; i < *n; i++)
    {
        while (p < end && isspace((unsigned char)*p))
            p++;
        if (p >= end)
        {
            munmap(map, len);
            return;
        }
        char *q;
        a[i] = strtof(p, &q);
        if (q == p)
        {
            munmap(map, len);
            return;
        }
        p = q;
    }
    *ioff = (int)(p - map);
    munmap(map, len);
    *ierr = 0;
}

/* Read *n raw floats from a binary rj sidecar (see src/KW/rjconvert.cpp:
   16-byte magic, int32 total count, float payload) starting at byte offset
   *ioff (first call passes 20 = past magic+count); advance *ioff.  The
   sidecar values were produced by strtof of the ASCII tokens. */

static const char rjb_magic[16] = "EPOS4RJTAB-BIN1";

void fastreadfb_(const char *fname, int *ioff, float *a, const int *n,
                 int *ierr)
{
    *ierr = 1;
    int fd = open(fname, O_RDONLY);
    if (fd < 0)
        return;
    struct stat st;
    if (fstat(fd, &st) != 0)
    {
        close(fd);
        return;
    }
    size_t len = (size_t)st.st_size;
    char *map = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (map == MAP_FAILED)
        return;

    *ierr = 2;
    int total;
    if (len < 20 || memcmp(map, rjb_magic, 16) != 0)
    {
        munmap(map, len);
        return;
    }
    memcpy(&total, map + 16, sizeof(int));
    if (len != 20 + (size_t)total * sizeof(float) || *ioff < 20 || (size_t)*ioff + (size_t)*n * sizeof(float) > len)
    {
        munmap(map, len);
        return;
    }
    memcpy(a, map + *ioff, (size_t)*n * sizeof(float));
    *ioff += *n * (int)sizeof(float);
    munmap(map, len);
    *ierr = 0;
}