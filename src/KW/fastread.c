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