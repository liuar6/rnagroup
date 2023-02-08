/* The MIT License (MIT)

   Copyright (c) 2023 Anrui Liu <liuar6@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   “Software”), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#include <stdlib.h>
#include <stdio.h>
#include "bioidx.h"
int main(){
    const char * n1 = "reg1";
    const char * n2 = "reg2";
    const char * n3 = "reg3";
    const char * n4 = "reg4";
    bioidx_t *bidx = bioidx_init();
    bioidx_itr_t *bitr ;
    bitr = bioidx_itr_init();
    bioidx_insert(bidx, 0, 200, 100000, (void *)n4);
    bioidx_insert(bidx, 0, 0, 100, (void *)n1);
    bioidx_insert(bidx, 0, 0, 100, (void *)n2);
    bioidx_insert(bidx, 0, 0, 1000000000, (void *)n3);
    const char *ret;
    bioidx_search(bidx, bitr, 0, 50, 1000);
    while ((ret = bioidx_itr_next(bitr)) != NULL) fprintf(stderr, "1:%s\n", ret);
    bioidx_search(bidx, bitr, 0, 150, 250);
    while ((ret = bioidx_itr_next(bitr)) != NULL) {
        bioidx_itr_remove(bitr);
        fprintf(stderr, "2:%s\n", ret);
    }
    bioidx_search(bidx, bitr, 0, 50, 1000);
    while ((ret = bioidx_itr_next(bitr)) != NULL) fprintf(stderr, "3:%s\n", ret);
    bioidx_destroy(bidx);
    bioidx_itr_destroy(bitr);
}