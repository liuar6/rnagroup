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
/* "vector.h" is a very simple implementation of sequence container in
  order to adopt a similar API like the well known "khash.h" */

#define VECTOR_VERSION 0.1.0

#define VECTOR_EXTENSION_FACTOR 0.5
#define VECTOR_MIN_CAPACITY 1


#define VEC_DECLARE(name, element_t) \
\
typedef struct vec_##name##_s{   \
    size_t size;    \
    size_t capacity;    \
    element_t *data;    \
}vec_##name##_t; \
vec_##name##_t *vec_init_##name(); \
int vec_add_##name(vec_##name##_t* v, element_t element); \
int vec_clear_##name(vec_##name##_t* v);  \
int vec_clear_##name(vec_##name##_t* v);  \
int vec_destroy_##name(vec_##name##_t* v);

#define VEC_INIT(name, element_t)   \
    \
typedef struct vec_##name##_s{   \
    size_t size;    \
    size_t capacity;    \
    element_t *data;    \
}vec_##name##_t; \
    \
static vec_##name##_t *vec_init_##name(){   \
    return(calloc(sizeof(vec_##name##_t), 1));   \
}   \
    \
static int vec_add_##name(vec_##name##_t* v, element_t element){   \
    if (v->size+1 > v->capacity){   \
        if (v->data == NULL) {  \
            v->data = malloc(VECTOR_MIN_CAPACITY * sizeof(element_t));    \
            v->capacity = VECTOR_MIN_CAPACITY;\
        }  \
        else {  \
            int extension=(int)(v->capacity) * VECTOR_EXTENSION_FACTOR; \
            if (extension < 1) extension = 1;   \
            element_t* new_data = realloc(v->data, (v->capacity + extension) * sizeof(element_t));    \
            if (new_data == NULL) return -1;    \
            v->data = new_data; \
            v->capacity = v->capacity + extension;  \
        }   \
    }   \
    v->data[v->size] = element; \
    v->size++;  \
    return 0;   \
}   \
static int vec_clear_##name(vec_##name##_t* v){  \
    v->size = 0;  \
    return 0; \
}   \
static int vec_destroy_##name(vec_##name##_t* v){  \
    free(v->data);  \
    free(v);    \
    return 0; \
}

#define vec_t(name) vec_##name##_t
#define vec_init(name) vec_init_##name()
#define vec_destroy(name, h) vec_destroy_##name(h)
#define vec_add(name, h, val) vec_add_##name(h, val)
#define vec_clear(name, h) vec_clear_##name(h)
#define vec_val(h, i) ((h)->data[(i)])