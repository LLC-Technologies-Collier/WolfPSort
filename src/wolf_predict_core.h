/*
 * Copyright (C) 2004-2006 Paul B. Horton (algorithm & design, APBC06)
 * Copyright 2026 Google LLC (clean-room portable C99 core library)
 * Licensed under the GNU General Public License v3 (GPLv3).
 */

#ifndef WOLF_PREDICT_CORE_H
#define WOLF_PREDICT_CORE_H

#include <stddef.h>

typedef struct {
    double *vals;
    double *ranks;
    size_t  n_u;
} WolfFeatureTable;

typedef struct {
    char   *id;
    char   *cls_name;
    size_t  cls_idx;
    double *features;
} WolfSample;

typedef struct {
    size_t  cls_idx;
    double  score;
} WolfClassScore;

typedef struct {
    size_t  sample_idx;
    double  similarity;
} WolfNeighbor;

typedef struct {
    size_t            n_feat;
    size_t            k_val;
    double           *w1;
    double           *w2;
    size_t            n_cls;
    char            **class_names;
    double           *utility;
    size_t            n_samples;
    WolfSample       *samples;
    WolfFeatureTable *tables;
} WolfModel;

WolfModel *wolf_model_load(const char *weight_path,
                           const char *training_path,
                           const char *utility_path,
                           char       *err_buf,
                           size_t      err_buf_sz);

void wolf_model_free(WolfModel *m);

/*
 * Predict localization scores and k nearest neighbors for a raw feature vector.
 * out_scores must have length m->n_cls; out_neighbors must have length k_eff.
 * Returns number of neighbors written to out_neighbors.
 */
size_t wolf_model_predict(const WolfModel *m,
                          const double    *raw_features,
                          const char      *query_id,
                          int              ignore_self,
                          WolfClassScore  *out_scores,
                          WolfNeighbor    *out_neighbors);

#endif /* WOLF_PREDICT_CORE_H */
