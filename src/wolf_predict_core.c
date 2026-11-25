/*
 * Copyright (C) 2004-2006 Paul B. Horton (algorithm & design, APBC06)
 * Copyright 2026 Google LLC (clean-room portable C99 core library)
 * Licensed under the GNU General Public License v3 (GPLv3).
 */

#include "wolf_predict_core.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WOLF_MAX_LINE 65536

typedef struct {
    size_t idx;
    double dist;
} WolfDistPair;

static char *xstrdup(const char *s) {
    size_t n = strlen(s) + 1;
    char  *p = (char *)malloc(n);
    if (!p) return NULL;
    memcpy(p, s, n);
    return p;
}

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static int cmp_dist_pair(const void *a, const void *b) {
    const WolfDistPair *na = (const WolfDistPair *)a;
    const WolfDistPair *nb = (const WolfDistPair *)b;
    if (na->dist < nb->dist) return -1;
    if (na->dist > nb->dist) return 1;
    return (na->idx > nb->idx) - (na->idx < nb->idx);
}

static int cmp_class_score(const void *a, const void *b) {
    const WolfClassScore *sa = (const WolfClassScore *)a;
    const WolfClassScore *sb = (const WolfClassScore *)b;
    if (sa->score > sb->score) return -1;
    if (sa->score < sb->score) return 1;
    return (sa->cls_idx > sb->cls_idx) - (sa->cls_idx < sb->cls_idx);
}

static void strip_comment_and_trim(char *line) {
    char *hash = strchr(line, '#');
    if (hash) *hash = '\0';
    size_t len = strlen(line);
    while (len > 0 && isspace((unsigned char)line[len - 1])) {
        line[--len] = '\0';
    }
}

static size_t split_ws(char *line, char ***out_tokens) {
    size_t cap = 64, count = 0;
    char **toks = (char **)malloc(cap * sizeof(char *));
    char  *p    = line;
    if (!toks) {
        *out_tokens = NULL;
        return 0;
    }
    while (*p) {
        while (*p && isspace((unsigned char)*p)) p++;
        if (!*p) break;
        char *start = p;
        while (*p && !isspace((unsigned char)*p)) p++;
        if (*p) *p++ = '\0';
        if (count == cap) {
            cap *= 2;
            char **nt = (char **)realloc(toks, cap * sizeof(char *));
            if (!nt) {
                free(toks);
                *out_tokens = NULL;
                return 0;
            }
            toks = nt;
        }
        toks[count++] = start;
    }
    *out_tokens = toks;
    return count;
}

/* Horton et al. (APBC06), Section 2.3.1: RLE empirical-CDF table per feature */
static WolfFeatureTable build_feature_table(const double *col, size_t n_samples) {
    double *sorted = (double *)malloc(n_samples * sizeof(double));
    memcpy(sorted, col, n_samples * sizeof(double));
    qsort(sorted, n_samples, sizeof(double), cmp_double);

    double *vals  = (double *)malloc(n_samples * sizeof(double));
    double *ranks = (double *)malloc(n_samples * sizeof(double));
    size_t  n_u = 0, i = 0;
    while (i < n_samples) {
        double v = sorted[i];
        size_t j = i + 1;
        while (j < n_samples && sorted[j] == v) j++;
        vals[n_u]  = v;
        ranks[n_u] = (0.5 * ((double)i + (double)j - 1.0) / (double)n_samples) * 100.0;
        n_u++;
        i = j;
    }
    free(sorted);
    WolfFeatureTable tbl = { vals, ranks, n_u };
    return tbl;
}

static double scale_value(const WolfFeatureTable *tbl, double x) {
    size_t lo = 0, hi = tbl->n_u;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (tbl->vals[mid] < x) lo = mid + 1;
        else hi = mid;
    }
    if (lo < tbl->n_u && tbl->vals[lo] == x) {
        return tbl->ranks[lo];
    }
    if (tbl->n_u <= 1) return tbl->ranks[0];

    /* Horton et al. (APBC06), Section 2.3.1, p. 3, para. 1, sentences 3-4:
     * linear interpolation (and extrapolation outside [min, max]). */
    size_t i0 = (lo == 0) ? 0 : (lo >= tbl->n_u ? tbl->n_u - 2 : lo - 1);
    size_t i1 = i0 + 1;
    double dv = tbl->vals[i1] - tbl->vals[i0];
    if (dv == 0.0) return tbl->ranks[i0];
    return tbl->ranks[i0] + ((x - tbl->vals[i0]) / dv) * (tbl->ranks[i1] - tbl->ranks[i0]);
}

void wolf_model_free(WolfModel *m) {
    if (!m) return;
    free(m->w1);
    free(m->w2);
    if (m->class_names) {
        for (size_t i = 0; i < m->n_cls; i++) free(m->class_names[i]);
        free(m->class_names);
    }
    free(m->utility);
    if (m->samples) {
        for (size_t i = 0; i < m->n_samples; i++) {
            free(m->samples[i].id);
            free(m->samples[i].cls_name);
            free(m->samples[i].features);
        }
        free(m->samples);
    }
    if (m->tables) {
        for (size_t f = 0; f < m->n_feat; f++) {
            free(m->tables[f].vals);
            free(m->tables[f].ranks);
        }
        free(m->tables);
    }
    free(m);
}

WolfModel *wolf_model_load(const char *weight_path,
                           const char *training_path,
                           const char *utility_path,
                           char       *err_buf,
                           size_t      err_buf_sz) {
    WolfModel *m = (WolfModel *)calloc(1, sizeof(WolfModel));
    if (!m) return NULL;

    char *line = (char *)malloc(WOLF_MAX_LINE);
    if (!line) {
        wolf_model_free(m);
        return NULL;
    }

    /* 1. Weights (.wolfw) */
    FILE *wf = fopen(weight_path, "r");
    if (!wf) {
        if (err_buf) snprintf(err_buf, err_buf_sz, "Cannot open weight file '%s'", weight_path);
        free(line);
        wolf_model_free(m);
        return NULL;
    }
    while (fgets(line, WOLF_MAX_LINE, wf)) {
        strip_comment_and_trim(line);
        if (!line[0]) continue;
        if (strncmp(line, "feature", 7) == 0 && isspace((unsigned char)line[7])) {
            char **toks = NULL;
            m->n_feat   = split_ws(line + 7, &toks);
            free(toks);
        } else if (strncmp(line, "w1", 2) == 0 && isspace((unsigned char)line[2])) {
            char **toks = NULL;
            size_t nt   = split_ws(line + 2, &toks);
            m->w1       = (double *)malloc(nt * sizeof(double));
            for (size_t i = 0; i < nt; i++) m->w1[i] = strtod(toks[i], NULL);
            free(toks);
        } else if (strncmp(line, "w2", 2) == 0 && isspace((unsigned char)line[2])) {
            char **toks = NULL;
            size_t nt   = split_ws(line + 2, &toks);
            m->w2       = (double *)malloc(nt * sizeof(double));
            for (size_t i = 0; i < nt; i++) m->w2[i] = strtod(toks[i], NULL);
            free(toks);
        } else if (strncmp(line, "k=", 2) == 0) {
            m->k_val = (size_t)strtoul(line + 2, NULL, 10);
        }
    }
    fclose(wf);
    if (!m->n_feat || !m->w1 || !m->w2 || !m->k_val) {
        if (err_buf) snprintf(err_buf, err_buf_sz, "Invalid weight file '%s'", weight_path);
        free(line);
        wolf_model_free(m);
        return NULL;
    }

    /* 2. Utility matrix (.wolfu) */
    FILE *uf = fopen(utility_path, "r");
    if (!uf) {
        if (err_buf) snprintf(err_buf, err_buf_sz, "Cannot open utility file '%s'", utility_path);
        free(line);
        wolf_model_free(m);
        return NULL;
    }
    size_t u_row = 0;
    while (fgets(line, WOLF_MAX_LINE, uf)) {
        strip_comment_and_trim(line);
        if (!line[0]) continue;
        char **toks = NULL;
        size_t nt   = split_ws(line, &toks);
        if (!m->class_names) {
            m->n_cls       = nt;
            m->class_names = (char **)calloc(m->n_cls, sizeof(char *));
            for (size_t i = 0; i < m->n_cls; i++) m->class_names[i] = xstrdup(toks[i]);
            m->utility = (double *)calloc(m->n_cls * m->n_cls, sizeof(double));
        } else if (u_row < m->n_cls) {
            for (size_t j = 0; j < nt && j < m->n_cls; j++) {
                m->utility[u_row * m->n_cls + j] = strtod(toks[j], NULL);
            }
            u_row++;
        }
        free(toks);
    }
    fclose(uf);
    if (!m->n_cls || u_row != m->n_cls) {
        if (err_buf) snprintf(err_buf, err_buf_sz, "Invalid utility matrix '%s'", utility_path);
        free(line);
        wolf_model_free(m);
        return NULL;
    }

    /* 3. Training data (.wolff) */
    FILE *tf = fopen(training_path, "r");
    if (!tf) {
        if (err_buf) snprintf(err_buf, err_buf_sz, "Cannot open training file '%s'", training_path);
        free(line);
        wolf_model_free(m);
        return NULL;
    }
    if (!fgets(line, WOLF_MAX_LINE, tf)) {
        fclose(tf);
        if (err_buf) snprintf(err_buf, err_buf_sz, "Empty training file '%s'", training_path);
        free(line);
        wolf_model_free(m);
        return NULL;
    }

    size_t cap = 1024;
    m->samples = (WolfSample *)calloc(cap, sizeof(WolfSample));
    while (fgets(line, WOLF_MAX_LINE, tf)) {
        strip_comment_and_trim(line);
        if (!line[0]) continue;
        char **toks = NULL;
        size_t nt   = split_ws(line, &toks);
        if (nt < m->n_feat + 2) {
            free(toks);
            continue;
        }
        size_t c_idx = m->n_cls;
        for (size_t c = 0; c < m->n_cls; c++) {
            if (strcmp(m->class_names[c], toks[1]) == 0) {
                c_idx = c;
                break;
            }
        }
        if (c_idx == m->n_cls) {
            if (err_buf) snprintf(err_buf, err_buf_sz, "Unknown class '%s'", toks[1]);
            free(toks);
            fclose(tf);
            free(line);
            wolf_model_free(m);
            return NULL;
        }
        if (m->n_samples == cap) {
            cap *= 2;
            m->samples = (WolfSample *)realloc(m->samples, cap * sizeof(WolfSample));
        }
        double *fvec = (double *)malloc(m->n_feat * sizeof(double));
        for (size_t f = 0; f < m->n_feat; f++) fvec[f] = strtod(toks[f + 2], NULL);
        m->samples[m->n_samples].id       = xstrdup(toks[0]);
        m->samples[m->n_samples].cls_name = xstrdup(toks[1]);
        m->samples[m->n_samples].cls_idx  = c_idx;
        m->samples[m->n_samples].features = fvec;
        m->n_samples++;
        free(toks);
    }
    fclose(tf);
    free(line);

    /* 4. Build UniformScaler tables & scale training features */
    m->tables       = (WolfFeatureTable *)calloc(m->n_feat, sizeof(WolfFeatureTable));
    double *col_buf = (double *)malloc(m->n_samples * sizeof(double));
    for (size_t f = 0; f < m->n_feat; f++) {
        for (size_t i = 0; i < m->n_samples; i++) col_buf[i] = m->samples[i].features[f];
        m->tables[f] = build_feature_table(col_buf, m->n_samples);
        for (size_t i = 0; i < m->n_samples; i++) {
            m->samples[i].features[f] = scale_value(&m->tables[f], m->samples[i].features[f]);
        }
    }
    free(col_buf);

    return m;
}

size_t wolf_model_predict(const WolfModel *m,
                          const double    *raw_features,
                          const char      *query_id,
                          int              ignore_self,
                          WolfClassScore  *out_scores,
                          WolfNeighbor    *out_neighbors) {
    double       *q_scaled = (double *)malloc(m->n_feat * sizeof(double));
    WolfDistPair *nbs      = (WolfDistPair *)malloc(m->n_samples * sizeof(WolfDistPair));
    double       *counts   = (double *)calloc(m->n_cls, sizeof(double));

    for (size_t f = 0; f < m->n_feat; f++) {
        q_scaled[f] = scale_value(&m->tables[f], raw_features[f]);
    }

    /* Horton et al. (APBC06), Section 2.3.2, p. 4, para. 1, Equation 1 */
    for (size_t i = 0; i < m->n_samples; i++) {
        double        l1   = 0.0;
        double        l2sq = 0.0;
        const double *s    = m->samples[i].features;
        for (size_t f = 0; f < m->n_feat; f++) {
            double d = s[f] - q_scaled[f];
            if (m->w1[f] != 0.0) l1 += m->w1[f] * fabs(d);
            if (m->w2[f] != 0.0) l2sq += m->w2[f] * d * d;
        }
        nbs[i].idx  = i;
        nbs[i].dist = l1 + sqrt(l2sq);
    }
    qsort(nbs, m->n_samples, sizeof(WolfDistPair), cmp_dist_pair);

    size_t k_eff = m->k_val < m->n_samples ? m->k_val : m->n_samples;
    size_t got   = 0;
    for (size_t i = 0; i < m->n_samples && got < k_eff; i++) {
        size_t idx = nbs[i].idx;
        if (ignore_self && query_id && strcmp(m->samples[idx].id, query_id) == 0) continue;
        out_neighbors[got].sample_idx = idx;
        out_neighbors[got].similarity = -nbs[i].dist;
        counts[m->samples[idx].cls_idx] += 1.0;
        got++;
    }

    /* Horton et al. (APBC06), Section 2.3.3, p. 4, para. 1 & Table 3 */
    for (size_t c = 0; c < m->n_cls; c++) {
        double s = 0.0;
        for (size_t j = 0; j < m->n_cls; j++) {
            s += m->utility[c * m->n_cls + j] * counts[j];
        }
        out_scores[c].cls_idx = c;
        out_scores[c].score   = s;
    }
    qsort(out_scores, m->n_cls, sizeof(WolfClassScore), cmp_class_score);

    free(counts);
    free(nbs);
    free(q_scaled);
    return got;
}
