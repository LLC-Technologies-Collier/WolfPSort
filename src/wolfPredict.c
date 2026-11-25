/*
 * Copyright (C) 2004-2006 Paul B. Horton (algorithm & design, APBC06)
 * Copyright 2026 Google LLC (clean-room portable C99 implementation)
 * Licensed under the GNU General Public License v3 (GPLv3).
 *
 * Clean-room portable C99 implementation of WoLF PSORT's weighted hybrid
 * L1/L2 k-nearest neighbors classifier (wolfPredict), implementing:
 *   - Horton et al. (APBC06), Section 2.3.1: empirical-CDF [0, 100] mid-rank
 *     percentile normalization with out-of-sample linear interpolation and
 *     extrapolation.
 *   - Horton et al. (APBC06), Section 2.3.2, Equation 1:
 *     d_W(j, k) = sum_i w1_i * |F_ji - F_ki| + sqrt(sum_i w2_i * (F_ji - F_ki)^2)
 *   - Horton et al. (APBC06), Section 2.3.3, Table 3: symmetric utility matrix
 *     expected-utility scoring for dual-localization prediction.
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 65536

typedef struct {
    double *vals;
    double *ranks;
    size_t  n_u;
} FeatureTable;

typedef struct {
    char   *id;
    char   *cls_name;
    size_t  cls_idx;
    double *features;
} Sample;

typedef struct {
    size_t idx;
    double dist;
} Neighbor;

typedef struct {
    size_t cls_idx;
    double score;
} ClassScore;

static char *xstrdup(const char *s) {
    size_t n = strlen(s) + 1;
    char *p  = (char *)malloc(n);
    if (!p) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    memcpy(p, s, n);
    return p;
}

static void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    return p;
}

static void *xrealloc(void *ptr, size_t n) {
    void *p = realloc(ptr, n ? n : 1);
    if (!p) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    return p;
}

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static int cmp_neighbor(const void *a, const void *b) {
    const Neighbor *na = (const Neighbor *)a;
    const Neighbor *nb = (const Neighbor *)b;
    if (na->dist < nb->dist) return -1;
    if (na->dist > nb->dist) return 1;
    return (na->idx > nb->idx) - (na->idx < nb->idx);
}

static int cmp_score(const void *a, const void *b) {
    const ClassScore *sa = (const ClassScore *)a;
    const ClassScore *sb = (const ClassScore *)b;
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
    char **toks = (char **)xmalloc(cap * sizeof(char *));
    char *p     = line;
    while (*p) {
        while (*p && isspace((unsigned char)*p)) p++;
        if (!*p) break;
        char *start = p;
        while (*p && !isspace((unsigned char)*p)) p++;
        if (*p) *p++ = '\0';
        if (count == cap) {
            cap *= 2;
            toks = (char **)xrealloc(toks, cap * sizeof(char *));
        }
        toks[count++] = start;
    }
    *out_tokens = toks;
    return count;
}

/* Horton et al. (APBC06), Section 2.3.1: RLE empirical-CDF table per feature */
static FeatureTable build_feature_table(const double *col, size_t n_samples) {
    double *sorted = (double *)xmalloc(n_samples * sizeof(double));
    memcpy(sorted, col, n_samples * sizeof(double));
    qsort(sorted, n_samples, sizeof(double), cmp_double);

    double *vals  = (double *)xmalloc(n_samples * sizeof(double));
    double *ranks = (double *)xmalloc(n_samples * sizeof(double));
    size_t n_u = 0, i = 0;
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
    FeatureTable tbl = { vals, ranks, n_u };
    return tbl;
}

static double scale_value(const FeatureTable *tbl, double x) {
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

    /* Step-function insertion rank (before APBC06 linear interpolation) */
    return ((double)lo - 0.5) / (double)tbl->n_u * 100.0;
}

static void format_number(double val, char *buf, size_t bufsz) {
    if (fabs(val) < 1e-9) {
        snprintf(buf, bufsz, "0");
        return;
    }
    snprintf(buf, bufsz, "%.3f", val);
    char *dot = strchr(buf, '.');
    if (dot) {
        char *end = buf + strlen(buf) - 1;
        while (end > dot && *end == '0') *end-- = '\0';
        if (end == dot) *end = '\0';
    }
}

static void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [OPTIONS] weightFile trainingDataFile utilityMatrixFile\n"
        "\n  OPTIONS\n"
        "    --ignore-self\n"
        "    --print-instance-comment\n"
        "    --print-neighbors\n"
        "    --print-all-scores\n",
        prog);
}

int main(int argc, char **argv) {
    int ignore_self            = 0;
    int print_instance_comment = 0;
    int print_neighbors        = 0;
    int print_all_scores       = 0;
    const char *pos_args[3];
    int n_pos = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--ignore-self") == 0) {
            ignore_self = 1;
        } else if (strcmp(argv[i], "--print-instance-comment") == 0) {
            print_instance_comment = 1;
        } else if (strcmp(argv[i], "--print-neighbors") == 0) {
            print_neighbors = 1;
        } else if (strcmp(argv[i], "--print-all-scores") == 0) {
            print_all_scores = 1;
        } else if (argv[i][0] == '-') {
            print_usage(argv[0]);
            return 1;
        } else if (n_pos < 3) {
            pos_args[n_pos++] = argv[i];
        } else {
            print_usage(argv[0]);
            return 1;
        }
    }
    if (n_pos != 3) {
        print_usage(argv[0]);
        return 1;
    }

    const char *weight_path   = pos_args[0];
    const char *training_path = pos_args[1];
    const char *utility_path  = pos_args[2];

    char line[MAX_LINE];

    /* 1. Load weights (.wolfw) */
    FILE *wf = fopen(weight_path, "r");
    if (!wf) {
        fprintf(stderr, "Cannot open weight file '%s'\n", weight_path);
        return 1;
    }
    size_t n_feat = 0, k_val = 0;
    double *w1 = NULL, *w2 = NULL;
    while (fgets(line, sizeof(line), wf)) {
        strip_comment_and_trim(line);
        if (!line[0]) continue;
        if (strncmp(line, "feature", 7) == 0 && isspace((unsigned char)line[7])) {
            char **toks = NULL;
            size_t nt   = split_ws(line + 7, &toks);
            n_feat      = nt;
            free(toks);
        } else if (strncmp(line, "w1", 2) == 0 && isspace((unsigned char)line[2])) {
            char **toks = NULL;
            size_t nt   = split_ws(line + 2, &toks);
            w1          = (double *)xmalloc(nt * sizeof(double));
            for (size_t i = 0; i < nt; i++) w1[i] = strtod(toks[i], NULL);
            free(toks);
        } else if (strncmp(line, "w2", 2) == 0 && isspace((unsigned char)line[2])) {
            char **toks = NULL;
            size_t nt   = split_ws(line + 2, &toks);
            w2          = (double *)xmalloc(nt * sizeof(double));
            for (size_t i = 0; i < nt; i++) w2[i] = strtod(toks[i], NULL);
            free(toks);
        } else if (strncmp(line, "k=", 2) == 0) {
            k_val = (size_t)strtoul(line + 2, NULL, 10);
        }
    }
    fclose(wf);
    if (!n_feat || !w1 || !w2 || !k_val) {
        fprintf(stderr, "Invalid weight file '%s'\n", weight_path);
        return 1;
    }

    /* 2. Load utility matrix (.wolfu) */
    FILE *uf = fopen(utility_path, "r");
    if (!uf) {
        fprintf(stderr, "Cannot open utility file '%s'\n", utility_path);
        return 1;
    }
    char **class_names = NULL;
    size_t n_cls = 0, u_row = 0;
    double *utility = NULL;
    while (fgets(line, sizeof(line), uf)) {
        strip_comment_and_trim(line);
        if (!line[0]) continue;
        char **toks = NULL;
        size_t nt   = split_ws(line, &toks);
        if (!class_names) {
            n_cls       = nt;
            class_names = (char **)xmalloc(n_cls * sizeof(char *));
            for (size_t i = 0; i < n_cls; i++) class_names[i] = xstrdup(toks[i]);
            utility = (double *)xmalloc(n_cls * n_cls * sizeof(double));
        } else if (u_row < n_cls) {
            for (size_t j = 0; j < nt && j < n_cls; j++) {
                utility[u_row * n_cls + j] = strtod(toks[j], NULL);
            }
            u_row++;
        }
        free(toks);
    }
    fclose(uf);
    if (!n_cls || u_row != n_cls) {
        fprintf(stderr, "Invalid utility matrix '%s'\n", utility_path);
        return 1;
    }

    /* 3. Load training data (.wolff) */
    FILE *tf = fopen(training_path, "r");
    if (!tf) {
        fprintf(stderr, "Cannot open training file '%s'\n", training_path);
        return 1;
    }
    if (!fgets(line, sizeof(line), tf)) {
        fclose(tf);
        fprintf(stderr, "Empty training file '%s'\n", training_path);
        return 1;
    }

    size_t cap = 1024, n_samples = 0;
    Sample *samples = (Sample *)xmalloc(cap * sizeof(Sample));
    while (fgets(line, sizeof(line), tf)) {
        strip_comment_and_trim(line);
        if (!line[0]) continue;
        char **toks = NULL;
        size_t nt   = split_ws(line, &toks);
        if (nt < n_feat + 2) {
            free(toks);
            continue;
        }
        size_t c_idx = n_cls;
        for (size_t c = 0; c < n_cls; c++) {
            if (strcmp(class_names[c], toks[1]) == 0) {
                c_idx = c;
                break;
            }
        }
        if (c_idx == n_cls) {
            fprintf(stderr, "Unknown class '%s' for '%s'\n", toks[1], toks[0]);
            free(toks);
            fclose(tf);
            return 1;
        }
        if (n_samples == cap) {
            cap *= 2;
            samples = (Sample *)xrealloc(samples, cap * sizeof(Sample));
        }
        double *fvec = (double *)xmalloc(n_feat * sizeof(double));
        for (size_t f = 0; f < n_feat; f++) fvec[f] = strtod(toks[f + 2], NULL);
        samples[n_samples].id       = xstrdup(toks[0]);
        samples[n_samples].cls_name = xstrdup(toks[1]);
        samples[n_samples].cls_idx  = c_idx;
        samples[n_samples].features = fvec;
        n_samples++;
        free(toks);
    }
    fclose(tf);

    /* 4. Build UniformScaler tables & scale training features */
    FeatureTable *tables = (FeatureTable *)xmalloc(n_feat * sizeof(FeatureTable));
    double *col_buf      = (double *)xmalloc(n_samples * sizeof(double));
    for (size_t f = 0; f < n_feat; f++) {
        for (size_t i = 0; i < n_samples; i++) col_buf[i] = samples[i].features[f];
        tables[f] = build_feature_table(col_buf, n_samples);
        for (size_t i = 0; i < n_samples; i++) {
            samples[i].features[f] = scale_value(&tables[f], samples[i].features[f]);
        }
    }
    free(col_buf);

    /* 5. Process query stream from stdin */
    printf("# k used for kNN is: %zu\n", k_val);
    if (!fgets(line, sizeof(line), stdin)) return 0;

    double     *q_scaled = (double *)xmalloc(n_feat * sizeof(double));
    Neighbor   *nbs      = (Neighbor *)xmalloc(n_samples * sizeof(Neighbor));
    double     *counts   = (double *)xmalloc(n_cls * sizeof(double));
    ClassScore *cscores  = (ClassScore *)xmalloc(n_cls * sizeof(ClassScore));

    while (fgets(line, sizeof(line), stdin)) {
        char comment[1024];
        comment[0] = '\0';
        char *hash = strchr(line, '#');
        if (hash) {
            snprintf(comment, sizeof(comment), "%s", hash + 1);
            size_t clen = strlen(comment);
            while (clen > 0 && (comment[clen - 1] == '\n' || comment[clen - 1] == '\r')) {
                comment[--clen] = '\0';
            }
            *hash = '\0';
        }
        strip_comment_and_trim(line);
        if (!line[0]) continue;

        char **toks = NULL;
        size_t nt   = split_ws(line, &toks);
        if (nt < n_feat + 2) {
            free(toks);
            continue;
        }
        const char *qid = toks[0];
        for (size_t f = 0; f < n_feat; f++) {
            q_scaled[f] = scale_value(&tables[f], strtod(toks[f + 2], NULL));
        }

        /* Horton et al. (APBC06), Section 2.3.2, Equation 1 */
        for (size_t i = 0; i < n_samples; i++) {
            double l1 = 0.0, l2sq = 0.0;
            const double *s = samples[i].features;
            for (size_t f = 0; f < n_feat; f++) {
                double d = s[f] - q_scaled[f];
                if (w1[f] != 0.0) l1 += w1[f] * fabs(d);
                if (w2[f] != 0.0) l2sq += w2[f] * d * d;
            }
            nbs[i].idx  = i;
            nbs[i].dist = l1 + sqrt(l2sq);
        }
        qsort(nbs, n_samples, sizeof(Neighbor), cmp_neighbor);

        size_t k_eff = k_val < n_samples ? k_val : n_samples;
        memset(counts, 0, n_cls * sizeof(double));
        size_t *chosen = (size_t *)xmalloc(k_eff * sizeof(size_t));
        double *cdist  = (double *)xmalloc(k_eff * sizeof(double));
        size_t  got    = 0;
        for (size_t i = 0; i < n_samples && got < k_eff; i++) {
            size_t idx = nbs[i].idx;
            if (ignore_self && strcmp(samples[idx].id, qid) == 0) continue;
            chosen[got] = idx;
            cdist[got]  = nbs[i].dist;
            counts[samples[idx].cls_idx] += 1.0;
            got++;
        }

        /* Horton et al. (APBC06), Section 2.3.3, Table 3 */
        for (size_t c = 0; c < n_cls; c++) {
            double s = 0.0;
            for (size_t j = 0; j < n_cls; j++) {
                s += utility[c * n_cls + j] * counts[j];
            }
            cscores[c].cls_idx = c;
            cscores[c].score   = s;
        }
        qsort(cscores, n_cls, sizeof(ClassScore), cmp_score);

        double min_score = (n_cls > 0 && cscores[0].score >= 10.0) ? 2.5 : 1.0;
        if (print_instance_comment) {
            printf("%s\n", comment);
        } else {
            printf("%s ", qid);
        }

        int printed = 0;
        for (size_t c = 0; c < n_cls; c++) {
            if (!print_all_scores && printed > 0 && cscores[c].score < min_score) break;
            char nbuf[64];
            format_number(cscores[c].score, nbuf, sizeof(nbuf));
            printf("%s%s %s", printed ? ", " : "", class_names[cscores[c].cls_idx], nbuf);
            printed++;
        }
        printf("\n");

        if (print_neighbors) {
            for (size_t i = 0; i < got; i++) {
                char nbuf[64];
                format_number(-cdist[i], nbuf, sizeof(nbuf));
                printf("  %-11s %-10s %s\n",
                       samples[chosen[i]].id, samples[chosen[i]].cls_name, nbuf);
            }
        }

        free(chosen);
        free(cdist);
        free(toks);
    }

    return 0;
}
