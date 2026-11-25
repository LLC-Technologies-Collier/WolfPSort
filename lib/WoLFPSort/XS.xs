/*
 * Copyright 2026 Google LLC
 * Licensed under the GNU General Public License v3 (GPLv3).
 */

#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "wolf_predict_core.h"
#include "wolf_predict_core.c"

typedef WolfModel * WoLFPSort__XS;

MODULE = WoLFPSort::XS  PACKAGE = WoLFPSort::XS

PROTOTYPES: DISABLE

WoLFPSort::XS
new(class, weight_path, training_path, utility_path)
    const char *class
    const char *weight_path
    const char *training_path
    const char *utility_path
  PREINIT:
    char err_buf[512];
    WolfModel *m;
  CODE:
    PERL_UNUSED_VAR(class);
    err_buf[0] = '\0';
    m = wolf_model_load(weight_path, training_path, utility_path, err_buf, sizeof(err_buf));
    if (!m) {
        croak("%s", err_buf[0] ? err_buf : "Failed to load WoLFPSort model");
    }
    RETVAL = m;
  OUTPUT:
    RETVAL

void
DESTROY(self)
    WoLFPSort::XS self
  CODE:
    wolf_model_free(self);

UV
k(self)
    WoLFPSort::XS self
  CODE:
    RETVAL = (UV)self->k_val;
  OUTPUT:
    RETVAL

UV
n_samples(self)
    WoLFPSort::XS self
  CODE:
    RETVAL = (UV)self->n_samples;
  OUTPUT:
    RETVAL

HV *
predict_vector(self, features_av, query_id_sv=&PL_sv_undef, ignore_self=0)
    WoLFPSort::XS self
    AV *features_av
    SV *query_id_sv
    int ignore_self
  PREINIT:
    SSize_t len;
    double *raw;
    WolfClassScore *scores;
    WolfNeighbor *nbs;
    size_t got, i;
    const char *qid;
    AV *scores_av;
    AV *nbs_av;
    HV *res_hv;
  CODE:
    len = av_len(features_av) + 1;
    if ((size_t)len != self->n_feat) {
        croak("Feature count mismatch: got %ld, expected %lu", (long)len, (unsigned long)self->n_feat);
    }
    Newx(raw, self->n_feat, double);
    for (i = 0; i < self->n_feat; i++) {
        SV **svp = av_fetch(features_av, (SSize_t)i, 0);
        raw[i] = (svp && *svp) ? SvNV(*svp) : 0.0;
    }
    qid = SvOK(query_id_sv) ? SvPV_nolen(query_id_sv) : NULL;

    Newx(scores, self->n_cls, WolfClassScore);
    Newx(nbs, self->k_val < self->n_samples ? self->k_val : self->n_samples, WolfNeighbor);
    got = wolf_model_predict(self, raw, qid, ignore_self, scores, nbs);

    scores_av = newAV();
    for (i = 0; i < self->n_cls; i++) {
        HV *entry = newHV();
        (void)hv_stores(entry, "class", newSVpv(self->class_names[scores[i].cls_idx], 0));
        (void)hv_stores(entry, "score", newSVnv(scores[i].score));
        av_push(scores_av, newRV_noinc((SV *)entry));
    }

    nbs_av = newAV();
    for (i = 0; i < got; i++) {
        HV *entry = newHV();
        const WolfSample *sp = &self->samples[nbs[i].sample_idx];
        (void)hv_stores(entry, "id", newSVpv(sp->id, 0));
        (void)hv_stores(entry, "class", newSVpv(sp->cls_name, 0));
        (void)hv_stores(entry, "similarity", newSVnv(nbs[i].similarity));
        av_push(nbs_av, newRV_noinc((SV *)entry));
    }

    Safefree(nbs);
    Safefree(scores);
    Safefree(raw);

    res_hv = newHV();
    (void)hv_stores(res_hv, "scores", newRV_noinc((SV *)scores_av));
    (void)hv_stores(res_hv, "neighbors", newRV_noinc((SV *)nbs_av));
    RETVAL = res_hv;
  OUTPUT:
    RETVAL
