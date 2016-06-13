#include "encint.h"

void od_encode_checkpoint(const daala_enc_ctx *enc, od_rollback_buffer *rbuf) {
  od_ec_enc_checkpoint(&rbuf->ec, &enc->ec);
  OD_COPY(&rbuf->adapt, &enc->state.adapt, 1);
}

void od_encode_rollback(daala_enc_ctx *enc, const od_rollback_buffer *rbuf) {
  od_ec_enc_rollback(&enc->ec, &rbuf->ec);
  OD_COPY(&enc->state.adapt, &rbuf->adapt, 1);
}
