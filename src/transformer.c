#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <assert.h>
#include "transformer.h"


// for testing dim as dipicted in http://jalammar.github.io/illustrated-transformer/
#define N_SEQ 2 // number of token size, i.e. sequence length
#define DIM_EMBEDDING 4 // embedding dimension
#define DIM_INTERNAL 3  // internal dimension

#define DEBUG 1

const char* fmt = "%8" SML_PRTREAL "g ";

// simple soft max in axis=1
void soft_max(Tensor a) {
  size_t i, j, rows=MatRows(a), cols=MatCols(a);

  // for numerical stability, max of row can be used per row instead
  float largest = MatMax(a);

  for (i=0; i != rows; ++i) {
    float sum_exp = 0.;
    for (j = 0; j != cols; ++j) {
      float e = sml_exp(MatGet0(a, i, j) - largest);
      sum_exp += e;
      MatSet0(a, i, j, e);
    }
    for (j = 0; j != cols; ++j) {
      float e = MatGet0(a, i, j)/sum_exp;
      MatSet0(a, i, j, e);
    }
#if DEBUG
    float sum = 0.;
    for (j = 0; j != cols; ++j) {
      sum = sum + MatGet0(a, i, j);
    }
    assert(sum == 1.);
#endif
  }
}

/*
 * self attention block
 *
 * input sa: self attention structure
 * input A: input Tensor (k x d) where k number of token, d is embedding dim
 */
Tensor self_attention(struct SelfAtten *sa, Tensor x) {
  // NOTE: skip mask
  // TODO, dim should be carried by sa
  struct Transformer *trfm = sa->trfm;

  Tensor q, k, v;
  q = MatDim(trfm->n_seq, trfm->dim_internal);
  k = MatDim(trfm->n_seq, trfm->dim_internal);
  v = MatDim(trfm->n_seq, trfm->dim_internal);

#if DEBUG
  MatWrite(stdout , x, fmt,"x:");
  MatWrite(stdout , sa->w_q, fmt,"w_q:");
  MatWrite(stdout , sa->w_k, fmt,"w_k:");
  MatWrite(stdout , sa->w_v, fmt,"w_v:");
#endif

  MatMul(q, x, sa->w_q);
  MatMul(k, x, sa->w_k);
  MatMul(v, x, sa->w_v);

#if DEBUG
  MatWrite(stdout , q, fmt,"q:");
  MatWrite(stdout , k, fmt,"k:");
  MatWrite(stdout , v, fmt,"v:");
#endif

  Tensor qk_t, qk_t_out, z, k_t;
  qk_t = MatDim(trfm->n_seq, trfm->n_seq);
  qk_t_out = MatDim(trfm->n_seq, trfm->n_seq);
  k_t = MatDim(MatCols(k), MatRows(k)); // k transpose

  // TODO: create API
  // z = MatDimLike(v);
  z = MatDim(MatRows(v), MatCols(v));

  // transpose
  MatTran(k_t, k);
#if DEBUG
  MatWrite(stdout , v, fmt,"k_t:");
#endif

  MatMul(qk_t, q, k_t);
#if DEBUG
  MatWrite(stdout , qk_t, fmt,"qk_t:");
#endif

  // TODO API for div by scalar
  // element wise division
  float d_sqrt = sqrt(sa->trfm->dim_embedding);
  // create a wrapper for element wise division
  Tensor d_sqrt_t;
  d_sqrt_t = MatDim(MatRows(qk_t), MatCols(qk_t));
  MatFill(d_sqrt_t, d_sqrt);
  MatDivE0(qk_t_out, qk_t, d_sqrt_t);

#if DEBUG
  MatWrite(stdout , qk_t_out, fmt,"qk_t_out:");
#endif

  // softmax probability
  soft_max(qk_t_out);

#if DEBUG
  MatWrite(stdout , qk_t_out, fmt,"qk_t_out softmax:");
#endif

  MatMul(z, qk_t_out, v);

  // free
  MatUnDim(q);
  MatUnDim(k);
  MatUnDim(v);

  MatUnDim(k_t);
  MatUnDim(qk_t);
  MatUnDim(qk_t_out);
  MatUnDim(d_sqrt_t);

#if DEBUG
  MatWrite(stdout , z, fmt,"z:");
#endif
  return z;
}

/*
 * multi head attention
 *
 * input: transformer struct
 * input: x input
 */ 
Tensor multi_head_attention(struct Transformer *trfm, Tensor x) {

  Tensor z_sa[N_MHA];
  for (int i = 0; i < N_MHA; i++) {
    z_sa[i] = self_attention(&trfm->sa[i], x);
#if DEBUG
    printf("%s:%d head=%d\n", __func__ , __LINE__, i);
    MatWrite(stdout , z_sa[i], fmt,"z_sa:");
#endif
  }

  Tensor z_sa_cat, z;
  z_sa_cat = MatDim(MatRows(z_sa[0]), MatCols(z_sa[0])*trfm->n_mha);
  z = MatDim(trfm->n_seq, trfm->dim_embedding);

  // matrix concat
  // axis=1: horizontal concat
  // maybe there is a better way of concat
  for (int i = 0; i < trfm->n_mha; i++) {
    for (int r = 0; r < MatRows(z_sa[0]); r++) {
      for (int c = 0; c < MatCols(z_sa[0]); c++) {
        float v = MatGet0(z_sa[i], r, c);
        MatSet0(z_sa_cat, r, c*i, v);
      }
    }
  }

#if DEBUG
  MatWrite(stdout , z_sa_cat, fmt,"z_sa_cat:");
  MatWrite(stdout , trfm->w_o, fmt,"w_o:");
#endif

  MatMul(z, z_sa_cat, trfm->w_o);

#if DEBUG
  MatWrite(stdout , z, fmt,"z:");
#endif

  // free mat
  MatUnDim(z_sa_cat);
  for (int i = 0; i < N_MHA; i++) {
    MatUnDim(z_sa[i]);
  }

  return z;
}

Tensor transformer_block(struct Transformer *trfm, Tensor x) {

  Tensor z;
  z = multi_head_attention(trfm, x);

  Tensor residual;
#if DEBUG
  MatWrite(stdout , x, fmt,"x:");
  MatWrite(stdout , z, fmt,"z:");
#endif
  MatAdd(z, z, x);

  // TODO: add layer norm
  // TODO: add FFN (feed forward neural network)

#if 0
  Tensor x_atten, x_ff1, x_ff2;

  // XXX
//  x_atten = layer_norm(residual);

  // XXX
//  x_ff1 = layer_ff(trfm->ff1, x_atten);

//  x_ff2 = layer_ff(trfm->ff2, x_ff1);


  Tensor residual2, x_out;
  MatAdd(x_atten, x_ff2, residual2);

//  x_out = layer_norm(residual2);

  // skip dropout

  return x_out;
#endif

  return z;
}

// instantiate transformer block pass the config struct
struct Transformer *init_transformer(void) {

  struct Transformer *trfm = malloc(sizeof(struct Transformer));
  memset(trfm, 0, sizeof(struct Transformer));
  trfm->dim_embedding = DIM_EMBEDDING;
  trfm->dim_internal = DIM_INTERNAL;
  trfm->n_seq = N_SEQ;
  trfm->n_mha = N_MHA;

  trfm->w_o = MatDim(trfm->dim_internal * trfm->n_mha, trfm->dim_embedding);
  // DEBUG FILL
  MatFill(trfm->w_o, 10.);

  struct SelfAtten *sa;
  for (int i = 0; i < N_MHA; i++) {
    sa = &trfm->sa[i];
    sa->trfm = trfm;  // set backpointer
    sa->w_q = MatDim(trfm->dim_embedding, trfm->dim_internal);
    sa->w_k = MatDim(trfm->dim_embedding, trfm->dim_internal);
    sa->w_v = MatDim(trfm->dim_embedding, trfm->dim_internal);

    // set values
    // DEBUG FILL
    MatFill(sa->w_q, 1.);
    MatFill(sa->w_k, 2.);
    MatFill(sa->w_v, 3.);
  }

  return trfm;
}

// load weight parameters by pre-trained model weight
int load_params(void) {

  return 0;
}


int test_attention(void) {

  printf("===================================\n");
  printf("Testing single self attention block\n");
  printf("===================================\n");

  struct Transformer *trfm = init_transformer();
  assert(trfm);

  Tensor x, z;
  // test input x
  x = MatDim(trfm->n_seq, trfm->dim_embedding);
  // set values, use MatSet
  // DEBUG FILL
  MatFill(x, 4);

  z = self_attention(&trfm->sa[0], x);

#if DEBUG
  MatWrite(stdout , z, fmt,"z:");
#endif

  MatUnDim(x);
  MatUnDim(z);

  // TODO: free trfm
  return 0;
}


int test_multiheadattention(void) {

  struct Transformer *trfm = init_transformer();
  assert(trfm);

  printf("===================================\n");
  printf("Testing multi head attention block\n");
  printf("===================================\n");

  Tensor x, z;
  // test input x
  x = MatDim(trfm->n_seq, trfm->dim_embedding);
  // set values, use MatSet
  // DEBUG FILL
  MatFill(x, 4);

  z = multi_head_attention(trfm, x);

#if DEBUG
  MatWrite(stdout , z, fmt,"z:");
#endif

  return 0;
}

int test_transformer_block(void) {

  printf("===================================\n");
  printf("Testing single transformer block\n");
  printf("===================================\n");

  struct Transformer *trfm = init_transformer();
  assert(trfm);

  Tensor x, z;
  // test input x
  x = MatDim(trfm->n_seq, trfm->dim_embedding);
  // set values, use MatSet
  // DEBUG FILL
  MatFill(x, 4);

  z = transformer_block(trfm, x);

#if DEBUG
  MatWrite(stdout , z, fmt,"z:");
#endif
  return 0;
}

int main(int argc, char **argv) {

  printf("Transformer simple\n");
//  test_attention();

//  test_multiheadattention();
  test_transformer_block();

  return 0;
}
