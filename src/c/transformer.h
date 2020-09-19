#include "sml.h"

typedef MATRIX Tensor;

struct Transformer;

/*
 * self attention block
 */
struct SelfAtten {
  Tensor w_q; // dim_embedding x dim_internal
  Tensor w_k; // ..
  Tensor w_v; // ..

  struct Transformer *trfm; // back pointer
};

// number of multi head attention
#define N_MHA 8

struct Transformer {
  //int n_mha; // number of multi head attention
  struct SelfAtten sa[N_MHA]; // self attention

  Tensor w_o; // weight for converting concatenated multi headed attention

  Tensor w_ff1; // weight for ffn1 (feed forward network)
  Tensor w_ff2;

  int dim_embedding; // token embedding dimension
  int dim_internal; // internal dimension of the token
  int n_seq; // number token, i.e. sequence length
  int n_mha; // number of multi-head attention
  int dim_hidden; // dim of hidden layer of ffn
};
