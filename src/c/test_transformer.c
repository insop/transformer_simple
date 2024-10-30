#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "transformer.h"

void test_self_attention() {
    struct Transformer *trfm = init_transformer();
    assert(trfm);

    Tensor x, z;
    x = MatDim(trfm->n_seq, trfm->dim_embedding);
    MatFill(x, 4);

    z = self_attention(&trfm->sa[0], x);

    assert(MatRows(z) == trfm->n_seq);
    assert(MatCols(z) == trfm->dim_internal);

    MatUnDim(x);
    MatUnDim(z);
}

void test_multi_head_attention() {
    struct Transformer *trfm = init_transformer();
    assert(trfm);

    Tensor x, z;
    x = MatDim(trfm->n_seq, trfm->dim_embedding);
    MatFill(x, 4);

    z = multi_head_attention(trfm, x);

    assert(MatRows(z) == trfm->n_seq);
    assert(MatCols(z) == trfm->dim_embedding);

    MatUnDim(x);
    MatUnDim(z);
}

void test_transformer_block() {
    struct Transformer *trfm = init_transformer();
    assert(trfm);

    Tensor x, z;
    x = MatDim(trfm->n_seq, trfm->dim_embedding);
    MatFill(x, 4);

    z = transformer_block(trfm, x);

    assert(MatRows(z) == trfm->n_seq);
    assert(MatCols(z) == trfm->dim_embedding);

    MatUnDim(x);
    MatUnDim(z);
}

int main() {
    test_self_attention();
    test_multi_head_attention();
    test_transformer_block();

    printf("All tests passed.\n");
    return 0;
}
