"""
classifier using transformer

"""
import torch
from torch import nn
import torch.nn.functional as F

import random, math

from util import d, here

from transformer_simple import TransformerBlock_naive


## for classify
class TransformerSimpleClassify(nn.Module):
    def __init__(self, n_seq, dim_emb, dim_internal, num_tokens, num_classes, max_pool=True, heads=8, depth=6, mask=False, ff_hidden_mult=4, dropout=0.0, dtype=torch.float32):
        super().__init__()

        self.num_tokens = num_tokens
        self.max_pool = max_pool
        self.depth = depth
        self.dim_emb = dim_emb
        self.dim_internal = dim_internal

        self.token_embedding = nn.Embedding(embedding_dim=dim_emb, num_embeddings=num_tokens).type(dtype)
        self.pos_embedding = nn.Embedding(embedding_dim=dim_emb, num_embeddings=n_seq).type(dtype)

        trfm_blocks = [TransformerBlock_naive(n_seq=n_seq, dim_emb=dim_emb, dim_internal=dim_emb, heads=heads, dtype=dtype) \
                       for _ in range(depth)]

        self.trfm_blocks = nn.Sequential(*trfm_blocks).type(dtype)

        self.toprobs = nn.Linear(dim_emb, num_classes).type(dtype)

        self.do = nn.Dropout(dropout).type(dtype)

    def forward(self, x):

        # x = [batch, n_seq]
        tokens = self.token_embedding(x)

        b, t, e = tokens.size()

        positions = self.pos_embedding(torch.arange(t, device=d()))[None, :, :].expand(b, t, e)
        x = tokens + positions

        x = self.do(x)

        x = self.trfm_blocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool in sequence direction

        x = self.toprobs(x)

        return F.log_softmax(x, dim=1)


