import torch
from torch import nn
import torch.nn.functional as F

import random, math

from util import d, here

class SelfAttention(nn.Module):
    def __init__(self, dim_emb, dim_internal, heads=8, mask=False, dropout=0.0):
        """
        A single self attention block

        :param dim_emb: embedding dimension
        :param dim_internal: dimension of internal representation, usually the same as dim_emb
        :param head: number of multi head
        :param mask

        """
        super().__init__()

        self.dim_emb = dim_emb
        self.dim_internal = dim_internal
        self.heads = heads
        self.mask = mask

        self.toqueries = nn.Linear(dim_emb, dim_internal)
        self.tokeys = nn.Linear(dim_emb, dim_internal)
        self.tovalues = nn.Linear(dim_emb, dim_internal)

        self.kSqrt_dim_emb = math.sqrt(self.dim_emb)

    def forward(self, x):
        # [batch size, seq, embedding]

        b, t, e = x.size()

        assert e == self.dim_emb, f'Input embedding ({e}) should match the layer embedding ({self.dim_emb})'

        queries = self.toqueries(x)
        keys = self.tokeys(x)
        values = self.tovalues(x)

        keys_transposed =  keys.transpose(-2, -1)
        dot = torch.matmul(queries, keys_transposed)  / self.kSqrt_dim_emb

        # softmax on row-wise element
        p_attn = F.softmax(dot, dim=2)

        z = torch.matmul(p_attn, values)

        return z


class MultiHeadAttention(nn.Module):
    def __init__(self, n_seq, dim_emb, dim_internal, heads=8, mask=False, dropout=0.0):
        """
        multi head attention block

        :param n_seq: number of token sequence
        :param dim_emb: embedding dimension
        :param dim_internal: dimension of internal representation, usually the same as dim_emb
        :param head: number of multi head
        :param mask

        """
        super().__init__()

        self.n_seq = n_seq
        self.dim_emb = dim_emb
        self.heads = heads
        self.mask = mask

        self.attentions = nn.ModuleList([SelfAttention(dim_emb, dim_internal, heads, mask, dropout) \
                                         for _ in range(heads)])
        self.w_o = nn.ModuleList([nn.Linear(dim_internal, dim_emb) \
                                         for _ in range(heads)])

    def forward(self, x):

        output = torch.zeros_like(x)

        for attention, w_o in zip(self.attentions, self.w_o):
            z = attention(x)
            output += w_o(z)

        return output


class TransformerBlock(nn.Module):
    def __init__(self, n_seq, dim_emb, dim_internal, heads=8, mask=False, ff_hidden_mult=4, dropout=0.0):
        """
        :ff_hidden_mult: number of multiples of embedding for total hidden size
        """
        super().__init__()

        self.mha = MultiHeadAttention(n_seq=n_seq, dim_emb=dim_emb, dim_internal=dim_internal, heads=heads, mask=mask, dropout=dropout)
        self.mask = mask

        self.norm1 = nn.LayerNorm(dim_emb)
        self.norm2 = nn.LayerNorm(dim_emb)

        self.ff = nn.Sequential(
            nn.Linear(dim_emb, ff_hidden_mult * dim_emb),
            nn.ReLU(),
            nn.Linear(ff_hidden_mult * dim_emb, dim_emb)
        )
        self.do = nn.Dropout(dropout)


    def forward(self, x):
        """
        """

        z = self.mha(x)

        # layer norm
        out = self.norm1(z + x)
        out = self.do(out)

        ff = self.ff(out)

        out = self.norm2(ff + out)
        out = self.do(out)


        return out

## for classify
class TransformerSimpleClassify(nn.Module):
    def __init__(self, n_seq, dim_emb, dim_internal, num_tokens, num_classes, max_pool=True, heads=8, depth=6, mask=False, ff_hidden_mult=4, dropout=0.0):
        super().__init__()

        self.num_tokens = num_tokens
        self.max_pool = max_pool
        self.depth = depth
        self.dim_emb = dim_emb
        self.dim_internal = dim_internal

        self.token_embedding = nn.Embedding(embedding_dim=dim_emb, num_embeddings=num_tokens)
        self.pos_embedding = nn.Embedding(embedding_dim=dim_emb, num_embeddings=n_seq)

        trfm_blocks = [TransformerBlock(n_seq=n_seq, dim_emb=dim_emb, dim_internal=dim_emb, heads=heads) \
                       for _ in range(depth)]

        self.trfm_blocks = nn.Sequential(*trfm_blocks)

        self.toprobs = nn.Linear(dim_emb, num_classes)

        self.do = nn.Dropout(dropout)

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


if __name__ == "__main__":

    """
    finish up the tf block
    put them together, and be able to train
    """

    # create model
    #model = transformer()

    sa = SelfAttention(4, 3)
    print(sa)
    mha = MultiHeadAttention(2, 4, 3)
    print(mha)
    model = TransformerBlock(2, 4, 3)
    print(model)

    model_tb = TransformerSimpleClassify(2, 4, 4, 10, 2)
    print(model_tb)

    """
    opt = torch.optim.Adam(lr=arg.lr, params=model.parameters())
    sch = torch.optim.lr_scheduler.LambdaLR(opt, lambda i: min(i / (arg.lr_warmup / arg.batch_size), 1.0))

    """

    n_batch  = 1
    x = torch.ones([n_batch, 2,4])
    # forward path
    model.train(False)

    out = model(x)
    print(out)
