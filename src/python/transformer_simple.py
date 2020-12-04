"""

modules for implementing basic transformer block

functions with `_naive` prefix means that it is a naive implementation


"""


import torch
from torch import nn
import torch.nn.functional as F

import random, math

from util import d, here

def init_weight_zavier(x):
    nn.init.xavier_uniform_(x.weight)
    if x.bias is not None:
        nn.init.constant_(x.bias, 0)

class SelfAttention_naive(nn.Module):
    def __init__(self, dim_emb, dim_internal, heads=8, mask=False, dropout=0.0, dtype=torch.float32):
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


        self.toqueries = nn.Linear(dim_emb, dim_internal).type(dtype)
        self.tokeys = nn.Linear(dim_emb, dim_internal).type(dtype)
        self.tovalues = nn.Linear(dim_emb, dim_internal).type(dtype)

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


class MultiHeadAttention_naive(nn.Module):
    def __init__(self, n_seq, dim_emb, dim_internal, heads=8, mask=False, dropout=0.0, dtype=torch.float32):
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

        self.attentions = nn.ModuleList([SelfAttention_naive(dim_emb, dim_internal, heads, mask, dropout, dtype=dtype) \
                                         for _ in range(heads)])
        self.w_o = nn.ModuleList([nn.Linear(dim_internal, dim_emb).type(dtype) \
                                         for _ in range(heads)])

    def forward(self, x):

        output = torch.zeros_like(x)

        for attention, w_o in zip(self.attentions, self.w_o):
            z = attention(x)
            output += w_o(z)

        return output


class TransformerBlock_naive(nn.Module):
    def __init__(self, n_seq, dim_emb, dim_internal, heads=8, mask=False, ff_hidden_mult=4, dropout=0.0, dtype=torch.float32):
        """
        :ff_hidden_mult: number of multiples of embedding for total hidden size
        """
        super().__init__()

        self.mha = MultiHeadAttention_naive(n_seq=n_seq, dim_emb=dim_emb, dim_internal=dim_internal, heads=heads, mask=mask, dropout=dropout, dtype=dtype)
        self.mask = mask

        self.norm1 = nn.LayerNorm(dim_emb).type(dtype)
        self.norm2 = nn.LayerNorm(dim_emb).type(dtype)

        self.ff = nn.Sequential(
            nn.Linear(dim_emb, ff_hidden_mult * dim_emb).type(dtype),
            nn.ReLU().type(dtype),
            nn.Linear(ff_hidden_mult * dim_emb, dim_emb).type(dtype)).type(dtype)
        self.do = nn.Dropout(dropout).type(dtype)

        if dtype == torch.float32 or dtype == torch.float64:
            init_weight_zavier(self.ff[0])  # 1st linear
            init_weight_zavier(self.ff[2])  # 2nd linear


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

if __name__ == "__main__":

    from classifier import TransformerSimpleClassify

    """
    finish up the tf block
    put them together, and be able to train
    """
    dtype = torch.float32

    # create model
    #model = transformer()

    sa = SelfAttention_naive(4, 3, dtype=dtype)
    print(sa)
    mha = MultiHeadAttention_naive(2, 4, 3, dtype=dtype)
    print(mha)
    model = TransformerBlock_naive(2, 4, 3, dtype=dtype)
    print(model)

    model_tb = TransformerSimpleClassify(2, 4, 4, 10, 2, dtype=dtype)
    print(model_tb)

    """
    opt = torch.optim.Adam(lr=arg.lr, params=model.parameters())
    sch = torch.optim.lr_scheduler.LambdaLR(opt, lambda i: min(i / (arg.lr_warmup / arg.batch_size), 1.0))

    """

    n_batch  = 1
    x = torch.ones([n_batch, 2,4], dtype=dtype)
    # forward path
    model.train(False)

    out = model(x)
    print(out)
