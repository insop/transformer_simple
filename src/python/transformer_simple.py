import torch
from torch import nn
import torch.nn.functional as F

import random, math

class SelfAttention(nn.Module):
    def __init__(self, emb, dim_internal, heads=8, mask=False, dropout=0.0):
        """
        :param emb: embedding dimension
        :param dim_internal: dimension of internal representation
        :param head: number of multi head
        :param mask

        """
        super().__init__()

        self.emb = emb
        self.dim_internal = dim_internal
        self.heads = heads
        self.mask = mask

        self.toqueries = nn.Linear(emb, dim_internal)
        self.tokeys = nn.Linear(emb, dim_internal)
        self.tovalues = nn.Linear(emb, dim_internal)

    def forward(self, x):
        # single batch

        # seq, embedding
        t, e = x.size()

        assert e == self.emb, f'Input embedding ({e}) should match the layer embedding ({self.emb})'

        queries = self.toqueries(x)
        keys = self.tokeys(x)
        values = self.tovalues(x)

        dot = torch.matmul(queries, keys.transpose(0, 1)) / math.sqrt(self.emb)

        # softmax on row-wise element
        p_attn = F.softmax(dot, dim=1)

        z = torch.matmul(p_attn, values)

        return z



class MultiHeadAttention(nn.Module):
    def __init__(self, n_seq, emb, dim_internal, heads=8, mask=False, dropout=0.0):
        """
        :param emb: embedding dimension
        :param dim_internal: dimension of internal representation
        :param head: number of multi head
        :param mask

        """
        super().__init__()

        self.n_seq = n_seq
        self.emb = emb
        self.heads = heads
        self.mask = mask

        self.toqueries = nn.Linear(emb, dim_internal)
        self.tokeys = nn.Linear(emb, dim_internal)

        self.attentions = nn.ModuleList([SelfAttention(emb, dim_internal, heads, mask, dropout) \
                                         for _ in range(heads)])
        self.w_o = nn.ModuleList([nn.Linear(dim_internal, emb) \
                                         for _ in range(heads)])

        self.layer_norm = nn.LayerNorm(self.n_seq, eps=1e-6)

    def forward(self, x):

        output = torch.zeros_like(x)

        for attention, w_o in zip(self.attentions, self.w_o):
            z = attention(x)
            output += w_o(z)

        #output = F.dropout(output, p=dropout)
        # residual
        #output = output + x

        #output = self.layer_norm(output)

        return output


class TransformerBlock(nn.Module):
    def __init__(self, n_seq, emb, dim_internal, heads=8, mask=False, ff_hidden_mult=4, dropout=0.0):
        """
        """
        super().__init__()

        self.mha = MultiHeadAttention(n_seq=n_seq, emb=emb, dim_internal=dim_internal, heads=heads, mask=mask, dropout=dropout)
        self.mask = mask

        self.norm1 = nn.LayerNorm(emb)
        self.norm2 = nn.LayerNorm(emb)

        self.ff = nn.Sequential(
            nn.Linear(emb, ff_hidden_mult * emb),
            nn.ReLU(),
            nn.Linear(ff_hidden_mult * emb, emb)
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

    """
    opt = torch.optim.Adam(lr=arg.lr, params=model.parameters())
    sch = torch.optim.lr_scheduler.LambdaLR(opt, lambda i: min(i / (arg.lr_warmup / arg.batch_size), 1.0))

    """

    x = torch.ones([2,4])
    # forward path
    model.train(False)

    out = model(x)
    print(out)
