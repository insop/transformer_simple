import torch
import torch.nn as nn
import torch.nn.functional as F
import unittest
from transformer_simple import SelfAttention_naive, MultiHeadAttention_naive, TransformerBlock_naive

class TestSelfAttentionNaive(unittest.TestCase):
    def setUp(self):
        self.dim_emb = 4
        self.dim_internal = 3
        self.heads = 8
        self.n_seq = 2
        self.dtype = torch.float32
        self.model = SelfAttention_naive(self.dim_emb, self.dim_internal, self.heads, dtype=self.dtype)
        self.x = torch.ones([1, self.n_seq, self.dim_emb], dtype=self.dtype)

    def test_forward(self):
        output = self.model(self.x)
        self.assertEqual(output.shape, (1, self.n_seq, self.dim_internal))

class TestMultiHeadAttentionNaive(unittest.TestCase):
    def setUp(self):
        self.dim_emb = 4
        self.dim_internal = 3
        self.heads = 8
        self.n_seq = 2
        self.dtype = torch.float32
        self.model = MultiHeadAttention_naive(self.n_seq, self.dim_emb, self.dim_internal, self.heads, dtype=self.dtype)
        self.x = torch.ones([1, self.n_seq, self.dim_emb], dtype=self.dtype)

    def test_forward(self):
        output = self.model(self.x)
        self.assertEqual(output.shape, (1, self.n_seq, self.dim_emb))

class TestTransformerBlockNaive(unittest.TestCase):
    def setUp(self):
        self.dim_emb = 4
        self.dim_internal = 3
        self.heads = 8
        self.n_seq = 2
        self.dtype = torch.float32
        self.model = TransformerBlock_naive(self.n_seq, self.dim_emb, self.dim_internal, self.heads, dtype=self.dtype)
        self.x = torch.ones([1, self.n_seq, self.dim_emb], dtype=self.dtype)

    def test_forward(self):
        output = self.model(self.x)
        self.assertEqual(output.shape, (1, self.n_seq, self.dim_emb))

if __name__ == '__main__':
    unittest.main()
