"""
bert base quantized model

"""

import argparse, json, os, requests, sys, time
from io import BytesIO
from os.path import join, isfile
from PIL import Image

import timeit

from mxnet.gluon.model_zoo import vision
import numpy as np

import tvm
from tvm import te
from tvm import rpc, autotvm, relay
from tvm.contrib import graph_runtime, utils, download
from tvm.contrib.debugger import debug_runtime
from tvm.relay import transform

# Make sure that TVM was compiled with RPC=1
assert tvm.runtime.enabled("rpc")

import logging
import numpy as np
import os
import random
import sys
import time
import torch

from argparse import Namespace
from torch.utils.data import (DataLoader, RandomSampler, SequentialSampler,
                              TensorDataset)
from tqdm import tqdm
from transformers import (BertConfig, BertForSequenceClassification, BertTokenizer,)
from transformers import glue_compute_metrics as compute_metrics
from transformers import glue_output_modes as output_modes
from transformers import glue_processors as processors
from transformers import glue_convert_examples_to_features as convert_examples_to_features

import logging
import numpy as np
import os
import random
import sys
import time
import torch

import transformers

from transformers import BertModel, BertTokenizer, BertConfig
import numpy

import torch

# Setup logging
logger = logging.getLogger(__name__)
logging.basicConfig(format = '%(asctime)s - %(levelname)s - %(name)s -   %(message)s',
                    datefmt = '%m/%d/%Y %H:%M:%S',
                    level = logging.WARN)

logging.getLogger("transformers.modeling_utils").setLevel(
   logging.WARN)  # Reduce logging

print(torch.__version__)


# -------- defines --------------

TIMEIT_MODEL = False
TIMEIT_NRUN = 10
N_TRIAL = 2000
BERT_SIZE = 'large'


"""MRPC data download

pwd
ls
wget https://gist.githubusercontent.com/W4ngatang/60c2bdb54d156a41194446737ce03e2e/raw/17b8dd0d724281ed7c3b2aeeda662b92809aadd5/download_glue_data.py
python download_glue_data.py --data_dir='glue_data' --tasks='MRPC'
ls glue_data/MRPC

wget https://download.pytorch.org/tutorial/MRPC.zip
unzip MRPC.zip
ls
pwd
"""

configs = Namespace()

# The output directory for the fine-tuned model.
# configs.output_dir = "/content/MRPC/"
configs.output_dir = "./MRPC/"

# The data directory for the MRPC task in the GLUE benchmark.
# configs.data_dir = "/content/glue_data/MRPC"
configs.data_dir = "./glue_data/MRPC"

# The model name or path for the pre-trained model.
configs.model_name_or_path = "bert-%s-uncased" % BERT_SIZE
log_tune_file = "bert-%s-quantized-tuning" % BERT_SIZE
# The maximum length of an input sequence
configs.max_seq_length = 128

# Prepare GLUE task.
configs.task_name = "MRPC".lower()
configs.processor = processors[configs.task_name]()
configs.output_mode = output_modes[configs.task_name]
configs.label_list = configs.processor.get_labels()
configs.model_type = "bert".lower()
configs.do_lower_case = True

# Set the device, batch size, topology, and caching flags.
configs.device = "cpu"
configs.per_gpu_eval_batch_size = 8
configs.n_gpu = 0
configs.local_rank = -1
configs.overwrite_cache = False


# Set random seed for reproducibility.
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
set_seed(42)


#model_org = BertForSequenceClassification.from_pretrained(configs.output_dir) # model in the output_dir is 'base'
model_org = BertForSequenceClassification.from_pretrained(configs.model_name_or_path)
model_org.to(configs.device)


enc = BertTokenizer.from_pretrained(configs.model_name_or_path)

# Tokenizing input text
text = "[CLS] Who was Jim Henson ? [SEP] Jim Henson was a puppeteer [SEP]"
tokenized_text = enc.tokenize(text)

# Masking one of the input tokens
masked_index = 8
tokenized_text[masked_index] = '[MASK]'
indexed_tokens = enc.convert_tokens_to_ids(tokenized_text)
segments_ids = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]

# Creating a dummy input
tokens_tensor = torch.tensor([indexed_tokens])
segments_tensors = torch.tensor([segments_ids])
dummy_input = [tokens_tensor, segments_tensors]

model_org.eval()
for p in model_org.parameters():
    p.requires_grad_(False)

transformers.__version__

# Creating the trace
traced_model_org = torch.jit.trace(model_org, [tokens_tensor, segments_tensors])
traced_model_org.eval()
for p in traced_model_org.parameters():
    p.requires_grad_(False)

print("traced_mode", traced_model_org)


tt_c = tokens_tensor.cpu()
st_c = segments_tensors.cpu()
res_pt = model_org(tt_c, st_c)


print("#"*100)
print("org model run")
print("#"*100)


#if TIMEIT_MODEL:
ti = timeit.timeit("model_org(tt_c, st_c)", globals=globals(), number=TIMEIT_NRUN)
print("Original pytorch bert run:", ti/TIMEIT_NRUN)

# ===========================
print("#"*100)
print("shape_list")
print("#"*100)

shape_list = [(i.debugName().split('.')[0], i.type().sizes()) for i in  list(traced_model_org.graph.inputs())[1:]]
#shape_list


"""

mod_bert, params_bert = tvm.relay.frontend.pytorch.from_pytorch(traced_model_org,
                        shape_list, default_dtype="float32")

print("mod_bert, params_bert", mod_bert, params_bert)


# ===========================


target = 'llvm'
ctx = tvm.cpu(0)


target_host = 'llvm'


tt_a = tvm.nd.array(tokens_tensor.numpy(), ctx)
st_a = tvm.nd.array(segments_tensors.numpy(), ctx)


tvm.relay.backend.compile_engine.get().clear() # just to be sure, see https://github.com/apache/incubator-tvm/pull/5724

with tvm.transform.PassContext(opt_level=3):
        graph, lib, params = tvm.relay.build(mod_bert,
                                     target=target,
                                     target_host=target_host,
                                     params=params_bert)
module_tvm = tvm.contrib.graph_runtime.create(graph, lib, ctx)

print("#"*100)
print("module_tvm is created")
print("#"*100)

module_tvm.set_input("input_ids", tt_a)
module_tvm.set_input("attention_mask", st_a)
module_tvm.set_input(**params)
# module_tvm.run()
o0 = module_tvm.get_output(0)
diff = numpy.abs((res_pt[0].cpu().numpy() - o0.asnumpy())).max()

print("#"*100)
print("Diff btwen mod_bert and tvm_module", diff)
print("#"*100)



print("#"*100)
print("model_tvm run")
print("#"*100)

if TIMEIT_MODEL:
    print("TVM module_tvm bert base:", timeit.timeit(model_tvm.run(), globals=globals()))


"""

# load quantized model

print("#"*100)
print("load quantized model")
print("#"*100)


quantized_model = torch.quantization.quantize_dynamic(
    model_org, {torch.nn.Linear}, dtype=torch.qint8
)

ti = timeit.timeit("quantized_model(tt_c, st_c)", globals=globals(), number=TIMEIT_NRUN)
print("Pytorch quantized bert run:", ti/TIMEIT_NRUN)

print("#"*100)
print("quantized model", quantized_model)
print("#"*100)

# Creating the trace
traced_quantized_model = torch.jit.trace(quantized_model, [tokens_tensor, segments_tensors])
traced_quantized_model.eval()
for p in traced_quantized_model.parameters():
    p.requires_grad_(False)

print("#"*100)
print("traced_quantized_model", traced_quantized_model)

q_mod_bert, q_params_bert = tvm.relay.frontend.pytorch.from_pytorch(traced_quantized_model,
                        shape_list, default_dtype="int8")

print("#"*100)
print("q_mod_bert,", q_mod_bert)

# model()
tt_c = tokens_tensor.cpu()
st_c = segments_tensors.cpu()
res_pt = quantized_model(tt_c, st_c)


print("quantized model")

#quantized_model(tt_c, st_c)


target = 'llvm'
ctx = tvm.cpu(0)

target_host = 'llvm'

tt_a = tvm.nd.array(tokens_tensor.numpy(), ctx)
st_a = tvm.nd.array(segments_tensors.numpy(), ctx)

tvm.relay.backend.compile_engine.get().clear() # just to be sure, see https://github.com/apache/incubator-tvm/pull/5724

with tvm.transform.PassContext(opt_level=3):
        q_graph, q_lib, q_params = tvm.relay.build(q_mod_bert,
                                     target=target,
                                     target_host=target_host,
                                     params=q_params_bert)

# TODO: save this image and load to test

q_module_tvm = tvm.contrib.graph_runtime.create(q_graph, q_lib, ctx)

q_module_tvm.set_input("input_ids", tt_a)
q_module_tvm.set_input("attention_mask", st_a)
q_module_tvm.set_input(**q_params_bert)

print("#"*100)
print("Start: quantized bert base tvm run")
print("#"*100)

q_module_tvm.run()
o0 = q_module_tvm.get_output(0)
diff = numpy.abs((res_pt[0].cpu().numpy() - o0.asnumpy())).max()


print("#"*100)
print("Done: quantized bert base tvm run")
print("#"*100)

ti = timeit.timeit("q_module_tvm.run()", globals=globals(), number=TIMEIT_NRUN)
print("TVM quantized bert run:", ti/TIMEIT_NRUN)

#import pdb;pdb.set_trace()

# save ISS
"""


# https://colab.research.google.com/github/d2l-ai/d2l-tvm-colab/blob/master/chapter_getting_started/from_mxnet.ipynb#scrollTo=AAMNEtJDb4LS

name = './bert-large-uncase-quantized-tvm'
graph_fn, mod_fn, params_fn = [name+ext for ext in ('.json','.tar','.params')]
q_lib.export_library(mod_fn)
with open(graph_fn, 'w') as f:
    f.write(q_graph)
with open(params_fn, 'wb') as f:
    f.write(tvm.relay.save_param_dict(q_params))


# load ISS

if True:
    name = './bert-large-uncase-quantized'
    graph_ld, mod_ld, params_ld = [name+ext for ext in ('.json','.tar','.params')]

    loaded_graph = open(graph_ld).read()
    loaded_mod = tvm.runtime.load_module(mod_ld)
    loaded_params = open(params_ld, "rb").read()

    module_loaded = tvm.contrib.graph_runtime.create(loaded_graph, loaded_mod, ctx)
    module_loaded.load_params(loaded_params)

# test loaded module
    module_loaded.set_input("input_ids", tt_a)
    module_loaded.set_input("attention_mask", st_a)
# module_loaded.set_input(**loaded_params)
    module_loaded.run()
    o0_loaded = module_loaded.get_output(0)
    diff = numpy.abs((res_pt[0].cpu().numpy() - o0_loaded.asnumpy())).max()
    print("loaded model diff", diff)
    ti = timeit.timeit("module_loaded.run()", globals=globals(), number=TIMEIT_NRUN)
    print("loaded model TVM quantized bert run:", ti/TIMEIT_NRUN)
"""
