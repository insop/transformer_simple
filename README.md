# Implement Transformer (TODO reference) with C

Implement Tranformer with C and Python for educational purpose. Code is written for readability.

# Build

```
cd lib/sml
make
cd ..

cd src/c
make

```

# Run

## C version
```
cd src/c
./transformer

```

## Python version
```
cd src/python

python ./experiments/classify.py  --random-seed=1234 --num-epochs=1 --tiny

```

# TODO:
- main Makefile for build library and c executable
- add config load
- add trained weight load
- add python code to generate test vector, and use that to test the C code

# Reference:
- SML: small math library, http://www.bios.unc.edu/distrib/bios235/sml/
- Transformer tutorial: http://jalammar.github.io/illustrated-transformer/
- Python Transformer implementation: http://peterbloem.nl/blog/transformers
