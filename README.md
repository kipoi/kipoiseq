# kipoiseq
<a href='https://circleci.com/gh/kipoi/kipoiseq'>
	<img alt='CircleCI' src='https://circleci.com/gh/kipoi/kipoiseq.svg?style=svg' style="max-height:20px;width:auto">
</a>
<a href=https://coveralls.io/github/kipoi/kipoiseq?branch=master>
	<img alt='Coverage status' src=https://coveralls.io/repos/github/kipoi/kipoiseq/badge.svg?branch=master style="max-height:20px;width:auto;">
</a>

Standard set of data-loaders for training and making predictions for DNA sequence-based models.

All dataloaders in `kipoiseq.dataloaders` decorated with `@kipoi_dataloader` (SeqIntervalDl and StringSeqIntervalDl) are compatible Kipoi models and can be directly used when specifying a new model in `model.yaml`:
```yaml
...
default_dataloader:
  defined_as: kipoiseq.dataloaders.SeqIntervalDl
  default_args:
    auto_resize_len: 1000 # override default args in SeqIntervalDl
    
dependencies:
  pip:
    - kipoiseq
...
```

## Installation

```bash
pip install kipoiseq
```

Optional dependencies:
```bash
pip install cyvcf2, pyranges
conda install cyvcf2, pyranges
```

## Getting started

```python
from kipoiseq.dataloaders import SeqIntervalDl

dl = SeqIntervalDl.init_example()  # use the provided example files
# your own files
dl = SeqIntervalDl("intervals.bed", "genome.fa")

len(dl)  # length of the dataset

dl[0]  # get one instance. # returns a dictionary: 
# dict(inputs=<one-hot-encoded-array>, 
#      targets=<additional columns in the bed file>, 
#      metadata=dict(ranges=GenomicRanges(chr=, start, end)...

all = dl.load_all()  # load the whole dataset

# load batches of data
it = dl.batch_iter(32, num_workers=8)  # load batches of data in parallel using 8 workers
# returns a dictionary with all three keys: inputs, targets, metadata

it = dl.batch_train_iter(32, num_workers=8)
# returns a tuple: (inputs, targets), can be used directly with keras' `model.fit_generator`
```

More info:
- Follow the getting-started [colab notebook](https://colab.research.google.com/github/kipoi/kipoiseq/blob/master/notebooks/getting-started.ipynb).
- See [docs](https://kipoi.org/kipoiseq/)

## How to write your own data-loaders
- Read the pytorch [Data Loading and Processing Tutorial](https://pytorch.org/tutorials/beginner/data_loading_tutorial.html) to become more familiar with transforms and dataloaders
- Read the code for `SeqIntervalDl` in [kipoiseq/dataloaders/sequence.py](https://github.com/kipoi/kipoiseq/blob/master/kipoiseq/dataloaders/sequence.py)
  - you can skip the `@kipoi_dataloader` and the long yaml doc-string. These are only required if you want to use dataloaders in Kipoi's model.yaml files.
- Explore the available transforms ([functional](http://kipoi.org/kipoiseq/transforms/functional/), [class-based](http://kipoi.org/kipoiseq/transforms/transforms/)) or extractors ([kipoiseq](https://github.com/kipoi/kipoiseq/blob/master/kipoiseq/extractors.py))
