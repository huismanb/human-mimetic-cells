#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# %%
import copy
import gc
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback
from typing import List, Tuple, Dict, Union, Optional
import warnings
import pandas as pd
# from . import asyn
import pickle
import torch
from anndata import AnnData
import scanpy as sc
import scvi
import seaborn as sns
import numpy as np
import wandb
from scipy.sparse import issparse
import matplotlib.pyplot as plt
from torch import nn
from torch.nn import functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)
from sklearn.metrics import confusion_matrix

sys.path.insert(0, "../")
import scgpt as scg
from scgpt.model import TransformerModel, AdversarialDiscriminator
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.loss import (
    masked_mse_loss,
    masked_relative_error,
    criterion_neg_log_bernoulli,
)
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor
from scgpt import SubsetsBatchSampler
from scgpt.utils import set_seed, category_str2int, eval_scib_metrics

sc.set_figure_params(figsize=(6, 6))
os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings('ignore')


# In[ ]:
 
os.environ['WANDB_API_KEY']='your_key
os.environ['WANDB_ENTITY']='your_name' 

os.chdir("/n/groups/cbdm_lab/brh040/analysis/scGPT/mimetics-only_version")


# In[ ]:


## Parameters can adjust:
# dataset_name, epochs, lr, batch_size, dropout, schedule_ratio, save_eval_interval

hyperparameter_defaults = dict(
    seed=0,
    dataset_name="zebrafish_mtec",
    do_train=True,
    load_model='/n/groups/cbdm_lab/brh040/analysis/scGPT/mimetics-only_version/save_human_scGPT/dev_human_mtec-Mar09-22-17_mimetics-only',
    mask_ratio=0.0,
    epochs=10,
    n_bins=51,
    MVC=False, # Masked value prediction for cell embedding
    ecs_thres=0.0, # Elastic cell similarity objective, 0.0 to 1.0, 0.0 to disable
    dab_weight=0.0,
    lr=1e-4,
    batch_size=4, 
    layer_size=128,
    nlayers=4,  # number of nn.TransformerEncoderLayer in nn.TransformerEncoder
    nhead=4,  # number of heads in nn.MultiheadAttention
    dropout=0.2,  # dropout probability
    schedule_ratio=0.9,  # ratio of epochs for learning rate schedule
    save_eval_interval=5,
    fast_transformer=True,
    pre_norm=False,
    amp=True,  # Automatic Mixed Precision
    include_zero_gene = True,
    freeze = False, #freeze
    DSBN = False,  # Domain-spec batchnorm
)


# In[ ]:


run = wandb.init(
    config=hyperparameter_defaults,
    project="Human_sc_mTEC_scGPT_predict_for_zebrafish_mimetic",
    reinit=True,
    settings=wandb.Settings(start_method="fork"),
)
config = wandb.config
print(config)

set_seed(config.seed)


# In[ ]:


## Parameters can adjust:
# max_seq_len, input_style, output_style, log_interval

# settings for input and preprocessing
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
mask_ratio = config.mask_ratio
mask_value = "auto"  # for masked values, now it should always be auto

include_zero_gene = config.include_zero_gene  # if True, include zero genes among hvgs in the training
max_seq_len = 8800 # Set max_seq_len to the number of genes in the sc dataset (or slightly larger)
n_bins = config.n_bins

# input/output representation
input_style = "binned"  # "normed_raw", "log1p", or "binned"
output_style = "binned"  # "normed_raw", "log1p", or "binned"

# settings for training
MLM = False  # whether to use masked language modeling, currently it is always on.
CLS = True  # celltype classification objective
ADV = False  # Adversarial training for batch correction
CCE = False  # Contrastive cell embedding objective
MVC = config.MVC  # Masked value prediction for cell embedding
ECS = config.ecs_thres > 0  # Elastic cell similarity objective
DAB = False  # Domain adaptation by reverse backpropagation, set to 2 for separate optimizer
INPUT_BATCH_LABELS = False  # TODO: have these help MLM and MVC, while not to classifier
input_emb_style = "continuous"  # "category" or "continuous" or "scaling"
cell_emb_style = "cls"  # "avg-pool" or "w-pool" or "cls"
adv_E_delay_epochs = 0  # delay adversarial training on encoder for a few epochs
adv_D_delay_epochs = 0
mvc_decoder_style = "inner product"
ecs_threshold = config.ecs_thres
dab_weight = config.dab_weight

explicit_zero_prob = MLM and include_zero_gene  # whether explicit bernoulli for zeros
do_sample_in_train = False and explicit_zero_prob  # sample the bernoulli in training

per_seq_batch_sample = False

# settings for optimizer
lr = config.lr  # TODO: test learning rate ratio between two tasks
lr_ADV = 1e-3  # learning rate for discriminator, used when ADV is True
batch_size = config.batch_size
eval_batch_size = config.batch_size
epochs = config.epochs
schedule_interval = 1

# settings for the model
fast_transformer = config.fast_transformer
fast_transformer_backend = "flash"  # "linear" or "flash"
embsize = config.layer_size  # embedding dimension
d_hid = config.layer_size  # dimension of the feedforward network in TransformerEncoder
nlayers = config.nlayers  # number of TransformerEncoderLayer in TransformerEncoder
nhead = config.nhead  # number of heads in nn.MultiheadAttention
dropout = config.dropout  # dropout probability

# logging
log_interval = 100  # iterations
save_eval_interval = config.save_eval_interval  # epochs
do_eval_scib_metrics = True


# In[ ]:


# %% validate settings
assert input_style in ["normed_raw", "log1p", "binned"]
assert output_style in ["normed_raw", "log1p", "binned"]
assert input_emb_style in ["category", "continuous", "scaling"]
if input_style == "binned":
    if input_emb_style == "scaling":
        raise ValueError("input_emb_style `scaling` is not supported for binned input.")
elif input_style == "log1p" or input_style == "normed_raw":
    if input_emb_style == "category":
        raise ValueError(
            "input_emb_style `category` is not supported for log1p or normed_raw input."
        )

if input_emb_style == "category":
    mask_value = n_bins + 1
    pad_value = n_bins  # for padding gene expr values
    n_input_bins = n_bins + 2
else:
    mask_value = -1
    pad_value = -2
    n_input_bins = n_bins

if ADV and DAB:
    raise ValueError("ADV and DAB cannot be both True.")
DAB_separate_optim = True if DAB > 1 else False


# In[ ]:


# change the name of save_dir

dataset_name = config.dataset_name
save_dir = Path(f"./save_human_scGPT/dev_{dataset_name}-{time.strftime('%b%d-%H-%M')}_predict_for_zebrafish_mimetic/")
save_dir.mkdir(parents=True, exist_ok=True)
print(f"save to {save_dir}")
logger = scg.logger
scg.utils.add_file_handler(logger, save_dir / "run.log")


# In[ ]:


## "batch_id" 0 for ref dataset for fine-tuning, and 1 for query
## Change data_is_raw accordingly

if dataset_name == "zebrafish_mtec":
    data_dir = Path("/n/groups/cbdm_lab/brh040/analysis/scGPT/mimetics-only_version/zebrafish")
    
    adata = sc.read(data_dir / "zebrafish_mec_aftercut5-humanlabels_mimeticsonly-20240306_cleaned.h5ad")	

    adata_test = adata.copy() 

    adata_test.obs["batch_id"] = adata_test.obs["str_batch"] = "1"     

    print(np.sum(adata.var.index == adata_test.var.index))
    adata_test.var.set_index(adata.var.index, inplace=True)

    data_is_raw = True
    filter_gene_by_counts = False
    adata_test_raw = adata_test.copy()
    adata = adata_test.copy()

# make the batch category column
batch_id_labels = adata.obs["str_batch"].astype("int8")
adata.obs["batch_id"] = batch_id_labels # make batch_id to int8
adata.var["gene_name"] = adata.var.index.tolist() # make it to non-category


# 8719

# In[ ]:


adata_test#.obs


# In[ ]:


adata_test.obs#['seurat_clusters']


# In[ ]:


np.unique(adata_test.obs['seurat_clusters'])


# In[ ]:


if config.load_model is not None:
    model_dir = Path(config.load_model)
    model_config_file = model_dir / "args.json"
    model_file = model_dir / "model_best.pt"
    vocab_file = model_dir / "vocab.json"

    vocab = GeneVocab.from_file(vocab_file)
    shutil.copy(vocab_file, save_dir / "vocab.json")
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)

    adata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in adata.var["gene_name"]
    ]
    gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
    logger.info(
        f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
        f"in vocabulary of size {len(vocab)}."
    )

    # Only take the genes that are in the vocabulary
    adata = adata[:, adata.var["id_in_vocab"] >= 0]

    # model
    # The last five model args will be overide by the model_config_file
    with open(model_config_file, "r") as f:
        model_configs = json.load(f)
    logger.info(
        f"Resume model from {model_file}, the model args will override the "
        f"config {model_config_file}."
    )
    embsize = model_configs["embsize"]
    nhead = model_configs["nheads"]
    d_hid = model_configs["d_hid"]
    nlayers = model_configs["nlayers"]
    n_layers_cls = model_configs["n_layers_cls"]


# In[ ]:


# set up the preprocessor, use the args to config the workflow
preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=filter_gene_by_counts,  # step 1
    filter_cell_by_counts=False,  # step 2
    normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
    result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
    log1p=data_is_raw,  # 4. whether to log1p the normalized data
    result_log1p_key="X_log1p",
    subset_hvg=False,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
    binning=n_bins,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)

adata_test = adata[adata.obs["str_batch"] == "1"]
preprocessor(adata_test, batch_key=None)


# scGPT - INFO - Filtering cells by counts ...
# scGPT - INFO - Normalizing total counts ...
# scGPT - INFO - Log1p transforming ...
# scGPT - INFO - Binning data ...

# In[ ]:


genes = adata.var["gene_name"].tolist()
if config.load_model is None:
    vocab = Vocab(
        VocabPybind(genes + special_tokens, None)
    )  # bidirectional lookup [gene <-> int]
vocab.set_default_index(vocab["<pad>"])
gene_ids = np.array(vocab(genes), dtype=int)


# In[ ]:


def prepare_data(sort_seq_batch=False) -> Tuple[Dict[str, torch.Tensor]]:
    masked_values_train = random_mask_value(
        tokenized_train["values"],
        mask_ratio=mask_ratio,
        mask_value=mask_value,
        pad_value=pad_value,
    )
    masked_values_valid = random_mask_value(
        tokenized_valid["values"],
        mask_ratio=mask_ratio,
        mask_value=mask_value,
        pad_value=pad_value,
    )
    print(
        f"random masking at epoch {epoch:3d}, ratio of masked values in train: ",
        f"{(masked_values_train == mask_value).sum() / (masked_values_train - pad_value).count_nonzero():.4f}",
    )

    input_gene_ids_train, input_gene_ids_valid = (
        tokenized_train["genes"],
        tokenized_valid["genes"],
    )
    input_values_train, input_values_valid = masked_values_train, masked_values_valid
    target_values_train, target_values_valid = (
        tokenized_train["values"],
        tokenized_valid["values"],
    )

    tensor_batch_labels_train = torch.from_numpy(train_batch_labels).long()
    tensor_batch_labels_valid = torch.from_numpy(valid_batch_labels).long()

    tensor_celltype_labels_train = torch.from_numpy(train_celltype_labels).long()
    tensor_celltype_labels_valid = torch.from_numpy(valid_celltype_labels).long()

    if sort_seq_batch:  # TODO: update to random pick seq source in each traning batch
        train_sort_ids = np.argsort(train_batch_labels)
        input_gene_ids_train = input_gene_ids_train[train_sort_ids]
        input_values_train = input_values_train[train_sort_ids]
        target_values_train = target_values_train[train_sort_ids]
        tensor_batch_labels_train = tensor_batch_labels_train[train_sort_ids]
        tensor_celltype_labels_train = tensor_celltype_labels_train[train_sort_ids]

        valid_sort_ids = np.argsort(valid_batch_labels)
        input_gene_ids_valid = input_gene_ids_valid[valid_sort_ids]
        input_values_valid = input_values_valid[valid_sort_ids]
        target_values_valid = target_values_valid[valid_sort_ids]
        tensor_batch_labels_valid = tensor_batch_labels_valid[valid_sort_ids]
        tensor_celltype_labels_valid = tensor_celltype_labels_valid[valid_sort_ids]

    train_data_pt = {
        "gene_ids": input_gene_ids_train,
        "values": input_values_train,
        "target_values": target_values_train,
        "batch_labels": tensor_batch_labels_train,
        "celltype_labels": tensor_celltype_labels_train,
    }
    valid_data_pt = {
        "gene_ids": input_gene_ids_valid,
        "values": input_values_valid,
        "target_values": target_values_valid,
        "batch_labels": tensor_batch_labels_valid,
        "celltype_labels": tensor_celltype_labels_valid,
    }

    return train_data_pt, valid_data_pt


# dataset
class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}


# data_loader
def prepare_dataloader(
    data_pt: Dict[str, torch.Tensor],
    batch_size: int,
    shuffle: bool = False,
    intra_domain_shuffle: bool = False,
    drop_last: bool = False,
    num_workers: int = 0,
) -> DataLoader:
    if num_workers == 0:
        num_workers = min(len(os.sched_getaffinity(0)), batch_size // 2)

    dataset = SeqDataset(data_pt)

    if per_seq_batch_sample:
        # find the indices of samples in each seq batch
        subsets = []
        batch_labels_array = data_pt["batch_labels"].numpy()
        for batch_label in np.unique(batch_labels_array):
            batch_indices = np.where(batch_labels_array == batch_label)[0].tolist()
            subsets.append(batch_indices)
        data_loader = DataLoader(
            dataset=dataset,
            batch_sampler=SubsetsBatchSampler(
                subsets,
                batch_size,
                intra_subset_shuffle=intra_domain_shuffle,
                inter_subset_shuffle=shuffle,
                drop_last=drop_last,
            ),
            num_workers=num_workers,
            pin_memory=True,
        )
        return data_loader

    data_loader = DataLoader(
        dataset=dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=drop_last,
        num_workers=num_workers,
        pin_memory=True,
    )
    return data_loader


# In[ ]:


## All model parameters are not frozen

num_types = 20
num_batch_types = 2

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

ntokens = len(vocab)  # size of vocabulary
model = TransformerModel(
    ntokens,
    embsize,
    nhead,
    d_hid,
    nlayers,
    nlayers_cls=3, # Used for decoder for classification task.
    n_cls=num_types if CLS else 1, # Cell type num
    vocab=vocab,
    dropout=dropout,
    pad_token=pad_token,
    pad_value=pad_value,
    do_mvc=MVC,
    do_dab=DAB,
    use_batch_labels=INPUT_BATCH_LABELS,
    num_batch_labels=num_batch_types,
    domain_spec_batchnorm=config.DSBN,
    input_emb_style=input_emb_style,
    n_input_bins=n_input_bins,
    cell_emb_style=cell_emb_style,
    mvc_decoder_style=mvc_decoder_style,
    ecs_threshold=ecs_threshold,
    explicit_zero_prob=explicit_zero_prob,
    use_fast_transformer=fast_transformer,
    fast_transformer_backend=fast_transformer_backend,
    pre_norm=config.pre_norm,
)
if config.load_model is not None:
    try:
        model.load_state_dict(torch.load(model_file))
        logger.info(f"Loading all model params from {model_file}")
    except:
        # only load params that are in the model and match the size
        model_dict = model.state_dict()
        pretrained_dict = torch.load(model_file)
        pretrained_dict = {
            k: v
            for k, v in pretrained_dict.items()
            if k in model_dict and v.shape == model_dict[k].shape
        }
        for k, v in pretrained_dict.items():
            logger.info(f"Loading params {k} with shape {v.shape}")
        model_dict.update(pretrained_dict)
        model.load_state_dict(model_dict)

param_count = sum(dict((p.data_ptr(), p.numel()) for p in model.parameters() if p.requires_grad).values())

logger.info(f"Total Params {(param_count )}")

wandb.log(
        {
            "info/param_count": param_count,
        },
)

model.to(device)
wandb.watch(model)

if ADV:
    discriminator = AdversarialDiscriminator(
        d_model=embsize,
        n_cls=num_batch_types,
    ).to(device)


# In[ ]:


criterion = masked_mse_loss
criterion_cls = nn.CrossEntropyLoss()
criterion_dab = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(
    model.parameters(), lr=lr, eps=1e-4 if config.amp else 1e-8
)
scheduler = torch.optim.lr_scheduler.StepLR(
    optimizer, schedule_interval, gamma=config.schedule_ratio
)
if DAB_separate_optim:
    optimizer_dab = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler_dab = torch.optim.lr_scheduler.StepLR(
        optimizer_dab, schedule_interval, gamma=config.schedule_ratio
    )
if ADV:
    criterion_adv = nn.CrossEntropyLoss()  # consider using label smoothing
    optimizer_E = torch.optim.Adam(model.parameters(), lr=lr_ADV)
    scheduler_E = torch.optim.lr_scheduler.StepLR(
        optimizer_E, schedule_interval, gamma=config.schedule_ratio
    )
    optimizer_D = torch.optim.Adam(discriminator.parameters(), lr=lr_ADV)
    scheduler_D = torch.optim.lr_scheduler.StepLR(
        optimizer_D, schedule_interval, gamma=config.schedule_ratio
    )

scaler = torch.cuda.amp.GradScaler(enabled=config.amp)


# In[ ]:


def evaluate_raw(model: nn.Module, loader: DataLoader, return_raw: bool = False) -> float:
    """
    Evaluate the model on the evaluation data.
    """
    model.eval()
    total_num = 0
    predictions = []
    pred_scores_list = []
    i = 0
    with torch.no_grad():
        for batch_data in loader:
            input_gene_ids = batch_data["gene_ids"].to(device)
            input_values = batch_data["values"].to(device)
            target_values = batch_data["target_values"].to(device)
            batch_labels = batch_data["batch_labels"].to(device)
            #celltype_labels = batch_data["celltype_labels"].to(device)

            src_key_padding_mask = input_gene_ids.eq(vocab[pad_token])
            with torch.cuda.amp.autocast(enabled=config.amp):
                output_dict = model(
                    input_gene_ids,
                    input_values,
                    src_key_padding_mask=src_key_padding_mask,
                    batch_labels=batch_labels if INPUT_BATCH_LABELS or config.DSBN else None,
                    CLS=CLS,  # evaluation does not need CLS or CCE
                    CCE=False,
                    MVC=False,
                    ECS=False,
                    do_sample=do_sample_in_train,
                    #generative_training = False,
                )
                output_values = output_dict["cls_output"]
                #loss = criterion_cls(output_values, celltype_labels)

                if DAB:
                    loss_dab = criterion_dab(output_dict["dab_output"], batch_labels)

            total_num += len(input_gene_ids)

            pred_scores = output_values.cpu().numpy()
            pred_scores_list.append(pred_scores)
            
            preds = output_values.argmax(1).cpu().numpy()
            predictions.append(preds)
            print(i)
            i += 1
            
    if return_raw:
        return np.concatenate(predictions, axis=0), np.concatenate(pred_scores_list, axis=0)

    return


# In[ ]:


# %% inference
def test_raw_mimetic(model: nn.Module, adata: DataLoader) -> float:
    all_counts = (
        adata.layers[input_layer_key].A
        if issparse(adata.layers[input_layer_key])
        else adata.layers[input_layer_key]
    )

    batch_ids = adata.obs["batch_id"].tolist()
    batch_ids = np.array(batch_ids)

    tokenized_test = tokenize_and_pad_batch(
        all_counts,
        gene_ids,
        max_len=max_seq_len,
        vocab=vocab,
        pad_token=pad_token,
        pad_value=pad_value,
        append_cls=True,  # append <cls> token at the beginning
        include_zero_gene=include_zero_gene,
    )

    input_values_test = random_mask_value(
        tokenized_test["values"],
        mask_ratio=mask_ratio,
        mask_value=mask_value,
        pad_value=pad_value,
    )

    test_data_pt = {
        "gene_ids": tokenized_test["genes"],
        "values": input_values_test,
        "target_values": tokenized_test["values"],
        "batch_labels": torch.from_numpy(batch_ids).long(),
    }

    test_loader = DataLoader(
        dataset=SeqDataset(test_data_pt),
        batch_size=eval_batch_size,
        shuffle=False,
        drop_last=False,
        num_workers=min(len(os.sched_getaffinity(0)), eval_batch_size // 2),
        pin_memory=True,
    )

    model.eval()
    predictions_raw, prediction_scores = evaluate_raw(
        model,
        loader=test_loader,
        return_raw=True,
    )

    softmax_scores = []

    for i in range(len(prediction_scores)):
        exp_prediction_scores = np.exp(prediction_scores[i])
        softmaxscore = exp_prediction_scores/np.sum(exp_prediction_scores)
        softmax_scores.append(softmaxscore)

    predictions = [ori_pred if np.max(softmax_scores[i]) > 0.5 else -1 for i,ori_pred in enumerate(list(predictions_raw))]
    predictions_new = predictions

    return predictions_new, softmax_scores


# In[ ]:


input_layer_key = {  # the values of this map coorespond to the keys in preprocessing
    "normed_raw": "X_normed",
    "log1p": "X_normed",
    "binned": "X_binned",
}[input_style]


# In[ ]:


predictions, softmax_scores = test_raw_mimetic(model, adata_test)


# In[ ]:

with open(save_dir / "predictions.npy", "wb") as f:
    np.save(f, np.array(predictions))

with open(save_dir / "softmax_scores.npy", "wb") as f:
    np.save(f, np.array(softmax_scores))


# In[ ]:


data_dir_original = Path("/n/groups/cbdm_lab/brh040/analysis/scGPT/mimetics-only_version/human_data")

adata_ori = sc.read(data_dir_original / "meclo_ht2-ht6_seurat_table_20240309_labelled-mimeticsonly_preprocessed.h5ad")

adata_test_ori = adata_ori[adata_ori.obs["orig.ident"]=="HT6",:].copy() 
adata_ref_ori = adata_ori[adata_ori.obs["orig.ident"].isin(["HT2", "HT3", "HT4", "HT5"]),:].copy() 

adata_test_ori.obs["celltype"] = adata_test_ori.obs["celltype"].astype("category")
adata_ref_ori.obs["celltype"] = adata_ref_ori.obs["celltype"].astype("category")

adata_ref_ori.obs["batch_id"] = adata_ref_ori.obs["str_batch"] = "0"
adata_test_ori.obs["batch_id"] = adata_test_ori.obs["str_batch"] = "1"     

print(np.sum(adata_ori.var.index == adata_ref_ori.var.index))
print(np.sum(adata_ori.var.index == adata_test_ori.var.index))
adata_ref_ori.var.set_index(adata_ori.var.index, inplace=True)
adata_test_ori.var.set_index(adata_ori.var.index, inplace=True)
print(np.sum(adata_ref_ori.var.index == adata_test_ori.var.index))

adata_ori = adata_ref_ori.concatenate(adata_test_ori, batch_key="str_batch")

# make the batch category column
batch_id_labels = adata_ori.obs["str_batch"].astype("category").cat.codes.values
adata_ori.obs["batch_id"] = batch_id_labels # make batch_id to int8
celltype_id_labels = adata_ori.obs["celltype"].astype("category").cat.codes.values # int8
celltypes = adata_ori.obs["celltype"].unique()
num_types = len(np.unique(celltype_id_labels))
id2type = dict(enumerate(adata_ori.obs["celltype"].astype("category").cat.categories))
adata_ori.obs["celltype_id"] = celltype_id_labels
adata_ori.var["gene_name"] = adata_ori.var.index.tolist() # make it to non-category



# In[ ]:


celltypes_new = adata_ori.obs["celltype"].unique()
celltypes_new = list(celltypes_new.categories)
celltypes_new.append('Novel')
cell_types_correct_order = [celltypes_new[-1]]
cell_types_correct_order.extend(celltypes_new[:-1])


# In[ ]:


id2type


# In[ ]:


id2type[-1] = 'Novel'
adata_test_raw.obs["predictions"] = [id2type[p] for p in predictions]


# In[ ]:


palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] 
palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"]
palette_ = {c: palette_[i] for i, c in enumerate(cell_types_correct_order)}
palette_['Novel'] = "black"


# In[ ]:


with plt.rc_context({"figure.figsize": (10, 5), "figure.dpi": (300)}):
    sc.pl.umap(
        adata_test_raw,
        color=["predictions"],
        palette=palette_,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(save_dir / "results_with_novel_predictions.png", dpi=300)


# In[ ]:


id2type_new = id2type.copy()
id2type_new


# In[ ]:


adata_test_raw.obs["predictions_mimetic"] = [id2type_new[p] for p in predictions]


# In[ ]:


with plt.rc_context({"figure.figsize": (10, 5), "figure.dpi": (300)}):
    sc.pl.umap(
        adata_test_raw,
        color=["predictions_mimetic"],
        palette=palette_,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(save_dir / "results_with_novel_predictions_mimetic.png", dpi=300)


# In[ ]:


def get_cell_embedding(
    model: nn.Module,
    adata: AnnData,
    include_types: List[str] = ["cls"],
) -> Optional[Dict]:
    """evaluate the model on test dataset of adata_t"""
    model.eval()

    # copy adata_t to avoid reuse previously computed results stored in adata_t
    adata_t = adata.copy()
    adata_t.obs["predictions"] = [id2type[p] for p in predictions]
    adata_t.obs["predictions_mimetic"] = [id2type_new[p] for p in predictions]

    all_counts = (
        adata_t.layers[input_layer_key].A
        if issparse(adata_t.layers[input_layer_key])
        else adata_t.layers[input_layer_key]
    )

    batch_ids = adata_t.obs["batch_id"].tolist()
    batch_ids = np.array(batch_ids)

    # Evaluate cls cell embeddings
    if "cls" in include_types:
        logger.info("Evaluating cls cell embeddings")
        tokenized_all = tokenize_and_pad_batch(
            all_counts,
            gene_ids,
            max_len=max_seq_len,
            vocab=vocab,
            pad_token=pad_token,
            pad_value=pad_value,
            append_cls=True,  # append <cls> token at the beginning
            include_zero_gene=include_zero_gene, # check include_zero_gene=include_zero_gene,
        )
        all_gene_ids, all_values = tokenized_all["genes"], tokenized_all["values"]
        src_key_padding_mask = all_gene_ids.eq(vocab[pad_token])
        with torch.no_grad(), torch.cuda.amp.autocast(enabled=config.amp):
            cell_embeddings = model.encode_batch(
                all_gene_ids,
                all_values.float(),
                src_key_padding_mask=src_key_padding_mask,
                batch_size=config.batch_size,
                batch_labels=torch.from_numpy(batch_ids).long() if config.DSBN else None,
                time_step=0,
                return_np=True,
            )
        cell_embeddings = cell_embeddings / np.linalg.norm(
            cell_embeddings, axis=1, keepdims=True
        )

        palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] 
        palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"]
        palette_ = {c: palette_[i] for i, c in enumerate(cell_types_correct_order)}

        adata_t.obsm["X_scGPT"] = cell_embeddings
        sc.pp.neighbors(adata_t, use_rep="X_scGPT")
        sc.tl.umap(adata_t, min_dist=0.3)
            
        with plt.rc_context({"figure.figsize": (20, 10), "figure.dpi": (300)}):
            sc.pl.umap(
                adata_t,
                color=["predictions"],
                palette=palette_,
                show=False,
            )
            plt.savefig(save_dir / "results_with_novel_cell_embed.png", dpi=300)
            
        with plt.rc_context({"figure.figsize": (20, 10), "figure.dpi": (300)}):
            sc.pl.umap(
                adata_t,
                color=["predictions_mimetic"],
                palette=palette_,
                show=False,
            )
            plt.savefig(save_dir / "results_with_novel_cell_embed_mimetic.png", dpi=300)

    if len(include_types) == 1:
        return adata_t


# In[ ]:


adata_test_cell_emb = get_cell_embedding(model, adata_test, include_types=["cls"])


# In[ ]:


palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] 
palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"]
palette_ = {c: palette_[i] for i, c in enumerate(cell_types_correct_order)}
palette_['Novel'] = "black"


# In[ ]:


with plt.rc_context({"figure.figsize": (10, 5), "figure.dpi": (300)}):
    sc.pl.umap(
        adata_test_cell_emb,
        color=["predictions"],
        palette=palette_,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(save_dir / "results_with_novel_cell_embed_predictions.png", dpi=300)


# In[ ]:


with plt.rc_context({"figure.figsize": (10, 5), "figure.dpi": (300)}):
    sc.pl.umap(
        adata_test_cell_emb,
        color=["predictions_mimetic"],
        palette=palette_,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(save_dir / "results_with_novel_cell_embed_predictions_mimetic.png", dpi=300)


# In[ ]:


with plt.rc_context({"figure.figsize": (10, 8), "figure.dpi": (300)}):
    sc.pl.umap(
        adata_test_cell_emb,
        color=["predictions"],
        palette=palette_,
        legend_loc='on data', 
        legend_fontsize=10,
        show=False,
    )
    plt.savefig(save_dir / "results_with_novel_cell_embed_predictions_v2.png", dpi=300)


# In[ ]:


with plt.rc_context({"figure.figsize": (10, 8), "figure.dpi": (300)}):
    sc.pl.umap(
        adata_test_cell_emb,
        color=["predictions_mimetic"],
        palette=palette_,
        legend_loc='on data', 
        legend_fontsize=10,
        show=False,
    )
    plt.savefig(save_dir / "results_with_novel_cell_embed_predictions_mimetic_v2.png", dpi=300)


# In[ ]:


adata_test_cell_emb.write(save_dir / "adata_test_cell_emb_large_v2_zebrafish_mimetic.h5ad")


# ## Export data for Seurat

# In[ ]:


import scanpy as sc
import numpy as np

import anndata2ri


# In[ ]:


import pandas as pd


# In[ ]:


### EDIT FOLDER
adata_emb = adata_test_cell_emb 
adata_emb


# In[ ]:


adata_emb.obs


# In[ ]:


adata_emb.obsm['X_scGPT']


# In[ ]:


adata_emb.obsm['X_scGPT'].shape


# In[ ]:


emb_df= pd.DataFrame(adata_emb.obsm['X_scGPT'], index=adata_emb.obs.index.tolist())
emb_df


# In[ ]:


### EDIT FOLDER
emb_df.to_csv(save_dir / 'zebrafish_embeddings_mimetics-only_version.csv', index=True)


# In[ ]:


### EDIT FOLDER
adata_emb.obs.to_csv(save_dir / 'zebrafish_obs_info_mimetics-only_version.csv', index=True)


# In[ ]:




