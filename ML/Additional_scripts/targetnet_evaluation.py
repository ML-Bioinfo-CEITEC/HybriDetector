# to run the TargetNet evaluation you need to install it first
# run: pip install git+https://github.com/katarinagresova/TargetNet.git
#
# you also need to donwload a pretrained TargetNet model into Models
# wget https://github.com/katarinagresova/TargetNet/blob/master/pretrained_models/TargetNet.pt?raw=true -O Models/TargetNet.pt
#
# and then you need to change the array containing paths to dataset with tree columns: mirna, mrna, label (names can be different, order is important)
# 
# script will create a folder `targetnet_metrics` with predictions for all datasets in the array

import sys

import torch
import pandas as pd
from pathlib import Path

from targetnet.train import Trainer
from targetnet.utils import set_seeds, set_output
from targetnet.model.model_utils import get_model
from targetnet.data import get_dataset_from_configs
from targetnet.data import reverse

# basicaly just a copy of the original TargetNet implementation, but I can pass config as a dict instead of a file
class ModelConfig():
    def __init__(self, cfg, idx="model_config"):
        """ model configurations """
        self.idx = idx
        self.type = None
        self.num_channels = None
        self.num_blocks = None
        self.stem_kernel_size = None
        self.block_kernel_size = None
        self.pool_size = None

        for key, value in cfg.items():
            if key == "skip_connection":                self.skip_connection = value
            elif key == "num_channels":                 self.num_channels = value
            elif key == "num_blocks":                   self.num_blocks = value
            elif key == "stem_kernel_size":             self.stem_kernel_size = value
            elif key == "block_kernel_size":            self.block_kernel_size = value
            elif key == "pool_size":                    self.pool_size = value
            else: sys.exit("# ERROR: invalid key [%s] in model-config file" % key)

# basicaly just a copy of the original TargetNet implementation, but I can pass config as a dict instead of a file
class DataConfig():
    def __init__(self, cfg, idx="data_config"):
        """ data configurations """
        self.idx = idx
        self.with_esa = True
        self.path = {}

        for key, value in cfg.items():
            if "with_esa" in key:             self.with_esa = value
            elif "path" in key:               self.path[key.split("_")[0]] = value
            else: sys.exit("# ERROR: invalid key [%s] in data-config file" % key)

def prepare_dataset(dataset):

    df = pd.read_csv(dataset, sep='\t')
    df.columns = ['mirna_seq', 'mrna_seq', 'label']
    df['mirna_id'] = df['mirna_seq']
    df['mrna_id'] = df['mrna_seq']
    # TargetNet is reversing the sequence internaly, but we don't need to reverse it - so we reverse it here and then it is reversed back
    df['mrna_seq'] = df['mrna_seq'].apply(lambda x: reverse(x))
    df = df.reindex(columns=['mirna_id', 'mirna_seq', 'mrna_id', 'mrna_seq', 'label'])

    dataset = Path(dataset)
    targetnet_dataset = dataset.with_name(dataset.stem + '_targetnet' + dataset.suffix)

    df.to_csv(targetnet_dataset, sep='\t', index=False)
    return targetnet_dataset


if __name__ == "__main__":

    DATASETS = [
        "../Datasets/test_set_1_1_CLASH2013_paper.tsv",
        "../Datasets/test_set_1_10_CLASH2013_paper.tsv",
        "../Datasets/test_set_1_100_CLASH2013_paper.tsv"
    ]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    set_seeds(2020)

    model_cfg = ModelConfig({
        'skip_connection': True,
        'num_channels': [16, 16, 32],
        'num_blocks': [2, 1, 1],
        'stem_kernel_size': 5,
        'block_kernel_size': 3,
        'pool_size': 3
    })

    dset_cfg = {}
    for index, dataset in enumerate(DATASETS):

        targetnet_dataset = prepare_dataset(dataset)
        dset_cfg['dataset-' + str(index) + '_path'] = targetnet_dataset

    data_cfg = DataConfig(dset_cfg)

    output, save_prefix = set_output({'output_path': 'targetnet_metrics/'}, "evaluate_model_log")

    model, params = get_model(model_cfg, with_esa=True)

    trainer = Trainer(model)
    trainer.load_model("../Models/TargetNet.pt", output)
    trainer.set_device(device)

    dataset_idxs, datasets, iterators = data_cfg.path.keys(), [], []
    for idx in dataset_idxs:
        dataset = get_dataset_from_configs(data_cfg, idx)
        iterator = torch.utils.data.DataLoader(dataset, 64, shuffle=False, pin_memory=True, num_workers=2)

        ### validation
        for B, batch in enumerate(iterator):
            trainer.evaluate(batch, device)

        ### save outputs
        trainer.aggregate(dataset.set_labels)
        trainer.save_outputs(idx, save_prefix)

        ### format predictions
        preds = pd.read_csv(save_prefix + "/%s_outputs.txt" % (idx), sep='\t')
        df = pd.read_csv(data_cfg.path[idx], sep='\t')
        df = pd.concat([df, preds], axis=1)
        df.columns = ['mirna_id', 'mirna', 'mrna_id', 'mrna', 'label', 'set_idx', 'prediction']
        df[['mirna', 'mrna', 'label', 'prediction']].to_csv(save_prefix + '/targetnet_score_' + idx + '.tsv', sep='\t', index=False)
        
        ### cleanup - remove unformated predictions
        Path(save_prefix + "/%s_outputs.txt" % (idx)).unlink()

    ### cleanup - delete log file
    Path(save_prefix + "/evaluate_model_log.txt").unlink()