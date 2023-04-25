import argparse, sys, os
import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
import fnmatch
from sklearn import svm

## parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description='Drug response prediction on new drugs')

parser.add_argument('-i', '--input', required=True, help='path to input PRISM prediction results (csv file)')
parser.add_argument('-smiles', '--input_smiles', required=True, help='path to input SMILES encoding of molecules (txt file)')
parser.add_argument('-o', '--output', default='./', help='path to output directory, default=\'./\'')

args = parser.parse_args()

# check input
if not os.path.exists(args.input):
    sys.exit('The input path does not exist.')
if fnmatch.fnmatch(args.input, "*.csv") == False:
    sys.exit('The input file is not a csv file.')

if not os.path.exists(args.input_smiles):
    sys.exit('The input path does not exist.')
if fnmatch.fnmatch(args.input_smiles, "*.txt") == False:
    sys.exit('The input file is not a txt file.')


# check output
if not os.path.isdir(args.output):
    sys.exit('The output directory does not exist.')


class new_drug_prediction:
    def _init__(self):
        self.load_prediction()
        self.prepare_fingerprints()
        self.load_input_mol()
        self.sensitivity_prediction()
        self.compute_output()
        
    def load_prediction(self):
        self.cadrres_pred = pd.read_csv(args.input, index_col=0, header=[0,1])
        mapping_df = pd.read_csv('/scDrug/data/drugID_smiles_map.csv', index_col = 0)
        self.train_smiles = mapping_df.loc[[self.cadrres_pred.columns.get_level_values(0), 'smiles'].values]
    
    def RDKfp_convert(self, smiles_ls):
        mol_rdkit = list(map(Chem.MolFromSmiles,smiles_ls))
        fps = [list(map(int, list(Chem.RDKFingerprint(x).ToBitString()))) for x in mol_rdkit]
        fps = np.array(fps)
        return fps
    
    def load_input_mol(self):
        with open(args.input_smiles) as f:
            self.input_smiles = f.read().strip().split("\n")

    def prepare_fingerprints(self):
        self.X = self.RDKfp_convert(self.train_smiles)
        self.input_mol_fps = self.RDKfp_convert(self.input_mol)
    
    def sensitivity_prediction(self):
        pred_auc = np.empty((len(self.input_smiles), len(self.cadrres_pred.index)))
        pred_auc[:] = np.nan
        self.pred_auc_df = pd.DataFrame(pred_auc, columns= self.cadrres_pred.index.tolist(), index = self.input_smiles)

        for cluster in self.cadrres_pred.index:
            y = self.cadrres_pred.loc[cluster]
            model = svm.SVR()
            model.fit(self.X, y)
            cluster_pred = model.predict(self.input_mol_fps)
            self.pred_auc_df[cluster] = cluster_pred
        self.pred_auc_df.to_csv(os.path.join(args.output, "new_drug_prediction.csv"), index = True)
    
    def compute_output(self):
        self.pred_auc_output = self.pred_auc_df.reset_index().melt(id_vars=["index"], var_name = 'cluster', 
                                                                   value_vars=self.pred_auc_df.columns.tolist(), 
                                                                   value_name = "AUC prediction")
        self.pred_auc_output['classification'] = ['potnetial' if pred > 0.6 else ('inactive' if pred < 0.2 else 'unclear') 
                                                  for pred in self.pred_auc_output['AUC prediction']]
        self.pred_auc_output['rank'] = list(map(int, self.pred_auc_output.groupby("cluster")["AUC prediction"].rank(ascending = False)))
        self.pred_auc_output = self.pred_auc_output.set_index(["cluster", 'index'])
        self.pred_auc_output.to_csv(os.path.join(args.output, "drug_level_prediction.csv"), index = True)

job = new_drug_prediction()