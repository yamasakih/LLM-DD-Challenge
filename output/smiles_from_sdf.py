import os
from rdkit import Chem
from rdkit.Chem import AllChem

input_file = '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf'
output_file = '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.txt'

suppl = Chem.SDMolSupplier(input_file)

with open(output_file, 'w') as f:
    for mol in suppl:
        if mol is not None:
            ic50_range = mol.GetProp('IC50_range_nM')
            if ic50_range == '<100':
                smiles = Chem.MolToSmiles(mol)
                f.write(smiles + '\n')