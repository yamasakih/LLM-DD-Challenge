## はじめに
いくつかのツールを利用可能な状態で Goal (目標)を箇条書きで提示するとステップバイステップで実行するアプリである [AutoGPT](https://github.com/Significant-Gravitas/Auto-GPT) を用いて、[ChemCrow](https://arxiv.org/abs/2304.05376) 同様に、創薬を行うことができるか試みました。

結論を先に述べますと、ある程度ツールを用意することで Docking score を指標とする Reinvent を用いて化学構造を生成することができました。
[ChemCrow](https://arxiv.org/abs/2304.05376)では化学構造の生成や Docking は行っていなかったためその点については新規性があるかなと思います。

## AutoGPT
[AutoGPT](https://github.com/Significant-Gravitas/Auto-GPT) では `ai_goals` セクションを含む設定ファイル `ai_setting.yaml` を指定することでそこに書かれている Goal を実行するために何をすべきかを AI が自発的に提案し、ユーザが許可をすると実行をします。

以下に例をあげます。

```
ai_goals:
- What's is scaffold hopping? Search it and output to scaffold.txt
- Shutdown when task is done.
ai_role: Medicinal Chemist and Chemoinformatician
ai_name: Hoshino AI
api_budget: 3.0
```

`ai_role` として `Medicinal Chemist and Chemoinformatician` を指定していますのでベースである ChatGPT-4 はメディシナルケミストかつケモインフォマティシャンとして振る舞います。
`ai_name`, `api_budget` はオプショナルであり解析には今回は影響ありません。

上記の設定で実行したところ AutoGPT の思考などを含めた出力が得られます。
各セクションは以下のような意味を含んでいます。(具体例はこの後の実作業を確認してください)

- THOUGHTS: AutoGPT が Goal を達成するために何をすべきなどを考えた思考が表示されます
- REASONING: その考えに至った理由が表示されます
- PLAN: Goal を達成するために必要なタスクが箇条書きに表示されます
- CRITICISM: 制約等がある場合表示されます
- NEXT ACTION: AutoGPT が利用可能なアクション (例えばWeb検索で情報を集めるなど) からタスクを達成するために選択したものとその引数が表示されます。

ユーザは AutoGPT の行動を許可する `y` を入力したりネクストアクションを指定する `feedback` といった行動を取ることができますが、今回は LLM ができる限り自発的に行動だけで創薬を行うことをコンセプトにするということで `feedback` を用いていません。

## 手順
実際に AutoGPT に依頼した `ai_golas` を設定し実行してもらった `ai_setting.yaml` とそれに対して AutoGPT が行ったアクションを時系列順に示していきます。
マニュアルで入力に用意したデータやログすべてなどは別途リポジトリにアップロードしています。

### 1. 既知化合物から活性値が高い化学構造だけを抽出する
[CBLB_inhibitors_vsF.sdf](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/input/CBLB_inhibitors_vsF.sdf) から活性値が高いリガンドだけに注目するために `IC50_range_nM` タグが "<100" であるリガンドだけを抽出する `smiles_from_sdf.py` を AutoGPT に実装させました。

[CBLB_inhibitors_vsF.sdf](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/input/CBLB_inhibitors_vsF.sdf) はコンテストにて周知されているものをあらかじめダウンロードしておきました。

[ai_settings_230603_1.yaml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_1.yaml)

```yaml
ai_goals:
- Implement a Python script named smiles_from_sdf.py using RDKit to output the SMILES key of the chemical structure with IC50_range_nM tag "<100" (note that it is a string) from /workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf to /workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.txt using RDKit.
- Run smiles_from_sdf.py
- Shutdown when task is done.
ai_name: Hoshino AI
ai_role: Medicinal Chemist and Chemoinformatician
api_budget: 3.0
```

以下にこの `ai_goals` に対して AutoGPT が行ったアクションの[ログ](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/logs/activity_230603_1.log)を示します。
なお、前半部分に AutoGPT の新機能などを説明する文章が必ず出力される仕様になっていますがコンテストには不要と判断し削除しています。また、DEBUGモードでは AutoGPT 自身がプロンプトを作成し ChatGPT-4 API に利用しています。それらを確認したい場合は [activity_230603_1.log](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/logs/activity_230603_1.log) を直接参照してください。これ以降のログも同様です。

```txt
2023-06-03 15:13:09,686 INFO GPT4 Only Mode:  ENABLED
2023-06-03 15:13:09,817 INFO Using AI Settings File: ai_settings_230603.yaml
2023-06-03 15:13:14,318 INFO  Auto-GPT-Chemoinformatics-plugin: 0.1.0 - This is the Auto-GPT chemoinformatics plugin.
2023-06-03 15:13:14,461 INFO Name : Hoshino AI
2023-06-03 15:13:14,602 INFO Role : Medicinal Chemist and Chemoinformatician
2023-06-03 15:13:14,766 INFO Goals: ['Implement a Python script named smiles_from_sdf.py using RDKit to output the SMILES key of the chemical structure with IC50_range_nM tag "<100" (note that it is a string) from /workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf to /workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.txt using RDKit.', 'Run smiles_from_sdf.py', 'Shutdown when task is done.']
2023-06-03 15:13:15,276 INFO API Budget: $3.0
2023-06-03 15:13:15,396 INFO Hoshino AI has been created with the following details:
2023-06-03 15:13:15,650 INFO Name: Hoshino AI
2023-06-03 15:13:15,763 INFO Role: Medicinal Chemist and Chemoinformatician
2023-06-03 15:13:15,929 INFO Goals: 
2023-06-03 15:13:16,005 INFO - Implement a Python script named smiles_from_sdf.py using RDKit to output the SMILES key of the chemical structure with IC50_range_nM tag "<100" (note that it is a string) from /workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf to /workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.txt using RDKit.
2023-06-03 15:13:16,573 INFO - Run smiles_from_sdf.py
2023-06-03 15:13:16,669 INFO - Shutdown when task is done.
2023-06-03 15:13:16,925 INFO Using memory of type: LocalCache
2023-06-03 15:13:17,076 INFO Using Browser: chrome
2023-06-03 15:14:08,074 INFO  THOUGHTS: I need to create a Python script to extract the SMILES key of chemical structures with IC50_range_nM tag '<100' from the given SDF file.
2023-06-03 15:14:08,535 INFO REASONING: To achieve the goal, I should start by writing a Python script using RDKit to process the SDF file and extract the required information.
2023-06-03 15:14:09,035 INFO PLAN: 
2023-06-03 15:14:09,106 INFO -  Write the Python script
2023-06-03 15:14:09,243 INFO -  Save the script to a file
2023-06-03 15:14:09,427 INFO -  Execute the script
2023-06-03 15:14:09,572 INFO -  Shutdown when the task is done
2023-06-03 15:14:09,830 INFO CRITICISM: None
2023-06-03 15:14:09,959 INFO NEXT ACTION:  COMMAND = write_to_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/smiles_from_sdf.py', 'text': "import os\nfrom rdkit import Chem\nfrom rdkit.Chem import AllChem\n\ninput_file = '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf'\noutput_file = '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.txt'\n\nsuppl = Chem.SDMolSupplier(input_file)\n\nwith open(output_file, 'w') as f:\n    for mol in suppl:\n        if mol is not None:\n            ic50_range = mol.GetProp('IC50_range_nM')\n            if ic50_range == '<100':\n                smiles = Chem.MolToSmiles(mol)\n                f.write(smiles + '\\n')"}
2023-06-03 15:14:10,503 INFO  Enter 'y' to authorise command, 'y -N' to run N continuous commands, 's' to run self-feedback commands, 'n' to exit program, or enter feedback for ...
2023-06-03 15:14:10,506 INFO  Asking user via keyboard...
2023-06-03 15:14:46,372 INFO -=-=-=-=-=-=-= COMMAND AUTHORISED BY USER -=-=-=-=-=-=-= 
2023-06-03 15:14:46,603 INFO SYSTEM:  Command write_to_file returned: File written to successfully.
2023-06-03 15:15:09,148 INFO  THOUGHTS: I have written the Python script to extract the SMILES key of chemical structures with IC50_range_nM tag '<100' from the given SDF file. Now, I need to execute the script.
2023-06-03 15:15:09,611 INFO REASONING: To achieve the goal, I should execute the Python script I just wrote.
2023-06-03 15:15:09,957 INFO PLAN: 
2023-06-03 15:15:10,011 INFO -  Execute the script
2023-06-03 15:15:10,127 INFO -  Shutdown when the task is done
2023-06-03 15:15:10,334 INFO CRITICISM: None
2023-06-03 15:15:10,434 INFO NEXT ACTION:  COMMAND = execute_python_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/smiles_from_sdf.py'}
2023-06-03 15:15:10,597 INFO [36mAUTHORISED COMMANDS LEFT: [0m2 
2023-06-03 15:15:10,718 INFO  Executing file '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/smiles_from_sdf.py'
2023-06-03 15:15:11,276 INFO SYSTEM:  Command execute_python_file returned: 
2023-06-03 15:15:26,615 INFO  THOUGHTS: I have executed the Python script to extract the SMILES key of chemical structures with IC50_range_nM tag '<100' from the given SDF file. The task is complete.
2023-06-03 15:15:27,087 INFO REASONING: Since the task is complete, I should shut down.
2023-06-03 15:15:27,354 INFO PLAN: 
2023-06-03 15:15:27,424 INFO -  Shutdown
2023-06-03 15:15:27,495 INFO CRITICISM: None
2023-06-03 15:15:27,574 INFO NEXT ACTION:  COMMAND = task_complete  ARGUMENTS = {'reason': 'Task complete'}
2023-06-03 15:15:27,850 INFO [36mAUTHORISED COMMANDS LEFT: [0m1 
2023-06-03 15:15:27,967 INFO Shutting down...
 Task complete
```

AutoGPT は RDKit を用いて SDF のタグを調べ `<100` である化学構造の SMILES を保存するスクリプト `smiles_from_sdf.py` を実装し、それを利用することで Goal である「`IC50_range_nM` タグが "<100" であるリガンドだけを抽出する」 ことを達成しています。

以下にコード [smiles_from_sdf.py](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/smiles_from_sdf.py) と成果物である [CBLB_inhibitors_vsF.txt](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/CBLB_inhibitors_vsF.txt) を示しておきます。

```py
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
```

```txt
C[C@@H](Cc1nncn1C)c1cccc(-c2[nH]nc3c(CN4CC[C@H](F)C4)cc(C(F)(F)F)nc23)c1
Cn1cnnc1CC1(c2cccc(NC(=O)N3CC4(CCNCC4)c4cccnc43)c2)COC1
C[C@H](Cc1nncn1C)c1cccc(NC(=O)N2CCOc3cccnc32)c1
C[C@H](Cc1nncn1C)c1cccc(NC(=O)N2CC(C)(C)Oc3cccnc32)c1
Cn1cnnc1CC1(c2cccc(NC(=O)N3CC4(CNC4)Oc4cccnc43)c2)COC1
CN1CC2(C1)CN(C(=O)Nc1cccc(C3(Cc4nncn4C)COC3)c1)c1ncccc1O2
Cn1cnnc1CC1(c2cccc(NC(=O)N3C[C@@]4(CCNC4)Oc4cccnc43)c2)COC1
Cn1cnnc1C(F)(F)[C@](C)(F)c1cccc(NC(=O)N2CC3(CNC3)Oc3cccnc32)c1
Cc1cc2c(c(C3CC3)n1)CN(c1cccc(C3(Cc4nncn4C)COC3)c1)C2=O
Cn1cnnc1CC1(c2cccc(N3Cc4c(cc(C5CC5)nc4C4CC4)C3=O)c2)COC1
...(以下略)
```

### 2. アポ構造のタンパクとリガンドを用意し、さらに Docking のための target protein preparation を行う
続けて Docking のために 8gcy.pdb からアポ構造のタンパクとリガンドのみを含む PDB を出力させました。
こちらは PDB などの理解をしている必要があり現状の AutoGPT では一から実装するの不可能のため [Chemoinformatics plugin](https://github.com/yamasakih/Auto-GPT-Plugin-Chemoinformatics) と言う名前でプラグインを作成し私が実装した関数を利用するように指示しました。

[ai_setting_230603_2.yml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_2.yaml)

```yaml
ai_goals:
- Make protein apo structure and ligand-only PDBs from 8gcy.pdb by using the Chemoinformatics plugin respectively.
- Read target_preparation_example.jsonc to know how to output setting json to run target protein preparation
- Update the setting json file to use the two PDBs you just created protein apo pdb and ligand Z3N.pdb, and output the file as 8gcy_preparation_setting.json.
- Run target protein preparation by using the Chemoinformatics plugin with 8gcy_preparation_setting.json
- Shutdown when task is done.
ai_name: Hoshino AI
ai_role: Medicinal Chemist and Chemoinformatician
api_budget: 3.0
```

Docking は [DockStream](https://github.com/MolecularAI/DockStream) と呼ばれるソフトウェアを用いました。こちらのソフトは Glide, AutodockVina などケモインフォマティクスの分野では広く用いられている Docking ソフトに対応し、JSON 形式で設定ファイルを記載しておくと Docking を実行することができるソフトウェアです。
こちらを採用した理由として、事前調査により LLM は JSON フォーマットの理解度が高く、 JSON の設定ファイルを考えて出力できることができるため、あらかじめ設定ファイルの例として [target_preparation_example.jsonc](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/input/target_preparation_example.jsonc) を私が用意しそちらを読み、今回のデータに適したように設定を書き換えて出力せよと指定するだけで十分 LLM は JSON を出力し Docking を行うことができると考えたからです。なお Docking tool として無償利用可能な AutoDock Vina を採用しています。

DockStream および DockStream と Reinvent の組み合わせなどはこちらの[記事](https://iwatobipen.wordpress.com/2021/11/02/reinforcement-learning-with-docking-score-rdkit-reinvent-chemoinforamtics/)や公式サイトを参考にしました。

最後に作成した JSON を利用し Docking の前処理である target protein preparation を行うように指示しました。
こちらも Chemoinformatics plugin として用意はしていますが DockStream を subprocess で呼び出すだけの簡単な[関数](https://github.com/yamasakih/Auto-GPT-Plugin-Chemoinformatics/blob/f2e84ceb6d1d908ad3de3e0416e3c4ca6d4a2fd6/src/auto_gpt_plugin/chemoinformatics/chemoinformatics.py#L82)にしています。

これらをまとめたを [ai_setting_230603_2.yml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_2.yaml)
 を実行した結果が以下です。 [JSONファイル](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/8gcy_preparation_setting.json) を今回のデータにフィットするように出力し target protein preparation が正常に終了しています。

また、 preparation を実行されたタンパクデータが [8gcy_receptor.pdbqt](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/8gcy_receptor.pdbqt) として保存されました。

```sh
2023-06-03 15:25:28,162 INFO GPT4 Only Mode:  ENABLED
2023-06-03 15:25:28,335 INFO Using AI Settings File: ai_settings_230603-2.yaml
2023-06-03 15:25:32,387 INFO  Auto-GPT-Chemoinformatics-plugin: 0.1.0 - This is the Auto-GPT chemoinformatics plugin.
2023-06-03 15:25:32,507 INFO Name : Hoshino AI
2023-06-03 15:25:32,683 INFO Role : Medicinal Chemist and Chemoinformatician
2023-06-03 15:25:32,864 INFO Goals: ['Make protein apo structure and ligand-only PDBs from 8gcy.pdb by using the Chemoinformatics plugin respectively.', 'Read target_preparation_example.jsonc to know how to output setting json to run target protein preparation', 'Update the setting json file to use the two PDBs you just created protein apo pdb and ligand Z3N.pdb, and output the file as 8gcy_preparation_setting.json.', 'Run target protein preparation by using the Chemoinformatics plugin with 8gcy_preparation_setting.json', 'Shutdown when task is done.']
2023-06-03 15:25:33,692 INFO API Budget: $3.0
2023-06-03 15:25:33,820 INFO Hoshino AI has been created with the following details:
2023-06-03 15:25:34,005 INFO Name: Hoshino AI
2023-06-03 15:25:34,091 INFO Role: Medicinal Chemist and Chemoinformatician
2023-06-03 15:25:34,239 INFO Goals: 
2023-06-03 15:25:34,313 INFO - Make protein apo structure and ligand-only PDBs from 8gcy.pdb by using the Chemoinformatics plugin respectively.
2023-06-03 15:25:34,660 INFO - Read target_preparation_example.jsonc to know how to output setting json to run target protein preparation
2023-06-03 15:25:35,085 INFO - Update the setting json file to use the two PDBs you just created protein apo pdb and ligand Z3N.pdb, and output the file as 8gcy_preparation_setting.json.
2023-06-03 15:25:35,521 INFO - Run target protein preparation by using the Chemoinformatics plugin with 8gcy_preparation_setting.json
2023-06-03 15:25:35,824 INFO - Shutdown when task is done.
2023-06-03 15:25:36,027 INFO Using memory of type: LocalCache
2023-06-03 15:25:36,224 INFO Using Browser: chrome
2023-06-03 15:25:53,953 INFO  THOUGHTS: I need to make protein apo structure and ligand-only PDBs from 8gcy.pdb.
2023-06-03 15:25:54,313 INFO REASONING: The first step is to separate the protein and ligand from the given PDB file.
2023-06-03 15:25:54,651 INFO PLAN: 
2023-06-03 15:25:54,727 INFO -  Use make_apo_protein_pdb to create protein apo structure
2023-06-03 15:25:54,943 INFO -  Use make_only_ligand_compound_pdb to create ligand-only PDB
2023-06-03 15:25:55,137 INFO CRITICISM: None
2023-06-03 15:25:55,231 INFO NEXT ACTION:  COMMAND = make_apo_protein_pdb  ARGUMENTS = {'pdb': '8gcy.pdb'}
2023-06-03 15:25:55,451 INFO  Enter 'y' to authorise command, 'y -N' to run N continuous commands, 's' to run self-feedback commands, 'n' to exit program, or enter feedback for ...
2023-06-03 15:25:55,452 INFO  Asking user via keyboard...
2023-06-03 15:26:09,328 INFO -=-=-=-=-=-=-= COMMAND AUTHORISED BY USER -=-=-=-=-=-=-= 
2023-06-03 15:26:09,785 INFO SYSTEM:  Command make_apo_protein_pdb returned: 8gcy_apo.pdb
2023-06-03 15:26:35,194 INFO  THOUGHTS: Now I need to create the ligand-only PDB.
2023-06-03 15:26:35,469 INFO REASONING: I have created the protein apo structure, and now I need to create the ligand-only PDB to complete the first goal.
2023-06-03 15:26:35,897 INFO PLAN: 
2023-06-03 15:26:35,964 INFO -  Use make_only_ligand_compound_pdb to create ligand-only PDB
2023-06-03 15:26:36,183 INFO CRITICISM: None
2023-06-03 15:26:36,236 INFO NEXT ACTION:  COMMAND = make_only_ligand_compound_pdb  ARGUMENTS = {'pdb': '8gcy.pdb'}
2023-06-03 15:26:36,550 INFO AUTHORISED COMMANDS LEFT: 9 
2023-06-03 15:26:36,803 INFO SYSTEM:  Command make_only_ligand_compound_pdb returned: ['ZN.pdb', 'SO4.pdb', 'Z3N.pdb', 'EDO.pdb', 'UNX.pdb']
2023-06-03 15:26:54,979 INFO  THOUGHTS: I need to read target_preparation_example.jsonc to know how to output setting json to run target protein preparation.
2023-06-03 15:26:55,351 INFO REASONING: Reading the example file will help me understand the format and required parameters for creating the setting json file.
2023-06-03 15:26:55,800 INFO PLAN: 
2023-06-03 15:26:55,883 INFO -  Use read_file to read target_preparation_example.jsonc
2023-06-03 15:26:56,065 INFO CRITICISM: None
2023-06-03 15:26:56,164 INFO NEXT ACTION:  COMMAND = read_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/target_preparation_example.jsonc'}
2023-06-03 15:26:56,381 INFO AUTHORISED COMMANDS LEFT: 8 
2023-06-03 15:26:56,555 INFO SYSTEM:  Command read_file returned: {
2023-06-03 15:27:42,585 INFO  THOUGHTS: I need to update the setting json file to use the two PDBs I just created, protein apo pdb and ligand Z3N.pdb, and output the file as 8gcy_preparation_setting.json.
2023-06-03 15:27:43,050 INFO REASONING: Updating the setting json file with the correct input paths will allow me to run target protein preparation with the created PDBs.
2023-06-03 15:27:43,457 INFO PLAN: 
2023-06-03 15:27:43,530 INFO -  Update the setting json file
2023-06-03 15:27:43,711 INFO -  Write the updated setting json to 8gcy_preparation_setting.json
2023-06-03 15:27:43,911 INFO CRITICISM: None
2023-06-03 15:27:44,019 INFO NEXT ACTION:  COMMAND = write_to_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/8gcy_preparation_setting.json', 'text': '{\n    "target_preparation": {\n        "header": {\n            "logging": {\n                "logfile": "8gcy_target_prep.log"\n            }\n        },\n        "input_path": "8gcy_apo.pdb",\n        "fixer": {\n            "enabled": true,\n            "standardize": true,\n            "remove_heterogens": true,\n            "fix_missing_heavy_atoms": true,\n            "fix_missing_hydrogens": true,\n            "fix_missing_loops": false,\n            "add_water_box": false,\n            "fixed_pdb_path": "8gcy_fixed_target.pdb"\n        },\n        "runs": [\n            {\n                "backend": "AutoDockVina",\n                "output": {\n                    "receptor_path": "8gcy_receptor.pdbqt"\n                },\n                "parameters": {\n                    "pH": 7.4,\n                    "extract_box": {\n                        "reference_ligand_path": "Z3N.pdb",\n                        "reference_ligand_format": "PDB"\n                    }\n                }\n            }\n        ]\n    }\n}'}
2023-06-03 15:27:44,711 INFO AUTHORISED COMMANDS LEFT: 7 
2023-06-03 15:27:44,848 INFO SYSTEM:  Command write_to_file returned: File written to successfully.
2023-06-03 15:28:07,866 INFO  THOUGHTS: I need to run target protein preparation by using the Chemoinformatics plugin with 8gcy_preparation_setting.json.
2023-06-03 15:28:08,193 INFO REASONING: Running target protein preparation will process the input PDBs according to the settings in the json file.
2023-06-03 15:28:08,543 INFO PLAN: 
2023-06-03 15:28:08,605 INFO -  Use target_protein_preparation with 8gcy_preparation_setting.json
2023-06-03 15:28:08,792 INFO CRITICISM: None
2023-06-03 15:28:08,902 INFO NEXT ACTION:  COMMAND = target_protein_preparation  ARGUMENTS = {'setting_json': '8gcy_preparation_setting.json'}
2023-06-03 15:28:09,164 INFO AUTHORISED COMMANDS LEFT: 6 
2023-06-03 15:28:25,052 INFO SYSTEM:  Command target_protein_preparation returned: ('8gcy_receptor.pdbqt', '8gcy_target_prep.log')
2023-06-03 15:28:48,716 INFO  THOUGHTS: I have completed the target protein preparation. Now, I should shut down.
2023-06-03 15:28:49,064 INFO REASONING: All the tasks have been completed, and there is no further action required.
2023-06-03 15:28:49,417 INFO PLAN: 
2023-06-03 15:28:49,465 INFO -  Shutdown
2023-06-03 15:28:49,578 INFO CRITICISM: None
2023-06-03 15:28:49,668 INFO NEXT ACTION:  COMMAND = task_complete  ARGUMENTS = {'reason': 'All tasks completed'}
2023-06-03 15:28:49,927 INFO AUTHORISED COMMANDS LEFT: 5 
2023-06-03 15:28:50,027 INFO Shutting down...
```

### 3. Docking の設定ファイルを作成する
続けて DockStream では `target protein preparation` 同様に Docking を行うために設定ファイルも必要であるため JSON フォーマットの設定ファイルを AutoGPT に作成するよう指示しました。
こちらも私があらかじめ [docking_example.jsonc](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/input/docking_example.jsonc) というファイルを用意しております。

[ai_settings_230603_3.yaml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_3.yaml)

```yaml
ai_goals:
- Read docking_example.jsonc to know how to output setting json to run docking ligand compound and target protein
- Update the setting json file to use the 8gcy_receptor.pdbqt created in target protein preparation and that log file 8gcy_target_prep.log. So you must assign `--center_x`, `--center_y`, `--center_z`, `--size_x`, `--size_y`, `--size_z` parameters in docking setting json with reference to target protein preparation log file. Each size parameter should be rounded up or otherwise margined after subtracting the min value from the max value. Output the file as docking_setting.json.
- Shutdown when task is done.
ai_name: Hoshino AI
ai_role: Medicinal Chemist and Chemoinformatician
api_budget: 3.0
```

docking の設定ファイルで工夫したところは target protein のポケットの指定です。
2 の `target protein preparation` では以下のように[ログ](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/8gcy_target_prep.log)にポケットの座標や大きさが出力されています。

```txt
2023-06-03 15:28:10 - DockStream version used: 1.0.0
2023-06-03 15:28:11 - Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
2023-06-03 15:28:18 - Started preparation run number 0.
2023-06-03 15:28:18 - In specify_cavity
2023-06-03 15:28:19 - In _export_as_pdb2pdbqt
2023-06-03 15:28:19 - In _export_as_pdb2pdbqt argument=['/tmp/tmp7p9jjwnm.pdb', '-opdbqt', '-O8gcy_receptor.pdbqt', '-xr', '-p', 7.4, '--partialcharge', 'gasteiger']
2023-06-03 15:28:24 - Ligand from file Z3N.pdb has the following dimensions:
X coordinates: min=-4.67, max=8.59, mean=2.21
Y coordinates: min=4.02, max=12.84, mean=8.41
Z coordinates: min=13.81, max=23.18, mean=18.57
2023-06-03 15:28:24 - In elif AUTODOCKVINA
2023-06-03 15:28:24 - Wrote AutoDock Vina target to file 8gcy_receptor.pdbqt.
2023-06-03 15:28:24 - Completed target preparation run number 0.
```

こちらの座標を AutoGPT が読みポケットの座標や大きさの指定を自動で判断し設定ファイルでアサインするようにしました。(少し `ai_goals` の文章が長くなってしまっていますが)

これらをまとめたを [ai_settings_230603_3.yaml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_3.yaml) を実行した結果が以下です。 ポケットの大きさを判断し [JSONファイル](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/docking_setting.json) を今回のデータにフィットするように出力し、正常に終了しています。

```txt
2023-06-03 15:51:04,300 INFO SYSTEM:  Command read_file returned: 2023-05-30 13:26:14 - DockStream version used: 1.0.0
2023-06-03 15:52:16,286 INFO GPT4 Only Mode:  ENABLED
2023-06-03 15:52:16,443 INFO Using AI Settings File: ai_settings_230603_3.yaml
2023-06-03 15:52:19,899 INFO  Auto-GPT-Chemoinformatics-plugin: 0.1.0 - This is the Auto-GPT chemoinformatics plugin.
2023-06-03 15:52:20,007 INFO Name : Hoshino AI
2023-06-03 15:52:20,172 INFO Role : Medicinal Chemist and Chemoinformatician
2023-06-03 15:52:20,373 INFO Goals: ['Read docking_example.jsonc to know how to output setting json to run docking ligand compound and target protein', 'Update the setting json file to use the 8gcy_receptor.pdbqt created in target protein preparation and that log file 8gcy_target_prep.log. So you must assign `--center_x`, `--center_y`, `--center_z`, `--size_x`, `--size_y`, `--size_z` parameters in docking setting json with reference to target protein preparation log file. Each size parameter should be rounded up or otherwise margined after subtracting the min value from the max value. Output the file as docking_setting.json.', 'Shutdown when task is done.']
2023-06-03 15:52:21,033 INFO API Budget: $3.0
2023-06-03 15:52:21,143 INFO Hoshino AI has been created with the following details:
2023-06-03 15:52:21,393 INFO Name: Hoshino AI
2023-06-03 15:52:21,504 INFO Role: Medicinal Chemist and Chemoinformatician
2023-06-03 15:52:21,657 INFO Goals: 
2023-06-03 15:52:21,720 INFO - Read docking_example.jsonc to know how to output setting json to run docking ligand compound and target protein
2023-06-03 15:52:22,089 INFO - Update the setting json file to use the 8gcy_receptor.pdbqt created in target protein preparation and that log file 8gcy_target_prep.log. So you must assign `--center_x`, `--center_y`, `--center_z`, `--size_x`, `--size_y`, `--size_z` parameters in docking setting json with reference to target protein preparation log file. Each size parameter should be rounded up or otherwise margined after subtracting the min value from the max value. Output the file as docking_setting.json.
2023-06-03 15:52:22,667 INFO - Shutdown when task is done.
2023-06-03 15:52:22,865 INFO Using memory of type: LocalCache
2023-06-03 15:52:23,007 INFO Using Browser: chrome
2023-06-03 15:52:41,516 INFO  THOUGHTS: I need to read the docking_example.jsonc file to understand how to output the setting json for running docking ligand compound and target protein.
2023-06-03 15:52:41,932 INFO REASONING: Reading the example file will provide me with the necessary information to create the docking setting json.
2023-06-03 15:52:42,274 INFO PLAN: 
2023-06-03 15:52:42,352 INFO -  Read docking_example.jsonc
2023-06-03 15:52:42,501 INFO -  Update the setting json file
2023-06-03 15:52:42,692 INFO -  Output the file as docking_setting.json
2023-06-03 15:52:42,836 INFO -  Shutdown when task is done
2023-06-03 15:52:42,996 INFO CRITICISM: None
2023-06-03 15:52:43,112 INFO NEXT ACTION:  COMMAND = read_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/docking_example.jsonc'}
2023-06-03 15:52:43,361 INFO  Enter 'y' to authorise command, 'y -N' to run N continuous commands, 's' to run self-feedback commands, 'n' to exit program, or enter feedback for ...
2023-06-03 15:52:43,363 INFO  Asking user via keyboard...
2023-06-03 15:53:04,935 INFO -=-=-=-=-=-=-= COMMAND AUTHORISED BY USER -=-=-=-=-=-=-= 
2023-06-03 15:53:05,141 INFO SYSTEM:  Command read_file returned: {
2023-06-03 15:53:23,575 INFO  THOUGHTS: I need to read the 8gcy_target_prep.log file to get the necessary parameters for the docking setting json.
2023-06-03 15:53:24,034 INFO REASONING: Reading the log file will provide me with the information needed to update the docking setting json with the correct parameters.
2023-06-03 15:53:24,449 INFO PLAN: 
2023-06-03 15:53:24,511 INFO -  Read 8gcy_target_prep.log
2023-06-03 15:53:24,659 INFO -  Update the setting json file
2023-06-03 15:53:24,875 INFO -  Output the file as docking_setting.json
2023-06-03 15:53:25,076 INFO -  Shutdown when task is done
2023-06-03 15:53:25,274 INFO CRITICISM: None
2023-06-03 15:53:25,362 INFO NEXT ACTION:  COMMAND = read_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/8gcy_target_prep.log'}
2023-06-03 15:53:25,590 INFO AUTHORISED COMMANDS LEFT: 4 
2023-06-03 15:53:25,735 INFO SYSTEM:  Command read_file returned: 2023-06-03 15:28:10 - DockStream version used: 1.0.0
2023-06-03 15:54:40,169 INFO  THOUGHTS: I will update the docking setting json with the parameters from the 8gcy_target_prep.log file.
2023-06-03 15:54:40,558 INFO REASONING: The log file contains the necessary information to update the docking setting json with the correct parameters.
2023-06-03 15:54:40,889 INFO PLAN: 
2023-06-03 15:54:40,971 INFO -  Update the setting json file
2023-06-03 15:54:41,206 INFO -  Output the file as docking_setting.json
2023-06-03 15:54:41,338 INFO -  Shutdown when task is done
2023-06-03 15:54:41,553 INFO CRITICISM: None
2023-06-03 15:54:41,709 INFO NEXT ACTION:  COMMAND = write_to_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/docking_setting.json', 'text': '{\n    "docking": {\n        "header": {\n            "logging": {\n                "logfile": "/workspace/Auto-GPT/tutorial/ADV_docking.log"\n            }\n        },\n        "ligand_preparation": {\n            "embedding_pools": [\n                {\n                    "pool_id": "RDkit",\n                    "type": "RDkit",\n                    "parameters": {\n                        "prefix_execution": ""\n                    },\n                    "input": {\n                        "standardize_smiles": false,\n                        "type": "console"\n                    },\n                    "output": {\n                        "conformer_path": "/workspace/Auto-GPT/tutorial/ADV_embedded_ligands.sdf",\n                        "format": "sdf"\n                    }\n                }\n            ]\n        },\n        "docking_runs": [\n            {\n                "backend": "AutoDockVina",\n                "run_id": "AutoDockVina",\n                "input_pools": [\n                    "RDkit"\n                ],\n                "parameters": {\n                    "binary_location": "/home/vscode/micromamba/envs/DockStream/bin/",\n                    "parallelization": {\n                        "number_cores": 6\n                    },\n                    "seed": 42,\n                    "receptor_pdbqt_path": [\n                        "/workspace/Auto-GPT/tutorial/8gcy_receptor.pdbqt"\n                    ],\n                    "number_poses": 10,\n                    "search_space": {\n                        "--center_x": 2.21,\n                        "--center_y": 8.41,\n                        "--center_z": 18.57,\n                        "--size_x": 14,\n                        "--size_y": 9,\n                        "--size_z": 10\n                    }\n                },\n                "output": {\n                    "poses": {\n                        "poses_path": "/workspace/Auto-GPT/tutorial/ADV_ligands_docked.sdf",\n                        "overwrite": false\n                    },\n                    "scores": {\n                        "scores_path": "/workspace/Auto-GPT/tutorial/ADV_scores.csv",\n                        "overwrite": false\n                    }\n                }\n            }\n        ]\n    }\n}'}
2023-06-03 15:54:42,397 INFO AUTHORISED COMMANDS LEFT: 3 
2023-06-03 15:54:42,616 INFO SYSTEM:  Command write_to_file returned: File written to successfully.
2023-06-03 15:55:04,193 INFO  THOUGHTS: I have updated the docking setting json file. Now, I will shut down.
2023-06-03 15:55:04,546 INFO REASONING: The task is complete, and there is no further action required.
2023-06-03 15:55:04,833 INFO PLAN: 
2023-06-03 15:55:04,888 INFO -  Shutdown
2023-06-03 15:55:04,974 INFO CRITICISM: None
2023-06-03 15:55:05,089 INFO NEXT ACTION:  COMMAND = task_complete  ARGUMENTS = {'reason': 'Task completed successfully'}
2023-06-03 15:55:05,392 INFO AUTHORISED COMMANDS LEFT: 2 
2023-06-03 15:55:05,525 INFO Shutting down...
```

### 4. REINVNET の設定ファイルを作成する
はじめは既知化合物を参照しつつ、 AutoDocking Vina による docking score を指標とした REINVENT の転移学習を行いたいと考えました。2, 3 同様 REINVENT の設定ファイルを作成します。今回も JSON フォーマットの設定ファイルを AutoGPT に作成するよう指示しました。
こちらも私があらかじめ [reinforcement_learning_setting.jsonc](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/input/reinforcement_learning_setting.jsonc) というファイルを用意しております。
設定ファイルを作成する際、既知構造を参考に新規化学構造を REINVENT で発生させるよう `smiles` という設定項目に 1 にて `CBLB_inhibitors_vsF.sdf` から高活性の化学構造にフィルターした `CBLB_inhibitors_vsF.txt` の化学構造10個をアサインするように指示しました。
自動で `CBLB_inhibitors_vsF.txt` を読み込み10個を選ぶようにしたかったのですが AutoGPT はそのような行動が苦手のようで `ai_goals` の指示の仕方を変えてみてもうまくワークしなかったため、少し残念ですが、化学構造10個を直接 `ai_goals` に記載しました。

[ai_setting_230603_4.yaml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_4.yaml)

```yaml
ai_goals:
- Read reinforcement_learning_setting.jsonc to know how to output setting json.
- Update the reinforcement_learning_setting.jsonc file to use `/workspace/Auto-GPT/autogpt/auto_gpt_workspace/docking_setting.json` to configuration_path. Update the reinforcement_learning_setting.jsonc file to use "C[C@@H](Cc1nncn1C)c1cccc(-c2[nH]nc3c(CN4CC[C@H](F)C4)cc(C(F)(F)F)nc23)c1", "Cn1cnnc1CC1(c2cccc(NC(=O)N3CC4(CCNCC4)c4cccnc43)c2)COC1" and "C[C@H](Cc1nncn1C)c1cccc(NC(=O)N2CCOc3cccnc32)c1" SMILES keys to smiles. Output the file as reinforcement_learning_setting.json.
- Shutdown when task is done.
ai_name: Hoshino AI
ai_role: Medicinal Chemist and Chemoinformatician
api_budget: 3.0
```

こちらを実行したところ AutoGPT は正常に Docking score を指標とする REINVENT の転移学習のための[設定ファイル](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/reinforcement_learning_setting.json)を作成しました。
[ログ](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/logs/activity_230603_4.log)を以下に示します。

```txt
2023-06-03 16:56:22,053 INFO GPT4 Only Mode:  ENABLED
2023-06-03 16:56:22,229 INFO Using AI Settings File: ai_settings_230603_4.yaml
2023-06-03 16:56:26,130 INFO  Auto-GPT-Chemoinformatics-plugin: 0.1.0 - This is the Auto-GPT chemoinformatics plugin.
2023-06-03 16:56:26,242 INFO Name : Hoshino AI
2023-06-03 16:56:26,403 INFO Role : Medicinal Chemist and Chemoinformatician
2023-06-03 16:56:26,578 INFO Goals: ['Read reinforcement_learning_setting.jsonc to know how to output setting json.', 'Update the reinforcement_learning_setting.jsonc file to use `/workspace/Auto-GPT/autogpt/auto_gpt_workspace/docking_setting.json` to configuration_path. Update the reinforcement_learning_setting.jsonc file to use "C[C@@H](Cc1nncn1C)c1cccc(-c2[nH]nc3c(CN4CC[C@H](F)C4)cc(C(F)(F)F)nc23)c1", "Cn1cnnc1CC1(c2cccc(NC(=O)N3CC4(CCNCC4)c4cccnc43)c2)COC1" and "C[C@H](Cc1nncn1C)c1cccc(NC(=O)N2CCOc3cccnc32)c1" SMILES keys to smiles. Output the file as reinforcement_learning_setting.json.', 'Shutdown when task is done.']
2023-06-03 16:56:27,126 INFO API Budget: $3.0
2023-06-03 16:56:27,240 INFO Hoshino AI has been created with the following details:
2023-06-03 16:56:27,505 INFO Name: Hoshino AI
2023-06-03 16:56:27,590 INFO Role: Medicinal Chemist and Chemoinformatician
2023-06-03 16:56:27,763 INFO Goals: 
2023-06-03 16:56:27,816 INFO - Read reinforcement_learning_setting.jsonc to know how to output setting json.
2023-06-03 16:56:28,049 INFO - Update the reinforcement_learning_setting.jsonc file to use `/workspace/Auto-GPT/autogpt/auto_gpt_workspace/docking_setting.json` to configuration_path. Update the reinforcement_learning_setting.jsonc file to use "C[C@@H](Cc1nncn1C)c1cccc(-c2[nH]nc3c(CN4CC[C@H](F)C4)cc(C(F)(F)F)nc23)c1", "Cn1cnnc1CC1(c2cccc(NC(=O)N3CC4(CCNCC4)c4cccnc43)c2)COC1" and "C[C@H](Cc1nncn1C)c1cccc(NC(=O)N2CCOc3cccnc32)c1" SMILES keys to smiles. Output the file as reinforcement_learning_setting.json.
2023-06-03 16:56:28,586 INFO - Shutdown when task is done.
2023-06-03 16:56:28,771 INFO Using memory of type: LocalCache
2023-06-03 16:56:28,938 INFO Using Browser: chrome
2023-06-03 16:56:48,860 INFO  THOUGHTS: I need to read the reinforcement_learning_setting.jsonc file and update it with the required information.
2023-06-03 16:56:49,278 INFO REASONING: The first step is to read the existing file and then update it with the new configuration path and SMILES keys.
2023-06-03 16:56:49,761 INFO PLAN: 
2023-06-03 16:56:49,827 INFO -  Read reinforcement_learning_setting.jsonc
2023-06-03 16:56:49,915 INFO -  Update the file with the new configuration path and SMILES keys
2023-06-03 16:56:50,218 INFO -  Output the updated file as reinforcement_learning_setting.json
2023-06-03 16:56:50,431 INFO -  Shutdown when the task is done
2023-06-03 16:56:50,697 INFO CRITICISM: None
2023-06-03 16:56:50,807 INFO NEXT ACTION:  COMMAND = read_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/reinforcement_learning_setting.jsonc'}
2023-06-03 16:56:51,065 INFO  Enter 'y' to authorise command, 'y -N' to run N continuous commands, 's' to run self-feedback commands, 'n' to exit program, or enter feedback for ...
2023-06-03 16:56:51,068 INFO  Asking user via keyboard...
2023-06-03 16:56:57,301 INFO  Asking user via keyboard...
2023-06-03 16:57:00,458 INFO -=-=-=-=-=-=-= COMMAND AUTHORISED BY USER -=-=-=-=-=-=-= 
2023-06-03 16:57:00,688 INFO SYSTEM:  Command read_file returned: {
2023-06-03 16:58:13,904 INFO  THOUGHTS: I need to update the reinforcement_learning_setting.jsonc file with the new configuration_path and SMILES keys.
2023-06-03 16:58:14,230 INFO REASONING: The task requires me to update the reinforcement_learning_setting.jsonc file with the new configuration_path and SMILES keys.
2023-06-03 16:58:14,540 INFO PLAN: 
2023-06-03 16:58:14,603 INFO -  Update the reinforcement_learning_setting.jsonc
2023-06-03 16:58:14,758 INFO -  Output the updated file as reinforcement_learning_setting.json
2023-06-03 16:58:14,996 INFO -  Shutdown when the task is done
2023-06-03 16:58:15,157 INFO CRITICISM: None
2023-06-03 16:58:15,342 INFO NEXT ACTION:  COMMAND = write_to_file  ARGUMENTS = {'filename': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/reinforcement_learning_setting.json', 'text': '{"logging": {"job_id": "demo", "job_name": "Reinforcement learning demo", "logging_frequency": 5, "logging_path": "/workspace/Auto-GPT/autogpt/auto_gpt_workspace/progress.log", "recipient": "local", "result_folder": "/workspace/Auto-GPT/autogpt/auto_gpt_workspace/results", "sender": "http://0.0.0.1"}, "parameters": {"diversity_filter": {"minscore": 0.4, "minsimilarity": 0.4, "name": "IdenticalMurckoScaffold", "nbmax": 25}, "inception": {"memory_size": 100, "sample_size": 10, "smiles": ["C[C@@H](Cc1nncn1C)c1cccc(-c2[nH]nc3c(CN4CC[C@H](F)C4)cc(C(F)(F)F)nc23)c1", "Cn1cnnc1CC1(c2cccc(NC(=O)N3CC4(CCNCC4)c4cccnc43)c2)COC1", "C[C@H](Cc1nncn1C)c1cccc(NC(=O)N2CCOc3cccnc32)c1"]}, "reinforcement_learning": {"prior": "/workspace/Auto-GPT/ReinventCommunity/notebooks/models/random.prior.new", "agent": "/workspace/Auto-GPT/ReinventCommunity/notebooks/models/random.prior.new", "n_steps": 50, "sigma": 128, "learning_rate": 0.0001, "reset": 0, "reset_score_cutoff": 0.5, "batch_size": 16, "margin_threshold": 50}, "scoring_function": {"name": "custom_sum", "parallel": false, "parameters": [{"component_type": "dockstream", "name": "dockstream", "weight": 1, "specific_parameters": {"configuration_path": "/workspace/Auto-GPT/autogpt/auto_gpt_workspace/docking_setting.json", "docker_script_path": "/workspace/Auto-GPT/DockStream/docker.py", "environment_path": "/home/vscode/micromamba/envs/DockStream/bin/python", "transformation": {"transformation_type": "reverse_sigmoid", "high": -5, "k": 0.25, "low": -12}}}]}, "run_type": "reinforcement_learning", "version": 3, "model_type": "default"}'}
2023-06-03 16:58:15,967 INFO AUTHORISED COMMANDS LEFT: 4 
2023-06-03 16:58:16,125 INFO SYSTEM:  Command write_to_file returned: File written to successfully.
2023-06-03 16:58:34,712 INFO  THOUGHTS: I have updated the reinforcement_learning_setting.json file. Now, I need to shut down.
2023-06-03 16:58:34,983 INFO REASONING: The task is complete, and I need to shut down as instructed.
2023-06-03 16:58:35,289 INFO PLAN: 
2023-06-03 16:58:35,367 INFO -  Shutdown
2023-06-03 16:58:35,480 INFO CRITICISM: None
2023-06-03 16:58:35,587 INFO NEXT ACTION:  COMMAND = task_complete  ARGUMENTS = {'reason': 'Task completed successfully'}
2023-06-03 16:58:35,811 INFO AUTHORISED COMMANDS LEFT: 3 
2023-06-03 16:58:35,944 INFO Shutting down...
```

### 5. Docking score を指標とした REINVNET を実行する
3, 4 の工程で作成した設定ファイルを用いて Docking score を指標とした REINVENT を実行するように AutoGPT に指示しました。

[https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_5.yaml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_5.yaml)

```yaml
ai_goals:
- Run run reinvent with docking by using the Chemoinformatics plugin with /workspace/Auto-GPT/autogpt/auto_gpt_workspace/reinforcement_learning_setting.json
- Shutdown when task is done.
ai_name: Hoshino AI
ai_role: Medicinal Chemist and Chemoinformatician
api_budget: 3.0
```

約15時間で計算が終了しました。[ログ](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/logs/activity_230603_5.log)を以下に示します。

```txt
2023-06-03 17:42:28,924 INFO GPT4 Only Mode:  ENABLED
2023-06-03 17:42:29,084 INFO Using AI Settings File: ai_settings_230603_5.yaml
2023-06-03 17:42:34,498 INFO  Auto-GPT-Chemoinformatics-plugin: 0.1.0 - This is the Auto-GPT chemoinformatics plugin.
2023-06-03 17:42:34,733 INFO Name : Hoshino AI
2023-06-03 17:42:34,896 INFO Role : Medicinal Chemist and Chemoinformatician
2023-06-03 17:42:35,099 INFO Goals: ['Run run reinvent with docking by using the Chemoinformatics plugin with /workspace/Auto-GPT/autogpt/auto_gpt_workspace/reinforcement_learning_setting.json', 'Shutdown when task is done.']
2023-06-03 17:42:35,486 INFO API Budget: $3.0
2023-06-03 17:42:35,571 INFO Hoshino AI has been created with the following details:
2023-06-03 17:42:35,825 INFO Name: Hoshino AI
2023-06-03 17:42:35,928 INFO Role: Medicinal Chemist and Chemoinformatician
2023-06-03 17:42:36,097 INFO Goals: 
2023-06-03 17:42:36,170 INFO - Run run reinvent with docking by using the Chemoinformatics plugin with /workspace/Auto-GPT/autogpt/auto_gpt_workspace/reinforcement_learning_setting.json
2023-06-03 17:42:36,468 INFO - Shutdown when task is done.
2023-06-03 17:42:36,659 INFO Using memory of type: LocalCache
2023-06-03 17:42:36,827 INFO Using Browser: chrome
2023-06-03 17:42:53,884 INFO  THOUGHTS: I need to run reinvent with docking using the Chemoinformatics plugin.
2023-06-03 17:42:54,210 INFO REASONING: The goal is to run reinvent with docking using the provided reinforcement learning setting JSON file.
2023-06-03 17:42:54,640 INFO PLAN: 
2023-06-03 17:42:54,724 INFO -  Run reinvent with docking using the provided JSON file
2023-06-03 17:42:55,033 INFO -  Shutdown when the task is done
2023-06-03 17:42:55,250 INFO CRITICISM: None
2023-06-03 17:42:55,382 INFO NEXT ACTION:  COMMAND = run_reinvent_with_docking  ARGUMENTS = {'setting_json': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/reinforcement_learning_setting.json'}
2023-06-03 17:42:55,626 INFO  Enter 'y' to authorise command, 'y -N' to run N continuous commands, 's' to run self-feedback commands, 'n' to exit program, or enter feedback for ...
2023-06-03 17:42:55,629 INFO  Asking user via keyboard...
2023-06-03 17:43:10,343 INFO -=-=-=-=-=-=-= COMMAND AUTHORISED BY USER -=-=-=-=-=-=-= 
2023-06-04 08:40:58,477 INFO SYSTEM:  Command run_reinvent_with_docking returned: None
2023-06-04 08:41:16,295 INFO  THOUGHTS: I have completed running reinvent with docking using the Chemoinformatics plugin.
2023-06-04 08:41:16,633 INFO REASONING: The goal was to run reinvent with docking using the provided reinforcement learning setting JSON file, and I have done that.
2023-06-04 08:41:17,164 INFO PLAN: 
2023-06-04 08:41:17,243 INFO -  Shutdown since the task is done
2023-06-04 08:41:17,404 INFO CRITICISM: None
2023-06-04 08:41:17,495 INFO NEXT ACTION:  COMMAND = task_complete  ARGUMENTS = {'reason': 'Completed running reinvent with docking using the Chemoinformatics plugin.'}
2023-06-04 08:41:17,879 INFO  Enter 'y' to authorise command, 'y -N' to run N continuous commands, 's' to run self-feedback commands, 'n' to exit program, or enter feedback for ...
2023-06-04 08:41:17,881 INFO  Asking user via keyboard...
2023-06-04 10:49:59,330 INFO -=-=-=-=-=-=-= COMMAND AUTHORISED BY USER -=-=-=-=-=-=-= 
2023-06-04 10:49:59,508 INFO Shutting down...
```

REINVENT にて発生された構造として多様な代表構造が [scaffold_memory.csv](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/scaffold_memory.csv) に保存されました。
以下に一部を示します。

```csv
Step,Scaffold,SMILES,dockstream,raw_dockstream,total_score,ID
49.0,c1ccc2nc(COc3ccc4cc(CCc5nn[nH]n5)ccc4c3)ccc2c1,Fc1ccc2ccc(COc3ccc4cc(CCc5nn[nH]n5)ccc4c3)nc2c1,0.7884759306907654,-10.100000381469727,0.7884759306907654,Reinforcement learning demo_0
44.0,O=C(c1ccc2ccccc2n1)N1CCCCC1,Cc1ccc2nc(C(=O)N3CCC(O)(C(F)(F)F)CC3)ccc2c1,0.711891233921051,-9.600000381469727,0.711891233921051,Reinforcement learning demo_1
30.0,O=C(c1cc2cnoc2s1)N1CCN(S(=O)(=O)c2ccccc2)CC1,Cc1ccc(C)c(S(=O)(=O)N2CCN(C(=O)c3cc4c(C)noc4s3)CC2)c1,0.711891233921051,-9.600000381469727,0.711891233921051,Reinforcement learning demo_2
34.0,O=c1oc2c(c3[nH]c4ccccc4c13)CCc1ccccc1-2,O=c1oc2c(c3[nH]c4ccccc4c13)CCc1ccccc1-2,0.6947353482246399,-9.5,0.6947353482246399,Reinforcement learning demo_3
13.0,N=c1nc[nH]c2c1C(=O)N(c1ccc(C3CCCCC3)cc1)CCO2,N=c1nc[nH]c2c1C(=O)N(c1ccc(C3CCC(CC(N)=O)CC3)cc1)CCO2,0.6947353482246399,-9.5,0.6947353482246399,Reinforcement learning demo_4
33.0,O=C(Cc1ccccc1)NCCc1ccc(-c2ccccc2)cc1,O=C(Cc1ccccc1F)NC(Cc1ccc(-c2cc(F)cc(F)c2F)cc1)C(=O)O,0.6947353482246399,-9.5,0.6947353482246399,Reinforcement learning demo_5
1.0,O=C(c1cc2c(=O)[nH]nc(-c3ccccc3)c2[nH]1)c1c(-c2ccccc2)[nH]c2ccccc12,O=C(c1cc2c(=O)[nH]nc(-c3ccc(F)cc3)c2[nH]1)c1c(-c2ccccc2)[nH]c2ccc(Br)cc12,0.6770213842391968,-9.399999618530273,0.6770213842391968,Reinforcement learning demo_6
11.0,O=C(NCC1COc2ccccc2O1)c1cccc(S(=O)(=O)N2CCOCC2)c1,O=C(NCC1COc2ccccc2O1)c1cccc(S(=O)(=O)N2CCOCC2)c1,0.6587843298912048,-9.300000190734863,0.6587843298912048,Reinforcement learning demo_7
35.0,O=C(N=c1cco[nH]1)c1ccc2c(=O)n3c(nc2c1)CCCCC3,Cc1cc(=NC(=O)c2ccc3c(=O)n4c(nc3c2)CCCCC4)[nH]o1,0.6587843298912048,-9.300000190734863,0.6587843298912048,Reinforcement learning demo_8
(以下略)
```

### 6. `scaffold_memory.csv` 既知化合物に含まれる化学構造を取り除く
5 にて発生された化学構造から既知化学構造を除くように指示しました。
スクリプトを一から作成するのも AutoGPT に指示をしてみたりもしたのですがかなり正解に近いコードを出力するものの、正常に実行できるコードを出力できなかったため、今回も私が Plugin として[関数](https://github.com/yamasakih/Auto-GPT-Plugin-Chemoinformatics/blob/f2e84ceb6d1d908ad3de3e0416e3c4ca6d4a2fd6/src/auto_gpt_plugin/chemoinformatics/chemoinformatics.py#L266)を用意し実行するように指示しました。

[ai_settings_230603_6.yaml](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/ai_settings/ai_settings_230603_6.yaml)

```yaml
ai_goals:
- compare smiles (/workspace/Auto-GPT/autogpt/auto_gpt_workspace/results/scaffold_memory.csv) and with sdf (/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf) by using the Chemoinformatics plugin respectively. Output result to comapre_result.csv.
- Shutdown when task is done.
ai_name: Hoshino AI
ai_role: Medicinal Chemist and Chemoinformatician
api_budget: 3.0
```

こちらで AutoGPT を実行したところ Docking score 上位10化合物を [comapre_result.csv](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/output/compare_result.csv) として出力したので、最終成果物として提出いたしました。
(すいません、ログの保存が正常にできなかったため、後日再度同じ実行を[行ったもの](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/logs/activity_230603_6.log)を掲載しています)


```txt
2023-06-04 16:54:12,508 INFO GPT4 Only Mode:  ENABLED
2023-06-04 16:54:12,650 INFO Using AI Settings File: ai_settings_230603_6.yaml
2023-06-04 16:54:16,405 INFO  Auto-GPT-Chemoinformatics-plugin: 0.1.0 - This is the Auto-GPT chemoinformatics plugin.
2023-06-04 16:54:16,521 INFO Name : Hoshino AI
2023-06-04 16:54:16,644 INFO Role : Medicinal Chemist and Chemoinformatician
2023-06-04 16:54:16,850 INFO Goals: ['compare smiles (/workspace/Auto-GPT/autogpt/auto_gpt_workspace/results/scaffold_memory.csv) and with sdf (/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf) by using the Chemoinformatics plugin respectively. Output result to comapre_result.csv.', 'Shutdown when task is done.']
2023-06-04 16:54:17,305 INFO API Budget: $3.0
2023-06-04 16:54:17,390 INFO Hoshino AI has been created with the following details:
2023-06-04 16:54:17,647 INFO Name: Hoshino AI
2023-06-04 16:54:17,824 INFO Role: Medicinal Chemist and Chemoinformatician
2023-06-04 16:54:17,998 INFO Goals: 
2023-06-04 16:54:18,061 INFO - compare smiles (/workspace/Auto-GPT/autogpt/auto_gpt_workspace/results/scaffold_memory.csv) and with sdf (/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf) by using the Chemoinformatics plugin respectively. Output result to comapre_result.csv.
2023-06-04 16:54:18,424 INFO - Shutdown when task is done.
2023-06-04 16:54:18,649 INFO Using memory of type: LocalCache
2023-06-04 16:54:18,818 INFO Using Browser: chrome
2023-06-04 16:54:45,857 INFO  THOUGHTS: I need to compare smiles in a CSV file with an SDF file and output the result to a new CSV file.
2023-06-04 16:54:46,261 INFO REASONING: The task requires me to compare the smiles and sdf files and output the result. I should use the 'compare_smiles_with_sdf' command to achieve this.
2023-06-04 16:54:46,724 INFO PLAN: 
2023-06-04 16:54:46,806 INFO -  Use the 'compare_smiles_with_sdf' command
2023-06-04 16:54:46,987 INFO -  Provide the required input files and output file
2023-06-04 16:54:47,243 INFO CRITICISM: None
2023-06-04 16:54:47,372 INFO NEXT ACTION:  COMMAND = compare_smiles_with_sdf  ARGUMENTS = {'csv_file': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/results/scaffold_memory.csv', 'sdf_file': '/workspace/Auto-GPT/autogpt/auto_gpt_workspace/CBLB_inhibitors_vsF.sdf', 'output_file': 'compare_result.csv'}
2023-06-04 16:54:47,701 INFO  Enter 'y' to authorise command, 'y -N' to run N continuous commands, 's' to run self-feedback commands, 'n' to exit program, or enter feedback for ...
2023-06-04 16:54:47,702 INFO  Asking user via keyboard...
2023-06-04 16:55:59,479 INFO -=-=-=-=-=-=-= COMMAND AUTHORISED BY USER -=-=-=-=-=-=-= 
2023-06-04 16:56:02,537 INFO SYSTEM:  Command compare_smiles_with_sdf returned: compare_result.csv
2023-06-04 16:56:18,847 INFO  THOUGHTS: I have completed the task of comparing smiles and sdf files and outputting the result to a CSV file.
2023-06-04 16:56:19,346 INFO REASONING: Since the task is complete, I should use the 'task_complete' command to shut down.
2023-06-04 16:56:19,655 INFO PLAN: 
2023-06-04 16:56:19,747 INFO -  Use the 'task_complete' command
2023-06-04 16:56:19,901 INFO -  Provide a reason for shutting down
2023-06-04 16:56:20,114 INFO CRITICISM: None
2023-06-04 16:56:20,225 INFO NEXT ACTION:  COMMAND = task_complete  ARGUMENTS = {'reason': 'Task completed successfully'}
2023-06-04 16:56:20,541 INFO AUTHORISED COMMANDS LEFT: 9 
2023-06-04 16:56:20,672 INFO Shutting down...
```

また、行ったアクションも含めてどのようにするべきか考えてほしいと AutoGPT に指示もしたのですが SMILES 文字列が長いことから ChatGPT にリクエストすることのできる文字数制限以上になってしまい、まったくアクションを行うことができなかったです。(本来ならPAINSフィルターなどを行うべきなのですが残念ですが AutoGPT が提案できなかったということで行っていません)
今後は文字数制限にどう対応するかといったことがこのようなツールを利用する際の課題になると思います。

最後に選出した10化合物の画像を表示しておきます。

![selected_compounds.png](https://github.com/yamasakih/LLM-DD-Challenge/blob/master/img/selected_compounds.png)
## まとめ
AutoGPT を使い標的タンパクに結合しそうな化学構造を選出することができました。

AutoGPT は ChemCrow 同様人間がツールを用意することでそれを武器として使うかどうかを判断したり JSON を出力すると言ったことを判断できていました。
一方で簡単なコードは実装できるものの、多少複雑なコードはまだ100%正常に動作するといったことはできず、現状は人がツールを用意し、 LLM に利用してもらうというアプローチが LLM を最も有効活用できると思いました。

もちろん今後の LLM の発展次第ではツールの用意も LLM 自身が可能になってくる可能性があると思います。そうなると爆発的にできることが増えると期待ができます。
