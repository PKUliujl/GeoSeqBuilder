# GeoSeqBuilder
GeoSeqBuilder, an user-friendly tool package based on deep learning, is used for protein backbone-fixed sequence design and simultaneously gives highly accurate side-chain conformations.

![Alt text](https://github.com/PKUliujl/GeoSeqBuilder/blob/main/image/flow.png)

The dependent packages for GeoSeqBuilder are listed in the [requirement.txt](https://github.com/PKUliujl/GeoSeqBuilder/blob/main/requirement.txt).


INSTALLATION
======================
 Users can download it by `git clone https://github.com/PKUliujl/GeoSeqBuilder.git` (without pretrained model's parameters in the [params](https://github.com/PKUliujl/GeoSeqBuilder/blob/geoseqbuilder/params) directory  due to the large size), alternatively, 
it is available to access [our disk](https://disk.pku.edu.cn:443/link/449F22FE2A06CD29D3C6DB182F4C38C2) with 
Access Code `MBYd` to get these files;


Before using GeoSeqBuilder, please check whether the dependent packages are available in your environment (see [requirement.txt](https://github.com/PKUliujl/GeoSeqBuilder/blob/main/requirement.txt)). If not, using `pip` or `conda` to install them.


USAGE
======================
```
To better use GeoSeqBuilder toolkit for backbone-fixed protein sequence design and side-chain conformation prediction, please add some of
these parameters

optional arguments:
  -h, --help            show this help message and exit
  --purpose {0,1,2}     0 is used only for sequence design, 1 is chosen only for side-chain conformation prediction, and 2 involves
                        performing both tasks simultaneously. default: 0
  --inputPATH INPUTPATH, -iP INPUTPATH
                        The directory path should contain the pdb file
  --chainID CHAINID     The given chain for design or repacking
  --inputfile INPUTFILE, -i INPUTFILE
                        A pdb file under inputPATH, eg. T1024_A.pdb.
  --ST ST               Sampling temperature(scaling factor) belongs to (0,1). Larger values increase the number of interation steps for
                        convegence. Only used when --purpose is 0 or 2. default: 0.1
  --n N                 The number of designed sequences. Only activated when --purpose is 0 or 2. default: 10
  --SM SM               Starting mode with different sequence initializations (how many residues are randomly masked). default = 100(%)
                        masked, i.e. random initialization)
  --Fixed FIXED         Those residues are expected to remain unchanged with names following the format
                        {residuenumber_chainid_resname,residuenumber_chainid_resname}, e.g. 1_A_W,10_A_C Please note that the residuenumber
                        is renumbered begining at 0.
  --outputPATH OUTPUTPATH, -oP OUTPUTPATH
                        the directory path of the outputfile. default: inputPATH
  --verbose {0,1}       Display the intermediate sequences
  --noCYS {0,1}         Do NOT generate Cys as far as possible
  
```

EXAMPLE
=====================
For repacking, 
```
python run_GeoSeqBuilder.py --purpose 1 --inputPATH examples/ --inputfile 3mpc_A.pdb --outputPATH ./ --chainID A   
```

For design,  
```
python run_GeoSeqBuilder.py --purpose 0 --inputPATH examples/ --inputfile 3mpc_A.pdb --outputPATH ./ --chainID A
```
or,
```
python run_GeoSeqBuilder.py --purpose 2 --inputPATH examples/ --inputfile 3mpc_A.pdb --outputPATH ./ --chainID A
```



Reminder: Only the regular pdb format files are accepted as inputs. If many GPUs are available in your machine, please specify the param CUDA_VISIBLE_DEVICES which you use.

Feel free to contact me via email at liujl@stu.pku.edu.cn for other issues.  


CITATION
=====================
If you find GeoSeqBuilder useful in your research, please cite it:
```
@article{liu2024all,
  title={All-atom protein sequence design based on geometric deep learning},
  author={Liu, Jiale and Guo, Zheng and Zhang, Changsheng and Lai, Luhua},
  journal={bioRxiv},
  pages={2024--03},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```


