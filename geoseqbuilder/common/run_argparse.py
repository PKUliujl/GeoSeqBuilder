
import argparse
import os

def run_inputparameters():
    description='To better use GeoSeqBuilder toolkit for backbone-fixed protein sequence design and side-chain conformation prediction, please add some of these parameters'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--purpose',type=int,choices=[0,1,2],help='0 is used only for sequence design, 1 is chosen only for side-chain conformation prediction, and 2 involves performing both tasks simultaneously. default: 0',default=0)
    parser.add_argument('--inputPATH','-iP',type=str,help='The directory path should contain the pdb file', required=True)
    parser.add_argument('--chainID',type=str,help='The given chain for design or repacking',required=True)
    parser.add_argument('--inputfile','-i',type=str,help='A pdb file under inputPATH, eg. T1024_A.pdb.', required=True)
    parser.add_argument('--ST',type=float,help='Sampling temperature(scaling factor) belongs to (0,1). Larger values increase the number of interation steps for convegence. Only used when --purpose is 0 or 2. default: 0.1', default = 0.1 )
    parser.add_argument('--n',type= int,help='The number of  designed sequences. Only activated when --purpose is 0 or 2. default: 10', default = 10)
    parser.add_argument('--SM',type=float,help='Starting mode with different sequence initializations (how many residues are randomly masked). default = 100(%%) masked, i.e. random initialization)', default = 100)
    parser.add_argument('--Fixed',type= str,help='Those residues are expected to remain unchanged with names following the format {residuenumber_chainid_resname,residuenumber_chainid_resname}, e.g. 1_A_W,10_A_C Please note that the residuenumber is renumbered begining at 0.',default = None)
    parser.add_argument('--outputPATH','-oP',type=str,help='the directory path of the outputfile. default: inputPATH')
    parser.add_argument('--verbose',type = int, choices=[0,1],default = 0, help = 'Display the intermediate sequences')
    parser.add_argument('--noCYS',type=int,choices={0,1},help='Do NOT generate Cys as far as possible', default = 0)

    args = parser.parse_args()
    
    if args.purpose != 1:
        if args.Fixed is not None  and args.inputfile is None:
            parser.error('With --Fixed given, only a corresponding pdbfie is required')
    if args.ST < 0 or args.ST >1:
        parser.error('--ST, Sampling temperature(scaling factor) belongs to (0,1)')
    if args.SM <0 or args.SM>100:
        parser.error('--SM, random masking rate for sequence initialization belings to (0,100)')


    '''
    if args.seqfile:
        args.purpose = 1

    if args.outputPATH is None:
        args.outputPATH = args.inputPATH

    if args.outputfile is None and args.purpose==0:
        args.outputfile =  args.pdbname +'_' + args.chainID + '_repacked.pdb'

    if args.sequencesPATH:
        if args.purpose == 0:
            parser.error('sequencesPATH can be specified only for design purposes')
        elif args.outputfile is not None or args.seqfile is not None:
            parser.error('With the sequencesPATH specified, a seqfile and outputfile are not compatible')

    if args.outputfile is None and args.purpose==1:
        args.outputfile =  args.pdbname +'_' + args.chainID + '_design.pdb'
    '''
    return  args


if __name__=='__main__':
    run_inputparameters()
