
import torch
from geoseqbuilder.common.run_argparse import *
from geoseqbuilder.sampling import samPling

def GeoSeqBuilder(args,models):
    samPling(args,models)
    return

if __name__ == '__main__':
    args = run_inputparameters()
    #device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  
    if args.purpose == 0:
        model_se = torch.load('geoseqbuilder/params/Se.pt')#.to(device)
        GeoSeqBuilder(args,[model_se,''])
    elif args.purpose == 1:
        model_sc = torch.load('geoseqbuilder/params/Sc.pt')#.to(device)
        GeoSeqBuilder(args,['',model_sc])
    else:
        model_se = torch.load('geoseqbuilder/params/Se.pt')
        model_sc = torch.load('geoseqbuilder/params/Sc.pt')
        GeoSeqBuilder(args,[model_se,model_sc])
    print(' '*13+'##############################################\n\
             ##      Thanks for using GeoSeqBuilder      ##\n\
             ##    More details see mdl.ipc.pku.edu.cn   ##\n\
             ##############################################\n')
