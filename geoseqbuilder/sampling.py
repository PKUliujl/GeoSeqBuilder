
from geoseqbuilder import Val_builder,builder,batch_fea,batch_tri,fetchseq,Rotamer_number,AA,AA_reverse
from geoseqbuilder.distributionR.samplingR import samplingr
import torch,os
import numpy as np
import datetime,random
import torch.nn.functional as F
import tensorflow as tf
config = tf.compat.v1.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.1 
session = tf.compat.v1.Session(config=config)


standard_aa_names = {
                   "ALA":0,
                   "CYS":1,
                   "ASP":2,
                   "GLU":3,
                   "PHE":4,
                   "GLY":5,
                   "HIS":6,
                   "ILE":7,
                   "LYS":8,
                   "LEU":9,
                   "MET":10,
                   "ASN":11,
                   "PRO":12,
                   "GLN":13,
                   "ARG":14,
                   "SER":15,
                   "THR":16,
                   "VAL":17,
                   "TRP":18,
                   "TYR":19,
                   }
standard_aa_names_r = dict( [val,key ] for key,val in standard_aa_names.items() )



def samPling(args,models):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    Val_builder(os.path.join( args.inputPATH, args.inputfile), args.chainID, args.inputfile[:-4]+'_Val.pdb')
    seq_alpha,seqnumber = fetchseq( os.path.join( args.inputPATH, args.inputfile), args.chainID )
    protein_graph = batch_fea(  args.inputfile[:-4]+'_Val.pdb', args.chainID , len(seq_alpha) )
    ##features
    y = torch.tensor(seqnumber).to(device)

    x = torch.from_numpy(protein_graph['Node_features']).float().to(device)
    edge_index = torch.from_numpy(protein_graph['edge_index']).T.long().to(device)
    edge_attr = torch.from_numpy(protein_graph['edge_attrs']).float().to(device)

    if args.purpose == 0:
        triang = batch_tri( args.inputfile[:-4]+'_Val.pdb', args.chainID, len(seq_alpha), protein_graph['edge_index'])
        model_se = models[0].to(device).eval()
    elif args.purpose == 1:
        model_sc = models[1]
        modeldis = [ model_sc['model%d'%i].to(device).eval() for i in [1,2,3]]
        modelft = [ model_sc['model_ft%d'%i].to(device).eval() for i in [1,2,3]]
    else:
        triang = batch_tri( args.inputfile[:-4]+'_Val.pdb', args.chainID, len(seq_alpha), protein_graph['edge_index'])
        model_se = models[0].to(device).eval()
        model_sc = models[1]
        modeldis = [ model_sc['model%d'%i].to(device).eval() for i in [1,2,3]]
        modelft = [ model_sc['model_ft%d'%i].to(device).eval() for i in [1,2,3]]
    os.system(f'rm {args.inputfile[:-4]}_Val.pdb') 
    print('\n'+'-'*50)
    with torch.no_grad():
        if args.purpose == 1:
            ## from y_initial
            y_initial = y
            condition = torch.scatter( torch.zeros(y_initial.size(0),20).to(device), 1, y.unsqueeze(1),1 ).float()
            distalP1 = modeldis[0](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = True)
            distalP2 = modeldis[1](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = True)
            distalP3 = modeldis[2](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = True)
            #distalP = (distalP1 + distalP2 + distalP3 )/3
            rotamerP1 = modelft[0](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = distalP1)
            rotamerP2 = modelft[1](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = distalP2)
            rotamerP3 = modelft[2](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = distalP3)
            rotamerP = (rotamerP1 + rotamerP2 + rotamerP3)/3
            ##write rotamer
            rotamerP = rotamerP[:,:48,:].argmax(1)
            rotamers = []
            rotamers_padding = []
            j = 0
            while j < x.size(0):
                for rn in range(  Rotamer_number[ y_initial[j].item() ]) :
                    sampling_interval = rotamerP[j,rn].item()
                    res = standard_aa_names_r[ y_initial[j].item() ]
                    #print(seq[j],res,rn,sampling_interval)
                    angle = samplingr( res, rn, sampling_interval)
                    rotamers.append( str(angle) )
                    rotamers_padding.append( angle)
                for rn in range(  4- Rotamer_number[ y_initial[j].item() ]):
                    rotamers_padding.append(0.)
                j+=1
            rotamers = [ float(i) for i in rotamers]
            print('Side-chain repacking modeling ...')
            print('Writing to ' + os.path.abspath(os.path.join( args.outputPATH ,args.inputfile[:-4]+'_repacked.pdb')))
            builder(os.path.join( args.inputPATH, args.inputfile), args.chainID ,rotamers, os.path.join( args.outputPATH ,args.inputfile[:-4]+'_repacked.pdb'))
            print('Done!')
        elif args.purpose == 0 or args.purpose == 2:
            average_acc = 0
            final_seqs = open(os.path.join( args.outputPATH ,args.inputfile[:-4])+ '_' + str(datetime.date.today() )+'.fa','w')
            for output_iter in range(args.n):
                num = range(20)
                np.random.seed( output_iter)
                if args.SM == 100:
                    y_initial = torch.from_numpy( np.random.choice( num,size = x.size(0), replace=True)  ).long().to(device)
                elif args.SM == 0:
                    y_initial = y
                else:
                    y_ = y.clone()
                    mask_pos = np.random.choice([i for i in range(x.size(0))],size = int(x.size(0)*args.SM/100),replace=False)
                    y_initial = y_
                    y_initial[mask_pos] = 0
                condition = torch.scatter( torch.zeros(y_initial.size(0),20).to(device), 1, y_initial.unsqueeze(1),1 ).float()
                if args.Fixed:
                    fix = args.Fixed.split(',')
                    fixid = [i.split('_') for i in fix]
                    fixid_list = [ int(i) for (i,j,k) in fixid]
                    for (i,j,k) in fixid:
                        if j != args.chainID:
                            raise ValueError('--Fixed, %s dose not mathed the designed chain %s'%(j,args.chainID))
                        resid,fixres = int(i),k
                        condition[resid][:20] = F.one_hot( torch.tensor(AA[fixres],device=device)  ,20)
                aa = ''
                proba = 0
                sim_old_new = []
                seq = []
                if args.ST > 0.5:
                    ks = [3] * int(condition.size(0)*1.5 ) + [1]*int(condition.size(0)/2 )
                else:
                    ks = [3] * int(condition.size(0)/2 ) + [1]*int(condition.size(0)/2 )
                for (i,k) in enumerate(ks):
                    pred1,rotamerP,confidenceP = model_se(x = x, edge_index = edge_index, edge_attr = edge_attr.float(), condition = condition, triangle_nodes = torch.from_numpy(triang['triangleindex']).T.long().to(device), triangle_edges= torch.from_numpy(triang['triangle2line']).T.long().to(device))
                    proba_softmax = F.softmax(pred1/args.ST, 1)
                    non_optimal = ( condition[:,:20].argmax(1) - proba_softmax.argmax(1)).nonzero().squeeze(dim=1).tolist()
                    if args.Fixed:
                        non_optimal = [i for i in non_optimal if i not in fixid_list]
                    muteres = random.sample( non_optimal, min(len(non_optimal), k))
                    step_acc = sum( condition[:,:20].argmax(1) == pred1.argmax(1)).item()/len(y)
                    sim_old_new.append(step_acc)
                    seq.append(''.join( AA_reverse[index.item()] for index in condition[:,:20].argmax(1) ) )
                    if (len(sim_old_new)>10 and sim_old_new[-1] - sim_old_new[-10]<1e-6 and step_acc==1) or len(non_optimal)==0:
                        break
                    for j in muteres:
                        if args.noCYS==1:
                            cys_sample = torch.multinomial(proba_softmax[j],num_samples=1)
                            cys_num = 0
                            while cys_sample.item()==1 and cys_num<5:
                                cys_sample = torch.multinomial(proba_softmax[j],num_samples=1)
                                cys_num += 1
                            condition[j] = F.one_hot( cys_sample, 20).float()
                        else:
                            condition[j] = F.one_hot( torch.multinomial(proba_softmax[j],num_samples=1), 20).float()
                    logP = torch.log( torch.gather(F.softmax(pred1,1),index = condition[:,:20].argmax(1).unsqueeze(1),dim =1) )
                    logits = torch.gather(pred1,dim=1,index=condition[:,:20].argmax(1).unsqueeze(1)).mean().item()
                    if args.verbose == 1:
                        print( seq[-1], 'convergent identity: %.4f,'%step_acc , 'iterative steps: %d,'%i ,'#(unconvergent residues): %d,'%len(non_optimal) , 'logP: %.4f,'%logP.mean().item() , 'logits: %.4f'%logits )
                seq = ''.join( AA_reverse[index.item()] for index in condition[:,:20].argmax(1) )
                acc = sum( condition[:,:20].argmax(1) == y).item()/len(y)
                print('Convergent sequence of Epoch %d,'%output_iter+' '*(3-len(str(output_iter))), seq+',', 'relative to native identity: %.4f,'%acc, 'logP: %.4f,'%(torch.log(F.softmax(pred1.detach(),1).max(1)[0]+1e-10).mean().item()), 'logits: %.4f'%(pred1.detach().max(1)[0].mean().item()) )
                final_seqs.write('>%d'%output_iter +' | convergent identity: %.4f'%step_acc +' | iterative steps: %d'%i + ' | #(unconvergent residues): %d'%len(non_optimal)+ ' | relative to native identity: %.4f'%acc + ' | logP: %.4f'%torch.log(F.softmax(pred1,1).max(1)[0]+1e-10).mean().item() + ' | logits: %.4f'%pred1.max(1)[0].mean().item() + '\n')
                final_seqs.write(seq +'\n')
                average_acc += acc
                if args.purpose == 2:
                    alpha = [AA[i] for i in seq]
                    alpha = torch.tensor(alpha).to(device)
                    condition = torch.scatter( torch.zeros(alpha.size(0),20).to(device), 1, alpha.unsqueeze(1),1 ).float()
                    distalP1 = modeldis[0](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = True)
                    distalP2 = modeldis[1](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = True)
                    distalP3 = modeldis[2](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = True)
                #distalP = (distalP1 + distalP2 + distalP3 )/3
                    rotamerP1 = modelft[0](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = distalP1)
                    rotamerP2 = modelft[1](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = distalP2)
                    rotamerP3 = modelft[2](x = x, edge_index = edge_index, edge_attr = edge_attr, condition=condition, distalP = distalP3)
                    rotamerP = (rotamerP1 + rotamerP2 + rotamerP3)/3
                    ##write rotamer
                    rotamerP = rotamerP[:,:48,:].argmax(1)
                    rotamers = []
                    rotamers_padding = []
                    j = 0
                    while j < x.size(0):
                        for rn in range(  Rotamer_number[ alpha[j].item() ]) :
                            sampling_interval = rotamerP[j,rn].item()
                            res = standard_aa_names_r[ alpha[j].item() ]
                            #print(seq[j],res,rn,sampling_interval)
                            angle = samplingr( res, rn, sampling_interval)
                            rotamers.append( str(angle) )
                            rotamers_padding.append( angle)
                        for rn in range(  4- Rotamer_number[ alpha[j].item() ]):
                            rotamers_padding.append(0.)
                        j+=1
                    rotamers = [ float(i) for i in rotamers]
                    print(f'Side-chain packing modeling for sequence {output_iter} ...')
                    print('Writing all-atom model to ' + os.path.abspath(os.path.join( args.outputPATH ,args.inputfile[:-4]+'_%d.pdb'%output_iter)))
                    builder(os.path.join( args.inputPATH, args.inputfile), args.chainID ,rotamers, os.path.join( args.outputPATH ,args.inputfile[:-4]+'_%d.pdb'%output_iter),seq_tobe_designed=seq )
            final_seqs.close()
            print("Average sequence recovery: {:.4f}".format(average_acc/(output_iter+1) ) )
            print('Writing sequences to '+os.path.abspath(os.path.join( args.outputPATH ,args.inputfile[:-4])+ '_' + str(datetime.date.today() )+'.fa'))
            print('Done!')
            return
        
