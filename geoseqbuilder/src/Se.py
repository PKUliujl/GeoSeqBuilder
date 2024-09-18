
import torch
from torch_geometric.nn import *
import torch.nn.functional as F
from torch_scatter import scatter
from geoseqbuilder.src.bases import GCNconv


class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()

        self.conv1_1 = GCNconv(in_channels = 26 , out_channels=50, edge_length=345)  # 160
        self.BN1_1 = torch.nn.BatchNorm1d(50)

        self.conv2_1 = GCNconv(in_channels = 26 , out_channels=50, edge_length=345)  # 160
        self.BN2_1 = torch.nn.BatchNorm1d(50)
        self.conv2_2 = GCNconv(50,50)
        self.BN2_2 = torch.nn.BatchNorm1d(50)

        self.conv3_1 = GCNconv(in_channels = 26 , out_channels=50, edge_length=345)
        self.BN3_1 = torch.nn.BatchNorm1d(50)
        self.conv3_2 = GCNconv(50,50)
        self.BN3_2 = torch.nn.BatchNorm1d(50)
        self.conv3_3 = GCNconv(50,50)
        self.BN3_3 = torch.nn.BatchNorm1d(50)


        self.conv4_1 = GCNconv(in_channels = 26 , out_channels=50, edge_length=345)  # 160
        self.BN4_1 = torch.nn.BatchNorm1d(50)
        self.conv4_2 = GCNconv(50,50)
        self.BN4_2 = torch.nn.BatchNorm1d(50)
        self.conv4_3 = GCNconv(50,50)
        self.BN4_3 = torch.nn.BatchNorm1d(50)
        self.conv4_4 = GCNconv(50,50)
        self.BN4_4 = torch.nn.BatchNorm1d(50)

        self.conv5_1 = GCNconv(in_channels = 26 , out_channels=50, edge_length=345)  # 160
        self.BN5_1 = torch.nn.BatchNorm1d(50)
        self.conv5_2 = GCNconv(50,50)
        self.BN5_2 = torch.nn.BatchNorm1d(50)
        self.conv5_3 = GCNconv(50,50)
        self.BN5_3 = torch.nn.BatchNorm1d(50)
        self.conv5_4 = GCNconv(50,50)
        self.BN5_4 = torch.nn.BatchNorm1d(50)
        self.conv5_5 = GCNconv(50,50)
        self.BN5_5 = torch.nn.BatchNorm1d(50)


        self.lin = torch.nn.Linear(50*5 + 20+ 10*3 ,20) #380  +8  +20 prior
        self.BN = torch.nn.BatchNorm1d(20)
        self.drop = torch.nn.Dropout(0.5)
        
        self.linR = torch.nn.Linear(50*5 + 20+ 10*3,49*4)
        self.BNR = torch.nn.BatchNorm1d(49*4)
 
        self.linPA = torch.nn.Linear(25 +20+ 10*3 , 2)
        self.linbondL = torch.nn.Linear(50 +20 + 10*3, 1)
        
        self.edge_mlp = torch.nn.Sequential(
                torch.nn.Linear(50*5, 50),
                torch.nn.BatchNorm1d(50),
                torch.nn.ELU(),
                torch.nn.Linear(50, 20*20 ),
                )

        self.edge_mlp2 = torch.nn.Sequential(
                torch.nn.Linear(50*5, 50),
                torch.nn.BatchNorm1d(50),
                torch.nn.ELU(),
                torch.nn.Linear(50, 20*10 ),
                )



    def forward(self,x,edge_index,edge_attr,condition,triangle_nodes, triangle_edges):  # distanc
        #x = torch.cat((x,prior),1)
        #x = torch.cat((x,prior),1)
        x1_1,edge_attr1_1 = self.conv1_1(x,edge_index,edge_attr)
        x1_1 = self.BN1_1(x1_1)
        x1_1 = F.elu(x1_1)

        x2_1,edge_attr2_1 = self.conv2_1(x,edge_index,edge_attr)
        x2_1 = self.BN2_1(x2_1)
        x2_1 = F.elu(x2_1)
        #x2_C = x1+x2   # torch.cat((x1,x2),1)
        #edge_attr2_C = edge_attr + edge_attr2
        x2_2, edge_attr2_2 = self.conv2_2(x2_1,edge_index,edge_attr2_1)
        x2_2 = F.elu( self.BN2_2(x2_2) ) # + x1 # self.lin(x1)
        #x2_C = x2_1+x2_2

        x3_1,edge_attr3_1 = self.conv3_1(x,edge_index,edge_attr)
        x3_1 = self.BN3_1(x3_1)
        x3_1 = F.elu(x3_1)
        #x2_C = x1+x2   # torch.cat((x1,x2),1)
        #edge_attr2_C = edge_attr + edge_attr2
        x3_2, edge_attr3_2 = self.conv3_2(x3_1,edge_index,edge_attr3_1)
        x3_2 = F.elu( self.BN3_2(x3_2) ) # + x1 # self.lin(x1)
        #x3_2 = x3_1+x3_2
        x3_3, edge_attr3_3 = self.conv3_3(x3_2,edge_index,edge_attr3_2)
        x3_3 = F.elu( self.BN3_3(x3_3) ) # + x1 # self.lin(x1)
        #x3_C = x3_1 + x3_2 + x3_3

        x4_1,edge_attr4_1 = self.conv4_1(x,edge_index,edge_attr)
        x4_1 = self.BN4_1(x4_1)
        x4_1 = F.elu(x4_1)
        #x2_C = x1+x2   # torch.cat((x1,x2),1)
        #edge_attr2_C = edge_attr + edge_attr2
        x4_2, edge_attr4_2 = self.conv4_2(x4_1,edge_index,edge_attr4_1)
        x4_2 = F.elu( self.BN4_2(x4_2) ) # + x1 # self.lin(x1)
        #x3_2 = x3_1+x3_2
        x4_3, edge_attr4_3 = self.conv4_3(x4_2,edge_index,edge_attr4_2)
        x4_3 = F.elu( self.BN4_3(x4_3) ) # + x1 # self.lin(x1)
        #x3_C = x3_1 + x3_2 + x3_3
        x4_4, edge_attr4_4 = self.conv4_4(x4_3,edge_index,edge_attr4_3)
        x4_4 = self.BN4_4(x4_4)
        x4_4 = F.elu(x4_4)
        #x4_C = x4_1+x4_2+x4_3+x4_4
        #edge_attr4_C = edge_attr + edge_attr2 + edge_attr3 + edge_attr4

        x5_1,edge_attr5_1 = self.conv5_1(x,edge_index,edge_attr)
        x5_1 = self.BN5_1(x5_1)
        x5_1 = F.elu(x5_1)
        #x2_C = x1+x2   # torch.cat((x1,x2),1)
        #edge_attr2_C = edge_attr + edge_attr2
        x5_2, edge_attr5_2 = self.conv5_2(x5_1,edge_index,edge_attr5_1)
        x5_2 = F.elu( self.BN5_2(x5_2) ) # + x1 # self.lin(x1)
        #x3_2 = x3_1+x3_2
        x5_3, edge_attr5_3 = self.conv5_3(x5_2,edge_index,edge_attr5_2)
        x5_3 = F.elu( self.BN5_3(x5_3) ) # + x1 # self.lin(x1)
        #x3_C = x3_1 + x3_2 + x3_3
        x5_4, edge_attr5_4 = self.conv5_4(x5_3,edge_index,edge_attr5_3)
        x5_4 = F.elu(self.BN5_4(x5_4))
        x5_5,edge_attr5_5 = self.conv5_5(x5_4,edge_index,edge_attr5_4)
        x5_5 = F.elu(self.BN5_5(x5_5))

        #x5_C = x5_1+x5_2+x5_3+x5_4+x5_5
        #edge_attr4_C = edge_attr + edge_attr2 + edge_attr3 + edge_attr4


        row,col = edge_index
        edge_weight = self.edge_mlp( torch.cat((edge_attr1_1, edge_attr2_2, edge_attr3_3, edge_attr4_4, edge_attr5_5),1)  ).view(-1,20,20)
        masked_node_update = torch.matmul( condition[col].unsqueeze(1), edge_weight).squeeze(1)
        condition_out = scatter( masked_node_update,row,dim = 0,reduce = 'sum')
        
        ## triangular interactions
        index_i,index_j,index_k = triangle_nodes
        ed_ij,ed_kj,ed_ik,ed_ki = triangle_edges

        edge_weight = self.edge_mlp2( torch.cat((edge_attr1_1, edge_attr2_2, edge_attr3_3, edge_attr4_4, edge_attr5_5),1)  ).view(-1,20,10)

        c_ij = torch.matmul( condition[index_i].unsqueeze(1),edge_weight[ed_ij] ).squeeze(1)
        c_kj = torch.matmul( condition[index_k].unsqueeze(1),edge_weight[ed_kj] ).squeeze(1)
        c_ik = torch.matmul( condition[index_i].unsqueeze(1),edge_weight[ed_ik] ).squeeze(1)
        c_ki = torch.matmul( condition[index_k].unsqueeze(1),edge_weight[ed_ki] ).squeeze(1)
        triangle_node_update = torch.cat((c_ij,c_kj,(c_ik+c_ki)/2),1)
        condition_triangle_out = scatter( triangle_node_update, index_j ,dim = 0,reduce = 'sum')
        
        #print(x5_C.size(), condition_out.size(),condition_triangle_out.size())
        o = torch.cat((x1_1,x2_2,x3_3,x4_4,x5_5, condition_out,condition_triangle_out),1)  #prior
        
        
        #o =self.drop(o)  # x20_C x10_C

        o1 = self.lin(o) # x10_C
        o1 = self.BN(o1)

        o2 = self.linR(o)
        o2 = self.BNR(o2)
        ''' 
        o3 = torch.sigmoid( self.linPA(o) )*2 -1
        o4 = self.linbondL(o)
        '''
        return o1,o2.view(-1,49,4),o, torch.cat((edge_attr1_1, edge_attr2_2, edge_attr3_3, edge_attr4_4, edge_attr5_5),1)       #, o3, o4




class finetuning_Net(torch.nn.Module):
    def __init__(self):
        super(finetuning_Net, self).__init__()
        
        self.Net = Net()
        #self.Net.eval()

        self.lin1 = torch.nn.Linear( 50*5 + 20+ 10*3 +20, 20 )

        self.confidence = torch.nn.Sequential(
                    torch.nn.Linear(50*5 + 20+ 10*3,50),
                    torch.nn.BatchNorm1d(50),
                    torch.nn.ELU(),
                    torch.nn.Dropout(0.5),
                    torch.nn.Linear(50,1)
                    )
        
        self.edge_mlp = torch.nn.Sequential(
                torch.nn.Linear(50*5, 50),
                torch.nn.BatchNorm1d(50),
                torch.nn.ELU(),
                torch.nn.Linear(50, 20*20 ),
                )

        self.edge_mlp2 = torch.nn.Sequential(
                torch.nn.Linear(50*5, 50),
                torch.nn.BatchNorm1d(50),
                torch.nn.ELU(),
                torch.nn.Linear(50, 20*10 ),
                )

        self.MLP = torch.nn.Sequential(
                torch.nn.Linear(20 + 20 + 10*3+ 49*4, 100 ),
                    torch.nn.BatchNorm1d(100),
                    torch.nn.ELU(),
                    torch.nn.Dropout(0.5),
                    torch.nn.Linear(100, 49*4)
                    )
        
    def forward( self,x,edge_index,edge_attr,condition,triangle_nodes, triangle_edges ):
        seqP, rotamerP,intermediateP,intermediateEDGEP = self.Net( x,edge_index,edge_attr,condition,triangle_nodes, triangle_edges  ) 

        confidence = torch.clamp(self.confidence(intermediateP), max=1.0)
        
        condition = self.lin1( torch.cat( (condition,intermediateP), 1 ) )

        row,col = edge_index
        edge_weight = self.edge_mlp( intermediateEDGEP  ).view(-1,20,20)
        masked_node_update = torch.matmul( condition[col].unsqueeze(1), edge_weight).squeeze(1)
        condition_out = scatter( masked_node_update,row,dim = 0,reduce = 'sum')

        ## triangular interactions
        index_i,index_j,index_k = triangle_nodes
        ed_ij,ed_kj,ed_ik,ed_ki = triangle_edges

        edge_weight = self.edge_mlp2( intermediateEDGEP  ).view(-1,20,10)

        c_ij = torch.matmul( condition[index_i].unsqueeze(1),edge_weight[ed_ij] ).squeeze(1)
        c_kj = torch.matmul( condition[index_k].unsqueeze(1),edge_weight[ed_kj] ).squeeze(1)
        c_ik = torch.matmul( condition[index_i].unsqueeze(1),edge_weight[ed_ik] ).squeeze(1)
        c_ki = torch.matmul( condition[index_k].unsqueeze(1),edge_weight[ed_ki] ).squeeze(1)
        triangle_node_update = torch.cat((c_ij,c_kj,(c_ik+c_ki)/2),1)
        condition_triangle_out = scatter( triangle_node_update, index_j ,dim = 0,reduce = 'sum')

        #print(x5_C.size(), condition_out.size(),condition_triangle_out.size())
        o = torch.cat(( condition, condition_out,condition_triangle_out, rotamerP.view(-1,49*4) ),1)  #prior
        
        o = self.MLP(o)

        return seqP, o.view(-1, 49,4), confidence
    
