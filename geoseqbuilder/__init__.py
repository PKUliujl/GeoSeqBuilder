
#from geoseqbuilder.src.bases import GCNconv
from geoseqbuilder.src import Sc
from geoseqbuilder.src import Se
#from geoseqbuilder.sampling import samPling
from geoseqbuilder.builder.test import builder,Val_builder
from geoseqbuilder.Utils.pdb_processor import fetchseq,batch_fea,batch_tri,Rotamer_number,AA,AA_reverse

sequence_design = Se.finetuning_Net()
sidechains = Sc.Net()

