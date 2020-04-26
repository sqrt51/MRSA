import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle #required to display simulated trees
#require PyQt5 package for graphically displaying trees - not automatically installed

#### =========================================================== ####
## class to handle the global list of genomes across all hospitals ##
class Pool:
  def __init__(self, strains):
    self.strains = strains
    self.tree = Tree()        #initialise empty tree object
    self.d_var = {}           #dictionary to reference each genome in global tree
    for i,v in enumerate(strains):
        self.d_var["S{0}".format(i)] = self.tree.add_child(
          name="0",
          dist=self.strains[i].t_birth) #add root node to tree
        self.d_var["S{0}".format(i)].add_feature('t_birth',self.strains[i].t_birth)

  def tree(self):
    return(self.tree)

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'Pool(strains = '+str(self.strains)+')'

  def add_new_strain(self, mut_rate, index, t, k, len_seq):
    num_muts = np.random.poisson(mut_rate*(t - self.strains[index].t_birth))
    parent_strain = self.strains[index]
    if num_muts > 0:
        if k + num_muts <= len_seq:
            parent_strain.num_child += 1
            patient_strain = Strain(
              parent_strain.parent + [parent_strain.num_child],
              0,
              num_muts,
              t,
              parent_strain.gen_seq + [])   #add empty list so stored in memory as new object
            for i in range(k,k+num_muts):   #mutate num_mutations from parent strain
                patient_strain.gen_seq[i] = 1
            self.strains.append(patient_strain)
            self.d_var["S{0}".format(len(self.strains)-1)] = self.d_var["S{0}".format(index)].add_child(
              name=".".join(str(i) for i in patient_strain.parent),
              dist=t - self.strains[index].t_birth) #add new genome to global tree
            self.d_var["S{0}".format(len(self.strains)-1)].add_features(t_birth=t)
            self.d_var["S{0}".format(index)] = self.d_var["S{0}".format(index)].add_child(
              name=self.d_var["S{0}".format(index)].name,
              dist=0)             #duplicate partent genome to ensure bifurcating tree in output
            self.d_var["S{0}".format(index)].add_features(t_birth="t")
            return len(self.strains)-1
        else:
            return -1         #stop run
    elif num_muts == 0:       #no mutation event so genome observed is identical to parent
        return index

#### ========================================= ####
## class to set properties of individual genomes ##
class Strain:
  def __init__(self, parent, num_child, mut, t_birth, gen_seq):
    self.parent = parent
    self.num_child = num_child
    self.mut = mut
    self.t_birth = t_birth
    self.gen_seq = gen_seq

  def __repr__(self):         #set how to display complete output
    return  'Strain(t_birth = '+str("%.2f" % self.t_birth) \
      +', mut = '+str(self.mut) \
      +', num_child = '+str(self.num_child) \
      +', parent = '+str(self.parent) \
      +', gen_seq = '+str(self.gen_seq) \
      +')'

  def __str__(self):          #set how to display readable output to user
    return  '['+str("%.2f" % self.t_birth) \
      +','+str(self.mut) \
      +','+str(self.num_child) \
      +','+str(self.parent) \
      +']'
      ## add this back in penultimate line to get gen_seq in output: +','+str(self.gen_seq) \

#### ============================================= ####
## class to set properties of one hospital in system ##
class Hospital:
  def __init__(self, indices, tfin, mut_rate, inf_rate):
    self.indices = indices
    self.tfin = tfin
    self.mut_rate = mut_rate
    self.inf_rate = inf_rate
    self.tree = Tree()        #initialise empty tree object

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'Hospital(indicies = '+str(self.indices) \
      +', tfin = '+str(self.tfin) \
      +', mut_rate = '+str(self.mut_rate) \
      +', inf_rate = '+str(self.inf_rate) \
      +')'

  def next_inf_time(self,weights):
    return [-np.log(np.random.uniform())/(self.inf_rate*w) for w in weights]

  def add_index(self, index, t_obs):
    self.indices.append([index,t_obs])  #add index of genome in pool and time to list in hospital

#### ================================================================ ####
## class to calculate and execute the next event to occur in simulation ##
class System:
  def __init__(self, inf_rate, mut_rate, tfin, len_seq, alpha):
    self.tfin = tfin          #max number of time units to run simulation for
    self.inf_rate = inf_rate  #rate at which strains mutate
    self.mut_rate = mut_rate  #rate of number of muatations of a strain when it mutates
    strains = [[0,0]]         #add index of root genome
    self.len_seq = len_seq    #length of genetic sequence of ones and zeros
    self.alpha = alpha        #weighting coefficient for dependence on time difference
    self.pool = Pool([Strain([],0,0,0,[0]*self.len_seq)]) #add root genome to pool
    self.H = Hospital(strains, self.tfin, self.mut_rate, self.inf_rate) #one hospital object

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'System(inf_rate = '+str(self.inf_rate) \
           +', mut_rate = '+str(self.mut_rate) \
           +', tfin = '+str(self.tfin) \
           +', len_seq = '+str(self.len_seq) \
           +', pool = '+str(self.pool) \
           +', H = '+str(self.H) \
           +')'

  def mutation(self):
      ## calculate the time of the next event to occur
    inf_times = self.H.next_inf_time([1])
    inf_time = min(inf_times)
    x = 0                     #index for where to mutate genetic sequence from
    while inf_time < self.tfin: #run simulation up to tfin time steps
      hap_ind = self.H.indices[inf_times.index(inf_time)][0] #select index of genome in pool to mutate
      temp_ind = self.pool.add_new_strain(
        self.mut_rate,
        hap_ind,
        inf_time,
        x,
        self.len_seq)         #return index of new genome that mutated from the current hospital
      if temp_ind == -1:      #if x + num_mut > len_seq then stop as overflow length of gen_seq
        inf_time = self.tfin
      else:
        self.H.add_index(temp_ind,inf_time)     #add index of new genome to hospital list
        weights = [0]*len(self.H.indices)       #empty list of global weights
        for i,m in enumerate(self.H.indices):
          weights[i] = np.exp(self.alpha*(m[1]-inf_time)) #weight = exp(alpha(t_obs - current time))
        s = sum(weights)      #total sum of weights
        norm_w = [float(i)/s for i in weights]  #dividing by s normalises weights to sum to 1
        ## calculte new infection times updating with current time
        inf_times = [inf_time + t for t in self.H.next_inf_time(norm_w)]
        inf_time = min(inf_times)
        x = x + self.pool.strains[temp_ind].mut #update index of where to mutate genetic sequence from

    ## print output of simulation 
    fin_strain = []           #list of all genomes in global pool
    for i,strain in enumerate(self.pool.strains):
      fin_strain.append("{i}, {strain}\n".format(i="%03d" % i,strain=strain))
    fin_strain = "".join(fin_strain)      #print as string
    indices = []              #list of genomes in each hospital with time of observation
    for j in self.H.indices:
      indices.append([j[0],"%.2f" % j[1]])#convert time to 2dp
    hosp_strains = []
    hosp_strains.append("Hospital: {hosp}\n".format(hosp=indices))
    hosp_strains = "".join(hosp_strains)  #print as string
    return print(fin_strain,hosp_strains,sep="\n")

# ##### NOTES ON HOW TO CONTROL INPUTS #####
# ## len_seq is the length of the sequence of 1s and 0s
# ## tfin is the time to go up to if the len_seq is not completely used up
# ## mut_rate=inf_rate is the rates for the poisson rv and exponential rv which decides the number of mutations and the time between mutations respectively

np.random.seed(11)           #set seed to reproduce results as in report

sys = System(inf_rate=5,      #infection rate
             mut_rate=1,      #mutation rate
             tfin=10,         #max total time units to run simulation to
             len_seq=700,     #length of the genome sequence where difference occur
             alpha=0.5)       #weighting coefficient for dependence on time difference
sys.mutation()

## ###### UNCOMMENT THIS TO PRINT LIST OF STRAIN GENETIC SEQUENCES ######
##for i,strain in enumerate(sys.pool.strains):
##  print(strain.gen_seq)
## ###### -------------------------------------------------------- ######
##print(sys.pool.tree.get_ascii(show_internal=True,attributes=['name','hospital'])) #Print tree in ascii characters in terminal
##print(sys.pool.tree.write(format=3)) #Print Newick Format

#### format and sytle setting for ete3 tree output ####
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False
ns = NodeStyle()
ns["size"] = 0
for node in sys.pool.tree.traverse():
  node.set_style(ns)
  
## display trees in ete3 tree viewer
sys.pool.tree.show(tree_style=ts)
