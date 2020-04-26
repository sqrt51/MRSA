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
        name=".".join(str(j) for j in self.strains[i].parent),
        dist=self.strains[i].t_birth) #add root node to tree with attributes
      self.d_var["S{0}".format(i)].add_features(t_birth=self.strains[i].t_birth,alive='T')

  def tree(self):
    return(self.tree)

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'Pool(strains = '+str(self.strains)+')'

  def add_new_strain(self, mut_rate, index, t, hosp):
    num_muts = np.random.poisson(mut_rate*(t - self.strains[index].t_birth))
    parent_strain = self.strains[index]
    if num_muts > 0:          #new strain observed different from parent add to pool and return index
      parent_strain.num_child += 1
      patient_strain = Strain(
        parent_strain.parent + [parent_strain.num_child],
        0,
        num_muts,
        t)
      self.strains.append(patient_strain)
      self.d_var["S{0}".format(len(self.strains)-1)] = self.d_var["S{0}".format(index)].add_child(
        name=".".join(str(i) for i in patient_strain.parent),
        dist=t - self.strains[index].t_birth)
      self.d_var["S{0}".format(len(self.strains)-1)].add_features(
        hospital=hosp,
        t_birth=t,
        alive='T')          #add new strain to tree with atributes
      self.d_var["S{0}".format(index)] = self.d_var["S{0}".format(index)].add_child(
        name=self.d_var["S{0}".format(index)].name,
        dist=0)             #duplicate partent strain to ensure bifurcating tree in output
      self.d_var["S{0}".format(index)].add_features(hospital=hosp,t_birth="t",alive='T')
      return len(self.strains)-1
    elif num_muts == 0:     #no mutation event so strain observed is identical to parent one
      return index

#### ========================================= ####
## class to set properties of individual genomes ##
class Strain:
  def __init__(self, parent, num_child, mut, t_birth):
    self.parent = parent
    self.num_child = num_child
    self.mut = mut
    self.t_birth = t_birth

  def __repr__(self):         #set how to display complete output
    return  'Strain(t_birth = '+str("%.2f" % self.t_birth) \
      +', mut = '+str(self.mut) \
      +', num_child = '+str(self.num_child) \
      +', parent = '+str(self.parent) \
      +')'

  def __str__(self):          #set how to display readable output to user
    return  '['+str("%.2f" % self.t_birth) \
      +','+str(self.mut) \
      +','+str(self.num_child) \
      +','+str(self.parent) \
      +']'

#### ============================================== ####
## class to set properties of each hospital in system ##
class Hospital:
  def __init__(self, indices, tfin, mut_rate, inf_rate, tfr_rate):
    self.indices = indices
    self.tfin = tfin
    self.mut_rate = mut_rate
    self.inf_rate = inf_rate
    self.tfr_rate = tfr_rate
    self.d_hosp = {}          #dictionary to reference each genome in tree
    self.tree = Tree()        #initialise empty tree object

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'Hospital(indicies = '+str(self.indices) \
      +', tfin = '+str(self.tfin) \
      +', mut_rate = '+str(self.mut_rate) \
      +', inf_rate = '+str(self.inf_rate) \
      +', tfr_rate = '+str(self.tfr_rate) \
      +')'

  def next_inf_time(self):
    return -np.log(np.random.uniform())/self.inf_rate

  def next_tfr_time(self):
    return [-np.log(np.random.uniform())/i if i not in [0] else np.inf for i in self.tfr_rate]

  def add_index(self, index, t, parent):#replace an index in hospital with the new index
    del_ind = np.random.choice(len(self.indices[0]),p=[i/sum(self.indices[1]) for i in self.indices[1]])
    del_str = self.indices[0][del_ind]  #det_str index of strain in pool, del_ind index in hospital
    if self.indices[1][del_ind] > 1:    #more than one of genome in hospital
      self.indices[1][del_ind] -= 1     #number of strain counts from 1, dictionary counts from 0
      self.d_hosp["H{0}_{1}".format(del_str,self.indices[1][del_ind])].add_feature('alive',"F")
    elif self.indices[1][del_ind] == 1: #only one of genome in hospital
      self.d_hosp["H{0}_{1}".format(del_str,0)].add_feature('alive',"F")
      del self.indices[0][del_ind]
      del self.indices[1][del_ind]
    if index in self.indices[0]:        #if index of genome alread in hospital
      ind = self.indices[0].index(index)#index is of strain in pool, ind is index in hospital
      self.d_hosp["H{0}_{1}".format(index,self.indices[1][ind])] = self.d_hosp["H{0}_{1}".format(index,0)].add_child(
        name=self.d_hosp["H{0}_{1}".format(index,0)].name,
        dist=t-self.d_hosp["H{0}_{1}".format(index,0)].dist)
      self.indices[1][ind] += 1         #number of strain counts from 1, dictionary counts from 0
    else:                               #if index of genome new to hospital
      self.indices[0].append(index)
      self.indices[1].append(1)
      self.d_hosp["H{0}_{1}".format(index,0)] = self.tree.add_child(
        name=".".join(str(j) for j in parent),
        dist=t)
    return del_str                      #index of stain that one copy of was remvoed

#### ================================================================ ####
## class to calculate and execute the next event to occur in simulation ##
class System:
  def __init__(self, n, inf_rate, mut_rate, tfr_rate, tfin, pool, strains):
    self.n = n                #number of hospitals
    self.tfin = tfin          #max number of time units to run simulation for
    self.inf_rate = inf_rate  #rate at which strains mutate
    self.mut_rate = mut_rate  #rate of number of muatations of a strain when it mutates
    self.tfr_rate = tfr_rate  #rate of transfers out of a Hospital
    self.pool = Pool(pool)    #initialise pool specified in function variables
    self.H = []               #initialise empty list to hold hospital objects
    for i in range(n):
      self.H.append(Hospital(
        strains[i],           #add the lists of strain indices to the corresponding hospital
        self.tfin,
        self.mut_rate,
        self.inf_rate,
        self.tfr_rate[i]))    #add transfer rates out from that hospital
    for i,h in enumerate(self.H):
      for j,s in enumerate(h.indices[0]):
        prop = h.indices[1][j]
        while prop >=1:       #add each copy of strain to corresponding hospital tree
          h.d_hosp["H{0}_{1}".format(s,h.indices[1][j]-prop)] = h.tree.add_child(
            name=".".join(str(j) for j in self.pool.strains[s].parent),
            dist=self.pool.strains[s].t_birth)
          prop -= 1

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'System(n = '+str(self.n) \
      +', inf_rate = '+str(self.inf_rate) \
      +', mut_rate = '+str(self.mut_rate) \
      +', tfr_rate = '+str(self.tfr_rate) \
      +', tfin = '+str(self.tfin) \
      +', pool = '+str(self.pool) \
      +', strains = '+str(self.strains) \
      +', H = '+str(self.H) \
      +')'

  def transfers(self):
    transfers = []            #list of transfers that occur in simulation
    ## calculate the time of the next event to occur
    inf_times = [t.next_inf_time() for t in self.H]
    tfr_times = [t.next_tfr_time() for t in self.H]
    inf_time = min(inf_times)
    tfr_time = min(min(i) for i in tfr_times)
    min_time = min(inf_time,tfr_time)
    while min_time < self.tfin: #run simulation up to tfin time steps

      ## infection event occur next      
      if inf_time < tfr_time:
        hosp_ind = inf_times.index(inf_time)#select the index of the current min of times
        next_hosp = self.H[hosp_ind]        #select the hospital corresponding to that index
        hap_ind = next_hosp.indices[0][np.random.choice(
          len(next_hosp.indices[0]),
          p=[i/sum(next_hosp.indices[1]) for i in next_hosp.indices[1]]
          )] #randomly with weighting select an index of strain from hospital with min time
        temp_ind = self.pool.add_new_strain(#return index of new strain that is the child of hap_ind
          self.mut_rate,
          hap_ind,
          inf_time,
          hosp_ind) 
        del_ind = self.H[hosp_ind].add_index(#add new strain to list and delete uniformly a random strain
          temp_ind,
          inf_time,
          self.pool.strains[hap_ind].parent)
        ## calculte new infection and transfer times updating with current time
        inf_times = [inf_time + t.next_inf_time() for t in self.H]
        tfr_times = [[inf_time + i for i in t.next_tfr_time()] for t in self.H]
        tfr_time = min(min(i) for i in tfr_times)
        inf_time = min(inf_times)

      ## transfer event occurs next
      elif tfr_time < inf_time:
        hosp_ind = [min(i) for i in tfr_times].index(tfr_time) #select the index of the current min of times
        next_hosp = self.H[hosp_ind]        #select the hospital corresponding to that index
        hap_ind = next_hosp.indices[0][np.random.choice(
          len(next_hosp.indices[0]),
          p=[i/sum(next_hosp.indices[1]) for i in next_hosp.indices[1]]
          )] #randomly with weighting select an index of strain from hospital with min time
        new_ind = tfr_times[hosp_ind].index(tfr_time) #index of hospital to transfer to
        transfers.append(["%.2f" % tfr_time, hap_ind, hosp_ind, new_ind]) #add transfer to list for output
        del_ind = self.H[new_ind].add_index(#add new strain to list and delete uniformly a random strain
          hap_ind,
          tfr_time,
          self.pool.strains[hap_ind].parent)
        ## calculte new infection and transfer times updating with current time
        tfr_times = [[tfr_time + i for i in t.next_tfr_time()] for t in self.H]
        inf_times = [tfr_time + t.next_inf_time() for t in self.H]
        inf_time = min(inf_times)
        tfr_time = min(min(i) for i in tfr_times)

      ## if next infection time = next transfer time stop code       
      else:
        print("Clashing inf and tfr times")
        break
      
      min_time = min(inf_time,tfr_time)

    ## print output of simulation        
    fin_strain = []           #list of all genomes in global pool
    for i,strain in enumerate(self.pool.strains):
      fin_strain.append("{i}, {strain}\n".format(i="%03d" % i,strain=strain))
    fin_strain = "".join(fin_strain)        #print as string
    hosp_strains = []         #list of genomes in each hospital with time of observation
    for i,hosp in enumerate(self.H):
      hosp_strains.append("Hospital {i}: {hosp}\n".format(i="%02d" % i, hosp=hosp.indices))
    hosp_strains = "".join(hosp_strains)    #print as string
    pool_strains = []
    for i,strain in enumerate(self.pool.strains): #iterate over all strains to get total still alive
      tot = 0
      for H in self.H:
        if i in H.indices[0]:
          prop = H.indices[1][H.indices[0].index(i)] #number of that strain
          for j in range(prop):#extend branch length of alive strains to tfin
            H.d_hosp["H{0}_{1}".format(i,j)].add_child(
              name=H.d_hosp["H{0}_{1}".format(i,j)].name,
              dist=self.tfin - H.d_hosp["H{0}_{1}".format(i,j)].get_distance(H.tree))
          tot += prop         #store running total
      pool_strains.append(tot)
    for i,num in enumerate(pool_strains):
      if num >= 1:            #extend branch length of alive strains to tfin
        self.pool.d_var["S{0}".format(i)].add_child(
          name=self.pool.d_var["S{0}".format(i)].name,
          dist=self.tfin - self.pool.strains[i].t_birth)
      elif num == 0:          #update attribut for strains no longer in global pool
        self.pool.d_var["S{0}".format(i)].add_feature('alive',"F")
    return print(fin_strain,hosp_strains,pool_strains,"","Transfers:",*transfers, sep="\n")

np.random.seed(3)             #set seed to reproduce results as in report

n = 3                         #number of hospital nodes in system
tfr_rt = [[0,1,1.5],          #Q rate matrix with 0 on diagonal
          [0.1,0,1],
          [2,1,0]]
pool = [Strain(["A"],0,0,0),Strain(["B"],0,0,0),Strain(["C"],0,0,0)] #pool of all strains at initialisation
strains = [[[0,1],[2,3]],[[0],[5]],[[2],[5]]] #list of length n where [strain index],[number of each strain]
sys = System(n,
             inf_rate=1.5,    #infection rate
             mut_rate=1,      #mutation rate
             tfr_rate=tfr_rt, #transfer rate matrix
             tfin=5,          #max total time units to run simulation to
             pool=pool,       #define pool at system initialisation
             strains=strains) #define list of number of each strain are in each hospiatl
sys.transfers()

#### format and sytle setting for ete3 tree output ####
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False

## set branch style setting to use for different hospital attributes
colour=["red","blue","green","black"]
nstyle = [NodeStyle(),NodeStyle(),NodeStyle(),NodeStyle()]
for i,val in enumerate(nstyle):
  nstyle[i]["fgcolor"] = colour[i]
  nstyle[i]["size"] = 0
  nstyle[i]["vt_line_color"] = "black"
  nstyle[i]["hz_line_color"] = colour[i]
  nstyle[i]["vt_line_width"] = 1
  nstyle[i]["hz_line_width"] = 1
  nstyle[i]["vt_line_type"] = 0
  nstyle[i]["hz_line_type"] = i%3     # 0 solid, 1 dashed, 2 dotted

for node in sys.pool.tree.traverse(): #set basic style to all nodes
  node.set_style(nstyle[3])
for i in range(n):                    #set different styles for each hospital
  for node in sys.pool.tree.iter_search_nodes(hospital = i):
    node.set_style(nstyle[i])
for i in range(3):                    #set the basic style to each of the local trees
    for node in sys.H[i].tree.traverse():
        node.set_style(nstyle[3])

## display trees in ete3 tree viewer
sys.pool.tree.show(tree_style=ts)     #global phylogenetic tree
for i in range(n):                    #local phylogenetic trees for each hospital
    sys.H[i].tree.show(tree_style=ts)
    
