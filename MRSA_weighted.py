import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle #required to display simulated trees
#require PyQt5 package for graphically displaying trees - not automatically installed

#### =========================================================== ####
## class to handle the global list of genomes across all hospitals ##
class Pool:
  def __init__(self, strains, hinit):
    self.strains = strains
    self.tree = Tree()        #initialise empty tree object
    self.d_var = {}           #dictionary to reference each genome in global tree
    self.d_var["S{0}".format(0)] = self.tree.add_child(
      name="0",
      dist=self.strains[0].t_birth) #add root node to tree
    self.d_var["S{0}".format(0)].add_features(
      t_birth=self.strains[0].t_birth,
      hospital=hinit)         #add attributes to the root node

  def tree(self):
    return(self.tree)

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'Pool(strains = '+str(self.strains)+')'

  def add_new_strain(self, mut_rate, index, t, hosp, k, len_seq):
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
          dist=t - self.strains[index].t_birth)
        self.d_var["S{0}".format(len(self.strains)-1)].add_features(
          hospital=hosp,
          t_birth=t)          #add new genome to global tree with atributes
        self.d_var["S{0}".format(index)] = self.d_var["S{0}".format(index)].add_child(
          name=self.d_var["S{0}".format(index)].name,
          dist=0)             #duplicate partent genome to ensure bifurcating tree in output
        self.d_var["S{0}".format(index)].add_features(hospital=hosp,t_birth=t)
        return len(self.strains)-1
      else:
        return -1             #stop run
    elif num_muts == 0:       #no mutation event so genome observed is identical to parent
      self.d_var["S{0}".format(index)] = self.d_var["S{0}".format(index)].add_child(
        name=self.d_var["S{0}".format(index)].name,
        dist=t - self.strains[index].t_birth) #duplicate partent strain
      self.d_var["S{0}".format(index)].add_features(hospital=hosp,t_birth=t)
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

#### ===================================================================== ####
## class to set properties of each hospital in system and handle local trees ##
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
    return 'Hospital(indices = '+str(self.indices) \
      +', tfin = '+str(self.tfin) \
      +', mut_rate = '+str(self.mut_rate) \
      +', inf_rate = '+str(self.inf_rate) \
      +', tfr_rate = '+str(self.tfr_rate) \
      +')'

  def next_inf_time(self,s,k,weights):
    if self.indices!=[]:      #dividing by s normalises weights to sum to 1
      return [-np.log(np.random.uniform())/(self.inf_rate*w/s) for w in weights] 
    else:                     #if the list of strains in the hospital is empty
      return [np.inf]         # return infinity so not chosen as mimimum time

  def next_tfr_time(self):
    if self.indices!=[]:
      return [-np.log(np.random.uniform())/i if i not in [0] else np.inf for i in self.tfr_rate]
    else:                     #if the list of strains in the hospital is empty
      return [np.inf]         # return infinity so not chosen as mimimum time

  def add_index(self, index, t_obs):
    self.indices.append([index,t_obs])  #add index of genome in pool and time to list in hospital

#### ================================================================ ####
## class to calculate and execute the next event to occur in simulation ##
class System:
  def __init__(self, n, inf_rate, mut_rate, tfr_rate, tfin, hinit, len_seq, alpha):
    self.n = n                #number of hospitals
    self.tfin = tfin          #max number of time units to run simulation for
    self.inf_rate = inf_rate  #rate at which strains mutate
    self.mut_rate = mut_rate  #rate of number of muatations of a strain when it mutates
    self.tfr_rate = tfr_rate  #rate of transfers out of a Hospital
    self.hinit = hinit        #set hospital index that root genome starts in
    strains = [[] for i in range(self.n)] #generate empty lists for each hospital
    strains[self.hinit] = [[0,0]]         #add index of root genome to hinit hospital
    self.len_seq = len_seq    #length of genetic sequence of ones and zeros
    self.alpha = alpha        #weighting coefficient for dependence on time difference
    self.pool = Pool([Strain([],0,0,0,[0]*self.len_seq)],hinit) #add root genome to pool
    self.H = []               #initialise empty list to hold hospital objects
    for i in range(n):
      self.H.append(Hospital(
        strains[i],           #add the lists of genome indices to the corresponding hospital
        self.tfin,
        self.mut_rate,
        self.inf_rate,
        self.tfr_rate[i]))    #add transfer rates out from that hospital
    self.H[hinit].d_hosp["H{0}-{1}".format(hinit,0)] = self.H[hinit].tree.add_child(
      name="0",
      dist=0)                 #add root node to hospital tree

  def __repr__(self):         #set how to display complete output
    return self.__str__()

  def __str__(self):          #set how to display readable output to user
    return 'System(n = '+str(self.n) \
      +', inf_rate = '+str(self.inf_rate) \
      +', mut_rate = '+str(self.mut_rate) \
      +', tfr_rate = '+str(self.tfr_rate) \
      +', tfin = '+str(self.tfin) \
      +', len_seq = '+str(self.len_seq) \
      +', pool = '+str(self.pool) \
      +', H = '+str(self.H) \
      +')'

  def transfers(self):
    transfers = []            #list of transfers that occur in simulation
    global_wt = [[] for i in range(self.n)] #empty list of weights for each hospital
    global_wt[self.hinit] = [1] #weight = 1 for initial parent genome
    ## calculate the time of the next event to occur
    inf_times = [[t for t in h.next_inf_time(1,1,global_wt[i])] for i,h in enumerate(self.H)]
    tfr_times = [[t for t in h.next_tfr_time()] for h in self.H]
    inf_time = min(min(i) for i in inf_times)
    tfr_time = min(min(i) for i in tfr_times)
    min_time = min(inf_time,tfr_time)
    x = 0                     #index for where to mutate genetic sequence from
    while min_time < self.tfin: #run simulation up to tfin time steps
      
      ## infection event occur next
      if inf_time < tfr_time:
        for h,t in enumerate(inf_times):
          if inf_time in t:   #select index of hospital and genome from position of inf_time in list
            hosp_ind = h      #index of the hospital with the minimum time
            hosp = self.H[h]  #hospital object with the current minimum time
            hap_ind = hosp.indices[t.index(inf_time)][0] #select index of genome in pool to mutate
        temp_ind = self.pool.add_new_strain(  
          self.mut_rate,
          hap_ind,
          inf_time,
          hosp_ind,
          x,
          self.len_seq)       #return index of new genome that mutated from the current hospital
        if temp_ind == -1:    #if x + num_mut > len_seq then stop as overflow length of gen_seq
          inf_time = self.tfin
          tfr_time = self.tfin
        else:
          self.H[hosp_ind].add_index(temp_ind,inf_time) #add index of new genome to hospital list
          self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,hap_ind)] = self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,hap_ind)].add_child(
            name=self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,hap_ind)].name,
            dist=0)           #add copy of parent strain as child to tree to make bifurcating
          if hap_ind != temp_ind: #if new genome is not old genome with zero mutations
            self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,temp_ind)] = self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,hap_ind)].add_child(
              name=".".join(str(i) for i in self.pool.strains[temp_ind].parent),
              dist=inf_time - self.H[hosp_ind].indices[inf_times[hosp_ind].index(inf_time)][1])
          global_wt = []      #empty list of global weights
          for i,h in enumerate(self.H):
            global_wt.append([0]*len(h.indices)) #empty list of weights of genomes across all hospitals
            for j,n in enumerate(h.indices):
              global_wt[i][j] = np.exp(self.alpha*(n[1]-inf_time)) #weight = exp(alpha(t_obs - current time)) for each genome
          wt_unlist = [i for s in global_wt for i in s] #reduce level of list (unlist once) to sum over
          s = sum(wt_unlist)  #total sum of weights
          k = len(wt_unlist)  #total number of genomes across all hospitals (not necessarily unique)
          ## calculte new infection times updating with current time
          inf_times = [[inf_time + t for t in h.next_inf_time(s,k,global_wt[i])] for i,h in enumerate(self.H)]
          inf_time = min(min(i) for i in inf_times)
          x = x + self.pool.strains[temp_ind].mut #update index of where to mutate genetic sequence from
          
      ## transfer event occurs next
      elif tfr_time < inf_time:
        for h,t in enumerate(tfr_times):
          if tfr_time in t:
            hosp_ind = h      #index of the hospital with the minimum time
            hosp = self.H[h]  #hospital object with the current minimum time
        weights = [0]*len(hosp.indices)     #empty list of weights of genomes across current hospital
        for i,m in enumerate(hosp.indices):
          weights[i] = np.exp(m[1]-tfr_time)#weight = exp(alpha(t_obs - current time)) for each genome in current hospital
        s_wt = sum(weights)   #sum of weights in current hospital
        norm_w = [float(i)/s_wt for i in weights] #normalised weights
        hap_ind = hosp.indices[np.random.choice(len(hosp.indices),p=norm_w)][0] #select a genome with weighting to transfer
        next_hosp = tfr_times[hosp_ind].index(tfr_time) #index of new hospital to transfer to
        transfers.append(["%.2f" % tfr_time, hap_ind, hosp_ind, next_hosp]) #add transfer to list for output
        if hap_ind in [i[0] for i in self.H[next_hosp].indices]: #if genome is already present in hospital transferred to
          strain_ind = max(i for i, val in enumerate(self.H[next_hosp].indices) if val[0] == hap_ind)
          self.H[next_hosp].d_hosp["H{0}-{1}".format(next_hosp,hap_ind)] = self.H[next_hosp].d_hosp["H{0}-{1}".format(next_hosp,hap_ind)].add_child(
            name=self.H[next_hosp].d_hosp["H{0}-{1}".format(next_hosp,hap_ind)].name,
            dist=tfr_time-self.H[next_hosp].indices[strain_ind][1]) #branch length to most recently observed identical genome
        else:                 #genome transferred is not already in that hospital
          self.H[next_hosp].d_hosp["H{0}-{1}".format(next_hosp,hap_ind)] = self.H[next_hosp].tree.add_child(
            name=self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,hap_ind)].name,
            dist=min_time)    #add new child to root node for transfer in as not directly related to decendents in current tree
        self.H[next_hosp].add_index(hap_ind,tfr_time) #add index and time of event to hospital list
        global_wt = []        #empty list of global weights
        for i,h in enumerate(self.H):
          global_wt.append([0]*len(h.indices)) #empty list of weights of genomes across all hospitals
          for j,n in enumerate(h.indices):
            global_wt[i][j] = np.exp(self.alpha*(n[1]-tfr_time)) #weight = exp(alpha(t_obs - current time)) for each genome
        wt_unlist = [i for s in global_wt for i in s] #reduce level of list (unlist once) to sum over
        s = sum(wt_unlist)    #total sum of weights
        k = len(wt_unlist)    #total number of genomes across all hospitals (not necessarily unique)
        ## calculate the updated time of the next event to occur including just transferred genome
        inf_times = [[tfr_time + t for t in h.next_inf_time(s,k,global_wt[i])] for i,h in enumerate(self.H)]
        tfr_times = [[tfr_time + t for t in h.next_tfr_time()] for h in self.H]
        inf_time = min(min(i) for i in inf_times)
        tfr_time = min(min(i) for i in tfr_times)

      ## if next infection time = next transfer time stop code
      else:
        print("Clashing inf and tfr times. Rerun")
        break
      
      min_time = min(inf_time,tfr_time)

    ## print output of simulation  
    fin_strain = []           #list of all genomes in global pool
    for i,strain in enumerate(self.pool.strains):
      fin_strain.append("{i}, {strain}\n".format(i="%03d" % i,strain=strain))
    fin_strain = "".join(fin_strain)        #print as string
    hosp_strains = []         #list of genomes in each hospital with time of observation
    for i,hosp in enumerate(self.H):
      indices = []
      for j in hosp.indices:
        indices.append([j[0],"%.2f" % j[1]])#convert time to 2dp
      hosp_strains.append("Hospital {i}: {h}\n".format(i="%02d" % i, h=indices))
    hosp_strains = "".join(hosp_strains)    #print as string
    return print(fin_strain,hosp_strains,"Transfers:",*transfers, sep="\n")

np.random.seed(343)           #set seed to reproduce results as in report

n = 3                         #number of hospital nodes in system
tfr_rt = [[0,1,1.5],          #Q rate matrix with 0 on diagonal
          [0.1,0,1],
          [2,1,0]]
sys = System(n,                     
             inf_rate=1.1,    #infection rate
             mut_rate=1,      #mutation rate
             tfr_rate=tfr_rt, #transfer rate matrix
             tfin=10,         #max total time units to run simulation to
             hinit=1,         #set hospital root genome starts in
             len_seq=1000,    #length of the genome sequence where difference occur
             alpha=1)         #weighting coefficient for dependence on time difference
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
