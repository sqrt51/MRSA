import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle #also need PyQt5 package istalled for graphically displaying trees which is not automatically installed as a dependency

class Pool:
  def __init__(self, strains):
    self.strains = strains
    self.tree = Tree()
    self.d_var = {}
    for i,v in enumerate(strains):
        self.d_var["S{0}".format(i)] = self.tree.add_child(name=".".join(str(j) for j in self.strains[i].parent),dist=self.strains[i].t_birth)
        self.d_var["S{0}".format(i)].add_features(t_birth=self.strains[i].t_birth, alive='T')

  def tree(self):
    return(self.tree)

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'Pool(strains = '+str(self.strains)+')'

  def add_new_strain(self, mut_rate, index, t, hosp, k, len_seq):
    num_muts = np.random.poisson(mut_rate*(t - self.strains[index].t_birth))
    parent_strain = self.strains[index]
    if num_muts > 0:
        if k + num_muts <= len_seq:
            parent_strain.num_child += 1
            patient_strain = Strain(parent_strain.parent + [parent_strain.num_child], 0, num_muts, t, parent_strain.gen_seq + []) #add empty list so python dosen't edit parent gen_seq and generate copy that refers to diferent object to source
            for i in range(k,k+num_muts): #mutate num_mutations from parent strain
                patient_strain.gen_seq[i] = 1
            self.strains.append(patient_strain)
            self.d_var["S{0}".format(len(self.strains)-1)] = self.d_var["S{0}".format(index)].add_child(name=".".join(str(i) for i in patient_strain.parent),dist=t - self.strains[index].t_birth)
            self.d_var["S{0}".format(len(self.strains)-1)].add_features(hospital=hosp, t_birth=t, alive='T')
            self.d_var["S{0}".format(index)] = self.d_var["S{0}".format(index)].add_child(name=self.d_var["S{0}".format(index)].name,dist=0) #duplicate partent strain to ensure bifurcating tree in output
            self.d_var["S{0}".format(index)].add_features(hospital=hosp,t_birth=t,alive='T')
            return len(self.strains)-1
        else:
            return -1 #stop run
    elif num_muts == 0: #no mutation event so strain observed is identical to parent one
        return index

  def del_strain(self, index, time): ##unused at the moment in the script
    self.d_var["S{0}".format(index)].add_feature('alive',"F")
    self.d_var["S{0}".format(index)].add_child(name="X",dist=time - self.d_var["S{0}".format(index)].dist)

class Strain:
  def __init__(self, parent, num_child, mut, t_birth, gen_seq):
    self.parent = parent
    self.num_child = num_child
    self.mut = mut
    self.t_birth = t_birth
    self.gen_seq = gen_seq

  def __repr__(self):
    return  'Strain(t_birth = '+str(self.t_birth) \
      +', mut = '+str(self.mut) \
      +', num_child = '+str(self.num_child) \
      +', parent = '+str(self.parent) \
      +', gen_seq = '+str(self.gen_seq) \
      +')'

  def __str__(self):
    return  '['+str("%.2f" % self.t_birth) \
      +','+str(self.mut) \
      +','+str(self.num_child) \
      +','+str(self.parent) \
      +']'


class Hospital:
  def __init__(self, indices, tfin, mut_rate, inf_rate, tfr_rate):
    self.indices = indices
    self.tfin = tfin
    self.mut_rate = mut_rate
    self.inf_rate = inf_rate
    self.tfr_rate = tfr_rate
    self.d_hosp = {}
    self.tree = Tree()

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'Hospital(indicies = '+str(self.indices)+', tfin = '+str(self.tfin)+', mut_rate = '+str(self.mut_rate)+', inf_rate = '+str(self.inf_rate)+', tfr_rate = '+str(self.tfr_rate)+')'

  def next_inf_time(self,s,k,weights):
    if self.indices!=[]:
        return [-np.log(np.random.uniform())/(k*self.inf_rate*w/s) for w in weights] # dividing by s normalises weights to sum to 1
    else:
        return [np.inf]

  def next_tfr_time(self):
    if self.indices!=[]:
        return [-np.log(np.random.uniform())/i if i not in [0] else np.inf for i in self.tfr_rate]
    else:
        return [np.inf]

  # def add_inf(self, t, pool):
    # index = np.random.choice(len(self.indices[0]),p=[i/sum(self.indices[1]) for i in self.indices[1]])
    # index = self.indices[0][index]
    # add_inf = pool.add_new_strain(self.mut_rate, index, t)
    # del_index = np.random.choice(len(self.indices[0]),p=[i/sum(self.indices[1]) for i in self.indices[1]])
    # if self.indices[1][del_index] > 1:
    #     self.indices[1][del_index] -= 1
    # elif self.indices[1][del_index] == 1:
    #     del self.indices[0][del_index]
    #     del self.indices[1][del_index]
    # self.indices[0].append(add_inf)
    # self.indices[1].append(1)

  def add_index(self, index, t_obs):
    #del_ind = np.random.choice(len(self.indices))
    #del_str = self.indices[del_ind]
    self.indices.append([index,t_obs])
    #return del_str


class System:
  def __init__(self, n, inf_rate, mut_rate, tfr_rate, tfin, hinit, len_seq):
    self.n = n #number of hospitals
    self.tfin = tfin
    self.inf_rate = inf_rate #rate at which strains mutate
    self.mut_rate = mut_rate #rate of number of muatations of a strain when it mutates
    self.tfr_rate = tfr_rate #rate of transfers out of a Hospital
    #self.p = p #threshold probability for transfer
    self.hinit = hinit
    strains = [[] for i in range(self.n)]
    strains[self.hinit] = [[0,0]]
    #self.k = k #current global counter value for where in gen_seq to start mutating values
    self.len_seq = len_seq #length of genetic sequence of ones and zeros
    self.pool = Pool([Strain([],0,0,0,[0]*self.len_seq)])
    self.H = []
    for i in range(n):
      self.H.append(Hospital(strains[i], self.tfin, self.mut_rate, self.inf_rate, self.tfr_rate[i]))
    for i,h in enumerate(self.H):
        for j,s in enumerate(h.indices):
            h.d_hosp["H{0}-{1}".format(i,s[0])] = h.tree.add_child(name=0,dist=0)

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'System(n = '+str(self.n)+', inf_rate = '+str(self.inf_rate)+', mut_rate = '+str(self.mut_rate)+', tfr_rate = '+str(self.tfr_rate)+', tfin = '+str(self.tfin)+', len_seq = '+str(self.len_seq)+', pool = '+str(self.pool)+', H = '+str(self.H)+')'
    #', strains = '+str(self.strains)+

  def transfers(self):
    transfers = []
    global_wt = []
    for i,h in enumerate(self.H):
      global_wt.append([0]*len(h.indices)) # generate empty list of weights for each strain in each hospital
      for j,n in enumerate(h.indices):
        global_wt[i][j] = np.exp(n[1]) # fill weights list with exp(t_obs - current time) for each strain in list
    wt_unlist = [i for s in global_wt for i in s] # unlist list of list
    s = sum(wt_unlist) # reduce level of list (unlist once) to sum over
    k = len(wt_unlist) # total number of strains across all hospitals (not necessarily unique)
    inf_times = [[t for t in h.next_inf_time(s,k,global_wt[i])] for i,h in enumerate(self.H)]
    tfr_times = [[t for t in h.next_tfr_time()] for h in self.H]
    inf_time = min(min(i) for i in inf_times)
    tfr_time = min(min(i) for i in tfr_times)
    min_time = min(inf_time,tfr_time)
    x = 0
    while min_time < self.tfin:
##        print(inf_times)
##        print(tfr_times)
##        print(inf_time)
##        print(tfr_time)
        if inf_time < tfr_time:
            for h,t in enumerate(inf_times):
              if inf_time in t:
                hosp_ind = h # index of the hospital with the minimum time
                hosp = self.H[h]  # select the hospital with the current min of times
                hap_ind = hosp.indices[t.index(inf_time)][0] #select from the indicies attribute of Hospital corresponding to the index of a strain in pool
            temp_ind = self.pool.add_new_strain(self.mut_rate, hap_ind, inf_time, hosp_ind, x, self.len_seq) #return index of new strain that is the child from one of the strains in the current incidence of H
##            print(hosp_ind)
##            print(hosp)
##            print(hap_ind)
##            print(temp_ind)
            if temp_ind == -1: #if k + num_mut > len_seq then stop
                inf_time = self.tfin
                tfr_time = self.tfin
            else:
                del_ind = self.H[hosp_ind].add_index(temp_ind,inf_time) #add that strain to the list of strains in that H and deletes a random strain from the list uniformly
                self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,temp_ind)] = self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,hap_ind)].add_child(name=".".join(str(i) for i in self.pool.strains[temp_ind].parent),dist=inf_time - self.H[hosp_ind].indices[inf_times[hosp_ind].index(inf_time)][1])
                # self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,del_ind)] = self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,del_ind)].add_child(name="X",dist=t - self.pool.strains[hap_ind].t_birth)
                # self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,del_ind)].add_feature("alive","F")
                # self.H[hosp_ind].d_hosp["H{0}-{1}".format(hosp_ind,temp_ind)] = self.tree.add_child(name=".".join(str(i) for i in self.pool.strains[temp_ind].parent),dist=min_time - self.pool.strains[hap_ind].t_birth)
                # if del_ind >= 0:
                #   self.pool.del_strain(del_ind, inf_time)
                #new_ind = np.random.choice(range(len(self.H)),p=self.p[hosp_ind]) #randomly pick a new hospital index with weighted probability where a transfer to the same hospital is used as a convention for no transfer taking place
                # if new_ind != hosp_ind:
                #     transfers.append([inf_time, hosp_ind, new_ind])
                #     self.H[new_ind].add_index(temp_ind)
                global_wt = []
                for i,h in enumerate(self.H):
                  global_wt.append([0]*len(h.indices)) # generate empty list of weights for each strain in each hospital
                  for j,n in enumerate(h.indices):
                    global_wt[i][j] = np.exp(n[1]-inf_time) # fill weights list with exp(t_obs - current time) for each strain in list
##                print(global_wt)
                wt_unlist = [i for s in global_wt for i in s] # unlist list of list
                s = sum(wt_unlist) # reduce level of list (unlist once) to sum over
                k = len(wt_unlist) # total number of strains across all hospitals (not necessarily unique)
                inf_times = [[inf_time + t for t in h.next_inf_time(s,k,global_wt[i])] for i,h in enumerate(self.H)]
                inf_time = min(min(i) for i in inf_times)
                x = x + self.pool.strains[temp_ind].mut
        elif tfr_time < inf_time:
            for h,t in enumerate(tfr_times):
              if tfr_time in t:
                hosp_ind = h # index of the hospital with the minimum time
                hosp = self.H[h]  # select the hospital with the current min of times
            weights = [0]*len(hosp.indices) # generate empty list of weights for each strain in the hospital
            for i,m in enumerate(hosp.indices):
                weights[i] = np.exp(m[1]-tfr_time) # fill weights list with exp(t_obs - current time) for each strain in list
            s_wt = sum(weights)
            norm_w = [float(i)/s_wt for i in weights]
##            print(norm_w)
            hap_ind = hosp.indices[np.random.choice(len(hosp.indices),p=norm_w)][0] 
            next_hosp = tfr_times[hosp_ind].index(tfr_time) # index of new hospital to transfer to
            transfers.append(["%.2f" % tfr_time, hap_ind, hosp_ind, next_hosp])
            if hap_ind in [i[0] for i in self.H[next_hosp].indices]:
                strain_ind = max(i for i, val in enumerate(self.H[next_hosp].indices) if val[0] == hap_ind)
                self.H[next_hosp].d_hosp["H{0}-{1}".format(next_hosp,hap_ind)] = self.H[next_hosp].d_hosp["H{0}-{1}".format(next_hosp,hap_ind)].add_child(name=".".join(str(i) for i in self.pool.strains[hap_ind].parent),dist=tfr_time-self.H[next_hosp].indices[strain_ind][1])
            else:
                self.H[next_hosp].d_hosp["H{0}-{1}".format(next_hosp,hap_ind)] = self.H[next_hosp].tree.add_child(name=".".join(str(i) for i in self.pool.strains[hap_ind].parent),dist=min_time)
            self.H[next_hosp].add_index(hap_ind,tfr_time)
            # self.H[new_ind].d_hosp["H{0}-{1}".format(new_ind,del_ind)] = self.H[new_ind].d_var["H{0}-{1}".format(new_ind,del_ind)].add_child(name="X",dist=t - self.pool.strains[hap_ind].t_birth)
            # self.H[new_ind].d_hosp["H{0}-{1}".format(new_ind,del_ind)].add_feature("alive","F")
            # if del_ind >= 0:
            # self.pool.del_strain(del_ind, tfr_time)
            global_wt = []
            for i,h in enumerate(self.H):
              global_wt.append([0]*len(h.indices)) # generate empty list of weights for each strain in each hospital
              for j,n in enumerate(h.indices):
                global_wt[i][j] = np.exp(n[1]-tfr_time) # fill weights list with exp(t_obs - current time) for each strain in list
            wt_unlist = [i for s in global_wt for i in s] # unlist list of list
            s = sum(wt_unlist) # reduce level of list (unlist once) to sum over
            k = len(wt_unlist) # total number of strains across all hospitals (not necessarily unique)
            inf_times = [[tfr_time + t for t in h.next_inf_time(s,k,global_wt[i])] for i,h in enumerate(self.H)]
            tfr_times = [[tfr_time + t for t in h.next_tfr_time()] for h in self.H]
            inf_time = min(min(i) for i in inf_times)
            tfr_time = min(min(i) for i in tfr_times)
##            tfr_times = [[tfr_time + i for i in t.next_tfr_time()] for t in self.H]
##            inf_times = [tfr_time + t.next_inf_time() for t in self.H]
##            inf_time = min(inf_times)
##            tfr_time = min(min(i) for i in tfr_times)
        else:
            print("Clashing inf and tfr times. Rerun")
            break
        min_time = min(inf_time,tfr_time)
    fin_strain = []
    for i,strain in enumerate(self.pool.strains):
      fin_strain.append("{i}, {strain}\n".format(i="%03d" % i,strain=strain))
    fin_strain = "".join(fin_strain)
    hosp_strains = []
    indices = []
    for i,hosp in enumerate(self.H):
      indices = []
      for j in hosp.indices:
        indices.append([j[0],"%.2f" % j[1]])
      hosp_strains.append("Hospital {i}: {h}\n".format(i="%02d" % i, h=indices))
    hosp_strains = "".join(hosp_strains)
    return print(fin_strain,hosp_strains,"Transfers:",*transfers, sep="\n")

x = [Strain(["A"],0,0,0,[0]*20),Strain(["B"],0,0,0,[0]*20),Strain(["C"],0,0,0,[0]*20),Strain(["D"],0,0,0,[0]*20)] #pool of all strains at initialisation
n = 4
strains = [[0,1],[0],[2],[3]] #list of length n=#hosp where [strain index] for each hospital
# p = [[0,0,1,0],
#      [0,0,0,1],
#      [0,0,1,0],
#      [0,0,0,1]]
tfr_rate = [[0,1,0.5],
            [0,0,1],
            [1,0,0]]
tfr_rate_r = [[]]
sys = System(n=3,inf_rate=0.2,mut_rate=1,tfr_rate=tfr_rate,tfin=5,hinit=1,len_seq=1000)
# print(sys)
sys.transfers()
# print(sys.pool.tree.get_ascii(show_internal=True,attributes=['name','hospital']))
# print(sys.pool.tree.write(format=3))
# print(sys.pool.d_var)
ts = TreeStyle()
#ts.optimal_scale_level = "full"
ts.show_leaf_name = True
ts.show_branch_length = False

colour=["red","blue","green","black"]
nstyle = [NodeStyle(),NodeStyle(),NodeStyle(),NodeStyle()]
for i,val in enumerate(nstyle):
  nstyle[i]["fgcolor"] = colour[i]
  nstyle[i]["size"] = 0

for i in range(n):
  for node in sys.pool.tree.iter_search_nodes(hospital = i):
    node.set_style(nstyle[i])
for i in range(3):
    for node in sys.H[i].tree.traverse():
        node.set_style(nstyle[0])
sys.pool.tree.show(tree_style=ts)
sys.H[1].tree.show(tree_style=ts)
#sys.pool.tree.render("test3.png", w=180, units="mm")
#print(sys.pool.d_var["S7"].t_birth)
