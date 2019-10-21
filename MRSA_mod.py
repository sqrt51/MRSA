import numpy as np

class Pool:
  def __init__(self, strains):
    self.strains = strains

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'Pool(strains = '+str(self.strains)+')'

  def add_new_strain(self, mut_rate, index, t):
    num_muts = np.random.poisson(mut_rate*(t - self.strains[index].t_birth))
    parent_strain = self.strains[index]
    parent_strain.num_child += 1
    patient_strain = Strain(parent_strain.parent + [parent_strain.num_child], 0, num_muts, t)
    self.strains.append(patient_strain)
    return len(self.strains)-1


class Strain:
  def __init__(self, parent, num_child, mut, t_birth):
    self.parent = parent
    self.num_child = num_child
    self.mut = mut
    self.t_birth = t_birth

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return  'Strain(t_birth = '+str("%.2f" % self.t_birth) \
      +', mut = '+str(self.mut) \
      +', num_child = '+str(self.num_child) \
      +', parent = '+str(self.parent) \
      +')'


class Hospital:
  def __init__(self, indices, tfin, mut_rate, inf_rate, pop, obs_prob = 0.2):
    self.indices = indices
    self.tfin = tfin
    self.mut_rate = mut_rate
    self.inf_rate = inf_rate
    self.pop = pop
    self.obs_prob = obs_prob

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'Hospital(indicies = '+str(self.indices)+', tfin = '+str(self.tfin)+', mut_rate = '+str(self.mut_rate)+', inf_rate = '+str(self.inf_rate)+', obs_prob = '+str(self.obs_prob)+')'

  def next_time(self):
    return -np.log(np.random.uniform())/self.inf_rate

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

  def add_index(self, index):
    del_ind = np.random.choice(len(self.indices[0]),p=[i/sum(self.indices[1]) for i in self.indices[1]])
    if self.indices[1][del_ind] > 1:
        self.indices[1][del_ind] -= 1
    elif self.indices[1][del_ind] == 1:
        del self.indices[0][del_ind]
        del self.indices[1][del_ind]
    self.indices[0].append(index)
    self.indices[1].append(1)


class System:
  def __init__(self, n, inf_rate, mut_rate, tfin, pool, strains, p, pop):
    self.n = n #number of hospitals
    self.tfin = tfin
    self.inf_rate = inf_rate
    self.mut_rate = mut_rate
    self.p = p #threshold probability for transfer
    self.pool = Pool(pool)
    self.pop = pop #population size of total number of cases at each site (same at each site)
    self.H = []
    for i in range(n):
      self.H.append(Hospital(strains[i], self.tfin, self.mut_rate, self.inf_rate, self.pop))

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'System(n = '+str(self.n)+', inf_rate = '+str(self.inf_rate)+', mut_rate = '+str(self.mut_rate)+', tfin = '+str(self.tfin)+', pool = '+str(self.pool)+', p = '+str(self.p)+', H = '+str(self.H)+')'
    #', strains = '+str(self.strains)+

  def transfers(self):
    transfers = []
    times = [t.next_time() for t in self.H]
    inf_time = min(times)
    while inf_time < self.tfin:
      hosp_ind = times.index(inf_time)  #select the index of the current min of times
      next_hosp = self.H[hosp_ind]      #select the hospital corresponding to that index
      hap_ind = next_hosp.indices[0][np.random.choice(len(next_hosp.indices[0]),p=[i/sum(next_hosp.indices[1]) for i in next_hosp.indices[1]])] #randomly with weighting select an entry from the indicies attribute of Hospital corresponding to the index of a strain in pool
      temp_ind = self.pool.add_new_strain(self.mut_rate, hap_ind, inf_time) #return index of new strain that is the child from one of the strains in the current incidence of H
      self.H[hosp_ind].add_index(temp_ind) #add that strain to the list of strains in that H
      new_ind = np.random.choice(range(len(self.H)),p=self.p[hosp_ind]) #randomly pick a new hospital index with weighted probability where a transfer to the same hospital is used as a convention for no transfer taking place
      if new_ind != hosp_ind:
        transfers.append([inf_time, hosp_ind, new_ind])
        self.H[new_ind].add_index(temp_ind)
      times = [inf_time + t.next_time() for t in self.H]
      inf_time = min(times)
    fin_strain = []
    for i,strain in enumerate(self.pool.strains):
      fin_strain.append("{i}, {strain}\n".format(i="%03d" % i,strain=strain))
    fin_strain = "".join(fin_strain)
    hosp_strains = []
    for i,hosp in enumerate(self.H):
      hosp_strains.append("Hospital {i}: {hosp}\n".format(i="%02d" % i, hosp=hosp.indices))
    hosp_strains = "".join(hosp_strains)
    return print(fin_strain,hosp_strains,"Transfers:",*transfers, sep="\n")

x = [Strain([1],0,0,0),Strain([2],0,0,0),Strain([3],0,0,0),Strain([4],0,0,0)] #pool of all strains at initialisation
n = 4
strains = [[[0,1],[9,1]],[[0],[10]],[[2],[10]],[[3],[10]]] #list of length n=#hosp where [strain index],[number of each strain]
p = [[0,0,1,0],
     [0,0,0,1],
     [0,0,1,0],
     [0,0,0,1]]
sys = System(n,1,1,5,x,strains,p,pop=1)
sys.transfers()
