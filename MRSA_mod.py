import numpy as np

class Pool:
  def __init__(self, strains):
    self.strains = strains

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'Pool(strains = '+str(self.strains)+')'
    
  def add_new_strain(self, mut_rate, index, t):
    num_muts = np.random.poisson(mut_rate*(t - self.strains[index].birth_time), 1)[0]
    parent_strain = self.strains[index]
    parent_strain.num_child += 1
    patient_strain = Strain(parent_strain.parent + [parent_strain.num_child], 0, num_muts, t)
    self.strains.append(patient_strain)
    return len(self.strains)-1
 
 
class Strain:
  def __init__(self, parent, num_child, mut, birth_time):
    self.parent = parent
    self.num_child = num_child
    self.mut = mut
    self.birth_time = birth_time

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return  'Strain(parent = '+str(self.parent)+', num_child = '+str(self.num_child)+', mut = '+str(self.mut)+', birth_time = '+str(self.birth_time)+')'
 
 
class Hospital:
  def __init__(self, indices, tfin, mut_rate, inf_rate, obs_prob = 0.2):
    self.indices = indices
    self.tfin = tfin
    self.mut_rate = mut_rate
    self.inf_rate = inf_rate
    self.obs_prob = obs_prob

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return 'Hospital(indicies = '+str(self.indices)+', tfin = '+str(self.tfin)+', mut_rate = '+str(self.mut_rate)+', inf_rate = '+str(self.inf_rate)+', obs_prob = '+str(self.obs_prob)+')'
  
  def next_time(self):
    return -np.log(np.random.uniform())/self.inf_rate
  
  def add_inf(self, t, pool):
    index = np.random.randint(len(self.indices))
    index = self.indices[index]
    self.indices.append(pool.add_new_strain(self.mut_rate, index, t))
    
  def add_index(self, index):
    self.indices.append(index)


class System:
  def __init__(self, n, inf_rate, mut_rate, tfin, pool, strains, p):
    self.n = n #number of hospitals
    self.tfin = tfin
    self.inf_rate = inf_rate
    self.mut_rate = mut_rate
    self.p = p #threshold probability for transfer
    self.pool = Pool(pool)
    self.H = []
    for i in range(n):
      self.H.append(Hospital(strains[i], self.tfin, self.mut_rate, self.inf_rate))

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
      hap_ind = next_hosp.indices[np.random.randint(len(next_hosp.indices))] #randomly select an entry from the indicies attribute of Hospital corresponding to the index of a strain in pool
      temp_ind = self.pool.add_new_strain(self.mut_rate, hap_ind, inf_time) #return index of new strain that is the child from one of the strains in the current incidence of H
      self.H[hosp_ind].add_index(temp_ind) #add that strain to the list of strains in that H
      trans_prob = np.random.uniform(0,1)
      if trans_prob < self.p:
        transfers.append([inf_time, hosp_ind])
        new_ind = np.random.randint(len(self.H))
        while new_ind == hosp_ind:      #ensure the strain is not transfered to the same hospital
          new_ind = np.random.randint(len(self.H))
        self.H[new_ind].add_index(temp_ind)
      times = [inf_time + t.next_time() for t in self.H]
      inf_time = min(times)
    
    return self.pool.strains, [hosp.indices for hosp in self.H], transfers
  
  
