from numpy import random
from inspect import Parameter
from matplotlib import pyplot as plt

class Walker:
  def __init__(self,initpos1,initpos2,deg,L):
    self.position = [(initpos1,initpos2)] #contains position history of walker. could not save history to save memory.
    self.degreeofinfection = [float(deg)]
    #print(float(deg))
    self.SIR = [getstate(float(deg))]
    self.L = L
  #contains the walker's infection status and the walker's position
  def updatepos(self,xmove,ymove):
    if abs(self.position[-1][0] + xmove) <= L and abs(self.position[-1][1] + ymove) <= L:
      self.position.append((self.position[-1][0] + xmove, self.position[-1][1] + ymove))
    else:
      #this only works properly when h << L
      newposx = self.position[-1][0] + xmove
      newposy = self.position[-1][1] + ymove
      if abs(self.position[-1][0] + xmove) > L:
        newposx = ((L - abs(L - abs(self.position[-1][0] + xmove))))*abs(newposx)/newposx
      if  abs(self.position[-1][1] + ymove) > L:
        newposy = ((L - abs(L - abs(self.position[-1][1] + ymove))))*abs(newposy)/newposy
      self.position.append((newposx, newposy))
  def changedeg(self,change):
    #if the walker isn't getting infected at this time step, we change the degree of infection by 1/tau1
    self.degreeofinfection.append(self.degreeofinfection[-1] - change)
    self.SIR.append(getstate(self.degreeofinfection[-1]))
  def becomeinfected(self):
    self.degreeofinfection.append(1.0)
    self.SIR.append(getstate(self.degreeofinfection[-1]))
    #if the walker becomes infected, this sets the degree of infection accordingly

def getstate(deg):
  if deg <= 0:
    return 0 #susceptible
  elif (1.0 - tau1/tau2) <= deg and (1.0 >= deg):
    return 1 #infectious
  else:
    return 2 #recovered

class RandomWalk:
  def __init__(self,L,h,nwalkers,p,tau1,tau2,starting):
    self.L = L
    self.h = h
    self.nwalkers = nwalkers
    self.walkers = [] #list of all Walker objects
    self.p = p
    self.degchange = 1.0/tau1
    self.starting = starting
    self.initializewalkers()
    self.ninfected = [self.starting]
    self.newinf = [0,0]
  def paramstring(self,string):
    self.paramstring = string
  def initializewalkers(self):
  #initializes the first element of poshistory with the positions of all of the walkers
  #use this to randomly initialize the position of a walker on the graph
    for i in range(self.nwalkers):
      self.walkers.append(Walker(random.randint(-L,L+1),random.randint(-L,L+1),i//(self.nwalkers - self.starting),L))
#given grid size (xmax,ymax); range of motion h; and position of a walker (origx,  origy), randomly return updated positions for the walkers.
#Also updates degrees of infection.
  def update(self):
    for w in self.walkers:
      w.updatepos(random.randint(-self.h,self.h+1),random.randint(-self.h,self.h+1))
    walkerdict = {}
    for w in self.walkers:
      if w.position[-1] in walkerdict:
        walkerdict[w.position[-1]].append(w)
      else:
        walkerdict[w.position[-1]] = [w]
    self.ninfected.append(0)
    self.newinf.append(0)
    for pos in walkerdict:
      infected = []
      susceptible = []
      recovered = []
      for w in walkerdict[pos]:
        if getstate(w.degreeofinfection[-1]) == 0:
          susceptible.append(w)
        elif getstate(w.degreeofinfection[-1]) == 1:
          infected.append(w)
          w.changedeg(self.degchange)
        else:
          recovered.append(w)
          w.changedeg(self.degchange)
      self.ninfected[-1] += len(infected)
      for s in susceptible:
        if random.binomial(1, 1-(1.0-self.p)**len(infected), size=None) == 1:
          s.becomeinfected()
          self.newinf[-1] += 1
        else:
          s.changedeg(self.degchange)
    #based on the current distribution of walkers, update the degree of infection of each walker
    #if walkers linger together for a time, this function will treat each time step as providing equal opportunity for infection
    #when there are multiple infectious walkers at a node, each independently has a p chance of infecting a susceptible walker
    #at the initialization step, no infection can occur.
      


#Below, we use the above classes to run simulations. Each tuple represents a different selection of parameters for the random walk.


tuples = [(50,5000,20,2,float(10),float(20),float(0.4),2)],(50,1000,50,2,float(3),float(20),float(0.4),2),(50,1000,50,2,float(10),float(20),float(0.2),2), (50,1000,50,2,float(10),float(20),float(0.6),2)]
walks = []
for t in tuples:
  TIMESTEPS, NWALKERS, L, H, tau1, tau2, P, starting = t
  walks.append(RandomWalk(L,H, NWALKERS, P,tau1,tau2,starting))
  walks[-1].paramstring("nwalkers = " + str(NWALKERS) + ", h = " + str(H) + ", tau1 = " + str(tau1) + ", tau2 = " + str(tau2) + ", P = " + str(P) + ", starting = " + str(starting)
)
  for t in range(TIMESTEPS):
    walks[-1].update()



#The function below allows us to visualize a random walk.


plt.rcParams["figure.autolayout"] = True
def visualize_walk(walk):
  fig, axs = plt.subplots(3,3)
  plt.axis(xmin = -50,xmax = 50,ymin = -50,ymax = 50)
  times = [i*5 for i in range(1,10)] #times to take snapshots at
  plt.xlim(-50, 50)
  plt.ylim(-50, 50)
  for i in range(9):
    time = times[i]
    j = i % 3
    k = i // 3
    for walker in walk.walkers:
      x = [walker.position[time][0]]
      y = [walker.position[time][1]]
      theplt = axs[j,k]
      theplt.grid()
      theplt.title.set_text("Walkers at time " + str(time))
      dimen = 40
      theplt.set_xlim([-dimen,dimen])
      theplt.set_ylim([-dimen,dimen])
      if walker.SIR[time] == 1:
        theplt.plot(x, y, marker=".", color="red",markersize = 5)
      elif walker.SIR[time] == 0:
        theplt.plot(x, y, marker=".", color="black",markersize = 5)
      else:
        theplt.plot(x, y, marker=".", color="blue",markersize = 5)

  fig.set_size_inches(10,10)
  paramstring = "nwalkers = " + str(NWALKERS) + ", h = " + str(H) + ", tau1 = " + str(tau1) + ", tau2 = " + str(tau2) + ", P = " + str(P) + ", starting = " + str(starting)

  fig.suptitle("With parameters "+ paramstring + ". Black dots represent susceptible walkers; red dots infectious; blue dots recovered.", y=1.02)
  fig.show()

visualize_walk(walks[0])

