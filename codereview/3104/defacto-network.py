import types
import pickle
import math
import time
import random
import os
import networkx as nx
import numpy as np
from scipy import sparse,io
import scipy.io as sio
import matplotlib.pyplot as plt
import randht

ntask = 100;
nodes = 5000;
timesteps = 1000;
wins = np.zeros((nodes,timesteps+10));

scheduled = dict();
timestamp = time.strftime('%X-%y-%m-%d');
ofile = open('scheduled' + timestamp + '.txt','wb');
scheds = np.zeros((ntask,ntask));

for i in range(ntask):
  pl = randht.randht(1,'powerlaw',1.5);
  effective = int(pl[0]*100);
  while effective >= nodes:
		pl = randht.randht(1,'powerlaw',1.5);
  		effective = int(pl[0]*100);

  G = nx.generators.random_graphs.watts_strogatz_graph(effective, int(math.floor(np.log(effective)) ** 1.2), 0.1);
  correspondence = range(nodes);
  random.shuffle(correspondence);
  GG=nx.Graph()
  for ed in G.edges():
    if not GG.has_edge(ed[1],ed[0]):
    	GG.add_edge(correspondence[ed[0]],correspondence[ed[1]]);
  
  winsize = 1 + int(math.floor(GG.number_of_nodes()*0.1));
  ofile.write(str(winsize) + ',')
  runs = int(timesteps/winsize);
  for s in range(runs-1):
    chance = np.random.random_sample();
    if chance < 0.1:#5*GG.number_of_nodes()/nodes: 
      ofile.write(str(s+1) + ',');
      if i in scheduled:
	temp = scheduled[i];
	temp.append(s);
	scheduled[i] = temp;
      else:
	temp = list();
	temp.append(s);
	scheduled[i] = temp;
      centrality = nx.betweenness_centrality(GG);
      for node in GG:
	#print centrality[node];	
      	wins[node,(s)*winsize + int(centrality[node]*200)-1] = i+1;
      	for ad in nx.edges_iter(GG,node):
	  offset = np.random.randint(1,3,[1,1]);
	  wins[ad,(s)*winsize + int(centrality[node]*200) + offset] = i+1;
  ofile.seek(-1, os.SEEK_CUR);
  ofile.write(',-');  
  ofile.write('\n');
  print i;



sampled = np.random.randint(0,nodes-1,[1,100]);
sio.savemat('raster-' + timestamp + '.mat',{'wins':wins[sampled,:]})


