import matplotlib.pyplot as plt
import os
dir = "runtimes"
log_files = [file for file in os.listdir(dir) if file.endswith(".txt")]

L = []
nbparticules = [500, 1000, 1500, 2000, 2500]

for file in log_files:
    Lx = []
    Ly = []
    path = f"{dir}/{file}"
    with open(path,'r') as f:
        lines = f.readlines()#[:-1] # remove last line

    for index,line in enumerate(lines):
        time = line.strip()
        Lx.append(nbparticules[index])
        Ly.append(float(time))
    
    name = file.split(".")[0].replace("_nbody_brute_force","")
    plt.plot(Lx,Ly,label = name)
    plt.legend()
    plt.xlabel("Nombre de particules")
    plt.ylabel("Temps d'execution (s)")
    
graph_name = "graphOMP"
plt.savefig(dir+f"/{graph_name}.png")