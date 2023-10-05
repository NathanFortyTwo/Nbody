import matplotlib.pyplot as plt


log_files = ["brute_force","barnes_hut"]

L =[]
for file in log_files:
    Lx = []
    Ly = []
    path = f"results/{file}.txt"
    with open(path,'r') as f:
        lines = f.readlines()[:-1] # remove last line
        
    for line in lines:
        nbproc, time = line.split(':')
        Lx.append(int(nbproc))
        Ly.append(float(time))

    L.append([Lx,Ly])

for couple in L:
    Lx, Ly = couple
    plt.plot(Lx,Ly)

plt.savefig("results/graph.png")

