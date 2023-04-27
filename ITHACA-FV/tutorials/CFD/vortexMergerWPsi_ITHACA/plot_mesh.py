import pyvista
import numpy as np   
import ezyrb
from tqdm import tqdm

import matplotlib.tri as mtri
import matplotlib.pyplot as plt
from ezyrb import POD, RBF, Database
from ezyrb import ReducedOrderModel as ROM
#%matplotlib inline

Nt = 301
Nc = 65536

#for n in tqdm(range(Nt), desc="Importing data"):

filename = "/scratch/gbuccino/git/ITHACA-FV/tutorials/CFD/vortexMergerWPsi_ITHACA/ITHACAoutput/Offline/off.foam"
reader = pyvista.OpenFOAMReader(filename)
reader.set_active_time_value(300)
mesh = reader.read()
p = pyvista.Plotter(notebook=False, off_screen=True)
p.add_mesh(mesh, scalars = "W",  cmap="jet", clim = [0,9.4e-1],
    show_scalar_bar = False)
p.camera_position = 'xy'
p.screenshot("WFOM_24_mesh.png")


ROM = np.load("wAANN_w_24.npy")
tmp1 = np.reshape(ROM, (256,256)) 

plt.imshow(tmp1[:,:], cmap = plt.cm.jet , origin = 'lower')
plt.savefig('ciao.png')
#mesh["internalMesh"]["w_rom"] = mesh["internalMesh"]["W"].copy()   #qui creiamo una nuova variabile per copia e poi la sovrascriviamo con i valori nelle celle
#WW = mesh["internalMesh"]["W"] 
#mesh["internalMesh"]["W"] =  ROM
#mesh["internalMesh"].point_data["w_rom"] = mesh["internalMesh"].point_data["W"].copy() #qui creiamo una nuova variabile per copia e poi la sovrascriviamo con i valori nei punti 
#mesh["internalMesh"].point_data["w_rom"] = ROM
#print(mesh["internalMesh"]["w_rom"])
""" plt.figure()
p= pyvista.Plotter(notebook=False, off_screen=True)
#mesh.set_active_scalars("w_rom")
p.add_mesh(mesh, scalars = "W",  cmap="jet", clim = [0,9.4e-1],
    show_scalar_bar = True)
p.camera_position = 'xy'
p.screenshot("WROM_24_mesh.png")

plt.figure()
p = pyvista.Plotter(notebook=False, off_screen=True)
p.add_mesh(mesh, scalars = "Psi_z",  cmap="jet", clim = [0,4.8e-1],
    show_scalar_bar = True)
p.camera_position = 'xy'
p.screenshot("PsiFOM_24_mesh.png")
ROM = np.load("psiAANN_psi_24.npy")
#mesh["internalMesh"]["w_rom"] = mesh["internalMesh"]["W"].copy()   #qui creiamo una nuova variabile per copia e poi la sovrascriviamo con i valori nelle celle
mesh["internalMesh"]["Psi_z"] = ROM
#mesh["internalMesh"].point_data["w_rom"] = mesh["internalMesh"].point_data["W"].copy() #qui creiamo una nuova variabile per copia e poi la sovrascriviamo con i valori nei punti 
#mesh["internalMesh"].point_data["w_rom"] = 0
#print(mesh["internalMesh"]["w_rom"])
plt.figure()
p = pyvista.Plotter(notebook=False, off_screen=True)
#mesh.set_active_scalars("w_rom")
p.add_mesh(mesh, scalars = "Psi_z",  cmap="jet", clim = [0,4.8e-1],
    show_scalar_bar = False)
p.camera_position = 'xy'
p.screenshot("PsiROM_24_mesh.png")




#np.save("snap_test1_2.npy", snap_test1_2)
#np.save("param_test1_2.npy", param_test1_2)

 """