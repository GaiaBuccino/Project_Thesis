import pyvista
from tqdm import tqdm
import numpy as np
from ezyrb import POD, RBF, Database, GPR, ANN, AE
from ezyrb import ReducedOrderModel as ROM
import matplotlib.pyplot as plt
import torch.nn as nn
from time import time
import sys

# lancia con python3 dataFOM.py nome_file_csv w/psi a seconda della variabile richiesta 

file_name = sys.argv[1]
which_filed = sys.argv[2]

PERFORM_ROM = True
Ntime = 301
Ntrain = 251
Nc = 65536
snapshots_w = np.zeros((Ntime, Nc))
snapshots_psi= np.zeros((Ntime, Nc))
params = np.zeros((Ntime,1))


# Import dataset
for n in tqdm(range(Ntime), desc="Importing data"):
    filename = "/scratch/gbuccino/git/ITHACA-FV/tutorials/CFD/vortexMergerWPsi_ITHACA/ITHACAoutput/Offline/off.foam"
    reader = pyvista.OpenFOAMReader(filename)
    reader.set_active_time_value(n)
    mesh = reader.read()
    w = np.array(mesh["internalMesh"]["W"])
    psi_scalar = np.array(mesh["internalMesh"]["Psi_z"])
    snapshots_w[n, :] = w
    snapshots_psi[n, :] = psi_scalar
    params[n] = 0.08*n

np.save('snapshots_w.npy', snapshots_w)
np.save('params.npy', params)

param_extr = params[Ntrain+1:, :]


if which_filed == 'w':
    inputs = snapshots_w[0:Ntrain,:]
    snapshots = snapshots_w
    inputs_extr = snapshots_w[Ntrain+1:,:]
    param = params[0:Ntrain, :]
    N_snap = inputs.shape[0]
    dimensions = [14]
elif which_filed == 'psi':
    inputs = snapshots_psi[2:Ntrain, :]
    snapshots = snapshots_psi
    inputs_extr = snapshots_psi[Ntrain+1:,:]
    param = params[2:Ntrain, :]    
    N_snap = inputs.shape[0]
    dimensions = [6]
else:
    raise ValueError


parameters = param

print(f"N_snap: {N_snap}")

# splitting indeces for train and test
mask_train = [True, False]*(N_snap//2)
mask_test = [False, True]*(N_snap//2)
mask_train += [True]
mask_test += [False]

mask_train = np.array(mask_train)
mask_test = np.array(mask_test)
print(f"snap_training: {sum(mask_train)}")
print(f"snap_testing: {sum(mask_test)}")

params_training = parameters[mask_train, :]
snapshots_training = inputs[mask_train, :]
params_testing = parameters[mask_test, :]
snapshots_testing = inputs[mask_test, :]

if PERFORM_ROM:
    import pandas as pd
    data_res = []

    # input
    input_dim = inputs.shape[1]
    db_training = Database(params_training, snapshots_training)
    db_testing = Database(params_testing, snapshots_testing)

    method = {"RBF": RBF(),
              "GPR": GPR(),
              "ANN": ANN([16, 64, 64], nn.Tanh(), [30000, 1e-12])}


    for dim in dimensions:
    
        for key, interp in method.items():

            pod = POD('svd', rank=dim)
            rbf = interp
            rom = ROM(db_training, pod, rbf)
            start = time()
            rom.fit()
            end = time()

            train_error = rom.test_error(db_training)
            test_error = rom.test_error(db_testing)
            data_res.append(
                [f"POD + {key}", dim, 100*train_error, 100*test_error, end - start])

            print(f"POD + {key} == dimension: {dim}")
            print(f"    ROM train error : {train_error:.5%}")
            print(f"    ROM test error : {test_error:.5%} ")
            print()

            if key == 'ANN':
                loss = interp.loss_trend
                epochs = [i for i in range(len(loss))]

                plt.semilogy(epochs[::500],loss[::500],'--o')
                plt.ylabel(f"NN loss for {which_filed.title()}")
                plt.xlabel("epochs")
                plt.legend(['train loss','test loss'])
                plt.savefig(f"{file_name}.pdf", format='pdf',bbox_inches='tight',pad_inches = 0)
                plt.close()
    
            data_res.append(["-", "-", "-", "-", "-"])

            #usare come errore np.linalg.norm()/stessa_cosa(soluz reale)
            #computing extrapolation error
            xx = np.zeros(len(param_extr))
            error_extr = np.zeros(len(param_extr))

            for t in range(len(param_extr)):
                start = time()
                pred_solution = rom.predict(param_extr[t])
                end = time()
                error_extr[t] = 100 * np.linalg.norm(pred_solution - snapshots[Ntrain+t])/np.linalg.norm(snapshots[Ntrain+t])
                xx[t] = params[Ntrain+t]
            print(f'ti sto per dire il tempo di previsione extrapolation per {key} for the field {which_filed}:', end - start)
            plt.figure()
            plt.plot(xx, error_extr, 'r')
            plt.title(f"{key} extrapolation error")
            plt.ylabel(f"{which_filed} error %")
            plt.xlabel("time")

            plt.savefig(f"{key}_Extrapolation_error_{which_filed}.pdf", format='pdf',bbox_inches='tight',pad_inches = 0)
            plt.close()

                    
            tim = [10, 60, 120]
            for idx in tim:
                p = params_testing[idx]
                rom_sol = interp.predict(params_testing)
                rom_sol = (pod.inverse_transform(rom_sol.T).T)[idx].reshape((256, 256))
                real_sol = snapshots_testing[idx].reshape((256, 256))

                plt.imshow(rom_sol,cmap=plt.cm.jet,origin='lower')
                plt.title(f'time {p[0]} (s)')
                plt.colorbar()
                plt.savefig(f'rom_{which_filed}_time_{p}_method_{key}.pdf')
                plt.close()

                plt.imshow(real_sol,cmap=plt.cm.jet,origin='lower')
                plt.title(f'time {p[0]} (s)')
                plt.colorbar()
                plt.savefig(f'fom_{which_filed}_time_{p}.pdf')
                plt.close()

                plt.imshow(rom_sol - real_sol,cmap=plt.cm.jet,origin='lower')
                plt.title(f'time {p[0]} (s)')
                plt.colorbar()
                plt.savefig(f'error_{which_filed}_time_{p}_{key}.pdf')
                plt.close()
                
                

            tim_extr = [5, 25, 45]
            error_extr = np.zeros(len(param_extr))
            xxxx = np.zeros(len(param_extr))

            for idxx in tim_extr:
                p = param_extr[idxx]
                rom_sol = interp.predict(param_extr)
                rom_sol = (pod.inverse_transform(rom_sol.T).T)[idxx].reshape((256, 256))
                real_sol = inputs_extr[idxx].reshape((256, 256))

                plt.imshow(rom_sol,cmap=plt.cm.jet,origin='lower')
                plt.title(f'time {p[0]} (s)')
                plt.colorbar()
                plt.savefig(f'rom_{which_filed}_time_{p}_method_{key}.pdf')
                plt.close()

                plt.imshow(real_sol,cmap=plt.cm.jet,origin='lower')
                plt.title(f'time {p[0]} (s)')
                plt.colorbar()
                plt.savefig(f'fom_{which_filed}_time_{p}.pdf')
                plt.close()

                plt.imshow(rom_sol - real_sol,cmap=plt.cm.jet,origin='lower')
                plt.title(f'time {p[0]} (s)')
                plt.colorbar()
                plt.savefig(f'error_{which_filed}_time_{p}_{key}.pdf')
                plt.close()

                pred_solution = rom.predict(p)
                error_extr[idxx] = 100 * np.linalg.norm(pred_solution - snapshots[Ntrain+idxx])/np.linalg.norm(snapshots[Ntrain+idxx])
                xxxx[idxx] = idxx

                plt.figure()


                plt.plot(xxxx, error_extr, 'o')
                plt.title(f"{key} error extrapolated")
                plt.ylabel(f"{which_filed} error %")
                plt.xlabel("extrapolated instants")

                plt.savefig(f"ROM_extrapolation_error_{key}_{which_filed}_24.pdf", format='pdf',bbox_inches='tight',pad_inches = 0)
                plt.close()

            #computing error of ROM approx over the time interval
            error_ROM = np.zeros(len(params_testing))
            xxx = np.zeros(len(params_testing))
            for t in range(len(params_testing)):

                pred_sol = rom.predict(params_testing[t])
                error_ROM[t] = 100 * np.linalg.norm(pred_sol - snapshots_testing[t])/np.linalg.norm(snapshots_testing[t])
                xxx[t] = params_testing[t]


                plt.figure()


                plt.plot(xxx, error_ROM, 'r')
                plt.title(f"{key} error")
                plt.ylabel(f"{which_filed} error %")
                plt.xlabel("time")

                plt.savefig(f"ROM_error_{key}_{which_filed}.pdf", format='pdf',bbox_inches='tight',pad_inches = 0)
                plt.close()

                
        print("===============")

            

    df_res = pd.DataFrame(data_res, columns=['method', 'dim', 'train %', 'test %', 'train time [s]'])
    df_res.to_csv(f'{file_name}.csv')












