import matplotlib.pyplot as plt
import numpy as np   

loss_trial = np.load("1_70_65536_Tanh_trainLossWWW0_20_unosi_unono.npy")
loss_test = np.load("1_70_65536_Tanh_testLossWWW0_20_unosi_unono.npy")

x = np.arange(0, len(loss_test))

plt.figure()

plt.semilogy(x,loss_trial,'b', marker = 'o')
plt.plot(x,loss_test,'r')
plt.ylabel("NN loss for W")
plt.xlabel("epochs")
plt.legend(['train loss','test loss'])
plt.savefig("plot_loss_WWW_70n_0_20.pdf", format='pdf',bbox_inches='tight',pad_inches = 0)

""" print(len(loss_test))
print(len(loss_trial))
print(loss_test)
print(loss_trial) """