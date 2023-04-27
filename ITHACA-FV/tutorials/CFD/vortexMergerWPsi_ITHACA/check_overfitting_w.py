import numpy as np
import torch
from torch import nn
from sklearn.model_selection import train_test_split
import os
from torch.utils.data import Dataset
torch.manual_seed(0)

torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

snap_trial_w = np.load("snap_trial_w0_20.npy")
param_trial = np.load("param_trial0_20.npy")
snap_test_w = np.array(np.load("snap_test_w0_20.npy"))
param_test = np.array(np.load("param_test0_20.npy"))

class NN(nn.Module):
    def __init__(self, in_features, nb_classes,
        hidden_sizes, act=nn.ReLU):
        super(NN, self).__init__()
        self.actName = act
        self.act = act()
        self.layers = []
        self.hidden_sizes = hidden_sizes
        self.in_features = in_features
        self.nb_classes = nb_classes
        self.layers.append(nn.Linear(in_features, hidden_sizes[0]))
        self.layers.append(act())
        for i, k in enumerate(hidden_sizes[:-1]):
            self.layers.append(nn.Linear(hidden_sizes[i], hidden_sizes[i+1]))
            self.layers.append(act())
        self.layers.append(nn.Linear(hidden_sizes[-1], nb_classes))
        self.net = nn.Sequential(*self.layers)

    def forward(self, x):
        x = self.net.forward(x)
        return x

    def hash(self):
        name = str(self.in_features)
        for i in self.hidden_sizes:
            name+= "_"+str(i)
        name+= "_"+str(self.nb_classes)
        name+= ("_"+str(self.act))[:-2]
        return name

class NNDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.from_numpy(X).type(torch.float32)
        self.y = torch.from_numpy(y).type(torch.float32)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return [self.X[idx], self.y[idx]]

###############################################################################

def train(**kwargs):
    '''
    Training of the neural network
    '''
    input_data = kwargs.get('input_data')
    output_data = kwargs.get('output_data')
    device = kwargs.get('device', torch.device('cpu'))
    Net = kwargs.get('nn')
    epochs = kwargs.get('epochs', 50000)
    learning_rate = kwargs.get('learning_rate', 1e-5)
    batch_size = kwargs.get('batch_size', 4)

    lossplottrain = []
    lossplottest = []

    X_train, X_test, Y_train, Y_test = param_trial, param_test, snap_trial_w, snap_test_w

    trainset = NNDataset(X_train, Y_train)
    testset = NNDataset(X_test, Y_test)
    trainloader = torch.utils.data.DataLoader(trainset,
            batch_size=batch_size, shuffle=True)
    testloader = torch.utils.data.DataLoader(testset,
            batch_size=X_test.shape[0], shuffle=True)


    #modelfile = Net.hash()+".pt"
    optimizer = torch.optim.Adam(Net.parameters(), lr=learning_rate)
    #if os.path.isfile(modelfile):
       # print('Reading existing training for inp->out mapping {}'.format(
         #   modelfile))
       # Net = torch.jit.load(modelfile)
    #else:
    print('Training the NN with data of vectors inp and out')
    for t in range(epochs):
        batch_losses = []
        for inputs, labels in trainloader:
            inputs, labels =  inputs.to(device), labels.to(device)
            optimizer.zero_grad()
            # forward + backward + optimize
            outputs = Net(inputs)
            mse_loss = torch.nn.MSELoss()
            loss = mse_loss(outputs, labels)
            loss.backward()
            optimizer.step()
            batch_losses.append(loss.item())
        loss = np.mean(batch_losses)

        # evaluate accuracy on test set
        batch_test_losses = []
        Net.eval()
        for inputs_test, labels_test in testloader:
            inputs_test, labels_test =  inputs_test.to(device), labels_test.to(device)
            outputs_test = Net(inputs_test)
            mse_loss = torch.nn.MSELoss()
            test_loss = mse_loss(outputs_test, labels_test)
            batch_test_losses.append(test_loss.item())
        test_loss = np.mean(batch_test_losses)
        lossplottrain.append(loss)
        lossplottest.append(test_loss)
        if t % 100 == 0:
            print(t, "Loss on train", loss)
            print(t, "Loss on test", test_loss)
        Net.train()

    # Save the model and scaling
    m = torch.jit.script(Net)
    np.save(Net.hash()+"_trainLossWWW0_20_unosi_unono.npy", lossplottrain)
    np.save(Net.hash()+"_testLossWWW0_20_unosi_unono.npy",lossplottest)
    m.save(Net.hash()+".pt")


if __name__ == '__main__':
    
    input_data = np.load("param_trial0_20.npy")
    target_data = np.load("snap_trial_w0_20.npy")

    # Create Net
    
    layers = [80]
    Net = NN(input_data.shape[1], target_data.shape[1], layers, torch.nn.Tanh)
    # Train the NN
    train(input_data=input_data,
          output_data=target_data,
          nn=Net, learning_rate=1e-5, epochs=50000)