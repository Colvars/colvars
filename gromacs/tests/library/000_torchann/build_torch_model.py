import torch

class MyModel(torch.nn.Module):
    def __init__(self):
        super().__init__()
    def forward(self, x):
        return x

print ("PyTorch Version:", torch.__version__)

model = MyModel()
scripted_cv_filename = f'./identity.pt'
torch.jit.script(model).save(scripted_cv_filename)

