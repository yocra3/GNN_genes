# Docker inheritance
FROM  graphcore/pytorch-geometric:3.3.0

# Install python packages
RUN pip3 install pandas

# Install pytorch version compatible with cuda (version 1.1)
RUN pip3 install --upgrade torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
RUN pip3 install torchinfo

# Install tensorboard (version 1.2)
RUN pip3 install torch tensorboard --index-url https://download.pytorch.org/whl/cu118



