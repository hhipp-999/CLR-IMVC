addpath("dataset")
load("EYaleB10.mat")
for iv = 1 : length(X)
    Z{iv} = X{iv}';
end
Y=gt;