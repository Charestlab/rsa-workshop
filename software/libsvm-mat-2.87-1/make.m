% This make.m is used under Windows

mex -largeArrayDims -O -c svm.cpp
mex -largeArrayDims -O -c svm_model_matlab.c
mex -largeArrayDims -O svmtrain.c svm.obj svm_model_matlab.obj
mex -largeArrayDims -O svmpredict.c svm.obj svm_model_matlab.obj
mex -largeArrayDims -O read_sparse.c
