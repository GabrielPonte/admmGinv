using LinearAlgebra, MAT, DelimitedFiles, SparseArrays

function getMatlabInstance(instanceName::String,matrixName::String)
    myInst = string("..\\Instances\\",instanceName,".mat");
    file = matopen(string(myInst));
    A = read(file, string(matrixName));
    close(file);
    A = Matrix(A);
    return A
end

# function getGinvInstance(instanceName::String)
#     A = getMatlabInstance(instanceName,"A");
#     R = getMatlabInstance(instanceName,"R");
#     C = getMatlabInstance(instanceName,"C");
#     return A,R,C
# end