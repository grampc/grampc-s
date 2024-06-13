% This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
%
% GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
%
% Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
% All rights reserved.
%
% GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

function writeGaussianProcessData(fileName, inputData, outputData, outputNoiseVariance)
    dimInput = size(inputData);
    dimOutput = size(outputData);

    if(dimInput(2) ~= dimOutput(2) || dimOutput(1) ~= 1)
        error('Dimensions of Matrices do not match.')
    end

    writematrix([dimInput outputNoiseVariance], fileName, 'WriteMode', 'overwrite', 'Delimiter', ' ')
    writematrix(inputData, fileName, 'WriteMode', 'append', 'Delimiter', ' ')
    writematrix(outputData, fileName, 'WriteMode', 'append', 'Delimiter', ' ')
end