function [CNN] = build_network(NumInputs, NumOutputs, net_type, neurons, filterSize, numHiddenUnits)

switch net_type
    case 0

        layers = [
            sequenceInputLayer(NumInputs,'Name','Sequence Input Layer',"Normalization","zscore")% prende i dati in input
        
            lstmLayer(numHiddenUnits) % long-short term layer (legami tra le variabili in input)
            dropoutLayer(0.2)
%             reluLayer('Name','Rectified Linear Unit Layer') % Applica la funzione ReLu
            
%             lstmLayer(numHiddenUnits,'OutputMode','last')
%             dropoutLayer(0.2)

            fullyConnectedLayer(NumOutputs,'Name','Fully Connected Layer')% Restituisce l'output (PM10 giorno successivo) 
            regressionLayer('Name','Regression Layer') % Minimizza il HMSE (half mean squared error)
            ];
        CNN = layerGraph(layers);% creo la rete neurale
    case 1
        layers = [
            sequenceInputLayer(NumInputs,"Name","sequence","Normalization","rescale-zero-one")
%             sequenceInputLayer(NumInputs,"Name","sequence")
        
            convolution1dLayer(filterSize,neurons,"Name","conv1d","Padding",'same')
            batchNormalizationLayer("Name","batchnorm")
            reluLayer("Name","Rectified Linear Unit Layer")
        
%             convolution1dLayer(filterSize,neurons,"Name","conv1d2","Padding",'same')
%             batchNormalizationLayer("Name","batchnorm2")
%             reluLayer("Name","Rectified Linear Unit Layer 2")
        
            fullyConnectedLayer(NumOutputs,"Name","fc")
%             sigmoidLayer('Name',"sigmoidoutput")
            regressionLayer("Name","regressionoutput")
            ];
        CNN = layerGraph(layers);% creo la rete neurale
    case 2
        layers = [
            sequenceInputLayer(NumInputs,'Name','Sequence Input Layer',"Normalization","rescale-zero-one")% prende i dati in input
        
            convolution1dLayer(filterSize,neurons,"Name","conv1d","Padding",'same')
            batchNormalizationLayer("Name","batchnorm")
            reluLayer("Name","Rectified Linear Unit Layer 1")
        
%             convolution1dLayer(N_conv,neurons,"Name","conv1d2","Padding",'same')
%             batchNormalizationLayer("Name","batchnorm2")
%             reluLayer("Name","Rectified Linear Unit Layer 2")
        
            lstmLayer(numHiddenUnits,'Name','Long-Short Term Layer') % long-short term layer (legami tra le variabili in input)
            reluLayer('Name','Rectified Linear Unit Layer 3') % Applica la funzione ReLu
        
            fullyConnectedLayer(NumOutputs,'Name','Fully Connected Layer')% Restituisce l'output (PM10 giorno successivo) 
            regressionLayer('Name','Regression Layer') % Minimizza il HMSE (half mean squared error)
            ];
        CNN = layerGraph(layers);% creo la rete neurale
end

end