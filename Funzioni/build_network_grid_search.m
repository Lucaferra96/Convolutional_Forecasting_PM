function [ANN] = build_network_grid_search(NumInputs, NumOutputs, net_type, neurons, filterSize, numHiddenUnits)

switch net_type
    case 0

        layers1 = [
            sequenceInputLayer(NumInputs,'Name','Sequence Input Layer',"Normalization","zscore")% prende i dati in input

            lstmLayer(numHiddenUnits) % long-short term layer (legami tra le variabili in input)
            dropoutLayer(0.2)

            fullyConnectedLayer(NumOutputs,'Name','Fully Connected Layer')% Restituisce l'output (PM10 giorno successivo) 
            regressionLayer('Name','Regression Layer') % Minimizza il HMSE (half mean squared error)
            ];
        % layers2 = [
        %     sequenceInputLayer(NumInputs, 'Name', 'Sequence Input Layer', "Normalization", "zscore")
        %     lstmLayer(numHiddenUnits, 'OutputMode', 'sequence')
        %     lstmLayer(numHiddenUnits, 'OutputMode', 'sequence')
        %     dropoutLayer(0.2)
        %     fullyConnectedLayer(NumOutputs, 'Name', 'Fully Connected Layer')
        %     regressionLayer('Name', 'Regression Layer')
        %     ];
        % layers3 = [
        %     sequenceInputLayer(NumInputs, 'Name', 'Sequence Input Layer', "Normalization", "zscore")
        %     bilstmLayer(numHiddenUnits, 'OutputMode', 'sequence')
        %     bilstmLayer(numHiddenUnits, 'OutputMode', 'sequence')
        %     dropoutLayer(0.2)
        %     fullyConnectedLayer(NumOutputs, 'Name', 'Fully Connected Layer')
        %     regressionLayer('Name', 'Regression Layer')
        %     ];

        ANN = {
            % Architecture 1: Shallow CNN
            layerGraph(layers1)
            
            % % Architecture 2: Deeper CNN
            % layerGraph(layers2)
            % 
            % % Architecture 3: 1D CNN with Pooling and Dropout
            % layerGraph(layers3)
            };

    case 1
        layers1 = [
            sequenceInputLayer(NumInputs,"Name","sequence","Normalization","rescale-zero-one")        
            convolution1dLayer(filterSize,neurons,"Name","conv1d","Padding",'same')
            batchNormalizationLayer("Name","batchnorm")
            reluLayer("Name","Rectified Linear Unit Layer")
            fullyConnectedLayer(NumOutputs,"Name","fc")
            regressionLayer("Name","regressionoutput")
            ];
        % layers2 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     convolution1dLayer(filterSize, neurons, "Name", "conv1d_1", "Padding", 'same')
        %     batchNormalizationLayer("Name", "batchnorm_1")
        %     reluLayer("Name", "ReLU_1")
        %     convolution1dLayer(filterSize, neurons, "Name", "conv1d_2", "Padding", 'same')
        %     batchNormalizationLayer("Name", "batchnorm_2")
        %     reluLayer("Name", "ReLU_2")
        %     convolution1dLayer(filterSize, neurons, "Name", "conv1d_3", "Padding", 'same')
        %     batchNormalizationLayer("Name", "batchnorm_3")
        %     reluLayer("Name", "ReLU_3")
        %     fullyConnectedLayer(NumOutputs, "Name", "fc")
        %     regressionLayer("Name", "regressionoutput")
        %     ];
        % layers3 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     convolution1dLayer(filterSize, neurons, "Name", "conv1d", "Padding", 'same')
        %     batchNormalizationLayer("Name", "batchnorm")
        %     reluLayer("Name", "ReLU")
        %     maxPooling1dLayer(1, "Stride", 1, "Name", "maxpooling")
        %     dropoutLayer(0.5, "Name", "dropout")
        %     fullyConnectedLayer(NumOutputs, "Name", "fc")
        %     regressionLayer("Name", "regressionoutput")
        %     ];
        % layers4 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     convolution1dLayer(filterSize, neurons, "Name", "conv1d", "Padding", 'same')
        %     batchNormalizationLayer("Name", "batchnorm")
        %     reluLayer("Name", "ReLU")
        %     fullyConnectedLayer(64, "Name", "fc_1")
        %     reluLayer("Name", "ReLU_1")
        %     fullyConnectedLayer(32, "Name", "fc_2")
        %     reluLayer("Name", "ReLU_2")
        %     fullyConnectedLayer(NumOutputs, "Name", "fc_3")
        %     regressionLayer("Name", "regressionoutput")
        %     ];
        % layers5 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     convolution1dLayer(filterSize, neurons, "Name", "conv1d_1", "Padding", 'same')
        %     batchNormalizationLayer("Name", "batchnorm_1")
        %     reluLayer("Name", "ReLU_1")
        %     convolution1dLayer(filterSize, neurons, "Name", "conv1d_2", "Padding", 'same')
        %     batchNormalizationLayer("Name", "batchnorm_2")
        %     reluLayer("Name", "ReLU_2")
        %     %additionLayer(2, "Name", "addition")
        %     fullyConnectedLayer(NumOutputs, "Name", "fc")
        %     regressionLayer("Name", "regressionoutput")
        %     ];

        ANN = {
            % Architecture 1: Shallow CNN
            layerGraph(layers1)
            
            % % Architecture 2: Deeper CNN
            % layerGraph(layers2)
            % 
            % % Architecture 3: 1D CNN with Pooling and Dropout
            % layerGraph(layers3)
            % 
            % % Architecture 4: 1D CNN with Multiple Fully Connected Layers
            % layerGraph(layers4)
            
            % Architecture 5: 1D CNN with Residual Connections
            % layerGraph(layers5)
            };

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
        ANN = layerGraph(layers);% creo la rete neurale

        case 3

        % layers1 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     fullyConnectedLayer(neurons(1), 'Name', 'FullyConnectedLayer')
        %     fullyConnectedLayer(neurons(2), 'Name', 'FullyConnectedLayer2')
        %     regressionLayer('Name', 'RegressionLayer')
        %     ];
        % layers2 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     fullyConnectedLayer(neurons, 'Name', 'FullyConnectedLayer')
        %     reluLayer('Name', 'ReLULayer')
        %     dropoutLayer(0.2, 'Name', 'DropoutLayer')
        %     fullyConnectedLayer(neurons, 'Name', 'FullyConnectedLayer2')
        %     reluLayer('Name', 'ReLULayer2')
        %     fullyConnectedLayer(NumOutputs, 'Name', 'OutputLayer')
        %     regressionLayer('Name', 'RegressionLayer')
        %     ];
        % layers3 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     fullyConnectedLayer(neurons, 'Name', 'FullyConnectedLayer')
        %     batchNormalizationLayer('Name', 'BatchNormalizationLayer')
        %     reluLayer('Name', 'ReLULayer')
        %     fullyConnectedLayer(neurons, 'Name', 'FullyConnectedLayer2')
        %     batchNormalizationLayer('Name', 'BatchNormalizationLayer2')
        %     reluLayer('Name', 'ReLULayer2')
        %     fullyConnectedLayer(NumOutputs, 'Name', 'OutputLayer')
        %     regressionLayer('Name', 'RegressionLayer')
        %     ];
        % layers4 = [
        %     sequenceInputLayer(NumInputs, "Name", "sequence", "Normalization", "rescale-zero-one")
        %     fullyConnectedLayer(neurons, 'Name', 'FullyConnectedLayer')
        %     reluLayer('Name', 'ReLULayer')
        %     fullyConnectedLayer(neurons, 'Name', 'FullyConnectedLayer2')
        %     reluLayer('Name', 'ReLULayer2')
        %     fullyConnectedLayer(neurons, 'Name', 'FullyConnectedLayer3')
        %     reluLayer('Name', 'ReLULayer3')
        %     fullyConnectedLayer(NumOutputs, 'Name', 'OutputLayer')
        %     regressionLayer('Name', 'RegressionLayer')
        %     ];
        

        % ANN = {
        %     % Architecture 1: Basic Feed Forward Neural Network
        %     layerGraph(layers1)
            
            % % Architecture 2: FNN with Dropouts
            % layerGraph(layers2)
            % 
            % % Architecture 3: FNN with Batch Normalization
            % layerGraph(layers3)
            % 
            % % Architecture 4: Deeper FNN
            % layerGraph(layers4)
            };

end

end