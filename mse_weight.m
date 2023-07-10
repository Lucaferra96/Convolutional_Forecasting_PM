function [loss,gradients] = mse_weight(net,X,T)

Y = forward(net,X);

% Weights
W = flip(1:10);
W = W/sum(W);

loss = 


% Compute loss.
loss = crossentropy(Y,T);

% Compute gradients.
gradients = dlgradient(loss,net.Learnables);


end