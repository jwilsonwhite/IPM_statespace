function kxy = mkkern_Tomales(x,y,Par,isJuv,Model_str,Meta,Y,Si)

% adapted from Will's code adapted from Easterling's code, used for Tomales
% oyster IPM

%SURVIVAL PART OF KERNEL
% Determine the model type & get mortality rates 
        switch Model_str.aM_space
            case 'Con'
                aM = exp(Par.aM(Y));
            case 'Grad'
                aM = exp(Par.aM(Y) + Meta.Distances(Si).*Par.aM(end));
            case 'Unim'
                aM = exp(Par.aM(Y) + Par.aM(end-1)*(Meta.Distances(Si)-Par.aM(end)).^2);
        end % end aM_space switch
        
        switch Model_str.jM_space
            case 'Con'
                jM = exp(Par.jM(Y));
            case 'Grad'
                jM = exp(Par.jM(Y) + Meta.Distances(Si).*Par.jM(end));
            case 'Unim'
                jM = exp(Par.jM(Y) + Par.jM(end-1)*(Meta.Distances(Si)-Par.jM(end)).^2);
        end % end jM_space switch
               
       
m = repmat((1-isJuv).*aM + isJuv.*jM,[1 length(x)]); %this is a matrix size x, 
                                                 %mortality for each size 
                                              
p1 = exp(-m); % convert mortality rate to survivorship, iterate over time steps

%GROWTH PART OF KERNEL
%add variability in growth to k
Linf = Par.Growth(Si,1);
k = Par.Growth(Si,3);

%growth
pmean1=Linf - (Linf - x).*exp(-k); % (do not add in x0 for the one-step growth)

%add variability around von Bertalanffy growth
psig1 = pmean1*Par.Growth(Si,4); 

%evaluate growth part of kernel
p2 = normpdf(y, pmean1, psig1);

p1 = max(0,p1);    %to make sure no negatives
p2 = max(0,p2);    %to make sure no negatives

kxy = p1.*p2;
%add fecundity to kernel if a closed population

