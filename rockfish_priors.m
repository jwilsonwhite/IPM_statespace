function rockfish_priors

% Estimate priors on F from stock assessments for blue & gopher rockfish
% (Key et al. 2005, 2007)


% Data from p. 7 of blue rockfish assessment
Blue_90s = [0.155 + 0.083 + 0;...
            0.121 + 0.058 + 0];
        
Blue_00s = [0.089 + 0.039 + 0.001;...
            0.062 + 0.015 + 0;...
            0.043 + 0.016 + 0;...
            0.051 + 0.013 + 0;...
            0.069 + 0.06 + 0;...
            0.045 + 0.009 + 0;...
            0.046 + 0.012 + 0;...
            0.087 + 0.011 + 0];
        
% Data from p. 5 of gopher rockfish assessment
Gopher_90s = [0.061 + 0.137;...
              0.085 + 0.14;...
              0.1 + 0.16;...
              0.095 + 0.191;...
              0.064 + 0.178;...
              0.098 + 0.091;...
              0.087 + 0.056;...
              0.064 + 0.049;...
              0.046 + 0.043];
          
Gopher_00s = [0.031 + 0.049;...
              0.035 + 0.079;...
              0.026 + 0.058;...
              0.01 + 0.102;...
              0.012 + 0.028];
          
% Blue calculation:
M = 0.14;  % natural mortality
H = Blue_90s; %harvest fraction
F0 = 0.1; % initial guess at F


for h = 1:length(H)
    Fx = @(F) abs((1-exp(-(F+M)))*(F/(F+M))-H(h)); % function to be minimized
F_Blue_90s(h) = fminsearch(@(f) Fx(f), F0); % returns the value of F
end
        
H = Blue_00s; %harvest fraction

for h = 1:length(H)
    Fx = @(F) abs((1-exp(-(F+M)))*(F/(F+M))-H(h)); % function to be minimized
F_Blue_00s(h) = fminsearch(@(f) Fx(f), F0); % returns the value of F
end

% Gophers:
M = 0.2;  % natural mortality
H = Gopher_90s; %harvest fraction
F0 = 0.1; % initial guess at F

for h = 1:length(H)
    Fx = @(F) abs((1-exp(-(F+M)))*(F/(F+M))-H(h)); % function to be minimized
F_Gopher_90s(h) = fminsearch(@(f) Fx(f), F0); % returns the value of F
end
        
H = Gopher_00s; %harvest fraction

for h = 1:length(H)
    Fx = @(F) abs((1-exp(-(F+M)))*(F/(F+M))-H(h)); % function to be minimized
F_Gopher_00s(h) = fminsearch(@(f) Fx(f), F0); % returns the value of F
end
      
mean(F_Blue_90s)
std(F_Blue_90s)

mean(F_Blue_00s)
std(F_Blue_00s)

mean(F_Gopher_90s)
std(F_Gopher_90s)

mean(F_Gopher_00s)
std(F_Gopher_00s)


          
              
              