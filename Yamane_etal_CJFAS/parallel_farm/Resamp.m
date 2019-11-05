function Out = Resamp(D,n)
Props = linspace(0,1,1000); % a range of proportions of D to try

D_tmp = repmat(D(:)',[length(Props),1]);
               P_tmp = repmat(Props(:),[1,length(D)]);
               % Multiply all the Ds by each possible set of proportions, then round to
               % get back to integers:
               X = round(D_tmp.*P_tmp);
               
               Sums = sum(X,2); % sample size using each prop
               Which = abs(Sums-n);
               Index = find(Which == min(Which),1);
               Prop_answer = Props(Index);
               
               % Prop_answer is the value to use
               
               Out = round(D*Prop_answer); % this is the dataset with the desired sample size
               %Out = poissrnd(Out); % add in Poisson variability.


end % end loop over reps
