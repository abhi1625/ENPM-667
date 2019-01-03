function sdot = cart_EOM(s,t,params,A,B,K)
%% Linear Model
sdot = (A-B*K)*s;
end

