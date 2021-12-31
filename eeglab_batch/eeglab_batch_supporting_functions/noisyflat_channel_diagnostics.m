function S = noisyflat_channel_diagnostics(S,EEG)

% (inverse) variance over data points/trials per channels
varv = std(reshape(EEG.data,size(EEG.data,1),[]),[],2);
invvar = 1./varv;

% (inverse) variance over channels
varv_chan = std(varv);
invvar_chan = std(invvar);

% outputs
S.prep.chandiag.var_chan = varv_chan;
S.prep.chandiag.invvar_chan = invvar_chan;
