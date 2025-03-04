% Wave detection for theta band
disp('Wave Detection for theta band ...')
threshold = 99;
xgp = arrayfun(@(s) s.xgptheta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
parameters.thetarhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
thetaWaves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.thetarhoThres);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xftheta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgptheta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wttheta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    thetaWaves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.thetarhoThres);
end

% Wave detection for beta band
disp('Wave Detection for beta band ...')
threshold = 99.9;
xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
parameters.betarhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);

if isfield(IntanBehaviour,'cueHitTrace')
    xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    betaWaves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.betarhoThres);
end
% Waves.wavesHit = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,0.5);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    betaWaves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.betarhoThres);
%     Waves.wavesMiss = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,0.5);
end
if isfield(IntanBehaviour,'missTrace')
    xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.missTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.missTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.missTrace, 'UniformOutput', false);
    betaWaves.wavesFA = detectWaves(xf,xgp,wt,IntanBehaviour.missTrace,parameters,parameters.betarhoThres);
end
if isfield(IntanBehaviour,'hitTrace')
    xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.hitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.hitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.hitTrace, 'UniformOutput', false);
    betaWaves.wavesHitReward = detectWaves(xf,xgp,wt,IntanBehaviour.hitTrace,parameters,parameters.betarhoThres);
end
if isfield(IntanBehaviour,'MIHitTrace')
    xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    betaWaves.wavesMIHit = detectWaves(xf,xgp,wt,IntanBehaviour.MIHitTrace,parameters,parameters.betarhoThres);
end
if isfield(IntanBehaviour,'MIFATrace')
    xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    betaWaves.wavesMIFA = detectWaves(xf,xgp,wt,IntanBehaviour.MIFATrace,parameters,parameters.betarhoThres);
end

% Wave detection for gamma band
disp('Wave Detection for gamma band ...')
threshold = 99.9;
xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
parameters.gammarhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);

if isfield(IntanBehaviour,'cueHitTrace')
    xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    gammaWaves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.gammarhoThres);
end
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    gammaWaves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.gammarhoThres);
end
if isfield(IntanBehaviour,'missTrace')
    xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.missTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.missTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.missTrace, 'UniformOutput', false);
    gammaWaves.wavesFA = detectWaves(xf,xgp,wt,IntanBehaviour.missTrace,parameters,parameters.gammarhoThres);
end
if isfield(IntanBehaviour,'hitTrace')
    xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.hitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.hitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.hitTrace, 'UniformOutput', false);
    gammaWaves.wavesHitReward = detectWaves(xf,xgp,wt,IntanBehaviour.hitTrace,parameters,parameters.gammarhoThres);
end
if isfield(IntanBehaviour,'MIHitTrace')
    xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    gammaWaves.wavesMIHit = detectWaves(xf,xgp,wt,IntanBehaviour.MIHitTrace,parameters,parameters.gammarhoThres);
end
if isfield(IntanBehaviour,'MIFATrace')
    xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    gammaWaves.wavesMIFA = detectWaves(xf,xgp,wt,IntanBehaviour.MIFATrace,parameters,parameters.gammarhoThres);
end