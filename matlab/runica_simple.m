function [weights, sphere] = runica_simple(data, extended)
% function [weights, sphere, meanvar, bias, signs, lrates, activations] = runica_simple(data, extended, pca, stop, maxstep, blockint)
% RUNICA_SIMPLE  Infomax ICA (logistic) with optional extended ICA (tanh) and optional PCA.
%
% Usage
%   [W,S,meanvar,bias,signs,lrates,acts] = runica_simple(data, extended, pca, stop, maxstep, block)
%
% Input
%   data      : (chans x samples)
%   extended  : 0 off, >0 update signs every N blocks, <0 fix #subgauss = -N
%   pca       : 0 off, N keep N PCs
%   stop      : stopping threshold on weight change
%   maxstep   : max passes
%   block     : block size
%
% Notes
%   Unmixing is W*S (if sphering on). If pca on, returned W already includes PCA, and S = eye(original_chans).

if nargin < 1, error('runica_simple: data required'); end
if nargin < 2, extended = 0; end
if nargin < 3, pca = uint32(0); end
if nargin < 4, stop = 0; end
if nargin < 5, maxstep = uint32(0); end
if nargin < 6, blockint = uint32(0); end
[chans, frames] = size(data);
if chans < 2 || frames < chans, error('runica_simple: data too small or rank limited'); end

% Set defaults for unspecified parameters
if stop == 0
    stop = (chans > 32) * 1e-7 + (chans <= 32) * 1e-6;
end
if maxstep == 0
    maxstep = uint32(512);
end
if blockint == 0
    block = ceil(min(5*log(frames), 0.3*frames));
else
    block = double(blockint);
end

% Fixed internal parameters
sphering      = 'on';
weights_init  = 0;
lrate_init    = 0.00065/log(chans);
anneal        = 0.90;
annealdeg     = 60;
useBias       = true;
momentum      = 0;
verbose       = true;
logfile       = [];
rndreset      = false;

% Extended ICA internals
DEFAULT_EXTANNEAL   = 0.98;
extmomentum         = 0.5;
signsbias           = 0.02;
MAX_KURTSIZE        = 6000;
MIN_KURTSIZE        = 2000;
SIGNCOUNT_THRESHOLD = 25;
SIGNCOUNT_STEP      = 2;

% Logging
fid = [];
if ~isempty(logfile)
    fid = fopen(logfile,'w');
    if fid < 0, error('runica_simple: cannot open logfile'); end
end

% Random seed
resetSeed = logical(rndreset);
warning('off','MATLAB:RandStream:ActivatingLegacyGenerators');
%if resetSeed, rand('state',sum(100*clock)); else, rand('state',0); end %#ok<RAND>
warning('on','MATLAB:RandStream:ActivatingLegacyGenerators');

% Remove channel means
rowmeans = mean(data,2);
data = data - rowmeans;

% Optional PCA
urchans = chans;
eigenvectors = [];
ncomps = chans;
if pca && pca ~= 0
    ncomps = double(pca);
    if ncomps < 0, ncomps = size(data,1) + ncomps; end
    if ncomps < 1 || ncomps > chans, error('runica_simple: invalid pca'); end

    % fprintf('Reducing to %d PCs\n', ncomps);
    C = (data*data.') / frames;
    [V,D] = eig(C);
    [d,ix] = sort(diag(D),'descend');
    V = V(:,ix);
    eigenvectors = double(real(V(:,1:ncomps)));
    newdata = eigenvectors.' * data;
    chans = ncomps;
else
    newdata = data;
end

% Sphering
if strcmpi(sphering,'on')
    % fprintf('Computing sphering matrix\n');
    sphere = real(2.0 * inv(sqrtm(double(cov(newdata.')))));
    newdata = sphere * newdata;
elseif strcmpi(sphering,'off')
    sphere = eye(chans);
else
    error('runica_simple: sphering must be on or off');
end

% Init weights
if isequal(weights_init,0)
    weights = eye(ncomps,chans);
else
    weights = weights_init;
    if ~isequal(size(weights),[ncomps,chans]), error('runica_simple: weights size mismatch'); end
end

% Extended ICA configuration
ext_on = 0;
extblocks = 0;
nsub = 1;
kurtsize = min(MAX_KURTSIZE, frames);
if extended ~= 0
    ext_on = 1;
    extblocks = fix(extended);
    anneal = DEFAULT_EXTANNEAL;
    if extblocks < 0
        nsub = -extblocks;
        extblocks = -1; % marker for fixed signs
    else
        if kurtsize < MIN_KURTSIZE
            % fprintf('Warning: kurtosis estimates may be noisy for %d points\n', kurtsize);
        end
    end
end
if nsub > ncomps, error('runica_simple: too many subgaussian components'); end

% Training constants
MAX_WEIGHT         = 1e8;
DEFAULT_BLOWUP     = 1e9;
DEFAULT_BLOWUP_FAC = 0.8;
DEFAULT_RESTART_FAC= 0.9;
MIN_LRATE          = 1e-6;

lrate     = lrate_init;
nochange  = stop;

lastt   = fix((frames/block - 1) * block + 1);

BI      = double(block) * eye(ncomps);
onesrow = ones(1,block);

bias    = zeros(ncomps,1);
lrates  = zeros(1,maxstep);

% Signs init for extended ICA
signvec = repmat(1,[1 ncomps]); %ones(1,ncomps);
signvec(1:nsub) = -1;
signs = diag(signvec);
oldsigns = zeros(size(signs));
signcount = 0;
old_kk = zeros(1,ncomps);

startweights = weights;
prevweights = weights;
prevwtchange = zeros(size(weights));
oldweights = weights;
olddelta = zeros(1,chans*ncomps);
oldchange = Inf;

% fprintf('Beginning ICA training: ncomps=%f, extended=%f, bias=%f\n', ncomps, ext_on, useBias);

step = 0;
while step < maxstep
    timeperm = randperm(frames);

    wts_blowup = false;
    blockno = 1;

    for t = 1:block:lastt
        Xb = double(newdata(:, timeperm(t:t+block-1)));
        if useBias
            U = weights*Xb + bias*onesrow;
        else
            U = weights*Xb;
        end

        if ext_on
            Y = tanh(U);
            weights = weights + lrate * (BI - signs*Y*U.' - U*U.') * weights;
            if useBias
                bias = bias + lrate * sum((-2*Y).',1).';
            end
        else
            Y = 1 ./ (1 + exp(-U));
            weights = weights + lrate * (BI + (1 - 2*Y)*U.') * weights;
            if useBias
                bias = bias + lrate * sum((1 - 2*Y).',1).';
            end
        end

        if momentum > 0
            weights = weights + momentum * prevwtchange;
            prevwtchange = weights - prevweights;
            prevweights = weights;
        end

        if max(abs(weights(:))) > MAX_WEIGHT
            wts_blowup = true;
            break
        end

        % Extended sign updates by kurtosis
        if ext_on && extblocks > 0 && rem(blockno, extblocks) == 0
            if kurtsize < frames
                rp = fix(rand(1,kurtsize) * frames);
                rp(rp == 0) = 1;
                A = weights * double(newdata(:, rp));
            else
                A = weights * double(newdata);
            end
            m2 = mean(A.'.^2).^2;
            m4 = mean(A.'.^4);
            kk = (m4 ./ m2) - 3.0;
            kk = extmomentum*old_kk + (1-extmomentum)*kk;
            old_kk = kk;

            signs(1:size(signs,1),1:size(signs,2)) = diag(sign(kk + signsbias));
            if isequal(signs, oldsigns)
                signcount = signcount + 1;
            else
                signcount = 0;
            end
            oldsigns = signs;
            if signcount >= SIGNCOUNT_THRESHOLD
                extblocks = fix(extblocks * SIGNCOUNT_STEP);
                signcount = 0;
            end
        end

        blockno = blockno + 1;
    end

    if wts_blowup
        lrate = lrate * DEFAULT_RESTART_FAC;
        if lrate < MIN_LRATE
            error('runica_simple: quitting, learning rate too small after blowup');
        end
        % fprintf('Weights blew up, restarting with lrate=%g\n', lrate);
        weights = startweights;
        oldweights = startweights;
        prevweights = startweights;
        prevwtchange = zeros(size(weights));
        bias = zeros(ncomps,1);
        step = 0;
        continue
    end

    % Step metrics
    step = step + 1;
    lrates(step) = lrate;

    dW = weights - oldweights;
    delta = reshape(dW, 1, chans*ncomps);
    change = delta*delta.';

    if step > 2 && isfinite(oldchange) && oldchange > 0
        angledelta = acos( (delta*olddelta.') / sqrt(change*oldchange) );
    else
        angledelta = 0;
    end
    deg = 180/pi * angledelta;

    fprintf('step %g  lrate %g  wchange %.9g  angledelta %.3g deg\n', step, lrate, change, deg);

    % Anneal
    if deg > annealdeg
        lrate = lrate * anneal;
        olddelta = delta;
        oldchange = change;
    elseif step == 1
        olddelta = delta;
        oldchange = change;
    end

    % Stop and blowup guard
    if step > 2 && change < nochange
        break
    elseif change > DEFAULT_BLOWUP
        lrate = lrate * DEFAULT_BLOWUP_FAC;
    end

    oldweights = weights;
end

lrates = lrates(1:step);

% Activations for output
activations = weights * double(newdata);

% If PCA was used, compose matrices back to original channel space
if ~isempty(eigenvectors)
    if strcmpi(sphering,'on')
        weights = weights * sphere * eigenvectors.';
        sphere = eye(urchans);
    else
        weights = weights * eigenvectors.';
        sphere = eye(urchans);
    end
end

% Component variance ranking
if size(weights,1) == size(weights,2)
    winv = inv(weights*sphere);
else
    winv = pinv(weights*sphere);
end
meanvar = sum(winv.^2,1) .* sum(activations.'.^2,1) / ((chans*frames) - 1);

% Sort descending
[meanvar, idx] = sort(meanvar, 'descend');
weights = weights(idx,:);
bias    = bias(idx);
if ext_on
    signs = diag(signs);
    signs = signs(idx);
else
    signs = ones(ncomps,1);
end
activations = activations(idx,:);

if ~isempty(fid), fclose(fid); end

end
