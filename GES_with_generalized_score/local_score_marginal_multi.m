
function [score] = local_score_marginal_multi(Data, Xi,PAi,parameters)
% calculate the local score by negative marginal likelihood
% based on a regression model in RKHS
% for variables with multi-variate dimensions

% parameters.d_label: index of each variable

T = size(Data,1);
X = Data(:,parameters.dlabel{Xi});
dX = size(X,2);

% set the kernel for X
GX = sum((X.*X),2);
Q = repmat(GX,1,T);
R = repmat(GX',T,1);
dists = Q + R - 2*(X*X');
dists = dists-tril(dists);
dists=reshape(dists,T^2,1);
widthX = sqrt(0.5*median(dists(dists>0)));
widthX = widthX*2.5; % kernel width
theta = 1/(widthX^2);
H =  eye(T) - ones(T,T)/T;
Kx = kernel(X, X, [theta,1]);
Kx = H * Kx * H;

Thresh = 1E-5;
[eig_Kx, eix] = eigdec((Kx+Kx')/2, min(400, floor(T/4))); % /2
IIx = find(eig_Kx > max(eig_Kx) * Thresh);
eig_Kx = eig_Kx(IIx);
eix = eix(:,IIx);

if(~isempty(PAi))    
    % set the kernel for PA
    PA_all = [];
    widthPA_all = [];
    for m = 1:length(PAi)
        PA = Data(:,parameters.dlabel{PAi(m)});
        PA_all = [PA_all, PA];
        G = sum(PA.*PA,2);
        Q = repmat(G,1,T);
        R = repmat(G',T,1);
        dists = Q + R - 2*PA*PA';
        dists = dists-tril(dists);
        dists=reshape(dists,T^2,1);
        widthPA = sqrt(0.5*median(dists(dists>0)));
        widthPA_all = [widthPA_all, widthPA*ones(1,length(parameters.dlabel{PAi(m)}))];
    end
    widthPA_all = widthPA_all*2.5; % kernel width
    
    covfunc = {'covSum', {'covSEard','covNoise'}};
    logtheta0 = [log(widthPA_all'); 0; log(sqrt(0.1))];
    [logtheta, fvals, iter] = minimize(logtheta0, 'gpr_multi_new', -300, covfunc, PA_all, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    
    [nlml dnlml] = gpr_multi_new(logtheta, covfunc, PA_all, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    
else
    
    covfunc = {'covSum', {'covSEard','covNoise'}};
    PA = zeros(T,1);
    logtheta0 = [100; 0; log(sqrt(0.1))];
    [logtheta, fvals, iter] = minimize(logtheta0, 'gpr_multi_new', -300, covfunc, PA, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    [nlml dnlml] = gpr_multi_new(logtheta, covfunc, PA, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    
end

score = nlml; % negative log-likelihood











