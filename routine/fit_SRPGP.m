function [ fit ] = fit_SRPGP( P, Y, regr, corr, lob, upb )
    
    H = regr((1:length(Y))'/length(Y));
    data = struct('regr',regr, 'corr',corr, 'P', P, 'Y',Y, 'H',H, 'HH',H'*H, 'YY',Y'*Y, 'HY',H'*Y);
    data.fftY = fft(Y)/sqrt(length(Y));    data.fftH = fft(full(H))/sqrt(length(Y));
    P = P(:);    index = (P == round(P));    
    data.intP = P(index);    data.iterP = data.intP;
    [data.Yb, data.Ys] = segment_mean(Y, data.intP);
    
    q = size(H,2);
    data.Hb = ones(max(data.intP), length(data.intP), q); % mean of G
    data.Hs = ones(max(data.intP), length(data.intP), q); % tail of G
    if q > 1
        for j = 1 : q
            [data.Hb(:,:,j), data.Hs(:,:,j)] = segment_mean(H(:,j),data.intP);
        end
    end
    data.Hb = permute(data.Hb, [1,3,2]);
    data.Hs = permute(data.Hs, [1,3,2]);
    
    % init parameters
    u = linspace(0,1,5);    [x1, x2] = meshgrid(u(2:end-1));
    thetas = lob + (upb-lob).*[x1(:),x2(:)];
    objs = nan(size(thetas,1),1);    fits = [];
    for i = 1 : size(thetas,1)
        [objs(i), fit] = likelihood_at_integer(thetas(i,:), data);
        fits = [fits; fit];
    end
    [~, index] = min(objs);
    theta0 = thetas(index,:);
    [~, index] = max(sum([fits(:).likelihood],2)); % optimal period
    data.iterP = fit.P(index); % figure; plot(sum([fits(:).likelihood],2));
    [theta, ~, fit, ~] = boxmin(@likelihood_at_integer, theta0, lob, upb, data);
    
    data.iterP = data.intP;
    [~, fit] = likelihood_at_integer(theta, data); % decimal?

    fit.theta0 = theta0;    fit.thetahat = theta;
    [~, index] = max(fit.likelihood);    fit.period = fit.P(index);
    fit.regr = regr;    fit.corr = corr;
    fit.trend = fit.H*fit.alpha;
end

function [obj, fit] = likelihood_at_integer(para, data) % likelihood at integer
    delta2 = para(1).^2;     theta = para(2:end); % Initialize
    n = length(data.Y);    q = size(data.H,2);    rho = 1;
    likelihood = nan(length(data.P),1);
    sigmas = nan(size(likelihood));    iterNum = nan(length(data.iterP),1);
    betas = nan(n, length(data.P));    alphas = nan(q, length(data.P));
    deltaTild2 = delta2 + 1/rho; % according to the definition of K   
    threshold = quantile(abs(data.fftY), 0.9);    lambda = threshold * rho;
    penality = nan(length(data.P),1);
    for i = 1 : length(data.iterP)
        p = data.iterP(i);    k = floor(n/p);    p1 = n - p*k;
        Yb = data.Yb(1:p,data.intP==p);    Ys = data.Ys(1:p1,data.intP==p);
        Hb = data.Hb(1:p,:,data.intP==p);    Hs = data.Hs(1:p1,:,data.intP==p);
        HH = data.HH-Hs'*Hs;    HY = data.HY-Hs'*Ys;
        R0.c = data.corr(0, theta, p, (1:p)', 1);    R0.eig = abs(fft(R0.c));
        Rdtild.c = k*R0.c/deltaTild2 + eye(p,1);    Rdtild.eig = abs(fft(Rdtild.c));
        Rd.c = k*R0.c/delta2 + eye(p,1);    Rd.eig = real(fft(Rd.c));

        fftHb = fft(Hb);    fftYb = fft(Yb);    Hbl = ifft(Rd.eig .\ fftHb);
        SHH = (HH+k*Hb'*(Hbl-Hb))/delta2;     SHY = (HY+k*(Hbl-Hb)'*Yb)/delta2;
        R0Rd_.eig = (k/delta2) * R0.eig ./ Rd.eig;    R0Rd_.c = ifft(R0Rd_.eig);
        if p1 ==0
            alpha = SHH \ SHY; % init alpha
        else
            R0Rd_Hb = ifft(R0Rd_.eig .* fftHb);    R0Rd_Yb = ifft(R0Rd_.eig .* fftYb);
            rpi = R0.c - ifft(R0Rd_.eig .* R0.eig) + delta2*eye(p,1);
            Gd = Hs - R0Rd_Hb(1:p1,:);    Yd = Ys - R0Rd_Yb(1:p1);
            S = toep_chol(rpi(1:p1));
            Gdl = S \ Gd;    GGd = Gdl'*Gdl;
            Ydl = S \ Yd;    GYd = Gdl'*Ydl;
            alpha = (SHH + GGd) \ (SHY + GYd);
            
            R0Rdtild_.eig = R0.eig ./ Rdtild.eig;    R0Rdtild_.c = ifft(R0Rdtild_.eig);
            rpi_tild = R0.c - ifft(R0Rdtild_.eig .* R0.eig)*k/deltaTild2 + deltaTild2*eye(p,1);
            Stild = toep_chol(rpi_tild(1:p1));
        end
        
        % compute Z without F*beta
        Zb = Yb - Hb*alpha;    Zt = Ys - Hs*alpha;
        Zb_hat = ifft(R0Rd_.eig .* fft(Zb)); % filtering
        if p1>0 % figure(100); plot(Zb_hat)
            Zd = Zt - Zb_hat(1:p1); % Z dot
            Idoc.c = eye(p,1) - R0Rd_.c; % Idot matrix
            IdR0.eig = fft(Idoc.c) .* R0.eig; % IdR.c = ifft(IdR.eig);

            Zdl = S \ Zd;    Pi_Zd = S' \ Zdl;    Pi_Zd(p,:) = 0;
            Zb_hat = Zb_hat + ifft(IdR0.eig .* fft(Pi_Zd));
            Z = [repmat(Zb_hat, k,1); Zb_hat(1:p1)];
        else
            Z = repmat(Zb_hat, k,1);
        end
        
        beta = zeros(size(data.Y));    beta0 = ones(size(data.Y));    TempMat = zeros(p, k); 
        v = data.fftY-data.fftH*alpha-fft(Z)/sqrt(n);    u = zeros(size(data.Y));
        count = 0;
        while max(abs(beta-beta0)) > 1e-6 % figure(200); plot(data.Y); hold on;  
            beta0 = beta;
            ifftV_U = ifft(v-u)*sqrt(n);
            res = (data.Y - data.H*alpha - ifftV_U)/rho;    TempMat(1:k*p) = res(1:k*p);
            rb = mean(TempMat,2);    rt = res(end-p1+1:end);
            Pi_rd = zeros(p, 1);    R0Rdtild_Pi_rd = zeros(p, 1);
            if p1>0
                R0Rdtild_rb = ifft(R0Rdtild_.eig .* fft(rb));
                rd = rt - R0Rdtild_rb(1:p1)*k/deltaTild2; % r dot
                rdl = Stild \ rd;    Pi_rd = Stild' \ rdl;    Pi_rd(p,:) = 0;
                R0Rdtild_Pi_rd = ifft(R0Rdtild_.eig .* fft(Pi_rd));
            end
            Rd_rb = ifft(fft(rb)./Rdtild.eig);
            K_r = [(res(1:k*p) + repmat(Rd_rb-rb-R0Rdtild_Pi_rd, k, 1))/deltaTild2; Pi_rd(1:p1)];
            beta = fft(K_r)/sqrt(n) + v - u;
            ifft_beta = K_r + ifftV_U;
            Y = data.Y - K_r - ifftV_U;    TempMat(1:k*p) = Y(1:k*p); 
            Yb = mean(TempMat,2);    Ys = Y(end-p1+1:end);    fftYb = fft(Yb);  
            HY = data.H'*Y - Hs'*Ys;    SHY = (HY+k*(Hbl-Hb)'*Yb)/delta2;
            if p1==0 %disp(alpha');
                alpha = SHH \ SHY; % init alpha
            else
                R0Rd_Yb = ifft(R0Rd_.eig .* fftYb);
                Yd = Ys - R0Rd_Yb(1:p1);
                Ydl = S \ Yd;    GYd = Gdl'*Ydl;
                alpha = (SHH + GGd) \ (SHY + GYd);
            end
            v = soft(beta + u, threshold);
            u = u + rho*(beta - v);
 
            count = count + 1;
        end % figure(100); hold on; plot(ifft(beta)*sqrt(n))
        iterNum(i) = count;
        
        Z = data.Y - data.H*alpha - ifft_beta;
        TempMat(1:k*p) = Z(1:k*p);    Zb = mean(TempMat,2);    Zt = Z(1+k*p:end);
        Zb_hat = ifft(R0Rd_.eig .* fft(Zb));
        if p1>0 % figure(100); plot(Zb_hat)
            Zd = Zt - Zb_hat(1:p1); % Z dot
            Idoc.c = eye(p,1) - R0Rd_.c; % Idot matrix
            IdR0.eig = fft(Idoc.c) .* R0.eig; % IdR.c = ifft(IdR.eig);

            Zdl = S \ Zd;    Pi_Zd = S' \ Zdl;    Pi_Zd(p,:) = 0;
            Zb_hat = Zb_hat + ifft(IdR0.eig .* fft(Pi_Zd));
        end
        
        if p1>0 % compute Z
            Zs(:,p==data.P) = [repmat(Zb_hat, k,1); Zb_hat(1:p1)];
        else 
            Zs(:,p==data.P) = repmat(Zb_hat, k,1);
        end

        ZtZ = sum(sum(TempMat.*TempMat));
        if p1 == 0
            sigma2 = (ZtZ+k*Zb'*(ifft(Rd.eig.\fft(Zb))-Zb))/delta2/n;
            likelihood(p==data.P) = (n*log(sigma2*delta2) + sum(log(Rd.eig)) + n + n*log(2*pi))/2;
        else
            sigma2 = ((ZtZ+k*Zb'*(ifft(Rd.eig.\fft(Zb))-Zb))/delta2 + Zdl'*Zdl)/n;
            likelihood(p==data.P) = (n*log(sigma2) + k*p*log(delta2) + sum(log(Rd.eig)) + sum(2*log(diag(S))) + n + n*log(2*pi))/2;
        end
        likelihood(p==data.P) = likelihood(p==data.P) + lambda*sum(abs(beta))/sigma2;
        penality(p==data.P) = lambda*sum(abs(beta));
        sigmas(p==data.P) = sqrt(sigma2);
        betas(:,p==data.P) = beta;    alphas(:,p==data.P) = alpha;      
    end % figure; plot(-likelihood);  figure; plot(sigmas);  figure; plot(power)
%     disp(sum(iterNum))

    [~, index] = min(likelihood);
    obj = likelihood(index);
    if  nargout > 1
        fit = struct('P',data.P, 'Y',data.Y, 'alpha',alphas(:,index), 'beta',betas(:,index), ...
            'sigma',sigmas(index), 'Z',Zs(:,index), 'H',data.H, 'likelihood',-likelihood);
    end
end

function y = soft(x, T)
%   x : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%   y : output of soft thresholding
    y = max(1 - T./abs(x), 0) .* x;
end

function [Yb, Ys] = segment_mean(Y, P)
    n = length(Y);    Y(end+max(P)) = 0; 
    Yseg = nan(max(P), length(P)); % sum of Y
    Ys = zeros(max(P), length(P)); % tail of Y
    for i = length(P) : -1 : 1
        if ~isnan(Yseg(1,i)), continue; end
        p = P(i);    p1 = n - p*floor(n/p);
        Yseg(1:p,i) = sum(reshape(Y(1:p*ceil(n/p)), p, ceil(n/p)), 2);
        Ys(1:p1,i)  = Y(n-p1+1:n);
        for j = 1 : i-1
            if ~isnan(Yseg(1,j)), continue; end
            q = P(j);    q1 = n - q*floor(n/q);
            if mod(p,q)==0
                Yseg(1:q,j) = sum(reshape(Yseg(1:p,i), q, p/q), 2);
                Ys(1:q1,j)  = Y(n-q1+1:n);
            end
        end
    end
    Yb = (Yseg - Ys) ./ floor(n./P');
end

function  [t, f, fit, perf] = boxmin(objfunc, t0, lo, up, data)
%BOXMIN  Minimize with positive box constraints

    % Initialize
    [t, f, fit, itdata] = start(objfunc, t0, lo, up, data);
    if  ~isinf(f)
      % Iterate
      p = length(t);
      if  p <= 2,  kmax = 2; else,  kmax = min(p,4); end
      for  k = 1 : kmax
        th = t;
        [t, f, fit, itdata] = explore(objfunc, t, f, fit, itdata, data);
        [t, f, fit, itdata] = move(objfunc, th, t, f, fit, itdata, data);
      end
    end
    perf = struct('nv',itdata.nv, 'perf',itdata.perf(:,1:itdata.nv));
end

function  [t, f, fit, itdata] = start(objfunc, t0, lo, up, data)
% Get starting point and iteration dataameters

    % Initialize
    t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
    D = 2 .^ ([1:p]'/(p+2));
    ee = find(up == lo);  % Equality constraints
    if  ~isempty(ee)
      D(ee) = ones(length(ee),1);   t(ee) = up(ee); 
    end
    ng = find(t < lo | up < t);  % Free starting values
    if  ~isempty(ng)
      t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  % Starting point
    end
    ne = find(D ~= 1);

    % Check starting point and initialize performance info
    [f  fit] = objfunc(t,data);   nv = 1;
    itdata = struct('D',D, 'ne',ne, 'lo',lo, 'up',up, ...
      'perf',zeros(p+2,200*p), 'nv',1);
    itdata.perf(:,1) = [t; f; 1];
    if  isinf(f)    % Bad dataameter region
      return
    end

    if  length(ng) > 1  % Try to improve starting guess
      d0 = 16;  d1 = 2;   q = length(ng);
      th = t;   fh = f;   jdom = ng(1);  
      for  k = 1 : q
        j = ng(k);    fk = fh;  tk = th;
        DD = ones(p,1);  DD(ng) = repmat(1/d1,q,1);  DD(j) = 1/d0;
        alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
        v = DD .^ alpha;   tk = th;
        for  rept = 1 : 4
          tt = tk .* v; 
          [ff  fitt] = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 1];
          if  ff <= fk 
            tk = tt;  fk = ff;
            if  ff <= f
              t = tt;  f = ff;  fit = fitt; jdom = j;
            end
          else
            itdata.perf(end,nv) = -1;   break
          end
        end
      end % improve

      % Update Delta  
      if  jdom > 1
        D([1 jdom]) = D([jdom 1]); 
        itdata.D = D;
      end
    end % free variables

    itdata.nv = nv;
end

function  [t, f, fit, itdata] = explore(objfunc, t, f, fit, itdata, data)
% Explore step

    nv = itdata.nv;   ne = itdata.ne;
    for  k = 1 : length(ne)
      j = ne(k);   tt = t;   DD = itdata.D(j);
      if  t(j) == itdata.up(j)
        atbd = 1;   tt(j) = t(j) / sqrt(DD);
      elseif  t(j) == itdata.lo(j)
        atbd = 1;  tt(j) = t(j) * sqrt(DD);
      else
        atbd = 0;  tt(j) = min(itdata.up(j), t(j)*DD);
      end
      [ff  fitt] = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 2];
      if  ff < f
        t = tt;  f = ff;  fit = fitt;
      else
        itdata.perf(end,nv) = -2;
        if  ~atbd  % try decrease
          tt(j) = max(itdata.lo(j), t(j)/DD);
          [ff  fitt] = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 2];
          if  ff < f
            t = tt;  f = ff;  fit = fitt;
          else
            itdata.perf(end,nv) = -2;
          end
        end
      end
    end % k

    itdata.nv = nv;
end

function  [t, f, fit, itdata] = move(objfunc, th, t, f, fit, itdata, data)
% Pattern move

    nv = itdata.nv;   ne = itdata.ne;   p = length(t);
    v = t ./ th;
    if  all(v == 1)
      itdata.D = itdata.D([2:p 1]).^.2;
      return
    end

    % Proper move
    rept = 1;
    while  rept
      tt = min(itdata.up, max(itdata.lo, t .* v));  
      [ff  fitt] = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 3];
      if  ff < f
        t = tt;  f = ff;  fit = fitt;
        v = v .^ 2;
      else
        itdata.perf(end,nv) = -3;
        rept = 0;
      end
      if  any(tt == itdata.lo | tt == itdata.up), rept = 0; end
    end

    itdata.nv = nv;
    itdata.D = itdata.D([2:p 1]).^.25;
end