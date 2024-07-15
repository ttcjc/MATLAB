%% Mie Scattering v1.0
% ----
% Lorem ipsum
% ----
% Usage: mieOut = mieScattering(m, x)


%% Changelog

% v1.0  -  Initial Commit


%% Main Function

function mieOut = mieScattering(m, x)
    
    % Avoid Negative Size Parameters
    mieOut = [];
    
    if x > 0
        nmax = round(2 + x + 4 * x^(1 / 3));
        
        n1 = nmax - 1;
        
        n = (1:nmax);cn = 2 * n + 1;
        
        c1n = n .* (n + 2) ./ (n + 1);
        
        c2n = cn ./ n ./ (n + 1);
        
        x2 = x*x;
        
        f = mieCoeffs(m,x); % [an; bn; cn; dn]
        
        anp = (real(f(1,:)));
        
        anpp = (imag(f(1,:)));
        
        bnp = (real(f(2,:)));
        
        bnpp = (imag(f(2,:)));
        
        g1(1:4,nmax) = [0; 0; 0; 0];
        
        g1(1,1:n1) = anp(2:nmax);
        
        g1(2,1:n1) = anpp(2:nmax);
        
        g1(3,1:n1) = bnp(2:nmax);
        
        g1(4,1:n1) = bnpp(2:nmax);
        
        dn = cn .* (anp + bnp);
        
        q = sum(dn);
        
        qext = 2 * q / x2;
        
        en = cn .* (anp .* anp + anpp .* anpp + bnp .* bnp + bnpp .* bnpp);
        
        q = sum(en);
        
        qsca = 2 * q / x2;
        
        qabs = qext - qsca;
        
        fn = (f(1,:) - f(2,:)) .* cn;
        
        gn = (-1).^n;
        
        f(3,:) = fn .* gn;
        
        q = sum(f(3,:));
        
        qb = q * q' / x2;
        
        asy1 = c1n .* (anp .* g1(1,:) + anpp .* g1(2,:) + bnp .* g1(3,:) + bnpp .* g1(4,:));
        
        asy2 = c2n .* (anp .* bnp + anpp .* bnpp);
        
        asy = 4 / x2 * sum(asy1 + asy2) / qsca;
        
        qratio = qb / qsca;
        
        mieOut = [qext qsca qabs qb asy qratio];
    end
    
end


%% Local Functions

function coeffsOut = mieCoeffs(m, x)

    nmax = round(2 + x + 4 * x^(1 / 3));
    
    n = (1:nmax);
    
    nu  =  (n + 0.5);
    
    z = m .* x;
    
    m2 = m .* m;
    
    sqx =  sqrt(0.5 * pi ./ x); sqz =  sqrt(0.5 * pi ./ z);
    
    bx  =  besselj(nu, x) .* sqx;
    
    bz  =  besselj(nu, z) .* sqz;
    
    yx  =  bessely(nu, x) .* sqx;
    
    hx  =  bx + 1i * yx;
    
    b1x = [sin(x) / x, bx(1:nmax - 1)];
    
    b1z = [sin(z) / z, bz(1:nmax - 1)];
    
    y1x = [ - cos(x) / x, yx(1:nmax - 1)];
    
    h1x =  b1x + 1i * y1x;
    
    ax  =  x .* b1x - n .* bx;
    
    az  =  z .* b1z - n .* bz;
    
    ahx =  x .* h1x - n .* hx;
    
    an  =  (m2 .* bz .* ax - bx .* az) ./ (m2 .* bz .* ahx - hx .* az);
    
    bn  =  (bz .* ax - bx .* az) ./ (bz .* ahx - hx .* az);
    
    cn  =  (bx .* ahx - hx .* ax) ./ (bz .* ahx - hx .* az);
    
    dn  =  m .* (bx .* ahx - hx .* ax) ./ (m2 .* bz .* ahx - hx .* az);
    
    coeffsOut = [an; bn; cn; dn];
    
end