function rr = geometry_trpv3_v2(n)
tt = linspace(0, 2*pi, n);
ecce = 0.8;    % ellipse eccentricity
para = 0.6;    % ellipse parameter
theta = -pi/2;   % angle origin
r = 0.45;    % ellipse radius correction factor

rr1 = r*para./(1-ecce.*cos(tt - theta));
rr2 = 0.5;
rr3 = 0.7 *r*para./(1-ecce.*cos(tt - theta + 3*pi/5));
rr = max(rr1, rr3);
rr = max(rr, rr2);
rr = [rr rr rr];
rr = movmean(rr, round(n/10));
rr = rr(n+1: 2*n);
rr0 = 0.5*(rr(1)+rr(end));
rr(1) = rr0;
rr(end) = rr0;
end