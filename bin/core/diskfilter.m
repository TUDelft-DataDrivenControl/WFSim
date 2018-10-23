function h = diskfilter(varargin)

% Check the number of input arguments.
narginchk(1,3);

p2    = varargin{2};
rad   = p2;
crad  = ceil(rad-0.5);
[x,y] = meshgrid(-crad:crad,-crad:crad);
maxxy = max(abs(x),abs(y));
minxy = min(abs(x),abs(y));
m1 = (rad^2 <  (maxxy+0.5).^2 + (minxy-0.5).^2).*(minxy-0.5) + ...
    (rad^2 >= (maxxy+0.5).^2 + (minxy-0.5).^2).* ...
    sqrt(rad^2 - (maxxy + 0.5).^2);
m2 = (rad^2 >  (maxxy-0.5).^2 + (minxy+0.5).^2).*(minxy+0.5) + ...
    (rad^2 <= (maxxy-0.5).^2 + (minxy+0.5).^2).* ...
    sqrt(rad^2 - (maxxy - 0.5).^2);
sgrid = (rad^2*(0.5*(asin(m2/rad) - asin(m1/rad)) + ...
    0.25*(sin(2*asin(m2/rad)) - sin(2*asin(m1/rad)))) - ...
    (maxxy-0.5).*(m2-m1) + (m1-minxy+0.5)) ...
    .*((((rad^2 < (maxxy+0.5).^2 + (minxy+0.5).^2) & ...
    (rad^2 > (maxxy-0.5).^2 + (minxy-0.5).^2)) | ...
    ((minxy==0)&(maxxy-0.5 < rad)&(maxxy+0.5>=rad))));
sgrid = sgrid + ((maxxy+0.5).^2 + (minxy+0.5).^2 < rad^2);
sgrid(crad+1,crad+1) = min(pi*rad^2,pi/2);
if ((crad>0) && (rad > crad-0.5) && (rad^2 < (crad-0.5)^2+0.25))
    m1  = sqrt(rad^2 - (crad - 0.5).^2);
    m1n = m1/rad;
    sg0 = 2*(rad^2*(0.5*asin(m1n) + 0.25*sin(2*asin(m1n)))-m1*(crad-0.5));
    sgrid(2*crad+1,crad+1) = sg0;
    sgrid(crad+1,2*crad+1) = sg0;
    sgrid(crad+1,1)        = sg0;
    sgrid(1,crad+1)        = sg0;
    sgrid(2*crad,crad+1)   = sgrid(2*crad,crad+1) - sg0;
    sgrid(crad+1,2*crad)   = sgrid(crad+1,2*crad) - sg0;
    sgrid(crad+1,2)        = sgrid(crad+1,2)      - sg0;
    sgrid(2,crad+1)        = sgrid(2,crad+1)      - sg0;
end
sgrid(crad+1,crad+1) = min(sgrid(crad+1,crad+1),1);
h = sgrid/sum(sgrid(:));

end