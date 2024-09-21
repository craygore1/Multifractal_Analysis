% calculates multifractal spectrum for 2d binary image
% uses a circular support
% objectcolor = color of object of interest ** 0 for black ** ** 1 for white **
% qvals = vector of desired q values
% outputs dq vs q plot and spectrum

function [Dq,myalpha,falpha] = multifractalcoordinate(bwimage,objectcolor,qvals,plots)

sz = size(bwimage);

if objectcolor
    mylog = bwimage == 1;
else
    mylog = bwimage == 0;
end

centerx = floor(sz(2)/2); %setting ellipse axes
centery = floor(sz(1)/2); 
radiusx = centerx;
radiusy = centery;
fullradius = radiusx+radiusy;
[imcols, imrows] = meshgrid(1:sz(2), 1:sz(1));
myellipse = (imrows - centery).^2 ./ radiusy^2 ...
    + (imcols - centerx).^2 ./ radiusx^2 <= 1;
mylog = mylog & myellipse;

xvec = round(linspace(-sz(2)/2,sz(2)/2,sz(2)));
yvec = round(linspace(-sz(1)/2,sz(1)/2,sz(1)));
[Xim,Yim] = meshgrid(xvec,yvec);

[theta,rho] = cart2pol(Xim,Yim); %polarcoordinates of all points

myind = find(mylog == objectcolor); %indices of points of interest
thetaind = theta(myind);
rhoind = rho(myind);

rhoind = rescale(rhoind,0,1);

%%%%%%% Box Counting Loop %%%%%%%%
maxboxes = 5; % maximum number of boxes = 4^maxboxes
M = zeros(4^maxboxes,maxboxes+1);

currentsize = 0;

for i = 0:maxboxes
    currentsize = currentsize+1;
    counter = 0;

    d1 = pi.*ones(2^i,1);
    d2 = -pi.*ones((2^i)-1,1);
    A = diag(d1) + diag(d2,-1);

    areas = (pi/(2^i)).*ones(2^i,1);
    rhorange = [0; sqrt(A\areas)];

    thetarange = linspace(-pi,pi,(2^i)+1);

    for j = 1:length(thetarange)-1
        for k = 1:length(rhorange)-1
            counter = counter + 1;
            temp1 = (thetarange(j) <= thetaind) & (thetaind <= thetarange(j+1));
            temp2 = (rhorange(k) <= rhoind) & (rhoind <= rhorange(k+1));
            temp3 = temp1 & temp2;
            M(counter,currentsize) = sum(temp3); 
        end
    end
end

%%%%%% Calculate Spectrum %%%%%%%
q = qvals;
prbM = M./sum(M);%converting pixel counts to proportions

truesz = (2*pi)./(2.^(0:maxboxes));

% X = flip(log2(recsz));
X = log2(truesz);

yD = zeros(length(q),length(prbM(1,:)));
yalph = zeros(length(q),length(prbM(1,:)));
yf = zeros(length(q),length(prbM(1,:)));
mu = zeros(length(prbM(:,1)),length(prbM(1,:)),length(q));
yind = 0;


for k = q %Calculates y values
    yind = yind + 1;
    for a = 1:length(prbM(1,:))
        if k == 1
            yD(yind,a) = sum(nonzeros(prbM(:,a)).*log2(nonzeros(prbM(:,a))));
        else
            yD(yind,a) = log2(sum(nonzeros(prbM(:,a)).^k));
            mu(1:length(nonzeros(prbM(:,a))),a,yind) = (nonzeros(prbM(:,a)).^k)./(sum(nonzeros(prbM(:,a)).^k));
            yalph(yind,a) = sum(nonzeros(mu(:,a,yind)).*log2(nonzeros(prbM(:,a))));
            yf(yind,a) = sum(nonzeros(mu(:,a,yind)).*log2(nonzeros(mu(:,a,yind))));
        end
    end
end

Dq = zeros(length(q),1);
tauq = zeros(length(q),1);
myalpha = zeros(length(q),1);
falpha = zeros(length(q),1);

for currq = 1:length(q)
    if q(currq) == 1
        Dq(currq) = abs(fit(X',yD(currq,:)','poly1').p1);
    else
        tauq(currq) = (fit(X',yD(currq,:)','poly1').p1); %Calculating linear fits
        Dq(currq) = tauq(currq)./(q(currq)-1);
        myalpha(currq) = fit(X',yalph(currq,:)','poly1').p1;
        falpha(currq) = fit(X',yf(currq,:)','poly1').p1;
    end
end

h = q(end)-q(end-1);

alphaleg = zeros(length(tauq),1);
alphaleg(1) = (tauq(2) - tauq(1))/h;
alphaleg(end) = (tauq(end) - tauq(end-1))/h;

for step = 2:length(alphaleg)-1
    alphaleg(step) = (tauq(step+1) - tauq(step-1))/(2*h);
end

fleg = q'.*alphaleg - tauq;

myalpha = alphaleg;
falpha = fleg;

%%%%%%% Plots %%%%%%%
if plots
    figure
    plot(q,Dq,LineWidth=1.25)
    %ylim([0 2])
    xlabel('q')
    ylabel('Dq')
    title('Dq vs q')
    
    figure
    scatter(myalpha,falpha,'.b')
    %xlim([0.6 1.8])
    %ylim([0 1])
    xlabel('alpha')
    ylabel('f(alpha)')
    title('Spectrum')
end


