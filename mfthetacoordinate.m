% calculates multifractal spectrum for 2d binary image
% calculates using theta angles mapped to a line segment
% objectcolor = color of object of interest ** 0 for black ** ** 1 for white **
% qvals = vector of desired q values
% outputs dq vs q plot and spectrum
% faster (and more accurate) than mfellipsebinarizedv2

function [Dq,myalpha,falpha] = mfthetacoordinate(bwimage,objectcolor,qvals,plots)

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

%Yim = (centerx/centery).*Yim;

[theta,rho] = cart2pol(Xim,Yim); %polarcoordinates of all points

myind = mylog == 1; %indices of points of interest
thetaind = theta(myind);

%%%%%%% Box Counting Loop %%%%%%%%
numslices = 2.^(0:6); 
sectang = 2.*pi./numslices;
M = zeros(max(numslices),length(numslices));

temptheta = linspace(-pi/2,pi/2,max(numslices)/2 + 1);

scaledtheta = atan2((centerx/centery).*tan(temptheta),1);
mydiff = diff(scaledtheta);

thetarange = zeros(max(numslices) + 1,1);
thetarange(1) = -pi;
thetarange(end) = pi;

for ii = 2:length(mydiff)
    thetarange(ii) = thetarange(ii-1) + mydiff(ii-1);
    thetarange(end-ii+1) = thetarange(end-ii+2) - mydiff(ii-1);
end

for j = 1:length(thetarange)-1
        temp1 = (thetarange(j) <= thetaind) & (thetaind < thetarange(j+1));
        M(j,1) = sum(temp1); 
end

for col = 2:length(M(1,:))
    a = 2^(col-1);
    for ii = 1:a:length(M(:,1))
        M(ii,col) = sum(M(linspace(ii,ii+a-1,a),col-1));
    end
end

%%%%%% Calculate Spectrum %%%%%%%
%% Set up regression
q = qvals; %vector of q values
prbM = M./sum(M(:,1));%converting pixel counts to proportions
X = log2(sectang);

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
        Dq(currq) = tauq(currq)./(1-q(currq));
        myalpha(currq) = fit(X',yalph(currq,:)','poly1').p1;
        falpha(currq) = fit(X',yf(currq,:)','poly1').p1;
    end
end

h = 0.1;

alphaleg = zeros(length(tauq),1);
alphaleg(1) = (tauq(2) - tauq(1))/h;
alphaleg(end) = (tauq(end) - tauq(end-1))/h;

for step = 2:length(alphaleg)-1
    alphaleg(step) = (tauq(step+1) - tauq(step-1))/(2*h);
end

fleg = q'.*alphaleg - tauq;

myalpha = -alphaleg;
falpha = -fleg;

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
