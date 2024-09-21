% calculates multifractal spectrum for 2d binary image inputted as a
% logical array
% calculates on rectangular support
% objectcolor = color of object of interest ** 0 for black ** ** 1 for white **
% qvals = vector of desired q values
% outputs dq vs q plot and spectrum
% plots := 1 for yes 0 for no

function [Dq,myalpha,falpha] = mfrectanglebinarized(bwmatrix,objectcolor,qvals,plots)

sz = size(bwmatrix);

if objectcolor
    nagalog = bwmatrix == 1;
else
    nagalog = bwmatrix == 0;
end

%imshow(nagalog)
%% Box Counting 
recsz = 2.^(0:1:4); %vector of rectangle sizes

M = zeros(max(recsz)^2,length(recsz));

currentsize = 0;

for x = recsz % Rectangular box code
    counter = 0;
    currentsize = currentsize + 1;
    pxsz = floor(sz./x);

    numboxR = floor(sz(1) / pxsz(1));
    rowvec = [pxsz(1) * ones(1, numboxR), rem(sz(1), pxsz(1))];

    numboxC = floor(sz(2) / pxsz(2));
    colvec = [pxsz(2) * ones(1, numboxC), rem(sz(2), pxsz(2))];

    logcell = mat2cell(nagalog,rowvec,colvec);

    for u = 1:numboxR
        for v = 1:numboxC
            counter = counter + 1;
            temp = cell2mat(logcell(u,v));
            M(counter,currentsize) = sum(sum(temp));
        end
    end

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = qvals;
prbM = M./sum(M);%converting pixel counts to proportions

truesz = sz(1)./recsz;

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Dq vs q and Spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%
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