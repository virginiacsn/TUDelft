function [ellipse, xr, yr, phi, Ew, Kel] = ellipse_fit(x,y)

%This function will calculate an ellipse using PCA on the data set
%consisting of x and y coordinates. The Eigenvalues of the PCA are scaled
%by normalizing the dataset by the Eigenvalues andcalculating the 95%
%confidence interval of the langths of the vectors to the different data
%points

%Calculating means of x and y
xm = mean(x);
ym = mean(y);

%putting x and y in a vector (first y and then x)
Y = y - ym;
X = x - xm;

XY = [x, y];
%  = FAb(:,2);% - mean(FAb(:,2));
% y = FAb(:,3);% - mean(FAb(:,3));
% x =[-10 5 -5 10];
% y = [-10 -5 5 10];

%Optie1
%Construct M
M = [2*X.*Y Y.^2 2*X 2*Y ones(size(X))];
% M = [X.*Y Y.^2 X Y ones(size(X))];
% Multiply (-X.^2) by pseudoinverse(M)
e = M\(-X.^2);

% Optie 2
% D = [X.*X X.*Y Y.*Y X Y ones(size(X))]; % design matrix
% S = D' * D; % scatter matrix
% C(6, 6) = 0; C(1, 3) = 2; C(2, 2) = -1; C(3, 1) = 2; % constraint matrix
% [gevec, geval] = eig(inv(S) * C); % solve eigensystem
% [PosR, PosC] = find(geval > 0 & ~isinf(geval)); % find positive eigenvalue
% e = gevec(:, PosC); % corresponding eigenvector

% %Optie 3
% centroid = mean(XY);   % the centroid of the data set
% 
% D1 = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
%       (XY(:,2)-centroid(2)).^2];
% D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
% S1 = D1'*D1;
% S2 = D1'*D2;
% S3 = D2'*D2;
% T = -inv(S3)*S2';
% M = S1 + S2*T;
% M = [M(3,:)./2; -M(2,:); M(1,:)./2];
% [evec,eval] = eig(M);
% cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
% A1 = evec(:,find(cond>0));
% A = [A1; T*A1];
% A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
% A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
% A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
%      A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
% A(4) = A4;  A(5) = A5;  A(6) = A6;
% e = A/norm(A);
% disp(e)

%Optie 4 (Taubin)
% centroid = mean(XY);   % the centroid of the data set
% 
% Z = [(XY(:,1)-centroid(1)).^2, 2*(XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
%      (XY(:,2)-centroid(2)).^2, 2*XY(:,1)-centroid(1), 2*XY(:,2)-centroid(2), ones(size(XY,1),1)];
% M = Z'*Z/size(XY,1);
% 
% P = [M(1,1)-M(1,6)^2, M(1,2)-M(1,6)*M(2,6), M(1,3)-M(1,6)*M(3,6), M(1,4), M(1,5);
%      M(1,2)-M(1,6)*M(2,6), M(2,2)-M(2,6)^2, M(2,3)-M(2,6)*M(3,6), M(2,4), M(2,5);
%      M(1,3)-M(1,6)*M(3,6), M(2,3)-M(2,6)*M(3,6), M(3,3)-M(3,6)^2, M(3,4), M(3,5);
%      M(1,4), M(2,4), M(3,4), M(4,4), M(4,5);
%      M(1,5), M(2,5), M(3,5), M(4,5), M(5,5)];
% 
% Q = [4*M(1,6), 2*M(2,6), 0, 0, 0;
%      2*M(2,6), M(1,6)+M(3,6), 2*M(2,6), 0, 0;
%      0, 2*M(2,6), 4*M(3,6), 0, 0;
%      0, 0, 0, 1, 0;
%      0, 0, 0, 0, 1];
% 
% [V,D] = eig(P,Q);
% 
% [Dsort,ID] = sort(diag(D));
% 
% A = V(:,ID(1));
% A = [A; -A(1:3)'*M(1:3,6)];
% A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
% A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
% A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
%      A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
% A(4) = A4;  A(5) = A5;  A(6) = A6;
% e = A/norm(A);

%Extract parameters from vector e
a = 1;
b = e(1);
c = e(2);
d = e(3);
f = e(4);
g = e(5);

% a = e(1);
% b = e(2);
% c = e(3);
% d = e(4);
% f = e(5);
% g = e(6);
% % %Use Formulas from Mathworld to find semimajor_axis, semiminor_axis, x0, y0
%, and phi

% delta = b^2-4*a*c;
delta = b^2-a*c;

x0 = (c*d - b*f)/delta;
y0 = (a*f - b*d)/delta;

xr = xm + x0;
yr = ym + y0;

if a < c && b == 0
    an = 0;
%     disp('1')
elseif a > c && b == 0
    an = pi/2;
%     disp('2')
elseif a < c && b ~= 0
%     an = - 0.5 * acot((c-a)/(2*b));
    an = 0.5 * acot((a-c)/(2*b));
%     disp('3')
elseif a > c && b ~= 0
%     an = pi/2 - 0.5 * acot((c-a)/(2*b));
an = pi/2 + 0.5 * acot((a-c)/(2*b));
%     disp('4')
end
% phi =  pi/2 - 0.5 * acot((c-a)/(2*b));

if an < 0
    phi = an + pi;
else
    phi = an;
end

nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
s = sqrt(1 + (4*b^2)/(a-c)^2);

a_prime = sqrt(nom/(delta* ( (c-a)*s -(c+a))));

b_prime = sqrt(nom/(delta* ( (a-c)*s -(c+a))));

% a_prime = sqrt(nom/(delta* sqrt((a-c)^2 +4*b^2) - (a+c)));
% b_prime = sqrt(nom/(delta* -sqrt((a-c)^2 +4*b^2) - (a+c)));

semimajor_axis = max(a_prime, b_prime);
semiminor_axis = min(a_prime, b_prime);

Ew = [semimajor_axis; semiminor_axis];

n = 0:0.0175:2*pi;
Kel = [cos(phi), -sin(phi); sin(phi), cos(phi)];

% ellipse = [cos(phi), -sin(phi); sin(phi), cos(phi)]*[semimajor_axis*cos(n);semiminor_axis*sin(n)];
ellipse = [cos(phi), -sin(phi); sin(phi), cos(phi)]*[semimajor_axis*cos(n);semiminor_axis*sin(n)];

% if (a_prime < b_prime)
%     phi = pi/2 - phi;
% end
% n = 0:0.001:2*pi;
% ellipse = [cos(phi), -sin(phi); sin(phi), cos(phi)]*[semimajor_axis*cos(n);semiminor_axis*sin(n)];
% figure(q)
% plot(ellipse(1,:) + mean(FAb(:,2)),ellipse(2,:) + mean(FAb(:,3)),'m'), hold on
% plot(mean(FAb(:,2)),mean(FAb(:,3)),'m+','markersize',15)
% 
% 
% 
% xm = mean(x);
% ym = mean(y);
% 
% %putting x and y in a vector (first y and then x)
% X = [y - ym, x - xm];
% 
% %Calculating Eigen Vectors (EV) and Eigen values (Ew)
% [EV,score,Ew,tsquared,explained] = pca(X);
% 
% %Using the eigenvector of the largest eigenvalue (always first) to
% %calculate the orientation of the ellipse
% % phi = atan((EV(1,1)/EV(2,1)));
% 
% an = atan((EV(1,1)/EV(2,1)));
% 
% if an < 0
%     phi = an + pi;
% else
%     phi = an;
% end
% 
% %Rotating the data point using the EV. Now the Ev are the x and y axis
% Xrot = EV' * X';
% 
% %Normalizing the coordinates of the datapoints to the Ew
% xn = Xrot(1,:) ./ max(Ew);
% yn = Xrot(2,:) ./ min(Ew);
% xn = xn.';
% yn = yn.';
% 
% %calculating the length of the vectors
% Ftot = sqrt(xn.^2 + yn.^2);
% 
% %Taking the 95% CI of the length of the normalized vectors and adding the
% %mean, to get the scaling vector
% S = 1.96*std(Ftot)/sqrt(length(Ftot)) + mean(Ftot);
% 
% %Calculating the ellipse
% n = 0:0.0175:2*pi;
% Kel = [cos(phi), -sin(phi); sin(phi), cos(phi)];
% 
% ellipse = [cos(phi), -sin(phi); sin(phi), cos(phi)]*[S*max(Ew)*cos(n);S*min(Ew)*sin(n)];
