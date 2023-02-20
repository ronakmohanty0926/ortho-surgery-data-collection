function [newCurve, newLabels] = resample(curve,label,nPts)
% CURVSPACE Evenly spaced points along an existing curve in 2D or 3D.
%   CURVSPACE(P,N) generates N points that interpolates a curve
%   (represented by a set of points) with an equal spacing. Each
%   row of P defines a point, which means that P should be a n x 2
%   (2D) or a n x 3 (3D) matrix.
%
%   (Example)
%   x = -2*pi:0.5:2*pi;
%   y = 10*sin(x);
%   z = linspace(0,10,length(x));
%   N = 50;
%   p = [x',y',z'];
%   q = curvspace(p,N);
%   figure;
%   plot3(p(:,1),p(:,2),p(:,3),'*b',q(:,1),q(:,2),q(:,3),'.r');
%   axis equal;
%   legend('Original Points','Interpolated Points');
%
%   See also LINSPACE.
%
%   22 Mar 2005, Yo Fukushima
%% initial settings %%
currentpt = curve(1,:); % current point
currentlbl = label(1,:);
indfirst = 2; % index of the most closest point in p from curpt
indlblfirst = 2;
len = size(curve,1); % length of p
newCurve = currentpt; % output point
newLabels = currentlbl;
k = 0;
%% distance between points in p %%
for k0 = 1:len-1
   dist_bet_pts(k0) = distance(curve(k0,:),curve(k0+1,:));
end
totaldist = sum(dist_bet_pts);
%% interval %%
intv = totaldist./(nPts-1);
%% iteration %%
for k = 1:nPts-1
   
   newpt = []; newlbl=[]; distsum = 0;
   ptnow = currentpt;
   lblnow = currentlbl;
   kk = 0;
   pttarget = curve(indfirst,:);
   lbltarget = label(indlblfirst,:);
   remainder = intv; % remainder of distance that should be accumulated
%    count = 0;
   while isempty(newpt)
      % calculate the distance from active point to the most
      % closest point in p
      disttmp = distance(ptnow,pttarget);
      distsum = distsum + disttmp;
      % if distance is enough, generate newpt. else, accumulate
      % distance
      if distsum >= intv
%          count = count + 1;
         newpt = interpintv(ptnow,pttarget,remainder);
         if (lblnow == lbltarget) 
             newlbl = lbltarget;
         elseif (lblnow ~= lbltarget)
             newlbl = lbltarget;
         end
      else
         remainder = remainder - disttmp;
         ptnow = pttarget;
         kk = kk + 1;
         if indfirst+kk > len
            newpt = curve(len,:);
            newlbl = label(len,1);
         else
            pttarget = curve(indfirst+kk,:);
            lbltarget = label(indlblfirst+kk,1);
%             if (label(k,1) == label(k+1,1)) newlbl = label(k,1);
%             elseif (label(k,1) ~= label(k+1,1))newlbl = label(count+1,1);
%             end
         end
      end
   end
   
   % add to the output points
   newCurve = [newCurve; newpt];
   newLabels = [newLabels; newlbl];
   
   % update currentpt and indfirst
   currentpt = newpt;
   currentlbl = newlbl;
   indfirst = indfirst + kk;
   indlblfirst = indlblfirst + kk;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SUBFUNCTIONS     %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function l = distance(x,y)
% DISTANCE Calculate the distance.
%   DISTANCE(X,Y) calculates the distance between two
%   points X and Y. X should be a 1 x 2 (2D) or a 1 x 3 (3D)
%   vector. Y should be n x 2 matrix (for 2D), or n x 3 matrix
%   (for 3D), where n is the number of points. When n > 1,
%   distance between X and all the points in Y are returned.
%
%   (Example)
%   x = [1 1 1];
%   y = [1+sqrt(3) 2 1];
%   l = distance(x,y)
%
% 11 Mar 2005, Yo Fukushima
%% calculate distance %%
if size(x,2) == 2
   l = sqrt((x(1)-y(:,1)).^2+(x(2)-y(:,2)).^2);
elseif size(x,2) == 3
   l = sqrt((x(1)-y(:,1)).^2+(x(2)-y(:,2)).^2+(x(3)-y(:,3)).^2);
else
   error('Number of dimensions should be 2 or 3.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newpt = interpintv(pt1,pt2,intv)
% Generate a point between pt1 and pt2 in such a way that
% the distance between pt1 and new point is intv.
% pt1 and pt2 should be 1x3 or 1x2 vector.
dirvec = pt2 - pt1;
dirvec = dirvec./norm(dirvec);
l = dirvec(1); m = dirvec(2);
newpt = [intv*l+pt1(1),intv*m+pt1(2)];
if length(pt1) == 3
   n = dirvec(3);
   newpt = [newpt,intv*n+pt1(3)];
end