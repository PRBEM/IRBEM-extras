function color = getcolor(i,cN)
% function color = getcolor(i,cN)
if (cN < 7),
   colors = [1 0 0
      0   0.8 0
      0   0   1
      0.7 0   0.7
      0   0.6 0.7
      0.7 0.6 0
   ];
   color = colors(i,:);
else
   j = (i-1)/(cN-1);
   g = 1-2*abs(j-0.5);
   if j < 0.5,
       r = 0;
       b = 2*abs(j-0.5);
   else
       b = 0;
       r = 2*abs(j-0.5);
   end       
   color = [r g b];
end

%    j = (i-1)/(cN-1);
%    a = rem(j*3,1);
%    b = rem(j*9,1);
%    color = 1/4+2/3*[j,1-b,a];
