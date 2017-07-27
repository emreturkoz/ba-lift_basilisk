figure(1)
clf
hold all
for i=20:40:260
a=sprintf('jet_%d.dat',i);%,l);
          dat = load(a);
          y=dat(:,1);
          x=dat(:,2);
         
         
          plot(x,-y,'.')
         
end
axis([0 1 -3 0])
pbaspect([1 3 1])
% 
% for i=0:0.1:1
%     
%     fx = (H0.*(1.-((r./R0).^2)).^1.25)*(2/pi)*atan(i/0.5) +0.08;
%     plot(r,fx)
% end