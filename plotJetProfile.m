figure(1)
clf
hold all

low_limit = 40;
increment = 40;
high_limit = 800;
subplot_count = (high_limit - low_limit)/increment + 1;

k = 1;



subplot_row_count = 1;
subplot_col_count = subplot_count;

if(mod(subplot_count,2) == 0 )
    subplot_row_count = 2;
    subplot_col_count = subplot_count/2;
end

length_scale = 24.6057; % in microns
solid_thickness = 3.0; % in microns
axial_radial_ratio = 8.0;

for i=low_limit:increment:high_limit
a=sprintf('twolayeraxi16/jet_%d.dat',i);%,l);
          dat = load(a);
          y=dat(:,1) - (solid_thickness/length_scale);
          x=dat(:,2);
         
          subplot(subplot_row_count,subplot_col_count,k)
          plot(length_scale.*x,-length_scale.*y,'b.')
          hold on
          plot(-length_scale.*x,-length_scale.*y,'b.')
          k = k+ 1;
          axis([-length_scale length_scale -axial_radial_ratio*length_scale 0])
          pbaspect([1 axial_radial_ratio/2 1])
          title_time = sprintf('t = %dτ_b',i);
          title(title_time)

end
% 
% for i=0:0.1:1
%     
%     fx = (H0.*(1.-((r./R0).^2)).^1.25)*(2/pi)*atan(i/0.5) +0.08;
%     plot(r,fx)
% end