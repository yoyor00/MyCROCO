%======================================================================
%
%     ---               Dune Test Case                ---     
%
%
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
% 
%  Patrick Marchesiello - 2012
%======================================================================
clear all
close all
%================== User defined parameters ===========================

fname='dune_his.nc';

%======================================================================
%
% Read data
%
j=2; 
% Remember you are plotting for time (in days output) = tndx - 1 
% For instance tndx=20 for 19 days simulated
tndx=20; 
%

nc=netcdf(fname);
% assuming two classes of sand, the second the finer
sandname = ncreadatt(fname,'/','Sd sand_2')*1e6; %um
time=nc{'scrum_time'}(:)./86400;
tndx=min(tndx,length(time));
disp([ 'tndx = ',num2str(tndx), ...
    ' - Time = ',num2str(time(tndx)),' days' ])
h=squeeze(nc{'h'}(j,:));
x=squeeze(nc{'x_rho'}(j,:));
zeta=squeeze(nc{'zeta'}(tndx,j,:));
L=length(zeta);
%
hmorph  =squeeze(nc{'hmorph'}(tndx,j,:));
bedthick=squeeze(nc{'bed_thick'}(tndx,:,j,:));
bedfrac =squeeze(nc{'bed_frac_sand_2'}(tndx,:,j,:));
NL=size(bedfrac,1);
%
u=squeeze(nc{'u'}(tndx,:,j,:));
N=size(u,1);
%
theta_s=nc.theta_s(:);
theta_b=nc.theta_b(:);
hc=nc.hc(:);
Vtrans=nc{'Vtransform'}(:);
close(nc);
%
zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',Vtrans);
zr=squeeze(zr);
x =reshape(x,1,L);
xr=repmat(x,[N 1]);
xu=0.5*(xr(:,1:end-1)+xr(:,2:end));
zu=0.5*(zr(:,1:end-1)+zr(:,2:end));
%
NL=NL+1;
bedthick=padarray(bedthick,1,'replicate','pre');
bedfrac =padarray(bedfrac, 1,'replicate','pre');
bedfrac(2:end-1,:)=bedfrac(3:end,:);
xbed=repmat(x,[NL 1]);
h2=repmat(hmorph,[NL 1]);
zbed=-h2;
zbed(2:end,:)=-h2(2:end,:)-cumsum(bedthick(2:end,:),1);
%
%  Plot bed stratigraphy
%
figure('position',[600 600 700 700])
subplot(2,1,1)
%subplot(3,1,1)
%hold on
%pcolor(xbed,zbed,100*bedfrac);
%%Cdata = test.CData;
%line(xr,-hmorph,'color','k','Linewidth',2)
%line(xr,-h,     'color','r','Linewidth',2,'linestyle','--')
%hold off
%caxis([40 60])
%%caxis([0 100])
%colorbar
%%colorbar('Ticks',[0,25,50,75,100])
%axis([-Inf Inf -8 0])
%ylabel('Depth [m]','Fontsize',15)
%title({'DUNE Test Case' [' Fine Sand Fraction ',num2str(sandname),'um - Day ', ...
%        num2str(time(tndx))] })
%%set(gca,'fontsize',15);
%set(gcf,'PaperPositionMode','auto');
%hold off
%
%%
%%  Plot bed stratigraphy different VERSION only boxes
%%
%subplot(3,1,2)
%hold on

plot(x,-hmorph,'color','k','Linewidth',2)
axis([-Inf Inf -8 0])

[foo, Sx]=size(x);
for ik=1:Sx-1
  for jk=1:NL-1
   line([ x(ik)-1 x(ik)-1]      ,[ zbed(jk,ik) zbed(jk+1,ik)   ],'Color','k');
   line([ x(ik)-1 x(ik+1)-1 ]   ,[ zbed(jk,ik) zbed(jk,ik)     ],'Color','k');
   line([ x(ik+1)-1  x(ik+1)-1 ],[ zbed(jk,ik) zbed(jk+1,ik)   ],'Color','k');
   line([ x(ik)-1 x(ik+1)-1 ]   ,[ zbed(jk+1,ik) zbed(jk+1,ik) ],'Color','k');
   patch( [ x(ik)-1  x(ik+1)-1  x(ik+1)-1 x(ik)-1 ] ,...
          [ zbed(jk,ik) zbed(jk,ik) zbed(jk+1,ik) zbed(jk+1,ik) ],...
          [ bedfrac(jk,ik)*100 bedfrac(jk,ik)*100 bedfrac(jk,ik)*100 bedfrac(jk,ik)*100 ])
  end
end
line(x,-h,'color','r','Linewidth',2,'linestyle','--')

hold off 
%caxis([0 100])
caxis([40 60])
colorbar
%colorbar('Ticks',[0,25,50,75,100]) 
ylabel('Depth [m]','Fontsize',15)
title({'DUNE Test Case' [' Fine Sand Fraction ',num2str(sandname),'um - Day ', ...
        num2str(time(tndx))] })
%%set(gca,'fontsize',15);
%set(gcf,'PaperPositionMode','auto');


hold off

%
%  Plot evolution bed
%
nc=netcdf(fname);
hmorph=squeeze(nc{'hmorph'}(:,j,:));
close(nc)
subplot(2,1,2)
line(xr,-hmorph( 1,:),'color','k','Linewidth',2) %topo init
hold on
for n=1:4
 it=2^(n-1)+1;
 if it<tndx
  line(xr,-hmorph(it,:),'color','k','Linewidth',2)
 end
end
line(xr,-hmorph(tndx,:),'color','r','Linewidth',2) %topo at tndx
hold off
grid on
axis([-Inf Inf -4.5 0])
title({'DUNE Test Case'  ['Bed evolution: 0 1 2 4 8 ',num2str(time(tndx)),' days']}) 
xlabel('Distance [m]','Fontsize',15)
ylabel('Depth [m]','Fontsize',15)

export_fig dune.png

return

